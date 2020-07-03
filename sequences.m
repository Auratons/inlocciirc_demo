addpath('functions/InLocCIIRC_utils/params');
addpath('functions/InLocCIIRC_utils/load_CIIRC_transformation');
addpath('functions/InLocCIIRC_utils/mkdirIfNonExistent');
addpath('functions/InLocCIIRC_utils/multiCameraPose');
addpath('functions/InLocCIIRC_utils/projectPointCloud');
addpath('functions/InLocCIIRC_utils/projectPointsUsingP');
addpath('functions/InLocCIIRC_utils/rotationDistance');
addpath('functions/InLocCIIRC_utils/R_to_numpy_array');
addpath('functions/InLocCIIRC_utils/T_to_numpy_array');
[ params ] = setupParams('holoLens1'); % NOTE: adjust

startIdx = 127; % the index of the first query to be considered in the sequence
k = 5; % the length of the sequence

%% extract HoloLens poses wrt initial unknown HoloLens CS
descriptionsTable = readtable(params.queryDescriptions.path); % decribes the reference poses

prevWarningState = warning();
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
rawHoloLensPosesTable = readtable(params.holoLens.poses.path);
warning(prevWarningState);

assert(size(descriptionsTable,1) == size(rawHoloLensPosesTable,1));
nQueries = size(descriptionsTable,1);

cameraPosesWrtHoloLensCS = zeros(nQueries,4,4);
queryInd = zeros(nQueries,1);
for i=1:nQueries
    queryId = descriptionsTable{i, 'id'};
    queryInd(i,:) = queryId;
    %space = descriptionsTable{i, 'space'}{1,1};
    %inMap = descriptionsTable{i, 'inMap'};
    t = [rawHoloLensPosesTable{i, 'Position_X'}; ...
                rawHoloLensPosesTable{i, 'Position_Y'}; ...
                rawHoloLensPosesTable{i, 'Position_Z'}];
    orientation = [rawHoloLensPosesTable{i, 'Orientation_W'}, ...
                    rawHoloLensPosesTable{i, 'Orientation_X'}, ...
                    rawHoloLensPosesTable{i, 'Orientation_Y'}, ...
                    rawHoloLensPosesTable{i, 'Orientation_Z'}];
    R = rotmat(quaternion(orientation), 'frame'); % what are the columns of R? 
        % Bases of WHAT wrt WHAT? (one of them is initial unknown HL CS, the other is HL camera CS)
        % -> it is most likely a rotation matrix from initial unknown HL CS to HL camera CS. i.e. the columns
        % are bases of initial unknown HL CS in HL camera CS coordinates

    % camera points to -z in HoloLens
    % see https://docs.microsoft.com/en-us/windows/mixed-reality/coordinate-systems-in-directx
    rFix = rotationMatrix([pi, 0.0, 0.0], 'ZYX'); % correct
    %rFix = eye(3); % incorrect

    %R = R';

    R1 = (rFix * R)';
    R2 = R' * rFix;
    Rd = R1-R2;
    eps = 1e-8;
    assert(all(Rd(:) < eps));

    cameraPositionWrtHoloLensCS = t';
    cameraOrientationWrtHoloLensCS = R2;

    pose = eye(4);
    pose(1:3,1:3) = cameraOrientationWrtHoloLensCS;
    pose(1:3,4) = cameraPositionWrtHoloLensCS;
    cameraPosesWrtHoloLensCS(i,:,:) = pose;
end

%% include only those in the sequence
queryInd = queryInd(startIdx:startIdx+k-1);
cameraPosesWrtHoloLensCS2 = zeros(k,4,4); % accounted for (possible) delay
% TODO: assert all query IDs are sorted and increasing by 1
for i=1:k
    queryId = queryInd(i);
    orientationDataIdx = queryId+params.HoloLensOrientationDelay;
    translationDataIdx = queryId+params.HoloLensTranslationDelay;
    %orientationDataIdx = queryId;
    %translationDataIdx = queryId;
    if (orientationDataIdx > nQueries || translationDataIdx > nQueries)
        error('No HoloLens pose data for query %d', queryId);
    end
    pose = eye(4);
    pose(1:3,1:3) = cameraPosesWrtHoloLensCS(orientationDataIdx,1:3,1:3);
    pose(1:3,4) = cameraPosesWrtHoloLensCS(translationDataIdx,1:3,4);
    cameraPosesWrtHoloLensCS2(i,:,:) = pose;
end
cameraPosesWrtHoloLensCS = cameraPosesWrtHoloLensCS2;

%% debug - we need to print: origin of each camera wrt Omega, bases of each camera wrt Omega,
for i=1:k
    queryId = queryInd(i);
    fprintf('query: %d\n', queryId);
    fprintf('camera origin wrt Omega: %s\n', T_to_numpy_array(cameraPosesWrtHoloLensCS(i,1:3,4)));
    fprintf('camera bases wrt Omega: %s\n', R_to_numpy_array(squeeze(cameraPosesWrtHoloLensCS(i,1:3,1:3))));
end

%% set up 2D-3D correspondences for the k queries
sensorSize = params.camera.sensor.size; % height, width
imageWidth = sensorSize(2);
imageHeight = sensorSize(1);
correspondences2D = [imageWidth/4, imageHeight/4; ...
                    imageWidth/4, imageHeight-imageHeight/4; ...
                    imageWidth-imageWidth/4, imageHeight/4; ...
                    imageWidth-imageWidth/4, imageHeight-imageHeight/4]';
correspondences3D = zeros(k,3,4);

j = 1;
for i=startIdx:startIdx+k-1
    queryId = queryInd(j);
    retrievedPosePath = fullfile(params.evaluation.retrieved.poses.dir, sprintf('%d.txt', queryId));
    retrievedPose = load_CIIRC_transformation(retrievedPosePath);
    retrievedT = -inv(retrievedPose(1:3,1:3))*retrievedPose(1:3,4); % wrt model
    retrievedR = retrievedPose(1:3,1:3); % modelBasesToEpsilonBases
    P = [params.camera.K*retrievedR, -params.camera.K*retrievedR*retrievedT]; 
    imageToModel = inv([P; 0,0,0,1]);
    toDeproject = ones(4,4);
    toDeproject(1:2,:) = correspondences2D;
    correspondences4D = imageToModel * toDeproject;
    correspondences3D(j,:,:) = correspondences4D(1:3,:);
    j = j + 1;
end

%% verification of correspondences - it works
%close all;
%figure;
%pointSize = 8.0;
%outputSize = params.camera.sensor.size;
%projectedPointCloud = projectPointCloud(params.pointCloud.path, params.camera.fl, retrievedR, ...
%                                    retrievedT, params.camera.sensor.size, outputSize, pointSize, ...
%                                    params.projectPointCloudPy.path); % TODO: use projectMesh instead, which can work in headless mode
%image(projectedPointCloud);
%axis image;
%
%hold on;
%scatter(correspondences2D(1,:), correspondences2D(2,:), 40, 'r', 'filled');
%reprojectedPts = projectPointsUsingP(squeeze(correspondences3D(end,:,:)), P);
%scatter(reprojectedPts(1,:), reprojectedPts(2,:), 20, 'g', 'filled');
%hold off;
%set(gcf, 'Position', get(0, 'Screensize'));

%% execute and collect results
%workingDir = tempname;
workingDir = '/Volumes/GoogleDrive/Můj disk/ARTwin/InLocCIIRC_dataset/evaluation/sequences'; % only for debugging; TODO: remove
    
inlierThreshold = 12.0; % TODO
numLoSteps = 10; % TODO
invertYZ = false; % TODO
pointsCentered = false;
undistortionNeeded = false; % TODO
estimatedPoses = multiCameraPose(workingDir, queryInd, cameraPosesWrtHoloLensCS, ...
                                    correspondences2D, correspondences3D, ...
                                    inlierThreshold, numLoSteps, ...
                                    invertYZ, pointsCentered, undistortionNeeded, params); % wrt model
mkdirIfNonExistent(params.evaluation.sequences.dir);

%% compare poses estimated by MultiCameraPose with reference poses
%% quantitative results
errors = struct();
for i=1:k
    queryId = queryInd(i);
    queryPoseFilename = sprintf('%d.txt', queryId);
    posePath = fullfile(params.dataset.query.dir, 'poses', queryPoseFilename);
    referenceP = load_CIIRC_transformation(posePath);
    referenceT = -inv(referenceP(1:3,1:3))*referenceP(1:3,4);
    referenceR = referenceP(1:3,1:3);

    estimatedT = estimatedPoses(i,1:3,4)';
    estimatedR = squeeze(estimatedPoses(i,1:3,1:3));

    errors(i).queryId = queryId;
    errors(i).translation = norm(estimatedT - referenceT);
    errors(i).orientation = rotationDistance(referenceR, estimatedR);
end
errorsTable = struct2table(errors);
errors = table2struct(sortrows(errorsTable, 'queryId'));
fprintf('id\ttranslation [m]\torientation [deg]\n');
for i=1:k
    fprintf('%d\t%0.2f\t%0.2f\n', errors(i).queryId, errors(i).translation, errors(i).orientation);
end
fprintf('Mean\t%0.2f\t%0.2f\n', mean([errors.translation]),  mean([errors.orientation]));
return;

%% qualitative results
close all;
for i=1:k
    figure;
    pointSize = 8.0;
    outputSize = params.camera.sensor.size;
    estimatedT = estimatedPoses(i,1:3,4)';
    estimatedR = squeeze(estimatedPoses(i,1:3,1:3));
    projectedPointCloud = projectPointCloud(params.pointCloud.path, params.camera.fl, estimatedR, ...
                                        estimatedT, params.camera.sensor.size, outputSize, pointSize, ...
                                        params.projectPointCloudPy.path);
    image(projectedPointCloud);
    axis image;
end
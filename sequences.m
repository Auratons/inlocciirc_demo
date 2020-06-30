addpath('functions/InLocCIIRC_utils/params');
addpath('functions/InLocCIIRC_utils/load_CIIRC_transformation');
addpath('functions/InLocCIIRC_utils/mkdirIfNonExistent');
addpath('functions/InLocCIIRC_utils/multiCameraPose');
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

cameraPoseWrtHoloLensCS = zeros(nQueries,4,4);
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
    rFix = rotationMatrix([pi, 0.0, 0.0], 'ZYX');
    R = rFix * R;

    cameraPositionWrtHoloLensCS = t';
    cameraOrientationWrtHoloLensCS = R';

    pose = eye(4);
    pose(1:3,1:3) = cameraOrientationWrtHoloLensCS;
    pose(1:3,4) = cameraPositionWrtHoloLensCS;
    cameraPoseWrtHoloLensCS(i,:,:) = pose;
end
queryInd = queryInd(startIdx:startIdx+k-1); % include only those in the sequence

%% convert k HL poses to be wrt 1st HL camera pose CS
cameraPoseWrtFirstCameraPose = zeros(k,4,4);
holoLensCSToFirstCameraPoseCS = inv(squeeze(cameraPoseWrtHoloLensCS(startIdx,:,:)));
j = 1;
for i=startIdx:startIdx+k-1
    cameraPoseWrtFirstCameraPose(j,:,:) = holoLensCSToFirstCameraPoseCS * squeeze(cameraPoseWrtHoloLensCS(i,:,:));
    j = j + 1;
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
    queryId = descriptionsTable{i, 'id'};
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

%% execute and collect results
%workingDir = tempname;
workingDir = '/Volumes/GoogleDrive/Můj disk/ARTwin/InLocCIIRC_dataset/evaluation/sequences'; % only for debugging; TODO: remove

inlierThreshold = 12.0; % TODO
numLoSteps = 10; % TODO
invertYZ = false; % TODO
pointsCentered = false;
undistortionNeeded = false; % TODO
estimatedPoses = multiCameraPose(workingDir, queryInd, cameraPoseWrtFirstCameraPose, ...
                                    correspondences2D, correspondences3D, ...
                                    inlierThreshold, numLoSteps, ...
                                    invertYZ, pointsCentered, undistortionNeeded, params); % wrt model
mkdirIfNonExistent(params.evaluation.sequences.dir);
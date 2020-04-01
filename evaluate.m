%% initialization
load(params.input.qlist_matname, 'query_imgnames_all');
densePV_matname = fullfile(params.output.dir, 'densePV_top10_shortlist.mat');
load(densePV_matname, 'ImgList');

if exist(params.evaluation.dir, 'dir') ~= 7
    mkdir(params.evaluation.dir)
end

%% visual evaluation
if exist(params.evaluation.query_vs_synth.dir, 'dir') ~= 7
    mkdir(params.evaluation.query_vs_synth.dir);
end

for i=1:size(query_imgnames_all,2)
    queryName = query_imgnames_all{i};
    queryImage = imread(fullfile(params.data.dir, params.data.q.dir, queryName));
    
    fun = @(x) strcmp(ImgList(x).queryname,queryName);
    tf = arrayfun(fun, 1:numel(ImgList));
    ImgListRecord = ImgList(find(tf));
    cutoutPath = ImgListRecord.topNname{1};
    synthName = strsplit(cutoutPath, '/');
    synthName = synthName{end};
    synthName = synthName(1:end-size(params.data.db.cutout.imgformat,2));
    synthPath = fullfile(params.output.synth.dir, queryName, [synthName, params.output.synth.matformat]);
    load(synthPath, 'RGBpersp');
    numRows = size(queryImage,1);
    numCols = size(queryImage,2);
    
    synthImage = RGBpersp;
    if isempty(synthImage)
        synthImage = zeros(numRows, numCols, 3, 'uint8');
    else
        synthImage = imresize(synthImage, [numRows numCols]);
    end

    imshowpair(queryImage, synthImage, 'montage');
    saveas(gcf, fullfile(params.evaluation.query_vs_synth.dir, queryName));
end

%% quantitative results
nQueries = size(query_imgnames_all,2);
errors = struct();
retrievedPoses = struct();
for i=1:nQueries
    queryName = query_imgnames_all{i};
    queryImage = imread(fullfile(params.data.dir, params.data.q.dir, queryName));
    
    fun = @(x) strcmp(ImgList(x).queryname,queryName);
    tf = arrayfun(fun, 1:numel(ImgList));
    ImgListRecord = ImgList(tf);
    
    cutoutPath = ImgListRecord.topNname{1};
    cutoutPath = strsplit(cutoutPath, '/');
    spaceName = cutoutPath{1};
    sweepId = cutoutPath{2};
    transPath = fullfile(params.data.dir, params.data.db.trans.dir, spaceName, 'transformations', sprintf('trans_%s.txt', sweepId));
    P1 = load_CIIRC_transformation(transPath);
    R1 = P1(1:3,1:3);
    T1 = P1(1:3,4);
    
    P2 = ImgListRecord.P{1};
    if any(isnan(P2(:)))
        T = nan(3,1);
        orientation = nan(3,1);
    else
        R2 = P2(1:3,1:3); % in fact this is K*R - are you sure? viz demo
        T2 = P2(1:3,4);

        P = P2*P1;
        T = -inv(P(:,1:3))*P(:,4);
        initialDirection = [0.0; 0.0; -1.0];
        rFix = [180.0, 0.0, 0.0];
        Rfix = rotationMatrix(deg2rad(rFix), 'XYZ');
        orientation = R2 * (R1 * (Rfix * initialDirection));
        % very ugly magic hack
        if orientation(3) < 0
            orientation = orientation .* [1.0; -1.0; 1.0];
        end
    end
    
    posesPath = fullfile(params.data.dir, params.data.q.dir, 'poses.csv');
    posesTable = readtable(posesPath);
    queryId = strsplit(queryName, '.');
    queryId = queryId{1};
    queryId = uint32(str2num(queryId));
    posesRow = posesTable(posesTable.id==queryId, :);
    
    referenceSpace = posesRow.space{1,1};

    referenceT = [posesRow.x; posesRow.y; posesRow.z];
    referenceOrientation = [posesRow.dirx; posesRow.diry; posesRow.dirz];
    
    errors(i).queryId = queryId;
    if strcmp(spaceName, referenceSpace)
        errors(i).translation = norm(T - referenceT);
    else
        errors(i).translation = 666;
    end
    errors(i).orientation = atan2d(norm(cross(orientation,referenceOrientation)),dot(orientation,referenceOrientation));
    errors(i).inMap = posesRow.inMap;
    
    retrievedPoses(i).id = queryId;
    retrievedPoses(i).x = T(1);
    retrievedPoses(i).y = T(2);
    retrievedPoses(i).z = T(3);
    retrievedPoses(i).dirx = orientation(1);
    retrievedPoses(i).diry = orientation(2);
    retrievedPoses(i).dirz = orientation(3);
    retrievedPoses(i).space = spaceName;
end

% errors
errorsTable = struct2table(errors);
errors = table2struct(sortrows(errorsTable, 'queryId'));
errorsFile = fopen(params.evaluation.errors.path, 'w');
fprintf(errorsFile, 'id,translation,orientation\n');
for i=1:nQueries
    fprintf(errorsFile, '%d,%0.2f,%0.2f\n', errors(i).queryId, errors(i).translation, errors(i).orientation);
end
fclose(errorsFile);

% retrievedPoses
retrievedPosesTable = struct2table(retrievedPoses);
retrievedPoses = table2struct(sortrows(retrievedPosesTable, 'id'));
retrievedPosesFile = fopen(params.evaluation.retrieved.poses.path, 'w');
fprintf(retrievedPosesFile, 'id x y z dirx diry dirz space\n');
for i=1:nQueries
    fprintf(retrievedPosesFile, '%d %g %g %g %g %g %g %s\n', retrievedPoses(i).id, ...
        retrievedPoses(i).x, retrievedPoses(i).y, retrievedPoses(i).z, ...
        retrievedPoses(i).dirx, retrievedPoses(i).diry, retrievedPoses(i).dirz, ...
        retrievedPoses(i).space);
end
fclose(retrievedPosesFile);

%% summary
summaryFile = fopen(params.evaluation.summary.path, 'w');
thresholds = [[0.25 10], [0.5 10], [1 10]];
scores = zeros(1, size(thresholds,2)/2);
inMapScores = scores;
offMapScores = scores;
fprintf(summaryFile, 'Conditions: ');
for i=1:2:size(thresholds,2)
    if i > 1
        fprintf(summaryFile, ' / ');
    end
    fprintf(summaryFile, '(%g [m], %g [deg])', thresholds(i), thresholds(i+1));
    
    count = 0;
    inMapCount = 0;
    offMapCount = 0;
    inMapSize = 0;
    offMapSize = 0;
    for j=1:size(errors,1)
        if errors(j).translation < thresholds(i) && errors(j).orientation < thresholds(i+1)
            count = count + 1;
            if errors(j).inMap
                inMapCount = inMapCount + 1;
            else
                offMapCount = offMapCount + 1;
            end
        end
        if errors(j).inMap
            inMapSize = inMapSize + 1;
        else
            offMapSize = offMapSize + 1;
        end
    end

    scores((i-1)/2+1) = count / size(errors,1) * 100;
    inMapScores((i-1)/2+1) = inMapCount / inMapSize * 100;
    offMapScores((i-1)/2+1) = offMapCount / offMapSize * 100;
end
fprintf(summaryFile, '\n');
for i=1:size(scores,2)
    if i > 1
        fprintf(summaryFile, ' / ');
    end
    fprintf(summaryFile, '%g [%%]', scores(i));
end
fprintf(summaryFile, '\n');

% inMap
for i=1:size(inMapScores,2)
    if i > 1
        fprintf(summaryFile, ' / ');
    end
    fprintf(summaryFile, '%0.2f [%%]', inMapScores(i));
end
fprintf(summaryFile, ' -- InMap\n');

% offMap
for i=1:size(offMapScores,2)
    if i > 1
        fprintf(summaryFile, ' / ');
    end
    fprintf(summaryFile, '%0.2f [%%]', offMapScores(i));
end
fprintf(summaryFile, ' -- OffMap\n');

fclose(summaryFile);
%% initialization
[ params ] = setup_project_ht_WUSTL;

load(params.input.qlist_matname, 'query_imgnames_all');
densePV_matname = fullfile(params.output.dir, 'densePV_top10_shortlist.mat');
load(densePV_matname, 'ImgList');

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
for i=1:nQueries
    queryName = query_imgnames_all{i};
    queryImage = imread(fullfile(params.data.dir, params.data.q.dir, queryName));
    
    fun = @(x) strcmp(ImgList(x).queryname,queryName);
    tf = arrayfun(fun, 1:numel(ImgList));
    ImgListRecord = ImgList(find(tf));
    
    cutoutPath = ImgListRecord.topNname{1};
    cutoutPath = strsplit(cutoutPath, '/');
    spaceName = cutoutPath{1};
    sweepId = cutoutPath{2};
    transPath = fullfile(params.data.dir, params.data.db.trans.dir, spaceName, 'transformations', sprintf('trans_%s.txt', sweepId));
    P1 = load_CIIRC_transformation(transPath);
    R1 = P1(1:3,1:3);
    T1 = P1(1:3,4);
    
    P2 = ImgListRecord.P{1};
    R2 = P2(1:3,1:3); % in fact this is K*R - are you sure? viz demo
    T2 = P2(1:3,4);
    
    dslevel = 8^-1;
    Iq = queryImage;
    fl = 3172 * dslevel;
    K = [fl, 0, size(Iq, 2)/2.0; 0, fl, size(Iq, 1)/2.0; 0, 0, 1];
    
    P = P2*P1;
    T = -inv(P(:,1:3))*P(:,4);
    initialDirection = [0.0; 0.0; 1.0]; % TODO: verify
    orientation = R2 * (R1 * initialDirection);
    
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
end

errorsTable = struct2table(errors);
errors = table2struct(sortrows(errorsTable, 'queryId'));

errorsFile = fopen(params.evaluation.errors.path, 'w');
fprintf(errorsFile, 'id,translation,orientation\n');
for i=1:nQueries
    fprintf(errorsFile, '%d,%0.2f,%0.2f\n', errors(i).queryId, errors(i).translation, errors(i).orientation);
end
fclose(errorsFile);

%% summary
summaryFile = fopen(params.evaluation.summary.path, 'w');
thresholds = [[0.25 10], [0.5 10], [1 10]];
scores = zeros(1, size(thresholds,2)/2);
fprintf(summaryFile, 'Conditions: ');
for i=1:2:size(thresholds,2)
    if i > 1
        fprintf(summaryFile, ' / ');
    end
    fprintf(summaryFile, '(%g [m], %g [deg])', thresholds(i), thresholds(i+1));
    
    count = 0;
    for j=1:size(errors,1)
        if errors(j).translation < thresholds(i) && errors(j).orientation < thresholds(i+1)
            count = count + 1;
        end
    end

    scores((i-1)/2+1) = count / size(errors,1) * 100;
end
fprintf(summaryFile, '\n');
for i=1:size(scores,2)
    if i > 1
        fprintf(summaryFile, ' / ');
    end
    fprintf(summaryFile, '%g [%%]', scores(i));
end
fclose(summaryFile);
fprintf('\n');
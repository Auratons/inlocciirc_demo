% user configuration
addpath('utils/');
addpath('tools/');
addpath('visual_inspection/');
%% initialization

[filepath, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'functions', 'yaml'));
addpath(fullfile(filepath, '..', 'functions', 'utils'));
addpath(fullfile(filepath, '..', 'functions', 'inLocCIIRC_utils', 'at_netvlad_function'));
addpath(fullfile(filepath, '..', 'functions', 'at_netvlad_function'));
addpath(fullfile(filepath, '..', 'functions', 'relja_matlab'));
addpath(fullfile(filepath, '..', 'functions', 'relja_matlab', 'matconvnet'));
addpath(fullfile(filepath, '..', 'functions', 'netvlad'));
addpath(fullfile(filepath, '..', 'functions', 'wustl_function'));
addpath(fullfile(filepath, '..', 'functions', 'yael_matlab_linux64_v438'));
addpath(fullfile(filepath, '..', 'functions', 'at_netvlad_function'));
addpath(fullfile(filepath, '..', 'functions', 'ht_pnp_function'));
run(fullfile(filepath, '..', 'functions', 'vlfeat', 'toolbox', 'vl_setup.m'));
run(fullfile(filepath, '..', 'functions', 'matconvnet', 'matlab', 'vl_setupnn.m'));

% To support debug in Matlab -desktop and running without an X server on the cluster.
% See dvc/scripts/inloc_pose_verification.sh where these variables are generated dynamically.
if ~exist('params_file', 'var')
    params_file = '/home/kremeto1/inloc/dvc/pipeline-artwin-conv5-pyrender/params.yaml';
    experiment_name = 'main';
end

parameters = ReadYaml(params_file).(experiment_name);
eval_params = parameters.evaluation;

if isfield(eval_params, 'input_query_mat_path')
    qlist = eval_params.input_query_mat_path;
else
    qlist = fullfile(eval_params.root_to_process, 'query_imgnames_all.mat');
end

params = struct();
params.output.dir = eval_params.root_to_process;
params.input.qlist.path = qlist;
params.evaluation.dir = fullfile(params.output.dir, 'evaluations');
params.evaluation.retrieved.poses.dir = fullfile(params.evaluation.dir, 'retrieved_poses');
params.evaluation.query_vs_synth.dir = fullfile(params.evaluation.dir, 'query_vs_synth');
params.evaluation.query_segments_vs_synth_segments.dir = fullfile(params.evaluation.dir, 'query_segments_vs_synth_segments');
params.input_candidate_pose_renders_path = fullfile(params.output.dir, 'candidate_renders');
params.evaluation.errors.path = fullfile(params.evaluation.dir, 'localization_errors.txt');
params.evaluation.retrieved.queries.path = fullfile(params.evaluation.dir, 'retrieved_queries.txt');
params.evaluation.summary.path = fullfile(params.evaluation.dir, 'summary.txt');
params.evaluation.overlays = fullfile(params.evaluation.dir,'all_results');

load(params.input.qlist.path);

value = load(params.input.qlist.path, 'filenames');
query_imgnames_all = value.filenames;

densePV_matname = fullfile(params.output.dir, 'densePV_top10_shortlist.mat');
load(densePV_matname, 'ImgList_rendered');
ImgList = ImgList_rendered;

mkdirIfNonExistent(params.evaluation.dir);
mkdirIfNonExistent(params.evaluation.retrieved.poses.dir);
mkdirIfNonExistent(params.evaluation.query_vs_synth.dir);
mkdirIfNonExistent(params.evaluation.query_segments_vs_synth_segments.dir);
mkdirIfNonExistent(params.evaluation.overlays);

%% quantitative results
nQueries = size(query_imgnames_all, 2);
query_eval = cell(1, numel(query_imgnames_all));
errors = struct();
retrievedQueries = struct();
inLocCIIRCLostCount = 0;
lostIds = [];
for i=1:numel(ImgList)
    fprintf('Processing %d/%d\n', i, numel(ImgList));

    render_pose_filename = sprintf('%s_params.json', erase(erase(erase(ImgList(i).render_path, "_color.png"), "_out.png"), "_neural"));
    queryPoseFilename = sprintf('%s_params.json', erase(ImgList(i).query_path, "_reference.png"));
    % queryPoseFilename = strrep(queryPoseFilename, 'joined-dataset-spheres', 'joined-dataset-pyrender-black_bg');

    % [P,ref_spaceName,fullName] = getReferencePose(i,ImgList,params);
    fid = fopen(queryPoseFilename);
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    val = jsondecode(str);
    P = val.camera_pose;

    if ~all(isnan(ImgList(i).P), 'all')
        fid = fopen(render_pose_filename);
        raw = fread(fid, inf);
        str = char(raw');
        fclose(fid);
        val = jsondecode(str);
        P_est.P = val.camera_pose;
        P_est.P(1:3, 2:3) = - P_est.P(1:3, 2:3);
    else
        P_est.P = ImgList(i).P;
    end

    % P is camera pose (extrinsic matrix) from _param.json files
    P_ref = {};
    P_ref.P = P;
    P_ref.R = P_ref.P(1:3, 1:3)';
    P_ref.t = P_ref.P(1:3, 4);
    P_ref.C = -P_ref.R * P_ref.t;

    P_est.R = P_est.P(1:3, 1:3)';
    P_est.t = P_est.P(1:3, 4);
    P_est.C = -P_est.R * P_est.t;

    % P is camera rotation and position
    % P_est = {};
    % P_est.P = ImgList(i).P;
    % [P_est.K, P_est.R, P_est.C] = P2KRC(P_est.P);
    % P_est.R = P_est.P(1:3, 1:3);
    % P_est.C = P_est.P(1:3, 4);
    % P_est.t = -P_est.R' * P_est.C;
    % est_spaceName = strsplit(ImgList(i).topNname{1},'/'); est_spaceName = est_spaceName{1};
    % est_mapName = strsplit(est_spaceName,'_'); est_mapName = est_mapName{1};
    % ref_mapName = strsplit(ref_spaceName,'_'); ref_mapName = ref_mapName{1};

    transform = eye(4);
    C_ref = P_ref.C;
    R_ref = P_ref.R;

    % query_eval{i}.pano_id = strsplit(fullName,'/');
    % query_eval{i}.pano_id = query_eval{i}.pano_id{1};
    query_eval{i}.id = i;
    query_eval{i}.C_ref = C_ref;
    query_eval{i}.R_ref = R_ref;
    query_eval{i}.C_ref_orig = P_ref.C;
    query_eval{i}.R_ref_orig = P_ref.R;
    query_eval{i}.C_est= P_est.C;
    query_eval{i}.R_est = P_est.R;
    query_eval{i}.spaceTransform = transform;

    estimate_params = sprintf('%s_params.json', erase(ImgList(i).render_path, "_color.png"));
    if exist(estimate_params, 'file') == 2
        fid = fopen(estimate_params);
        raw = fread(fid, inf);
        str = char(raw');
        fclose(fid);
        val = jsondecode(str);
        query_eval{i}.est_space = val.localized_scan;
        query_eval{i}.ref_space = val.source_scan;
    else
        query_eval{i}.est_space = "";
        query_eval{i}.ref_space = "";
    end

    errors(i).translation = norm(C_ref - P_est.C);
    errors(i).orientation = rotationDistance(R_ref, P_est.R);
    errors(i).queryId =  ImgList(i).query_path;
    errors(i).inMap = strcmp(query_eval{i}.ref_space, query_eval{i}.est_space);
    if isnan(errors(i).translation)
        lostIds = [lostIds i];
        inLocCIIRCLostCount = inLocCIIRCLostCount + 1;
    end

    visual_inspection = false;
    if errors(i).translation > 2 && visual_inspection
        %         getSynthView(params,ImgList,i,1,true,params.evaluation.dir,sprintf('err_%.2fm',errors(i).translation));
    end

    if visual_inspection  && ~isnan(errors(i).translation)
        prefix = sprintf('err_%.2fm_%.0fdeg', errors(i).translation, errors(i).orientation);
        fname = sprintf('%s_results_q_id_%d_best_db_%d.jpg', prefix, i, 1);
        if ~exist(fullfile(params.evaluation.overlays, fname))
            getSynthView(params, ImgList, i, 1, params.evaluation.overlays, prefix);
        end
    end
    % if visual_inspection && isnan(errors(i).translation)
    %     if ~exist(fullfile(params.evaluation.overlays,sprintf('%s_results_q_id_%d_best_db_%d.jpg',sprintf('err_NaN'),i,1)))
    %         getFailureView(params,ImgList,i,1,true,params.evaluation.overlays,sprintf('%s_results_q_id_%d_best_db_%d.jpg',sprintf('err_NaN'),i,1));
    %     end
    % end

    [~, name, ext] = fileparts(queryPoseFilename);
    retrievedPosePath = fullfile(params.evaluation.retrieved.poses.dir, strcat(name,ext));
    retrievedPoseFile = fopen(retrievedPosePath, 'w');
    fprintf(retrievedPoseFile, jsonencode(P_est.P));
    fclose(retrievedPoseFile);

    retrievedQueries(i).id = i;
    retrievedQueries(i).space = query_eval{i}.est_space;

    sprintf('DONE %s', ImgList(i).query_path);
end

save(fullfile(params.evaluation.dir,'query_eval.mat'),'query_eval', '-v7');

% errors
errors = struct2table(errors);
errors = table2struct(sortrows(errors, 'queryId'));
errorsFile = fopen(params.evaluation.errors.path, 'w');
fprintf(errorsFile, 'id,inMap,translation,orientation\n');
for i=1:nQueries
    inMapStr = 'No';
    if errors(i).inMap
        inMapStr = 'Yes';
    end
    fprintf(errorsFile, '%s,%s,%0.4f,%0.4f\n', errors(i).queryId, inMapStr, errors(i).translation, errors(i).orientation);
end
fclose(errorsFile);

meaningfulTranslationErrors = [errors(~isnan([errors.translation])).translation];
meaningfulOrientationErrors = [errors(~isnan([errors.orientation])).orientation];

% statistics of the errors
meanTranslation = mean(meaningfulTranslationErrors);
meanOrientation = mean(meaningfulOrientationErrors);
medianTranslation = median(meaningfulTranslationErrors);
medianOrientation = median(meaningfulOrientationErrors);
stdTranslation = std(meaningfulTranslationErrors);
stdOrientation = std(meaningfulOrientationErrors);

% retrievedQueries
retrievedQueries = struct2table(retrievedQueries);
retrievedQueries = table2struct(sortrows(retrievedQueries, 'id'));
retrievedQueriesFile = fopen(params.evaluation.retrieved.queries.path, 'w');
fprintf(retrievedQueriesFile, 'id space\n');
for i=1:nQueries
    fprintf(retrievedQueriesFile, '%d %s\n', retrievedQueries(i).id, retrievedQueries(i).space);
end
fclose(retrievedQueriesFile);

%%%%%%%%%%%
% SUMMARY %
%%%%%%%%%%%

summaryFile = fopen(params.evaluation.summary.path, 'w');
thresholds = [[0.25; 2], [0.5; 5], [2.5; 7.5] [5; 10], [7.5; 15], [10; 20], [15; 30], [20; 30]];
% thresholds =  logspace(0,2,32)/20;
% thresholds = [thresholds; 10*ones(1,size(thresholds,2))];

inMapSize = 0;
offMapSize = 0;
for j=1:length(errors)
    if errors(j).inMap
        inMapSize = inMapSize + 1;
    else
        offMapSize = offMapSize + 1;
    end
end

if inMapSize == 0
    inMapScores = zeros(1, size(thresholds, 2));
else
    inMapScores = 100.0 * ones(1, size(thresholds, 2)) / double(inMapSize);
end
if offMapSize == 0
    offMapScores = zeros(1, size(thresholds, 2));
else
    offMapScores = 100.0 * ones(1, size(thresholds, 2)) / double(offMapSize);
end
scores = 100.0 * ones(1, size(thresholds,2)) / double(length(errors));

for i=1:size(thresholds,2)
    count = 0.0;
    inMapCount = 0.0;
    offMapCount = 0.0;
    for j=1:length(errors)
        if errors(j).translation < thresholds(1,i) && errors(j).orientation < thresholds(2,i)
            if errors(j).inMap
                count = count + 1;
                inMapCount = inMapCount + 1;
            else
                offMapCount = offMapCount + 1;
            end
        end
    end

    % we want to include cases InLoc got lost, but not blacklisted queries (=no reference poses)
    scores(i) = count * scores(i);
    inMapScores(i) = inMapCount * inMapScores(i);
    offMapScores(i) = offMapCount * offMapScores(i);
end

fprintf(summaryFile, 'Thresholds            Percentage   InMap      OffMap\n');
for i=1:size(scores,2)
    fprintf(summaryFile, '(%5.2f m, %5.2f deg):   ', thresholds(1,i), thresholds(2,i));
    fprintf(summaryFile, '%5.1f %%  ', scores(i));
    fprintf(summaryFile, '%4.1f [%%]    ', inMapScores(i));
    fprintf(summaryFile, '%4.1f [%%]\n', offMapScores(i));
end
fprintf(summaryFile, '\n');
fprintf(summaryFile, '\nInLocCIIRC got completely lost %d out of %d times. Not included in the mean/median/std errors.\n', ...
    inLocCIIRCLostCount, nQueries);
fprintf(summaryFile, '\nInLocCIIRC selected a wrong map %d out of %d times.\n', ...
    offMapSize, nQueries);
fprintf(summaryFile, '\nErrors (InLocCIIRC poses wrt reference poses):\n');
fprintf(summaryFile, ' \ttranslation [m]\torientation [deg]\n');
fprintf(summaryFile, 'Mean\t%0.2f\t%0.2f\n', meanTranslation, meanOrientation);
fprintf(summaryFile, 'Median\t%0.2f\t%0.2f\n', medianTranslation, medianOrientation);
fprintf(summaryFile, 'Std\t%0.2f\t%0.2f\n', stdTranslation, stdOrientation);
fclose(summaryFile);
disp(fileread(params.evaluation.summary.path));

f = figure('visible', 'off');
grid on;
yyaxis left
plot(scores, thresholds(1,:), 'Marker', '.', 'MarkerSize', 20);
ylabel('Distance threshold [m]');
% hold on;
% plot3d([thresholds(2,:); scores], '-b', 'Marker', '.', 'MarkerSize', 20);

yyaxis right
plot(scores, thresholds(2,:), 'Marker', '.', 'MarkerSize', 20);
ylabel('Angular distance threshold [deg]');

xlabel('Correctly localised queries [%]');
legend('Distance', 'Angular distance');

path_parts = strsplit(params_file, filesep);
dset_name = path_parts(length(path_parts) - 1);
title(sprintf('InLoc at %s', dset_name{1}));
saveas(gcf, fullfile(params.evaluation.dir, 'correctly_localized_queries.jpg'));
close(f);


function distance = rotationDistance(R1, R2)
    R = R2*R1'; % residual rotation matrix
    distance = acos(0.5 * (trace(R)-1));
    distance = abs(rad2deg(distance));
end

function mkdirIfNonExistent(pathToDirectory)
    if exist(pathToDirectory, 'dir') ~= 7
        mkdir(pathToDirectory);
    end
end

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
    params_file = '/home/kremeto1/inloc/dvc/pipeline-artwin-conv5-spheres/params.yaml';
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
%params.output.dir = '/home/kremeto1/inloc/datasets/pipeline-artwin-conv5-spheres/';
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
params.output.gv_dense.dir = fullfile(params.output.dir, 'gv_dense');

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

%% quantitative results
nQueries = size(query_imgnames_all, 2);
query_eval = cell(1, numel(query_imgnames_all));
errors = struct();
retrievedQueries = struct();
inLocCIIRCLostCount = 0;
lostIds = [];
for i=1:numel(ImgList)
    fprintf('Processing %d/%d\n', i, numel(ImgList));

    queryPoseFilename = sprintf('%s_params.json', erase(ImgList(i).query_path, "_reference.png"));
    queryPoseFilename = strrep(queryPoseFilename, 'joined-dataset-spheres', 'joined-dataset-pyrender-black_bg');
    mkdirIfNonExistent(params.evaluation.dir);


    % [P,ref_spaceName,fullName] = getReferencePose(i,ImgList,params);
    fid = fopen(queryPoseFilename);
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    val = jsondecode(str);
    P = val.camera_pose;

    P_ref = {};
    P_ref.R = P(1:3, 1:3);
    P_ref.t = P(1:3, 4);
    P_ref.P = P;
    P_ref.C = -P_ref.R' * P_ref.t;

    P_est = {};
    P_est.P = ImgList(i).P;
    [P_est.K, P_est.R, P_est.C] = P2KRC(P_est.P);
    P_est.t = -P_est.R * P_est.C;
    % est_spaceName = strsplit(ImgList(i).topNname{1},'/'); est_spaceName = est_spaceName{1};
    % est_mapName = strsplit(est_spaceName,'_'); est_mapName = est_mapName{1};
    % ref_mapName = strsplit(ref_spaceName,'_'); ref_mapName = ref_mapName{1};

    C_ref = [];
    R_ref = [];
    if 0 %~strcmp(est_spaceName,ref_spaceName) &&  strcmp(est_mapName,ref_mapName)
        interesting = true;
        transform = [];
        E_h_12 = [1.000000000000 0.000265824958 -0.000320481340 -0.965982019901;
            -0.000265917915 0.999999940395 -0.000290183641 0.005340866279;
            0.000320404157 0.000290269381 1.000000119209 0.241866841912;
            0.000000000000 0.000000000000 0.000000000000 1.000000000000];

        E_l_21 = [ 0.999996125698 0.000008564073 0.002756817034 3.283028602600;
            -0.000006930272 0.999999821186 -0.000592759345 0.000593465462;
            -0.002756824018 0.000592722441 0.999996006489 1.970497488976;
            0.000000000000 0.000000000000 0.000000000000 1.000000000000];
        switch (est_spaceName)
            case params.dataset.db.space_names{1}
                transform = E_h_12;
            case params.dataset.db.space_names{2}
                transform = inv(E_h_12);
            case params.dataset.db.space_names{3}
                transform = inv(E_l_21);
            case params.dataset.db.space_names{4}
                transform = E_l_21;
            otherwise
                interesting = true
                error('ups')
                % getSynthView(params,ImgList,i,1,true);
        end

        C_ref = transform*[P_ref.C; 1];
        C_ref = C_ref(1:3);
        R_ref = P_ref.R*transform(1:3, 1:3);
    else
        transform = eye(4);
        C_ref = P_ref.C;
        R_ref = P_ref.R;
    end
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
    inLocCIIRCLostCount = inLocCIIRCLostCount + isnan(errors(i).translation);
    if isnan(errors(i).translation)
        lostIds = [lostIds i];
    end

    visual_inspection = true;
    if errors(i).translation > 2 && visual_inspection
        %         getSynthView(params,ImgList,i,1,true,params.evaluation.dir,sprintf('err_%.2fm',errors(i).translation));

    end

    mkdirIfNonExistent(fullfile(params.evaluation.dir,'all_results'));
    if visual_inspection  && ~isnan(errors(i).translation)
        if ~exist(fullfile(fullfile(params.evaluation.dir,'all_results'),sprintf('%s_results_q_id_%d_best_db_%d.jpg',sprintf('err_%.2fm_%.0fdeg',errors(i).translation,errors(i).orientation),i,1)))
            getSynthView(params,ImgList,i,1,true,fullfile(params.evaluation.dir,'all_results'),sprintf('err_%.2fm_%.0fdeg',errors(i).translation,errors(i).orientation));
        end
    end
    if visual_inspection && isnan(errors(i).translation)
        if ~exist(fullfile(fullfile(params.evaluation.dir,'all_results'),sprintf('%s_results_q_id_%d_best_db_%d.jpg',sprintf('err_NaN'),i,1)))
            getFailureView(params,ImgList,i,1,true,fullfile(params.evaluation.dir,'all_results'),sprintf('%s_results_q_id_%d_best_db_%d.jpg',sprintf('err_NaN'),i,1));
        end
    end



    [~, name] = fileparts(queryPoseFilename);
    retrievedPosePath = fullfile(params.evaluation.retrieved.poses.dir, name);
    mkdirIfNonExistent(fileparts(retrievedPosePath));
    retrievedPoseFile = fopen(retrievedPosePath, 'w');
    % P_str = P_to_str(P_est.P);
    fprintf(retrievedPoseFile, jsonencode(P_est.P));
    fclose(retrievedPoseFile);

    retrievedQueries(i).id = i;
    % retrievedQueries(i).space = ref_spaceName;

    sprintf('DONE %s', ImgList(i).query_path);
end
save(fullfile(params.evaluation.dir,'query_eval.mat'),'query_eval', '-v7');
% errors
errorsBak = errors;
errorsTable = struct2table(errors);
errors = table2struct(sortrows(errorsTable, 'queryId'));
errorsFile = fopen(params.evaluation.errors.path, 'w');
fprintf(errorsFile, 'id,inMap,translation,orientation\n');
for i=1:nQueries
    % inMapStr = 'No';
    % if errors(i).inMap
    %     inMapStr = 'Yes';
    % end
    fprintf(errorsFile, '%d,%s,%0.4f,%0.4f\n', errors(i).queryId, "Yes", errors(i).translation, errors(i).orientation);
end
fclose(errorsFile);
errors = errorsBak; % we cannot use the sorted. it would break compatibility with blacklistedQueries array!

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
retrievedQueriesTable = struct2table(retrievedQueries);
retrievedQueries = table2struct(sortrows(retrievedQueriesTable, 'id'));
retrievedQueriesFile = fopen(params.evaluation.retrieved.queries.path, 'w');
fprintf(retrievedQueriesFile, 'id space\n');
for i=1:nQueries
    fprintf(retrievedQueriesFile, '%d %s\n', retrievedQueries(i).id, ...
        "");
        %retrievedQueries(i).space);
end
fclose(retrievedQueriesFile);

%% summary
summaryFile = fopen(params.evaluation.summary.path, 'w');
% thresholds = [[0.05 10],[0.10 10],[0.15 10],[0.20 10],[0.25 10], [0.5 10], [0.75 10],[1 10],[2 10],[5 10]];
thresholds =  logspace(0,2,32)/50;
thresholds = [thresholds; 10*ones(1,size(thresholds,2))];
scores = zeros(1, size(thresholds,2));
inMapScores = scores;
offMapScores = scores;
fprintf(summaryFile, 'Conditions: ');
for i=1:size(thresholds,2)
    if i > 1
        fprintf(summaryFile, ' / ');
    end
    fprintf(summaryFile, '(%g [m], %g [deg])', thresholds(1,i), thresholds(2,i));

    count = 0;
    inMapCount = 0;
    offMapCount = 0;
    inMapSize = 0;
    offMapSize = 0;
    for j=1:length(errors)
        if errors(j).translation < thresholds(1,i) && errors(j).orientation < thresholds(2,i)
            count = count + 1;
            if errors(j).inMap
                inMapCount = inMapCount + 1;
            else
                offMapCount = offMapCount + 1;
            end
        end
%         if errors(j).inMap
%             inMapSize = inMapSize + 1;
%         else
%             offMapSize = offMapSize + 1;
%         end
    end

    % we want to include cases InLoc got lost, but not blacklisted queries (=no reference poses)
    nMeaningfulErrors = length(errors);
    scores(i) = count / nMeaningfulErrors * 100;
    inMapScores(i) = inMapCount / inMapSize * 100;
    offMapScores(i) = offMapCount / offMapSize * 100;
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
% for i=1:size(offMapScores,2)
%     if i > 1
%         fprintf(summaryFile, ' / ');
%     end
%     fprintf(summaryFile, '%0.2f [%%]', offMapScores(i));
% end
% fprintf(summaryFile, ' -- OffMap\n');
fprintf(summaryFile, '\nInLocCIIRC got completely lost %d out of %d times. Not included in the mean/median/std errors.\n', ...
    inLocCIIRCLostCount, nQueries);
fprintf(summaryFile, '\nInLocCIIRC selected a wrong map %d out of %d times.\n', ...
    offMapCount, nQueries);
fprintf(summaryFile, '\nErrors (InLocCIIRC poses wrt reference poses):\n');
fprintf(summaryFile, ' \ttranslation [m]\torientation [deg]\n');
fprintf(summaryFile, 'Mean\t%0.2f\t%0.2f\n', meanTranslation, meanOrientation);
fprintf(summaryFile, 'Median\t%0.2f\t%0.2f\n', medianTranslation, medianOrientation);
fprintf(summaryFile, 'Std\t%0.2f\t%0.2f\n', stdTranslation, stdOrientation);
fclose(summaryFile);
disp(fileread(params.evaluation.summary.path));

f = figure('visible', 'off');
thr_t = thresholds(1,:);
thr_t;
% plot3d([(0:size(scores,2))/2; 0 scores],'-b');
plot3d([thr_t;scores],'-b','Marker','.','MarkerSize',20); hold on;
grid on;
% ax = gca
xticks(gca,[0.1,0.15,0.2,(1:8)/4])
xtickangle(-75)
% xticklabels(gca,strsplit(num2str(thr_t)))
hold on;
xlabel('Distance threshold [m]');
ylabel('Correctly localised queries [%]');
title('InLoc SPRING at Broca');
saveas(gcf,fullfile(params.evaluation.dir,'correctly_localized_queries.jpg'));
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

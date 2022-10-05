[filepath, ~, ~] = fileparts(mfilename('fullpath'));
addpath(fullfile(filepath, '..', 'inLocCIIRC_dataset', 'inloc_build_file_lists'));
addpath(fullfile(filepath, '..', 'inLocCIIRC_dataset', 'inloc_compute_features'));
addpath(fullfile(filepath, '..', 'inLocCIIRC_dataset', 'inloc_compute_scores'));
addpath(fullfile(filepath));

% inloc_hw = getenv("INLOC_HW");
% if isempty(inloc_hw) || (~strcmp(inloc_hw, "GPU") && ~strcmp(inloc_hw, "CPU"))
%     fprintf('Please specify environment variable INLOC_HW to one of: "GPU", "CPU"\n');
%     fprintf('CPU mode will run on many cores (unsuitable for boruvka).\n');
%     fprintf('GPU mode will run on maximum of 4 cores, but with a GPU.\n');
%     fprintf('NOTE: You should first run InLocCIIRC on GPU, then run it on CPU.\n')
%     error("See above");
% end
% fprintf('InLocCIIRC is running in %s mode.\n', inloc_hw);
%
% delete(gcp('nocreate'));
% if strcmp(inloc_hw, "CPU")
%     if strcmp(environment(), 'laptop')
%         nWorkers = 8;
%     else
%         nWorkers = 45;
%     end
%     c = parcluster;
%     c.NumWorkers = nWorkers;
%     saveProfile(c);
%     p = parpool('local', nWorkers);
% end

inloc_build_file_lists(params_file, 'experiment_name', experiment_name);

inloc_compute_features(params_file, 'experiment_name', experiment_name);

inloc_compute_scores(params_file, 'experiment_name', experiment_name);

inloc_retrieval(params_file, 'experiment_name', experiment_name);

inloc_dense_pose_estimation(params_file, 'experiment_name', experiment_name);

inloc_neural_pose_verification(params_file, 'experiment_name', experiment_name);

% evaluate;

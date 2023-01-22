% Expects:
% params_file path to yaml
% experiment_name with name of dict in the yaml file

% To support debug in Matlab -desktop and running without an X server on the cluster.
% See dvc/scripts/inloc_pose_verification.sh where these variables are generated dynamically.
if ~exist('params_file', 'var')
    params_file = '/home/kremeto1/inloc/dvc/pipeline-grand-conv5-pyrender/params.yaml';
    experiment_name = 'main';
end

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


parameters = ReadYaml(params_file).(experiment_name);

workers_cnt = 16;

types = {'db', 'query'};

fprintf('\n%s >>> Running inloc_build_file_lists...\n', datestr(now,'HH:MM:SS.FFF'));
params = parameters.file_lists;

file_lists = struct();
for idx = 1:length(types)
    type = types{idx};
    output_mat_path = params.(sprintf('output_%s_mat_path', type));

    if exist(output_mat_path, 'file') ~= 2
        glob = get_with_default(params, sprintf('input_%s_glob', type), '*.jpg');
        files = dir(fullfile(params.(sprintf('input_%s_dir', type)), glob));
        nFiles = size(files, 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % nFiles = fix(nFiles / 10);%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % if strcmp(type, 'db')
        %     nFiles = 20;
        % else
        %     nFiles = 5;
        % end
        filenames = cell(1, nFiles);

        for i=1:nFiles
            filenames{1, i} = fullfile(files(i).folder, files(i).name);
        end

        create_parent_folder(output_mat_path);
        save(output_mat_path, 'filenames');

        file_lists.(sprintf('%s_filenames', type)) = filenames;
    else
        fprintf('SKIPPING BUILDING FILE LISTS, output "%s" already exists.\n', output_mat_path);
    end
end


fprintf('\n%s >>> Running inloc_compute_features...\n', datestr(now,'HH:MM:SS.FFF'));
params = parameters.features;

load(params.input_netvlad_pretrained, 'net');
net = relja_simplenn_tidy(net);
net = relja_cropToLayer(net, 'postL2');

for idx = 1:length(types)
    type = types{idx};
    output_path = params.(sprintf('output_%s_features_path', type));

    if exist(output_path, 'dir') ~= 7
        if ~isfield(file_lists, sprintf('%s_filenames', type))
            value = load(params.(sprintf('input_%s_mat_path', type)), 'filenames');
            file_lists.(sprintf('%s_filenames', type)) = value.filenames;
        end
        filenames = file_lists.(sprintf('%s_filenames', type));
        file_count = size(filenames, 2);
        features = struct();
        create_parent_folder(output_path);
        for i=1:file_count
            fprintf('%s Finding features for %s image #%d/%d\n\n', datestr(now,'HH:MM:SS.FFF'), type, i, file_count);
            img_path = filenames{i};
            img = imread(img_path);
            cnn = at_serialAllFeats_convfeat(net, img, 'useGPU', true);
            for l = setdiff([1, 2, 3, 4, 5, 6], [3, get_with_default(params, 'input_feature_layer', 6)])
                cnn{l} = [];
            end
            features.img_path = img_path;
            features.features = cnn;
            [~, name] = fileparts(img_path);
            save(fullfile(output_path, name + ".mat"), 'features', '-v7.3');
        end

    else
        fprintf('SKIPPING FEATURE COMPUTATION, output "%s" already exists.\n', output_path);
    end
end


fprintf('\n%s >>> Running inloc_compute_scores...\n', datestr(now,'HH:MM:SS.FFF'));
params = parameters.scores;

score = struct();
if exist(params.output_scores_mat_path, 'file') ~= 2
    imgs = dir(fullfile(params.input_db_features_path, '*.mat'));
    queries = dir(fullfile(params.input_query_features_path, '*.mat'));
    n_img = size(imgs, 1);
    n_query = size(queries, 1);

    coarse_feature_level = get_with_default(params, 'input_feature_layer', 6);
    score = struct('query_path', {}, 'scores', {}, 'db_score_paths', {});

    db_paths = cell(1, n_img);
    db_features = load(fullfile(imgs(1).folder, imgs(1).name));
    all_db_features = zeros(n_img, size(db_features.features.features{coarse_feature_level}.x(:), 1));
    for i=1:n_img
        db_features = load(fullfile(imgs(i).folder, imgs(i).name));
        all_db_features(i, :) = db_features.features.features{coarse_feature_level}.x(:)';
        db_paths{i} = db_features.features.img_path;
    end
    all_db_features = all_db_features';
    all_db_features = all_db_features ./ vecnorm(all_db_features);
    check_is_normalized(all_db_features);

    for i=1:n_query
        fprintf('%s Processing query %d/%d\n', datestr(now,'HH:MM:SS.FFF'), i, n_query);
        query_features = load(fullfile(queries(i).folder, queries(i).name));
        single_q_features = query_features.features.features{coarse_feature_level}.x(:)';
        single_q_features = single_q_features ./ vecnorm(single_q_features);
        check_is_normalized(single_q_features);
        single_q_features = repmat(single_q_features, n_img, 1)';
        similarityScores = dot(single_q_features, all_db_features);
        score(i).query_path = query_features.features.img_path;
        score(i).scores = single(similarityScores); % NOTE: this is not a probability distribution (and it does not have to be)
        score(i).db_score_paths = db_paths;
    end

    clearvars all_db_features

    create_parent_folder(params.output_scores_mat_path);
    save(params.output_scores_mat_path, 'score');
else
    fprintf('SKIPPING SCORE COMPUTATION, output "%s" already exists.\n', params.output_scores_mat_path);
end


fprintf('\n%s >>> Running inloc_retrieval...\n', datestr(now,'HH:MM:SS.FFF'));
params = parameters.retrieval;

ImgList = struct();
if exist(params.output_topN_mat_path, 'file') ~= 2
    shortlist_topN = get_with_default(params, 'input_topN', 100);
    ImgList = struct('query_path', {}, 'topN_db_paths', {}, 'topN_scores', {});

    if isempty(fieldnames(score))
        score = load(params.input_scores_mat_path, 'score').score;
    end

    for i = 1:size(score, 2)
        ImgList(i).query_path = score(i).query_path;
        [score_sort, score_idx] = sort(score(i).scores, 'descend');
        ImgList(i).topN_db_paths = score(i).db_score_paths(score_idx(1:shortlist_topN));
        ImgList(i).topN_scores = score_sort(1:shortlist_topN);
    end

    create_parent_folder(params.output_topN_mat_path);
    save('-v6', params.output_topN_mat_path, 'ImgList');
else
    fprintf('SKIPPING RETRIEVAL, output "%s" already exists.\n', params.output_topN_mat_path);
end


fprintf('\n%s >>> Running inloc_dense_pose_estimation...\n', datestr(now,'HH:MM:SS.FFF'));
params = parameters.pose_estimation;

pnp_topN = get_with_default(params, 'input_pnp_topN', 10);

% densePE (topN reranking -> pnp_topN pose candidate)
densePE_matname = params.output_candidate_mat_path;
if exist(densePE_matname, 'file') ~= 2

    fprintf('\n');
    delete(gcp('nocreate'));
    cluster = parcluster;
    cluster.NumWorkers = workers_cnt;
    saveProfile(cluster);
    parpool('local', workers_cnt);

    if isempty(fieldnames(ImgList))
        ImgList_retrieval = load(params.input_topN_mat_path, 'ImgList');
        ImgList_retrieval = ImgList_retrieval.ImgList;
    else
        ImgList_retrieval = ImgList;
    end

    ImgList = struct('query_path', {}, 'topN_db_paths', {}, 'topN_scores', {}, 'P', {}, 'K', {}, 'R', {}, 't', {});

    for ii = 1:length(ImgList_retrieval)
        curr_query_path = ImgList_retrieval(ii).query_path;
        curr_topN_db_paths = ImgList_retrieval(ii).topN_db_paths;

        ImgList(ii).query_path = curr_query_path;
        ImgList(ii).topN_db_paths = curr_topN_db_paths;

        % preload query feature
        [~, name] = fileparts(curr_query_path);
        query_features = load(fullfile(params.input_query_features_path, name + ".mat"));
        if curr_query_path == query_features.features.img_path
            cnnq = query_features.features.features;
        else
            error('Not the same!');
        end

        if exist(fullfile(params.output_gv_dense_dir, filename(curr_query_path)), 'dir') ~= 7
            mkdir(fullfile(params.output_gv_dense_dir, filename(curr_query_path)));
        end

        coarse_feature_level = get_with_default(params, 'input_feature_layer', 6);

        parfor (kk = 1:length(curr_topN_db_paths), workers_cnt)
            parfor_denseGV( cnnq, curr_query_path, curr_topN_db_paths{kk}, params );
        end

        for jj = 1:length(curr_topN_db_paths)
            this_gvresults = load(fullfile( ...
                params.output_gv_dense_dir, ...
                filename(curr_query_path), ...
                strcat(filename(curr_topN_db_paths{jj}), ".mat")));
            ImgList(ii).topN_scores(jj) = ImgList_retrieval(ii).topN_scores(jj) + size(this_gvresults.inls12, 2);
        end

        [sorted_score, idx] = sort(ImgList(ii).topN_scores, 'descend');
        ImgList(ii).topN_db_paths = ImgList(ii).topN_db_paths(idx);
        ImgList(ii).topN_scores = ImgList(ii).topN_scores(idx);

        fprintf('%s %s done.\n', datestr(now,'HH:MM:SS.FFF'), filename(curr_query_path));
    end

    % pnp list
    qlist = cell(1, length(ImgList) * pnp_topN);
    dblist = cell(1, length(ImgList) * pnp_topN);
    for ii = 1:length(ImgList)
        for jj = 1:pnp_topN
            qlist{pnp_topN * (ii - 1) + jj} = ImgList(ii).query_path;
            dblist{pnp_topN * (ii - 1) + jj} = ImgList(ii).topN_db_paths{jj};
        end
    end

    % dense pnp
    for ii = 1:length(ImgList)
        if exist(fullfile(params.output_pnp_dense_inlier_dir, filename(ImgList(ii).query_path)), 'dir') ~= 7
            mkdir(fullfile(params.output_pnp_dense_inlier_dir, filename(ImgList(ii).query_path)));
        end
        fprintf("Creating %s\n", fullfile(params.output_pnp_dense_inlier_dir, filename(ImgList(ii).query_path)));
    end
    parfor (ii = 1:length(qlist), workers_cnt)
        parfor_densePE( ...
            qlist{ii}, ...
            dblist{ii}, ...
            params ...
        );
    end

    % load top10 pnp
    for ii = 1:length(ImgList)
        ImgList(ii).P = cell(1, pnp_topN);
        for jj = 1:pnp_topN
            [~, dbbasename, ~] = fileparts(ImgList(ii).topN_db_paths{jj});
            this_densepe_matname = fullfile(params.output_pnp_dense_inlier_dir, filename(ImgList(ii).query_path), strcat(dbbasename, ".mat"));
            load(this_densepe_matname, 'P');
            ImgList(ii).P{jj} = P;
            [K, R, t] = P2KRC(P);
            ImgList(ii).K{jj} = K;
            ImgList(ii).R{jj} = R;
            ImgList(ii).t{jj} = t;
        end
    end

    create_parent_folder(densePE_matname);
    save('-v6', densePE_matname, 'ImgList');
    ImgList_densePE = ImgList;
else
    fprintf('SKIPPING POSE ESTIMATION, output "%s" already exists.\n', densePE_matname);
end


fprintf('\n%s >>> Running inloc_neural_pose_verification...\n', datestr(now,'HH:MM:SS.FFF'));
params = parameters.pose_verification;

PV_topN = get_with_default(params, 'input_pv_topN', 10);

neuralPV_matname = params.output_shortlist_mat_path;
candidate_renders = params.input_candidate_pose_renders_path;
if exist(candidate_renders, 'dir') == 7 && exist(neuralPV_matname, 'file') ~= 2
    %Load the poses from the .mat file
    if ~exist('ImgList_densePE', 'var')
        ImgList_densePE = load(params.input_candidate_mat_path, 'ImgList');
        ImgList_densePE = ImgList_densePE.ImgList;
    end
    ImgList_rendered = struct('query_path', {}, 'score', {}, 'P', {});

    for ii = 1:length(ImgList_densePE)
        qpath = ImgList_densePE(ii).query_path;
        [~, qstem, ext] = fileparts(filename(qpath));
        fprintf('\n%s Processing query %s%s\n', datestr(now,'HH:MM:SS.FFF'), qstem, ext) % debug
        % Load query image & output information
        Iq = imread(qpath);
        reference_size = size(Iq);
        reference_h = reference_size(1);
        reference_w = reference_size(2);

        this_P = nan(3, 4);
        this_score = 0;
        this_db = ImgList_densePE(ii).topN_db_paths(1);
        this_db = this_db{:};
        [~, dbstem, ~] = fileparts(this_db);
        % Just initial value for evaluation
        this_render = fullfile(candidate_renders, qstem, strcat(dbstem, "_color.png"));

        % Just desquerify if needed
        Is = imread(this_render);
        shape = size(Is);
        render_h = shape(1);
        render_w = shape(2);
        if (render_h ~= render_w && reference_h == reference_w) || (render_h == render_w && reference_h == reference_w && reference_h > render_h)
            square_size = reference_h;
            offset_h = fix((square_size - render_h) / 2);
            offset_w = fix((square_size - render_w) / 2);
            Iq = Iq(offset_h+1:offset_h + render_h, offset_w+1:offset_w + render_w, :);
        end

        % normalization
        Iq_norm = normalize(Iq);

        for jj = 1:PV_topN
            % Strip the rank of the rendering from the name of the image
            P = ImgList_densePE(ii).P(jj);
            P = P{:};
            db_path = ImgList_densePE(ii).topN_db_paths(jj);
            db_path = db_path{:};
            [~, dbstem, ~] = fileparts(db_path);

            candidate_render = fullfile(candidate_renders, qstem, strcat(dbstem, "_color.png"));
            if exist(candidate_render, 'file') ~= 2
                continue;
            end
            Is = imread(candidate_render);

            % desquarify images from NRIW rendering
            shape = size(Is);
            if (shape(1) == shape(2) && shape(1) > reference_h)
                square_size = shape(1);
                offset_h = fix((square_size - reference_h) / 2);
                offset_w = fix((square_size - reference_w) / 2);
                Is = Is(offset_h+1:offset_h + reference_h, offset_w+1:offset_w + reference_w, :);
            end

            % normalization
            I_synth = normalize(Is);

            % compute DSIFT
            [fq, dq] = vl_phow(im2single(Iq_norm), 'sizes', 8, 'step', 4);
            [fsynth, dsynth] = vl_phow(im2single(I_synth), 'sizes', 8, 'step', 4);
            dq = relja_rootsift(single(dq));
            dsynth = relja_rootsift(single(dsynth));

            % error
            err = sqrt(sum((dq - dsynth).^2, 1));
            score = quantile(err, 0.5)^-1;
            fprintf('%s Rendering top %d --> sim = %.4f\n', datestr(now,'HH:MM:SS.FFF'), jj, score);

            if score > this_score
                this_score = score;
                this_P = P;
                this_render = candidate_render;
                this_db = db_path;
            end
        end

        ImgList_rendered(ii).query_path = ImgList_densePE(ii).query_path;
        ImgList_rendered(ii).score = this_score;
        ImgList_rendered(ii).P = this_P;
        ImgList_rendered(ii).render_path = this_render;
        ImgList_rendered(ii).top_db_path = this_db;
    end
    save('-v6', neuralPV_matname, 'ImgList_rendered');
else
    fprintf('SKIPPING POSE VERIFICATION, output "%s" already exists.\n', neuralPV_matname);
end

function check_is_normalized(feat)
    tol = 1e-5;
    if ~all(abs(vecnorm(feat) - 1.0) < tol)
        fprintf('norm: %f\n', vecnorm(feat));
        error('Features are not normalized!');
    end
end

function name = filename(pth)
    [~, name, ext] = fileparts(pth);
    name = [name, ext];
end

function [normalized] = normalize(image)
    normalized = double(rgb2gray(image));
    if sum(normalized, 'all') ~= 0
        normalized = normalized - mean(reshape(normalized, 1, []));
        normalized = normalized ./ std(reshape(normalized, 1, []));
    end
end

function [value] = get_with_default(structure, field_name, default_value)
    % get Equivalent of Python's dict.get(name, default).
    if ~isfield(structure, field_name)
        value = default_value;
    else  % Yummy string field name access syntax:
        value = structure.(field_name);
    end
end

function create_parent_folder(filename)
    [parent_folder, ~, ~] = fileparts(filename);
    if exist(parent_folder, 'dir') ~= 7
        mkdir(parent_folder);
    end
end

function parfor_denseGV( cnnq, qname, dbname, params)
    [~, dbbasename] = fileparts(filename(dbname));
    [~, qbasename] = fileparts(filename(qname));
    coarselayerlevel = get_with_default(params, 'input_feature_layer', 6);
    finelayerlevel = 3;

    this_densegv_matname = fullfile(params.output_gv_dense_dir, filename(qname), strcat(filename(dbname), ".mat"));

    if exist(this_densegv_matname, 'file') ~= 2

        % load input feature
        features = load(fullfile(params.input_db_features_path, dbbasename + ".mat"));
        cnndb = features.features.features;

        % coarse-to-fine matching
        cnnfeat1size = size(cnnq{finelayerlevel}.x);
        cnnfeat2size = size(cnndb{finelayerlevel}.x);
        [match12,f1,f2,cnnfeat1,cnnfeat2] = at_coarse2fine_matching(cnnq,cnndb,coarselayerlevel,finelayerlevel);
        [inls12] = at_denseransac(f1,f2,match12,2);

        save('-v6', this_densegv_matname, 'cnnfeat1size', 'cnnfeat2size', 'f1', 'f2', 'inls12', 'match12');
        fprintf('denseGV: %s vs %s DONE. \n', qbasename, dbbasename);
    else
        fprintf('denseGV: %s vs %s ALREADY EXISTS. \n', qbasename, dbbasename);
    end
end

function parfor_densePE( qname, dbname, params )
    [~, dbbasename] = fileparts(filename(dbname));
    [~, qbasename] = fileparts(filename(qname));

    this_densepe_matname = fullfile(params.output_pnp_dense_inlier_dir, filename(qname), strcat(dbbasename, ".mat"));

    if exist(this_densepe_matname, 'file') ~= 2
        %geometric verification results
        this_densegv_matname = fullfile(params.output_gv_dense_dir, filename(qname), strcat(filename(dbname), ".mat"));
        this_gvresults = load(this_densegv_matname);
        tent_xq2d = this_gvresults.f1(:, this_gvresults.inls12(1, :));
        tent_xdb2d = this_gvresults.f2(:, this_gvresults.inls12(2, :));

        %depth information
        this_db_matname = fullfile(params.input_cutout_matfiles_path, strcat("cutout_", erase(filename(dbname), "_reference"), ".mat"));
        if exist(this_db_matname, 'file') ~= 2
            this_db_matname = fullfile(params.input_cutout_matfiles_path, strcat(filename(dbname), ".mat"));
        end
        load(this_db_matname, 'XYZcut');
        %load transformation matrix (local to global)
        transformation_txtname = fullfile(params.input_transforms_path, strcat("cutout_", erase(filename(dbname), "_reference"), ".mat"));
        if exist(transformation_txtname, 'file') ~= 2
            transformation_txtname = fullfile(params.input_transforms_path, strcat(filename(dbname), ".mat"));
        end
        load(transformation_txtname, 'R', 'position', 'calibration_mat');
        % P_after = eye(4);
        % P_after(1:3, 1:3) = R;
        % P_after(1:3, 4) = - R * position';
        %[ ~, P_after ] = load_WUSTL_transformation(transformation_txtname);
        %Feature upsampling
        Iqsize = size(imread(qname));
        Idbsize = size(XYZcut);
        tent_xq2d = at_featureupsample(tent_xq2d, this_gvresults.cnnfeat1size, Iqsize);
        tent_xdb2d = at_featureupsample(tent_xdb2d, this_gvresults.cnnfeat2size, Idbsize);
        %query ray
        tent_ray2d = calibration_mat(1:3, 1:3)^-1 * [tent_xq2d; ones(1, size(tent_xq2d, 2))];
        %DB 3d points
        indx = sub2ind(size(XYZcut(:,:,1)),tent_xdb2d(2,:),tent_xdb2d(1,:));
        X = XYZcut(:,:,1);Y = XYZcut(:,:,2);Z = XYZcut(:,:,3);
        tent_xdb3d = [X(indx); Y(indx); Z(indx)];
        % XYZcut is already expressed in world coordinate system, we don't need to transform here.
        % tent_xdb3d = bsxfun(@plus, P_after(1:3, 1:3)*tent_xdb3d, P_after(1:3, 4));
        %Select keypoint correspond to 3D
        idx_3d = all(~isnan(tent_xdb3d), 1);
        tent_xq2d = tent_xq2d(:, idx_3d);
        tent_xdb2d = tent_xdb2d(:, idx_3d);
        tent_ray2d = tent_ray2d(:, idx_3d);
        tent_xdb3d = tent_xdb3d(:, idx_3d);


        tentatives_2d = [tent_xq2d; tent_xdb2d];
        tentatives_3d = [tent_ray2d; tent_xdb3d];


        %solver
        if size(tentatives_2d, 2) < 3
            P = nan(3, 4);
            inls = false(1, size(tentatives_2d, 2));
        else
            [ P, inls ] = ht_lo_ransac_p3p( tent_ray2d, tent_xdb3d, 1.0*pi/180);
            if isempty(P)
                P = nan(3, 4);
            end
        end

        save('-v6', this_densepe_matname, 'P', 'inls', 'tentatives_2d', 'tentatives_3d');
        fprintf('densePE: %s vs %s DONE. \n', qbasename, dbbasename);
    else
        fprintf('densePE: %s vs %s ALREADY EXISTS. \n', qbasename, dbbasename);
    end
end

% Expects:
% params_file path to yaml
% experiment_name with name of dict in the yaml file

% To support debug in Matlab -desktop and running without an X server on the cluster.
% See dvc/scripts/inloc_pose_verification.sh where these variables are generated dynamically.
if ~exist('params_file', 'var')
    params_file = '/home/kremeto1/inloc/dvc/pipeline-inloc/params.yaml';
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
        if exist(this_render, 'file') ~= 2
            this_render = fullfile(candidate_renders, qstem, strcat(dbstem, "_out.png"));
        end

        % Just desquerify if needed
        if exist(this_render, 'file') == 2  % Try it, 6DOF could not be found
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
                candidate_render = fullfile(candidate_renders, qstem, strcat(dbstem, "_out.png"));
            end
            if exist(candidate_render, 'file') ~= 2
                continue;
            end
            Is = imread(candidate_render);

            % desquarify images from NRIW rendering
            shape = size(Is);
            if (shape(1) == shape(2) && reference_w ~= reference_h)
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

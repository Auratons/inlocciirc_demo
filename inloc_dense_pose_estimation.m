function inloc_dense_pose_estimation(varargin)
    % Note: It first rerank top100 original shortlist (ImgList_retrieval) in accordance
    % with the number of dense matching inliers. It then computes query
    % candidate poses by using top10 database images.

    [filepath, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'functions', 'inLocCIIRC_utils'));
    inloc_add_abs_fn_path('yaml');
    inloc_add_abs_fn_path('wustl_function');

    params = inloc_parse_inputs(varargin{:}).pose_estimation;

    % densePE (topN reranking -> pnp_topN pose candidate)
    densePE_matname = params.output_candidate_mat_path;
    if exist(densePE_matname, 'file') ~= 2
        pnp_topN = get_with_default(params, 'input_pnp_topN', 10);
        ImgList_retrieval = load(params.input_topN_mat_path, 'ImgList');
        ImgList_retrieval = ImgList_retrieval.ImgList;
        query_features = load(params.input_query_features_mat_path, 'features');
        query_features = query_features.features;
        db_features = load(params.input_db_features_mat_path, 'features');
        db_features = db_features.features;

        ImgList = struct('query_path', {}, 'topN_db_paths', {}, 'topN_scores', {}, 'P', {});

        % load(params.input_db_depth_mat_path, 'db_depthnames_all');
        for ii = 1:1:1%length(ImgList_retrieval)
            curr_query_path = ImgList_retrieval(ii).query_path;
            curr_topN_db_paths = ImgList_retrieval(ii).topN_db_paths;

            ImgList(ii).query_path = curr_query_path;
            ImgList(ii).topN_db_paths = ImgList_retrieval(ii).topN_db_paths;
            % ImgList(ii).topN_scores = {};

            % preload query feature
            idx = find(strcmp(convertCharsToStrings(curr_query_path), {query_features.img_path}));
            cnnq = query_features(idx).features;

            if exist(fullfile(params.output_gv_dense_dir, filename(curr_query_path)), 'dir') ~= 7
                mkdir(fullfile(params.output_gv_dense_dir, filename(qname)));
            end

            parfor kk = 1:1:length(ImgList_retrieval(ii).topN_db_paths)
                parfor_denseGV( cnnq, curr_query_path, curr_topN_db_paths{kk}, params );
                fprintf('dense matching: %s vs %s DONE. \n', curr_query_path, curr_topN_db_paths{kk});
            end

            for jj = 1:1:length(ImgList_retrieval(ii).topN_db_paths)
                this_gvresults = load(fullfile( ...
                    params.output_gv_dense_dir, ...
                    filename(curr_query_path), ...
                    strcat(filename(curr_topN_db_paths{jj}), ".mat")));
                ImgList(ii).topN_scores(jj) = ImgList_retrieval(ii).topN_scores(jj) + size(this_gvresults.inls12, 2);
            end

            [sorted_score, idx] = sort(ImgList(ii).topN_scores, 'descend');
            ImgList(ii).topN_db_paths = ImgList(ii).topN_db_paths(idx);
            ImgList(ii).topN_scores = ImgList(ii).topN_scores(idx);

            fprintf('%s done. \n', curr_query_path);
        end

        % pnp list
        qlist = cell(1, length(ImgList) * pnp_topN);
        dblist = cell(1, length(ImgList) * pnp_topN);
        for ii = 1:1:length(ImgList)
            for jj = 1:1:pnp_topN
                qlist{pnp_topN * (ii - 1) + jj} = ImgList(ii).query_path;
                dblist{pnp_topN * (ii - 1) + jj} = ImgList(ii).topN_db_paths{jj};
            end
        end

        % dense pnp
        for ii = 1:1:length(ImgList)
            if exist(fullfile(params.output_pnp_dense_inlier_dir, filename(ImgList(ii).query_path), 'dir') ~= 7
                mkdir(fullfile(params.output_pnp_dense_inlier_dir, filename(ImgList(ii).query_path));
            end
            fprintf("Creating %s\n", fullfile(params.output_pnp_dense_inlier_dir, filename(ImgList(ii).query_path)));
        end
        parfor ii = 1:1:length(qlist)
            parfor_densePE( ...
                qlist{ii}, ...
                dblist{ii}, ...
                ... %db_depthnames_all{ii}, ...
                params ...
            );
            fprintf('densePE: %s vs %s DONE. \n', qlist{ii}, dblist{ii});
        end

        % load top10 pnp
        for ii = 1:1:length(ImgList)
            ImgList(ii).P = cell(1, pnp_topN);
            for jj = 1:1:pnp_topN
                [~, dbbasename, ~] = fileparts(ImgList(ii).topN_db_paths{jj});
                this_densepe_matname = fullfile(params.output_pnp_dense_inlier_dir, filename(ImgList(ii).query_path), strcat(dbbasename, ".mat"));
                load(this_densepe_matname, 'P');
                ImgList(ii).P{jj} = P;
            end
        end

        create_parent_folder(densePE_matname);
        save('-v6', densePE_matname, 'ImgList');
    end
end

function name = filename(pth)
    [~, name, ext] = fileparts(pth);
    name = [name, ext];
end

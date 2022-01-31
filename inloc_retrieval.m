function inloc_retrieval(varargin)
    % Config fields:
    % input_scores_mat_path
    % input_topN
    % output_topN_mat_path
    % It loads localization score and output topN (default 100) database list for each query.

    [filepath, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'functions', 'inLocCIIRC_utils'));
    inloc_add_abs_fn_path('yaml');

    params = inloc_parse_inputs(varargin{:}).retrieval;

    if exist(params.output_topN_mat_path, 'file') ~= 2
        shortlist_topN = get_with_default(params, 'input_topN', 100);
        ImgList = struct('query_path', {}, 'topN_db_paths', {}, 'topN_scores', {});

        load(params.input_scores_mat_path, 'score');

        for i = 1:1:size(score, 2)
            ImgList(i).query_path = score(i).query_path;
            [score_sort, score_idx] = sort(score(i).scores, 'descend');
            ImgList(i).topN_db_paths = score(i).db_score_paths(score_idx(1:shortlist_topN));
            ImgList(i).topN_scores = score_sort(1:shortlist_topN);
        end

        create_parent_folder(params.output_topN_mat_path);
        save('-v6', params.output_topN_mat_path, 'ImgList');
    end
end

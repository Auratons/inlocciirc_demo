function inloc_dense_pose_verification(varargin)
    % densePV (top10 pose candidate -> pose verification)
    % Modified pose verification step to run with rendered images.
    % and compute similarity between original query and synthesized views.
    % Pose candidates are then re-scored by the similarity.

    [filepath, ~, ~] = fileparts(mfilename('fullpath'));
    addpath(fullfile(filepath, '..', 'functions', 'inLocCIIRC_utils'));
    inloc_add_abs_fn_path('yaml');
    inloc_add_abs_fn_path('utils');
    run(fullfile(filepath, '..', 'functions', 'vlfeat', 'toolbox', 'vl_setup.m'));

    params = inloc_parse_inputs(varargin{:}).pose_verification;

    PV_topN = get_with_default(params, 'input_pv_topN', 10);

    neuralPV_matname = params.output_shortlist_mat_path;
    if exist(neuralPV_matname, 'file') ~= 2
        %Load the poses from the .mat file
        load(params.input_candidate_mat_path, 'ImgList');
        ImgList_rendered = struct('query_path', {}, 'score', {}, 'P', {});

        for ii = 1:1:length(ImgList)
            qpath = ImgList(ii).query_path;
            fprintf('\nProcessing %s\n', qpath) % debug
            %Load query image & output information
            Iq = imread(qpath);
            reference_size = size(Iq);
            reference_h = reference_size(1);
            reference_w = reference_size(2);

            % normalization
            Iq_norm = normalize(Iq);

            this_P = nan(3, 4);
            this_score = 0;
            for jj = 1:1:PV_topN
                % Strip the rank of the rendering from the name of the image
                P = ImgList(ii).P(jj);
                db_path = ImgList(ii).topN_db_paths(jj);
                render_path = strcat(erase(db_path, "_reference.png"), "_color.png");
                I_synth = imread(render_path);
                % desquarify images from NRIW rendering
                shape = size(I_synth);
                if (shape(1) == shape(2))
                    square_size = shape(1);
                    offset_h = fix((square_size - reference_h) / 2);
                    offset_w = fix((square_size - reference_w) / 2);
                    I_synth = I_synth(offset_h+1:offset_h + reference_h, offset_w+1:offset_w + reference_w, :);
                end

                % normalization
                I_synth = normalize(I_synth);

                % compute DSIFT
                [fq, dq] = vl_phow(im2single(Iq_norm), 'sizes', 8, 'step', 4);
                [fsynth, dsynth] = vl_phow(im2single(I_synth), 'sizes', 8, 'step', 4);
                dq = relja_rootsift(single(dq));
                dsynth = relja_rootsift(single(dsynth));

                % error
                err = sqrt(sum((dq - dsynth).^2, 1));
                score = quantile(err, 0.5)^-1;
                fprintf('rendering top %d --> sim = %.4f\n', jj, score);

                if score > this_score
                    this_score = score;
                    this_P = P;
                end
            end

            ImgList_rendered(ii).query_path = ImgList(ii).query_path;
            ImgList_rendered(ii).score = this_score;
            ImgList_rendered(ii).P = this_P;
        end
        save('-v6', neuralPV_matname, 'ImgList_rendered');
    else
        fprintf('SKIPPING POSE VERIFICATION, output "%s" already exists.\n', neuralPV_matname);
    end
end

function [normalized] = normalize(image)
    normalized = double(rgb2gray(image));
    normalized = normalized - mean(reshape(normalized, 1, []));
    normalized = normalized ./ std(reshape(normalized, 1, []));
end

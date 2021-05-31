%Modified pose verification step to run with rendered images. 
%and compute similarity between original query and synthesized views. Pose
%candidates are then re-scored by the similarity. 

%% densePV (top10 pose candidate -> pose verification)
PV_topN = 10;
%Load the poses from the .mat file
load(params.input.topK_matfile, 'ImgList');
ImgList_rendered = struct('queryname', {}, 'score', {}, 'P', {});
neuralPV_matname = fullfile(params.output.dir, params.output.matname);

if exist(neuralPV_matname, 'file') ~= 2
    for ii = 1:1:length(ImgList)
        qname = ImgList(ii).queryname;
        fprintf('\nProcessing %s\n', qname) % debug
        if exist(fullfile(params.input.dir, qname), 'dir') == 7
            %Load query image & output information
            P_list = ImgList(ii).P(1:PV_topN);
            Iq = imread(fullfile(params.input.dir, qname, strcat('reference', params.input.imgformat)));
            reference_size = size(Iq);
            reference_h = reference_size(1);
            reference_w = reference_size(2);

            %normalization
            Iq_norm = double(rgb2gray(Iq));
            Iq_norm = Iq_norm - mean(reshape(Iq_norm, 1, []));
            Iq_norm = Iq_norm ./ std(reshape(Iq_norm, 1, []));

            list_rendered = dir(fullfile(params.input.dir, qname, strcat('*out*', params.input.imgformat)));

            this_P = nan(3, 4);
            this_score = 0;
            if length(list_rendered) ~= 0
                for jj = 1:1:length(list_rendered)
                    %Strip the rank of the rendering from the name of the image
                    img_name = list_rendered(jj).name;
                    img_index = img_name(1:end - length(params.input.imgformat) - length('_out'));
                    img_index = str2num(img_index) + 1;
                    P = P_list(img_index);
                    img_path = list_rendered(jj);
                    I_synth = imread(fullfile(img_path.folder, img_path.name));
                    % desquarify images from NRIW rendering
                    shape = size(I_synth);
                    if (shape(1) == shape(2))
                        square_size = shape(1);
                        offset_h = fix((square_size - reference_h) / 2);
                        offset_w = fix((square_size - reference_w) / 2);
                        I_synth = I_synth(offset_h+1:offset_h + reference_h, offset_w+1:offset_w + reference_w, :);
                    end

                    %normalization
                    I_synth = double(rgb2gray(I_synth));
                    I_synth = I_synth - mean(reshape(I_synth, 1, []));
                    I_synth = I_synth ./ std(reshape(I_synth, 1, []));
                    
                    %compute DSIFT
                    [fq, dq] = vl_phow(im2single(Iq_norm), 'sizes', 8, 'step', 4);
                    [fsynth, dsynth] = vl_phow(im2single(I_synth), 'sizes', 8, 'step', 4);
                    dq = relja_rootsift(single(dq)); 
                    dsynth = relja_rootsift(single(dsynth));
                    
                    %error
                    err = sqrt(sum((dq - dsynth).^2, 1));
                    score = quantile(err, 0.5)^-1;
                    fprintf('rendering %d --> sim = %.4f\n', img_index, score);

                    if score > this_score
                        this_score = score;
                        this_P = P;
                    end
                end
            end

            % errmap = err;
            % xuni = sort(unique(fsynth(1, :)), 'ascend');
            % yuni = sort(unique(fsynth(2, :)), 'ascend');
            % errmap = errmap(yuni, xuni);

            ImgList_rendered(ii).queryname = ImgList(ii).queryname;
            ImgList_rendered(ii).score = this_score;
            ImgList_rendered(ii).P = this_P;
        end
    end
    save('-v6', neuralPV_matname, 'ImgList_rendered');
else
    load(neuralPV_matname, 'ImgList_rendered');
end

if exist(params.output.dir, 'dir') ~= 7
    mkdir(params.output.dir);
end

% Save in submission format
fid = fopen(fullfile(params.output.dir, params.output.txtname), 'w');
compt = 0;
% load('~/InLoc_demo/outputs/densePE_top100_shortlist.mat', 'ImgList');
for ii = 1:1:length(ImgList_rendered)
    this_qname = ImgList_rendered(ii).queryname;
    if length(this_qname) > 0
        if length(ImgList_rendered(ii).P) == 0
            fprintf('Skipping submission format printing for %s\n', ImgList_rendered(ii).queryname)
            continue
        end
        compt = compt + 1;
        fprintf('Submission format printing for %s\n', ImgList_rendered(ii).queryname)
        this_P = ImgList_rendered(ii).P{1};
        % this_P = ImgList(ii).P{1};
        qtr = rot2qtr(this_P(1:3, 1:3));
        t = this_P(1:3, 4);
        
        if sum(isnan(qtr)) > 0 || sum(isinf(qtr)) > 0 || sum(isnan(t)) > 0 || sum(isinf(t)) > 0
            continue;
        end
        
        %1. image name
        fprintf(fid, '%s', this_qname);
        %2. pose
        fprintf(fid, ' %.10g', qtr);
        fprintf(fid, ' %.10g', t);
        
        fprintf(fid, '\n');
    end
end
fprintf('\nWritten %d / %d poses in submission format.\n', compt, length(ImgList))

fclose(fid);

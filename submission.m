% Submission format
params.input.list = '/home/kremeto1/inloc/datasets/pipeline-inloc-manual-pyrender/densePV_top10_shortlist_neural.mat';
params.output.dir = '/home/kremeto1/inloc/datasets/pipeline-inloc-manual-pyrender';
params.output.txtname = 'pyrender-neural-submission.txt';

ImgList = load(params.input.list);
ImgList_rendered = ImgList.ImgList_rendered;

fid = fopen(fullfile(params.output.dir, params.output.txtname), 'w');
compt = 0;
for ii = 1:1:length(ImgList_rendered)
    this_qname = ImgList_rendered(ii).query_path;
    if length(this_qname) > 0
        compt = compt + 1;
        this_P = ImgList_rendered(ii).P;
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
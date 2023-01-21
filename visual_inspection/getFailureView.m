function im_synth = getFailureView(params,ImgList,q_id,db_rank,showMontage,savePth,prefix)
    qpath =  ImgList(q_id).query_path;

    % cutoutName = ImgList(q_id).topNname;
    % cutoutName= cutoutName{db_rank};

    dbpath = ImgList(q_id).top_db_path;
    Iq = imread(qpath);

    im_synth = [];
    %dense features

    this_densegv_matname = fullfile(params.output.gv_dense.dir, filename(qpath), strcat(filename(dbpath), ".mat"));
    gvresults = load(this_densegv_matname);

    tent_xq2d = gvresults.f1(:, gvresults.match12(1, :));
    tent_xdb2d = gvresults.f2(:, gvresults.match12(2, :)); %matches

    inls_xq2d = gvresults.f1(:, gvresults.inls12(1, :));
    inls_xdb2d = gvresults.f2(:, gvresults.inls12(2, :)); %matches

    %Feature upsampling
    Idb = imread(dbpath);
    Idbsize = size(Idb);
    Iqsize = size(Iq);
    tent_xq2d = at_featureupsample(tent_xq2d,gvresults.cnnfeat1size,Iqsize);
    % without this, the features in query image would not match the cutout aspect ratio
    tent_xdb2d = at_featureupsample(tent_xdb2d,gvresults.cnnfeat2size,Idbsize);
    inls_xq2d = at_featureupsample(inls_xq2d,gvresults.cnnfeat1size,Iqsize);
    % without this, the features in query image would not match the cutout aspect ratio
    inls_xdb2d = at_featureupsample(inls_xdb2d,gvresults.cnnfeat2size,Idbsize);

    queryWidth = Iqsize(2);
    queryHeight = Iqsize(1);
    cutoutWidth = Idbsize(2);
    cutoutHeight = Idbsize(1);

    inls_xq2d = adjust_inliers_to_match_original_query(inls_xq2d, queryWidth, queryHeight, cutoutWidth, cutoutHeight);
    tent_xq2d = adjust_inliers_to_match_original_query(tent_xq2d, queryWidth, queryHeight, cutoutWidth, cutoutHeight);
    %     together_q = repmat(together, [1, 1, 3]);
    %     for feat_i = 1: size(tent_xq2d,2)
    %         color = [0 0 255];
    %         block = zeros(3,3,3);
    %         block(:,:,1:)
    %         together_q(tent_xq2d(1,feat_i)-1:tent_xq2d(1,feat_i)+1,tent_xq2d(2,feat_i)-1:tent_xq2d(2,feat_i)+1,:) = ;
    %     end
    f = figure('visible', 'off');
    imshow(rgb2gray([Iq Idb]));hold on;
    for c = 1:20:size(tent_xq2d,2)
        plot3d([tent_xq2d(:,c), tent_xdb2d(:,c) + [size(Iq,2) 0 ]'],'r-');
    end
    for c = 1:20:size(inls_xq2d,2)
        plot3d([inls_xq2d(:,c), inls_xdb2d(:,c) + [size(Iq,2) 0 ]'],'-g');
    end
    plot3d(tent_xq2d,'b.');
    plot3d(inls_xq2d,'g.');
    plot3d(tent_xdb2d + [size(Iq,2) 0 ]','b.');
    plot3d(inls_xdb2d + [size(Iq,2) 0 ]','g.');
    if showMontage
        %                 errmaps =load(fullfile(params.output.synth.dir, ImgList(q_id).queryname, sprintf('%d%s', db_rank, params.output.synth.matformat)),'errmaps');
        %                 errmaps= grs2rgb(errmaps.errmaps{1});
        if nargin >= 6
            if nargin >=7
                saveas(gcf,fullfile(savePth,sprintf('%s_results_q_id_%d_best_db_%d.jpg',prefix,q_id,db_rank)));
            else
                saveas(gcf,fullfile(savePth,sprintf('results_q_id_%d_best_db_%d.jpg',q_id,db_rank)));
            end
        end
    end
    close(f);
end

function name = filename(pth)
    [~, name, ext] = fileparts(pth);
    name = [name, ext];
end

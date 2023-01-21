function im_synth = getSynthView(params, ImgList, q_id, db_rank, showMontage, savePth, prefix)
    qpath =  ImgList(q_id).query_path;
    dbpath = ImgList(q_id).top_db_path;
    Iq = imread(qpath);
    RGBpersp = imread(ImgList(q_id).render_path);

    shape = size(RGBpersp); render_h = shape(1); render_w = shape(2);
    shape = size(Iq); iq_h = shape(1); iq_w = shape(2);
    if render_h ~= render_w && iq_h == iq_w
        square_size = iq_h;
        offset_h = fix((square_size - render_h) / 2);
        offset_w = fix((square_size - render_w) / 2);
        Iq = Iq(offset_h+1:offset_h + render_h, offset_w+1:offset_w + render_w, :);
    end

    if showMontage
        blend = Iq/2 + RGBpersp/2;
        comparison = {Iq, RGBpersp, blend, imread(dbpath)};
        f = figure('visible', 'off');
        montage(comparison,'Size',[1 4]);
        if nargin >= 6
            if nargin >=7
                if exist(savePth, 'dir') ~= 7
                    mkdir(savePth);
                end
                saveas(gcf, fullfile(savePth, sprintf('%s_results_q_id_%d_best_db_%d.jpg', prefix, q_id, db_rank)));
            else
            saveas(gcf,fullfile(savePth, sprintf('results_q_id_%d_best_db_%d.jpg', q_id, db_rank)));
            end
        end
        close(f);

    end
    im_synth = RGBpersp;
end

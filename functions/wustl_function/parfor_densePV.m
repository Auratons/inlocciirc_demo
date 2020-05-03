function parfor_densePV( qname, dbname, P, params )

spaceName = strsplit(dbname, '/');
spaceName = spaceName{1};
this_densePV_matname = fullfile(params.output.synth.dir, qname, buildCutoutName(dbname, params.output.synth.matformat));

if exist(this_densePV_matname, 'file') ~= 2
    if all(~isnan(P(:)))
        
        %load downsampled images
        Iq = imresize(imread(fullfile(params.data.dir, params.data.q.dir, qname)), params.data.q.dslevel);
        fl = params.data.q.fl * params.data.q.dslevel;
        R = P(1:3,1:3);
        t = P(1:3,4);

        meshPath = fullfile(params.data.models.dir, spaceName, 'mesh_rotated.obj');
        t = -inv(R)*t;
        rFix = [180.0, 0.0, 0.0];
        Rfix = rotationMatrix(deg2rad(rFix), 'XYZ');
        sensorSize = [size(Iq,2), size(Iq,1)];
        headless = ~strcmp(environment(), 'laptop');
        [RGBpersp, XYZpersp, depth] = projectMesh(meshPath, fl, inv(R)*Rfix, t, sensorSize, false, -1, params.input.projectMesh_py_path, headless);
        RGB_flag = all(~isnan(XYZpersp), 3);
        
        %compute DSIFT error
        if any(RGB_flag(:))
            %normalization
            Iq_norm = image_normalization( double(rgb2gray(Iq)), RGB_flag );
            I_synth = double(rgb2gray(RGBpersp));
            I_synth(~RGB_flag) = nan;
            I_synth = image_normalization( inpaint_nans(I_synth), RGB_flag );
            
            %compute DSIFT
            [fq, dq] = vl_phow(im2single(Iq_norm),'sizes',8,'step',4);
            [fsynth, dsynth] = vl_phow(im2single(I_synth),'sizes',8,'step',4);
            f_linind = sub2ind(size(I_synth), fsynth(2, :), fsynth(1, :));
            iseval = RGB_flag(f_linind);
            dq = relja_rootsift(single(dq)); dsynth = relja_rootsift(single(dsynth));
            
            %error
            err = sqrt(sum((dq(:, iseval) - dsynth(:, iseval)).^2, 1));
            score = quantile(err, 0.5)^-1;
            errmap = nan(size(I_synth));errmap(f_linind(iseval)) = err;
            xuni = sort(unique(fsynth(1, :)), 'ascend');yuni = sort(unique(fsynth(2, :)), 'ascend');
            errmap = errmap(yuni, xuni);
            
%             %debug
%             figure();set(gcf, 'Position', [0 0 1000, 300]);
%             ultimateSubplot( 3, 1, 1, 1, 0.01, 0.05 );
%             imshow(Iq);
%             ultimateSubplot( 3, 1, 2, 1, 0.01, 0.05 );
%             imshow(RGBpersp);
%             ultimateSubplot( 3, 1, 3, 1, 0.01, 0.05 );
%             imagesc(errmap);colormap('jet');axis image off;
%             keyboard;
            
        else
            score = 0;
            errmap = [];
        end
    else
        Iq = [];
        RGBpersp = [];
        RGB_flag = [];
        score = 0;
        errmap = 0;
    end
    
    if exist(fullfile(params.output.synth.dir, qname), 'dir') ~= 7
        mkdir(fullfile(params.output.synth.dir, qname));
    end
    save(this_densePV_matname, 'Iq', 'RGBpersp', 'RGB_flag', 'score', 'errmap');
end


end


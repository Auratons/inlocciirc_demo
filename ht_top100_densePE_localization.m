%Note: It first rerank top100 original shortlist (ImgList_original) in accordance
%with the number of dense matching inliers. TODO: and then?

shortlist_topN = 100;
topN_with_GV = 10;
mCombinations = 10;

%% densePE (top100 reranking -> top10 pose candidate)

densePE_matname = fullfile(params.output.dir, 'densePE_top100_shortlist.mat');
if exist(densePE_matname, 'file') ~= 2
    
    %dense feature extraction
    net = load(params.netvlad.dataset.pretrained);
    net = net.net;
    net= relja_simplenn_tidy(net);
    net= relja_cropToLayer(net, 'postL2');
    for ii = 1:1:length(ImgList_original)
        q_densefeat_matname = fullfile(params.input.feature.dir, params.dataset.query.dirname, [ImgList_original(ii).queryname, params.input.feature.q_matformat]);
        if exist(q_densefeat_matname, 'file') ~= 2
            queryImage = load_query_image_compatible_with_cutouts(fullfile(params.dataset.query.dir, ImgList_original(ii).queryname), ...
                                                                        params.dataset.db.cutout.size);
            cnn = at_serialAllFeats_convfeat(net, queryImage, 'useGPU', true);
            cnn{1} = [];
            cnn{2} = [];
            cnn{4} = [];
            cnn{6} = [];
            [feat_path, ~, ~] = fileparts(q_densefeat_matname);
            if exist(feat_path, 'dir')~=7; mkdir(feat_path); end
            save('-v6', q_densefeat_matname, 'cnn');
            fprintf('Dense feature extraction: %s done. \n', ImgList_original(ii).queryname);
        end
        
        for jj = 1:1:shortlist_topN
            db_densefeat_matname = fullfile(params.input.feature.dir, params.dataset.db.cutout.dirname, ...
                [ImgList_original(ii).topNname{jj}, params.input.feature.db_matformat]);
            if exist(db_densefeat_matname, 'file') ~= 2
                cutoutImage = imread(fullfile(params.dataset.db.cutouts.dir, ImgList_original(ii).topNname{jj}));
                cnn = at_serialAllFeats_convfeat(net, cutoutImage, 'useGPU', true);
                cnn{1} = [];
                cnn{2} = [];
                cnn{4} = [];
                cnn{6} = [];
                [feat_path, ~, ~] = fileparts(db_densefeat_matname);
                if exist(feat_path, 'dir')~=7; mkdir(feat_path); end
                save('-v6', db_densefeat_matname, 'cnn');
                fprintf('Dense feature extraction: %s done. \n', ImgList_original(ii).topNname{jj});
            end
        end
    end

    inloc_hw = getenv("INLOC_HW");
    if strcmp(inloc_hw, "GPU")
        exit(0);
    end
    
    %shortlist reranking
    ImgList = struct('queryname', {}, 'topNname', {}, 'topNscore', {}, 'Ps', {});
    for ii = 1:1:length(ImgList_original)
        ImgList(ii).queryname = ImgList_original(ii).queryname;
        ImgList(ii).topNname = ImgList_original(ii).topNname(1:shortlist_topN);
        
        %preload query feature
        qfname = fullfile(params.input.feature.dir, params.dataset.query.dirname, [ImgList(ii).queryname, params.input.feature.q_matformat]);
        cnnq = load(qfname, 'cnn');cnnq = cnnq.cnn;
        
        parfor kk = 1:1:shortlist_topN
            parfor_denseGV( cnnq, ImgList(ii).queryname, ImgList(ii).topNname{kk}, params );
            fprintf('dense matching: %s vs %s DONE. \n', ImgList(ii).queryname, ImgList(ii).topNname{kk});
        end
        
        for jj = 1:1:shortlist_topN
            cutoutPath = ImgList(ii).topNname{jj};
            this_gvresults = load(fullfile(params.output.gv_dense.dir, ImgList(ii).queryname, buildCutoutName(cutoutPath, params.output.gv_dense.matformat)));
            ImgList(ii).topNscore(jj) = ImgList_original(ii).topNscore(jj) + size(this_gvresults.inls12, 2);
        end
        
        [sorted_score, idx] = sort(ImgList(ii).topNscore, 'descend');
        ImgList(ii).topNname = ImgList(ii).topNname(idx);
        ImgList(ii).topNscore = ImgList(ii).topNscore(idx);

        fprintf('%s done. \n', ImgList(ii).queryname);
    end
    
    %% for each query, find top-mCombinations sequences of lengths params.sequence.length
    areQueriesFromHoloLensSequence = isfield(params, 'sequence') && isfield(params.sequence, 'length');
    if ~areQueriesFromHoloLensSequence
        desiredSequenceLength = 1;
    else
        desiredSequenceLength = params.sequence.length;
    end
    ImgListSequential = ImgList;

    % build queryInd (i-th query in ImgList does not mean i-th query in the whole sequence)
    queryInd = zeros(length(ImgList),1);
    for i=1:length(ImgList)
        queryIdx = ImgList(i).queryname;
        queryIdx = strsplit(queryIdx, '.');
        queryIdx = queryIdx{1};
        queryIdx = str2num(queryIdx);
        queryInd(i) = queryIdx;
    end
    [~,queryInd] = sort(queryInd);
    % TODO: assert difference between neighbouring IDs is 1

    for queryId=1:length(ImgList)
        i = queryInd(queryId); % position of queryId in ImgList
        % compute cumulative score for each combination

        % generate all combination indices
        if queryId-desiredSequenceLength+1 < 1
            actualSequenceLength = queryId;
        else
            actualSequenceLength = desiredSequenceLength;
        end
        permInd = permn([1:topN_with_GV], actualSequenceLength);

        permScores = zeros(size(permInd,1),1);
        for j=1:size(permInd)
            score = 0.0;
            permIndCol = 0;
            for k=queryId-actualSequenceLength+1:queryId
                permIndCol = permIndCol + 1;
                cutoutIdx = permInd(j,permIndCol);
                score = score + ImgList(k).topNscore(cutoutIdx);
            end
            permScores(j) = score;
        end

        % find indices of m sequences with the highest cumulative score
        [sorted_score, idx] = sort(permScores, 'descend');
        ImgListSequential(i).topNscore = sorted_score(idx(1:mCombinations))';

        topInd = permInd(idx(1:mCombinations),:);
        ImgListSequential(i).topNname = cell(actualSequenceLength,mCombinations);
        for j=1:size(topInd,1)
            permIndCol = 0;
            for k=queryId-actualSequenceLength+1:queryId
                permIndCol = permIndCol + 1;
                name = ImgList(i).topNname{topInd(j,permIndCol)};
                ImgListSequential(i).topNname{permIndCol,j} = name;
            end
        end
    end

    if areQueriesFromHoloLensSequence
        posesFromHoloLens = getPosesFromHoloLens(params.HoloLensOrientationDelay, params.HoloLensTranslationDelay, params);
        nQueries = length(ImgListSequential);
        assert(size(posesFromHoloLens,1) == nQueries);
    end

    qlist = cell(1, length(ImgListSequential)*mCombinations);
    dblist = cell(1, length(ImgListSequential)*mCombinations);
    dbind = cell(1, length(ImgListSequential)*mCombinations);
    posesFromHoloLensList = cell(1, length(ImgListSequential)*mCombinations);
    firstQueryInd = cell(1, length(ImgListSequential)*mCombinations);
    lastQueryInd = cell(1, length(ImgListSequential)*mCombinations);
    for ii = 1:length(ImgListSequential)
        lastQueryId = find(queryInd == ii); % the one for which we try to estimate pose
        for jj = 1:mCombinations
            idx = mCombinations*(ii-1)+jj;
            qlist{idx} = ImgListSequential(ii).queryname;
            dblist{idx} = ImgListSequential(ii).topNname(:,jj);
            dbind{idx} = jj;
            actualSequenceLength = size(ImgListSequential(ii).topNname, 1);
            firstQueryId = lastQueryId - actualSequenceLength + 1;
            if areQueriesFromHoloLensSequence
                posesFromHoloLensList{idx} = posesFromHoloLens(firstQueryId:lastQueryId,:,:);
            end
            firstQueryInd{idx} = firstQueryId;
            lastQueryInd{idx} = lastQueryId;
        end
    end

    %dense pnp
    parfor ii = 1:length(qlist)
    %for ii = 1:length(qlist)
        parfor_densePE(qlist{ii}, dblist{ii}, dbind{ii}, posesFromHoloLensList{ii}, firstQueryInd{ii}, lastQueryInd{ii}, params);
        fprintf('densePE: %s vs a cutout sequence DONE. \n', qlist{ii});
        fprintf('%d/%d done.\n', ii, length(qlist));
    end
    
    %load top-mCombinations poses
    for ii = 1:1:length(ImgListSequential)
        ImgListSequential(ii).Ps = cell(1, mCombinations);
        for jj = 1:1:mCombinations
            this_densepe_matname = fullfile(params.output.pnp_dense_inlier.dir, ImgListSequential(ii).queryname, ...
                                            sprintf('%d%s', jj, params.output.pnp_dense.matformat));
            load(this_densepe_matname, 'Ps');
            ImgListSequential(ii).Ps{jj} = Ps;
        end
    end
    
    if exist(params.output.dir, 'dir') ~= 7
        mkdir(params.output.dir);
    end
    ImgList = ImgListSequential;
    save('-v6', densePE_matname, 'ImgList');
else
    load(densePE_matname, 'ImgList');
end
ImgList_densePE = ImgList;

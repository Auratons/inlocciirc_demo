%Note: It first synthesize query views according to top10 pose candedates
%and compute similarity between original query and synthesized views. Pose
%candidates are then re-scored by the similarity.

%% densePV (top10 pose candidate -> pose verification)
PV_topN = 10;
densePV_matname = fullfile(params.output.dir, 'densePV_top10_shortlist.mat');
if exist(densePV_matname, 'file') ~= 2
    
    %synthesis list
    qlist = cell(1, PV_topN*length(ImgList_densePE));
    dblist = cell(1, PV_topN*length(ImgList_densePE));
    Plist = cell(1, PV_topN*length(ImgList_densePE));
    for ii = 1:1:length(ImgList_densePE)
        for jj = 1:1:PV_topN
            qlist{PV_topN*(ii-1)+jj} = ImgList_densePE(ii).queryname;
            dblist{PV_topN*(ii-1)+jj} = ImgList_densePE(ii).topNname{jj};
            Plist{PV_topN*(ii-1)+jj} = ImgList_densePE(ii).P{jj};
        end
    end
    %find unique scans
    dbscanlist = cell(size(dblist));
    dbscantranslist = cell(size(dblist));
    for ii = 1:1:length(dblist)
        this_floorid = strsplit(dblist{ii}, '/');this_floorid = this_floorid{1};
        info = parse_WUSTL_cutoutname( dblist{ii} );
        dbscanlist{ii} = strcat(this_floorid, params.dataset.db.scan.matformat);
        dbscantranslist{ii} = fullfile(this_floorid, 'transformations', ['trans_', info.scan_id, '.txt']);
    end
    [dbscanlist_uniq, sort_idx, uniq_idx] = unique(dbscanlist);
    dbscantranslist_uniq = cell(size(dbscanlist_uniq));
    qlist_uniq = cell(size(dbscanlist_uniq));
    dblist_uniq = cell(size(dbscanlist_uniq));
    Plist_uniq = cell(size(dbscanlist_uniq));
    for ii = 1:1:length(dbscanlist_uniq)
        idx = uniq_idx == ii;
        dbscantranslist_uniq{ii} = dbscantranslist(idx);
        qlist_uniq{ii} = qlist(idx);
        dblist_uniq{ii} = dblist(idx);
        Plist_uniq{ii} = Plist(idx);
    end
    
    %compute synthesized views and similarity

    % Because projectMesh in densePV requires up to 20 GB of RAM per one instance,
    % we need to limit the number of workers
    % TODO: optimize and leverage more workers
    poolobj = gcp('nocreate');
    delete(poolobj); % terminate any previous pool
    nWorkers = 4;
    c = parcluster;
    c.NumWorkers = nWorkers;
    saveProfile(c);
    p = parpool('local', nWorkers);

    for ii = 1:1:length(dbscanlist_uniq)
        this_dbscan = dbscanlist_uniq{ii};
        this_dbscantrans = dbscantranslist_uniq{ii};
        this_qlist = qlist_uniq{ii};
        this_dblist = dblist_uniq{ii};
        this_Plist = Plist_uniq{ii};
        
        %compute synthesized images and similarity scores
        parfor jj = 1:1:length(this_qlist)
            P = load_CIIRC_transformation(fullfile(params.dataset.db.trans.dir, this_dbscantrans{jj}));
            parfor_densePV( this_qlist{jj}, this_dblist{jj}, this_Plist{jj} * P, params );
            fprintf('densePV: %d / %d done. \n', jj, length(this_qlist));
        end
        fprintf('densePV: scan %s (%d / %d) done. \n', this_dbscan, ii, length(dbscanlist_uniq));
    end
    
    %load similarity score and reranking
    ImgList = struct('queryname', {}, 'topNname', {}, 'topNscore', {}, 'P', {});
    for ii = 1:1:length(ImgList_densePE)
        ImgList(ii).queryname = ImgList_densePE(ii).queryname;
        ImgList(ii).topNname = ImgList_densePE(ii).topNname(1:PV_topN);
        ImgList(ii).topNscore = zeros(1, PV_topN);
        ImgList(ii).P = ImgList_densePE(ii).P(1:PV_topN);
        for jj = 1:1:PV_topN
            cutoutPath = ImgList(ii).topNname{jj};
            load(fullfile(params.output.synth.dir, ImgList(ii).queryname, buildCutoutName(cutoutPath, params.output.synth.matformat)), 'score');
            ImgList(ii).topNscore(jj) = score;
        end
        
        %reranking
        [sorted_score, idx] = sort(ImgList(ii).topNscore, 'descend');
        ImgList(ii).topNname = ImgList(ii).topNname(idx);
        ImgList(ii).topNscore = ImgList(ii).topNscore(idx);
        ImgList(ii).P = ImgList(ii).P(idx);
    end
    
    if exist(params.output.dir, 'dir') ~= 7
        mkdir(params.output.dir);
    end
    save('-v6', densePV_matname, 'ImgList');
    
else
    load(densePV_matname, 'ImgList');
end
ImgList_densePV = ImgList;

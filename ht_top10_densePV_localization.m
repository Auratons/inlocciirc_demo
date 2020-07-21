%Note: It first synthesize query views according to top10 pose candedates
%and compute similarity between original query and synthesized views. Pose
%candidates are then re-scored by the similarity.

PV_topN = 10; % assuming this is not larger than mCombinations
densePV_matname = fullfile(params.output.dir, 'densePV_top10_shortlist.mat');
if exist(densePV_matname, 'file') ~= 2
    
    %synthesis list
    qlist = cell(1, PV_topN*length(ImgList_densePE));
    dblist = cell(1, PV_topN*length(ImgList_densePE));
    PsList = cell(1, PV_topN*length(ImgList_densePE));
    dbind = cell(1, PV_topN*length(ImgList_densePE));
    for ii = 1:1:length(ImgList_densePE)
        for jj = 1:1:PV_topN
            idx = PV_topN*(ii-1)+jj;
            qlist{idx} = ImgList_densePE(ii).queryname;
            dblist{idx} = ImgList_densePE(ii).topNname(:,jj);
            PsList{idx} = ImgList_densePE(ii).Ps{jj};
            dbind{idx} = jj;
        end
    end
    %find unique scans
    dbscanlist = cell(size(dblist));
    dbscantranslist = cell(size(dblist));
    for ii = 1:1:length(dblist)
        for j=1:size(dblist{ii},1)
            dbpath = dblist{ii}{j};
            this_floorid = strsplit(dbpath, '/');this_floorid = this_floorid{1};
            info = parse_WUSTL_cutoutname(dbpath);
            dbscanlist{ii} = strcat(this_floorid, params.dataset.db.scan.matformat);
            dbscantranslist{ii}{j} = fullfile(this_floorid, 'transformations', ['trans_', info.scan_id, '.txt']);
        end
    end
    [dbscanlist_uniq, sort_idx, uniq_idx] = unique(dbscanlist);
    dbscantranslist_uniq = cell(size(dbscanlist_uniq));
    qlist_uniq = cell(size(dbscanlist_uniq));
    dblist_uniq = cell(size(dbscanlist_uniq));
    PsList_uniq = cell(size(dbscanlist_uniq));
    dbind_uniq = cell(size(dbscanlist_uniq));
    for ii = 1:1:length(dbscanlist_uniq)
        idx = uniq_idx == ii;
        dbscantranslist_uniq{ii} = dbscantranslist(idx);
        qlist_uniq{ii} = qlist(idx);
        dblist_uniq{ii} = dblist(idx);
        PsList_uniq{ii} = PsList(idx);
        dbind_uniq{ii} = dbind(idx);
    end
    
    %compute synthesized views and similarity

    % Because projectMesh in densePV requires up to 20 GB of RAM per one instance,
    % we need to limit the number of workers
    % TODO: optimize and leverage more workers
    poolobj = gcp('nocreate');
    delete(poolobj); % terminate any previous pool
    if strcmp(environment(), 'laptop')
        nWorkers = 1;
    else
        nWorkers = 4;
    end
    c = parcluster;
    c.NumWorkers = nWorkers;
    saveProfile(c);
    p = parpool('local', nWorkers);

    for ii = 1:1:length(dbscanlist_uniq)
        this_dbscan = dbscanlist_uniq{ii};
        this_dbscantrans = dbscantranslist_uniq{ii};
        this_qlist = qlist_uniq{ii};
        this_dblist = dblist_uniq{ii};
        this_PsList = PsList_uniq{ii};
        this_dbind = dbind_uniq{ii};
        
        %compute synthesized images and similarity scores
        parfor jj = 1:1:length(this_qlist)
        %for jj = 1:1:length(this_qlist)
            actualPs = cell(1, size(this_dblist{jj},1));
            for k=1:length(actualPs)
                P = load_CIIRC_transformation(fullfile(params.dataset.db.trans.dir, this_dbscantrans{jj}{k}));
                actualPs{k} = this_PsList{jj}{k} * P;
            end
            parfor_densePV( this_qlist{jj}, this_dblist{jj}, this_dbind{jj}, actualPs, params );
            fprintf('densePV: %d / %d done. \n', jj, length(this_qlist));
        end
        fprintf('densePV: scan %s (%d / %d) done. \n', this_dbscan, ii, length(dbscanlist_uniq));
    end
    
    %load similarity score and reranking
    ImgList = struct('queryname', {}, 'topNname', {}, 'topNscore', {}, 'Ps', {});
    for ii = 1:1:length(ImgList_densePE)
        ImgList(ii).queryname = ImgList_densePE(ii).queryname;
        ImgList(ii).topNname = ImgList_densePE(ii).topNname(1:PV_topN);
        ImgList(ii).topNscore = zeros(1, PV_topN);
        ImgList(ii).Ps = ImgList_densePE(ii).Ps(1:PV_topN);
        for jj = 1:1:PV_topN
            dbnamesId = jj;
            load(fullfile(params.output.synth.dir, ImgList(ii).queryname, sprintf('%d%s', dbnamesId, params.output.synth.matformat)), 'scores');
            cumulativeScore = sum(cell2mat(scores));
            ImgList(ii).topNscore(jj) = cumulativeScore;
        end
        
        %reranking
        [sorted_score, idx] = sort(ImgList(ii).topNscore, 'descend');
        ImgList(ii).topNname = ImgList(ii).topNname(idx);
        ImgList(ii).topNscore = ImgList(ii).topNscore(idx);
        ImgList(ii).Ps = ImgList(ii).Ps(idx);
    end
    
    if exist(params.output.dir, 'dir') ~= 7
        mkdir(params.output.dir);
    end
    save('-v6', densePV_matname, 'ImgList');
    
else
    load(densePV_matname, 'ImgList');
end
ImgList_densePV = ImgList;

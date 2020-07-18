%Note: It loads localization score and output top100 database list for each query. 

%% Load query and database list
load(params.input.qlist.path);
load(params.input.dblist.path);

%% top100 retrieval
shortlist_topN = 100;
pnp_topN = 10;
top100_matname = fullfile(params.output.dir, 'original_top100_shortlist.mat');
if exist(top100_matname, 'file') ~= 2
    ImgList = struct('queryname', {}, 'topNname', {}, 'topNscore', {});
    
    %Load score
    load(params.input.scores.path, 'score');
    
    %shortlist format
    for i=1:size(query_imgnames_all,2)
        queryName = query_imgnames_all{i};
        ImgList(i).queryname = queryName;
        ii = find(strcmp({score.queryname}, queryName));
        [~, score_idx] = sort(score(ii).scores, 'descend');
        ImgList(i).topNname = cutout_imgnames_all(score_idx(1:shortlist_topN));
        ImgList(i).topNscore = score(ii).scores(score_idx(1:shortlist_topN));
    end
    
    if exist(params.output.dir, 'dir') ~= 7
        mkdir(params.output.dir);
    end
    save('-v6', top100_matname, 'ImgList');
else
    load(top100_matname, 'ImgList');
end
ImgList_original = ImgList;

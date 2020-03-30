startup;
[ params ] = setup_project_ht_WUSTL;


%1. retrieval
ht_retrieval;

%2. geometric verification
ht_top100_densePE_localization;

%3. pose verification
ht_top10_densePV_localization;

%4. evaluate
evaluate;

if ~strcmp(environment(), "laptop")
    exit(0); % avoid "MATLAB: management.cpp:671: find: Assertion `' failed."
end
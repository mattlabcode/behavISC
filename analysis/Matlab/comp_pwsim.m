% STEP 2b: Pairwise correlation 
% INPUT:
% OUTPUT: 

clc
clear all;
cd('C:\Users\mbain25\Documents\efflisb\Analysis')
addpath(genpath('matlab'))
rng_cfig = rng(257);


%(*)define file params
story      = 1;
n_sets     = 50;
smoothtype = 'lowess';  
wsize      = 5;
online     = 0;
split_by   = 'eengagement'; % dim to split by ('engagement' / 'eengagement'/'mental_sim'/'enjoyment'/'absorption'/'attention'           

%(*)define params
do_split   = 0;             % run analysis on median split data (1 -> by split; 0 -> no split)
split      = 1;             % which split to run analysis on
n_perm     = 500;           % iterations of ts perm
do_save    = 0;

% input
if online
    if do_split
        sets = load(['sets/sets_story_' num2str(story) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '_splitby_' split_by '_split_' num2str(split) '_online' '.mat']);
    else
        sets = load(['sets/sets_story_' num2str(story) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '_online' '.mat']);
    end
else
    if do_split
        sets = load(['sets/sets_story_' num2str(story) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '_splitby_' split_by '_split_' num2str(split) '.mat']);
    else
        sets = load(['sets/sets_story_' num2str(story) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '.mat']);
    end
end
sets = sets.sets;

% savefile
if online
    if do_split
        savefile = ['splits/splits_story_' num2str(story) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '_splitby_' split_by '_split_' num2str(split) '_online' '.mat'];
    else
        savefile = ['splits/splits_story_' num2str(story) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '_online' '.mat'];
    end
else
    if do_split
        savefile = ['splits/splits_story_' num2str(story) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '_splitby_' split_by '_split_' num2str(split) '.mat'];
    else
        savefile = ['splits/splits_story_' num2str(story) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '.mat'];
    end
end
    
% # SSs per set
n_ss = size(sets{1}, 1);

% initialize arrays to store correlations and distance metrics
[robs, z, prop] = deal(zeros(n_sets, 1));
rperm           = zeros(n_sets, n_perm);

% inialize output data struct
output = [];

% COMPUTE CORRELATIONS
for s = 1:n_sets
    % compute pw correlations for current set
    corr_tmp = sets{s}(triu(corr(sets{s}'), 1) ~= 0);
    
    % Fisher-transform correlations and get median 
    robs(s) = median(atanh(corr_tmp));
end

% COMPUTE SUMMARY VALS/STATS
robs_av = tanh(median(robs))

% STORE CORR AND STATS IN OUTPUT STRUCT -----------------------------------
output.robs    = robs;
output.robs_av = robs_av;

% SAVE --------------------------------------------------------------------
if do_save
    if exist(savefile, 'file')
        error(['file ' savefile ' exists already. Either delete or use new parameters!'])
    end
    
    save(savefile, '-struct', 'output')
end
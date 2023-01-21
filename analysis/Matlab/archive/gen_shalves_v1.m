% STEP 2: Split-Half Group Generation
% generate all distinct combinations of 2 groups of n SSs, where n = # SSs per set/2. 
% Collapse each group into single time course to give two time courses per set.
% INPUT: 
%   n_sets*1 cell array
% OUTPUT: 

clc
clear all;
cd('/Users/matthewbain/Documents/Science/Experiments/efflisisc/analysis/')
addpath(genpath('matlab'))
rng_cfig = rng(257);

% define file params
online     = 1;
story      = 1:2;
n_sets     = 50;
smoothtype = 'lowess';
%wsize      = [62 63];
wsize      = 10;
do_split   = 0;

if do_split
    split_by   = 'absorption';               
    split_trun = 0;         
    split      = 1;           
else
    split_by   = '';             
    split_trun = 0;
    split      = 0;           
end

trunc = 'no';

% define params
randsamp_thresh = 1000; % if # splits > thresh, randsamp subset of 'randsamp_thresh' combos

do_save = 0;
overwrite = 0;

% -------------------------------------------------------------------------
% loop through OL/IL; each story/wsize specified
for o = 1:length(online)
for st = 1:length(story)
for w = 1:length(wsize)
    
% input
datadir = 'sets';
if length(wsize) > 1
    datadir = 'sets/smoothing';
end
sets = load([datadir '/' 'sets' '_story_' num2str(story(st)) '_nsets_' num2str(n_sets) ...
    '_smoothtype_' smoothtype '_wsize_' num2str(wsize(w)) '_dosplit_' num2str(do_split) ...
    '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) ...
    '_trunc_' trunc '_online_' num2str(online(o)) '.mat']);
sets = sets.sets;

% savefile
outdir = 'splits';
if length(wsize) > 1
    outdir = 'splits/smoothing';
end
savefile = [outdir '/' 'splits' '_story_' num2str(story(st)) '_nsets_' num2str(n_sets) ...
    '_smoothtype_' smoothtype '_wsize_' num2str(wsize(w)) '_dosplit_' num2str(do_split) ...
    '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) ...
    '_trunc_' trunc '_online_' num2str(online(o)) '.mat'];
    
if do_save
    if exist(savefile, 'file')
        if overwrite 
            disp('overwriting existing savefile')
        else
            %continue
            error(['file ' savefile ' exists already. Either delete or use new parameters!'])
        end
    end
end
% -------------------------------------------------------------------------

% --- # SS per split ---
% # SS (adjust so divisible by 2)
n_ss = size(sets{1}, 1);
if mod(n_ss, 2) > 0
    n_ss = n_ss - 1;
end
splits_n = n_ss/2; % # SS per split

% all distinct combinations of 2 groups of SSs
bin    = 1:n_ss;
combos = combnk(bin, splits_n);

% randomly sample subset of combos if # combos > randsamp_thresh
if length(combos) > randsamp_thresh
    combos = combos(datasample(1:length(combos), randsamp_thresh, 'replace', false), :);
    %fprintf('using %d combos\n', randsamp_thresh)
end    

% eliminate half of combos that are mirror reflections (for efficiency - redundancy does not affect outcome of correlation) 
%***how does this not get messed up by taking random subset of combos (above section)??
for c = 1:length(combos)/2
    combos(ismember(combos, find(~ismember(bin, combos(c, :))), 'rows'), :) = [];
end

%fprintf('%d splits per set\n', length(combos))

% array to store split-half ts x 2 for each set*combo
splits = cell(length(sets), length(combos));

for s = 1:length(sets)    
    for c = 1:length(combos)
        % get combos for each group
        half1_combo = combos(c, :);
        half2_combo = bin(~ismember(bin, combos(c, :)));

        % mean ts for each group
        half1 = nanmean(sets{s}(half1_combo, :));
        half2 = nanmean(sets{s}(half2_combo, :));
        splits(s, c) = {[half1; half2]}; % store for current combo/set
    end
end

% save --------------------------------------------------------------------
if do_save
    save(savefile, '-mat', 'splits', '-v7.3')
end

end
end
end
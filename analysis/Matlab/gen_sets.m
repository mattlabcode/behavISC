% STEP 1: Set Generation
% resample N time courses, pooled into v versions, to generate s sets w/ 
% n supersubjects (SSs) each. 
% For N/v iterations, randomly sample (w/out replacement) v subjects 
% (1 from each version pool) and collapse ts to form set of N/n SSs.
% INPUT:  
%   N/v*v cell array containing time courses to correlate
% OUTPUT: 
%   s*1 cell array, where each cell contains set w/ N/v ts

clc
clear;
cd('/Users/matthewbain/Documents/Science/Experiments/efflisisc/analysis/')
addpath(genpath('matlab'))
rng_cfig = rng(257);

% data params
online = 1; 
do_split = 0; 

if do_split
    split_by   = 'absorption';               
    split_trun = 0;    
    split      = 1;              % which split to run analysis on
else
    split_by   = '';             
    split_trun = 0;
    split = 0;
end

% local params
scrambled  = 0:1;             % 1 for scrambled; 0 for not
story      = 1:2;
n_sets     = 50;   
smoothtype = 'lowess';      % 'movmean'/'median' -> rect window movavg; 'gaussian' -> gaussian window; 'lowess' -> local linear regression w/ 1st deg polynormial   
%wsize      = 1:1:100;            % size of smoothing window (in time points; 0 -> no smoothing)
wsize      = 10;

n_prtrials = 0;             % # trials to reject from start of time courses (practice trials) ***remove?
trunc      = 'no';          % 'no' -> use whole ts; 'nsig'/'sig' -> truncate SS ts to only n/sig sect. (concat); 

do_save    = 1; 
overwrite  = 1;             % overwrite any existing savefiles (NOT RECOMMENDED)

% -------------------------------------------------------------------------
% loop through OL/IL; each story/wsize specified
for sc = scrambled
for o = 1:length(online)
for st = 1:length(story)
for w = 1:length(wsize)

% input
ts = load(['timeseries' '/' 'ts' '_scrambled_' num2str(sc) '_dosplit_' num2str(do_split) '_splitby_' split_by '_splittrun_' num2str(split_trun) '_online_'  num2str(online(o)) '.mat']);

if do_split % load in data
    split_name = ['split' num2str(split)];
    tspool     = ts.ypool(story(st)).(split_name);
else
    ts     = ts.ts;
    tspool = ts(story(st)).ypool;
end

if ~strcmp(trunc, 'no')
    datadir = 'output_varana';
    if length(wsize) > 1
        datadir = 'output_varana/smoothing';
    end
    varana = load([datadir '/' 'stats_' 'varanalysis' '_story_' num2str(story(st)) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize(w)) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) '_online_' num2str(online(o)) '.mat']); % varana output
    tcsig = varana.tcsig_pav < .05;
end

% savefile
outdir = 'sets';
if length(wsize) > 1
    outdir = 'sets/smoothing';
end
savefile = [outdir '/' 'sets' '_story_' num2str(story(st)) '_scrambled_' num2str(sc) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize(w)) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) '_trunc_' trunc '_online_' num2str(online(o)) '.mat'];

if do_save
    if exist(savefile, 'file')
        %{
        if length(wsize) > 1
            fprintf('skipping wsize %d because savefile already written.\n', wsize)
            continue
        else
        %}
        if overwrite 
            disp('overwriting existing savefile')
        else
            error(['file ' savefile ' exists already. Either delete or use new parameters!'])
        end
        %end
    end
end

% -------------------------------------------------------------------------
% --- get dimensions ---
% number of versions
n_ver = size(tspool, 2);

% number of SSs (non-subj - empty - row = all nan)
tspool_empty = cellfun(@all, cellfun(@isnan, tspool, 'UniformOutput', false)); 
n_SS  = size(tspool, 1) - max(sum(tspool_empty)); % base on highest # of subjs in each v 

%fprintf('%d supersubjects per set\n', n_SS)

% generate sets -----------------------------------------------------------
sets = cell(n_sets, 1); % array to store sets
for i = 1:n_sets
    set = zeros(n_SS, length(tspool{1})); % array to store SS for current set 
    ssubj_ix = zeros(n_SS, n_ver);        % array to store new configuration of SSs
    
	% get random selection of v S from each vpool
    for j = 1:n_ver
        ssubj_ix(:, j) = randperm(n_SS); 
    end

    % --- generate current set --- 
    for ss = 1:n_SS
        % draw S from each vpool, create SS & store in set
        ss_curr = zeros(n_ver, length(tspool{1})); % array to store curr SS
        for v = 1:n_ver
            ss_curr(v, :) = tspool{ssubj_ix(ss, v), v}; 
        end
        set(ss, :) = nanmean(ss_curr); % store created SS in current set

        % --- smooth SS ---
        if wsize(w)
            % interpolate over any NaN values using linear regression
            set_tmp    = set(ss, :);
            t          = 1:length(set_tmp);
            set(ss, :) = interp1(t(~isnan(set_tmp)), set_tmp(~isnan(set_tmp)), t);

            % smooth
            set(ss, :) = smoothdata(set(ss, :), smoothtype, wsize(w));
        end
    end
    
    % reject practice trials from set***
	set = set(:, (n_prtrials + 1):end);
    
    % --- remove nonsig portions of ts (based on varana op) & concat. ---
    switch trunc
        case 'sig'
            set = set(:, tcsig);
        case 'nsig'
            set = set(:, ~tcsig);
    end
    
    
    % store current set
    sets(i, :) = {set}; 
end

% save --------------------------------------------------------------------
if do_save
    save(savefile, 'sets')
end
    
end
end
end
end
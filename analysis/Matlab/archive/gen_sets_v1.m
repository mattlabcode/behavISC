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
story      = 1:2;
n_sets     = 50;   
smoothtype = 'lowess';      % 'movmean'/'median' -> rect window movavg; 'gaussian' -> gaussian window; 'lowess' -> local linear regression w/ 1st deg polynormial   
%wsize      = 1:1:100;            % size of smoothing window (in time points; 0 -> no smoothing)
wsize      = 10;

n_prtrials = 0;             % # trials to reject from start of time courses (practice trials) ***remove?
trunc      = 'no';          % 'no' -> use whole ts; 'nsig'/'sig' -> truncate SS ts to only n/sig sect. (concat); 

do_save    = 0; 
overwrite  = 0;             % overwrite any existing savefiles (NOT RECOMMENDED)

% -------------------------------------------------------------------------
% loop through OL/IL; each story/wsize specified
for o = 1:length(online)
for st = 1:length(story)
for w = 1:length(wsize)

% input
ts = load(['timeseries' '/' 'ts' '_dosplit_' num2str(do_split) '_splitby_' split_by '_splittrun_' num2str(split_trun) '_online_'  num2str(online(o)) '.mat']);

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
savefile = [outdir '/' 'sets' '_story_' num2str(story(st)) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize(w)) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) '_trunc_' trunc '_online_' num2str(online(o)) '.mat'];

if do_save
    if exist(savefile, 'file')
        if length(wsize) > 1
            fprintf('skipping wsize %d because savefile already written.\n', wsize)
            continue
        else
            if overwrite 
                disp('overwriting existing savefile')
            else
                error(['file ' savefile ' exists already. Either delete or use new parameters!'])
            end
        end
    end
end

% -------------------------------------------------------------------------
% --- get dimensions ---
% number of versions
n_ver = size(tspool, 2);

% number of SSs (non-subj - empty - row = all nan)
tspool_empty = cellfun(@all, cellfun(@isnan, tspool, 'UniformOutput', false)); 
n_SS  = length(tspool) - max(sum(tspool_empty)); % base on highest # of subjs in each v 

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

%{
%% ***TEMP - compute mean r over all Ssets using split-half corr
f = waitbar(0, 'Please wait...');

% compute r using split-half method on all sets
for i = 1:length(sets)
    waitbar(i/length(sets), f, ['set ' num2str(i) ' / ' num2str(length(sets))]);
    [robs_all, p_all] = corr_ts(sets{i});
    robs_av(i) = mean(robs_all);
end

delete(f)

mean(robs_av)

%% ***TEMP - collapse all sets into 1 mean time course and plot
% define array to store mean ts for each set
ts_mean = [];

% average ts in each set
for s = 1:length(sets)
    ts_mean(s, :) = mean(sets{s}, 1);
end

% get mean ts over all sets (grand mean)
ts_gmean = mean(ts_mean, 1);

% plot
plot(ts_gmean)

%% ***TEMP - collapse all sets into 10 mean SSs and plot
n_SS = size(sets{1}, 1);

% initialize array to store sets by SS instead of set
sets_SS    = cell(n_SS, 1);
sets_SS(:) = {NaN(length(sets), length(sets{1}))};

% sets_SS row counter
i = 1;

for s = 1:length(sets)
    for ss = 1:n_SS
        sets_SS{ss}(i, :) = sets{s}(ss, :);
    end
	i = i + 1;
end

% get mean/var for each SS aggregate 
for ss = 1:n_SS
    sets_mean(ss, :) = mean(sets_SS{ss}, 1); 
    sets_var(ss, :)  = var(sets_SS{ss}, 1); 
end

plot(sets_mean')
%% ***TEMP - verify that all sets have same mean time course
sets_av = zeros(length(sets), length(sets{1}));

for i = 1:length(sets)
    sets_av(i, :) = nanmean(sets{i});
end

figure
for i = 1:size(sets_av, 1), hold on, plot(sets_av(i, :)); end

nnz(abs(diff(sets_av)) > .01)/numel(sets_av);

%% ***TEMP - for given Sset, plot all SS ts overlaid
set = sets{1};

% smooth more for visualization
for i = 1:size(set, 1)
    set(i, :) = smooth(set(i, :), 10);
end

figure
for i = 1:size(set, 1)
    plot(set(i, :), 'linewidth', 2)
    hold on
end

%% ***TEMP - for given Sset, plot mean ts w/ shaded error bars
plot_type = 'tsmean';
savefile = ['figures_ts' '/' plot_type '_story_' num2str(story) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '.eps'];

% Sset to use
set = smoothts(sets{1}, 'b', 5);

figure

% data to plot
x   = ts(story).x/60;
y   = mean(set); 
std = nanstd(set);
n   = sum(~isnan(set));
sem = std./sqrt(n);               

% ts w/ shaded error bars
fig = shadedErrorBar(x, y, sem); 
    
% axis params
%set(gca, 'ticklength', [0.001 0.001]);
%set(gca, 'xlim', [0 13.7]);
%set(gca, 'ylim', [-.3 .3]);

% labels
title(['Mean RT Time Course for Dual-Task: Story ' num2str(story)]);
xlabel('Story Time (min)');
ylabel('Mean RT (sec; norm)');

box off

print(gcf, savefile, '-deps', '-r300')

%% ***TEMP - for rand generated set and given smooth params get mean r obs 
% params
n_sets     = 30;
smoothtype = 'rect'; % 'rect' -> rectangular window movavg; 'gauss' -> gaussian window; 'lowess' -> local linear regression w/ 1st degree polynormial   
wsize       = 2;        % size of smoothing window (in time points; 0 -> no smoothing)

% initialize array to store observed correlations
robs_av = zeros(1, n_sets);

f = waitbar(0,'Please wait...');
for i = 1:n_sets
	waitbar(i/n_sets, f, ['Superset ' num2str(i) ' / ' num2str(n_sets)]);
    
    % rand set
    set = rand(min(sum(~cellfun(@isempty, tspool))), length(ts(story).ystacked))*.4 - .2;

    % SMOOTH
    for ss = 1:n_SS
        if wsize
            % interpolate over any NaN values using linear regression
            set_temp   = set(ss, :);
            t             = 1:length(set_temp);
            set(ss, :) = interp1(t(~isnan(set_temp)), set_temp(~isnan(set_temp)), t);

            switch smoothtype
                case 'rect'
                    % smooth using rectangular function-weighted moving average over wsize
                    set(ss, :) = smoothts(set(ss, :), 'b', wsize);
                case 'gauss'
                    set(ss, :) = smoothts(set(ss, :), 'g', wsize, wsize/3);
                case 'lowess'
                    set(ss, :) = smooth(set(ss, :), wsize, 'lowess');
            end
        end
    end

    % reject practice trials from set
	set = set(:, (n_prtrials + 1):end);
    
    % CORRELATE
    % compute r using split-half method on all Ssets and store
    [robs_all, p_all] = corr_ts(set);
    robs_av(i) = mean(robs_all);
end

delete(f)
    
r_obs = mean(mean((robs_av)))
%}
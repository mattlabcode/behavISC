% STEP 3: 
% Calc obs and perm corr dist on split-half ts for each set. 
% Compute distance metrics between obs and perm dists.

clc
clear all;
cd('/Users/matthewbain/Documents/Science/Experiments/efflisisc/analysis/')
addpath(genpath('matlab'))
rng_cfig = rng(257);

% define file params
online     = 1;
scrambled  = 0:1;
story      = 1:2;
n_sets     = 50;
smoothtype = 'lowess';  
%wsize      = 1:1:100;
wsize      = 10;
do_split   = 0;             

if do_split
    split_by   = 'absorption';               
    split_trun = 0;            
    split_num  = 1;            % ***replace var used later 'split' and change this back to 'split'
else
    split_by   = '';             
    split_trun = 0;
    split_num  = 0;           
end

trunc = 'no';

% define params
n_perm     = 500;       % it of ts perm (default = 500; 1 if running smoothing analysis for many window sizes); if 'nt', use # tps
nc         = 1000;          % # criteria (iterations) for computing ROC
nbins_null = 100;           % # bins for histograms

do_save    = 1;
overwrite  = 1;

% -------------------------------------------------------------------------
% loop through OL/IL; each story/wsize specified
for sc = scrambled
for o = 1:length(online)
for st = 1:length(story)
for w = 1:length(wsize)

% input
datadir = 'splits';
if length(wsize) > 1 
    datadir = 'splits/smoothing';
end
splits = load([datadir '/' 'splits' '_story_' num2str(story(st)) '_scrambled_' num2str(sc) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize(w)) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split_num) '_splittrun_' num2str(split_trun) '_trunc_' trunc '_online_' num2str(online(o)) '.mat']);
splits = splits.splits;

% savefile
outdir = 'correlations';
if length(wsize) > 1
    outdir = 'correlations/smoothing';
end
savefile = [outdir '/' 'corr' '_story_' num2str(story(st)) '_scrambled_' num2str(sc) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize(w)) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split_num) '_splittrun_' num2str(split_trun) '_trunc_' trunc '_nperm_' num2str(n_perm) '_online_' num2str(online(o)) '.mat'];

% check if savefile already exists
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
            if length(wsize) > 1
                continue
            else
                error(['file ' savefile ' exists already. Either delete or use new parameters!'])
            end
        end
        %end
    end
end

% -------------------------------------------------------------------------
% compute correlation distributions ---------------------------------------
% get split dimensions
n_splits = size(splits, 2);
n_t      = length(splits{1});

% get # of perm if using # tps (all possible cshifts, +/-ve, incl. 0)
if strcmp(n_perm, 'nt')
    n_perm = n_t;
end

% initialize arrays to store correlations and distance metrics
[robs, z, prop, p_ntest]        = deal(zeros(n_sets, n_splits));
rperm                           = zeros(n_splits, n_perm, n_sets); %***weird...make dim s,h,p (not h,p,s)
[auc, dprime]                   = deal(zeros(n_sets, 1));
[tpr, fpr]                      = deal(zeros(n_sets, nc));
[bincounts_obs, bincounts_perm] = deal([]);

output = []; % output data struct

%f = waitbar(0,'Please wait...');
tic
% go through all sets
for s = 1:n_sets
	%waitbar(s/size(splits, 1), f, ['set ' num2str(s) ' / ' num2str(size(splits, 1))]);
    
    % go through all splits
    for h = 1:size(splits, 2)    
        
        % ---obs dist--- 
        % get current split
        split = splits{s, h};
           
        % get observed corr
        r          = corr(split');
        robs(s, h) = atanh(r(2)); 
        
        % ---perm dist---
        for p = 1:n_perm
            
            % shift one ts halve (note: ARB which one)
            %shift = p;
            if length(wsize) > 1
                shift = randsample(1:n_t, 1); 
            else
            	%shift = randsample(round(wsize/2):n_t, 1); % wsize to avoid perm w/ smeared ts
            	shift = randsample(wsize:n_t, 1); % wsize to avoid perm w/ smeared ts
            end
            
            surr = circshift(split(1, :), shift);

            % get corresponding corr and store
            r              = corr([surr; split(2, :)]');
            rperm(h, p, s) = atanh(r(2));
        end           

        % test normality of perm dist (p = .05) (p_av = median(median(p_ntest, 2)))
        rperm_std = (rperm(h, :, s) - nanmean(rperm(h, :, s)))/nanstd(rperm(h, :, s)); % standardize pdist to test against std. N
        if ~isnan(rperm_std)
            [~, p_ntest(s, h)] = kstest(rperm_std);
        else
            p_ntest(s, h) = nan;
        end
            
        % get distance metric btwn robs and pdist for current h
        z(s, h)    = robs(s, h) - mean(rperm(h, :, s))/std(rperm(h, :, s));
        prop(s, h) = nnz(rperm(h, :, s) > robs(s, h))/n_perm;  
        
        % plot (s, h) robs value against rperm dist (as movie)
        %{
        hist(rperm(h, :, s)), xlim([-.4 .4]), ylim([0 150]) 
        line([robs(s, h) robs(s, h)], [0 100], 'color', [1 0 0], 'linestyle', '--')
        text(-.35, 130, {['z = ' num2str(z(s, h))], ['ntest p = ' num2str(p_ntest(s, h))] ['prop = ' num2str(prop(s, h))]}, 'FontSize', 9) 
        pause(5)
        %}
    end
    
    % collapse h pdists into one
    pdist = reshape(rperm(:, :, s), [1 n_splits*n_perm]);
    odist = robs(s, :);
    
    % compare pdist to rdist for current s & store distance metrics
    [auc(s), tpr(s, :), fpr(s, :), criteria] = ROC(odist, pdist, nc, 0);
    dprime(s)                                = sqrt(2)*norminv(auc(s));
    
    % get dist info for plotting ------------------------------------------
    if s == 1
        % set nbins for alt dist
        nbins_alt = round((nbins_null/range(pdist))*range(odist));
        
        % get normalized bin counts and edges for dists
        [bincounts_perm(s, :), binedges_perm] = histcounts(pdist, nbins_null, 'Normalization', 'probability');
        [bincounts_obs(s, :), binedges_obs]   = histcounts(odist, nbins_alt, 'Normalization', 'probability');
    else
        % get bin counts using edges from 1st set
        bincounts_perm(s, :) = histcounts(pdist, binedges_perm, 'Normalization', 'probability');
        bincounts_obs(s, :)  = histcounts(odist, binedges_obs, 'Normalization', 'probability');
    end
end
%delete(f)
%fprintf('correlations took %d seconds to run.\n', toc);
%fprintf('pdist normality test p_av = %f\n', median(median(p_ntest, 2)));

% ---compute summary vals/stats---
% get summary values
robs_av   = tanh(median(median(robs, 2))); % (s,h): avg across splits then sets
robs_var  = tanh(median(var(robs, 0, 2))); %***0 weird...
robs_ran  = [tanh(median(min(robs, [], 2))) tanh(median(max(robs, [], 2)))];
robs_n    = n_splits;

rperm_av  = tanh(median(median(median(rperm, 2), 1))); % (h,p,s): avg across perms, then splits, then sets  
rperm_var = tanh(median(median(var(rperm, 0, 2))));
rperm_ran = [tanh(median(median(min(rperm, [], 2)))) tanh(median(median(max(rperm, [], 2))))];
rperm_n   = n_perm; 

% convert z scores to p values (2-tailed) and average
z_av  = 10^mean(mean(log10(abs(z)), 2));
z_p   = 2*(1 - normcdf(abs(z)));
zp_av  = 10^mean(mean(log10(z_p), 2)); % ***not confident about this...
zp_av2 = 1 - normcdf(z_av);

% compute mean prop (p value), first log transforming all non-zero prop
%prop_av = 10^(sum(log10(prop(prop ~= 0))) / (n_splits*n_sets));
prop_av = median(median(prop, 2), 1); %***only correct option?...

% compute t-test between mean r for perm/obs dist ***no good. Assumes ind.
%[t, p, df, ci] = ttest_onmeans(robs_av, rperm_av, robs_var, rperm_var, robs_n, rperm_n, .05);

% store corr and stats in output struct -----------------------------------
output.robs           = robs;
output.robs_av        = robs_av;
output.robs_var       = robs_var;
output.robs_ran       = robs_ran;
output.robs_n         = robs_n;

output.pav_ntest      = median(median(p_ntest, 2));

output.rperm          = rperm;
output.rperm_av       = rperm_av;
output.rperm_var      = rperm_var;
output.rperm_ran      = rperm_ran;
output.rperm_n        = rperm_n;

output.z              = z;
output.z_av           = z_av;
output.zp_av          = zp_av;
output.prop_av        = prop_av;

output.auc            = auc;
output.auc_av         = mean(auc);
output.auc_var        = var(auc);
output.auc_ran        = [min(auc) max(auc)];
output.dprime         = dprime;
output.dprime_av      = mean(dprime);
output.dprime_var     = var(dprime);
output.dprime_ran     = [min(dprime) max(dprime)];

output.tpr            = tpr;
output.fpr            = fpr;
output.bincounts_obs  = bincounts_obs;
output.bincounts_perm = bincounts_perm;
output.binedges_obs   = binedges_obs;
output.binedges_perm  = binedges_perm;

output.rng_cfig       = rng_cfig;

% save --------------------------------------------------------------------
if do_save
    save(savefile, '-struct', 'output')
end

end
end
end
end
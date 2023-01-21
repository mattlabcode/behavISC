% STEP 3: 
% INPUT:
% OUTPUT: 

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
wsize      = 10;
do_split   = 0;      

if do_split
    split_by   = 'enjoyment';               
    split_trun = 0;         
    split      = 1;           
else
    split_by   = '';             
    split_trun = 0;
    split      = 0;           
end

% define params
alpha      = .05;
t_stamp    = 1;       % 1 -> show story times beneath each sig bar; 0 -> no times
ix_thresh  = 9;       % threshold size of bar to add num identifier above (in sec; 0 -> skip)
n_perm     = 500;     % iterations of ts perm (default = 500) ***make this param conditional on null_type == prop/z/tind
error_type = 'ci';    % 'sem' -> plot error bars as sem; 'ci' -> plot as ci
null_type  = 'srank'; % test for comparing obs ts to null: 'tzero'; 'srank'; 'prop'; 'z'; 'tind' 

bar_thresh = 1;

do_save    = 1;
overwrite  = 1;

% -------------------------------------------------------------------------
% loop through OL/IL; each story/wsize specified
for sc = scrambled
for o = 1:length(online)
for st = 1:length(story)
for w = 1:length(wsize)

%(*)define plotting params
vert_offset = [-.75 -.65]; % height of sig bars/ix
vert_offset_tstamp = [-.57 -.6]; % height of tstamps (alternate to save space)

% input
datadir = 'sets';
if length(wsize) > 1
    datadir = 'sets/smoothing';
end
sets = load([datadir '/' 'sets' '_story_' num2str(story(st)) '_scrambled_' num2str(sc) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize(w)) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) '_trunc_no' '_online_' num2str(online(o)) '.mat']);
sets = sets.sets;

ts = load(['timeseries' '/' 'ts' '_scrambled_' num2str(sc) '_dosplit_' num2str(do_split) '_splitby_' split_by '_splittrun_' num2str(split_trun) '_online_'  num2str(online(o)) '.mat']);
ts   = ts.ts;

% savefile
savefile_fig = ['figures_ts' '/' 'ts_varanalysis' '_story_' num2str(story(st)) '_scrambled_' num2str(sc) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize(w)) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) '_online_' num2str(online(o)) '.svg'];

outdir = 'output_varana';
if length(wsize) > 1
    outdir = 'output_varana/smoothing';
end
savefile_stats = [outdir '/' 'stats_varanalysis' '_story_' num2str(story(st)) '_scrambled_' num2str(sc) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize(w)) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) '_online_' num2str(online(o)) '.mat'];

% check if savefile already exists
if do_save
    if exist(savefile_fig, 'file') || exist(savefile_stats, 'file')
        if length(wsize) > 1
            fprintf('skipping wsize %d because savefile already written.\n', wsize)
            continue
        else
            if overwrite 
                disp('overwriting existing savefile')
            else
                error(['savefile exists already. Either delete or use new parameters!'])
            end
        end
    end
end

% -------------------------------------------------------------------------
% prepare comparison btwn obs & null ts -----------------------------------
% get set dimensions
n_ss = size(sets{1}, 1);
n_t  = size(sets{1}, 2);

% initialize arrays to store obs/perm ts and stats
tsobs                               = sets;
tsperm                              = zeros(n_ss, n_t, n_perm, n_sets); %***check these dimensions
[tsobs_av, tsobs_var, tsobs_sem]    = deal(zeros(n_sets, n_t));
[tsperm_av, tsperm_var, tsperm_sem] = deal(zeros(n_sets, n_t));

[tcsig_p, tcsig_stat]               = deal(zeros(n_sets, n_t));
tcsig_ci                            = zeros(2, n_t, n_sets);

output                              = [];                   % output struct

for s = 1:n_sets
    % ---obs ts---
    % calc stats on tsobs
    tsobs_av(s, :)  = nanmean(tsobs{s});
    tsobs_var(s, :) = nanvar(tsobs{s});
	tsobs_sem(s, :) = nanstd(tsobs{s})/sqrt(n_ss);

    % ----perm ts---    
    if any(strcmp(null_type, {'prop', 'z', 'tind'})) % only compute perm for null_types that require it 
        % generate n_perm permuted sets of SS; get av mean/var of each
        for p = 1:n_perm

            % get vector of shifts and apply to obs ts
            shifts = randi(n_t, [1 n_ss]); %***wsize:n_t??
            for c = 1:length(shifts)
                tsperm(c, :, p, s) = circshift(tsobs{s}(c, :), shifts(c));
            end
        end

        % get mean and var ts of each permuted set ***convert this to 2D (squeeze) and rotate (n_perm * t)
        tsperm_av_tmp  = mean(tsperm(:, :, :, s), 1);
        tsperm_var_tmp = var(tsperm(:, :, :, s), 1);

        % get mean tsperm mean and var for current set
        tsperm_av(s, :)  = mean(tsperm_av_tmp, 3);
        tsperm_var(s, :) = mean(tsperm_var_tmp, 3);
    end
    
    % compare obs and null ts (by tp) -------------------------------------
    switch null_type
        % --- Method 1: one-sample ttest comparing ts obs to zero ---
        case 'tzero'
            [~, p, ci, stats] = ttest(tsobs{s});
            
            tcsig_stat(s, :)  = stats.tstat;
            tcsig_p(s, :)     = p;
            tcsig_ci(:, :, s) = ci;
            
        % --- Method 2: signed-rank test of ts obs against zero ---
        case 'srank'
            %[p, stats] = deal(zeros(n_sets, n_t));
            for t = 1:n_t
                [p, ~, stats] = signrank(tsobs{s}(:, t), 0, 'tail', 'right');
                
                tcsig_stat(s, t)  = stats.signedrank;
                tcsig_p(s, t)     = p;
            end
            %tcsig_ci(:, :, s) = ci;
            
        % --- Method 3: compute prop of the n_perm av perm RT > av obs ---
        case 'prop'
            for t = 1:n_t    
                % proportion of perm > obs RT for current tp *plot to check: histogram(tsperm_av_tmp(1, t, :)) vs tsobs_av(s, t)
                tcsig_p(s, t) = (nnz(tsperm_av_tmp(1, t, :) > tsobs_av(s, t)))/n_perm;
             end

        % --- Method 4: get distance btwn av obs and n_perm av perm RT ---
        case 'z'
            for t = 1:n_t    
                % calc p from z comparing av obs RT to perm dist
                k(s, t)       = kstest(squeeze(tsperm_av_tmp(1, t, :))); % test normality 
                z             = (tsobs_av(s, t) - mean(tsperm_av_tmp(1, t, :))) / std(tsperm_av_tmp(1, t, :));
                tcsig_p(s, t) = normcdf(z); % one-tailed
            end
            
        % --- Method 5: ind ttest btwn obs and perm dist ---
        case 'tind'
            [n1, n2]          = deal(ones(1, n_t)*n_ss);
            [tval, p, ~, ci]  = ttest_onmeans(tsobs_av(s, :), tsperm_av(s, :), tsobs_var(s, :), tsperm_var(s, :), n1, n2, alpha);
            %[h, p, ci, stats]  = ttest2(tsobs_av(s, :), tsperm_av(s, :), 'vartype', 'unequal');
            
            tcsig_stat(s, :)  = tval;
            tcsig_p(s, :)     = p/2;
            tcsig_ci(:, :, s) = ci;
    end

end
  
% store ts and stats in output struct -------------------------------------
output.tsobs      = tsobs;
output.tsperm     = tsperm;
output.tsobs_av   = tsobs_av;
output.tsobs_var  = tsobs_var;
output.tsobs_sem  = tsobs_sem;
output.tsperm_av  = tsperm_av;
output.tsperm_var = tsperm_var;
output.tsperm_sem = tsperm_sem;

output.tcsig_stat = tcsig_stat;
output.tcsig_p    = tcsig_p;
output.tcsig_pav  = 10.^mean(log10(output.tcsig_p));
output.tcsig_ci   = tcsig_ci;

output.rng_cfig   = rng_cfig;

%%{
% -------------------------------------------------------------------------
% plot --------------------------------------------------------------------
% colour scheme
colours1 = [190 22 34; 7 162 228]/255; % R/B dark
colours2 = [242 209 212; 199 236 255]/255; % R/B light light

% data for plotting
x       = ts(story(st)).x;
yobs    = mean(output.tsobs_av); 
yperm   = mean(output.tsperm_av); 
varobs  = mean(output.tsobs_var); 
varperm = mean(output.tsperm_var); 
semobs  = sqrt(varobs/n_ss);
semperm = sqrt(varperm/n_ss);
ciobs   = abs(norminv(alpha/2))*semobs; % sem*z-crit for defined alpha 2-tailed
ciperm  = abs(norminv(alpha/2))*semperm;                     
p       = 10.^mean(log10(output.tcsig_p)); % average p values after log transforming   

output.ciobs = ciobs;

%x = 1:length(yobs); % adjustment for truncated ts

% --- get title tags ---
if online(o)
    online_tag = 'online';
elseif ~online(o)
    online_tag = 'lab';
end
if story(st) == 1
    story_tag = 'Arctic';
elseif story(st) == 2
    story_tag = 'Space';
end

% indices of all significant points 
if any(strcmp(null_type, {'prop', 'z', 'tind'})) % if null isn't 0
    sig_pts = find(p < alpha); 
else
    % check if confidence interval intersects 0
    sign_bounds = sign([yobs + ciobs; yobs - ciobs])';
    sig_pts = find((sign_bounds(:, 1) == sign_bounds(:, 2)));
end
    
% plot error bars
if strcmp(error_type, 'sem')
    if any(strcmp(null_type, {'prop', 'z', 'tind'})) % if perm not 0
        pl(1) = fill([x'/60; flipud(x'/60)], [yperm' - semperm'; flipud(yperm' + semperm')], colours2(2, :), 'linestyle', 'none'); hold on
    end
    pl(2) = fill([x'/60; flipud(x'/60)], [yobs' - semobs'; flipud(yobs' + semobs')], colours2(1, :), 'linestyle', 'none'); hold on
elseif strcmp(error_type, 'ci')
    %fill([x'/60; flipud(x'/60)], [yperm' - ciperm'; flipud(yperm' + ciperm')], colours2(2, :), 'linestyle', 'none'); hold on
    pl(1) = fill([x'/60; flipud(x'/60)], [yobs' - ciobs'; flipud(yobs' + ciobs')], colours2(1, :), 'linestyle', 'none'); hold on
end

% plot null ts (just 0 if not perm)
if any(strcmp(null_type, {'prop', 'z', 'tind'}))  
	f(2) = plot(x/60, yperm, 'color', colours1(2, :), 'LineWidth', 3); hold on
else
	f(2) = plot(x/60, zeros([1 length(x)]), 'color', colours1(2, :), 'LineWidth', 3); hold on
end

% plot obs ts (plot perm first so obs on top)
f(1) = plot(x/60, yobs, 'color', colours1(1, :), 'LineWidth', 3); hold on

        
% plot bars below sig sections --------------------------------------------
% initialize array to store start- and end-points of bars
line_ix = cell(1); 

% starting index of all bars
line_startix = find(diff(sig_pts) > 1) + 1;   

% ---draw first bar---
if ~isempty(line_startix) && line_startix(1) > 2

    % get starting and ending index
    line_start   = sig_pts(1);
    line_end     = sig_pts(line_startix(1) - 1);
    line_ix{end} = [ts(story(st)).x(line_start) x(line_end)];

    % draw bar
    if line_end - line_start > bar_thresh
        line([x(line_start) x(line_end)]/60, [vert_offset(1) vert_offset(1)], 'LineWidth', 4, 'Color', [0 0 0]);
    end
end

% ---draw all remaining bars---
for i = 1:length(line_startix)

    % draw last bar
    if i == length(line_startix)

        % get starting and ending index
        line_start       = sig_pts(line_startix(i));
        line_end         = sig_pts(end);
        line_ix{end + 1} = [x(line_start) x(line_end)];

        % draw bar
        if line_end - line_start > bar_thresh
            line([x(line_start) x(line_end)]/60, [vert_offset(1) vert_offset(1)], 'LineWidth', 4, 'Color', [0 0 0]);   
        end
        
    else
        % get starting and ending index
        line_start = sig_pts(line_startix(i));
        line_end   = sig_pts(line_startix(i + 1) - 1);

        if cellfun(@isempty, line_ix)
            line_ix{end}     = [x(line_start) x(line_end)];
        else
            line_ix{end + 1} = [x(line_start) x(line_end)];
        end

        % draw bar
        if line_end - line_start > bar_thresh
            line([x(line_start) x(line_end)]/60, [vert_offset(1) vert_offset(1)], 'LineWidth', 4, 'Color', [0 0 0]);
        end
    end
end

% add tstamps & ix for sig sections that exceed thresh len ----------------
if ~isempty(line_ix{1}) && ix_thresh > 0
    
    % start ix counter
    ix = 1; 
    
    for i = 1:length(line_ix)        
        
        % determine if section exceeds threshold length
        if line_ix{i}(2) - line_ix{i}(1) > ix_thresh
            
            % identify center of current bar
            
            % add ix
            text((line_ix{i}(1))/60, vert_offset(2), num2str(ix), 'FontSize', 9, 'fontname', 'Verdana') 
            
            % increment ix counter
            ix = ix + 1;
            
            % add timestamps
            if t_stamp
                
                % get timing of current sig section
                times = {datestr(seconds(line_ix{i}(1)), 'MM:SS') datestr(seconds(line_ix{i}(2)), 'MM:SS')};

                % write text
                text((line_ix{i}(1))/60, vert_offset_tstamp(1), [times{1} '-' times{2}], 'FontSize', 5)
                
                % flip vertical offset so next tstamp is staggered
                vert_offset_tstamp = fliplr(vert_offset_tstamp);
            end
        end
    end
end

% calculate prop of time pts that are sig (to 2 decimals)
prop_sig = length(sig_pts)/length(yobs)*100;
prop_sig = round(prop_sig, 2);
output.propsig = prop_sig;

% add text to state proportion
text(9.2, -.88, [num2str(prop_sig) '% RTs sig.'], 'FontSize', 12, 'fontname', 'Verdana')

% legend
%lgd = legend([f(1), f(2) pl(1)], 'Observed', 'Permutation', 'Confidence interval', 'Location', 'northwest', 'FontSize', 14);

% figure mods -------------------------------------------------------------
% labels
xlabel('Story time (min)')
ylabel('Normalized RT (s)')
%title([story_tag ' (' online_tag ')'], 'fontweight', 'normal')
title(story_tag, 'fontweight', 'normal')

% axis props
set(gca, 'xlim', [0 13.5])
set(gca, 'ylim', [-1 1])
set(gca, 'xtick', 0:3:12)
set(gca, 'ytick', -1:.5:1)

set(gca, 'fontname', 'Verdana') 
set(gca, 'FontSize', 16)
set(gca, 'linewidth', 2)

% figure props
box off

% legend props
%lgd.ItemTokenSize = [17, 12]; % size of lines in legend
%legend boxoff

hold off

% save --------------------------------------------------------------------
if do_save
    save(savefile_stats, '-struct', 'output')
    
    if length(wsize) == 1
        print(gcf, savefile_fig, '-dsvg', '-r300')
    end
end
pause(.1)
end
end
end
end
% Purpose: 
% INPUT:
% OUTPUT: 

clc
clear all;
cd('C:\Users\mbain25\Documents\efflisb\Analysis')
rng(257)

% define file params*****
story      = 1;
n_sets     = 1000;
smoothtype = 'rect'; % 'rect' -> rectangular window movavg; 'gauss' -> gaussian window; 'lowess' -> local linear regression w/ 1st degree polynormial   
span       = 10;

% define params*****

% read data
supsets = load(['Supersets/supsets_story_' num2str(story) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_span_' num2str(span) '.mat']);
supsets = supsets.supsets;
ts = load(['Group_Processing/ts_norm_1_smoothed_0.mat']);
ts = ts.ts;

% get savefile
savefile = ['figures_var' '/' 'varanalysis' '_story_' num2str(story) '_smoothtype_' smoothtype '_span_' num2str(span) '_pcrit_' num2str(p_crit) '.png'];
if exist(savefile, 'file')
    error(['file ' savefile ' exists already. Either delete or use new parameters!'])
end



savefile = ['Plotting' filesep 'ts_varanalysis' '_story_' num2str(s) '_smoothed_' num2str(span) '_pcrit_' num2str(p_crit) '.png'];
savefile_stats = ['Group_Processing' filesep 'stats_varanalysis' '_story_' num2str(s) '_smoothed_' num2str(span) '_pcrit_' num2str(p_crit) '.mat'];

% call function to permute each stack n_sim times to generate surrogates
[ts_perm, ts_perm_av] = permute_ts(ts(s).ystacked, size(ts(s).ystacked, 1));

% get mean ts
ts_av = nanmean(ts(s).ystacked);

% apply smoothing to permuted ts w/ same params as stacks
if span == 0
else
    for i = 1:size(ts_perm_av, 1)
        ts_perm_av(i, :) = smooth(ts_perm_av(i, :), span);
    end
end
% compute columnwise t-test btwn bootstrapped & original ts
[h, p, ci, stats] = ttest2(ts(s).ystacked, ts_perm_av, 'vartype', 'unequal');

% store stats
tscomp_stats(s).stats = stats;
tscomp_stats(s).p     = p;

% OVERLAY PERMUTED & ORIGINAL TS W/ SIG INDICATED BY COLOUR
% figure colour scheme
colour_scheme = {[0 0 0]/255, [0 80 191]/255, [255 17 0]/255};

% stats on stacks
y   = nanmean(ts(s).ystacked); 
std = nanstd(ts(s).ystacked);
n   = sum(~isnan(ts(s).ystacked));
sem = std./sqrt(n);               

%%{
% calculated CI for each point
tcrit  = tinv(p_crit, size(ts(s).ystacked, 1) - 1); % t-stat required for sig 
CI     = tcrit*sem;                                 % confidence interval size                
%}

% original ts w/ shaded error bars or confidence intervals*****
figure
fig = shadedErrorBar(ts(s).x/60, y, sem); % error bars
%fig = shadedErrorBar(ts(s).x/60, y, CI); % CIs

hold on
plot(ts(s).x/60, ts_av, 'color', colour_scheme{1}, 'LineWidth', 1)

% permuted ts
hold on
plot(ts(s).x/60, mean(ts_perm_av), 'color', colour_scheme{2}, 'LineWidth', 1)

% indicate significant differences in red
hold on
sig_pts = find(p < p_crit); 
scatter(ts(s).x(sig_pts)/60, ts_av(sig_pts), 'MarkerFaceColor', colour_scheme{3},'MarkerEdgeColor', colour_scheme{1})

% ADD HORIZONTAL LINES TO PLOT BELOW SIGNIFICANT SECTIONS
% initialize arrays to store start and endpoints of lines
line_ix(s).ix = cell(1); 

% starting index of all lines
line_startix = find(diff(sig_pts) > 1) + 1;   

% draw first line
if line_startix(1) > 2

    % get and store starting and ending index
    line_start         = sig_pts(1);
    line_end           = sig_pts(line_startix(1) - 1);
    line_ix(s).ix{end} = [ts(s).x(line_start) ts(s).x(line_end)];

    % draw line
    line([ts(s).x(line_start) ts(s).x(line_end)]/60, [-.12 -.12], 'LineWidth', 3, 'Color', colour_scheme{3});
end

% draw all remaining lines
for i = 1:length(line_startix)

    % draw last line
    if i == length(line_startix)

        % get and store starting and ending index
        line_start             = sig_pts(line_startix(i));
        line_end               = sig_pts(end);
        line_ix(s).ix{end + 1} = [ts(s).x(line_start) ts(s).x(line_end)];

        % draw line
        line([ts(s).x(line_start) ts(s).x(line_end)]/60, [-.12 -.12], 'LineWidth', 3, 'Color', colour_scheme{3});   

    else
        % get and store starting and ending index
        line_start                 = sig_pts(line_startix(i));
        line_end                   = sig_pts(line_startix(i + 1) - 1);

        if cellfun(@isempty, line_ix(s).ix)
            line_ix(s).ix{end}     = [ts(s).x(line_start) ts(s).x(line_end)];
        else
            line_ix(s).ix{end + 1} = [ts(s).x(line_start) ts(s).x(line_end)];
        end

        % draw line
        line([ts(s).x(line_start) ts(s).x(line_end)]/60, [-.12 -.12], 'LineWidth', 3, 'Color', colour_scheme{3});
    end
end

% labels
title({['Original RT Time Series vs. Permuted:'], ['Story ' num2str(s)]})
ylabel('RT (s, normalized)')
xlabel('Story Time (min)')

% axis params
%ylim([.5 1.1])
ylim([-.2 .2])

% add text to indicate sections of story
for i = 1:length(line_ix(s).ix)

    % get timing of current sig section
    time = datestr(seconds(line_ix(s).ix{i}(1)), 'MM:SS');

    % write text
    %text((line_ix(s).ix{i}(2))/60, .6, num2str(time))
    text((line_ix(s).ix{i}(2))/60, -.13, num2str(time), 'FontSize', 6)
end

box off
hold off

% save
save(savefile_stats, 'tscomp_stats')
print(gcf, savefile, '-dpng', '-r300')
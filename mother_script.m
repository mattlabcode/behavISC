clc
clear;
format compact

%cd('C:\Users\mbain25\Documents\efflisb\Analysis')    % experiment directory
cd('C:\Users\mbain25\Documents\efflisb\Analysis')
datadir = 'C:\Users\mbain25\Documents\efflisb\Data'; % data directory 

addpath(datadir);

addpath(genpath('Matlab'))                           % all scripts and functions in experiment directory

addpath(genpath('Timings'))                          % trial timings
addpath(genpath('Stories'))                          % stories

addpath(genpath('Group_Processing'))                 % output directories
addpath(genpath('Individual_Processing')) 
addpath(genpath('Plotting'))   

%% ------------------------------------------------------------------------
% RUN DIAGNOSTICS 
curr_file = 1; % file to start analysis from

run_diagnostics(datadir, curr_file)

%% ------------------------------------------------------------------------
% AGGREGATE AND SMOOTH TIME SERIES
%***store number of misses for each time point
%***store versions that have time point for each time point in cell array

norm = 1;  % normalized or non-normalized data
span = 0; % size of smoothing window (in time points; 0 = no smoothing)

datadir = 'E:\MBain Computer Dump\Users\mbain25\Documents\efflisb\Analysis\Individual_Processing';

process_tseries(datadir, norm, span)

%% ------------------------------------------------------------------------
% PLOT TIME SERIES 
clc
close all

norm     = 1;
span     = 10;
story    = 1;

groupdir = dir(fullfile('Group_Processing', '*.mat'));

plot_tseries(groupdir, norm, span, story)

%% ------------------------------------------------------------------------
% PLOT DIAGNOSTIC INFO
clc
close all

data = load('Group_Processing\av_processed', 'story_dx');
data = data.story_dx;

plot_dx(data)

%% ------------------------------------------------------------------------
% COMPUTE TIME SERIES CORRELATIONS***SCRAP
clc
close all

norm     = 1;
span     = 10;
story    = 1;
n_perm   = 100; % r computed on surrogate data (generated over n_sim its) n_boot times

groupdir = dir(fullfile('Group_Processing', '*.mat'));

compute_sim(groupdir, story, norm, span, n_perm)

%% ------------------------------------------------------------------------
% ANALYZE SYSTEMATIC VARIABILITY OF TS 
clc
close all

span   = 20;
p_crit = .05;

data   = load(['Group_Processing\ts_norm_1_smoothed_' num2str(span)]); 
ts     = data.ts;

analyze_var(ts, span, p_crit)

%% ------------------------------------------------------------------------
% RUN STATS ON UPPER/LOWERCASE RT AVERAGES 
%***also plot this

data     = load('Group_Processing\av_processed', 'story_dx');
story_dx = data.story_dx;

% run paired ttest
[h, p, ci, stats] = ttest(story_dx(1).avRT_upp, story_dx(1).avRT_low, 'Vartype', 'unequal');

%% ------------------------------------------------------------------------
% CALCULATE AVERAGE RT ON AVG AGGREGATE TS OVER ALTERNATING SEGMENTS***SCRAP
%***figure out why mean for last noise segment for first story is nan

norm     = 1;
smoothed = 1;
n_seg    = 14; % # of N/C segments (should be even integer)

groupdir = dir(fullfile('Group_Processing', '*.mat'));

comp_noiseclear(groupdir, norm, smoothed, n_seg)

%% 
% plot dual-task instructional image

savefile = ['figures_schematic' '/' 'barplot' '_dualtask_pred' '.svg'];

% set data to plot
x   = [1 2];
y   = [500 1000]';

% set plotting params
colours    = [190 22 34; 7 162 228]/255;
bar_width  = .4;

% plot bars
f(1) = bar(1, y(1), 'BarWidth', bar_width, 'FaceColor', colours(1, :), 'LineWidth', 2);  
hold on
f(2) = bar(2, y(2), 'BarWidth', bar_width, 'FaceColor', colours(2, :), 'LineWidth', 2);                

% add labels
set(gca, 'XTickLabel', {'Low', 'High'}, 'fontsize', 28)
ylabel('Mean RT (ms)', 'fontsize', 28)
xlabel('Cognitive Load', 'fontsize', 28)

% axis params
set(gca, 'linewidth', 4)
set(gca, 'FontSize', 28)
set(gca, 'fontname', 'Verdana') 
ylim([0 1000])

box off

hold off

print(gcf, savefile, '-dsvg', '-r300')

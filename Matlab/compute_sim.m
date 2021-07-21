function compute_sim(groupdir, story, norm, span, n_perm)
%% ------------------------------------------------------------------------ 
% LOAD IN REQUESTED DATA 
% define file based on specified parameters
file_ix = strfind({groupdir.name}, ['ts_' 'norm_' num2str(norm) '_smoothed_' num2str(span)]);
file = find(cellfun(@nnz, file_ix));

if isempty(file)
    % error if data w/ given parameters not stored 
    disp('Error: The requested file does not exist.')
    return
else
    % load corresponding data 
	file = groupdir(file).name;
    data = load(file);
    ts = data.ts;
end

%% ------------------------------------------------------------------------ 
%{
% METHOD 1: AVG ALL PAIRWISE CORRELATIONS BTWN SSs 
% savefile for correlation matrix and stats
savefile = ['Plotting' filesep 'corrmatrix_pairwise_' 'norm_' num2str(norm) '_smoothed_' num2str(span) '_story_' num2str(story) '.png'];
savefile_stats = ['Group_Processing' filesep 'stats_pwcorr_' 'norm_' num2str(norm) '_smoothed_' num2str(span) '_story_' num2str(story) '.mat'];

% store pw r/p values
[r, p] = corrcoef(ts(story).ystacked.', 'rows', 'pairwise');

% compute mean r/p value
mean_r = mean(r((triu(r, 1) ~= 0)));
mean_p = mean(p((triu(r, 1) ~= 0)));

% compute p for mean r***revisit this at some point
meanr_df = size(ts(story).ystacked, 1) - 1;
meanr_t  = mean_r*(sqrt(meanr_df/(1 - mean_r)));
meanr_p  = 2*(1 - tcdf(abs(meanr_t), meanr_df));

% store stats in output struct
pwcorr_stats.r = r;
pwcorr_stats.p = p;
pwcorr_stats.mean_r = mean_r;
pwcorr_stats.mean_p = mean_p;
pwcorr_stats.meanr_p = meanr_p;

fprintf('pairwise mean r value: %f, p = %f\n', mean_r, mean_p);

% set params for corr matrix
n = length(r);                              
M = r;                                                   % n x n correlation matrix
L = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10'}; % label for each bin

% plot correlation matrix
figure
imagesc(M);                 % plot the matrix
set(gca, 'XTick', 1:n);     % center x-axis ticks on bins
set(gca, 'YTick', 1:n);     % center y-axis ticks on bins
set(gca, 'XTickLabel', L);  % set x-axis labels
set(gca, 'YTickLabel', L);  % set y-axis labels
caxis([0 1]);               % set colourbar limits
title({'Pariwise Dual-Task RT Correlations:', ['Story ' num2str(story)']}, 'FontSize', 14); % set title
colormap('jet');            % set colorscheme
colorbar()                  % enable colorbar

%fcnCorrMatrixPlot(ts(story).ystacked.', L, 'hello')

save(savefile_stats, 'pwcorr_stats')
print(gcf, savefile, '-dpng', '-r300')
%}
%% ------------------------------------------------------------------------
% VALIDATION: PAIRWISE CORR OF POOLED TS 
%{
savefile = ['Group_Processing' filesep 'stats_pwcorrpools_' 'norm_' num2str(norm) '_smoothed_' num2str(span) '_story_' num2str(story) '.mat'];

for v = 1:size(ts(story).ypool, 2)
    % calculate pairwise r/p and store
    tspool = cell2mat(ts(story).ypool(:, v));
    [r, p] = corrcoef(tspool.', 'rows', 'pairwise');
    r_all(v) = {r};
    p_all(v) = {p};
    
    % get means
    mean_r(v) = mean(r((triu(r, 1) ~= 0)));
	mean_p(v) = mean(p((triu(p, 1) ~= 0)));

    fprintf('pairwise mean r value on version %d pool: %f\n', v, mean_r(v));
end

% save stats
pwcorrpools_stats.r_all = r_all;
pwcorrpools_stats.p_all = p_all;
pwcorrpools_stats.mean_r = mean_r;
pwcorrpools_stats.mean_p = mean_p;

save(savefile, 'pwcorrpools_stats')
%}
%% ------------------------------------------------------------------------ 
% PERMUTE RESAMPLED SUPERSETS TO COMPARE WITH ORIGINAL r VALUES
savefile_stats    = ['Group_Processing' filesep 'stats_shalfcorr_' 'norm_' num2str(norm) '_smoothed_' num2str(span) '_story_' num2str(story) '.mat'];
savefile_bootdist = ['Plotting' filesep 'permdist_shalfcorr_' 'norm_' num2str(norm) '_smoothed_' num2str(span) '_nperm_' num2str(n_perm) '_story_' num2str(story) '.png'];

[r_obs_all, p_all] = corr_ts(ts(story).ystacked);
r_obs_av = mean(r_obs_all);

% initialize array to store resampling stat
r_perm_av = zeros(1, n_perm);

% create wait bar to strack progress of resampling
f = waitbar(0,'Please wait...');

for i = 1:n_perm
    % update wait bar
    waitbar(i/n_perm, f, ['Permutation iteration ' num2str(i) ' / ' num2str(n_perm)]);
    
    % permute each ts over n_sim iterations*****simple perm or phase rand
    %***fix this function so that generates one surrogate for each ts in array
    [ts_perm, ts_perm_av] = permute_ts(ts(story).ystacked, 1);
    
    % initialize array to store surrogate superset
    superset_surr = zeros(size(ts(story).ystacked));
    
    for j = 1:size(ts(story).ystacked, 1) 
        superset_surr(j, :) = surrogate(ts(story).ystacked(j, :), 1, 'FT', 0, .5);
    end
    
    %%{
    % apply smoothing to permuted ts w/ same params as stacks
    for j = 1:size(ts_perm_av, 1)
        ts_perm_av(j, :) = smooth(ts_perm_av(j, :), span);
    end
    %}
    
    % compute r using independent groups method
    [r_all, p_all] = corr_ts(superset_surr);

    % store permutation stat
    r_perm_av(i) = mean(r_all);
end

% close wait bar
delete(f)

% p as proportion of times bootstrapped stat exceeded observed
p_boot_prop = (nnz(r_perm_av > r_obs_av))/n_perm;

% p as determined by z-score
z_r = (r_obs_av - mean(r_perm_av))/std(r_perm_av);
p_boot_z = 1 - normcdf(z_r);

fprintf('independent groups r value: %f, p = %f\n', mean(r_obs_all), p_boot_prop);


% SAVE
% save stats
groupcorr_stats.r_all_obs = r_all;
groupcorr_stats.mean_r_obs = r_obs_av;
groupcorr_stats.p_boot = p_boot_prop;

save(savefile_stats, 'groupcorr_stats')
%% ------------------------------------------------------------------------  
% PLOT PERMUTATION DISTRIBUTION
%null_dist   = r_perm_av;
%alt_dist    = r_obs_all;

% define values for curves to compare
null_dist   = randn(1, 1000)*10;
alt_dist    = randn(1, 1000)*10 + 30;
nbins_null  = 60;
colour_null = [137 117 232]/255;
colour_alt  = [43 204 197]/255;
p           = .001;

plot_permdist(null_dist, alt_dist, nbins_null, colour_null, colour_alt, p)

% labels
title({'Correlation vs Permutation Distribution:', [' Story ' num2str(story)]})
ylabel('Frequency')
xlabel('r (Pearson Correlation Coefficient)')

%print(gcf, savefile_bootdist, '-dpng', '-r300')
%%
% number of bins
nbins_null = 60;
nbins_alt  = round((nbins_null/range(null_dist))*range(alt_dist));

figure
p1 = histogram(null_dist, nbins_null, 'FaceColor', [137 117 232]/255);

% overlay distribution of observed r values
hold on
p2 = histogram(alt_dist, nbins_alt, 'FaceColor', [43 204 197]/255);

%***
%{
p1 = histfit(null_dist, nbins_null, 'normal');

% colour scheme
p1(1).FaceColor = [137 117 232]/255;
p1(2).Color     = p1(1).FaceColor*.5;

% adjust bar transparency
hpatch = get(p1(1), 'children');
set(hpatch, 'FaceAlpha', 0.1);

p2 = histfit(r_all_obs, nbins_alt, 'normal');

% colour scheme
p2(1).FaceColor = [43 204 197]/255;
p2(2).Color     = p2(1).FaceColor*.5;

% adjust bar transparency
hpatch = get(p2(1), 'children');
set(hpatch, 'FaceAlpha', 0.1);
%}

colormap([.5 .5 .5])

% labels
title({'Correlation vs Permutation Distribution:', [' Story ' num2str(story)]})
ylabel('Frequency')
xlabel('r (Pearson Correlation Coefficient)')

% add line for means
p3 = line([r_obs_av, r_obs_av], ylim, 'LineWidth', 2, 'Color', [255 53 53]/255);

% 95% confidence interval on r_obs
CI_r_boot = quantile(r_perm_av, [0.025 0.975]);

% add lines for 95% CIs
p4 = line([CI_r_boot(1), CI_r_boot(1)], ylim, 'LineWidth', 2, 'Color', [0 0 0], 'LineStyle', '--');
line([CI_r_boot(2), CI_r_boot(2)], ylim, 'LineWidth', 2, 'Color', [0 0 0], 'LineStyle', '--');

legend([p1 p2 p3 p4], 'Permutation Distribution', 'Observed Distribution', ['r (observed) = ' num2str(r_obs_av) ', p = ' num2str(p_boot_prop)], '95% CI', 'Location', 'northwest')

box off



print(gcf, savefile_bootdist, '-dpng', '-r300')

%% ------------------------------------------------------------------------
% PLOT ROC CURVE
% savefile
savefile_roc = ['Plotting' filesep 'roc_shalf_' 'norm_' num2str(norm) '_smoothed_' num2str(span) '_nperm_' num2str(n_perm) '_story_' num2str(story) '.png'];

% define values for curves to compare
null_dist = r_perm_av;
alt_dist  = r_obs_all;

% set upper/lower thresh as upper edge of null dist (FPR = 0) & lower edge of null (FPR = 100) - lower edge of alt (TPR = 100)
thresh_bounds = [max(null_dist), min(null_dist)];

% set step size and number of steps for incrementally adjusting threshold
step_size = .001;
steps     = ceil((thresh_bounds(1) - thresh_bounds(2))/step_size);

% set starting threshold
thresh = thresh_bounds(1);

% initialize arrays to store FP/TRP at each thresh
[tpr, fpr] = deal(zeros(1, steps));

for i = 1:steps
    % get TPR and FPR for current thresh
    tpr(i) = nnz(alt_dist > thresh)/length(alt_dist);
	fpr(i) = nnz(null_dist > thresh)/length(null_dist);
    
    % decrement threshold
    thresh = thresh - step_size;
end

% get d' ( = sqrt(Z)(AUC), Z = inverse of CDF of normal distribution)
auc = abs(trapz(fpr, tpr));
dprime = sqrt(2)*norminv(auc);

% plot TPR vs FPR
figure
plot(fpr, tpr, 'LineWidth', 1.5, 'Color', [0 0 0])
%roc = roc_curve(null_dist, alt_dist);

% set axis params
set(gca, 'TickLength', [0.01 0]);

% add labels
title('ROC Curve')
ylabel('TPR (True Positive Rate)')
xlabel('FPR (False Positive Rate)')

title({['ROC Curve - Correlation vs Permutation Dist:'], ['Story ' num2str(story)]})

% add d' prime to plot
annotation('textbox', [.7 .8 .1 .1], 'String', ['dprime: ' num2str(dprime)], 'FitBoxToText', 'on');

box off

% save
print(gcf, savefile_roc, '-dpng', '-r300')

% bootstrap correlation
%[R, Rpc, Rsd, Rpt, Z, Zpc, Zsd, Zpt] = bbcorr(ts(story).ystacked(1:2, :)', 5, 5000, 0.9, 1, 0)

%% ------------------------------------------------------------------------
%{
% METHOD 3: AVG LEAVE-ONE-OUT CORR ON STACKS 
n_subj = size(ts(story).ystacked, 1);

% track all leave-one-out correlations
r_all = zeros(1, n_subj);

for i = 1:n_subj
    % get supersubj i ts
    ts_i = ts(story).ystacked(i, :);
    
    % average over all ts except i
    ts_mean = nanmean(ts(story).ystacked(1:end ~= i, :));
    
    % compute corr between ts i and leave-one-out ts_mean
    r = (corrcoef(ts_i, ts_mean, 'rows', 'pairwise'));
    r_all(i) = r(2);
end

% calculate mean leave-one-out ISC 
r_mean = mean(r_all);

fprintf('leave-one-out r value: %f\n', r_mean);
%}
end

function [subjset] = subjset_gen(tspool)
% define # SSs based on ts version pool w/ fewest S
n_SS  = min(sum(~cellfun(@isempty, tspool)));

% number of timing versions
n_ver = size(tspool, 2);

% get random selection of n_SS S from each version pool
subjset_ix = randperm(n_SS); 

% preallocate array to store new SS
subjset = cell(n_ver, 1); 

for i = 1:n_ver
    % draw S from appropriate version pool to add to SS 
    subjset(i) = tspool(subjset_ix(i), i); 
end

subjset = nanmean(cell2mat(subjset));
end

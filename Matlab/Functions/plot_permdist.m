function plot_permdist(null_dist, alt_dist, nbins_null, colours)
% Function 
% INPUT:
%   null_dist   array 
% Author: Matthew Bain

% set colour of lines outlining hist
colour_histedges = [50 50 50]/255;

% number of bins
nbins_alt  = round((nbins_null/range(null_dist))*range(alt_dist));

figure
f(1) = histogram(null_dist, nbins_null, 'FaceColor', colours(2, :), 'EdgeColor', colour_histedges);

% overlay distribution of observed r values
hold on
f(2) = histogram(alt_dist, nbins_alt, 'FaceColor', colours(1, :), 'EdgeColor', colour_histedges);

colormap([.5 .5 .5])

% define ylim based on hist heights
ylim = [0 max(max(f(1).Values), max(f(2).Values))];

% add line for means
%l(1) = line([mean(alt_dist), mean(alt_dist)], ylim, 'LineWidth', 2, 'Color', [255 53 53]/255, 'LineWidth', 2);

% 95% confidence interval on r_obs
CI_r_boot = quantile(null_dist, [0.025 0.975]);

% add lines for 95% CIs
l(2) = line([CI_r_boot(1), CI_r_boot(1)], ylim, 'LineWidth', 2, 'Color', [80 80 80]/255, 'LineStyle', '--', 'LineWidth', 2);
line([CI_r_boot(2), CI_r_boot(2)], ylim, 'LineWidth', 2, 'Color', [80 80 80]/255, 'LineStyle', '--', 'LineWidth', 2);

lgd = legend([f(1) f(2)], 'Permutation', 'Case Judgement', 'Location', 'northwest');
%lgd = legend([f(1) f(2) l(1) l(2)], 'Permutation Distribution', 'Observed Distribution', ['r (observed) = ' num2str(mean(alt_dist)) ', p = ' num2str(p)], '95% Confidence Intervals', 'Location', 'northwest');

% legend params
set(lgd, 'FontSize', 28);
%lgd.ItemTokenSize = [15, 28]; % size of lines in legend
legend boxoff

box off

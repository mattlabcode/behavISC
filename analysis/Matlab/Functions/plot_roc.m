function [p, dprime] = plot_roc(null_dist, alt_dist, colour)
% Function 
% INPUT:
%   null_dist   array 
% Author: Matthew Bain

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
p = plot(fpr, tpr, 'Color', colour, 'LineWidth', 4);

% set axis params
%set(gca, 'TickLength', [0.01 0]);

box off
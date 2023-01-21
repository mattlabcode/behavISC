function [r_all, p_all] = corr_ts(ts)
% Correlate ts by enumerating all combinations of: 
% 2 distinct groups of n/2 ts from set of n and computing r on group averages.
% INPUT:
%   ts          array of ts to be correlated, arranged by row
% OUTPUT:
%   r_all       row vector containing all correlations over each combination
%   p_all       row vector containing all corresponding p values

% number of time series
size_ts = round(size(ts, 1));

% adjust useable time series so that number divisible by 2
if mod(size_ts, 2) > 0
    size_ts = size_ts - 1;
end

% get number of ts in each group
n_ts = size_ts / 2;

% get all distinct combinations of 2 groups of four from set of 8
pool = 1:size_ts;
combos = combnk(pool, n_ts);

% eliminate half of combos that are mirror reflections (for efficiency -
% but note that this redundancy will not affect outcome of correlation)
for c = 1:length(combos)/2
    % find index of corresponding mirror-image combo and eliminate
    combos(ismember(combos, find(~ismember(1:size_ts, combos(c, :))), 'rows'), :) = [];
end

% initialize struct to store bootstrapped stats
[r_all, p_all] = deal(zeros(1, length(combos)));

for c = 1:length(combos)
    % get current combination for each group
    ts1_combo = combos(c, :);
    ts2_combo = pool(~ismember(pool, combos(c, :)));
    
    % calculate mean ts for each group
    ts1 = nanmean(ts(ts1_combo, :));
    ts2 = nanmean(ts(ts2_combo, :));
    
    % correlate
    [r, p] = corrcoef(ts1, ts2, 'rows', 'pairwise');
    
    % store stats
    r_all(c) = r(2);
    p_all(c) = p(2);
end

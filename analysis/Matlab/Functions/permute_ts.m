function [ts_perm, ts_perm_av] = permute_ts(ts_stacked, n_sim)
% Function randomizes an array of ts by permuting data points.
% INPUT: 
%   ts             array of ts to be correlated, arranged by row
%   n_sim          number of permutations to generate new ts rows
% OUTPUT: 
%   surrogate_ts   permuted array of ts
% Author: Matthew Bain

% bootstrapping params
n = size(ts_stacked, 1);
N = length(ts_stacked);

% permute each row of ts n_sim times and store (all + avg)
ts_perm    = cell(n, 1);
ts_perm_av = zeros(n, N);

for i = 1:n        
    % store all permutations for current row
    ts_perm_curr = zeros(n_sim, N);
    
    for j = 1:n_sim
        % generate permutations for current stack
        ts_perm_curr(j, :) = ts_stacked(i, (randperm(N)));
    end

    % store current stack permutations and average 
    ts_perm{i}       = ts_perm_curr;
    ts_perm_av(i, :) = mean(ts_perm_curr);
end

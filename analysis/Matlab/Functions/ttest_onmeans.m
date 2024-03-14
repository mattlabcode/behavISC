function [tval, p, df, ci] = ttest_onmeans(mean1, mean2, var1, var2, n1, n2, alpha)
% [tval, p, df, ci] = ttest_onmeans(mean1, mean2, var1, var2, n1, n2, alpha)
% computes independent 2 sample t-test based on group means/variances
% assuming unequal variances (use Welch's t-test) and two-tailed test
% INPUT:  
%   mean1/mean2 - sample means (1*n, n = # groups to compare)
%   var1/var2   - sample variances
%   n1/n2       - sample sizes
%   alpha       - p-value cutoff
% OUTPUT: 
%   tval        - t-stats 
%   p           - p-value
%   df          - degrees of freedom
%   ci          - confidence intervals (2*n, rows = upper/lower bounds)

% check that all data vectors equal length
if ~isequal(length(mean1), length(mean2), length(var1), length(var2), length(n1), length(n2))
    error(message('data vectors must be equal length'));
end

% initialize arrays to store output
[tval, p, df] = deal(zeros(1, length(mean1)));
ci            = zeros(2, length(mean1));

% go through all groups
for i = 1:length(mean1)

    % compute df (Welch–Satterthwaite equation)
    df1   = (var1(i)/n1(i) + var2(i)/n2(i))^2;
    df2   = (((var1(i)/n1(i))^2) / (n1(i) - 1)) + (((var2(i)/n2(i))^2) / (n2(i) - 1));
    df(i) = df1/df2;
    
    % compute t
    ser     = sqrt(var1(i)/n1(i) + var2(i)/n2(i));
    xdiff   = mean1(i) - mean2(i);
    tval(i) = xdiff/ser;    

    % compute p-value (2-tailed)
    p(i) = 2 * tcdf(-abs(tval(i)), df(i));
    
    % get critical t value
    tcrit = tinv((1 - alpha/2), df(i));
    
    % compute confidence intervals
    ci(:, i) = [mean1(i) + tcrit*ser; mean1(i) - tcrit*ser];
end

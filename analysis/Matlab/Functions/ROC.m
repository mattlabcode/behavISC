function [AUC pHit pFA criteria] = ROC(x,y,nc,pflag)

% [AUC pHit pFA criteria] = ROC(x,y,nc,pflag) 
%
% x, y - vectors of firing rates, one for each condition, x is defined as hit
% nc - number of criteria
% pflag - plot 0 or 1
%
% AUC  - area under the curve
% pHit - proportion of hits for each criterion
% pFA  - proportion of false alarms for each criterion
% criteria - criteria used to calculate pHit and pFA
%
% References:
% http://www.mbfys.ru.nl/~robvdw/DGCN22/PRACTICUM_2011/LABS_2011/ALTERNATIVE_LABS/Lesson_10.html
% Britten et al. (1992) J Neurosci. 12(12):4745-65.
%
% Description: The script calculates the receiver operator curve for an ideal observer.
% --------------------------------------------------------------------------------------------------
% B. Herrmann, Email: bjoern.herrmann@outlook.com, 2015-02-24

% define criteria space
mi  = min([x(:);y(:)]);
ma  = max([x(:);y(:)]);
r   = (ma - mi)*1.1;
low = mi-(r-r/1.1)/2;
criteria = linspace(low,low+r,nc);

% calculate hits and false alarms
[pHit pFA] = deal(zeros([length(criteria) 1]));
for ii = 1 : length(criteria)
	pHit(ii) = sum(x > criteria(ii))/length(x);
	pFA(ii)  = sum(y > criteria(ii))/length(y);
end

% calculate area under the curve
AUC = -trapz(pFA,pHit);

% do plotting
if pflag
	figure, set(gcf,'Color',[1 1 1]), hold on
	set(gca,'FontSize',14)
	plot(pFA,pHit,'-r','LineWidth',2);
	%plot(pFA,pHit,'sr','LineWidth',2,'MarkerSize',5);
	plot([0,1],[0,1],'--k','LineWidth',2);
	xlabel('p(y>=crit)');
	ylabel('p(x>=crit)');
	axis([0 1 0 1])
	axis square
end

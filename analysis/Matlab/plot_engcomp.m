% PURPOSE: compare correlation values for median splits.

clc
%close all
clear all;
cd('C:\Users\mbain25\Documents\efflisb\Analysis')
addpath(genpath('matlab'))
rng_cfig = rng(257);

% define file params
online     = 1;
story      = 1;
n_sets     = 50;
smoothtype = 'lowess';
wsize      = 10;
split_by   = 'absorption';               
split_trun = 1;         
n_perm     = 500;

% define params
n = [50 50]; % size of group 1 & 2 (= # Ss in each set of SSs)

%colours    = [219 115 122; 102 195 235]/255; % R/B light
%colours    = ([190 22 34; 7 162 228] + 20)/255; % R/B dark
colours   = [106 219 214; 189 178 243]/255; % G/P light
%sig = {[1 2]};
yli = [-.05 .4];
filetype   = 'svg';
linewidth  = 1.5;

do_save    = 0;
overwrite  = 0;

% -------------------------------------------------------------------------
% input
[r, err] = deal(zeros(1));
for s = 1:2 %***replace w/ below if run splits again w/ 'trunc'
    input = load(['correlations/archive - b4 updating surr 280720' '/' 'corr' '_story_' num2str(story) '_nsets_' num2str(n_sets) ...
        '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '_dosplit_1' ...
        '_splitby_' split_by '_split_' num2str(s) '_splittrun_' num2str(split_trun) ...
        '_nperm_' num2str(n_perm) '_online_' num2str(online) '.mat']);

    %{
    input = load(['correlations' '/' 'corr' '_story_' num2str(story) '_nsets_' num2str(n_sets) ...
        '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '_dosplit_1' ...
        '_splitby_' split_by '_split_' num2str(s) '_splittrun_' num2str(split_trun) ...
        '_trunc_no' '_nperm_' num2str(n_perm) '_online_' num2str(online) '.mat']);
    %}  
        
    % st1:low st1:high st2:low st2:high
    r(end + 1)   = input.robs_av;
    err(end + 1) = input.robs_var; 
end

r(1)   = []; % remove zero row
err(1) = [];

% savefile
savefile = ['figures_engcomp' '/' 'corrcomp' '_story_' num2str(1) '_nsets_' num2str(n_sets) ...
    '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '_splitby_' split_by ...
    '_splittrun_' num2str(split_trun) '_nperm_' num2str(n_perm) '_online_' num2str(online) ['.' filetype]];

% check if savefile already exists
if do_save
    if exist(savefile, 'file')
        if overwrite 
            disp('overwriting existing savefile')
        else
            error(['file ' savefile ' exists already. Either delete or use new parameters!'])
        end
    end
end

% -------------------------------------------------------------------------
% compute sig of correlation difference -----------------------------------
rz = atanh(r);
szab = sqrt(1/(n(1) - 3) + 1/(n(2) - 3));
z = abs(rz(2) - rz(1)) / szab; % high > low
p = 1 - normcdf(z, 0, 1); % 1-tailed
disp(p)

% plot --------------------------------------------------------------------
% set data
x   = [1 2];
err = sqrt(err);                          

% plot bars
f = bar(diag(r), 'stacked', 'BarWidth', .4, 'LineWidth', linewidth); hold on   
%'BarWidth', bar_width, 'FaceColor', colours2(1, :), 'LineWidth'); hold on        
set(f(1), 'facecolor', colours(2, :))
set(f(2), 'facecolor', colours(1, :))

% plot error bars
er = errorbar(r, err, 'linewidth', linewidth); hold on

% plot significance bar
if p > .05
    sig = nan;
else
    sig = p;
end
sigstar([1 2], sig); hold on
  
% error line props
er.Color     = [0 0 0];                            
er.LineStyle = 'none';  

%{
% legend
lgd = legend([f(1) f(2)], 'Low', 'High', 'Location', 'northwest'); hold on

% legend props
lgd = set(lgd, 'FontSize', 16);
lgd.ItemTokenSize = [.3, .3]; % size of boxes in legend
legend boxoff
%}

% add text to state r
text(x(1), r(1) + .05, num2str(r(1)), 'FontSize', 10, 'fontname', 'Verdana');
text(x(2), r(2) + .05, num2str(r(2)), 'FontSize', 10, 'fontname', 'Verdana');

% add labels
title(['Story ' num2str(story) ': ' split_by], 'fontweight', 'normal')
set(gca, 'XTickLabel', {'Low', 'High'});
ylabel('Correlation')

% figure props
box off 

% axis props
ylim(yli)
xtickangle(45)
yticks(-2.:.1:.5)
set(gca, 'linewidth', linewidth)
set(gca, 'fontsize', 16, 'fontname', 'Verdana') 
set(gca, 'TickLength', [0.005 0]);

% save --------------------------------------------------------------------
if do_save
	set(gcf,'PaperUnits', 'inches', 'PaperPosition', [0 0 2.85 4.5]*.9)
    print(gcf, savefile, ['-d' filetype])
end
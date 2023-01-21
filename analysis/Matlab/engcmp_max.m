clc
clear all;
cd('C:\Users\mbain25\Documents\efflisb\Analysis')
addpath(genpath('matlab'))
rng_cfig = rng(257);

% define file params
online     = [0 1]; % what data to include in scatter
story      = [1 2];

% define params
smoothtype = 'lowess';
wsize      = 10;

dim        = {'enjoyment'}; % dimension of engagement to perform analysis on             
canon_meth = 'indts';  % 'indts' -> use individual subj ts; 'setts' -> use smoothed SS ts in sets

colours    = [190 22 34; 7 162 228]/255; % R/B dark
%colours = [242 209 212; 199 236 255]/255; % R/B light light

font         = 'Verdana';
linewidth    = 2;
%paperposition = [0 0 3 4.5];

if length(dim) > 1
    fontsize = 9;
else
    fontsize = 16;
end

txt_on       = 1; % turn labels on or off

do_save   = 0;
overwrite = 0;

% -------------------------------------------------------------------------
% loop through each dimension specified
figure
pl = 1; % plot counter
for d = 1:length(dim)
    
% --- input ---
% ts in v pools w/ run labels & all ts
ts_il = load(['timeseries' '/' 'ts' '_dosplit_0' '_splitby_' '_splittrun_0' '_online_'  num2str(0) '.mat']); 
ts_ol = load(['timeseries' '/' 'ts' '_dosplit_0' '_splitby_' '_splittrun_0' '_online_'  num2str(1) '.mat']); 
ts_il = ts_il.ts;
ts_ol = ts_ol.ts;
%tspool_il = ts_il(story(s)).ypool; %***SCRAP
%tspool_ol = ts_ol(story(s)).ypool;

% group-level processing metrics
dxav_il = load(['Group_Processing' '/' 'av_processed' '_logtrans_0' '_online_' num2str(0) '.mat']); 
dxav_ol = load(['Group_Processing' '/' 'av_processed' '_logtrans_0' '_online_' num2str(1) '.mat']); 
dxav_il = dxav_il.dxav;
dxav_ol = dxav_ol.dxav;
    
% --- output ---
% get savefile tags
if length(dim) > 1
    dim_string = 'both';
else 
    dim_string = dim{1};
end
online_string = num2str(online);
story_string = num2str(story);
online_string(online_string == ' ') = []; % remove whitespace
story_string(story_string == ' ') = []; % remove whitespace

% savefile
savefile = ['figures_engcorr' '/' 'line' '_online_' online_string ...
    '_story_' story_string '_dim_' dim_string '_canonmeth_' canon_meth ...
    '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '.svg'];

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
% rank order ts/ratings w/in each v pool based on eng ratings -------------
% concat final eng ratings/corr values if online/inlab; story 1/2 combined 
set_resp      = [];
corr_ssranked = [];
for o = 1:length(online)
for s = 1:length(story)

if online(o) == 0
    ts = ts_il;
    dxav = dxav_il;
elseif online(o) == 1
    ts = ts_ol;
    dxav = dxav_ol;
end

% get dimensions
n_ss  = min(sum(ts(story(s)).run ~= 0)); % base # SSs on largest # Ss all pools have
n_ver = size(ts(story(s)).ypool, 2);

% truncate ypool & run ix
ts(story(s)).ypool = ts(story(s)).ypool(1:n_ss, :); 
ts(story(s)).run   = ts(story(s)).run(1:n_ss, :); 

tspool_ranked = cell([n_ss n_ver]);
[runs_ranked, resppool_ranked] = deal(zeros([n_ss n_ver]));

for v = 1:n_ver
    
    % file indices for current version
    curr_ix = ts(story(s)).run(:, v);
    
    % ranked ratings & indices of sorted ratings for curr v
    [resppool_ranked(:, v), sort_ix] = sort(dxav(story(s)).(lower(dim{d}))(curr_ix));
    runs_ranked(:, v) = curr_ix(sort_ix);
    
    % sort current v pool ts according to ranked ratings    
    tspool_ranked(:, v) = ts(story(s)).ypool(sort_ix, v); 
end

% --- get set w/ SS ts/ratings ordered from least to most engaged ---
n_tp         = length(tspool_ranked{1});
set_ts       = zeros([n_ss n_tp]);
set_resp_tmp = zeros([n_ss 1]);

for ss = 1:n_ss
    set_ts(ss, :) = nanmean(cell2mat(tspool_ranked(ss, :)'));
    set_resp_tmp(ss) = mean(resppool_ranked(ss, :));
     
    % smooth ranked SS ts
    if wsize
        % interpolate over any NaN values using linear regression
        ss_curr = set_ts(ss, :);
        t          = 1:n_tp;
        set_ts(ss, :) = interp1(t(~isnan(ss_curr)), ss_curr(~isnan(ss_curr)), t);

        % smooth
        set_ts(ss, :) = smoothdata(set_ts(ss, :), smoothtype, wsize);
    end
end
    
% --- get unique canonical EL tc for each SS (ts agg - SS ts) ---
ts_canon = zeros(size(set_ts));
corr_ssranked_tmp = zeros([n_ss 1]);

for ss = 1:n_ss

    % get indices wrt ypool of subj that make up current ss
    ss_poolix = ts(story(s)).run == runs_ranked(ss, :);

    % get av EL tc using all ts other than curr SS
    ts_canon(ss, :) = nanmean(cell2mat(ts(story(s)).ypool(~ss_poolix)));            

    % smooth
    if wsize
        ts_canon(ss, :) = smoothdata(ts_canon(ss, :), smoothtype, wsize);
    end

    % --- correlate each ranked SS ts w/ corresponding canon tc ---
    corr_ssranked_tmp(ss) = corr(set_ts(ss, :)', ts_canon(ss, :)');
end
    
% concatenate il/ol data
set_resp = vertcat(set_resp, set_resp_tmp);
corr_ssranked = vertcat(corr_ssranked, corr_ssranked_tmp); 

end
end

% -------------------------------------------------------------------------
% plot -------------------------------------------------------------------- 
subplot(1, length(dim), pl); hold on
pl = pl + 1; 

% compute spearman correlation
[r, p] = corr(set_resp, corr_ssranked, 'type', 'spearman', 'tail', 'right');
df     = length(set_resp)*2 - 2; 
fprintf('corr(mean %s rating for rank-ordered SSs, corr btwn SS ts & can. EL tc): r(%d) = %f, p = %f\n', lower(dim{d}), df, r, p)

% plot SS rating by corr w/ canonical EL tc (25 = circle sz)
fig = scatter(set_resp, corr_ssranked, 30, [.15 .15 .6], 'filled', 'MarkerFaceAlpha', .8); hold on

% --- plot line of best fit (degree 1 polynomial) ---
coef = polyfit(set_resp, corr_ssranked, 1); 
x1   = linspace(min(set_resp), max(set_resp), 1000);
y1   = polyval(coef, x1);
fit  = plot(x1, y1, 'color', colours(1, :), 'LineWidth', 1.6, 'linestyle', '--'); hold on

% --- figure mods ---
% add text
if txt_on
    % add text stating correlation to plot
    text(3.1, .35, ['{\itr}' ' = ' num2str(r, '%.3f')], 'FontSize', 14, 'fontname', 'Verdana')
    
    %dims = [.2 .6 .3 .3]; % coord of lower L corner of box wrt lower L corner of fig; width; height of annotation
    %annotation('textbox', dims, 'string', ['r(' num2str(df) '): ' num2str(r, '%.3f') ', p = ' num2str(p, '%.3f')], 'FitBoxToText', 'on', 'FontSize', 9, 'fontname', 'Verdana', 'EdgeColor', [1 1 1])
end

% labels
title([dim{d}], 'FontWeight','Normal')
ylabel('Corr. w/ canonical time course')
xlabel('Average rating (1 - 7)')

% axis props
set(gca, 'ylim', [-.2 .4])
%axis square

set(gca, 'TickLength', [.0075 0]);


set(gca, 'linewidth', linewidth)
set(gca, 'FontSize', fontsize)
set(gca, 'fontname', font) 

% figure props
box off
hold off

end

% save --------------------------------------------------------------------
if do_save
    print(gcf, savefile, '-dsvg', '-r300')
end
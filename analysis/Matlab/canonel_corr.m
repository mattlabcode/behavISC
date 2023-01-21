clc
clear all;
cd('/Users/matthewbain/Documents/Science/Experiments/efflisisc/analysis/')
addpath(genpath('matlab'))
rng_cfig = rng(257);

% define file params
test       = 'datasets'; % 'scrambling' -> nonsig corr btwn OL IT and OL Sc; 'datasets' -> sig corr btwn IL IT and OL IT
n_sets     = 50;
smoothtype = 'lowess';
wsize      = 10;
do_split   = 0;             
if do_split
    split_by   = 'enjoyment';               
    split_trun = 1;         
    split      = 1;           
else
    split_by   = '';             
    split_trun = 0;
    split      = 0;           
end

% define params
n_clear = 100; % 63
n_noise = 100; % 36

plot_type = 'incon2'; % con1; con2; incon1; incon2

colours1 = [106 219 214; 189 178 243]/255; % G/P light (low/high)
colours2 = [109 215 232; 70 240 189]/255; % B/G light (st1/st2)

font = 'Verdana';
fontsize = 16;
linewidth = 3;
%paperposition = [0 0 3 4.5];

do_save = 1;
overwrite = 1;

% -------------------------------------------------------------------------
% savefiles
savefile_stats = ['output_canonelcorr' '/' 'relana_' plot_type '_' test '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(min(wsize)) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) '.mat'];
savefile_fig = ['figures_canonelcorr' '/' 'relana_' plot_type '_' test '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(min(wsize)) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) '.svg'];

% check if savefile already exists
if do_save
    if exist(savefile_fig, 'file')
        if overwrite 
            disp('overwriting existing savefile')
        else
            error(['file ' savefile_fig ' exists already. Either delete or use new parameters!'])
        end
    end
end

% --- input ---
% initializations
output = [];
[r_con1, r_con2, r_con, r_incon, z1, z2, z_comb, p1, p2, p_comb] = deal(zeros(1, length(wsize)));

% ***store loadfile in variable & mod in loop (online 1/0)?
if strcmp(test, 'scrambling')
    tsagg_1 = load(['timeseries' '/' 'ts_scrambled_0' '_dosplit_' num2str(do_split) '_splitby_' split_by '_splittrun_' num2str(split_trun) '_online_'  num2str(1) '.mat']);
    tsagg_2 = load(['timeseries' '/' 'ts_scrambled_1' '_dosplit_' num2str(do_split) '_splitby_' split_by '_splittrun_' num2str(split_trun) '_online_'  num2str(1) '.mat']);
elseif strcmp(test, 'datasets')
    tsagg_1 = load(['timeseries' '/' 'ts_scrambled_0' '_dosplit_' num2str(do_split) '_splitby_' split_by '_splittrun_' num2str(split_trun) '_online_'  num2str(1) '.mat']);
    tsagg_2 = load(['timeseries' '/' 'ts' '_dosplit_' num2str(do_split) '_splitby_' split_by '_splittrun_' num2str(split_trun) '_online_'  num2str(0) '.mat']);
    %tsagg_2 = load(['timeseries' '/' 'ts' '_dosplit_' num2str(do_split) '_splitby_' split_by '_splittrun_' num2str(split_trun) '_online_'  num2str(1) '.mat']);
end
tsagg_1 = tsagg_1.ts;
tsagg_2 = tsagg_2.ts;

% -------------------------------------------------------------------------
% compute tc corr for diff pairings of st and IL/OL -----------------------
% get # timepoints for correlated ts (based on limiting story)
n_tp = length(tsagg_1(1).yagg);

% store mean ts & smooth
[ts_clear, ts_noise] = deal([]);
[ts_clear_sem, ts_noise_sem] = deal([]); %*****TEMP
for b = 1:length(tsagg_1)
    ts_clear(b, :) = smoothdata(nanmean(tsagg_1(b).yagg(:, 1:n_tp)), smoothtype, wsize); 
    ts_noise(b, :) = smoothdata(nanmean(tsagg_2(b).yagg(:, 1:n_tp)), smoothtype, wsize); 

    %*****TEMP STORE VAR 4 PLOTTING - NOTE: var here not great bc want var on smoothed SSs (need 2 plot from analyze_var if doesn't look good)
    ts_clear_sem(b, :) = nanstd(tsagg_1(b).yagg(:, 1:n_tp)) ./ sqrt(sum(~isnan(tsagg_1(b).yagg(:, 1:n_tp)))); 
    ts_noise_sem(b, :) = nanstd(tsagg_2(b).yagg(:, 1:n_tp)) ./ sqrt(sum(~isnan(tsagg_2(b).yagg(:, 1:n_tp)))); 

    %ciobs   = abs(norminv(alpha/2))*semobs; % sem*z-crit for defined alpha 2-tailed
end

% calculate individual congruent correlations
r_con1 = corr(ts_clear(1, :)', ts_noise(1, :)');
r_con2 = corr(ts_clear(2, :)', ts_noise(2, :)');

% ***TEMP r_incon1/2
r_incon1 = corr(ts_clear(1, :)', ts_noise(2, :)');
r_incon2 = corr(ts_clear(2, :)', ts_noise(1, :)');

% average congruent correlations (IL/OL matching st1/2) & incongruent
r_con = mean([corr(ts_clear(1, :)', ts_noise(1, :)') corr(ts_clear(2, :)', ts_noise(2, :)')]);
r_incon = mean([corr(ts_clear(1, :)', ts_noise(2, :)') corr(ts_clear(2, :)', ts_noise(1, :)')]);

n = mean(n_clear, n_noise); % both rcon and incon combine ol/il samples %***should this be mean? or sum???

% --- compute sig of correlation differences ---
szab = sqrt(1/(n - 3) + 1/(n - 3)); %***note: df for both groups avg bc con & incon corr based on IL+OL

% story 1
rz1 = atanh([r_incon r_con1]);
z1 = abs(rz1(2) - rz1(1)) / szab; % high > low
p1 = 1 - normcdf(z1, 0, 1); % 1-tailed

% story 2
rz2 = atanh([r_incon r_con2]);
z2 = abs(rz2(2) - rz2(1)) / szab; % high > low
p2 = 1 - normcdf(z2, 0, 1); % 1-tailed

% combined
rz_comb = atanh([r_incon r_con]);
z_comb = abs(rz_comb(2) - rz_comb(1)) / szab; % high > low
p_comb = 1 - normcdf(z_comb, 0, 1); % 1-tailed

% PLOT IL/OL TIME COURSES OVERLAID ----------------------------------------
colours1 = [46 0 214; 0 135 0]/255; % B/G dark (IL/OL)
colours2 = [206 200 255; 170 228 180]/255; % B/G light

x = tsagg_1(1).x;

figure

% get story vector based on congruent or incongruent
if strcmp(plot_type, 'con1')
    st = [1 1];
    r_txt = r_con1;
elseif strcmp(plot_type, 'con2')
    st = [2 2];
    r_txt = r_con2;
elseif strcmp(plot_type, 'incon1')
    st = [1 2];
    r_txt = r_incon1;
elseif strcmp(plot_type, 'incon2')
    st = [2 1];
    r_txt = r_incon2;
end

% plot error bars
pl(1) = fill([x'/60; flipud(x'/60)], [ts_clear(st(1), :)' - ts_clear_sem(st(1), :)'; flipud(ts_clear(st(1), :)' + ts_clear_sem(st(1), :)')], colours2(1, :), 'linestyle', 'none', 'facealpha', 1); hold on
pl(2) = fill([x'/60; flipud(x'/60)], [ts_noise(st(2), :)' - ts_noise_sem(st(2), :)'; flipud(ts_noise(st(2), :)' + ts_noise_sem(st(2), :)')], colours2(2, :), 'linestyle', 'none', 'facealpha', .6); hold on

% plot obs ts (plot perm first so obs on top)
f(1) = plot(x/60, ts_clear(st(1), :), 'color', colours1(1, :), 'LineWidth', 2); hold on
f(2) = plot(x/60, ts_noise(st(2), :), 'color', colours1(2, :), 'LineWidth', 2); hold on

% figure mods -------------------------------------------------------------
% legend
if strcmp(test, 'scrambling')
    lgd = legend([f(1), f(2)], 'intact', 'scrambled', 'Location', 'northeast', 'FontSize', 14);
elseif strcmp(test, 'datasets')
    lgd = legend([f(1), f(2)], 'lab', 'online', 'Location', 'northeast', 'FontSize', 14);
end

% labels
xlabel('Story time (min)')
ylabel('Normalized RT (s)')

% title
stories = {'Arctic', 'Space'};
title(['corr(' stories{st(1)} ', ' stories{st(2)} ')'], 'fontweight', 'normal')

% add text to state correlation
text(1, 1, ['r = ' num2str(r_txt)], 'FontSize', 16)

% axis props
set(gca, 'xlim', [0 13.5])
set(gca, 'ylim', [-1 1.2])
set(gca, 'xtick', 0:3:12)
set(gca, 'ytick', -1:.5:1)

set(gca, 'fontname', 'Verdana') 
set(gca, 'FontSize', 16)
set(gca, 'linewidth', 2)

% figure props
box off

% legend props
lgd.ItemTokenSize = [17, 12]; % size of lines in legend
legend 
legend boxoff

hold off

% store stats -------------------------------------------------------------
output.r_con1         = r_con1;
output.r_con2         = r_con2;
output.r_con          = r_con;
output.r_incon        = r_incon;
output.z1             = z1;
output.z2             = z2;
output.z_comb         = z_comb;
output.p1             = p1;
output.p2             = p2;
output.p_comb         = p_comb;

% output stats
fprintf('r congruent st 1: %f; r congruent st 2: %f\n', output.r_con1, output.r_con2)
fprintf('r congruent: %f; r incongruent: %f\n', output.r_con, output.r_incon)

fprintf('r congruent 1 > r incongruent: z = %f, p = %f\n', output.z1, output.p1)
fprintf('r congruent 2 > r incongruent: z = %f, p = %f\n', output.z2, output.p2)
fprintf('r congruent comb > r incongruent: z = %f, p = %f\n', output.z_comb, output.p_comb)

% save --------------------------------------------------------------------
if do_save
    print(gcf, savefile_fig, '-dsvg', '-r300')
    save(savefile_stats, '-struct', 'output')
end
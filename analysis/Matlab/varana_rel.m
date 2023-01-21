clc
clear all;
cd('C:\Users\mbain25\Documents\efflisb\Analysis')
addpath(genpath('matlab'))
rng_cfig = rng(257);

% define file params
n_sets     = 50;
smoothtype = 'lowess';
wsize      = 1:100;
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
n_il       = 63;
n_ol       = 112;

plot_type = 'r2'; % 'r' (con & incon overlaid) or 'p' (rcon > rincon); 'con' (both con. curves b4 avg) (for r/p/z: 1, 2, or comb)

colours1 = [106 219 214; 189 178 243]/255; % G/P light (low/high)
colours2 = [109 215 232; 70 240 189]/255; % B/G light (st1/st2)

font = 'Verdana';
fontsize = 16;
linewidth = 3;
%paperposition = [0 0 3 4.5];

do_plot    = 1;

do_save    = 0;
overwrite  = 0;

% -------------------------------------------------------------------------
% savefile
savefile_stats = ['output_relana' '/' 'relana' '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype ...
    '_wsize_' num2str(min(wsize)) num2str(mean(diff(wsize))) num2str(max(wsize)) ...
    '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) ...
    '_splittrun_' num2str(split_trun) '.mat'];
savefile_fig = ['figures_smoothana' '/' 'relana_' plot_type ...
    '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype ...
    '_wsize_' num2str(min(wsize)) num2str(mean(diff(wsize))) num2str(max(wsize)) ...
    '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) ...
    '_splittrun_' num2str(split_trun) '.svg'];

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
[r_con1, r_con2, r_con, r_incon, z1, z2, z_comb, p1, p2, p_comb, dicecoef_con, dicecoef_incon] = deal(zeros(1, length(wsize)));

% skip loop for storing data if stats already saved
if ~exist(savefile_stats, 'file')
    output = load(savefile_stats);
else
    progbar = waitbar(0, 'Please wait...');
    for w = 10%1:length(wsize)
    waitbar(w/length(wsize), progbar, ['wsize ' num2str(w) ' / ' num2str(length(wsize))]);

    % ***store loadfile in variable & mod in loop (online 1/0)?
    tsagg_il = load(['timeseries' '/' 'ts' '_dosplit_' num2str(do_split) ... % ts
        '_splitby_' split_by '_splittrun_' num2str(split_trun) '_online_'  num2str(0) '.mat']);
    tsagg_ol = load(['timeseries' '/' 'ts' '_dosplit_' num2str(do_split) ...
        '_splitby_' split_by '_splittrun_' num2str(split_trun) '_online_'  num2str(1) '.mat']);
    tsagg_il = tsagg_il.ts;
    tsagg_ol = tsagg_ol.ts;

    for s = 1:length(tsagg_il)
        varana(s).inlab = load(['output_varana/smoothing' '/' 'stats_' 'varanalysis' '_story_' num2str(s) ... % varana output
            '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize(w)) ...
            '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) ...
            '_splittrun_' num2str(split_trun) '_online_'  num2str(0) '.mat']);
        varana(s).online = load(['output_varana/smoothing' '/' 'stats_' 'varanalysis' '_story_' num2str(s) ...
            '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(w) ...
            '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) ...
            '_splittrun_' num2str(split_trun) '_online_'  num2str(1) '.mat']);
        tcsig(s).inlab = varana(s).inlab.tcsig_pav < .05;
        tcsig(s).online = varana(s).online.tcsig_pav < .05;
    end

    % -------------------------------------------------------------------------
    % compute tc corr for diff pairings of st and IL/OL -----------------------
    % get # timepoints for correlated ts (based on limiting story)
    n_tp = length(tsagg_il(1).yagg);

    % store mean ts & smooth
    [ts_il, ts_ol] = deal([]);
    [ts_il_sem, ts_ol_sem] = deal([]); %*****TEMP
    for b = 1:length(tsagg_il)
        ts_il(b, :) = smoothdata(nanmean(tsagg_il(b).yagg(:, 1:n_tp)), smoothtype, wsize(w)); 
        ts_ol(b, :) = smoothdata(nanmean(tsagg_ol(b).yagg(:, 1:n_tp)), smoothtype, wsize(w)); 
        
        %*****TEMP STORE VAR 4 PLOTTING - NOTE: var here not great bc want var on smoothed SSs (need 2 plot from analyze_var if doesn't look good)
        ts_il_sem(b, :) = nanstd(tsagg_il(b).yagg(:, 1:n_tp)) ./ sqrt(sum(~isnan(tsagg_il(b).yagg(:, 1:n_tp)))); 
        ts_ol_sem(b, :) = nanstd(tsagg_ol(b).yagg(:, 1:n_tp)) ./ sqrt(sum(~isnan(tsagg_ol(b).yagg(:, 1:n_tp)))); 
        
        %ciobs   = abs(norminv(alpha/2))*semobs; % sem*z-crit for defined alpha 2-tailed
    end
    
    % calculate individual congruent correlations
    r_con1(w) = corr(ts_il(1, :)', ts_ol(1, :)');
    r_con2(w) = corr(ts_il(2, :)', ts_ol(2, :)');
    
    % ***TEMP r_incon1/2
    r_incon1(w) = corr(ts_il(1, :)', ts_ol(2, :)');
    r_incon2(w) = corr(ts_il(2, :)', ts_ol(1, :)');
    
    % average congruent correlations (IL/OL matching st1/2) & incongruent
    r_con(w) = mean([corr(ts_il(1, :)', ts_ol(1, :)') corr(ts_il(2, :)', ts_ol(2, :)')]);
    r_incon(w) = mean([corr(ts_il(1, :)', ts_ol(2, :)') corr(ts_il(2, :)', ts_ol(1, :)')]);
    
    n = mean(n_il, n_ol); % both rcon and incon combine ol/il samples %***should this be mean? or sum???

    % --- compute sig of correlation differences ---
	szab = sqrt(1/(n - 3) + 1/(n - 3)); %***note: df for both groups avg bc con & incon corr based on IL+OL

    % story 1
    rz1 = atanh([r_incon(w) r_con1(w)]);
    z1(w) = abs(rz1(2) - rz1(1)) / szab; % high > low
    p1(w) = 1 - normcdf(z1(w), 0, 1); % 1-tailed
    
    % story 2
    rz2 = atanh([r_incon(w) r_con2(w)]);
    z2(w) = abs(rz2(2) - rz2(1)) / szab; % high > low
    p2(w) = 1 - normcdf(z2(w), 0, 1); % 1-tailed
    
    % combined
    rz_comb = atanh([r_incon(w) r_con(w)]);
    z_comb(w) = abs(rz_comb(2) - rz_comb(1)) / szab; % high > low
    p_comb(w) = 1 - normcdf(z_comb(w), 0, 1); % 1-tailed
    
    % compute overlap btwn sig sections ---------------------------------------
    dicecoef(1).con = (2*nnz(tcsig(1).online == tcsig(1).inlab)) / (2*length(tcsig(1).online)); % limit ts length by st1 ***put in loop, apply n_tp to all & use n_tp for len
    dicecoef(2).con = (2*nnz(tcsig(2).online(1:n_tp) == tcsig(2).inlab(1:n_tp))) / (2*length(tcsig(1).online));
    dicecoef(1).incon = (2*nnz(tcsig(1).online == tcsig(2).inlab(1:n_tp))) / (2*length(tcsig(1).online));
    dicecoef(2).incon = (2*nnz(tcsig(2).online(1:n_tp) == tcsig(1).inlab)) / (2*length(tcsig(1).online));

    dicecoef_con(w) = mean([dicecoef(1).con, dicecoef(2).con]);
    dicecoef_incon(w) = mean([dicecoef(1).incon, dicecoef(2).incon]);

    %%
    %*****PLOT IL/OL TIME COURSES OVERLAID
    if w == 10        
        st = [2 1]; % 1 = Arctic; 2 = Space
        
        colours1 = [46 0 214; 0 135 0]/255; % B/G dark (IL/OL)
        colours2 = [206 200 255; 170 228 180]/255; % B/G light
                
        x = tsagg_il(1).x;
        
        figure
        %plot(1:n_tp, ts_il(2, :)), hold on
        %5plot(1:n_tp, ts_ol(2, :))
        
        % plot error bars
        pl(1) = fill([x'/60; flipud(x'/60)], [ts_il(st(1), :)' - ts_il_sem(st(1), :)'; flipud(ts_il(st(1), :)' + ts_il_sem(st(1), :)')], colours2(1, :), 'linestyle', 'none', 'facealpha', 1); hold on
        pl(2) = fill([x'/60; flipud(x'/60)], [ts_ol(st(2), :)' - ts_ol_sem(st(2), :)'; flipud(ts_ol(st(2), :)' + ts_ol_sem(st(2), :)')], colours2(2, :), 'linestyle', 'none', 'facealpha', .6); hold on

       
        % plot obs ts (plot perm first so obs on top)
        f(1) = plot(x/60, ts_il(st(1), :), 'color', colours1(1, :), 'LineWidth', 2); hold on
        f(2) = plot(x/60, ts_ol(st(2), :), 'color', colours1(2, :), 'LineWidth', 2); hold on
        
        
        % figure mods -------------------------------------------------------------
        % legend
        %lgd = legend([f(1), f(2)], 'Lab', 'Online', 'Location', 'northwest', 'FontSize', 14);
        
        % labels
        xlabel('Story time (min)')
        ylabel('Normalized RT (s)')
        title('corr(Space, Arctic)', 'fontweight', 'normal')

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
        %lgd.ItemTokenSize = [17, 12]; % size of lines in legend
        %legend boxoff

        hold off
        
        savefile_fig = [num2str(st) '.svg'];
        print(gcf, savefile_fig, '-dsvg', '-r300')
        
    end
    
    end
    delete(progbar)
    
    % identify 'dynamic window'
    dynwind = find(z_comb == max(z_comb));
    %dynwind = find(z == max(z));
    
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
    output.dicecoef_con   = dicecoef_con;
    output.dicecoef_incon = dicecoef_incon;
    output.dynwind        = dynwind;
end

% output stats for wsize of 10
fprintf('r congruent st 1 (wsize = 10): %f; r congruent st 2: %f\n', output.r_con1(10), output.r_con2(10))
fprintf('r congruent (wsize = 10): %f; r incongruent: %f\n', output.r_con(10), output.r_incon(10))

fprintf('r congruent 1 > r incongruent (wsize = 10): z = %f, p = %f\n', output.z1(10), output.p1(10))
fprintf('r congruent 2 > r incongruent (wsize = 10): z = %f, p = %f\n', output.z2(10), output.p2(10))
fprintf('r congruent comb > r incongruent (wsize = 10): z = %f, p = %f\n', output.z_comb(10), output.p_comb(10))

fprintf('dynamic window 1: %d samples\n', find(output.z1 == max(output.z1)))
fprintf('dynamic window 2: %d samples\n', find(output.z2 == max(output.z2)))
fprintf('dynamic window comb: %d samples\n', output.dynwind)

fprintf('dicecoef congruent (wsize = 10): %f; dicecoef incongruent: %f\n', output.dicecoef_con(10), output.dicecoef_incon(10))


% -------------------------------------------------------------------------
% plot --------------------------------------------------------------------
if do_plot

% get x vector
%isequal(tsagg_il(1).x(1:n_tp), tsagg_il(2).x(1:n_tp), tsagg_ol(1).x(1:n_tp), tsagg_ol(2).x(1:n_tp)) % verify that all tp arrays equivalent
%x = tsagg_il(1).x(1:n_tp); % all tp arrays equivalent so set x to any
x = 1:length(r_con);

% plot
figure
if any(strcmp(plot_type, {'r_comb', 'r1', 'r2'}))
    
    % story tag
    if strcmp(plot_type(2), '1')
        story_tag = 'Arctic';
    elseif strcmp(plot_type(2), '2')
        story_tag = 'Space';
    end
    
    
    if strcmp(plot_type, 'r_comb')
        f(1) = plot(x, output.r_con, 'color', colours1(1, :), 'LineWidth', linewidth); hold on 
    elseif strcmp(plot_type, 'r1')
        f(1) = plot(x, output.r_con1, 'color', colours1(1, :), 'LineWidth', linewidth); hold on 
    elseif strcmp(plot_type, 'r2')
        f(1) = plot(x, output.r_con2, 'color', colours1(1, :), 'LineWidth', linewidth); hold on 
    end

    f(2) = plot(x, output.r_incon, 'color', colours1(2, :), 'LineWidth', linewidth); hold on

    if strcmp(plot_type(2), '1')
    % legend
    lgd = legend([f(1), f(2)], 'Congruent', 'Incongruent', 'Location', 'southeast');
    
    % legend props
    lgd.ItemTokenSize = [15, 28]; % size of lines in legend
    legend boxoff
    end
    
    ylabel('Correlation ({\itr})')
    title(story_tag, 'fontweight', 'normal');
    
    set(gca, 'xlim', [0 100])
    set(gca, 'ylim', [-.06 .8])
    xticks(0:20:100)
    yticks(0:.1:.8)
    
elseif any(strcmp(plot_type, {'p_comb', 'p1', 'p2'}))
    if strcmp(plot_type, 'p_comb')
        plot(x, output.p_comb, 'color', [1 0 0], 'LineWidth', linewidth); hold on
    elseif strcmp(plot_type, 'p1')
        plot(x, output.p1, 'color', [1 0 0], 'LineWidth', linewidth); hold on
    elseif strcmp(plot_type, 'p2')
        plot(x, output.p2, 'color', [1 0 0], 'LineWidth', linewidth); hold on
    end
    
    % add vertical line at dynamic window
    %line([output.dynwind, output.dynwind], [-.1 .1], 'color', [0 0 0], 'linewidth', 1.5, 'linestyle', '--')
    
    ylabel('{\itp} ({\itr}_{con} > {\itr}_{incon})')
    %%{
    set(gca, 'xlim', [0 100])
    %set(gca, 'ylim', [0 .4])
    xticks(0:20:100)
    %yticks(0:.1:.4)
    %}
    
elseif any(strcmp(plot_type, {'z_comb', 'z1', 'z2'}))
    if strcmp(plot_type, 'z_comb')
        plot(x, output.z_comb, 'color', [1 0 0], 'LineWidth', linewidth); hold on
    elseif strcmp(plot_type, 'z1')
        plot(x, output.z1, 'color', [1 0 0], 'LineWidth', linewidth); hold on
    elseif strcmp(plot_type, 'z2')
        plot(x, output.z2, 'color', [1 0 0], 'LineWidth', linewidth); hold on
    end
    
    % add vertical line at dynamic window
    %line([output.dynwind, output.dynwind], [-.1 .1], 'color', [0 0 0], 'linewidth', 1.5, 'linestyle', '--')
    
    ylabel('{\itz} ({\itr}_{\itcon} > {\itr}_{\itincon})')
    %%{
    set(gca, 'xlim', [0 100])
    set(gca, 'ylim', [0 4])
    xticks(0:20:100)
    yticks(0:1:4)
    %}
    
elseif strcmp(plot_type, 'con')
    f(1) = plot(x, output.r_con1, 'color', colours2(1, :), 'LineWidth', linewidth); hold on
    f(2) = plot(x, output.r_con2, 'color', colours2(2, :), 'LineWidth', linewidth); hold on

    % legend
    lgd = legend([f(1), f(2)], 'Arctic', 'Space', 'Location', 'northwest');
    
    % legend props
    lgd.ItemTokenSize = [15, 28]; % size of lines in legend
    legend boxoff
    
    ylabel('Correlation ({\itr})')
    
    set(gca, 'xlim', [0 100])
    set(gca, 'ylim', [0 .8])
    xticks(0:20:100)
    yticks(0:.1:.8)
end

% --- figure mods ---
% labels
xlabel('Window size (samples)')
%ylabel(plot_type)

% axis props
%set(gca, 'xlim', [0 13.5])
%set(gca, 'xtick', 0:3:12)

set(gca, 'fontname', 'Verdana')
set(gca, 'FontSize', fontsize)
set(gca, 'linewidth', linewidth)

% figure props
axis square
box off
hold off

end

% save --------------------------------------------------------------------
if do_save
    print(gcf, savefile_fig, '-dsvg', '-r300')
    save(savefile_stats, '-struct', 'output')
end
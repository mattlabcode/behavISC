% STEP 4: Plot Perm Distributions and ROC Curves

clc
clear all;
cd('/Users/matthewbain/Documents/Science/Experiments/efflisisc/analysis/')
addpath(genpath('matlab'))
rng_cfig = rng(257);

% define file params
online     = 1;
scrambled  = 1;
story      = 1;
n_sets     = 50;
smoothtype = 'lowess';
wsize      = 10;
do_split   = 0; 

if do_split 
    split_by   = 'enjoyment';               
    split_trun = 0;         
    split      = 1;           
else
    split_by   = '';             
    split_trun = 0;
    split      = 0;           
end

n_perm     = 1; %500

% define params
plot_type  = 'barplot'; % 'permdist' -> perm & obs dists; 'roc_curve' - roc; 'barplot' -> bars cmp mean r for p/o dist; 'corrcause' -> bars cmp r for sig/nsig

if strcmp(plot_type, 'barplot') || strcmp(plot_type, 'corrcause')
    err_type = 'std';     % 'var'; 'std'; 'sem'; 'ci'
    if online
        n = 112;    % # subjects (= # pairs)
    elseif ~online
        n = 63;
    end
end

colours   = [7 162 228; 190 22 34]/255; % B/R dark
colours2  = [137 117 232; 43 204 197]/255; % P/G dark

fontname  = 'Verdana';
fontsize  = 16;
linewidth = 2;

do_save   = 1;

% -------------------------------------------------------------------------
% --- input ---
if strcmp(plot_type, 'corrcause')
    
    % ISC values for truncated ts (sig/nsig sections)
    input_sig = load(['correlations/corr_story_' num2str(story) '_scrambled_' num2str(scrambled) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) '_trunc_sig' '_nperm_' num2str(n_perm) '_online_' num2str(online) '.mat']);
    input_nsig = load(['correlations/corr_story_' num2str(story) '_scrambled_' num2str(scrambled) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) '_trunc_nsig' '_nperm_' num2str(n_perm) '_online_' num2str(online) '.mat']);
else
    input = load(['correlations/corr_story_' num2str(story) '_scrambled_' num2str(scrambled) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) '_trunc_no' '_nperm_' num2str(n_perm) '_online_' num2str(online) '.mat']);
end
    
% savefile
outdir = 'figures_sim';
savefile = [outdir '/' plot_type '_story_' num2str(story) '_scrambled_' num2str(scrambled) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize) '_dosplit_' num2str(do_split) '_splitby_' split_by '_split_' num2str(split) '_splittrun_' num2str(split_trun) '_nperm_' num2str(n_perm) '_online_' num2str(online) '.svg'];
    
% check if savefile already exists
if do_save
    if exist(savefile, 'file')
        choice = questdlg('Warning: Savefile already exists. Overwrite?', 'Warning', 'no', 'yes', 'proceed without saving', 'no');
        switch choice
            case 'no'
                error(['file ' savefile ' exists already. Either delete or use new parameters!'])
            case 'yes'
            case 'proceed without saving'
        end
    end
end

% -------------------------------------------------------------------------
% plot --------------------------------------------------------------------
% add title tag for online/in-lab
if online
    online_tag = 'online';
elseif ~online
    online_tag = 'lab';
end

if story == 1
    story_tag = 'Arctic';
elseif story == 2
    story_tag = 'Space';
end

switch plot_type
    case 'corrcause'
       
        % set data to plot
        x   = [.3 1.7];
        y   = [input_nsig.robs_av input_sig.robs_av]';
        
        switch err_type
            case 'var'
                err = [(input_nsig.robs_var) (input_sig.robs_var)];
            case 'std'
                err = [sqrt(input_nsig.robs_var) sqrt(input_sig.robs_var)];                            
            case 'sem'
                err = [sqrt(input_nsig.robs_var/input_nsig.rperm_n) sqrt(input_sig.robs_var/input_sig.rperm_n)]; 
            case 'ci'
                err = 1.96*[sqrt(input_nsig.robs_var/input_nsig.rperm_n) sqrt(input_sig.robs_var/input_sig.rperm_n)]; 
        end
        
        % plot bars
        for b = 1:2
            bar(x(b), y(b), 'BarWidth', .6, 'FaceColor', colours2(b, :), 'facealpha', .6, 'LineWidth', linewidth/2); hold on               
        end
        
        % plot error bars
        er = errorbar(x, y, err, err, 'LineWidth', linewidth/2); hold on

        % error line props
        er.Color     = [0 0 0];                            
        er.LineStyle = 'none';  
        
        % --- compute sig of correlation difference ---
        rz   = atanh([y(2) y(1)]);
        szab = sqrt(1/(n - 3) + 1/(n - 3));
        z    = abs(rz(2) - rz(1)) / szab;   % sig > nsig
        p    = 1 - normcdf(z, 0, 1);     % 1-tailed
        fprintf('ISC sig sections > ISC nsig sections: z = %f, p = %f\n', z, p)
        
        % plot significance bar
        if p > .05
            sig = nan;
        else
            sig = p;
        end
        sigstar({[x(1), x(2)]}, sig)
        
        % add labels
        set(gca, 'XTickLabel', {'Nonsig.', 'Sig.'})
        ylabel('Median correlation')
        title(story_tag, 'fontweight', 'normal')
        
        % axis props
        xticks(x)
        
        xlim_curr = get(gca, 'XLim');
        %xlim([xlim_curr(1)+.4 xlim_curr(2)-.4])
        
        ylim([-.2 1.2])
        yticks(-.2:.2:1)
        
        xtickangle(45)
        set(gca, 'TickLength', [.005 0]);
        
        set(gca, 'fontname', fontname)
        set(gca, 'FontSize', fontsize)
        set(gca, 'linewidth', linewidth)
        
        % figure props
        set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.5 4.5])
        box off 
        hold off
    
    case 'permdist'
        % get data for plotting (perm and obs distribution)
        nulldist_bincounts = mean(input.bincounts_perm);
        altdist_bincounts  = mean(input.bincounts_obs);
        nulldist_binedges  = input.binedges_perm;
        altdist_binedges   = input.binedges_obs;
        nulldist_ser       = std(input.bincounts_perm)/sqrt(length(input.bincounts_perm));
        altdist_ser        = std(input.bincounts_obs)/sqrt(length(input.bincounts_obs));
               
        % plotting params
        colour_histedges = [20 20 20]/255;
                
        % plot null/alt dists
        figure
        f(1) = histogram('BinEdges', nulldist_binedges, ...
            'BinCounts', nulldist_bincounts, 'FaceColor', colours(1, :), ...
            'EdgeColor', colour_histedges); hold on
        f(2) = histogram('BinEdges', altdist_binedges, ...
            'BinCounts', altdist_bincounts, 'FaceColor', colours(2, :), ...
            'EdgeColor', colour_histedges); hold on
                        
        %{
        % get adjusted binedges so error bars/fitted line centered on midpoints
        nulldist_binedges = nulldist_binedges(1:end - 1) + mode(diff(nulldist_binedges))/2;
        altdist_binedges  = altdist_binedges(1:end - 1) + mode(diff(altdist_binedges))/2;
            
        % plot line outlining normal dist over hists
        plot(nulldist_binedges, nulldist_bincounts, 'color', colours(2, :), 'linewidth', 2); hold on
        plot(altdist_binedges, altdist_bincounts, 'color', colours(1, :), 'linewidth', 2); hold on
                    
        % add error bars to hist bars
        %errorbar(nulldist_binedges , nulldist_bincounts, nulldist_ser, 'LineStyle', 'none', 'Color', colours(2, :)); hold on
        %errorbar(altdist_binedges , altdist_bincounts, altdist_ser, 'LineStyle', 'none', 'Color', colours(1, :));    
        %}
  
        % legend
        lgd = legend([f(1) f(2)], 'Permutation', 'Observed', 'Location', 'northwest');
        
        % legend props
        %set(lgd, 'FontSize', 28);
        lgd.ItemTokenSize = [17, 12]; % size of items in legend
        legend boxoff
        
        % labels
        xlabel('Correlation')
        ylabel('Frequency')        
        title([story_tag ' (' online_tag ')'], 'fontweight', 'normal')
                
        % axis props
        %axis square
        set(gca, 'xlim', [-.6 .64])
        set(gca, 'ylim', [0 .06])
        set(gca, 'XTick', [-.5 0 .5])
        set(gca, 'YTick', 0:.02:.06)
        
        set(gca, 'ticklength', [.005 0])
        
        set(gca, 'fontname', fontname) 
        set(gca, 'FontSize', fontsize)
        set(gca, 'linewidth', linewidth)
        
        % figure props
        box off
        hold off
        
    case 'roc_curve'
        % set data for plotting
        %null_dist = nulldist_bincounts; % using av dist
        %alt_dist  = altdist_bincounts;
        fpr = mean(input.fpr);
        tpr = mean(input.tpr);

        %{
        % plot curve
        auc    = ROC(alt_dist, null_dist, length(input.fpr), 1);
        dprime = sqrt(2)*norminv(auc);
        %}
                
        % plot
        plot(fpr, tpr, 'LineWidth', linewidth + 1, 'color', [1 0 0]); hold on % curve        
        plot([0, 1], [0, 1], '--k', 'LineWidth', linewidth, 'color', [0 0 0]); hold on % straight line w/ slope 1 through origin
        
        % add AUC & d' to plot (AUC doesn't show zero before decimal)
        auc_av = strrep(sprintf('%.3f', input.auc_av), '0.', '.');
        text(.65, .127, {['{\itd}'' = ' sprintf('%.2f', input.dprime_av)], ['AUC = ' auc_av]}, 'fontname', fontname, 'FontSize', fontsize - 2);
        %lgd = legend({['AUC = ' num2str(input.auc_av)], ['dprime = ' num2str(input.dprime_av)]}, 'Location', 'southeast');        
        %legend boxoff
        
        fprintf('d'' = %f (range: %f - %f)\n', input.dprime_av, input.dprime_ran(1), input.dprime_ran(2))
        fprintf('AUC = %f (range: %f - %f)\n', input.auc_av, input.auc_ran(1), input.auc_ran(2))
        
        % labels
        xlabel('False positive rate')
        ylabel('True positive rate')
        
        % axis props
        axis([0 1 0 1])
        axis square
        
        set(gca, 'fontname', fontname) 
        set(gca, 'FontSize', fontsize)
        set(gca, 'linewidth', linewidth)
                
        % figure props
        box off
        hold off
        
    case 'barplot'
        % set data to plot
        x   = [.3 1.7];
        y   = [input.rperm_av input.robs_av]';
        
        switch err_type
            case 'var'
                err = [(input.rperm_var) (input.robs_var)];
            case 'std'
                err = [sqrt(input.rperm_var) sqrt(input.robs_var)];                            
            case 'sem'
                err = [sqrt(input.rperm_var/input.rperm_n) sqrt(input.robs_var/input.robs_n)]; 
            case 'ci'
                err = 1.96*[sqrt(input.rperm_var/input.rperm_n) sqrt(input.robs_var/input.robs_n)]; 
        end
        
        % plot bars
        for b = 1:2
            bar(x(b), y(b), 'BarWidth', .6, 'FaceColor', colours(b, :), 'facealpha', .6, 'LineWidth', linewidth/2); hold on               
        end
        
        % plot error bars
        er = errorbar(x, y, err, err, 'LineWidth', linewidth/2); hold on

        % error line props
        er.Color     = [0 0 0];                            
        er.LineStyle = 'none';  
        
        % add text to state r
        %text(x(1) - bar_width, y(1) + .05, num2str(y(1)), 'FontSize', 17, 'fontname', 'Verdana');
        %text(x(2) - bar_width, y(2) + .05, num2str(y(2)), 'FontSize', 17, 'fontname', 'Verdana');
        
        % display stats
        fprintf('robs = %f (range: %f - %f), p = %f\n', input.robs_av, input.robs_ran(1), input.robs_ran(2), input.prop_av)
        
        % plot significance bar
        p = input.prop_av;
        if p > .05
            sig = nan;
        else
            sig = p;
        end
        sigstar({[x(1), x(2)]}, sig)
        
        % add labels
        set(gca, 'XTickLabel', {'Permutation', 'Observed'})
        set(gca, 'xticklabel', {[]})
        ylabel('Median correlation')
        %title({['Mean Correlation for Observed vs Permutation'], ['Correlation Distribution']})
        
        % axis props
        xticks([])
        ylim([-.2 .6])
        yticks(-.2:.1:.6)
        
        xtickangle(45)
        set(gca, 'TickLength', [.005 0]);
        
        set(gca, 'fontname', fontname)
        set(gca, 'FontSize', fontsize)
        set(gca, 'linewidth', linewidth)
        
        % figure props
        set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 3.5 4.5])
        box off 
        hold off
end

% save --------------------------------------------------------------------
if do_save
    if exist('choice')
        if strcmp(choice, 'proceed without saving')
            disp('file not saved')
        else
            print(gcf, savefile, '-dsvg', '-r300')
        end
    else
        print(gcf, savefile, '-dsvg', '-r300')
    end
end
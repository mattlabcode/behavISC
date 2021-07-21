clc
clear all;
cd('/Users/matthewbain/Documents/Science/Experiments/efflish/analysis/')
addpath(genpath('matlab'))
rng_cfig = rng(257);

% define file params
online = 0;       

% define params
plot_type = 'ratings_bystory';   % 'RTs_bycase'; 'ratings_bystory'; 'RTs_bystory'; 'hists_RTbycase'; 'ratings_bysplit'; 'ratings_bysplit2'

if strcmp(plot_type, 'ratings_bystory')  
    dim        = 'Absorption'; % (*)modify this only              
    
    do_split   = 0;
    split_trun = 0;               
    story      = 0;
elseif strcmp(plot_type, 'RTs_bycase')
    story      = 2;
    
    do_split   = 0;
    dim        = '';
    split_trun = 0;
    
elseif strcmp(plot_type, 'ratings_bysplit') || strcmp(plot_type, 'ratings_bysplit2')
    story      = 2;
    dim        = 'Absorption'; 
    split_trun = 1;           
    
    do_split   = 1;
    
    % read in ts (contains runs in each version pool split)
    ts = load(['timeseries' '/' 'ts' '_dosplit_0' ... 
    '_splitby_' '' '_splittrun_0' '_online_'  num2str(online) '.mat']);
    ts = ts.ts;
    
    colours_split   = [106 219 214; 189 178 243]/255; % G/P light
else
    story      = 0; 
    do_split   = 0;
    dim        = '';
    split_trun = 0;
end

%colours  = [219 115 122; 102 195 235]/255; % R/B light (obs/perm)

if strcmp(plot_type, 'RTs_bycase')
    colours   = [106 219 214; 189 178 243]/255; % G/P light (misc.)
else
    colours  = [109 215 232; 70 240 189]/255; % B/G light (st1/st2)
end
%colours  = [190 22 34; 7 162 228]/255; % R/B dark
%colours  = [149 235 131; 190 75 150]/255; % G/P mid (st1/st2?)
%colours  = [55 0 255; 0 135 0]/255; % G/B dark
%colours  = [55 0 250; 0 135 0]/255; % G/B dark

font = 'Verdana';
fontsize = 16;
linewidth = 2;
paperposition = [0 0 3 4.5];

plot_spread = 1;

do_save   = 0;
overwrite = 0;

% -------------------------------------------------------------------------
% input
dxav = load(['Group_Processing' '/' 'av_processed' '_logtrans_0' '_online_' num2str(online) '.mat']); % group-level metrics
dxav = dxav.dxav; 

% savefile
savefile = ['figures_dx' '/' plot_type '_story_' num2str(story) ...
    '_dosplit_' num2str(do_split) '_splitby_' lower(dim) '_splittrun_' num2str(split_trun) ...
    '_online_' num2str(online) '.svg'];

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
% add title tag for online/in-lab
if online
    online_tag = 'online';
elseif ~online
    online_tag = 'lab';
end

% add story tag for splits & letter-case plots
if do_split || strcmp(plot_type, 'RTs_bycase')
    if story == 1
        story_tag = 'Arctic';
    elseif story == 2
        story_tag = 'Space';
    end
end

switch plot_type
	% ---------------------------------------------------------------------
    case 'comprehension_bystory'
        % set up plotting structure
        w     = [0.3 0.8];
        sides = [0 1];
        c     = {[0 0 0; colours(1, :); colours(1, :)], [0 0 0; colours(2, :); colours(2, :)]};
        
        % initializations
        xvals      = [];
        XTick      = [];
        XTickLabel = [];
        
        % get data
        for b = 1:length(dxav)     
            x(b, :) = dxav(b).comp_score/10;
            
            % replace rejected runs w/ nan
            x(b, find(dxav(b).rejected)) = nan;
        end
        
        % reject (remove) all subj w/ at least 1 rejected run
        x(:, any(isnan(x), 1)) = [];
        
        % --- plot ---
        figure
        for b = 1:length(dxav)
            xvals{b} = plot_vboxplotPlus(x(b, :), w, c{b}, sides(b), plot_spread); hold on;

            % store tick positions
            XTick(b) = mean(w);

            % shift position of next boxplot
            if plot_spread
                w = w + .5;
            else
                w = w + .5;
            end
        end
        
        % connect corresponding points for each participant
        if plot_spread
            for s = 1:length(x)
                plot([xvals{1}(s) xvals{2}(s)], [x(1, s), x(2, s)], 'color', [0 0 0])
            end
        end

        % labels
        set(gca, 'XTick', XTick, 'XTickLabel', {'Arctic', 'Space'})
        ylabel('Correct answers (/10)'); 
        %title(['Comprehension' ' (' online_tag ')'], 'fontweight', 'normal');
        online_tag(1) = upper(online_tag(1));
        title(online_tag, 'fontweight', 'normal');
        
        % axis props
        set(gca, 'xlim', [0.2 1.4]);
        set(gca, 'ylim', [5 10]);
        yticks(5:1:10)
        
        xtickangle(45)
        
        set(gca, 'fontname', font) 
        set(gca, 'FontSize', fontsize)
        set(gca, 'linewidth', linewidth)

        hold off        
    % ---------------------------------------------------------------------
    case 'ratings_bystory'
        % set up plotting structure
        w     = [0.3 0.8];
        sides = [0 1];
        c     = {[0 0 0; colours(1, :); colours(1, :)], [0 0 0; colours(2, :); colours(2, :)]};

        % initializations
        xvals      = [];
        XTick      = [];
        XTickLabel = [];
        
        % get data
        for b = 1:length(dxav)     
            x(b, :) = dxav(b).(lower(dim));
            
            % replace rejected runs w/ nan
            x(b, find(dxav(b).rejected)) = nan;
        end
        
        % reject (remove) all subj w/ at least 1 rejected run
        x(:, any(isnan(x), 1)) = [];
        
        % --- plot ---
        figure
        for b = 1:length(dxav)
            xvals{b} = plot_vboxplotPlus(x(b, :), w, c{b}, sides(b), plot_spread); hold on;

            % store plot label parameters
            XTick(b) = mean(w);

            % shift position of next boxplot
            if plot_spread
                w = w + .5;
            else
                w = w + .5;
            end
        end
        
        % --- compare mean ratings story 1 v 2 (paired) ---
        [h, p, ~, stats] = ttest(x(1, :), x(2, :), 'tail', 'both');
        stats.fxsize = abs(mean(x(1, :)) - mean(x(2, :))); 
        fprintf('%s story 1 vs story 2: t(%d) = %f, p = %f\n', lower(dim), stats.df, stats.tstat, p)
                
        % also store descriptive stats & effect size
        for b = 1:length(dxav)
            stats.desc(b).data = x(b, :);
            stats.desc(b).mean = mean(x(b, :));
            stats.desc(b).var = var(x(b, :));
        end
        
        % plot significance bar
        if p > .05
            sig = nan;
        else
            sig = p;
        end
        sigstar({[XTick(1), XTick(2)]}, sig)

        % connect corresponding points for each participant
        if plot_spread
            for s = 1:length(x)
                plot([xvals{1}(s) xvals{2}(s)], [x(1, s), x(2, s)], 'color', [0 0 0])
            end
        end

        % labels
        set(gca, 'XTick', XTick, 'XTickLabel', {'Arctic', 'Space'})
        ylabel('Rating (1-7)');
        % title([dim ' (' online_tag ')'], 'fontweight', 'normal');
        title(dim, 'fontweight', 'normal');
        
        % axis props
        xtickangle(45)
        set(gca, 'xlim', [0.2 1.4]);
        set(gca, 'ylim', [1 7.5]);
        set(gca, 'YTick', 1:1:7)
        
        set(gca, 'fontname', font) 
        set(gca, 'fontsize', fontsize)
        set(gca, 'linewidth', linewidth)
        
        hold off
        
    % ---------------------------------------------------------------------
	% --- plot mean RT by story ---
    case 'RTs_bystory'        
        % set up plotting structure
        w     = [0.3 0.8];
        sides = [0 1];
        c     = {[0 0 0; colours(1, :); colours(1, :)], [0 0 0; colours(2, :); colours(2, :)]};
        
        % initializations
        xvals      = [];
        XTick      = [];
        XTickLabel = [];
        
        % get data
        for b = 1:length(dxav)    
            x(b, :) = dxav(b).av_RT;
            
            % replace rejected runs w/ nan
            x(b, find(dxav(b).rejected)) = nan;
        end
        
        % reject (remove) all subj w/ at least 1 rejected run
        x(:, any(isnan(x), 1)) = [];
        
        % --- plot ---
        figure
        for b = 1:length(dxav)
            xvals{b} = plot_vboxplotPlus(x(b, :), w, c{b}, sides(b)); hold on;

            % store plot label parameters to add after loop
            XTick(b) = mean(w);

            % shift position of next boxplot
            w = w + 1;        
        end

        % connect corresponding points for each participant
        for s = 1:length(x)
            plot([xvals{1}(s) xvals{2}(s)], [x(1, s) x(2, s)], 'color', [0 0 0])
        end

        % --- compare mean RTs story 1 v story 2(paired) ---
        [h, p, ~, stats] = ttest(x(1, :), x(2, :), 'tail', 'both');
        stats.fxsize = abs(mean(x(1, :)) - mean(x(2, :))); 
        fprintf('mean RT story 1 v story 2: t(%d) = %f, p = %f\n', stats.df, stats.tstat, p)
        
        % also store descriptive stats & effect size
        for b = 1:length(dxav)
            stats.desc(b).data = x(b, :);
            stats.desc(b).mean = mean(x(b, :));
            stats.desc(b).var = var(x(b, :));
        end
        
        % plot significance bar
        if p > .05
            sig = nan;
        else
            sig = p;
        end
        sigstar({[XTick(1), XTick(2)]}, sig)
        
        % labels
        set(gca, 'XTick', XTick, 'XTickLabel', {'Arctic', 'Space'})
        ylabel('Average RT (s)');
        title(['Response time' ' (' online_tag ')'], 'fontweight', 'normal');
        
        % axis props
        set(gca, 'xlim', [0.2 2]); 
        set(gca, 'ylim', [0 1.8]); 
        set(gca, 'YTick', 0:.4:1.6)
        %set(gca, 'ylim', [0 1.5]); 
        
        xtickangle(45)
        
        set(gca, 'fontname', font) 
        set(gca, 'fontsize', fontsize)
        set(gca, 'linewidth', linewidth)
        
        hold off        

    % ---------------------------------------------------------------------
	% plot mean RT by case for given story --------------------------------
    case 'RTs_bycase'
        % set up plotting structure
        w     = [0.3 0.8];
        sides = [0 1];
        c     = {[0 0 0; colours(1, :); colours(1, :)], [0 0 0; colours(2, :); colours(2, :)]};
                
        plot_spread = 0;
        
        % initializations
        xvals      = [];
        XTick      = [];
        XTickLabel = [];
        
        % get data
        x(1, :) = dxav(story).avRT_low;
        x(2, :) = dxav(story).avRT_upp;
        x(:, find(dxav(story).rejected)) = []; % reject subjects
        
        % --- plot ---
        figure
        for b = 1:length(dxav)
            xvals{b} = plot_vboxplotPlus(x(b, :), w, c{b}, sides(b), plot_spread); hold on;

            % store plot label parameters to add after loop
            XTick(b) = mean(w);

            % shift position of next boxplot
            if plot_spread
                w = w + 1;
            else
                w = w + .5;
            end
        end

        % connect corresponding points for each participant
        if plot_spread
            for s = 1:length(x)
                plot([xvals{1}(s) xvals{2}(s)], [x(1, s) x(2, s)], 'color', [0 0 0])
            end
        end
        
        % --- compare mean RTs upper- v lower-case (paired) ---
        [h, p, ~, stats] = ttest(x(1, :), x(2, :), 'tail', 'right');
        stats.fxsize = mean(x(1, :)) - mean(x(2, :)); 
        fprintf('story %d RTs lower > upper-case: t(%d) = %f, p = %f\n', story, stats.df, stats.tstat, p)
        
        % also store descriptive stats & effect size
        for b = 1:length(dxav)
            stats.desc(b).data = x(b, :);
            stats.desc(b).mean = mean(x(b, :));
            stats.desc(b).var = var(x(b, :));
        end
        
        % plot significance bar
        if p > .05
            sig = nan;
        else
            sig = p;
        end
        sigstar({[XTick(1), XTick(2)]}, sig)
        
        % labels
        set(gca, 'XTick', XTick, 'XTickLabel', {'Lower', 'Upper'})
        ylabel('Average RT (s)');
        %title([story_tag ' (' online_tag ')'], 'fontweight', 'normal');
        title(story_tag, 'fontweight', 'normal');
        
        % axis props
        set(gca, 'xlim', [0.2 (XTick(2) + .45)]); 
        %set(gca, 'xlim', [0.2 2]); 
        set(gca, 'ylim', [0 1.8]); 
        set(gca, 'YTick', 0:.3:1.8)
        
        xtickangle(45)
        
        set(gca, 'fontname', font) 
        set(gca, 'fontsize', fontsize)
        set(gca, 'linewidth', linewidth)
        
        hold off        

    % ---------------------------------------------------------------------
	% plot dists of avg RTs for each story by case ------------------------
    case 'hists_RTbycase'
        for s = 1:length(dxav)
            % lower-case plot
            subplot(2, 2, s + (s - 1))
            hist([dxav(s).avRT_low])

            % labels
            xlabel('Response Time (s)'); 
            ylabel('Frequency');
            title(['Lower-case RTs: Story ' num2str(s)], 'FontWeight', 'normal')

            % axis props
            box off
            set(gca, 'fontsize', 9)
            set(gca, 'fontname', 'Verdana')
            set(gca, 'linewidth', 1)
            
            % upper-case plot
            subplot(2, 2, s*2)
            hist([dxav(s).avRT_upp])

            % labels
            xlabel('Mean Response Time (s)'); 
            ylabel('Frequency');
            title(['Upper-case RTs: Story ' num2str(s)], 'FontWeight', 'normal')

            % axis props
            box off
            set(gca, 'fontsize', 9)
            set(gca, 'fontname', 'Verdana')
            set(gca, 'linewidth', 1)
            
            % figure props
            colormap([.5 .5 .5])
            
            % --- compare mean RTs upper- v lower-case (paired) ---
            [h, p, ~, stats] = ttest(dxav(s).avRT_low, dxav(s).avRT_upp, 'tail', 'both');
            fprintf('story %d RTs lower > upper-case: t(%d) = %f, p = %f\n', s, stats.df, stats.tstat, p)     
        end 

        hold off
        
    % ---------------------------------------------------------------------
    case 'ratings_bysplit' 
        % --- split behavioural responses by version for comparison ---
        % initialize arrays to store resp for each split
        for b = 1:length(dxav)
            [split(b).RTs, split(b).resp] = deal([]);
        end
        
        % loop through all versions
        for v = 1:size(ts(story).ypool, 2)

            % get fids for current version pool
            curr_fids = ts(story).run(:, v);
            curr_fids = curr_fids(curr_fids ~= 0); % remove zeros

            % if not even # of fids, remove final
            if mod(length(curr_fids), 2) ~= 0
                curr_fids = curr_fids(1:end - 1);
            end

            % get indices of sorted responses for relevant dimension
            curr_scores = dxav(story).(lower(dim))(curr_fids);
            [~, sort_ix] = sort(curr_scores);
            
            for b = 1:1:length(dxav)
                % get run indices of each split
                if b == 1
                    split(b).ix = curr_fids(sort_ix(1:length(sort_ix)/2 - split_trun)); % ***put split_ix into single array (start loop here)
                elseif b == 2
                    split(b).ix = curr_fids(sort_ix(length(sort_ix)/2 + 1 + split_trun:end));
                end
                
                % concatenate data for current version split runs onto split data array
                split(b).RTs = [split(b).RTs, dxav(story).av_RT(split(b).ix)]; % ***just a curiosity (Bjorn idea; probably remove)
                split(b).resp = [split(b).resp, dxav(story).(lower(dim))(split(b).ix)];
            end
        end

        % --- plot ratings by split ---
        figure
        % plotting params
        w     = [0.3 0.8]; % x axis position of boxes 
        c     = {[0 0 0; colours_split(1, :); colours_split(1, :)], [0 0 0; colours_split(2, :); colours_split(2, :)]};
        sides = [0 1];

        % plotting initialization
        XTickLabel = [];

        % plot
        for b = 1:length(dxav) 

            % data for plotting
            x(b, :) = split(b).resp;
            
            % plot
            plot_vboxplotPlus(x(b, :), w, c{b}, sides(1)); hold on;

            % store plot label parameters
            XTick(b) = mean(w);

            % shift position of next boxplot %***need this??
            w = w + .5;        
        end

        % --- compare split RTs/ratings ---
        % compare mean RTs high > low (independent)
        [~, p, ~, stats] = ttest2(split(2).RTs, split(1).RTs, 'tail', 'right');
        fprintf('story %d RTs high vs low splits: t(%d) = %f, p = %f\n', story, stats.df, stats.tstat, p)

        % compare mean ratings high > low (independent)
        [~, p, ~, stats] = ttest2(split(2).resp, split(1).resp, 'tail', 'right');
        stats.fxsize = mean(split(2).resp) - mean(split(1).resp); 
        fprintf('story %d ratings high vs low splits: t(%d) = %f, p = %f; effect size: %f\n', story, stats.df, stats.tstat, p, stats.fxsize)

        % also store descriptive stats & effect size
        for b = 1:length(dxav)
            stats.desc(b).data = split(b).resp;
            stats.desc(b).mean = mean(split(b).resp);
            stats.desc(b).var = var(split(b).resp);
        end

        % plot significance bar
        if p > .05
            sig = nan;
        else
            sig = p;
        end
        sigstar({[XTick(1), XTick(2)]}, sig)

        % labels
        set(gca, 'XTickLabel', {'Low', 'High'})
        ylabel('Rating (1-7)');
        title([story_tag ': ' dim ' (' online_tag ')'], 'FontWeight', 'Normal');

        % axis props
        xtickangle(45)
        set(gca, 'xlim', [0.2 1.4]);
        set(gca, 'ylim', [1 7.5])
        set(gca, 'XTick', XTick)
        set(gca, 'YTick', 1:1:7)
        
        set(gca, 'fontname', font)
        set(gca, 'fontsize', fontsize) 
        set(gca, 'linewidth', linewidth)

        hold off
    
    % ---------------------------------------------------------------------
    case 'ratings_bysplit2' 
        % --- split behavioural responses by version for comparison ---
        % initialize arrays to store resp for each split
        for b = 1:length(dxav)
            [split(b).RTs, split(b).resp] = deal([]);
        end
        
        % loop through all versions
        for v = 1:size(ts(story).ypool, 2)

            % get fids for current version pool
            curr_fids = ts(story).run(:, v);
            curr_fids = curr_fids(curr_fids ~= 0); % remove zeros

            % if not even # of fids, remove final
            if mod(length(curr_fids), 2) ~= 0
                curr_fids = curr_fids(1:end - 1);
            end

            % get indices of sorted responses for relevant dimension
            curr_scores = dxav(story).(lower(dim))(curr_fids);
            [~, sort_ix] = sort(curr_scores);
            
            for b = 1:1:length(dxav)
                % get run indices of each split
                if b == 1
                    split(b).ix = curr_fids(sort_ix(1:length(sort_ix)/2 - split_trun)); % ***put split_ix into single array (start loop here)
                elseif b == 2
                    split(b).ix = curr_fids(sort_ix(length(sort_ix)/2 + 1 + split_trun:end));
                end
                
                % concatenate data for current version split runs onto split data array
                split(b).RTs = [split(b).RTs, dxav(story).av_RT(split(b).ix)]; % ***just a curiosity (Bjorn idea; probably remove)
                split(b).resp = [split(b).resp, dxav(story).(lower(dim))(split(b).ix)];
            end
        end
        
        % --- plot boxplots ---
        % get data for plotting
        data = [];
        for b = 1:2
            data(:, b) = split(b).resp;
        end
        
        figure
        x = 1:size(data, 2);
        boxplot(data, [1 2], 'widths', .4, 'symbol', '', 'Whisker', 3); hold on
        
        h = findobj(gca, 'Tag', 'Box');

        % add colours to boxplots ***reversed for some reason??
        patch(h(1).XData, h(1).YData, colours_split(2, :), 'FaceAlpha', .5);
        patch(h(2).XData, h(2).YData, colours_split(1, :), 'FaceAlpha', .5);

        % define colours for plotting scatter
        colours1 = [];
        colours_split = [.85 .85 .85; .85 .85 .85];
        for j = 1:length(x)
            colours1 = vertcat(colours1, repmat(colours_split(j, :), length(data), 1));
        end

        % plot scatter
        x1 = repmat(x, length(data), 1); 
        scatter(x1(:), data(:), 20, colours1, 'filled', 'MarkerFaceAlpha', .6, ...
            'jitter', 'on', 'jitterAmount', 0.1, 'MarkerEdgeColor', [.4 .4 .4], 'LineWidth', .3);

        % --- compare split RTs/ratings ---
        % compare mean RTs high > low (independent)
        [~, p, ~, stats] = ttest2(split(2).RTs, split(1).RTs, 'tail', 'right');
        fprintf('story %d RTs high vs low splits: t(%d) = %f, p = %f\n', story, stats.df, stats.tstat, p)

        % compare mean ratings high > low (independent)
        [~, p, ~, stats] = ttest2(split(2).resp, split(1).resp, 'tail', 'right');
        stats.fxsize = mean(split(2).resp) - mean(split(1).resp); 
        fprintf('story %d ratings high vs low splits: t(%d) = %f, p = %f; effect size: %f\n', story, stats.df, stats.tstat, p, stats.fxsize)

        % also store descriptive stats & effect size
        for b = 1:length(dxav)
            stats.desc(b).data = split(b).resp;
            stats.desc(b).mean = mean(split(b).resp);
            stats.desc(b).var = var(split(b).resp);
        end

        % plot significance bar
        if p > .05
            sig = nan;
        else
            sig = p;
        end
        sigstar([1 2], sig)

        % labels
        set(gca,'XTickLabel',{'Low', 'High'})
        ylabel('Rating (1-7)');
        title([story_tag ': ' dim ' (' online_tag ')'], 'FontWeight','Normal');

        % figure props
        set(findobj(gcf, 'LineStyle', '--'), 'LineStyle', '-')
        set(findobj(gca, 'type', 'line'), 'linew', 2, 'color', [0 0 0])

        % axis props
        xtickangle(0)
        set(gca, 'xlim', [0.6 2.4]);
        set(gca, 'ylim', [1 7.5])
        yticks(1:1:7)
        
        set(gca, 'fontname', font)
        set(gca, 'fontsize', fontsize) 
        set(gca, 'linewidth', linewidth)
        
        box off
        
        paperposition = [0 0 2.5 4.5]*.9;

        hold off
end

% save --------------------------------------------------------------------
if do_save
    % set and save new figure size
    set(gcf,'PaperUnits', 'inches', 'PaperPosition', paperposition)
    %set(gcf,'PaperUnits', 'inches', 'PaperPosition', [0 0 5 7])
    print(gcf, savefile, '-dsvg', '-r300')    
end
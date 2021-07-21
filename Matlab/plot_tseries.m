function plot_tseries(groupdir, norm, span, story)
%% -----------------------------------------------------------------------
% PARAMS 
n_blocks = 2;
sf       = 2; % trial sampling rate

%% ------------------------------------------------------------------------
% LOAD REQUESTED DATA 
% define file based on specified parameters
file_ix = strfind({groupdir.name}, ['ts_' 'norm_' num2str(norm) '_smoothed_' num2str(span)]);
file = find(cellfun(@nnz, file_ix));

% error message is data not stored w/ given parameters 
if isempty(file)
    disp('Error: The requested file does not exist.')
    return
else
    % load corresponding data 
	file = groupdir(file).name;
    data = load(file);
    ts = data.ts;
end

%% ------------------------------------------------------------------------
% SET PLOTTING PARAMETERS
if norm
    % set plotting parameters accordingly
    ylim = [-.3 .3];
    ylab = 'Mean RT (norm)';
else
    ylim = [0 1.5];
    ylab = 'Mean RT (s)';
end 

%% ------------------------------------------------------------------------
% PLOT MEAN AGGREGATE TS
savefile = ['Plotting' filesep 'ts_avagg_' 'norm_' num2str(norm) '_smoothed_' num2str(span) '.png'];
savefile_fft = ['Plotting' filesep 'fft_avagg' 'norm_' num2str(norm) '_smoothed_' num2str(span) '.png'];

figure
for b = 1:n_blocks
    subplot(2, 1, b);

    % calculate stats on SSs
    y   = nanmean(ts(b).ystacked); 
    std = nanstd(ts(b).ystacked);
    n   = sum(~isnan(ts(b).ystacked));
    
    % with shaded error bars
    %fig = shadedErrorBar(ts(b).x/60, y, ts(b).err);
	fig = shadedErrorBar(ts(b).x/60, y, std./sqrt(n));
    
    % labels
    title(['Story ' num2str(b)], 'FontSize', 12);
    ylabel(ylab, 'FontSize', 12);
    
    % axis parameters
    set(gca, 'ticklength', [0.001 0.001])
    set(gca, 'xlim', [0 13.7]);
    set(gca, 'ylim', ylim);
    
    box off
        
    % add sine wave to plot with period of smoothing window*****
    %{
    hold on   
 
    t = ts(b).x(1)/60:0.01:ts(b).x(end)/60;
    a = .075*sin(1/(span*SR)*60*2*pi*t);
    
    plot(t, a, 'color', 'r', 'LineWidth', .3);
	hold off
    %}
end

% plot labels
xlabel('Story Time (min)', 'FontSize', 12);
suptitle('Mean Aggregate Case-Judgement Response Times');

print(gcf, savefile, '-dpng', '-r300')

% plot power spectra*****
%{
figure
for b = 1:n_blocks
    subplot(2, 1, b);

    plot_fft(nanmean(ts(b).ystacked), sf); 
    
    title(['Mean Aggregated RTs: Power Spectrum: Story ' num2str(b)])
    xlabel('Frequency')
    ylabel('Power')
end

print(gcf, savefile_fft, '-dpng', '-r300')
%}

%% ------------------------------------------------------------------------
% PLOT AVG TS FOR EACH VERSION W/ STEMS TO DEPICT SAMPLING
%*****
%{
% define which versions to plot
n_versions = size(ts(1).ypool, 2);
versions = 1:1;

savefile = ['Plotting' filesep 'ts_versionsamp_' 'norm_' num2str(norm) '_smoothed_' num2str(span) '_versions_' num2str(versions(1)) '-' num2str(versions(end)) '.png'];

x = ts(1).x;

figure
for v = versions
    % get indices of all trials for current version
    y = nanmean(cell2mat(ts(1).ypool(:, v)));
    y(~isnan(y)) = 1;
    
    stem(x(1:40), y(1:40), 'linewidth', 1.5)
    hold on
end

set(gca, 'xlim', [x(1) 150]);
set(gca, 'ylim', [0 2], 'ytick', []);

% set labels
legend({'v1', 'v2', 'v3', 'v4', 'v5', 'v6'})
title('RT Sampling by Timing Version', 'fontsize', 14) 
xlabel('Story Time (s)')
ylabel('Target')

box off

hold off

print(gcf, savefile, '-dpng', '-r300')
%}
%% ------------------------------------------------------------------------
% PLOT TS OVERLAID FOR EACH POOL
%*****
%{
savefile = ['Plotting' filesep 'ts_poolall_' 'norm_' num2str(norm) '_smoothed_' num2str(span) '_story_' num2str(story) '.png'];

x = ts(story).x/60;

figure
for v = 1:size(ts(story).ypool, 2)
    subplot(3, 2, v)
    
    % loop through ts in each pool
    for i = 1:size(cell2mat(ts(story).ypool(:, v)), 1)
        
        % define data and smooth for plotting
        y = cell2mat(ts(story).ypool(i, v));
        y = smooth(y, 10);
        %y(isnan(y)) = interp1(x(~isnan(y)), y(~isnan(y)), x(isnan(y)));
        
        plot(x, y, 'linewidth', 1)

        hold on
    end
    
	set(gca, 'ylim', [-.4 .4]);
    
    % add title
    title(['Version ' num2str(v)]) 

    % add labels
    if v == 3 || v == 4
        ylabel('Mean Response Time (normalized)')
    end
    
    if v == 5 || v == 6
        xlabel('Story Time (min)')
    end
    
    %legend({'s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9', 's10'}, 'FontSize', 5)
end

% set title
suptitle(['Pooled Time Series by Timing Version: Story ' num2str(story)] ) 

hold off

print(gcf, savefile, '-dpng', '-r300')
%}
%% ------------------------------------------------------------------------
% PLOT ALL SSs OVERLAID 
%{
savefile = ['Plotting' filesep 'ts_aggall_' 'norm_' num2str(norm) '_smoothed_' num2str(span) '_story_' num2str(story) '.png'];

x = ts(story).x/60;

figure
for s = 1:size(ts(story).ystacked, 1)
    plot(x, smoothts(ts(story).ystacked(s, :), 'b', 40), 'linewidth', 1.2)
	%plot(x, ts(story).ystacked(s, :), 'linewidth', 1.2)
    hold on
end

set(gca, 'ylim', ylim);

% set labels
legend({'ss1', 'ss2', 'ss3', 'ss4', 'ss5', 'ss6', 'ss7', 'ss8', 'ss9', 'ss10'})
title(['Aggregate Response Times: Story ' num2str(story)], 'fontsize', 14) 
xlabel('Story Time (min)')
ylabel('Mean Response Time (normalized)')

box off

hold off

print(gcf, savefile, '-dpng', '-r300')
%}
end
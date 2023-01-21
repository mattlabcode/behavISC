clc
clear all;
cd('/Users/matthewbain/Documents/Science/Experiments/efflisisc/analysis/')
addpath(genpath('matlab'))
rng_cfig = rng(257);

% define file params
online     = 1;
scrambled  = 0;
n_sets     = 50;
smoothtype = 'lowess';
n_perm     = 1; %***this has to be 1, right? For efficiency?

% define params
p_crit    = .05;
ver_ofset = .005;   % amount to separate sig stars from points in plots
plot_type = 'corr'; % 'corr' -> wsize*r; 'propsig' -> wsize*%sig; 'z' -> wsize*z
wsize     = 1:1:100;

%colours   = [7 162 228; 190 22 34]/255; % B/R dark
colours  = ([109 215 232; 70 240 189] - 50)/255; % B/G light (st1/st2)

font = 'Verdana';
fontsize = 16;
linewidth = 3;

do_save   = 1;
overwrite = 1;

% -------------------------------------------------------------------------
% savefile
savefile_stats = ['output_smoothana' '/' 'smoothana_' plot_type '_scrambled_' num2str(scrambled) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(min(wsize)) num2str(mean(diff(wsize))) num2str(max(wsize)) '_nperm_' num2str(n_perm) '_online_' num2str(online) '.mat'];
savefile_fig = ['figures_smoothana' '/' 'smoothana_' plot_type '_scrambled_' num2str(scrambled) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(min(wsize)) num2str(mean(diff(wsize))) num2str(max(wsize)) '_nperm_' num2str(n_perm) '_online_' num2str(online) '.svg'];

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
output = [];
[r, p, z, propsig] = deal(zeros([2, length(wsize)]));   % matrices to store input

% skip loop if stats already saved
if exist(savefile_stats, 'file')
    output = load(savefile_stats);
else
    for s = 1:2
        for w = 1:length(wsize)
            if strcmp(plot_type, 'corr') || strcmp(plot_type, 'z')
                input = load(['correlations/smoothing' '/' 'corr' '_story_' num2str(s) '_scrambled_' num2str(scrambled) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize(w)) '_dosplit_' num2str(0) '_splitby_' '' '_split_' num2str(0) '_splittrun_' num2str(0) '_trunc_no' '_nperm_' num2str(n_perm) '_online_' num2str(online) '.mat']);

                r(s, w) = input.robs_av;
                p(s, w) = input.prop_av; 
                z(s, w) = input.zp_av;
                
                % store stats
                output.r = r;
                output.p = p;
                output.z = z;

            elseif strcmp(plot_type, 'propsig') 
                input = load(['output_varana/smoothing' '/' 'stats_varanalysis' '_story_' num2str(s) '_scrambled_' num2str(scrambled) '_nsets_' num2str(n_sets) '_smoothtype_' smoothtype '_wsize_' num2str(wsize(w)) '_dosplit_' num2str(0) '_splitby_' '' '_split_' num2str(0) '_splittrun_' num2str(0) '_online_' num2str(online) '.mat']);

                propsig(s, w) = input.propsig; 
                
                % store stats
                output.propsig = propsig;
            end
        end
    end
end

% -------------------------------------------------------------------------
% plot --------------------------------------------------------------------
% add title tag for online/in-lab
if online
    online_tag = 'Online';
elseif ~online
    online_tag = 'Lab';
end

if strcmp(plot_type, 'corr')
    
    % get data
    r = output.r;
	p = output.p;
    
    % plot
    figure
    for s = 1:2
        f(s) = plot(wsize, r(s, :), 'color', colours(s, :), 'LineWidth', linewidth); hold on

        %{
        % add sig stars over r values
        sig_ix = p(s, :) < p_crit;
        x_sig  = wsize(sig_ix);
        y_sig  = r(s, sig_ix);
        text(x_sig, y_sig + ver_ofset, "*"), hold on; 
        %}
    end

    % labels
    ylabel('Median observed ISC')

    % axis props
    ylim([0 .8])
    yticks([0:.1:.8]);

elseif strcmp(plot_type, 'propsig')
    
    % get data
    propsig = output.propsig;
    
    % plot
    figure
    for s = 1:2
        f(s) = plot(wsize, propsig(s, :), 'color', colours(s, :), 'LineWidth', linewidth); hold on
    end

    % labels
    ylabel('Percent significant')

    % axis props
    ylim([0 100])
    yticks(0:20:100);

elseif strcmp(plot_type, 'z')
    z = output.z;
    
    % plot
    figure
    for s = 1:2
        f(s) = plot(wsize, z(s, :), 'color', colours(s, :), 'LineWidth', linewidth); hold on
    end

    % labels
    ylabel('z-score')

    % axis props
    %yticks(.1:.1:.8);
end

% legend
lgd = legend([f(1) f(2)], 'Arctic', 'Space', 'Location', 'northwest');

% legend props
%set(lgd, 'FontSize', 20, 'LineWidth', 4); %***this doesn't actually adjust size?!!!
lgd.ItemTokenSize = [15, 28]; % size of lines in legend
legend boxoff

% labels
xlabel('Window size (samples)')
title(online_tag, 'fontweight', 'normal')

% axis props
xticks(0:20:wsize(end));
%set(gca, 'TickLength', [0.005 0]);

set(gca,'fontname', font)
set(gca, 'FontSize', fontsize)
set(gca, 'linewidth', linewidth) 

% figure props
axis square
box off 

if do_save
	save(savefile_stats, '-struct', 'output')
    print(gcf, savefile_fig, '-dsvg', '-r300')
end
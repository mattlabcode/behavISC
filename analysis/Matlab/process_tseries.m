% STEP 2: 
clc
clear all;
cd('/Users/matthewbain/Documents/Science/Experiments/efflisisc/analysis/')
addpath(genpath('matlab'))
rng_cfig = rng(257);

% define file params
online     = 1;                   % read in files for online (1) or in-lab study (0)

% define params
scrambled = 1; % 1 if experiment involved scrambling; 0 otherwise

n_prtrials = 5;                   % # practice trials to reject from start of each ts
sf         = 2;                   % RT sampling rate
do_split   = 0;                   % perform median split (1 -> do split; 0 -> no split)
if do_split % (*) only modify these
    split_by   = 'enjoyment';     % dim to split by ('engagement' / 'eengagement'/'mental_sim'/'enjoyment'/'absorption'/'attention'           
    split_trun = 0;               % # of SSs to truncate split by (0 = all SSs incl in split)
else
    split_by   = '';             
    split_trun = 0;
end

do_save = 1;
overwrite = 1;

% -------------------------------------------------------------------------
% input
if online % processed log files
    logs = dir(fullfile('Individual_Processing/online/processed', '*.mat')); % processed logs directory
else
	logs = dir(fullfile('Individual_Processing', '*.mat'));
end
dxav = load(['Group_Processing' '/' 'av_processed' '_logtrans_0' '_online_' num2str(online) '.mat']); % group-level metrics ***make logtrans file param
dxav = dxav.dxav; 
    
% -------------------------------------------------------------------------
% get timing info ---------------------------------------------------------
timing_dir = dir(fullfile('Timings', '*.mat'));
for run = 1:length(timing_dir)
	trial_info(run) = load(['Timings' '/' timing_dir(run).name]);
end

% aggregate time series ---------------------------------------------------
% get # stories and versions
n_blocks = length(timing_dir);
n_ver    = length([trial_info(1).timings]);

% loop through intact/scrambled
for sc = 0:1
    if scrambled == 0 && sc == 1
        continue
    end 
    
    % savefile
    if scrambled == 0
        savefile = ['timeseries' '/' 'ts' '_dosplit_' num2str(do_split) '_splitby_' split_by '_splittrun_' num2str(split_trun) '_online_' num2str(online) '.mat'];
    elseif scrambled == 1
        savefile = ['timeseries' '/' 'ts' '_scrambled_' num2str(sc) '_dosplit_' num2str(do_split) '_splitby_' split_by '_splittrun_' num2str(split_trun) '_online_' num2str(online) '.mat'];
    end
        
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

    % set up data structure ***maybe don't intialize everything; consider high-D matrices
    for b = 1:n_blocks
        % all possible time points
        ts(b).x = (min([trial_info(b).timings{:}])) : sf : (max([trial_info(b).timings{:}]));

        % initialize other output fields
        ts(b).yagg = NaN(1, length(ts(b).x));
        
        if scrambled == 0
            ts(b).ypool = repmat({NaN(1, length(ts(b).x))}, [ceil(length(logs)/2/n_ver), n_ver]);
        elseif scrambled == 1
            ts(b).ypool = repmat({NaN(1, length(ts(b).x))}, [ceil(length(logs)/2/n_ver/2), n_ver]);
        end
        
        [ts(b).ystacked, ts(b).mean, ts(b).std, ts(b).len, ts(b).err] = deal([]);

        % initialize separate matrices to store pooled ts for splits (1/2 data)
        if do_split
            ypool(b).split1 = cell(floor(ceil(length(logs)/2/n_ver)/2), n_ver);
            ypool(b).split2 = cell(floor(ceil(length(logs)/2/n_ver)/2), n_ver);
        end
    end

    % add individual ts to data structure -------------------------------------
    for run = 1:length(logs)

        % --- get run info ---
        % load run data
        filename = logs(run).name;
        load([logs(run).folder '/' logs(run).name])

        % skip current participant if story not intact/scrambled
        if scrambled == 1
            if sc == 0
                if scram == story
                    continue
                end
            elseif sc == 1
                if scram ~= story
                    continue
                end
            end  
        end
        
        % get subj #
        if online
            subj_num = str2double(filename(3:(length(filename) - 23))); 
        else
            subj_num = str2double(filename(3:(length(filename) - 16))); 
        end

        % get current story (online story value is numerical)
        if ~online
            story = str2double(fname(6)); 
        end

        % define individual RT data
        y_ind = secs_norm;

        % --- add run ts info to aggregate ---
        % LOOP by Sc*****

        %if ~nnz(find(ismember(rej_files, run)))
        if ~dxav(story).rejected(subj_num) 

            % locate x indices for ind ts and add ts to new agg row
            if online
                timings = trial_info(story).timings{version};
            else
                timings = info.timings;
            end
            x_ix = find(ismember(ts(story).x, timings));
            ts(story).yagg(size(ts(story).yagg, 1), x_ix) = y_ind;

            % start new row of y agg unless on last subj
            if run ~= length(logs) && run ~= (length(logs) - 1)
                ts(story).yagg(size(ts(story).yagg, 1) + 1, :) =  NaN;
            end

            % get current version (online version value not in struct)
            if ~online
                version = info.version;
            end

            % get position of first unfilled row of version pool
            new_row = find(~any(~isnan(cell2mat(ts(story).ypool(:, version))), 2), 1);

            %***figure out if if statement necessary
            if isempty(new_row)
            else
                % start new row for appropriate version pool
                ts(story).ypool{length(ts(story).ypool(:, version)), version} = NaN(1, length(ts(story).x));   

                % store runs corresponding to elements in pooled ts ***don't need to get subj # here bc do above now (start of loop)
                if online
                    subj_num = filename(3:(length(filename) - 23)); 
                else
                    subj_num = filename(3:(length(filename) - 16)); 
                end
                ts(story).run(new_row, version) = str2double(subj_num);

                % add individual ts to version pool 
                ts(story).ypool{new_row, version}(x_ix) = y_ind;
            end
        end
    end

    % perform median split ----------------------------------------------------
    if do_split
        for s = 1:n_blocks

            % --- split time courses by version for correlation comparison ---
            for v = 1:n_ver

                % get fids for current version pool
                curr_fids = ts(s).run(:, v);
                curr_fids = curr_fids(curr_fids ~= 0); % remove zeros

                % if not even # of fids, remove final
                if mod(length(curr_fids), 2) ~= 0
                    curr_fids = curr_fids(1:end - 1);
                end

                % get ts indices of each split ***put split_ix into single array
                curr_scores = dxav(s).(split_by)(curr_fids); % scores on dim performing split on
                [~, sort_ix] = sort(curr_scores);
                split1_ix = sort_ix(1:length(sort_ix)/2 - split_trun);
                split2_ix = sort_ix(length(sort_ix)/2 + 1 + split_trun:end);

                % store split data in separate version pools
                ypool(s).split1(1:length(curr_scores)/2 - split_trun, v) = ts(s).ypool(split1_ix, v);
                ypool(s).split2(1:length(curr_scores)/2 - split_trun, v) = ts(s).ypool(split2_ix, v);
            end
        end 
    end

    % ---eliminate practice trials & stack ts --- ***clean up section (remove uneeded stats)
    for b = 1:n_blocks
        % reject practice trials from x and yagg
        ts(b).x(:, 1:n_prtrials)    = [];
        ts(b).yagg(:, 1:n_prtrials) = []; 

        for v = 1:n_ver
            % reject practice trials from pooled ts
            tspool = cell2mat(ts(b).ypool(:, v));
            tspool(:, 1:n_prtrials) = [];
            ts(b).ypool(:, v) = num2cell(tspool, 2);

            % do same for splits
            if do_split
                tspool1 = cell2mat(ypool(b).split1(:, v));
                tspool1(:, 1:n_prtrials) = [];
                ypool(b).split1(1:size(tspool1, 1), v) = num2cell(tspool1, 2);

                tspool2 = cell2mat(ypool(b).split2(:, v));
                tspool2(:, 1:n_prtrials) = [];
                ypool(b).split2(1:size(tspool2, 1), v) = num2cell(tspool2, 2);

                % replace all nan rows of ypool with empty arrays
                del_rows = find(~any(~isnan(cell2mat(ts(b).ypool(:, v))), 2));
                if ~isempty(del_rows)
                    ts(b).ypool(del_rows, v) = cell(1, 1);
                end
            end
        end

        % run stats on aggregate y values
        ts(b).std  = nanstd(ts(b).yagg);
        ts(b).len  = sum(~isnan(ts(b).yagg));
        ts(b).mean = nanmean(ts(b).yagg);    
        ts(b).err = nanstd(ts(b).yagg)./sqrt(sum(~isnan(ts(b).yagg))); %***simplify by replacing with std and len
        ts(b).err(find(ts(b).err == 0)) = NaN;

        % sequentially stack rows of pooled ts into aggregate ts w/ n_ver ts
        n_ss = min(sum(~cellfun(@isempty, ts(b).ypool))); % base # SSs on limiting vpool
        for s = 1:n_ss
            ts(b).ystacked(s, :) = nanmean(cell2mat(ts(b).ypool(s, :).'));
        end
        ts(b).ystacked(find(ts(b).err == 0)) = NaN; % ensure 0s removed for plotting err
    end

    % save --------------------------------------------------------------------
    if do_save
        if do_split
            save(savefile, 'ypool')
        else
            save(savefile, 'ts')
        end
    end
end
% ***archive this?
% STEP 1b:  
clc
clear all;
cd('/Users/matthewbain/Documents/Science/Experiments/efflisisc/analysis/')
addpath(genpath('matlab'))
rng_cfig = rng(257);

%(*)define params
cond_ids      = {'arctic' 'swim'};
logtrans      = 0;
do_save       = 0;
overwrite     = 0;
overwrite_ind = 0;

% input 
%[~, ~, raw] = xlsread('/Users/matthewbain/Documents/Science/Experiments/efflisisc/data/efflisisc_online.csv');
raw = readtable('/Users/matthewbain/Documents/Science/Experiments/efflisisc/data/efflisisc_online.xls');

% savefile
savefile_av  = ['Group_Processing' '/' 'av_processed' '_logtrans_' num2str(logtrans) '_online_1' '.mat']; % metrics averaged across subj
%savefile_rej = ['Group_Processing' '/' 'rejsubj_online_1' '.mat'];     % vector of blocks to reject

% check if savefile already exists
if do_save
    if exist(savefile_av, 'file')
        if overwrite 
            disp('overwriting existing savefile')
        else
            error(['file ' savefile_av ' exists already. Either delete or use new parameters!'])
        end
    end
end

% get timing info
timing_dir = dir(fullfile('Timings', '*.mat'));
n_blocks   = length(timing_dir);
for b = 1:n_blocks
	trial_info = load(['Timings' '/' timing_dir(b).name]);
    timings(b, :) = trial_info.timings;
end

% get header indices ------------------------------------------------------
% experiment section
%{
section_col  = find(strcmp(raw(1, :), 'designation'));

% subj ids
id_col       = find(strcmp(raw(1, :), 'identifier'));   % ordered by time of completion

% headphone check ('headphone-test')
hpass_col    = find(strcmp(raw(1, :), 'correct'));      % 1 -> correct; 0 -> incorrect - *hp/dt/comp*

% dual-task ('x-dual-task', x = 'swim'/'arctic')
story_col    = find(strcmp(raw(1, :), 'STORY'));        % 'arctic' or 'swim'
dtstim_col   = find(strcmp(raw(1, :), 'stimulus'));     % letter/comp q/engagement q id (EEx/MSx/Ex/Ax) - *dt/comp/eng*
dtresp_col   = find(strcmp(raw(1, :), 'key_press'));    % NA (miss) / 88 ('x') / 77 ('m')
rt_col       = find(strcmp(raw(1, :), 'rt'));           % in ms
dtcorr_col   = find(strcmp(raw(1, :), 'correct'));      
ttime_col    = find(strcmp(raw(1, :), 'ONSET'));        % theoretical onset of trials
tonset_col   = find(strcmp(raw(1, :), 'ONSET_ACTUAL')); % theoretical onset of trials

% comprehension qs ('x-comprehension')
compstim_col = find(strcmp(raw(1, :), 'stimulus'));
compresp_col = find(strcmp(raw(1, :), 'responses'));    % option selected
compcorr_col = find(strcmp(raw(1, :), 'correct'));

% engagement qs ('x-engagement')
engstim_col  = find(strcmp(raw(1, :), 'stimulus'));
engresp_col  = find(strcmp(raw(1, :), 'response'));     % 1-7

% demographic survey ('survey')
survq_col    = find(strcmp(raw(1, :), 'surveyQuestion'));
survresp_col = find(strcmp(raw(1, :), 'QRESP'));
%}

% define logfile fields based on headers
log_hpfield    = 'hpass';
log_dtfields   = {'dtstim', 'dtresp', 'rt', 'dtcorr', 'ttime', 'tonset', 'story', 'compresp', 'engresp', 'survq', 'survresp'}; 
raw_dtfields   = {'stimulus', 'key_press', 'rt', 'correct', 'ONSET', 'ONSET_ACTUAL', 'STORY', 'responses', 'response', 'surveyQuestion', 'QRESP'};

log_compfields = {'compstim', 'compresp', 'compcorr'}; 
raw_compfields = {'stimulus', 'responses', 'correct'};

log_engfields  = {'engstim', 'engresp'}; 
raw_engfields  = {'stimulus', 'response'}; 

log_survfields = {'survq', 'survresp'}; 
raw_survfields = {'surveyQuestion', 'QRESP'}; 

% mat to track # subjects rejected for failing each criterion
rej_fields = {'all', 'hp', 'comp', 'dt', 'misses'};
for f = 1:length(rej_fields)
    nrej.(rej_fields{f}) = 0;
end
rej_subj = []; % blocks (files) to reject

dxav = [];

% process individual data -------------------------------------------------
% subj list
%subjs = unique(raw(2:end, id_col));
subjs = unique(raw.identifier);

for s = 1:length(subjs)
    % struct to store start and end ix (col) of sections for each story (rows)
    sect_ix = [];
    
    % get 1st and last row for current subj
    %subj_range      = find(strcmp(raw(:, id_col), subjs(s)));
    subj_range      = find(strcmp(raw.identifier, subjs(s)));
    sect_ix(1).subj = [subj_range(1) subj_range(end)];
    
    % extract current subject raw data
    raw_curr = raw(sect_ix(1).subj(1):sect_ix(1).subj(2), :);
        
    % identify headphone check section
    %hp_range          = find(strcmp(raw_curr(:, section_col), 'headphone-test'));
    hp_range          = find(strcmp(raw_curr.designation, 'headphone-test'));
    
    % identify survey section
    surv_range        = find(strcmp(raw_curr.designation, 'survey'));
    sect_ix(1).surv   = [surv_range(1) surv_range(end)];
    
    for b = 1:n_blocks    
        % savefile/logfile for current block/subject
        logfile  = ['Individual_Processing' '/' 'online/' 'logs/' 'el' sprintf('%02d', s) 'b' num2str(b) '_online' '.mat'];
        savefile = ['Individual_Processing' '/' 'online/' 'processed/' 'el' sprintf('%02d', s) 'b' num2str(b) '_processed_online' '.mat'];
        
        % structs for store subj info and diagnostics
        log = [];
        dx  = [];
        
        % extract indices for each section --------------------------------
        % identify dual-task sections
        dt_range      = find(strcmp(raw_curr.designation, [cond_ids{b} '-dual-task']));
        sect_ix(b).dt = [dt_range(1) dt_range(end)];

        % identify comprehension sections
        comp_range      = find(strcmp(raw_curr.designation, [cond_ids{b} '-comprehension']));
        sect_ix(b).comp = [comp_range(1) comp_range(end)];

        % identify engagement sections
        eng_range       = find(strcmp(raw_curr.designation, [cond_ids{b} '-engagement']));
        sect_ix(b).eng  = [eng_range(1) eng_range(end)];

        % extract data to log ---------------------------------------------
        % headphone check
        %log.(log_hpfield) = raw_curr(hp_range, eval([log_hpfield '_col']));
        log.(log_hpfield) = raw_curr.correct(hp_range);
                    
        % dual-task
        for f = 1:length(log_dtfields)
            %log.(log_dtfields{f}) = raw_curr(sect_ix(b).dt(1):sect_ix(b).dt(2), eval([log_dtfields{f} '_col']));
            log.(log_dtfields{f}) = raw_curr.(raw_dtfields{f})(sect_ix(b).dt(1):sect_ix(b).dt(2));
        end
       
        % comprehension
        for f = 1:length(log_compfields)
            %log.(log_compfields{f}) = raw_curr(sect_ix(b).comp(1):sect_ix(b).comp(2), eval([log_compfields{f} '_col']));
            log.(log_compfields{f}) = raw_curr.(raw_compfields{f})(sect_ix(b).comp(1):sect_ix(b).comp(2));
        end
        
        % engagement
        for f = 1:length(log_engfields)
            %log.(log_engfields{f}) = raw_curr(sect_ix(b).eng(1):sect_ix(b).eng(2), eval([log_engfields{f} '_col']));
            log.(log_engfields{f}) = raw_curr.(raw_engfields{f})(sect_ix(b).eng(1):sect_ix(b).eng(2));
        end
        
        % survey
        for f = 1:length(log_survfields)
            %log.(log_survfields{f}) = raw_curr(sect_ix(1).surv(1):sect_ix(1).surv(2), eval([log_survfields{f} '_col']));
            log.(log_survfields{f}) = raw_curr.(raw_survfields{f})(sect_ix(1).surv(1):sect_ix(1).surv(2));
        end
            
        % compute dx ------------------------------------------------------
        % story id
        if strcmp(log.story{1}, 'arctic')
            dx.story = 1;
        elseif strcmp(log.story{1}, 'swim')
            dx.story = 2;
        end
        
        % version id (compare timing vectors to subj trial times)
        %ttimes = cell2mat(log.ttime');
        ttimes = log.ttime';
        for m = 1:length(timings)
            if length(ttimes) ~= length(timings{dx.story, m})
            elseif ttimes/1000 == timings{dx.story, m}
                dx.version = m;
            end
        end
        
        % --- demographics ---
        %dx.age = log.survresp{strcmp(log.survq, 'age')};
        dx.age = log.survresp(strcmp(log.survq, 'age'));
        dx.gender = {log.survresp{strcmp(log.survq, 'gender')}};
        dx.firstlang = {log.survresp{strcmp(log.survq, 'first-language')}};
        
        % --- quality check dx ---
        % headphone check score
        dx.hpcorr = nnz(cell2mat(log.hpass));
        
        % comprehension score
        dx.comp_score = nnz(cell2mat(log.compcorr)) * 10;
        
        % replace dtresp 'NA' w/ nan
        for i = 1:length(log.dtresp)
            if strcmp(log.dtresp{i}, 'NA')
                log.dtresp{i} = nan;
            end
        end
        
        % # missed trials
        dx.n_miss = nnz(isnan(cell2mat(log.dtresp)));

        % dual-task score (77 = 'm' = 'upper')
        dx.nrejtt_err = nnz( ~((cell2mat(log.dtresp) == 77) == cell2mat(isstrprop(log.dtstim, 'upper'))) );
        dx.dtscore = nnz( (cell2mat(log.dtresp) == 77) == (cell2mat(isstrprop(log.dtstim, 'upper'))) ) / length(log.dtstim);
        
        % --- reverse DT responses if subj mixed up letters ---
        dx.dt_rev = 0; % default value for reversal
        if dx.dtscore < .5 % if subj perf < chance, check score if resp reversed
            dtscore_rev = nnz( ~(cell2mat(log.dtresp) == 77) == (cell2mat(isstrprop(log.dtstim, 'upper'))) ) / length(log.dtstim);
            
            % reverse resp if reversed score > chance
            if dtscore_rev > .5 
                dx.dtscore = dtscore_rev;
                dx.dt_rev = 1;
            end
        end
        
        % average latency of trial onsets (ms)
        dx.av_latency = mean(cell2mat(log.tonset) - cell2mat(log.ttime));
        
        % proportion of dt targets that were uppercase
        dx.prop_upper = nnz(cell2mat(isstrprop((log.dtstim), 'upper'))) / length(log.dtstim);

        % --- RT dx --- 
        % replace RT 'NA' w/ nan and convert to matrix
        dx.secs = log.rt;
        for i = 1:length(dx.secs)
            if strcmp(dx.secs{i}, 'NA')
                dx.secs{i} = nan;
            end
        end
        dx.secs = cell2mat(dx.secs)/1000;
        
        % log transform RTs to make normally distributed
        if logtrans
            dx.secs = log10(dx.secs);
        end

        % letter case indices
        upper_ix = find(isstrprop([log.dtstim{:}], 'upper'));
        lower_ix = find(isstrprop([log.dtstim{:}], 'lower'));
        
        % RT descriptives by case        
        avRT_upp = nanmean(dx.secs(upper_ix));
        avRT_low = nanmean(dx.secs(lower_ix));
        sdRT_upp = nanstd(dx.secs(upper_ix));
        sdRT_low = nanstd(dx.secs(lower_ix));

        % reject trials (1.96 sd > av for corresponding case)     
        upp_rejix = find(abs(dx.secs(upper_ix) - avRT_upp) > 1.959*sdRT_upp);    
        low_rejix = find(abs(dx.secs(lower_ix) - avRT_low) > 1.959*sdRT_low); 
        
        % store # outliers
        dx.nrejtt_uppol = nnz(upp_rejix);
        dx.nrejtt_lowol = nnz(low_rejix);
        
        dx.secs(upp_rejix) = nan;
        dx.secs(low_rejix) = nan;
        
        dx.n_rej = nnz(upp_rejix) + nnz(low_rejix);
        
        % RT descriptives by case post-rejection
        dx.avRT_upp = nanmean(dx.secs(upper_ix));
        dx.avRT_low = nanmean(dx.secs(lower_ix));
        dx.sdRT_upp = nanstd(dx.secs(upper_ix));
        dx.sdRT_low = nanstd(dx.secs(lower_ix));
        
        dx.av_RT = nanmean(dx.secs);
        dx.sd_RT = nanstd(dx.secs);
        dx.min_RT = nanmin(dx.secs);
        dx.max_RT = nanmax(dx.secs);
        
        % normalize RTs based on upper/lowercase performance
        secs_norm = dx.secs;
        secs_norm(upper_ix) = (secs_norm(upper_ix) - dx.avRT_upp)/dx.sdRT_upp;
        secs_norm(lower_ix) = (secs_norm(lower_ix) - dx.avRT_low)/dx.sdRT_low;
        dx.secs_norm = secs_norm; % store RTs   
        
        % ---questionnaire dx---
        % extract non-numerical portions of narrative ids
        for m = 1:length(log.engstim) 
            eng_ids{m} = log.engstim{m}((find(isletter(log.engstim{m}))));
        end
    
        % store positions of different narrative qs for each dimension
        [unique_ids, ~, ixb] = unique(eng_ids, 'stable');

        % convert eng resp to mat
        engagement_mat = cell2mat(log.engresp);
        
        % calculate score for each dim
        for n = 1:length(unique_ids)
            curr_id = unique_ids{n}; % get new field name
            engagement_score.(curr_id) = mean(engagement_mat(find(ixb == n)));
        end
        dx.narrative_score = engagement_score;

        % calculate absorption
        dx.absorption = mean(engagement_mat(find([ixb == find(strcmp(unique_ids, 'EE')) | ...
            ixb == find(strcmp(unique_ids, 'A')) | ixb == find(strcmp(unique_ids, 'MS'))])));

        % store other metrics
        dx.engagement  = mean(engagement_mat); % aggregated engagement score
        dx.enjoyment   = engagement_score.E;
        dx.attention   = engagement_score.A;
        dx.eengagement = engagement_score.EE;
        dx.mental_sim  = engagement_score.MS;
            
        % reject subj -----------------------------------------------------
        if dx.comp_score < 70 || dx.dtscore < .85 || dx.n_miss > 14
            dx.rejected = 1;
            nrej.all = nrej.all + 1;
            if dx.hpcorr < 5 
                nrej.hp = nrej.hp + 1;
            end
            if dx.comp_score < 70
                nrej.comp = nrej.comp + 1;
            end
            if dx.dtscore < .85
                nrej.dt = nrej.dt + 1;
            end
            if dx.n_miss > 14
                nrej.misses = nrej.misses + 1;
            end
        else
            dx.rejected = 0;
        end
        
        %{
        %***SCRAP (just stored in dx)
        % get current file number and add to rejected subject vector ***just store this in dx
        if dx.rejected
            curr_fid = s*2 - 2 + b;
            rej_subj = [rej_subj curr_fid];
        end
        %}
        
        % store group-level dx --------------------------------------------
        % define all diagnostic metric names
        fields = fieldnames(dx);
        %fields(ismember(fields, {'secs' 'secs_norm'})) = []; % remove non-dxav fields
        
        % concatenate individual dx
        if size(dxav, 2) < b % if current block not yet in dxav
            for f = 1:length(fields)
                if ismember(fields{f}, {'secs', 'secs_norm'})
                    dxav(b).(fields{f}) = {dx.(fields{f})};
                else
                    dxav(b).(fields{f}) = [dx.(fields{f})];
                end
                
            end
        else
            for f = 1:length(fields)
                
                if ismember(fields{f}, {'secs', 'secs_norm'})
                    dxav(b).(fields{f}) = [dxav(b).(fields{f}), {dx.(fields{f})}];
                else
                    dxav(b).(fields{f}) = [dxav(b).(fields{f}), dx.(fields{f})];
                end
            end
        end
        
        % save ------------------------------------------------------------
        if do_save
            if exist(savefile) || exist(logfile)
                if overwrite_ind
                    fprintf('Overwriting %s.\n', savefile)
                    save(logfile, '-struct', 'log') 
                    save(savefile, '-struct', 'dx')
                else
                    disp('File not saved: block already exists.')
                end
            else
                save(logfile, '-struct', 'log') 
                save(savefile, '-struct', 'dx')
            end
        end
    end
    
end

% --- demo calculations ---
% final sample: array of 1/0 if subj rejected (i.e., both blocks rej)
for i = 1:length(dxav(1).rejected) 
    subj_incl(i) = ~isequal([dxav(1).rejected(i); dxav(2).rejected(i)], [1 1]');
end
    
n_female = nnz(strcmp(dxav(1).gender(subj_incl), 'Female'));
age_min  = min(dxav(1).age(subj_incl));
age_max  = max(dxav(1).age(subj_incl));
age_mean = nanmean(dxav(1).age(subj_incl));

%fprintf('proportion blocks rejected = %f\n', length(rej_subj) / (s*2))
fprintf('%d/%d complete subjects rejected.\n', nnz(any([dxav(1).rejected; dxav(2).rejected], 1)), s);
fprintf('%d/%d subjects had at least one run.\n', nnz(subj_incl), s);

% save group-level dx
if do_save
    save(savefile_av, 'dxav')
    %save(savefile_rej, 'rej_subj')
end
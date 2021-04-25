function trial_data = get_trial_data(expmt_list, expmt)
% Collect stimulation amplitude and pulse widths used for all trials
% pertaining to specific animal

data = expmt_list{expmt};

N_channels = length(data.elec_list);
N_trial = data.trial_list(end,end);
excluded_trials = find_skipped_trials(expmt_list, expmt);

trial_list = [];
stim_list = [];
selec_list = [];
pw_list = [];
for i=1:6
    N_ses_trials = length(nonzeros(data.stim_hist(i, :)));
    trial_list = [trial_list; (data.trial_list(i,1):data.trial_list(i,2))'];
    stim_list = [stim_list; nonzeros(data.stim_hist(i,:)')]; 
    selectGrid = sum(reshape(data.selectivityGrid(i,1:N_ses_trials,1:N_channels), N_ses_trials, N_channels),2);
    selec_list = [selec_list; selectGrid];
    pw_list = [pw_list; (linspace(data.pulseWidth(i),data.pulseWidth(i),N_ses_trials).*1000)'];
end

% Filling missing trial data
for trial=1:length(N_trial)
    if ismember(trial,excluded_trials)
        trial_list = [trial_list; trial];
        stim_list = [stim_list; nan];
        selec_list = [selec_list; nan];
        pw_list = [pw_list; nan];
    end  
end

    trial_data(:,1) = trial_list;
    trial_data(:,2) = stim_list;
    trial_data(:,3) = selec_list;
    trial_data(:,4) = pw_list;
    [~,idx] = sort(trial_data(:,1)); % sort just the first column
    trial_data = trial_data(idx,:);   % sort the whole matrix using the sort indices
    
end

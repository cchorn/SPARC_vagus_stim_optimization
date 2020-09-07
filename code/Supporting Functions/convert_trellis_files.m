%% Plot data from trellis file
% Made as script

% Written by: Jonathan Shulgach
% Updated: 5/8/20

% Changeable parameters
% ==============================================================
expmt       = 6; % Need to define animal and data file manually

% Set local path to trellis files
source_path = 'C:\Users\jonat\Documents\RNEL\source';

% Set new path to cloud repository
box_path    = 'C:\Users\jonat\Box\Data_for_Dylan';
% ==============================================================

STA_t_min   = -2;
STA_t_max   = 498;
fs          = 30e3;
N_channels  = 32;
plot_fig    = true;
chan_list   = 1:N_channels;
exp_data    = expmt_list{expmt};
stim_port   = exp_data.stim_port;
cohort      = exp_data.cohort;
cohort_folder = strrep(exp_data.cohort,'-','_');


N_trials = exp_data.trial_list(6,2); %last trial
skip_trials = find_skipped_trials(expmt_list, expmt);



for trial=1:N_trials
    
    if ~ismember(trial,skip_trials)==1
        
        ses_trial = find_trial(expmt_list, expmt, trial);
        session = find_session(expmt_list, expmt, trial);
        stim = expmt_list{expmt}.stim_hist(session,ses_trial);
        entity = exp_data.entity(session,:);  % necessary for reading trellis data files, specified as "raw", or "hi-res"
        
        % Channel 1 on headstage was used to stim nerve cuff electrode contacts 1-2, channel 9 on headstage for contacts 3-4
        if str2double(exp_data.cuff_list(session,1))==3
            stim_chan = 9;
        else
            stim_chan = 1;
        end
        if strcmpi(stim_port,'a')
            stimChan = stim_chan;
        elseif strcmpi(stim_port,'b')
            stimChan = 128 + stim_chan;
        elseif strcmpi(stim_port,'c')
            stimChan = 256 + stim_chan;
        elseif strcmpi(stim_port,'d')
            stimChan = 384 + stim_chan;
        end
        
        tic
        file_name = sprintf('datafile%04d.nev', trial);
        pw = num2str(exp_data.pulseWidth(session)*1000);
        % Path to trellis files should be on local machine for faster
        % loading
        %rec_filePath = fullfile([source_path,'\',cohort], file_name);
        rec_filePath = fullfile([source_path,'\F',cohort,'\trellis'], file_name);
        % Neuroshare function to read stimulation events file, rec_filePath contains directory to file starting
        %with "C:/...", stimChan has channel ID to read from
        %stimTimes = read_stimEvents(rec_filePath, stimChan);
        
        % save file path on box
        matfile = sprintf('F%s_trial_%d_data',cohort_folder,trial);
        matfilePath = fullfile([box_path,'\',cohort,'\.mat data files\',matfile]);
        
        
        % Neuroshare function to read all channels on continuous recording file. rec_filePath contains directory
        % to file starting with "C:/...", entity type for reading data based om sampling rate (ex: 'hi-res' = 2kHz,
        % 'raw' = 30Kz), elec_list is number array (1:32)       
        temp = read_continuousData(rec_filePath, char(entity), chan_list);
        % method 2, more advisable
        variable = struct();
        variable.(matfile) = temp;
        
        save(matfilePath, '-struct', 'variable')  % EDITED
        toc
    end
end

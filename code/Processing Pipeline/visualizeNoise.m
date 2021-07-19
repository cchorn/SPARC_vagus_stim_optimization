function [baseline] = visualizeNoise(expmt_list, expmt, source_path)
% Calculates standard deviation and mean potentials for each trial in one
% animal then plots the data

%   Written by: Dylan Beam on 3/3/2020
%   Last edited by: Jonathan Shulgach on 5/11/2020

% Inputs:
% =======================================================================
%   source_path: the path where all the raw data being processed is being stored
%   animal_num: the number of the animal to process (only one)
%   excluded_trials: an array of trial number(s) that should not be processed
%
% Outputs:
% =======================================================================
%   baseline: (struct) struct containing channel baseline noise info


% Changeable Parameters
% ===================================================================
cloud_path = true;
% ===================================================================

%determine constants
exp_data = expmt_list{expmt}; %extract meta data for the desired animal
stim_port = exp_data.stim_port; %extract the stim port value
%cohort = sprintf('F%s', exp_data.cohort); %the animal identifier
cohort = exp_data.cohort;
cohort_path = strrep(exp_data.cohort,'-','_');
rec_path = fullfile(source_path, cohort, 'Trellis'); %path to the folder with data
N_trials = expmt_list{expmt,1}.trial_list(6,2); %last trial number to be analyzed
break_point = expmt_list{expmt,1}.trial_list(3,2); %the trial number where stimulation switches from 1:2 to 3:4
fs=30e3; %sampling frequency
STA_bin = [-2 400]; %window for STA calculation where stim event is 0
STA_t_min = STA_bin(1)*fs/1000; %convert ms to s, then multiply by sample rate
STA_t_max = STA_bin(2)*fs/1000;
window_size = STA_t_max - STA_t_min; %the number of sampled points in the appropriate window size
buffer = 3;% 1.5; %time (ms) around the stimulation event to be removed
FailedTrials = zeros(32,N_trials); %indicates which trials were not successfully processed
excluded_trials = find_skipped_trials(expmt_list, expmt);
% FILTER, RANDOMIZE AND TAKE STATS FROM THE DATA

for chan_num = 1:32
    
    %loop through all trials
    for trial = 1:N_trials
        fprintf("Trial %d ", trial)
        try
            
            %exclude false trials
            if trial ~= excluded_trials
            else
                continue
            end
            
            %find the pulse width for this trial
            pw = find_pw(exp_data, trial);
            
            %name the file to use
            file_name = sprintf('datafile%04d.nev', trial);
            
            %determine which cuff pair is used, then define stim_chan as 1 or 9
            if trial > break_point
                stim_chan = 9;
            else
                stim_chan = 1;
            end
            
            %redefine stimChan based on the stimport
            switch stim_port
                case 'a'
                    stimChan = stim_chan;
                case 'b'
                    stimChan = 128 + stim_chan;
                case 'c'
                    stimChan = 256 + stim_chan;
                case 'd'
                    stimChan = 384 + stim_chan;
            end
            
            %extract the entity
            entity = expmt_list{expmt}.entity(1);
            
            mat_file_name = sprintf('F%s_trial_%d_data',cohort_path,trial);
            matfilePath = fullfile([source_path,'\',cohort,'\.mat data files\', mat_file_name,'.mat']);
            
            %create full path name for the recording file being read
            if cloud_path==true
                rec_filePath = fullfile([source_path,'\',cohort], file_name);
            else
                rec_filePath = fullfile([source_path,'\F',cohort,'\Trellis'], file_name);
            end
            %rec_filePath = fullfile(rec_path,rec_file);
            
            %Gather the timing of stim events, then convert into the index based on
            %the sampling frequency
            stimTimes = read_stimEvents(rec_filePath, stimChan);
            stimIdx = ceil(stimTimes{1}*fs);
            
            %adjust stim indices length if necessary
            if length(stimIdx) > 120 %this should only be happening when pulsewidth is 0.85 mS or more
                stimIdx = stimIdx(1:2:end);
            end
            
            
            %load in the data
            if cloud_path==true
                if exist('analogData','var')==0
                    temp = load(matfilePath);
                    analogData = temp.(mat_file_name);
                end
                analogData = analogData(chan_num,:);
            else
                analogData = read_continuousData(rec_filePath, char(entity), chan_num);
            end
            
            % As a precaution look for any missing data packets and replace
            % with NaN
            [newData, missing_packet_times, ~] = replace_missing_data_packets(analogData,8000);
            
            if ~isempty(missing_packet_times{1,1})
                a=1;
            end
            
            [blankedData, blankedOverlay] = blanking(newData, stimIdx, pw, fs, buffer);
            
            %apply filters
            [b,a] = butter(2,150/(30e3/2),'high');
            tmpChan = filtfilt(b,a,blankedData);
            % Implement 60Hz notch filter
            d = designfilt('bandstopiir','FilterOrder',2, ...
                'HalfPowerFrequency1',50,'HalfPowerFrequency2',70, ...
                'DesignMethod','butter','SampleRate',fs);
            tmpChan_notch = filtfilt(d,tmpChan);
            
            %create an array with stim events removedBand
            no_stim_data = no_stim(tmpChan_notch, stimIdx, pw, fs, buffer);
            
            %randomize index times within actual time range
            no_stim_Idx = int32((length(no_stim_data)-(window_size+1)) * rand(1,length(stimIdx)));
            
            %window-averages data, then calculates mean and STD
            temp = cell2mat(arrayfun(@(x) no_stim_data(x:(x+window_size)), no_stim_Idx,'UniformOutput',false)');
            
            % As a precaution remove any remaining artifacts in signal (no longer needed)
            %[snipData, ~] = remove_artifacts(temp,-300);
            
            avg_chanData = mean(snipData);
            baseline{expmt,1}.(['Chan_',num2str(chan_num)]).std_trials(trial) = std(avg_chanData);
            baseline{expmt,1}.(['Chan_',num2str(chan_num)]).mean_trials(trial) = mean(avg_chanData);
            
        catch ME
            
            sprintf('Error processing trial %d on channel %d: %s', trial, chan_num, ME.message)
            FailedTrials(chan_num, trial) = 1;
            
        end
        
    end
    
end

% PLOT MEAN AND STD VS TRIALS TO SEE IF THERE ARE ANY EXCESSIVELY NOISY TRIALS
create_baseline_figure(expmt_list, expmt, baseline)
%for chan_num = 1:32
%    figure(1)
%    subplot(4,8,chan_num)
%    plot(baseline{expmt,1}.(['Chan_',num2str(chan_num)]).std_trials)
    
%    figure(2)
%    subplot(4,8,chan_num)
%    plot(baseline{expmt,1}.(['Chan_',num2str(chan_num)]).mean_trials)
%end

end


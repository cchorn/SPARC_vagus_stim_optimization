%% This code is used to calculate the baseline noise for trial that failed in the initial trial

% Written by: Dylan Beam on 02/18/2020
% Last edited by: Dylan Beam on 02/28/2020

%% Define variables and load data

clc,clear

%this data will change by user, and animal being analyzed
rec_path = 'R:\users\dyb11\source\F19-19\Trellis'; %path to the folder with data
expmt = 4; %number corresponding to the animal selected (in above path)
num_trials = 45; %number of trials to analyze
break_point = 23; %the trial number where stimulation switches from 1:2 to 3:4
excluded_trials = [24]; %an array of all trials to be excluded

%should be constant across animals
fs=30e3; %sampling frequency
STA_bin = [-2 400]; %window for STA calculation where stim event is 0
STA_t_min = STA_bin(1)*fs/1000; %convert ms to s, then multiply by sample rate
STA_t_max = STA_bin(2)*fs/1000;
window_size = STA_t_max - STA_t_min; %the number of sampled points in the appropriate window size
buffer = 1.5; %time (ms) around the stimulation event to be removed

%load data
load('expmt_list_1_6_20.mat')
load('baseline.mat')
load('FailedTrials.mat')
exp_data = expmt_list{expmt};
stim_port = exp_data.stim_port;

%% Start processing the data

%identify channels that have failed trials
fail_chan = find(sum(FailedTrials'));

for i = 1:length(fail_chan)
    
    
        chan_num = fail_chan(i);

        trial_num = find(FailedTrials(chan_num,:));

        try 

        switch length(trial_num)
            case 1
                
                fprintf('Reprocessing trial %d on channel %d \n', trial_num(1), chan_num)
                
                %find the pulse width for this trial
                pw = find_pw(exp_data, trial_num(1));

                %name the file to use
                rec_file = sprintf('datafile%04d.nev', trial_num(1));

                %determine which cuff pair is used, then define stim_chan as 1 or 9 
                if trial_num(1) > break_point
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

                %create full path name for the recording file being read
                rec_filePath = fullfile(rec_path,rec_file);
                
                %Gather the timing of stim events, then convert into the index based on
                %the sampling frequency
                stimTimes = read_stimEvents(rec_filePath, stimChan);
                stimIdx = ceil(stimTimes{1}*fs);

                %adjust stim indices length if necessary
                if length(stimIdx) > 120 %this should only be happening when pulsewidth is 0.85 mS or more
                    stimIdx = stimIdx(1:2:end);
                end

                %load in the data
                analogData = read_continuousData(rec_filePath, char(entity), chan_num);
                
                %blank the data
                [blankedData, blankedOverlay] = blanking(analogData, stimIdx, pw, fs, buffer);

                %apply filters
                [b,a] = butter(2,150/(30e3/2),'high');
                tmpChan = filtfilt(b,a,blankedData);
                % Implement 60Hz notch filter
                d = designfilt('bandstopiir','FilterOrder',2, ...
                    'HalfPowerFrequency1',50,'HalfPowerFrequency2',70, ...
                    'DesignMethod','butter','SampleRate',fs);
                tmpChan_notch = filtfilt(d,tmpChan);

                %create an array of all data points before the first stim event
                pre_stim_data = tmpChan_notch(1:stimIdx(1)-(0.01*fs));

                %create an array with stim events removedBand 
                no_stim_data = no_stim(tmpChan_notch, stimIdx, pw, fs, buffer);

                %create an array with stim events removed and randomized
                no_stim_data_rand = no_stim_rand(tmpChan_notch, stimIdx, pw, fs, buffer);
                
                %randomize index times within actual time range
                no_stim_Idx = int32((length(no_stim_data)-(window_size+1)) * rand(1,length(stimIdx)));

                baseline{expmt,1}.(['Chan_',num2str(chan_num)]).mean_raw(trial_num(1)) = mean(tmpChan_notch);
                baseline{expmt,1}.(['Chan_',num2str(chan_num)]).std_raw(trial_num(1)) = std(tmpChan_notch); 

                %calculate mean and std for each trial at each channel (pres-stim data)
                baseline{expmt,1}.(['Chan_',num2str(chan_num)]).mean_pre(trial_num(1)) = mean(pre_stim_data);
                baseline{expmt,1}.(['Chan_',num2str(chan_num)]).std_pre(trial_num(1)) = std(pre_stim_data); 

                %calculate mean and std for each trial at each channel (pres-stim data)
                baseline{expmt,1}.(['Chan_',num2str(chan_num)]).mean_no(trial_num(1)) = mean(no_stim_data);
                baseline{expmt,1}.(['Chan_',num2str(chan_num)]).std_no(trial_num(1)) = std(no_stim_data);
                
                %calculate mean and std for each trial at each channel (no-stim
                %window averaged)
                snipData = cell2mat(arrayfun(@(x) no_stim_data(x:(x+window_size)), no_stim_Idx,'UniformOutput',false)');
                avg_chanData = mean(snipData);
                baseline{expmt,1}.(['Chan_',num2str(chan_num)]).std_win_avg(trial_num) = std(avg_chanData);
                baseline{expmt,1}.(['Chan_',num2str(chan_num)]).mean_win_avg(trial_num) = mean(avg_chanData);

                %calculate mean and std for each trial at each channel (no-stim
                %window averaged)
                snipDataRand = cell2mat(arrayfun(@(x) no_stim_data_rand(x:(x+window_size)), no_stim_Idx,'UniformOutput',false)');
                avg_chanDataRand = mean(snipDataRand);
                baseline{expmt,1}.(['Chan_',num2str(chan_num)]).std_win_avg_rand(trial_num) = std(avg_chanDataRand);
                baseline{expmt,1}.(['Chan_',num2str(chan_num)]).mean_win_avg_rand(trial_num) = mean(avg_chanDataRand);
            otherwise
                
                for j = 1:length(trial_num)
                    fprintf('Reprocessing trial %d on channel %d \n', trial_num(j), chan_num)
                    
                    %find the pulse width for this trial
                    pw = find_pw(exp_data, trial_num(j));

                    %name the file to use
                    rec_file = sprintf('datafile%04d.nev', trial_num(j));

                    %determine which cuff pair is used, then define stim_chan as 1 or 9 
                    if trial_num(j) > break_point
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

                    %create full path name for the recording file being read
                    rec_filePath = fullfile(rec_path,rec_file);
                    
                    %Gather the timing of stim events, then convert into the index based on
                    %the sampling frequency
                    stimTimes = read_stimEvents(rec_filePath, stimChan);
                    stimIdx = ceil(stimTimes{1}*fs);

                    %adjust stim indices length if necessary
                    if length(stimIdx) > 120 %this should only be happening when pulsewidth is 0.85 mS or more
                        stimIdx = stimIdx(1:2:end);
                    end

                    %load in the data
                    analogData = read_continuousData(rec_filePath, char(entity), chan_num);
                    
                    %blank the data
                    [blankedData, blankedOverlay] = blanking(analogData, stimIdx, pw, fs, buffer);

                    %apply filters
                    [b,a] = butter(2,150/(30e3/2),'high');
                    tmpChan = filtfilt(b,a,blankedData);
                    % Implement 60Hz notch filter
                    d = designfilt('bandstopiir','FilterOrder',2, ...
                        'HalfPowerFrequency1',50,'HalfPowerFrequency2',70, ...
                        'DesignMethod','butter','SampleRate',fs);
                    tmpChan_notch = filtfilt(d,tmpChan);

                    %create an array of all data points before the first stim event
                    pre_stim_data = tmpChan_notch(1:stimIdx(1)-(0.01*fs));

                    %create an array with stim events removedBand 
                    no_stim_data = no_stim(tmpChan_notch, stimIdx, pw, fs, buffer);

                    %create an array with stim events removed and randomized
                    no_stim_data_rand = no_stim_rand(tmpChan_notch, stimIdx, pw, fs, buffer);
                    
                    %randomize index times within actual time range
                    no_stim_Idx = int32((length(no_stim_data)-(window_size+1)) * rand(1,length(stimIdx)));

                    baseline{expmt,1}.(['Chan_',num2str(chan_num)]).mean_raw(trial_num(j)) = mean(tmpChan_notch);
                    baseline{expmt,1}.(['Chan_',num2str(chan_num)]).std_raw(trial_num(j)) = std(tmpChan_notch); 

                    %calculate mean and std for each trial at each channel (pres-stim data)
                    baseline{expmt,1}.(['Chan_',num2str(chan_num)]).mean_pre(trial_num(j)) = mean(pre_stim_data);
                    baseline{expmt,1}.(['Chan_',num2str(chan_num)]).std_pre(trial_num(j)) = std(pre_stim_data); 

                    %calculate mean and std for each trial at each channel (pres-stim data)
                    baseline{expmt,1}.(['Chan_',num2str(chan_num)]).mean_no(trial_num(j)) = mean(no_stim_data);
                    baseline{expmt,1}.(['Chan_',num2str(chan_num)]).std_no(trial_num(j)) = std(no_stim_data);
                    
                    %calculate mean and std for each trial at each channel (no-stim
                    %window averaged)
                    snipData = cell2mat(arrayfun(@(x) no_stim_data(x:(x+window_size)), no_stim_Idx,'UniformOutput',false)');
                    avg_chanData = mean(snipData);
                    baseline{expmt,1}.(['Chan_',num2str(chan_num)]).std_win_avg(trial_num(j)) = std(avg_chanData);
                    baseline{expmt,1}.(['Chan_',num2str(chan_num)]).mean_win_avg(trial_num(j)) = mean(avg_chanData);

                    %calculate mean and std for each trial at each channel (no-stim
                    %window averaged)
                    snipDataRand = cell2mat(arrayfun(@(x) no_stim_data_rand(x:(x+window_size)), no_stim_Idx,'UniformOutput',false)');
                    avg_chanDataRand = mean(snipDataRand);
                    baseline{expmt,1}.(['Chan_',num2str(chan_num)]).std_win_avg_rand(trial_num(j)) = std(avg_chanDataRand);
                    baseline{expmt,1}.(['Chan_',num2str(chan_num)]).mean_win_avg_rand(trial_num(j)) = mean(avg_chanDataRand);
                    
                end
        end
        
        %edit FailedTrials to know the data is processed properly
        FailedTrials(chan_num,trial_num) = 0;
        
    catch ME
        switch length(trial_num)
            case 1
                fprintf('Error reprocessing trial %d has failed on channel %d: %s \n', trial_num(1), chan_num, ME.message)
            case 2
                fprintf('Error reprocessing trial %d has failed on channel %d: %s \n', trial_num(j), chan_num, ME.message)
        end
    end
    
end
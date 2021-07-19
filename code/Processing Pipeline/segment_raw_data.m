function data_table = segment_raw_data(expmt_list, expmt, source_path)
%   Neural Rehabilitation Engineering Lab
%   Jonathan Shulgach
%   Last Updated: 8/5/20
%
%   This function reads trellis files and stores them into a struct for
%   offline analysis, making it easier and faster to analyze data if any
%   parameters need to be changed (e: threshold criteria)
%
% Inputs:
% ===================================================================
% source_path: (char) the path where all the raw data being processed is being stored
%
% Outputs:
% ==================================================================
% data_table: (struct) experiment information struct

% Changeable Parameters
% ===================================================================
cloud_path = true;
save_data = true;
% ===================================================================

%define high level variables
prompt = 'y';                          % Used for the auto-advancement of loading new files, and keeping track of stim amp per trial
data_table = struct();
fs            = 30e3;                  % sampling frequency
STA_bin       = [-2 498];              %time window used for STA calculations where the stim event is at time 0
STA_t_min     = STA_bin(1)*fs/1000;    %convert ms to s, then multiply by sample rate
STA_t_max     = STA_bin(2)*fs/1000;
blanking_window_offset = 3.0;

%Place files to ingore here...
ignore = {...
    %    'datafile0053.nev',...
    };

%loops through all animals you would like to process
for expmt=expmt
    if strcmpi(prompt,'n')
        break;
    end
    
    % Grab animal header information
    exp_data = expmt_list{expmt};
    stim_port = exp_data.stim_port;
    %cohort = ['F',strrep(exp_data.cohort,'-','_')];
    cohort = exp_data.cohort;
    cohort_path = strrep(exp_data.cohort,'-','_');
    chan_list = 1:32;
    if exp_data.rec_chan==2
        chan_list = chan_list + 32;   % Needed in case the recording cable was using a "Y-type" cable and connected to 2nd port
    end
    new_cmd = 1; % reset cmd_hist
    
    for session=1:size(exp_data.trial_list, 1) % for each trial list session (usually 6)
        fprintf("------------------------------------------------\n")
        fprintf("Session %d\n", session);
        fprintf("------------------------------------------------\n")
        start_trial = exp_data.trial_list(session, 1);
        end_trial = exp_data.trial_list(session, 2);
        entity = exp_data.entity(session,:);  % necessary for reading trellis data files, specified as "raw", or "hi-res"
        first_file = 1;
        %cohort = strrep(exp_data.cohort,'-','_');
        
        %Check for lower bound
        if size(exp_data.stimAmp,2)==1 %
            prev_cathAmp = 0;
        else
            prev_cathAmp = exp_data.stimAmp(session, 2);
        end
        cathAmp = exp_data.stimAmp(session, 1);
        
        % Channel 1 on headstage was used to stim nerve cuff electrode contacts 1-2, channel 9 on headstage for contacts 3-4
        if str2double(exp_data.cuff_list(session,1))==3
            stim_chan = 9;
        else
            stim_chan = 1;
        end
        
        % Trellis assigns 128 unique electrode channel ID's to each port. Depending on which Grapevine system port the stimulator
        % was plugged in, the stim channel assignment changes
        if strcmpi(stim_port,'a')
            stimChan = stim_chan;
        elseif strcmpi(stim_port,'b')
            stimChan = 128 + stim_chan;
        elseif strcmpi(stim_port,'c')
            stimChan = 256 + stim_chan;
        elseif strcmpi(stim_port,'d')
            stimChan = 384 + stim_chan;
        end
        
        for trial=1:length(start_trial:end_trial) % for all trials from recording
            actual_trial = start_trial + trial - 1;
            cathAmp = round(cathAmp/20)*20; % Force stim input to be multiple of 20uA (grapevine stim adapter hardware limitation)
            
            % Check to see that stim is not multiple of 20.
            if mod(cathAmp,20)~=0
                fprintf("Error with stimulus, not increment of 20uA. Recommend restarting at \n", round(cathAmp/20)*20);
                break;
            end
            %name the file to use
            file_name = sprintf('datafile%04d.nev', actual_trial);
            
            mat_file_name = sprintf('F%s_trial_%d_data',cohort_path,actual_trial);
            matfilePath = fullfile([source_path,'\',cohort,'\.mat data files\', mat_file_name,'.mat']);
            
            %if first_file==1        %Only ask for filepath once, just need folder location
                %rec_path = [source_path, exp_data.cohort]
                if cloud_path==true
                    rec_filePath = fullfile([source_path,'\',cohort], file_name);
                else
                    rec_path = [source_path, '\F', exp_data.cohort, '\Trellis'];
                end
                %rec_path  = ['R:\data_raw\ferret\2019\', exp_data.date_dash,' - aF', exp_data.cohort, '\Trellis'];
                %if isequal(rec_path,0)
                %    disp('User selected Cancel');
                %    return;
                %else
                    %new_file_num = start_trial;
                    %file_name = sprintf('datafile%04d.nev', trial);
                    %file_name = ['datafile',num2str(new_file_num,'%04d'),'.nev'];
                    %rec_filePath = fullfile(rec_path,file_name);
                    %first_file = 0;
                %end
            %else
                %new_file_num = str2double(file_name(end-7:end-4))+1;             % parse only the numbers in the file name, increment by one
                %file_name = sprintf('datafile%04d.nev', trial);
                %file_name = [file_name(1:8),num2str(new_file_num,'%04d'),'.nev']; % replace rec_file with new file name
                %rec_filePath = fullfile(rec_path, file_name);                    % rebuild path to file
            %end
            
            % Place files to skip here
            if sum(contains(ignore,file_name))==0 % Logical check whether recording file is in ignore list; if so, skip
                
                fprintf("%s", file_name);
                fprintf(" | Cohort: %s\n", exp_data.cohort);
                
                fprintf("------------------------------------------------------------\n")
                fprintf("Reading recording data\n")
                fprintf("------------------------------------------------------------\n")
                
                % Neuroshare function to read stimulation events file, rec_filePath contains directory to file starting
                %with "C:/...", stimChan has channel ID to read from
                stimTimes = read_stimEvents(rec_filePath, stimChan);
                
                % Load data
                if cloud_path==true
                    % Neuroshare function takes a long time to read in data, so faster method is to load the .mat files 
                    %if exist('analogData','var')==0
                    temp = load(matfilePath);
                    analogData = temp.(mat_file_name);
                    %end
                else
                    % Neuroshare function to read all channels on continuous recording file. rec_filePath contains directory
                    % to file starting with "C:/...", entity type for reading data based om sampling rate (ex: 'hi-res' = 2kHz,
                    % 'raw' = 30Kz), elec_list is number array (1:32)
                    analogData = read_continuousData(rec_filePath, char(entity), chan_list);
                end

                for chan = 1:length(chan_list)
                    fprintf('Analyzing Channel %s ', num2str(chan))
                    % Grab stim indices and muliply with samplig rate to convert to sample indices
                    stimIdx = ceil(stimTimes{1}*fs);
                    
                    if length(stimIdx) > 120 %this should only be happening when pulsewidth is 0.85 mS or more
                        stimIdx = stimIdx(1:2:end); % grab the actual start of the dual phase stim events
                    end
                    
                    % check that first (or more) stim happens after the user specified "before stim event" time
                    % has passed. Using a while loop will keep comparing the first stim index element until false
                    while stimIdx(1) < STA_t_min         % check that first stim is included within time window
                        stimIdx = stimIdx(2:end);
                    end
                    
                    [newData, missing_packet_times, ~] = replace_missing_data_packets(analogData(chan,:),8000);
                    
                    %find the pulsewidth for the trial for blanking
                    pw = find_pw(exp_data, actual_trial);
                    
                    %blank the data
                    [blankedData, blankedOverlay] = blanking(newData, stimIdx, pw, fs, blanking_window_offset);
                    
                    % Finla check to remove remaining artifacts from signal
                    
                    %                     keyboard
                    %----apply bandpass filter----
                    % Create 2nd order butterworth filter, high pass at 150Hz, with output parameters
                    [b,a] = butter(2,150/(30e3/2),'high'); %7/30/19 - used to be 300 but lowered to catch unknown CAPs just in case
                    tmpChan = filtfilt(b,a,blankedData);
                    
                    % Implement 2nd order butterworth 60Hz notch filter on top of filtered signal, band pass between 50-70Hz. Create
                    % filter object "d", then pass filter objetc into "filtfilt"
                    d = designfilt('bandstopiir','FilterOrder',2, ...
                        'HalfPowerFrequency1',50,'HalfPowerFrequency2',70, ...
                        'DesignMethod','butter','SampleRate',fs);
                    tmpChan_notch = filtfilt(d,tmpChan);
                    
                    % Using a MATLAB function "cell2mat", and the stim indices list "new_stimIdx" as the function input, 1) collect the signal
                    % data from "tmpChan_notch" indexed at x and everything within the time window (between STA_t_min and STA_t_max), 2) the
                    % store the resulting cell output into a matrix of values
                    temp = cell2mat(arrayfun(@(x) tmpChan_notch((x-abs(STA_t_min)):(x+abs(STA_t_max))), stimIdx,'UniformOutput',false)');
                    
                    % As a precaution remove any remaining artifacts in signal
                    %[snipData, ~] = remove_artifacts(temp,200);
                    
                    % create a list of numbers between two values, with a specified number of values in between
                    x_time = linspace(STA_t_min,STA_t_max,size(temp,2))*1000/fs;
                    
                    % Just using "pw" for string placeholder in figures and text
                    if strcmp(num2str(exp_data.pulseWidth(session)),'0.1')
                        pw = '100';
                    elseif strcmp(num2str(exp_data.pulseWidth(session)),'0.4')
                        pw = '400';
                    elseif strcmp(num2str(exp_data.pulseWidth(session)),'0.5')
                        pw = '500';
                    elseif strcmp(num2str(exp_data.pulseWidth(session)),'1')
                        pw = '1000';
                    end
                    
                    % % Similar placeholder function for cuff contact pair string in figures and text
                    if session<4
                        pair = '1_2';
                    else
                        pair = '3_4';
                    end
                    
                    % If true, fill struct with relevant data relating to PW used,
                    if save_data == true
                        data_table.(['F',cohort_path]).(['Cuff_',pair]).(['PW_',pw,'_us']){chan,trial} = temp;
                        data_table.(['F',cohort_path]).Sampling_Frequency = fs;
                        data_table.(['F',cohort_path]).Filter_Settings = '2nd order butter, high pass at 150Hz, notch at 60Hz with 10Hz+- half power band';
                        data_table.(['F',cohort_path]).x_time = x_time;
                    else
                        data_table = 0;
                    end
                    fprintf('\n')
                end
            end
            
            
            % Beginning of querying trial run commands from header data. During actual acute experiments, prompt is given by user to either continue
            % with yes command ('y'), no to exit ('no'), restart stim amp ('r'), skip file due to data corruption or noise ('s'), or reject algorithm
            % new stim ('m') and override with new greater ('h') or smaller ('l') stim amp, calculated as half the difference between current and previous stim
            %
            % ex: (last stim: 1500uA)
            % >> m;
            % >> "higher or lower (h or l)?: "; >> h;
            % >>"Increasing Stimulation - New Stimulation: 2260uA"
            % >> "stim difference: 1500uA"
            
            %manual_prompt = input(text, 's');
            prompt = exp_data.cmd_hist(new_cmd);
            new_cmd = new_cmd + 1;
            if strcmpi(prompt,'n')
                break;
            elseif strcmpi(prompt,'r') %restart loop with same stim params
                fprintf("restart\n")
                if save_data==true
                    data_table.(['F',cohort_path]).stim_hist(session, trial) = cathAmp;
                end
            elseif strcmpi(prompt,'y')
                fprintf("continue\n")
                if save_data==true
                    data_table.(['F',cohort_path]).stim_hist(session, trial) = cathAmp;
                end
                prev_cathAmp = cathAmp;
                cathAmp = round(new_cathAmp/20)*20;
            elseif strcmpi(prompt, 's') %skip due to corrupted file or other reason
                fprintf("skip\n")
                manual_prompt = exp_data.cmd_hist(new_cmd);
                new_cmd = new_cmd + 1;
                %skipped file decision
                if strcmpi(manual_prompt,'h')
                    new_cathAmp = round(round(cathAmp + abs(prev_cathAmp - cathAmp)/2)/20)*20;
                    fprintf("Increasing Stimulation - New Stimulation: %duA", new_cathAmp);
                elseif strcmpi(manual_prompt,'l')
                    final_cathAmp = cathAmp;
                    new_cathAmp = round(round(cathAmp - abs(prev_cathAmp - cathAmp)/2)/20)*20;
                    fprintf("Decreasing Stimulation - New Stimulation: %duA", new_cathAmp);
                else
                    fprintf("Unknown input, continuing with recommended stimulus")
                end
                prev_cathAmp = cathAmp;
                cathAmp = round(new_cathAmp/20)*20;
                fprintf(" | Stim difference: %duA\n", abs(prev_cathAmp - cathAmp))
                % Update new file name
                %new_file_num = str2double(file_name(end-7:end-4))+1;
                %file_name = [file_name(1:8),num2str(new_file_num,'%04d'),'.nev'];
                %rec_filepath = fullfile(rec_path, file_name);
                % Second Decision
                if strcmpi(manual_prompt,'h')
                    new_cathAmp = round(round(cathAmp + abs(prev_cathAmp - cathAmp)/2)/20)*20;
                    fprintf("Increasing Stimulation - New Stimulation: %duA", new_cathAmp);
                elseif strcmpi(manual_prompt,'l')
                    final_cathAmp = cathAmp;
                    new_cathAmp = round(round(cathAmp - abs(prev_cathAmp - cathAmp)/2)/20)*20;
                    fprintf("Decreasing Stimulation - New Stimulation: %duA", new_cathAmp);
                else
                    fprintf("Unknown input, continuing with recommended stimulus")
                end
                if save_data==true
                    data_table.(['F',cohort_path]).stim_hist(session, trial) = cathAmp;
                end
                prev_cathAmp = cathAmp;
                cathAmp = round(new_cathAmp/20)*20;
                fprintf(" | Stim difference: %duA\n", abs(prev_cathAmp - cathAmp))
            elseif strcmpi(prompt,'m')
                fprintf("manual override detected\n")
                %text = 'higher or lower (h or l)?: '; % commented out for automatic stepping through commands
                %manual_prompt = input(text, 's');
                manual_prompt = exp_data.cmd_hist(new_cmd);
                new_cmd = new_cmd + 1;
                if strcmpi(manual_prompt,'h')
                    new_cathAmp = round(round(cathAmp + abs(prev_cathAmp - cathAmp)/2)/20)*20;
                    fprintf("Increasing Stimulation - New Stimulation: %duA", new_cathAmp);
                elseif strcmpi(manual_prompt,'l')
                    final_cathAmp = cathAmp;
                    new_cathAmp = round(round(cathAmp - abs(prev_cathAmp - cathAmp)/2)/20)*20;
                    fprintf("Decreasing Stimulation - New Stimulation: %duA", new_cathAmp);
                else
                    fprintf("Unknown input, continuing with recommended stimulus")
                end
                if save_data==true
                    data_table.(['F',cohort_path]).stim_hist(session, trial) = cathAmp;
                end
                prev_cathAmp = cathAmp;
                cathAmp = round(new_cathAmp/20)*20;
                fprintf(" | Stim difference: %duA\n", abs(prev_cathAmp - cathAmp))
            end
            % Lots of this code can be put into separate functions to considerably clean up presentation.
            
        end
    end
end
end

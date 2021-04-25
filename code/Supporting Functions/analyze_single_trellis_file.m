%% Plot data from trellis file
% Made as script
% Written by: Jonathan Shulgach

% Changeable Parameters
% ===================================================================
expmt       = 3; % Need to define animal and data file
source_path = 'C:\Users\jonat\Documents\RNEL\source'; % Set local data path
box_path    = 'C:\Users\jonat\Box\Data_for_Dylan';
trial       = 51; % Trial number which matches data file number
plot_snips         = false; % Plot raw data, STA, and RMS of trial
plot_cv_shades     = false;
plot_STA_all_chans = true;
plot_STA_chan      = false;
plot_RMS_all_chans = true;
plot_RMS_chan      = false;
use_trellis_file = false; % choose either to load data from box or cell data table
use_cloud_path  = false; % Choose whether to read from local or cloud data repository
rms_offset = 15; %in step counts for RMS, actually 13 milliseconds
verbose = true;
offset_1_chans = 22;%[7,24,27];
%offset_1_chans = [1,2,3,4,5,6,8,9,10,11,12,13,14,16,17,20,21,22,25,27,29,32];% double manual offset
% =====================================================================

blanking_window_offset = 3.0; % Agreed upon on 4/28/20
STA_t_min   = -2;
STA_t_max   = 498;
fs          = 30e3;
N_channels  = 32;
file_name   = sprintf('datafile%04d.nev', trial); % data file name
exp_data    = expmt_list{expmt};

% Custom array for displaying channel outputs in appropriate subplots to
% match Utah MEA chnanel configuration on screen
elec_chan_map = [...
    1,  9, 17, 25,...
    2, 10, 18, 26,...
    3, 11, 19, 27,...
    4, 12, 20, 28,...
    5, 13, 21, 29,...
    6, 14, 22, 30,...
    7, 15, 23, 31,...
    8, 16, 24, 32];

cohort = exp_data.cohort;
cohort_path = strrep(exp_data.cohort,'-','_');
session = find_session(expmt_list, expmt, trial);
ses_trial = find_trial(expmt_list, expmt, trial);
stim = expmt_list{expmt}.stim_hist(session,ses_trial);
std_gain = expmt_list{expmt}.std_gain;
skip_channels = expmt_list{expmt,1}.exclude_channels;

N_trials = exp_data.trial_list(6,2); %last trial
skip_trials = find_skipped_trials(expmt_list, expmt);
sub_size = numSubplots(N_channels);

windowSize=30;
stepSize=3;
pw = find_pw(exp_data, trial);

rms_offset = rms_offset*10;

if use_cloud_path==true
    rec_filePath = fullfile([box_path,'\',cohort], file_name);
else
    rec_filePath = fullfile([source_path,'\F',cohort,'\Trellis'], file_name);
end

if use_trellis_file==true
    stimTimes = getStimEvents(expmt_list, expmt, rec_filePath);
    tic
    analogData = getAnalogData(expmt_list, expmt, rec_filePath, N_channels, use_cloud_path);
    toc
end

total_chan_resp = 0;
if verbose==true
    fprintf("\n----------------------------------------------------\n")
    fprintf("Trial:  %d | Stim : %duA | PW: %sus\n", trial, exp_data.stim_hist(session,ses_trial), pw);
    elec_update_bytes = fprintf('\nElectrode: %d of %d \n', 1, N_channels);
    spike_graphic = [repmat('| ',1,N_channels),'|'];
    spike_graphic_bytes = fprintf(spike_graphic);
end
for chan = 1:N_channels
    
    Vrms_list = [];
    %analogData2 = avg_signal_data{chan,trial};
    
    if use_trellis_file == true
        
        % Grab stim indices and muliply with samplig rate to convert to sample indices
        stimIdx = ceil(stimTimes{1}*fs);
        %stimIdx = find(tmpChan_notch > 100);
        
        if length(stimIdx) > 120 %this should only be happening when pulsewidth is 0.85 mS or more
            stimIdx = stimIdx(1:2:end); % grab the actual start of the dual phase stim events
        end
        
        % check that first (or more) stim happens after the user specified "before stim event" time
        % has passed. Using a while loop will keep comparing the first stim index element until false
        while stimIdx(1) < STA_t_min         % check that first stim is included within time window
            stimIdx = stimIdx(2:end);
        end
        
        % There is still the possibility that stimulation artifacts below 8mV
        % (railing) will be seen. Must keep this here 
        [newData, missing_packet_times, ~] = replace_missing_data_packets(analogData(chan,:),800);
        
        %blank the data
        [blankedData, blankedOverlay] = blanking(newData, stimIdx, pw, fs, blanking_window_offset);
        
        %----apply bandpass filter----
        % Create 2nd order butterworth filter, high pass at 150Hz, with output parameters
        [b,c_region] = butter(2,150/(30e3/2),'high'); %7/30/19 - used to be 300 but lowered to catch unknown CAPs just in case
        tmpChan = filtfilt(b,c_region,blankedData);
        
        % Implement 2nd order butterworth 60Hz notch filter on top of filtered signal, band pass between 50-70Hz. Create
        % filter object "d", then pass filter objetc into "filtfilt"
        d = designfilt('bandstopiir','FilterOrder',2, ...
            'HalfPowerFrequency1',50,'HalfPowerFrequency2',70, ...
            'DesignMethod','butter','SampleRate',fs);
        
        tmpChan_notch = filtfilt(d,tmpChan);
        
        sig_time = linspace(0,60,length(tmpChan_notch));
        
        % Using a MATLAB function "cell2mat", and the stim indices list "new_stimIdx" as the function input, 1) collect the signal
        % data from "tmpChan_notch" indexed at x and everything within the time window (between STA_t_min and STA_t_max), 2) the
        % store the resulting cell output into a matrix of values
        snipData = cell2mat(arrayfun(@(x) tmpChan_notch((x-abs(STA_t_min*fs/1000)):(x+abs(STA_t_max*fs/1000))), stimIdx,'UniformOutput',false)');
        % Collect segmented data
        STA_data = mean(snipData);
    else
        % Collect segmented data from loaded dataset
        STA_data = avg_signal_data{chan, trial};
        
    end
    
    % Create time period list for plotting data
    STA_time = linspace(STA_t_min,STA_t_max,length(STA_data));
    
    % Get the rms of the signal
    [Vrms_list, Vrms_time] = getRmsData(STA_data, 1);
    
    % Get response threshold for channel
    N_data = length(Vrms_list);
    %if ismember(chan, offset_1_chans)
    %    manual_offset = 0.7;
    %else
    manual_offset = -0.1;%0.5;
    %end
    % Some assistance with noisy channels
    [threshold, ~] = getResponseThreshold(expmt_list, expmt, baseline, trial, chan, N_data);%, manual_offset);
    
    %Detect response behavior for signal input and collect response times for conduction velocities
    [crossed, max_val, max_idx, time_bins, resp_bins] = response_detection(exp_data, threshold, Vrms_list, Vrms_time, rms_offset);
    
    
    if crossed==true
        if verbose==true
            fprintf(repmat('\b',1,(elec_update_bytes+spike_graphic_bytes+1)))
            elec_update_bytes = fprintf('\nElectrode: %d of %d | Spike!\n', chan, N_channels);
            spike_graphic(2*chan) = '!';
            spike_graphic_bytes = fprintf(spike_graphic)-1;
        end
        %[max_val, max_idx] = max(Vrms_list(rms_offset:end-rms_offset));
        peak_time = (max_idx+STA_t_min)*1000/fs;
        total_chan_resp = total_chan_resp + 1;
    else
        if verbose==true
            fprintf(repmat('\b',1,(elec_update_bytes+spike_graphic_bytes)))
            elec_update_bytes = fprintf('Electrode: %d of %d \n', chan, N_channels);
            spike_graphic(2*chan) = '.';
            spike_graphic_bytes = fprintf(spike_graphic);
        end
        max_idx=0;
    end
    
    %RMS
    if plot_RMS_all_chans==true
        titleName = sprintf('Trial %d PW%dus %duA RMS Window %d us, step %.3g us', trial, pw*1000, stim, windowSize/0.03, stepSize/0.03);
        plot_RMS_grid(Vrms_list, Vrms_time, titleName, chan, threshold, rms_offset, 0.2);
    end
    % STA
    if plot_STA_all_chans==true
        titleName = sprintf('Trial %d PW%dus %duA STA', trial, pw*1000, stim);
        plot_STA_grid(STA_data, STA_time, titleName, chan, threshold, rms_offset);
    end
    if plot_snips==true
        titleName = sprintf('Trial %d PW%dus %duA raw data', trial, pw*1000, stim);
        plot_raw_grid(snipData', STA_time, titleName, chan, threshold, rms_offset)
    end
    if plot_STA_chan==true
        figure(22)
        %plot(x_time, snipData, '-k','LineWidth',0.1);
        plot(STA_time,STA_data, '-k','LineWidth',0.1);
        title(['Trial ', num2str(trial),' Chan ', num2str(chan),' ', num2str(stim),'mA PW',num2str(pw*1000),'us raw'])
        ylim([-5,5]);
        ylabel('uV')
        xlabel('time (msec)')
        xlim([STA_t_min,STA_t_max]);
        drawnow
    end
    if plot_RMS_chan==true
        fig_title = ['Chan ', num2str(chan), ' Trial ', num2str(trial),' RMS Window ', num2str(windowSize/0.03), ' us, step ',num2str(stepSize/0.03),' us'];
        if plot_cv_shades==true
            plot_RMS_figure(Vrms_list, Vrms_time, fig_title, threshold, rms_offset, resp_bins, time_bins);
        else
            plot_RMS_figure(Vrms_list, Vrms_time, fig_title, threshold, rms_offset);
        end
    end
    if plot_RMS_chan==true || plot_STA_chan==true
        d=input('Next chan?: ');
    end
    
    
end

fprintf('Total channels responding: %d\n', total_chan_resp);

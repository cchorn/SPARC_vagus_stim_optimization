function acuteStimBinarySearch(cathElec, anodElec, stimAmp, preQuietT_s, stimT_s, postQuietT_s, freq, phaseDur_ms, elec_list, stim_port, rec_chan, entity, STA_bin, thresh_bin, enable_stim)
%   Neural Rehabilitation Engineering Lab
%   Jonathan Shulgach
%   Last Updated: 3/2/20
%
%   DESCRIPTION
%   ===================================================================
%   Automated optimization bisectional algorithm that:
%   1) performs stimulation trials using Trellis API
%   2) performs STA or spike threshold analysis on recording data
%   3) Adjust stimulation parameter until difference of 60uA
%
%   INPUTS
%   ===================================================================
%   cathElec        :   (int) cathode electrode 1 to 8
%   anodElec        :   (int) anode electrode 1 to 8
%   cathAmp         :   (1x2 numeric) stimulation values for stim trials
%                           where first value is current stim, second is
%                           previous stim value to define search window
%   preQuietT_s     :   (numeric) pre stim duration in seconds
%   stimT_s         :   (numeric) stim duration in seconds
%   postQuietT_s    :   (numeric) post stim duration in seconds
%   freq            :   (int) stimulatio nfrequency
%   phaseDur_ms     :   (int) phase duration of stim pulse
%                             (see assumption 4)
%   elec_list       :   (numeric) list of electrode channels to record
%                           [1:32]
%   stim_port       :   (char) grapevine port that stimulation headstage is
%                           plugged into
%   stim_chan       :   (int) channel
%   entity          :   (char vector) type of recording data to load
%   STA_bin         :   (1x2 numeric) time window for STA analysis
%   thresh_bin      :   (1x2 numeric) time window for threshold detection
%   enable_stim     :   (log) enabele or skip stim event
%
%   OUTPUTS
%   ===================================================================
%
%   EXAMPLE
%   ===================================================================
%   Ex 1: (using stim channels 1-2, 3mA 60sec stim at 2Hz, 1ms pw, electrodes
%           1-32, grapevine stim port B, rec port 1, raw analog data, STA and
%           spike threshold time window 2ms prior, 400ms post, stim enabled)
%   autoStimSeq_Shulgach(1, 2, [3000], 0, 60, 0, 2, 1, [1:32], 'b', 1, 1, 'raw', [2 400], [2 400], true)
%
%   Ex 2: (using stim channels 3-4, 760uA stim bound with 100uA, 60sec stim at 2Hz, 500us pw, electrodes
%           1-32, grapevine stim port B, rec port 1, raw analog data, STA and
%           spike threshold time window 2ms prior, 400ms post, stim disabled)
%   autoStimSeq_Shulgach(3, 4, [760 100], 0, 60, 0, 2, 0.5, [1:32], 'b', 1, 'raw', [2 400], [2 400], false)
%
%   NOTES
%   ===================================================================
%   1) Figure can be saved to data directory by uncommenting save section
%


% ====================== Initialize Prameters ============================
first_file = true;
trial_count = 1;
stim_diff = 30; % uA
fs=30e3; %sampling frequency
STA_t_min = STA_bin(1)*fs/1000; %convert ms to s, then multiply by sample rate
STA_t_max = STA_bin(2)*fs/1000;
search_offset = 3000;
plot_raw_data_fig = true;
plot_STA_data_fig = true;

N_channels = 32;

%=========================================================================


% Determine which cuff pair is used, then define stim_chan as 1 or 9
if cathElec==3
    stim_chan = 9;
else
    stim_chan = 1;
end

% Redefine stimChan based on the stimport (port "A" is chan 1, port "B" is
% chan 129, "C" would be 257, "D" would be 385
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

%Check for lower bound, if empty replace with 0
if length(stimAmp)==1
    prev_cathAmp = 0;
else
    prev_cathAmp = stimAmp(2);
end
cathAmp = stimAmp(1);

%Recording channel configuration (In case wire is dual port)
if rec_chan==2
    elec_list = elec_list + 32;
end

% Perform quiet trial for baseline data and collect
[base_file, base_path] = uigetfile('*.nev','Please select baseline data');
baseline = getBaselineParams(fullfile(base_file,base_path),1:N_channels);

% Display experiment information on command window for logging history and
% user feedback
fprintf("============================================================\n")
fprintf("Experiment Begin: %s, at %s\n", date, datestr(now,'HH:MM:SS'))
fprintf("Cuff Channels: %d-%d | Starting Stim: %duA | Frequency: %dHz | Phase duration: %0.3fms\n",cathElec, anodElec, stimAmp(1), freq, phaseDur_ms)
fprintf("Electrodes: %s\n", num2str(elec_list))
fprintf("Stim port: %s | Stim channel: %d | recording channel: %d | entity: %s\n", stim_port, stim_chan, rec_chan, entity)
fprintf("STA time bin: %dms to %dms | Threshold time bin: %dms to %dms\n",-STA_bin(1), STA_bin(2), -thresh_bin(1), thresh_bin(2))
LogicalStr = {'false', 'true'};
fprintf("Stimulation: %s\n", LogicalStr{enable_stim + 1})

% Provide prompt for user to continue with experiment when ready
text = 'Continue (y or n)?: ';
prompt = input(text, 's');
if strcmpi(prompt,'y')
    try
        % Begin binary search algorithm
        while abs(prev_cathAmp - cathAmp) > stim_diff
            
            % Force input to be multiple of 20
            cathAmp = round(cathAmp/20)*20;
            if mod(cathAmp,2)==1
                fprintf("Error with stimulus, not increment of 20uA. Recommend restarting at \n", round(cathAmp/20)*20)
                break;
            end
            
            % Get new values for omnetics stimulator as 4 channels signal
            % output ([1x4] step stim output values)
            cathAmp_step = getStimSteps(cathAmp);
            
            % Execute stimulator control in Trellis
            fprintf("------------------------------------------------------------\n")
            fprintf("Performing stimulation trial %d | Trial stim: %duA\n", trial_count, cathAmp)
            if enable_stim==true
                stimseqDesign(cathElec, anodElec, cathAmp_step, preQuietT_s, stimT_s, postQuietT_s, freq, phaseDur_ms);
            end
            fprintf("Stimulation trial done\n")
            fprintf("------------------------------------------------------------\n")
            
            % Ask user for file directory once, then auto-loop through files
            if first_file==true
                %Only ask for filepath once, most recent file
                [rec_file, rec_path] = uigetfile('*.nev','Please select stim recording data');
                rec_filepath = fullfile(rec_path,rec_file);
                first_file = false;
            else
                % grab next Trellis file saving should have auto-increment enabled
                new_file_num = str2double(rec_file(end-7:end-4))+1;
                rec_file = sprintf('datafile%04d.nev', new_file_num);
                rec_filepath = fullfile(rec_path, rec_file);
            end
            
            % Current file being read
            disp(rec_file)
            
            tic
            fprintf("------------------------------------------------------------\n")
            fprintf("Reading recording data\n")
            fprintf("------------------------------------------------------------\n")
            
            % Create grid with nice arrangement of figures for all channels
            sub_size = numSubplots(length(N_channels));
            %subplot(sub_size(1),sub_size(2), elec);
            
            % Graphic update in command window showing what electrode is being analyzed and plotted
            chanEvalGraphic('init', trial_counter, stimAmp, phaseDur_ms, N_channels);
            
            for chan = N_channels
                
                % Neuroshare function to read all channels on continuous recording file. rec_filePath contains directory
                % to file starting with "C:/...", entity type for reading data based om sampling rate (ex: 'hi-res' = 2kHz,
                % 'raw' = 30Kz), elec is channel to read)
                analogData = read_continuousData(rec_filePath, entity, chan);
                
                % Neuroshare function to read stimulation events file, rec_filePath contains directory to file starting
                %with "C:/...", stimChan has channel ID to read from. Muliply with samplig rate to convert to sample indices
                stimTimes = read_stimEvents(rec_filePath, stimChan);
                stimIdx = ceil(stimTimes{1}*fs);
                
                %adjust stim indices length if necessary
                if length(stimIdx) > 120 %this should only be happening when pulsewidth is 0.85 mS or more
                    stimIdx = stimIdx(1:2:end);
                end
                
                % Create 2nd order butterworth filter, high pass at 150Hz, with output parameters
                [b,a] = butter(2,150/(30e3/2),'high'); %7/30/19 - used to be 300 but lowered to catch unknown CAPs just in case
                tmpChan = filtfilt(b,a,analogData(chan,:));
                
                % Using a MATLAB function "cell2mat", and the stim indices list "stimIdx" as the function input, 1) collect the signal
                % data from "tmpChan_notch" indexed at x and everything within the time window (between STA_t_min and STA_t_max), 2) the
                % store the resulting cell output into a matrix of values
                snipData = cell2mat(arrayfun(@(x) tmpChan_notch((x-abs(STA_t_min)):(x+abs(STA_t_max))), stimIdx,'UniformOutput',false)');
                
                % Colletc average of signal snippets
                avg_chanData = mean(snipData);
                
                % Evaluate recorded aignal and identify CAPs in STA signals
                base_mean = baseline.(['Chan_',num2str(chan)]).mean;
                base_std = baseline.(['Chan_',num2str(chan)]).std;
                
                % Calculate the threshold which the signal must cross to
                % verify stimulation response
                thresh = base_mean + std_scale_low*base_std;
                
                % Narrow down threshold crossing evaluation range with the
                % STA signal to avoid detection of stimulation events
                search_window = avg_canData(STA_t_min+search_offset: end-search_offset);
                [s_freq, s_val] = hist(sig_thresh.s_search_window,100);
                
                % Plot raw signal figures
                if plot_raw_data_fig
                    plot_raw_data(sipData, STA_t_min, STA_t_max, chan, search_window)
                end           
                
                % Calculate the total number of threshold crossings within STA signal
                crossings = sum(abs(search_window) > thresh);
                if crossings >= cross_dur && max(search_window) < 3*max(search_window)
                    [max_val, max_idx] = max(search_window);
                    peak_time = (max_idx+STA_t_min+search_offset)*1000/fs;
                    chanEvalGraphic('spike', trial_counter, stimAmp, phaseDur_ms, N_channels);                  
                    STA_chan_List(elec) = 1;
                else
                    chanEvalGraphic('none', trial_counter, stimAmp, phaseDur_ms, N_channels);
                    STA_chan_List(elec) = 0;
                end
                
                % Plot STA signal figures
                if plot_STA_data_fig
                    plot_STA_data(sipData, STA_t_min, STA_t_max, chan, search_window, thresh, max_val, peak_time)
                end
            end
            
            % Make decision as to whether increase or decrease stim value
            fprintf("\n----------------------------------------------------\n")
            fprintf("Total Number of Channels with Responses: %d\n", sum(STA_chan_List))
            fprintf("----------------------------------------------------\n")
            if sum(STA_chan_List) >= 30
                fprintf("WARNING: Abnormally large amount of spikes detected....Recommend Restart\n")
            end
            
            if sum(STA_chan_List)>=1 % spike seen, lower stim
                new_cathAmp = round(round(cathAmp - abs(prev_cathAmp - cathAmp)/2)/20)*20;
                fprintf("Decreasing Stimulation - New Stimulation: %duA", round(new_cathAmp/20)*20);
            else
                new_cathAmp = round(round(cathAmp + abs(prev_cathAmp - cathAmp)/2)/20)*20;
                fprintf("Increasing Stimulation - New Stimulation: %duA", round(new_cathAmp/20)*20);
            end
            
            %trial_count = trial_count + 1;
            fprintf(" | Stim difference: %duA\n", abs(cathAmp - round(new_cathAmp/20)*20))
            
            beep
            % Manual Check (can disable)
            text = 'Continue (y or n), manual override (m), or restart (r)?: ';
            prompt = input(text, 's');
            switch prompt
                case 'n'
                    break;
                case 'r'
                    %restart loop with same stim params
                case 'y'
                    prev_cathAmp = cathAmp;
                    cathAmp = round(new_cathAmp/20)*20;
                case 'm'
                    text = 'higher or lower (h or l)?: ';
                    manual_prompt = input(text, 's');
                    if strcmpi(manual_prompt,'h')
                        new_cathAmp = round(round(cathAmp + abs(prev_cathAmp - cathAmp)/2)/20)*20;
                        fprintf("Increasing Stimulation - New Stimulation: %duA", new_cathAmp);
                    elseif strcmpi(manual_prompt,'l')
                        new_cathAmp = round(round(cathAmp - abs(prev_cathAmp - cathAmp)/2)/20)*20;
                        fprintf("Decreasing Stimulation - New Stimulation: %duA", new_cathAmp);
                    else
                        fprintf("Unknown input, continuing with recommended stimulus")
                    end
                    prev_cathAmp = cathAmp;
                    cathAmp = round(new_cathAmp/20)*20;
                    fprintf(" | Stim difference: %duA\n", abs(prev_cathAmp - cathAmp))
            end
        end
    catch
        sprintf('Error running algorithm. Stopped with file %d at stim %d: %s', rec_file, cathAmp, ME.message)
        return
    end
end

fprintf("============================================================\n")
fprintf("Final stimulatim: %duA\n", round(prev_cathAmp/20)*20)
fprintf("done\n")

end

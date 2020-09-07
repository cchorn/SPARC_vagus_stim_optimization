function expmt_list = analyze_stim_data_final(varargin)
%   DESCRIPTION
%   ===================================================================
%   Calculates the false positive and true positive rate based on given
%   response data
%
%   Jonathan Shulgach
%   Last updated: 8/26/2020
%
%   INPUTS
%   ===================================================================
%   expmt_list      :   (struct) data file header for experiments
%   baseline        :   (string) file name for baseline data to be collected
%   avg_signal_data :   (MxN cell) Cell array with STA average data signals
%                           for M channels and N trials
%   expmt           :   (numeric) experiment number
%
%   OUTPUTS
%   ===================================================================
%   ratio_data         :   (1x1 struct) struct containing TP,FP,TN,FN rates and gain values tested for a
%                               range of RMS window and step sizes
%
%   EXAMPLE
%   ===================================================================
%   expmt_list = analyze_stim_data_final(expmt_list, baseline, avg_signal_data, expmt)

expmt_list = varargin{1};
baseline = varargin{2};
avg_signal_data = varargin{3};
expmt = varargin{4};

% Adjust these parameters depending on the animal
plot_RMS_chan   = false;
plot_RMS_all_chans = false;
plot_hist = false;

verbose = true; % enable to see trial channel responses to stim expressed on command line

STA_t_min   = -2;
STA_t_max   = 498;
fs = 30e3;
rms_offset = 13; %in step counts for RMS, actually 13 milliseconds taking into consideration the -2 startpoint

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

exp_data = expmt_list{expmt};
expmt_list{expmt}.selectivityGrid = []; % overwrite previous values
[N_channels, N_trials] = size(avg_signal_data);
sub_size = numSubplots(N_channels);

%conductVel = nan(1,length(N_channels));
%sub_size = numSubplots(N_channels);
selectivityList = zeros(1,N_channels);
expmt_list{expmt}.minthresh = nan(N_channels,6);
expmt_list{expmt}.maxstim = nan(N_channels,6);
expmt_list{expmt}.selectivityGrid = zeros(size(exp_data.stim_hist,1),size(exp_data.stim_hist,2),N_channels);
expmt_list{expmt}.cv = [];

% Update time from last run (helps for determining  what datasets need
% to be re-analyzed/updated)
expmt_list{expmt,1}.last_run = char(datetime('now'));


skip_trials = expmt_list{expmt,1}.exclude_trials;
skip_channels = expmt_list{expmt,1}.exclude_channels;
% Get period of data actually being evaluated (with respect
% to the offset) so we don't analyze portions of the signal
% affected by the stim events
rms_offset = (abs(STA_t_min) + rms_offset)*10;

if verbose==true
    cohort = exp_data.cohort;
    fprintf(['\nAnimal ',cohort,'\n'])
    std_gain = expmt_list{expmt}.std_gain;
    windowSize=30;
    stepSize=3;
    fprintf(['\nGain: ', num2str(std_gain),' | RMS Window: ', num2str(windowSize/0.03),', Step: ', num2str(stepSize/0.03)])
end



for trial=1:N_trials
    
    if ismember(trial,skip_trials)==false

        session = find_session(expmt_list, expmt, trial);
        ses_trial = find_trial(expmt_list, expmt, trial);
        pw = find_pw(exp_data, trial);
        
        % Graphic update in command window showing what electrode is being analyzed and plotted
        if verbose==true
            fprintf("\n----------------------------------------------------\n")
            fprintf("Trial:  %d | Stim : %duA | PW: %sus\n", trial, exp_data.stim_hist(session,ses_trial), pw);
            elec_update_bytes = fprintf('\nElectrode: %d of %d \n', 1, N_channels);
            spike_graphic = [repmat('| ',1,N_channels),'|'];
            spike_graphic_bytes = fprintf(spike_graphic);
        end
        
        for chan=1:N_channels
            if ismember(chan,skip_channels)==true
                selectivityList(chan) = 0;
            else
                
                % Collect segmented data
                STA_data = avg_signal_data{chan, trial};
                
                % Get the rms of the signal
                [Vrms_list, Vrms_time] = getRmsData(STA_data, 1);
                
                % Get response threshold for channel
                N_data = length(Vrms_list);
                [threshold, ~] = getResponseThreshold(expmt_list, expmt, baseline, trial, chan, N_data);
                
                %Detect response behavior for signal input and collect response times for conduction velocities
                [crossed, ~, ~, time_bins, resp_bins] = response_detection(exp_data, threshold, Vrms_list, Vrms_time, rms_offset);
                
                % Save time bin responses because they record CV in steps of 0.5m/s
                expmt_list{expmt}.cv(chan,trial).time_bins = time_bins;
                expmt_list{expmt}.cv(chan,trial).resp_bins = resp_bins;
                
                if crossed==true % ----- Response detected -----
                    if verbose==true
                        fprintf(repmat('\b',1,(elec_update_bytes+spike_graphic_bytes+1)))
                        elec_update_bytes = fprintf('\nElectrode: %d of %d | Spike!\n', chan, N_channels);
                        spike_graphic(2*chan) = '!';
                        spike_graphic_bytes = fprintf(spike_graphic)-1;
                    end
                    
                    %Update list that holds which channels have response
                    selectivityList(chan) = 1;
                    detected_response_data(chan,trial) = 1;
                    
                    % Collect cuff threshold data
                    if expmt_list{expmt}.stim_hist(session,ses_trial) < expmt_list{expmt}.minthresh(chan,session) || isnan(expmt_list{expmt}.minthresh(chan,session))
                        expmt_list{expmt}.minthresh(chan,session) = expmt_list{expmt}.stim_hist(session,ses_trial);
                    end
                    
                    %Collect max stim channel response data
                    if (expmt_list{expmt}.stim_hist(session,ses_trial) > expmt_list{expmt}.maxstim(chan,session)) || isnan(expmt_list{expmt}.maxstim(chan,session))
                        expmt_list{expmt}.maxstim(chan,session) = expmt_list{expmt}.stim_hist(session,ses_trial);
                    end
                    %peakTimes(chan) = (max_idx(1)+STA_t_min)*1000/fs;
                    
                else % ----- No response detected ----
                    
                    if verbose==true
                        fprintf(repmat('\b',1,(elec_update_bytes+spike_graphic_bytes)))
                        elec_update_bytes = fprintf('Electrode: %d of %d \n', chan, N_channels);
                        spike_graphic(2*chan) = '.';
                        spike_graphic_bytes = fprintf(spike_graphic);
                    end
                    selectivityList(chan) = 0;
                    %peakTimes(chan) = 0;
                    detected_response_data(chan,trial) = 0;
                end
                prompt_bits=0;
                if plot_RMS_chan==true
                    
                    fig_title = ['Chan ', num2str(chan), ' Trial ', num2str(trial),' RMS Window ', num2str(windowSize/0.03), ' us, step ',num2str(stepSize/0.03),' us'];
                    plot_RMS_figure(Vrms_list, Vrms_time, fig_title, threshold, rms_offset, resp_bins, time_bins);
                    %plot_RMS_and_data_figure(); % for showing both STA and RMS side-by-side
                    
                    if verbose==true
                        % 7th element is 3m/s, which is the max CV for vagus C-Fibers (Lee Fisher - 6/23/20)
                        fprintf("\nTotal bins with response: %d\n", sum(resp_bins));
                        t_start_idx = find(floor(10*Vrms_time)/10==floor(time_bins(7)),1); % get index of time bin mark WRT Vrms time
                        t_end_idx = find(floor(10*Vrms_time)/10==floor(Vrms_time(end-rms_offset)),1);
                        if sum(Vrms_list(t_start_idx:t_end_idx)>threshold)>0
                            fprintf("C-fiber: YES");
                        else
                            fprintf("C-fiber: NO");
                        end
                        if floor(time_bins(end)) < rms_offset/10
                            t_start = floor(rms_offset/10);
                        end
                        t_start_idx = find(floor(10*Vrms_time)/10==t_start,1); % get index of time bin mark WRT Vrms time
                        t_end_idx = find(floor(10*Vrms_time)/10==floor(time_bins(7)),1);
                        if sum(Vrms_list(t_start_idx:t_end_idx)>threshold)>0
                            fprintf("\t aDelta-fiber: YES\n");
                        else
                            fprintf("\t aDelta-fiber:  NO\n");
                        end
                        prompt_bits = prompt_bits + 16;
                    end
                    
                end
                if plot_RMS_all_chans == true
                    titleName = sprintf('Trial %d PW%sus RMS Window %d us, step %.3g us', trial, num2str(pw*1000), windowSize/0.03, stepSize/0.03);
                    plot_RMS_grid(Vrms_list, Vrms_time, titleName, chan, threshold, rms_offset);
                end
                
                
                %Histogram Plot
                if plot_hist==true
                    figure(7)
                    hist(Vrms_list(rms_offset:end),100);
                    h = findobj(gca,'Type','patch');
                    h.FaceColor = [0 0.4470 0.7410];
                    %hist(rms_threshold.rms_search_window,100)
                    hold on;
                    plot(thresh, 1:250, 'r.'); %RMS line plot
                    plot(thresh, 1:250, 'r.'); %RMS line plot
                    xlabel('Signal');
                    ylabel('Frequency');
                    title('Signal Histogram')
                    hold off;
                    drawnow
                end
                
                if plot_hist==true || plot_RMS_chan==true
                    a=input('Next channel: ');
                end
                
                if verbose==true && plot_RMS_chan ==true
                    fprintf(repmat('\b',1,61+prompt_bits))
                elseif verbose==true
                    %fprintf(repmat('\b',1,1))
                end
                %if save_fig==true
                %    sgtitle(sprintf('Cohort %s | Trial %d | Channel %d | Spike: %s',strrep(cohort,'_','-'),ct,chan,spike))
                %    save_path = [pwd, '\F19_19_ground_truth_eval\', 'Trial',num2str(trial),'_Chan',num2str(chan)];
                %    saveas(gcf, [save_path '.png']);
                %end
                %keyboard
            end
        end
        if plot_RMS_all_chans==true
            b=input('\nNext trial: ');
        end
        expmt_list{expmt}.selectivityGrid(session, ses_trial, :) = selectivityList;
    end
end


fprintf("\nDone!\n")
end
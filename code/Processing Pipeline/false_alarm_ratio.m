function [ratio_data, detected_response_data, incorrect_response_data] = false_alarm_ratio(varargin)
%   DESCRIPTION
%   ===================================================================
%   Calculates the false positive and true positive rate based on given
%   response data
%
%   Jonathan Shulgach
%   Last updated: 5/4/2020
%
%   INPUTS
%   ===================================================================
%   expmt_list      :   (struct) data file header for experiments
%   baseline        :   (string) file name for baseline data to be collected
%   avg_signal_data :   (MxN cell) Cell array with STA average data signals
%                           for M channels and N trials
%   ground_truth_data:  (MxN double) numeric array containing boolean values representing responses on channels
%   expmt           :   (numeric) experiment number
%
%   OUTPUTS
%   ===================================================================
%   ratio_data         :   (1x1 struct) struct containing TP,FP,TN,FN rates and gain values tested for a
%                               range of RMS window and step sizes
%
%   EXAMPLE
%   ===================================================================
%
%   >> ratio_data = false_alarm_ratio(expmt_list, baseline, avg_signal_data, ground_truth_dataset, expmt, excluded_trials)

expmt_list = varargin{1};
baseline = varargin{2};
avg_signal_data = varargin{3};
ground_truth_dataset = varargin{4};
expmt = varargin{5};

% Adjust these parameters depending on the animal
N_channels  = 32;
plot_fig1   = false;
save_fig    = false;
fs = 30e3;
gain_list=linspace(2,4,41);
rms_offset=150; %in step counts for RMS, actually 15 milliseconds
STA_t_min     = -2;
exp_data = expmt_list{expmt};
N_trials=length(ground_truth_dataset);
cohort = ['F',strrep(exp_data.cohort,'-','_')];
if exist('ct','var') == 0
    ct = 1;
end

% Find trial/datafile numbers not included in actual trial list
skip_trials = find_skipped_trials(expmt_list, expmt);

% Add additional trials to skip in case bad trials are detected from
% shivering or increased noise
if length(varargin)>5
    skip_trials = sort([skip_trials, varargin{6}]);
end


% Iterate over 2 window sizes (500us and 1000us) for RMS window
for windowSize=30%15:15:30
    fprintf(['Window Size:',num2str(windowSize)])
    
    % Iterate from 100us to 1000us step sizes for RMS window
    fprintf('\nStep Size:')
    for stepSize=3%:3:12
        fprintf([num2str(stepSize),', ']);
        % Refresh constants for every stepsize or window
        TP_rate = 0;
        TN_rate = 0;
        FP_rate = 0;
        FN_rate = 0;
        
        % Iterate through gain list
        for gain_ct=1:length(gain_list)
            TP = 0;
            FP = 0;
            TN = 0;
            FN = 0;
            fprintf(['\nGain: ', num2str(gain_list(gain_ct)),' | trial '])
            for trial=1:N_trials
                fprintf([num2str(trial),', ']);
                Vrms_list = [];
                for chan=1:N_channels
                    
                    % Check whether trial is included in the list of trials
                    % to skip, and if so... skip and replace all detected
                    % trial responses to zero
                    if ismember(trial,skip_trials)==1
                        detected_response_data(chan,trial) = 0;
                    else
                        
                        % Collect segmented data
                        STA = avg_signal_data{chan, trial};

                        % STA_t_max now defined by length of data instead
                        % of explicitely
                        STA_t_max = round(length(STA(abs(STA_t_min)*fs/1000:end))*1000/fs);
                        
                        % Create time period list for plotting data
                        x_time = linspace(STA_t_min,STA_t_max,length(STA));
                        
                        %calculate Vrms with one window/step size
                        window = 1:windowSize;
                        idx = 1;
                        while max(window) < length(STA)
                            Vrms_list(idx) = rms(STA(window));
                            window = window + stepSize;
                            %window = window + (windowSize-stepSize);
                            idx = idx + 1;
                        end
                        window = window(1):length(STA);
                        Vrms_list(idx) = rms(STA(window));
                        Vrms_time = linspace(STA_t_min,STA_t_max,length(Vrms_list));
                        
                        % Find threshold by obtaining mean and std
                        % from "baseline" struct, multiply std with gain
                        % and add with mean to establish threshold value
                        base_std = baseline{expmt,1}.(['Chan_',num2str(chan)]).std_trials(trial);
                        base_mean = baseline{expmt,1}.(['Chan_',num2str(chan)]).mean_trials(trial);
                        std_gain = gain_list(gain_ct);
                        thresh = base_mean + std_gain*base_std;
                        threshold = linspace(thresh,thresh,length(Vrms_list));
                        x_thresh = linspace(STA_t_min, STA_t_max,length(Vrms_list));
                        
                        
                        % For STA check all instances of threshold crossings and evaluate whether a
                        % crossing of 1 millisecond occurs denoting true
                        % response. For RMS just see if it crosses
                        crossed = false;
                        cross_dur = 1;%round(1000/(stepSize/0.03)); % convert milliseconds to elements
                        for i=rms_offset:(length(Vrms_list)-rms_offset)
                            if sum(Vrms_list(i:i+cross_dur)>thresh)>1% just cross %cross_dur
                                crossed = true;
                            end
                        end
                        if crossed==true
                        %if sum(Vrms_list(rms_thresh_offset:end-rms_thresh_offset) > threshold(rms_thresh_offset:end-rms_thresh_offset)) >= 1
                            detected_response_data(chan,trial) = 1;
                            [max_val, max_idx] = max(Vrms_list(rms_offset:end-rms_offset));
                            peak_time = (max_idx+STA_t_min)*1000/fs;
                        else
                            max_idx=0;
                            detected_response_data(chan,trial) = 0;
                        end
                        
                        if plot_fig1==true
                            figure(6)
                            sgtitle(sprintf('Cohort %s | Trial %d | Channel %d |          ',strrep(cohort,'_','-'),trial,chan))
                            subplot(1,2,1)
                            plot(x_time, STA, '-k','LineWidth',0.1);
                            title('Filtered signal')
                            ylim([-10,10]);
                            xlim([STA_t_min,STA_t_max]);
                            drawnow
                            titleName = sprintf('RMS Window %d us, step %.3g us', windowSize/0.03, stepSize/0.03);
                            subplot(1,2,2)
                            plot(Vrms_time, Vrms_list)
                            hold on
                            plot(Vrms_time,threshold,'-r','LineWidth', 1);
                            if max_idx~=0
                                plot((rms_offset+max_idx)/10,threshold(1),'or')
                            end
                            hold off
                            title(titleName)
                            axis tight
                            ylim([-1,5])
                            drawnow
                            if save_fig==true
                                sgtitle(sprintf('Cohort %s | Trial %d | Channel %d | Spike: %s',strrep(cohort,'_','-'),ct,chan,spike))
                                save_path = [pwd, '\F19_19_ground_truth_eval\', 'Trial',num2str(trial),'_Chan',num2str(chan)];
                                saveas(gcf, [save_path '.png']);
                            end
                            %keyboard
                        end
                        
                        % Evaluate True Positive and False Positive rates
                        if detected_response_data(chan,trial)==ground_truth_dataset(chan,trial)
                            if ground_truth_dataset(chan,trial)==1
                                TP = TP + 1;
                                incorrect_response_data(chan,trial) = 0;
                            else
                                TN = TN + 1;
                                incorrect_response_data(chan,trial) = 0;
                            end
                        else
                            if ground_truth_dataset(chan,trial)==0
                                FP = FP + 1;
                                incorrect_response_data(chan,trial) = 1;
                            else
                                FN = FN + 1;
                                incorrect_response_data(chan,trial) = 1;
                            end
                        end
                        
                    end
                end
            end
            
            % Update rates depending on gain value used
            TP_rate(gain_ct,1) = 100*TP/(TP + FN);
            TN_rate(gain_ct,1) = 100*TN/(TN + FP);
            FP_rate(gain_ct,1) = 100*FP/(FP + TN);
            FN_rate(gain_ct,1) = 100*FN/(FN + TP);
            
            
        end
        %plot(gain_list,[FP_rate,TP_rate,FN_rate,TN_rate])
        %legend('FP rate','TP rate','FN rate','TN rate')
        %xlabel('gain');
        %ylabel('rate (%)');
        
        
        % Save run information to struct
        ratio_data(ct).window=windowSize;
        ratio_data(ct).step=stepSize;
        ratio_data(ct).gain_list = gain_list;
        ratio_data(ct).TP_rate = TP_rate;
        ratio_data(ct).FP_rate = FP_rate;
        ratio_data(ct).FN_rate = FN_rate;
        ratio_data(ct).TN_rate = TN_rate;
        ct=ct+1;
        
    end
end

fprintf("\nDone!\n")

% Generate final figure with ROC curves
figure(7)
hold on
for i=1%:length(ratio_data)
    plot(ratio_data(i).FP_rate,ratio_data(i).TP_rate);
end
hold off
Legend=cell(length(ratio_data),1);
for iter=1:length(ratio_data)
    Legend{iter}=['window: ',num2str(ratio_data(iter).window/0.03),', step: ',num2str(ratio_data(iter).step/0.03)];
end
legend(Legend)
title(['F',exp_data.cohort,' ROC']);
xlabel('FPR');ylabel('TPR');
%xlim([0 10]);ylim([70 100]);
xlim([0 100]);ylim([0 100]);

end

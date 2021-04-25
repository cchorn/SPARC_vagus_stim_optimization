function plot_trial(expmt_list, expmt, trial, avg_signal_data, baseline)
% Plot signal behavior for selected animal and trial
%
% Inputs
% =================================================================
% trial           : (num) trial selector within valid trial range
% avg_signal_data : (NxM cell) cell array containing segmented channal data across all channels (N) and trials (M)
% baseline        : (struct) struct containing mean and std baseline information across all channels and trials 


STA_t_min  = -2;
windowSize = 30;
stepSize   = 3;
fs         = 30e3;
plot_fig1  = true;
save_fig   = true;
rms_offset = 150;
data = expmt_list{expmt};
N_channels = length(data.elec_list);
N_trial = data.trial_list(end,end);
std_gain = data.std_gain;
excluded_trials = find_skipped_trials(expmt_list, expmt);
cohort = ['F',data.cohort];
session = find_session(expmt_list, expmt, trial);
pw = num2str(data.pulseWidth(session)*1000);
ses_trial = find_trial(expmt_list, expmt, trial);
stim = expmt_list{expmt}.stim_hist(session,ses_trial);
if ~ismember(trial,excluded_trials)
    
    for chan=1:N_channels
        
        % Collect segmented data
        STA_data = avg_signal_data{chan, trial};
        [~,N_data] = size(STA_data);
        STA_t_max = round(length(STA_data(abs(STA_t_min)*fs/1000:end))*1000/fs);
        % Create time period list for plotting data
        x_time = linspace(STA_t_min,STA_t_max,N_data);

        %calculate Vrms with one window/step size
        window = 1:windowSize;
        idx = 1;
        while max(window) < N_data
            Vrms_list(idx) = rms(STA_data(window));
            window = window + stepSize;
            %window = window + (windowSize-stepSize);
            idx = idx + 1;
        end
        window = window(1):length(STA_data);
        Vrms_list(idx) = rms(STA_data(window));
        Vrms_time = linspace(STA_t_min,STA_t_max,length(Vrms_list));

                
        % Find threshold by obtaining mean and std
        % from "baseline" struct, multiply std with gain
        % and add with mean to establish threshold value
        base_std = baseline{expmt,1}.(['Chan_',num2str(chan)]).std_trials(trial);
        base_mean = baseline{expmt,1}.(['Chan_',num2str(chan)]).mean_trials(trial);
        thresh = base_mean + std_gain*base_std;
        x_thresh = linspace(STA_t_min, STA_t_max,length(Vrms_list));
        threshold = linspace(thresh,thresh,length(x_thresh));
        
        % Check all instances of threshold crossings and evaluate whether a
        % crossing of 1 millisecond occurs denoting true response
        crossed = false;
        cross_dur = round(1000/(stepSize/0.03)); % convert milliseconds to elements
        for i=rms_offset:(length(Vrms_list)-rms_offset)
            if sum(Vrms_list(i:i+cross_dur)>thresh)>1%cross_dur
                crossed = true;
            end
        end
        if crossed==true
            detected_response_data(chan,trial) = 1;
            [max_val, max_idx] = max(Vrms_list(rms_offset:end-rms_offset));
            peak_time = (max_idx+STA_t_min)*1000/fs;
        else
            max_idx=0;
            detected_response_data(chan,trial) = 0;
        end
        
        if plot_fig1==true
            figure(6)
            sgtitle(sprintf('%s | Trial %d | Chan %d | PW%sus | %duA', cohort, trial, chan, pw, stim));
            subplot(1,2,1)
            %plot(x_time,new_signal_data{chan,trial}');
            hold on
            plot(x_time, STA_data, '-k','LineWidth',0.1);            
            hold off
            title('Filtered signal')
            ylim([-30,30]);
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
            ylim([-1,3])
            drawnow
            %if save_fig==true
            %    sgtitle(sprintf('Cohort %s | Trial %d | Channel %d | Spike: %s',strrep(cohort,'_','-'),ct,chan,spike))
            %    save_path = [pwd, '\F19_19_ground_truth_eval\', 'Trial',num2str(trial),'_Chan',num2str(chan)];
            %    saveas(gcf, [save_path '.png']);
            %end
            %keyboard
        end
        a = input('next chan: '); 
    end
else
    fprintf('Trial doesnt exist for this animal')
end


end

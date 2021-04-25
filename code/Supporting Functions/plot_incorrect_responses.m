function plot_incorrect_responses(expmt_list, expmt, clean_signal_data, baseline, plotting_array)

plot_fig1=true;
exp_data = expmt_list{expmt};
cohort = ['F',strrep(exp_data.cohort,'-','_')];
STA_t_min = -2;
STA_t_max = 498;
windowSize=15;
stepSize=12;
fs=30e3;
std_gain = 2.7;
rms_thresh_offset=40; %in step counts for RMS, not milliseconds



for trial=1:size(plotting_array,2)
    
    for chan=1:size(plotting_array,1)
        if plotting_array(chan,trial)==1
            snipData = clean_signal_data{chan, trial};
            STA = mean(snipData);
            
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
            
            % Find threshold by obtaining mean and std
            % from "baseline" struct, multiply std with gain
            % and add with mean to establish threshold value
            base_std = baseline{expmt,1}.(['Chan_',num2str(chan)]).std_trials(trial);
            base_mean = baseline{expmt,1}.(['Chan_',num2str(chan)]).mean_trials(trial);
            thresh = base_mean + std_gain*base_std;
            threshold = linspace(thresh,thresh,length(Vrms_list));
            x_thresh = linspace(STA_t_min, STA_t_max,length(Vrms_list));
            
            
            %Detect response
            if sum(Vrms_list(rms_thresh_offset:end-rms_thresh_offset) > threshold(rms_thresh_offset:end-rms_thresh_offset)) >= 1
                [max_val, max_idx] = max(Vrms_list(rms_thresh_offset:end-rms_thresh_offset));
                peak_time = (max_idx+STA_t_min)*1000/fs;
            else
                max_idx=0;
            end
            
            if plot_fig1==true
                figure(6)
                sgtitle(sprintf('Cohort %s | Trial %d | Channel %d |          ',strrep(cohort,'_','-'),trial,chan))
                subplot(1,2,1)
                plot(x_time, snipData')
                hold on
                plot(x_time, STA, '-k','LineWidth',0.1);
                hold off
                title('Filtered signal')
                ylim([-100,100]);
                xlim([STA_t_min,STA_t_max]);
                drawnow
                titleName = sprintf('RMS Window %d us, step %.3g us', windowSize/0.03, stepSize/0.03);
                subplot(1,2,2)
                plot(Vrms_list)
                hold on
                plot(threshold,'-r','LineWidth', 1);
                if max_idx~=0
                    plot((rms_thresh_offset+max_idx),threshold(1),'or')
                end
                hold off
                title(titleName)
                axis tight
                %xlim([2,400])
                ylim([-1,4])
                drawnow
                %keyboard
                choice = input('Press Enter to continue... ','s');
            end
        end
    end
end

end


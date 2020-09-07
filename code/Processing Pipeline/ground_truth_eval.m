function ground_truth_dataset = ground_truth_eval(data_table, expmt_list, expmt)
%   DESCRIPTION
%   ===================================================================
%   Function that loads trellis data for user analysis and generates matrix
%   of boolean values for stim responses
%
%   Jonathan Shulgach
%   Last updated: 4/7/2020
%
%   INPUTS
%   ===================================================================
%   expmt_list      :   (struct) data file header for experiments
%   expmt           :   (numeric) experiment number
%
%   OUTPUTS
%   ===================================================================
%   excluded_trials :   (1xN numeric) list of trials left out as
%                                   datafile numbers
%   NOTES
%   ===================================================================
%   Commenting out the function lines allows the ground_truth_dataset
%   matrix to remain intact in case an error persist 

exp_data = expmt_list{expmt};
N_channels = 32;
N_trials = exp_data.trial_list(end, 2);
windowSize=30;
stepSize = 15;
skip_trials = find_skipped_trials(expmt_list, expmt);
exp_data = expmt_list{expmt};
fs            = 30e3;                  % sampling frequency
STA_bin       = [-2 498];              %time window used for STA calculations where the stim event is at time 0
STA_t_min     = STA_bin(1)*fs/1000;    %convert ms to s, then multiply by sample rate
STA_t_max     = STA_bin(2)*fs/1000;
if exist('ct','var') == 0
    ct = 1;
end

% Collect Data into one cell table
signal_data=[];
cohort = ['F',strrep(exp_data.cohort,'-','_')];
for i=1:6
    pair = strrep(exp_data.cuff_list(i,:),':','_');
    pw = num2str(exp_data.pulseWidth(i)*1000);
    signal_data = [signal_data, data_table.(cohort).(['Cuff_',pair,]).(['PW_',pw,'_us'])];
end

% Add cell columns of zeros for skipped trials to create the full
% 32xN_trials cell table
for skip_col=skip_trials
    if skip_col==1
        signal_data = [cell(32,1), signal_data];
    else
        signal_data = [signal_data(:,1:skip_col-1), cell(32,1),signal_data(:,skip_col:end)];
    end
end

[clean_signal_data, cleaned_trials] = remove_artifacts(signal_data,100);


for trial=1:N_trials
    for chan=1:N_channels
        % Check whether trial is included in the list of trials
        % to skip, and if so... skip and replace all detected
        % trial responses to zero        
        if ismember(trial, skip_trials)
            ground_truth_dataset(chan,ct) = 0;
        else
            
            % Plot signal snippets and STA data onto figure 1
            snipData = clean_signal_data{chan, trial};
            STA = mean(snipData);
            x_time = linspace(STA_t_min,STA_t_max,size(snipData,2))*1000/fs;
            figure(4)
            sgtitle(sprintf('Cohort %s | Trial %d | Channel %d |          ',strrep(cohort,'_','-'),trial,chan))
            subplot(1,2,1)
            plot(x_time, snipData')
            hold on
            plot(x_time, STA, '-k','LineWidth',0.1);
            hold off
            title('Filtered signal')
            xlim([STA_bin(1),STA_bin(2)])
            ylim([-100,100])
            drawnow
            
            %calculate Vrms with one window/step size
            window = 1:windowSize;
            idx = 1;
            while max(window) < length(STA)
                Vrms_list(idx) = rms(STA(window));
                window = window + stepSize;
                %window = window + (windowSize-stepSize);
                idx = idx + 1;
            end
            
            % Plot RMS data onto figure 2
            titleName = sprintf('RMS Window %d us, step %.3g us', windowSize/0.03, stepSize/0.03);
            subplot(1,2,2)
            plot(Vrms_list)
            hold on
            x_thresh = linspace(0,length(Vrms_list),length(Vrms_list));
            thresh_temp = linspace(0.86,0.86,length(Vrms_list));
            plot(x_thresh,thresh_temp,'r')
            hold off
            title(titleName)
            axis tight
            ylim([-1,4])
            drawnow
            
            % Prompt user for choice and save response to ground truth
            % dataset
            choice = input('Response or not? (y or n) ','s');
            if strcmpi(choice,'y')==true
                ground_truth_dataset(chan,trial) = 1;
                spike='YES';
            else
                ground_truth_dataset(chan,trial) = 0;
                spike='NO';
            end
            
            % Save figures for future analysis
            %sgtitle(sprintf('Cohort %s | Trial %d | Channel %d | Spike: %s',strrep(cohort,'_','-'),ct,chan,spike))
            %save_path = [pwd, '\F19_19_ground_truth_eval\', 'Trial',num2str(trial),'_Chan',num2str(chan)];
            %saveas(gcf, [save_path '.png']);
        end
    end
end

end

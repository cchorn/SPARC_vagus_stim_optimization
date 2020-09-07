function create_baseline_figure(expmt_list, expmt, baseline)
% Plot new baseline figure with only valid trials
%
% Written by: Jonathan Shulgach
% Updated: 5/8/20
%  
% Inputs
% =================================================================
% baseline        : (struct) struct containing mean and std baseline information across all channels and trials 



%determine constants
exp_data = expmt_list{expmt}; %extract meta data for the desired animal
%stim_port = exp_data.stim_port; %extract the stim port value
cohort = sprintf('F%s', exp_data.cohort); %the animal identifier
N_electrodes = 32;
%rec_path = fullfile(source_path, cohort, 'Trellis'); %path to the folder with data
N_trials = expmt_list{expmt,1}.trial_list(6,2); %last trial number to be analyzed
%break_point = expmt_list{expmt,1}.trial_list(3,2); %the trial number where stimulation switches from 1:2 to 3:4
fs=30e3; %sampling frequency
STA_bin = [-2 400]; %window for STA calculation where stim event is 0
STA_t_min = STA_bin(1)*fs/1000; %convert ms to s, then multiply by sample rate
STA_t_max = STA_bin(2)*fs/1000;
%window_size = STA_t_max - STA_t_min; %the number of sampled points in the appropriate window size
%buffer = 3;% 1.5; %time (ms) around the stimulation event to be removed
%FailedTrials = zeros(32,N_trials); %indicates which trials were not successfully processed
sub_size = numSubplots(N_electrodes);

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

for chan = 1:N_electrodes
    figure(17)
    sgtitle(sprintf('%s Baseline STD',cohort))
    grid_fig = subplot(sub_size(1),sub_size(2), elec_chan_map(chan));
    ct = 1;
    temp = baseline{expmt,1}.(['Chan_',num2str(chan)]).std_trials;
    for trial = 1:N_trials
        if temp(trial)~=0
            baseline_y(ct) = temp(trial);
            baseline_x(ct) = trial;
            ct = ct + 1;
        end
    end
    plot(baseline_x,baseline_y)
    ylim([0 1.5])
    if ~mod(chan,4)==0
        set(gca,'XTick',[])
    else
        xlabel('Trial');
    end
    if chan>4
        set(gca,'YTick',[])
    else
        ylabel('\muV')
    end
        
end
end
%
% Analysis script that goes through the entire processing pipeline in
% segregated processes, allowing for intermmediate savepoints and updating
% old data.
%
% Comment out the portions of the pipeline not relevant to your analysis. 

% =======================================================================ZAZZ
% ======== Step 1: Choose animal =========
% Animals are analyzed in the chronological order of the experiments. 
% 1 : F20-19, 2 : F21-19, 3 : F19-19, 4: F25-19, 6 : F26-19, 7 : F34-19
%
% Note: Because animal 1 (F20-19) has only 2PW tested with one cuff pair,
% analysis on this animal is ignored, and performed on 2-7--------
% Ex:  expmt = 3 : F22-19 
expmt = 2;


%% =======================================================================
% ===== Step 2: Configure Data Path ======
% Load the header file containing all the wxperiment information for each
% animal, and set the source path to the data files 
load('expmt_list.mat') % Will load the "expmt_list.m" header file containing all expeirment and stim params info
source_path = 'C:\Users\jonat\Box\Data_for_Dylan'; % Box cloud repo path where all animal folders containingtrellis recordings are saved to Box
%source_path = 'C:\Users\jonat\Documents\RNEL\source'; %local repo

%% =======================================================================
% ========= Step 3: Collect Baseline Data =============
% Acquire response thresholds for each channel and trial by analyzing the
% baseline data
if exist('baseline','var') == false
    %baseline = visualizeNoise(expmt_list, expmt, source_path);
    %create_baseline_figure(expmt_list, expmt, baseline)
end
% In case artifacts infiktrate the recordings, analyze the baseline data to
% check values are realistic (e.g. sanity check). If channels or trials are
% too noisy, remove them from the analysis
if exist('bad_trials','var') == false
    %fprintf("Checking baseline values\n====\n")
    %manual_trial_skip = [];
    %manual_chan_skip = []; % Caution: requires manual checking of channels. 
    %[good_trials, bad_trials, trial_results] = analyze_baseline_values(baseline, expmt);
    if exist('excluded_trials','var') == false
        % Excluded trials have now been saved to header file, so reqriting is
        % not necessary
    
        %excluded_trials = expmt_list{expmt}.exclude_trials;
        %excluded_trials = find_skipped_trials(expmt_list, expmt);
    end
    %expmt_list{expmt,1}.exclude_trials = sort([excluded_trials, bad_trials, manual_trial_skip]);
    %expmt_list{expmt,1}.exclude_channels = manual_chan_skip;
end

%% =======================================================================
% ================ Step 4: Load raw data and collect STA ===============
% Read through recording files, perform filtering and segmentation (no STA
% yet)
% Note: Skip this step if you already have STA data collected
if exist('data_table','var') == false
    % data_table = segment_raw_data(expmt_list, expmt, source_path);
end

% Collect the mean of segmented data for each channel and trial to get STA
% data table
if exist('avg_signal_data','var') == false
%    fprintf("collecting average values from segmented data\n====\n")
%    [new_signal_data, avg_signal_data] = convert_seg_to_table(expmt_list, expmt, data_table);
end


%% ======================================================================
% ========== Step 5: Analyze responses in Data ============
% The expmt_list will become poppulated with the response data across
% trials for each channel. If the RMS of the signal crosses the threshold,
% a response is recorded.
%
% Note: Steps 2-4 can be skipped if you saved the intermmediate baseline and 
% avg_signal_data files beforehand 
tic
expmt=2;
load('C:\Users\YATESLAB\Documents\MATLAB\F21_19_average_signal_data_8-6-20.mat')
load('C:\Users\YATESLAB\Documents\MATLAB\F21_19_baseline_4-28-20.mat')
expmt_list = analyze_stim_data_final(expmt_list, baseline, avg_signal_data, expmt, true, 50);
expmt=3;
load('C:\Users\YATESLAB\Documents\MATLAB\F22_19_baseline_4-28-20.mat')
load('C:\Users\YATESLAB\Documents\MATLAB\F22_19_average_signal_data_7-27-20.mat')
expmt_list = analyze_stim_data_final(expmt_list, baseline, avg_signal_data, expmt);
expmt=4;
load('C:\Users\YATESLAB\Documents\MATLAB\F19_19_baseline_4-28-20.mat')
load('C:\Users\YATESLAB\Documents\MATLAB\F19_19_average_signal_data_4-28-20.mat')
expmt_list = analyze_stim_data_final(expmt_list, baseline, avg_signal_data, expmt);
expmt=5;
load('C:\Users\YATESLAB\Documents\MATLAB\F25_19_average_signal_data_4-28-20.mat')
load('C:\Users\YATESLAB\Documents\MATLAB\F25_19_baseline_4-28-20.mat')
expmt_list = analyze_stim_data_final(expmt_list, baseline, avg_signal_data, expmt);
expmt=6;
load('C:\Users\YATESLAB\Documents\MATLAB\F26_19_baseline_4-28-20.mat')
load('C:\Users\YATESLAB\Documents\MATLAB\F26_19_average_signal_data_4-28-20.mat')
expmt_list = analyze_stim_data_final(expmt_list, baseline, avg_signal_data, expmt);
expmt=7;
load('C:\Users\YATESLAB\Documents\MATLAB\F34_19_baseline_4-28-20.mat')
load('C:\Users\YATESLAB\Documents\MATLAB\F34_19_average_signal_data_4-28-20.mat')
expmt_list = analyze_stim_data_final(expmt_list, baseline, avg_signal_data, expmt);
toc



%% ======================================================================
%  ========== Step 6: Analyze Response Activity ==========
% Count all responses and compare MEA performance across all PW
% combinations between cuff pairs

[SI_combos_list, overlap_combos_list] = SI_search(expmt_list, 2:7);

%% ======================================================================
% ========== Step 7: Analyze Sepective Response Behavior =========
% 
% Note: May need to set patrameters inside function to plot desired figures

% Generate plots to show the optimal stimulation parameter settings to get
% the most "selective" responses, using the "Selectivity Index" function
% output. % PW####.max_SI.max_SI_selec_chan contains 90% selectivity results
opt_params_list = get_opt_params(SI_combos_list, expmt_list, 2:7, 3, true);

plot_response_distribution(expmt_list, 2:7, opt_params_list);

% Generate plots to show selective boundaries with x-axis and y-axis as
% stim amplitude for each cuff pair
[avg_line_data, expmt_list] = get_selec_stim_bound(SI_combos_list, expmt_list);




%% ======================================================================
% ========== Step 8: Recruitement Curves ===============
% Generate cuff threshold RC curves with x-axis as PW duration and y-axis
% as stim amplitude
min_thresh_plot(expmt_list)


%% ======================================================================
% ========== Step 9:  Analyze Conduction Velocity =========
% Generate cv figures showing response behavior at max stim and cuff
% threshold
[cv_min_data, cv_max_data] = get_conduction_velocity(expmt_list, 2:7, true);
cv_table = organize_cv_table(expmt_list, cv_min_data);



%%======================================================================
%Other functions below are used to assist in the data collection/validation and are otherwise not essential for the processing pipeline
%%======================================================================
 
% ======= Ground Truth Datasets ======
%load('F22_19_ground_truth_dataset.mat')

%if exist('ground_truth_dataset','var') == 0
    % Collect ground truth from signals, with final cell table as
    % "ground_truth_dataset". This can be used with "false_alarm_ratio" to
    % find the optimal settings for threshold gain.
    
%end
%if exist('ratio_data','var') == 0
%    fprintf("conducting false alarm ration analysis\n====\n")
%    [ratio_data, detected_response_data, incorrect_response_data] = false_alarm_ratio(expmt_list, baseline, avg_signal_data, ground_truth_dataset, expmt, bad_trials);
%    plot_incorrect_responses(expmt_list, expmt, raw_signal_data, baseline, incorrect_response_data)
%end



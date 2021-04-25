function [threshold, x_threshold] = getResponseThreshold(varargin)
% Grabs the std value and uses the animal specific gain and optional manual
% offset to establish the threshold value for particular channel 
%
% Note: manual offset is only added if manual_std_ovveride is set to tru
% for specific channel and trial

% Written By: Jonathan Shulgach
% Last updated: 8/26/20


expmt_list = varargin{1};
expmt = varargin{2};
baseline = varargin{3};
trial = varargin{4};
chan = varargin{5};

std_gain = expmt_list{expmt}.std_gain;

% Find threshold by obtaining mean and std
% from "baseline" struct, multiply std with gain
% and add with mean to establish threshold value
base_std = baseline{expmt,1}.(['Chan_',num2str(chan)]).std_trials(trial);
base_mean = baseline{expmt,1}.(['Chan_',num2str(chan)]).mean_trials(trial);
threshold = base_mean + std_gain*base_std;

if length(varargin)>5
    N_data = varargin{6};
    % Create array with length of data with values of threshold (used for plotting)
    if length(varargin)>6 
        % If choosing to personally add a manual offset for this channel and trial
        threshold = threshold + varargin{7};
    elseif  expmt_list{expmt}.manual_std_override(chan,trial)==true
        % If predefined for this channel and trial
        threshold = threshold + expmt_list{expmt}.manual_std_offset(chan,trial);
    end
    x_threshold = linspace(threshold,threshold,N_data);
else
    x_threshold = [];
    
end

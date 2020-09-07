function [good_trials, bad_trials, trial_results] = analyze_baseline_values(baseline, expmt)
%   DESCRIPTION
%   ===================================================================
%   Function that checks baseline recordings from "visualizeNoise" and
%   makes suggestion of trials to further exclude from analysis
%
%   Written by: Jonathan Shulgach
%   Last updated: 4/7/2020
%
%   INPUTS
%   ===================================================================
%   baseline      :   (1xN cell) cell array holding baseline information
%                           across all trials and channels
%   expmt         :   (numeric) experiment number
%
%   OUTPUTS
%   ===================================================================
%   good_trials   :   (1xN numeric) list of trials with std below threshold
%   bad_trials    :   (1xN numeric) list of trials with std above threshold
%   trial_results :   (1xN cell) cell array of trials with lists of
%                         channels having std above threshold

good_trials = [];
bad_trials = [];
N_trials = size(baseline{expmt,1}.Chan_1.mean_trials,2); % Assuming there a Chan_1
for trial=1:N_trials
    bad = [];
    for chan_num=1:32
        
        % Collect std value for channel and evaluate whether it crosses
        % threshold
        temp = baseline{expmt,1}.(['Chan_',num2str(chan_num)]).std_trials(trial);
        if temp > 2 && ~ismember(temp,bad)
            bad = [bad,temp];
        end
    end
    
    % Store list of good and bad trials
    if isempty(bad)
        good_trials = [good_trials, trial];
    else
        bad_trials = [bad_trials, trial];
    end
    trial_results{trial} = bad;
end

end

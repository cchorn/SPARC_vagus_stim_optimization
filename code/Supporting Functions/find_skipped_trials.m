function excluded_trials = find_skipped_trials(expmt_list, expmt)
%   DESCRIPTION
%   ===================================================================
%   Function that find trial/datafile numbers not included in actual trial list
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

exp_data = expmt_list{expmt};
actual_trials=[];

% First create a list of all actual trials run
for session=1:size(exp_data.trial_list, 1) % for each trial list session (usually 6)
    start_t = exp_data.trial_list(session, 1);
    end_t = exp_data.trial_list(session, 2);
    actual_trials = [actual_trials, start_t:end_t];
end

% Find largest trial number run
N_trial = exp_data.trial_list(end, 2);

% The function "ismember" will output a boolean list of size(1xN_trial) with zeros 
% for any numbers missing in "actual_trials" indexed from 1 to N_trial. Negating this output 
% produces ones in the missing trials, multiplying this output by element with the same 
% 1:N_trial list and removing remaining zeros produces a list of only
% missing trials
excluded_trials = nonzeros((1:N_trial).*~ismember(1:N_trial,actual_trials))';

end

function pw = find_pw(meta_data, trial_num)
%Parses through the meta data to identify the pulsewidth for a specified
%trial number

%meta_data: a struct containing information about the specific experiment
%trial_num: the trial number the pulsewidth is being found for

%   Made by: Dylan Beam on 2/10/2020
%   Last edited by: Dylan Beam on 2/10/2020

%found out which binary search the trial is nested in
if((trial_num >= meta_data.trial_list(1,1)) && (trial_num <= meta_data.trial_list(1,2)))
    pw = meta_data.pulseWidth(1);
elseif ((trial_num >= meta_data.trial_list(2,1)) && (trial_num <= meta_data.trial_list(2,2)))
    pw = meta_data.pulseWidth(2);
elseif ((trial_num >= meta_data.trial_list(3,1)) && (trial_num <= meta_data.trial_list(3,2)))
    pw = meta_data.pulseWidth(3);
elseif ((trial_num >= meta_data.trial_list(4,1)) && (trial_num <= meta_data.trial_list(4,2)))
    pw = meta_data.pulseWidth(4);
elseif ((trial_num >= meta_data.trial_list(5,1)) && (trial_num <= meta_data.trial_list(5,2)))
    pw = meta_data.pulseWidth(5);
elseif ((trial_num >= meta_data.trial_list(6,1)) && (trial_num <= meta_data.trial_list(6,2)))
    pw = meta_data.pulseWidth(6);
else 
    error('Invalid trial number; this trial does not exist or was excluded for this animal')
end


end


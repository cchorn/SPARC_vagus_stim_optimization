function session_trial = find_trial(expmt_list, expmt, trial)
% This function parses through the data to find the trial number where 
% trial 1 is the first trial for the associated session, not the first 
% trial in the experiment. 
%
%   Jonathan Shulgach
%   Last updated: 4/14/2020
%
%

exp_data = expmt_list{expmt};
session_trial  =[];
for i=1:size(exp_data.trial_list,1)
    t = exp_data.trial_list(i,:);
    if trial>=t(1) && trial<=t(2)
        trial_list = t(1):t(2);
        session_trial = find(trial_list==trial);
        break
    else
        
        %error('Invalid trial number; this trial does not exist or was excluded for this animal')
    end
end
    
end

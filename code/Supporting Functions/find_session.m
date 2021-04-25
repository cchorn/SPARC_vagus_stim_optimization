function session = find_session(expmt_list, expmt, trial)
% This function parses through the data to find the session number 
% associated with the experiment and trial number.
%
%   Jonathan Shulgach
%   Last updated: 7/21/2020
%

exp_data = expmt_list{expmt};
session  =[];
for i=1:size(exp_data.trial_list,1)
    t = exp_data.trial_list(i,:);
    if trial>=t(1) && trial<=t(2)
        session = i;
        break
    else
        
        %error('Invalid trial number; this trial does not exist or was excluded for this animal')
    end
end

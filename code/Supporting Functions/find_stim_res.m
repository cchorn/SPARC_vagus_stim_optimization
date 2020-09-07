function stim_res_list = find_stim_res(expmt_list)
% stimulation difference condition was set at less than 20uA, but this
% program should verify for each channel what the minimum difference in
% stimulation amplitude was for each channel between threshold and the
% next closest value

for expmt=1:length(expmt_list) % for each day
    exp_data = expmt_list{expmt};
    if size(exp_data.trial_list,1) > 2 % avoid data with less than 2 pw points
        for i=1:6
            stim_list = unique(expmt_list{expmt,1}.stim_hist(i,:));
            for chan=1:length(exp_data.elec_list)
                %check value is not NaN
                if isnan(exp_data.minthresh(chan,i))
                    stim_res_list{expmt,1}.thresh_res(chan,i) = exp_data.minthresh(chan,i);
                else 
                    % Get list of stim vals and find range of min thresh index +/- one
                    min_thresh_idx = find(stim_list==exp_data.minthresh(chan,i));
                    if min_thresh_idx == 1
                        stim_res_list{expmt,1}.thresh_res(chan,i) = abs((stim_list(min_thresh_idx) - stim_list(min_thresh_idx+1)));
                    else
                        stim_res_list{expmt,1}.thresh_res(chan,i) = abs((stim_list(min_thresh_idx) - stim_list(min_thresh_idx-1)));
                    end
                end
            end
            stim_res_list{expmt,1}.cohort = exp_data.cohort;
            stim_res_list{expmt,1}.pulseWidth = exp_data.pulseWidth;
            stim_res_list{expmt,1}.smallest_stim_res(i) = min(stim_res_list{expmt,1}.thresh_res(:,i));
        
            histogram(stim_res_list{expmt,1}.thresh_res(:,i), 100)
            xlim([0 200]);
            title([exp_data.cohort, ' | stim resolution distribution'])
            xlabel('Stim amplitude');
            ylabel('Frequency');
        end
    end
end

end

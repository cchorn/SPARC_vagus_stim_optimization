function [new_signal_data, avg_signal_data] = convert_seg_to_table(expmt_list, expmt, data_table)
%   DESCRIPTION
%   ===================================================================
%   Function that creats a cell table of trial data instead of the struct
%   formatas well as a cell table of STA values
%
%   Written by: Jonathan Shulgach
%   Last updated: 8/5/2020
%
%   INPUTS
%   ===================================================================
%   expmt_list    :   (1x1 struct) header for experiment data
%   expmt         :   (numeric) experiment selection number
%   data_table    :   (1x1 struct) output from "segment_raw_data" function
%
%   OUTPUTS
%   ===================================================================
%   raw_signal_data   :   (MxN numeric) cell matrix with M channels and N
%                               trials
%   avg_signal_data   :   (MxN numeric) cell matrix with M channels and N
%                               trials


exp_data = expmt_list{expmt};

skip_trials = find_skipped_trials(expmt_list, expmt);

% Collect Data into one cell table
raw_signal_data=[];
cohort = ['F',strrep(exp_data.cohort,'-','_')];
for i=1:6
    pair = strrep(exp_data.cuff_list(i,:),':','_');
    pw = num2str(exp_data.pulseWidth(i)*1000);
    raw_signal_data = [raw_signal_data, data_table.(cohort).(['Cuff_',pair,]).(['PW_',pw,'_us'])];
end

% Add cell columns of zeros for skipped trials to create the full
% 32xN_trials cell table
for skip_col=skip_trials
    if skip_col==1
        raw_signal_data = [cell(32,1), raw_signal_data];
    else
        raw_signal_data = [raw_signal_data(:,1:skip_col-1), cell(32,1),raw_signal_data(:,skip_col:end)];
    end
end
new_signal_data = raw_signal_data;
%[new_signal_data, ~, ~] = replace_missing_data_packets(raw_signal_data, 1050);
%[new_signal_data, ~] = remove_artifacts(raw_signal_data,150);
% Get averages from data snippets and create new cell matrix
fprintf("Averaging Signal | Trial: ")
for trial=1:size(new_signal_data,2)
    fprintf([num2str(trial),'|'])
    for chan=1:size(new_signal_data,1)
        temp = new_signal_data{chan, trial};
        avg_signal_data{chan,trial} = mean(temp);
    end
end

end


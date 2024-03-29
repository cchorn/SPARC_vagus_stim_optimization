function [crossed, max_val, max_idx, time_bins, resp_bins] = response_detection(exp_data, threshold, V_data, V_time, rms_offset, bin_time, varargin)
% This function determines if a response to stimulation happened, and records a series of time bins where
% repsonses took place. Conduction velocity is measures in 0.5m/s increments, starting with the first bin
% (0-0.5m/s) between the end of the data and the first time marker found from:
%                { Response_Time = Cuff_Array_Distance/Conduction_Velocity }
%
% exp_data : header information
% thresh   : (int) threshold crossing value for detecting stim response
% V_data   : (1xN) signal data
% V_time   : (1xN) time indices of signal data
% offset   : (int) time indices offset to start analyzing threshold crossing (compensate for stim events)

% rms_offset format: 150elements which is 15ms
% Simple classification: For STA check any threshold crossinginstance and evaluate whether a
% crossing of 1 millisecond occurs denoting true response. For
% RMS just see if it crosses at least once
verbose=false;
scale=1; % keep scale at 1 for Vrms
fs=30e3;
crossed = false;
cross_dur = 1;%round(1000/(stepSize/0.03)); % convert milliseconds to elements
exp_data;
if ~isempty(varargin) && varargin{1}==true
    % extra offset towards end of signal in case stim artifacts exist (needed for F21-19 Trial49) 
    rms_end_offset = round(varargin{2});
else
    rms_end_offset = 0;
end
if sum(V_data(rms_offset:end-rms_offset-rms_end_offset)>threshold)>=cross_dur
    crossed = true;
    %bin_time = 20; %ms
    [max_val, max_idx] = find_peaks(V_data(rms_offset:end-rms_offset-rms_end_offset),threshold, fs, bin_time);
    %[max_val, max_idx] = max(V_data(rms_offset:end-rms_offset));
else
    max_val = 0;
    max_idx = 0;
end


% Crossing points: Check whether crossing occurs within time bins
L_nerve = exp_data.nerve_length;
cv = linspace(0.5,30,60); % initialize an array of 0-to-30 cv values, stepped at 0.5m/s

% for 0.5m/s conduction velocity steps, get the actual time indices
bin = L_nerve./cv;
%time_bins = round(1000*bin); % unique time bins rounded to ms
time_bins = 1000*bin;
% keep in mind we also need to add the range between the first bin value and the end
time_bins = [V_time(end-(rms_offset/scale)), time_bins];

% Find the first occurences of times in milliseconds bc the signal has more
% precision than milliseconds
temp_time = floor(10*V_time)/10; % round all times
[C,ia] = unique(temp_time); % find time indices and positions
% better to flip time indices because of 1) descending bin size and 2)
% avoiding the blanking period
t_idx = flip(ia(C>0));
t_val = flip(C(C>0)); % should match index with value (ex:1=1,2=2,...498=498)

% Vrms_list has values from 0-499ms but has ~5000 elements, so need to
% include 10 elements for every 1ms
resp_bins = linspace(0,0,length(time_bins));
if max_idx~=0
    actual_peak_time = V_time(rms_offset+max_idx-1);
    %actual_peak_time_idx = rms_offset/10 + V_time(max_idx);
else
    actual_peak_time = 0;
end
for t=1:length(time_bins)-1
    
    if t==52
        a=1;
    end
    
    % remember times are in descending order (to collect cv from 0 and higher),
    % so start time indexed ahead (higher t==lower cv)
    t_start = floor(10*time_bins(t+1))/10; % time bin range begin (ex:196)
    t_end = floor(10*time_bins(t))/10; % time bin range end   (ex:392)
    
    % First we need to match the time values with t_bin with the indices of
    % those same times in the data
    t_start_idx = t_val==t_start; % get index of time bin mark WRT Vrms time
    rms_t_start = t_idx(t_start_idx);
    t_end_idx = t_val==t_end;
    rms_t_end = t_idx(t_end_idx);
    
    % grab all Vrms values between time bins
    if sum(V_data(rms_t_start:rms_t_end)>threshold)>=1
        % Make sure peaks detected have a window around them where
        % no othe rpeaks will get detected. This is useful for peaks detected
        % between 15-30m/s where the spacing differes by 0.1-10millis.
        if sum((V_time(rms_t_start) < actual_peak_time).*(V_time(rms_t_end) >= actual_peak_time))==1
            resp_bins(t) = 1;
        else
            resp_bins(t) = 0;
        end
    else
        resp_bins(t) = 0;
    end
end
% Shorten bin data if it exceeds the rms offset, since we don't want to
% include false positive responses from stim events
%while time_bins(end) < V_time(rms_offset) % 15ms offset, consistent across all animals regardless of peak analysis offset
%    time_bins = time_bins(1:end-1);
%    resp_bins = resp_bins(1:end-1);
%end



end

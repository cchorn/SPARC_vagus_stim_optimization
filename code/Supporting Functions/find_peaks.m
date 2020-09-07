function [spike_amp_data, spike_amp_t, bin_ct] = find_peaks(data, threshold, fs, bin_time)
% This function searches for local peaks separated by time bins
% Written by: Jonathan Shulgach
% Last Updated: 8/11/20

all_bin_peaks=false; %enable to collect all threshold crossing peaks of the signal without binning

spike_amp_t = [];
spike_amp_data = [];
bin_size = bin_time*fs/10;

% Find all the points in the signal where the signal crosses threshold
% which amount to detecting a spike.
if threshold > 0
    aboveThreshold = (data > threshold);
else
    aboveThreshold = (data < threshold);
end

aboveThreshold = [false, aboveThreshold, false];  %pad with 0's at ends
edges = diff(aboveThreshold);
temp_rising = find(edges==1);     %rising/falling edges
temp_falling = find(edges==-1);

% Check indices of falling does not pass length of data
if ~isempty(temp_falling)
    while temp_falling(end) > length(data)
        temp_falling = temp_falling(1:end-1);
    end
end

if all_bin_peaks==false
    % check if any following spike times fall within the bin after current spike
    [rising, falling] = check_spike_overlap(temp_rising, temp_falling, bin_size);
end

% Find and store peaks between rising/falling periods
for i=1:length(falling)
    if threshold > 0
        [val,pos] = max(data(rising(i):falling(i)));
    else
        [val,pos] = min(data(rising(i):falling(i)));
    end
    
    temp_t = rising(i) + pos-1;
    spike_amp_t = [spike_amp_t, temp_t];
    spike_amp_data = [spike_amp_data, val];
end


t = 0;
bin_ct = [];
% Partition the number of spikes detected according to specific time bins
while t <= length(data) + bin_size
    temp_ct = length(spike_amp_t((spike_amp_t < (t+bin_size)) & (spike_amp_t > t)));
    if isempty(temp_ct)
        temp_ct = 0;
    end
    bin_ct = [bin_ct, temp_ct];
    t = t + bin_size;
end

% Double check we have the right amount of bins for the data
while length(bin_ct) > floor(length(data)/(fs))
    bin_ct = bin_ct(1:end-1);
end

end

function [updated_rising, updated_falling] = check_spike_overlap(temp_rising, temp_falling, bin_size)

updated_rising = temp_rising;
updated_falling = temp_falling;
done = false;

while done==false
    if length(updated_rising)>1 % If only one value left, stop
        small_intervals = diff(updated_falling) < bin_size;
        if sum(small_intervals) > 0 % If there are no remaining indices smaller than the time bin, stop
            [~,idx] = find(small_intervals==1);
            % Just remove the second index, then check interval spacing
            % again and again until all that's left in the array are spaced
            % intervals
            updated_falling(idx(1)+1) = nan; 
            updated_rising(idx(1)+1) = nan;
            updated_rising = updated_rising(~isnan(updated_rising));
            updated_falling = updated_falling(~isnan(updated_falling));
        else
            done=true;
        end
    else
        done=true;
    end
end


end

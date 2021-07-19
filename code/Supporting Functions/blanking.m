function [blankedData, blankedOverlay] = blanking(data, stimEvents, pw, fs, buffer)
%Blanks out stimulation events for the input data with the indicated list of
%stimulation events, pulse width, sampling freq and the desired buffer size

%   Written by: Dylan Beam on 3/19/2020
%   Last edited by: Dylan Beam on 3/23/2020
%   data: the data stream input for blanking
%   stimEvents: where stimulation events start in the data vector
%   pw: pulsewidth of the stimuli (ms)
%   fs: the sample frequency (Hz)
%   buffer: desired buffer to be blanked on either side of the stimuli (ms)


%identify the size of windows to be cut
window_size = (2*pw)/1000; %2x assumes same pulsewidth for the hyper and depolarization phases; 1000 converts ms to s

%transform window_size and buffer into sampling space
window_size_sampled = ceil(window_size*fs);
buffer_sampled = ceil(buffer*fs/1000);
blankedLength = window_size_sampled + 2*buffer_sampled; 

%create an array to index the start of the stim event and all the other
%data points within the window
window = [];
blankedOverlay = cell(length(stimEvents), 1);

%loop through all stimulation events
for i = 1:length(stimEvents)
        
    %identify the indices that should be blanked
    window(1) = stimEvents(i)-buffer_sampled;
    window(2) = stimEvents(i)+window_size_sampled+buffer_sampled;
    
    %identify the length of the window for linear interpolation
    windowLength = window(2) - window(1) + 1;
    
    %blank the data
    data(window(1):window(2)) = linspace(data(window(1)-1), data(window(2)+1), windowLength);
        
    %create a matrix of all the data so it can be super imposed on top of
    %one another
    if i == length(stimEvents)
        blankedOverlay{i,1} = data(window(1)-1:end);
    else
        blankedOverlay{i,1} = data(window(1)-1:stimEvents(i+1)-buffer_sampled-2);
    end
    
end

%create blanked data stream
blankedData = data;

end


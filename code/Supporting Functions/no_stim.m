function [data_transformed] = no_stim(data, stim_indices, pulse_width, fs, buffer)
%Eliminates stimulation artifacts from data

%data: array of sampled data to be processed
%stim_indices: an array of numbers referring to the indices in the raw data
%of which stimulation artifacts are initially detected
%pulse_width: the length of a pulse (in ms)
%fs: sampling frequency (in Hz)
%buffer: time (ms) around each side of the stimulation artifact that you would
%like to eliminate

%   Made by: Dylan Beam on 2/10/2020
%   Last edited by: Dylan Beam on 2/10/2020

%identify the size of windows to be cut
window_size = (2*pulse_width)/1000; %2x assumes same pulsewidth for the hyper and depolarization phases; 1000 converts ms to s

%transform window_size and buffer into sampling space
window_size_sampled = ceil(window_size*fs);
buffer_sampled = ceil(buffer*fs/1000);

%create an array to index the start of the stim event and all the other
%data points within the window
ID_w_window = [];
for i = 1:length(stim_indices)
    
    lengthID2 = length(ID_w_window);
    
    ID_w_window(lengthID2+1:lengthID2+1+(window_size_sampled+2*buffer_sampled)) = stim_indices(i)-buffer_sampled:stim_indices(i)+window_size_sampled+buffer_sampled;
    
end

%eliminate data within the appropriate window
data(ID_w_window) = [];

data_transformed = data;

end


function [clean_data, cleaned_trials] = remove_artifacts(varargin)
%   DESCRIPTION
%   ===================================================================
%   Function that removes artifacts in data by implementing linear interpolated blanking
%   during the period where the signal crosses a threshold and clean the
%   signals
%
%   Jonathan Shulgach
%   Last updated: 5/3/2020
%
%   INPUTS
%   ===================================================================
%   signal_data :  (MxN cell) Cell array with raw data signal snippets
%                           for M channels and N trials
%   thresh      :  (numeric) numeric threshold value for artifact detection
%   offset      :  (numeric) offset to begin threshold crossing detection 
%   fs          :  (numeric) sampling rate
%
%   OUTPUTS
%   ===================================================================
%   clean_data        :  (MxN cell) Cell array with cleaned average data signals
%                           for M channels and N trials
%
%   EXAMPLE
%   ===================================================================
%   
%   >> [clean_data, cleaned_trials] = remove_artifacts(raw_signal_data, 200)


%Set default parameterd
plot_fig=false;

signal_data = varargin{1};


cleaned_trials = [];
if iscell(signal_data)==1
    clean_data = cell(size(signal_data));
else
    temp = signal_data;
    signal_data = 0;
    signal_data = cell(1);
    signal_data{1,1} = temp;
    clean_data = [];
end

% Load signal data and paramters, otherwise set default values if no additional inputs 
if length(varargin) > 2
    thresh = abs(varargin{2});
    offset = abs(varargin{3});
elseif length(varargin) > 1
    thresh = abs(varargin{2}); % microvolts
    offset = 91;
else
    thresh = 200;
    offset = 91;
end

if length(varargin)>3
    fs = varargin{4};
else
    fs=30e3;
end

% Loop through each trial (column) and channel (row) in cell array 
fprintf("Searching for Artifacts | Trial: ")
for trial=1:size(signal_data,2)
    fprintf([num2str(trial),'|'])
    for chan=1:size(signal_data,1)
        
        temp = signal_data{chan, trial};

        
        % Each trial and channel run contains a data array of ~120 rows of
        % data representing ~120 segments of recording from stim events.
        % Each of these runs are checked to see if any artifacts exist
        for run=1:size(temp,1)
            
            % If a signal is detected crossing the threshold in the defined
            % window...
            if (sum(find(temp(run,offset:end-offset) > thresh | temp(run,offset:end-offset) < -thresh)) > 0)
                
                % Find the index of the maximum value and apply linear interpolated blanking
                peak_idx = find(temp(run,offset:end-offset) > thresh | temp(run,offset:end-offset) < -thresh);
                for i=1:length(peak_idx)
                    
                    %linspace is the best way to do this; draw a line
                    %betwen idx-30 and idx+90 from peak idx
                    temp(run,(offset+peak_idx(i)-fs/1000):(offset+peak_idx(i)+4*fs/1000)) = linspace(temp(run,(offset+peak_idx(i)-fs/1000)),temp(run,(offset+peak_idx(i)+4*fs/1000)),151);
                end
                
                % Update list of trials that needed cleaning as metadata
                if ~ismember(trial, cleaned_trials)
                    cleaned_trials = [cleaned_trials, trial];
                end
            end
        end
        
        % Extra plot for debugging
        if isempty(temp)==0 && plot_fig==true
            x_time = linspace(-2,498,size(temp,2));
            figure(12)
            subplot(1,2,1)
            plot(x_time, signal_data{chan, trial}')
            title('before')
            subplot(1,2,2)
            plot(x_time, temp')
            title('after')
            drawnow
        end
        if iscell(clean_data)==1
            clean_data{chan,trial} = temp;
        else
            clean_data = temp; 
        end
    end
end
fprintf('\nDone!\n')


end


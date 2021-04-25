function [clean_data, missing_packets, cleaned_trials] = replace_missing_data_packets(varargin)
% This function searches for missing data packets in a signal (represented as a vast amlitude of 8000 or greater)
% and replaces those amplitudes with NaN

% INPUTS
% =========================================================
%    signal_data   : (MxN num or cell array) array containing data to search through
%    thresh        : (num) threshold value representing detection of dropped data packets in signal
%
% OUTPUTS
% =========================================================
%    clean_data    : (MxN num or cell array) array containing cleaned data
%    missing_packets  : (1xN cell array) cell array holding list of time indices
%                           where data packets were lost for each trial
%    cleaned_trials: (1xN num) list of trials that needed to be cleaned of
%                           lost data packets
%
%
% Jonathan Shulgach
% Last updated: 4/29/2020

plot_fig=false;
cleaned_trials = [];
packets_lost = [];
t_start = [];
t_end = [];
fs=30e3;
signal_data = varargin{1};
if iscell(signal_data)==1
    clean_data = cell(size(signal_data));
    [N_chan, N_trial] = size(signal_data);
else
    temp_data = signal_data;
    [N_chan, N_data] = size(signal_data);
    N_seg = 1;
    N_trial = 1;
    signal_data = cell(1); % convert to cell
    for i=1:N_chan
        signal_data{i,1} = temp_data(i,:);
    end
    clean_data = [];
end

% Load signal data and paramters, otherwise set default values if no additional inputs
if length(varargin) > 1
    thresh = abs(varargin{2});
else
    thresh = 8000; % microvolts
end

%fprintf("Searching for missing data packets | Trial: ")
for trial=1:N_trial
    
    %fprintf([num2str(trial),'|'])
    
    for chan=1:N_chan
        
        temp_data = signal_data{chan, trial};
        [N_seg, N_data] = size(temp_data);
        t_max = round(N_data/fs);
        
        % In case there are multiple segments to look at, each trial and channel run contains a data array of ~120 rows of
        % data representing ~120 segments of recording from stim events. each of these runs are checked to see if any artifacts exist
        for j=1:N_seg
            for i=2:N_data
                % get time indices for packet loss
                if (temp_data(j,i) > thresh && temp_data(j,i-1) < thresh) || (temp_data(j,i) < -thresh && temp_data(j,i-1) > -thresh)
                    t_start = [t_start, i];
                    packets_lost = [packets_lost,i/fs];
                elseif (temp_data(j,i) < thresh && temp_data(j,i-1) > thresh) || (temp_data(j,i) > -thresh && temp_data(j,i-1) < -thresh)
                    t_end = [t_end, i];
                    packets_lost = [packets_lost,i/fs];
                elseif (temp_data(j,i) < -thresh && temp_data(j,i-1) < -thresh) || (temp_data(j,i) > thresh && temp_data(j,i-1) > thresh)
                    packets_lost = [packets_lost,i/fs];
                end
            end
            
            % Make sure time indices match length
            while length(t_start) > length(t_end)
                t_start = t_start(1:end-1);
            end
            while length(t_start) < length(t_end)
                t_end = t_end(1:end-1);
            end
            
            % Apply linear interpolated blanking
            for i=1:length(t_start)
                temp_data(j,t_start(i)-1:t_end(i)+1) = linspace(temp_data(j,t_start(i)-1),temp_data(j,t_end(i)+1),length(t_start(i)-1:t_end(i)+1));
            end
            
            % Extra plot for debugging
            if isempty(temp_data)==0 && plot_fig==true
                x_time = linspace(-2,t_max,N_data);
                figure(12)
                subplot(1,2,1)
                plot(x_time, signal_data{chan, trial}')
                xlim([-2 t_max]);
                ylim([-thresh 1000]);
                %axis tight
                title('before')
                subplot(1,2,2)
                plot(x_time, temp_data)
                xlim([-2 t_max]);
                ylim([-thresh 1000]);
                title('after')
                %axis tight
                drawnow
            end
            if iscell(clean_data)==1
                clean_data{chan,trial} = temp_data;
            else
                clean_data(chan,:) = temp_data;
            end
            
            
        end
        missing_packets{trial} = packets_lost;
        
    end
    %fprintf("\n")
end
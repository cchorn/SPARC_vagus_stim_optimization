function baseline = getBaselineParams(base_filePath, N_channels)

%STA_avg_enable = true;

% %randomize index times within actual time range
% new_stimIdx = int32(max(stimIdx)*rand(1,length(stimIdx)));

for chan = N_channels
    
    % load data from specified file path, with the entity data type, and ID of electode channel
    % NOTE: num_chan can be single number for single output (1), or array for multiple (1:32)
    analogData = read_continuousData(base_filePath, char(entity), chan);
    
    %----apply bandpass filter----
    % Create 2nd order butterworth filter, high pass at 150Hz, with output parameters
    [b,a] = butter(2,150/(30e3/2),'high'); %7/30/19 - used to be 300 but lowered to catch unknown CAPs just in case
    %tmpChan = filter(b,a,analogData);
    tmpChan = filtfilt(b,a,analogData);
    
    % get mean and baseline
    baseline{expmt,1}.(['Chan_',num2str(chan)]).mean = mean(tmpChan);
    baseline{expmt,1}.(['Chan_',num2str(chan)]).std = std(tmpChan);
    
end

end
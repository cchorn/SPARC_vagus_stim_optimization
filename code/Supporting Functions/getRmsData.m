function [Vrms_list, Vrms_time] = getRmsData(signalData, time_windowSize)

fs = 30e3;
STA_t_min   = -2; %ms
STA_t_max   = 498; %ms

windowSize=round(fs/(1000*time_windowSize)); %window size should be 30 elements if time_window is 1ms
stepSize=windowSize/10; % step size should be around 33us
Vrms_list = [];

% there might still seriously be large artifacts picked up that are below
% the 8000 railing and an actual response. 400-500uV is the usua response, 
% so anything double that should constitute an artifact 
[signalData, ~, ~] = replace_missing_data_packets(signalData,800); 

%calculate Vrms with one window/step size
window = 1:windowSize;
idx = 1;
while max(window) < length(signalData)
    Vrms_list(idx) = rms(signalData(window));
    window = window + stepSize;
    %window = window + (windowSize-stepSize);
    idx = idx + 1;
end
window = window(1):length(signalData);
Vrms_list(idx) = rms(signalData(window));
Vrms_time = linspace(STA_t_min,STA_t_max,length(Vrms_list));


end


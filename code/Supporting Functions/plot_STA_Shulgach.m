function [sig_thresh, false_peak_count] = plot_STA_Shulgach(analogData, stimTimes, neuralChan, STA_t_min, STA_t_max, rms_scale, rms_thresh_scale, std_thresh_scale, rec_filePath, search_offset, plot_data)

%   DESCRIPTION
%   ===================================================================
%   Plots stim triggered average for specified analog channels from
%   a given ns5 file
%
% The following signal prossessing is done for each channel:
% 1) Filter raw data by passing the signal through butterworth filter with 300Hz high pass at 30kHz
%
% 2) Acquire array of filtered signal snippets containing -2ms:400ms time periods
%        indexed by stim event (120 stim events in total, 120 data snippets of )
%
% 3) Calculate the average of the signal across all 120 data snippets, for
%        every data point in time to get average signal
%
% 4) Perform moving window RMS calculation by running a 10ms RMS window
%        over average signal data snippet, using mid-section method (taking
%        RMS of signal 5ms before, 5ms after). This is the data being
%        evaluated with baseline threshold
%
% 5) Remove abmormal RMS spikes by identifying peak RMS values and dampening surrounding data (5ms before, 5ms after) to 0.5
%
% 6) Calculate baseline threshold by taking RMS of RMS Avg data multiplied by a scale value, the standard deviation of RMS
%        Avg data multiplied by a scale value, add both values and create a
%        line spanning the 400ms period to acquire threshold baseline
%        NOTES:
%            a) The scaling parameters have been manually tuned to minimize total FP and FN
%            b) Adding standard deviation to baseline threshold accounts
%                   for broken electrode channels with larger peak-to-peak
%                   RMS amplitudes.
%
% 7) Plot data
%
%   INPUTS
%   ===================================================================
%   filePath   : (string) path to data file
%   stimChan   : (int) stim channel number
%   neuralChan : (int) eng data channel
%
%   OUTPUTS
%   ===================================================================
%   Figure per STA channel and snippet matrix from [-1 20] ms around
%   stim event
%
%   EXAMPLE
%   ===================================================================
%   plot_STA('R:\data_raw\ferret\2018\01162018 - aF43-17\Trellis\datafile0006.nev', 1,129)
%
%   ACN created 11/17
%   JS modified 6/13/19

if size(stimTimes,2) ~= 1
    error('multiple stim chans detected...')
end

false_peak_count = 0;
rms_window       = 5;     %ms
step             = 1;     %33us
damp_window      = 360;   %ms
fs               = 30e3;
rms_method       = 'mid';
numChans         = length(neuralChan);
stimIdx          = ceil(stimTimes{1}*fs);
std_scale_low    = 4.5;
std_scale_high   = 20;
fig1             = true;
fig2             = true;
fig3             = false;
fig4             = true;
fig5             = true;
fig6             = true;
fig7             = false;


chanSnips = cell(1,numChans);
Base_chanSnips = cell(1,numChans);

for i = 1:numChans
    while stimIdx(1) < (abs(STA_t_max) + abs(STA_t_min) + 1) % check that first stim is included within time window
        stimIdx = stimIdx(2:end);
    end
        %----apply bandpass filter----
    %old method
   [b,a] = butter(2,150/(30e3/2),'high'); %7/30/19 - used to be 300 but lowered to catch unknown CAPs just in case
   tmpChan = filter(b,a,analogData(i,:));
    
    % Implement 60Hz notch filter
    d = designfilt('bandstopiir','FilterOrder',2, ...
        'HalfPowerFrequency1',30,'HalfPowerFrequency2',90, ...
        'DesignMethod','butter','SampleRate',fs);
    
    tmpChan_notch = filtfilt(d,tmpChan);
    
    
    baseline_mean = mean(tmpChan_notch); % baseline mean
    baseline_std = std(tmpChan_notch); % baseline std
    
    %raw_chanSnips{i} = cell2mat(arrayfun(@(x) analogData(i,(x-abs(STA_t_min)):(x+abs(STA_t_max))), stimIdx,'UniformOutput',false)');
    %filt_chanSnips{i} = cell2mat(arrayfun(@(x) tmpChan((x-abs(STA_t_min)):(x+abs(STA_t_max))), stimIdx,'UniformOutput',false)');
    notch_chanSnips{i} = cell2mat(arrayfun(@(x) tmpChan_notch((x-abs(STA_t_min)):(x+abs(STA_t_max))), stimIdx,'UniformOutput',false)');
    
    %remove crazy signal spikes
    snipData = notch_chanSnips{i};
    snipData(snipData(:,:)>=60)=baseline_mean;
    snipData(snipData(:,:)<=-60)=baseline_mean;
    avg_chanData = mean(snipData);
    
    x_time = linspace(STA_t_min,STA_t_max,size(snipData,2))*1000/fs;
      
    

    
    %Moving window RMS calculation for notch filter
    if strcmpi(rms_method,'lag')
        temp = [zeros(1,rms_window*fs/1000), avg_chanData];
        for k=1:step:length(avg_chanData)
            rms_avg(k) = rms(temp(k:k+rms_window));
        end
    elseif strcmpi(rms_method,'mid')
        temp = [zeros(1,(rms_window*fs/1000)/2-1), avg_chanData, zeros(1,(rms_window*fs/1000)/2-1)];
        for k=1:step:length(avg_chanData)
            rms_avg(k) = rms(temp(k:k+(rms_window*fs/1000)/2-1));
        end
    elseif strcmpi(rms_method,'lead')
        for k=1:step:(length(avg_chanData)-rms_window*fs/1000)
            rms_avg(k) = rms((avg_chanData:(avg_chanData+rms_window*fs/1000)));
        end
    end
    %rms_avg(rms_avg==0) = []; % remove zeros
    rms_x_time = linspace(STA_t_min,STA_t_max,size(rms_avg,2))*1000/fs;
    
    %Moving window mean calculation
%         if strcmpi(rms_method,'lag')
%             temp = [zeros(1,rms_window*fs/1000), avg_chanData];
%             for k=1:step:length(avg_chanData)
%                 s_avg(k) = mean(temp(k:k+rms_window));
%             end
%         elseif strcmpi(rms_method,'mid')
%             temp = [zeros(1,(rms_window*fs/1000)/2-1), avg_chanData, zeros(1,(rms_window*fs/1000)/2-1)];
%             for k=1:step:length(avg_chanData)
%                 s_avg(k) = mean(temp(k:k+(rms_window*fs/1000)/2-1));
%             end
%         elseif strcmpi(rms_method,'lead')
%             for k=1:step:(length(avg_chanData)-rms_window*fs/1000)
%                 s_avg(k) = mean((avg_chanData:(avg_chanData+rms_window*fs/1000)));
%             end
%         end

    s_avg = movmean(avg_chanData, (rms_window*fs/1000));


    s_x_time = linspace(STA_t_min,STA_t_max,size(s_avg,2))*1000/fs;

    %     %remove crazy mean spikes
    %     [peak_val, peak_idx] = findpeaks(s_avg);
    %     damp_window = round(damp_window/2);
    %     for j=1:length(peak_val)
    %         if peak_val(j) > 4
    %             if ((length(s_avg) - peak_idx(j)) <= damp_window)
    %                 s_avg(peak_idx(j):end) = 0.1;
    %             elseif peak_idx(j) <= damp_window
    %                 s_avg(1:peak_idx(j)) = 0.1;
    %             else
    %                 s_avg((peak_idx(j)-damp_window):(peak_idx(j)+damp_window)) = 0.1;
    %             end
    %             %false_peak_count = false_peak_count + 1;
    %         end
    %     end
    
    
    sig_thresh.raw_avg = avg_chanData;
    sig_thresh.raw_mean = mean(avg_chanData(STA_t_min+search_offset:end-search_offset));
    sig_thresh.raw_std = std(avg_chanData(STA_t_min+search_offset:end-search_offset));
    sig_thresh.raw_thresh_low = sig_thresh.raw_mean + std_scale_low*sig_thresh.raw_std; % low_thresh = mean + A*std_offset
    sig_thresh.raw_thresh_high = sig_thresh.raw_mean + std_scale_high*sig_thresh.raw_std; % high_thresh = mean + B*std_offset
    sig_thresh.raw_search_window = avg_chanData(STA_t_min+search_offset:end-search_offset);
    [raw_freq, raw_val] = hist(sig_thresh.raw_search_window,100);
    sig_thresh.raw_freq = raw_freq;
    sig_thresh.raw_val = raw_val;
    
    sig_thresh.s_avg = s_avg;
    sig_thresh.s_mean = mean(s_avg(STA_t_min+search_offset:end-search_offset));
    sig_thresh.s_std = std(s_avg(STA_t_min+search_offset:end-search_offset));
    sig_thresh.s_thresh_low = sig_thresh.s_mean + std_scale_low*sig_thresh.s_std; % low_thresh = mean + A*std_offset
    sig_thresh.s_thresh_high = sig_thresh.s_mean + std_scale_high*sig_thresh.s_std; % high_thresh = mean + B*std_offset
    sig_thresh.s_search_window = s_avg(STA_t_min+search_offset:end-search_offset);
    [s_freq, s_val] = hist(sig_thresh.s_search_window,100);
    sig_thresh.s_freq = s_freq;
    sig_thresh.s_val = s_val;
    
    sig_thresh.rms_avg = rms_avg;
    sig_thresh.rms_mean = mean(rms_avg(STA_t_min+search_offset:end-search_offset));
    sig_thresh.rms_std = std(rms_avg(STA_t_min+search_offset:end-search_offset));
    sig_thresh.rms_thresh_low = sig_thresh.rms_mean + std_scale_low*sig_thresh.rms_std; % low_thresh = mean + A*std_offset
    sig_thresh.rms_thresh_high = sig_thresh.rms_mean + std_scale_high*sig_thresh.rms_std; % high_thresh = mean + B*std_offset
    sig_thresh.rms_search_window = rms_avg(STA_t_min+search_offset:end-search_offset);
    [rms_freq, rms_val] = hist(sig_thresh.rms_search_window,100);
    sig_thresh.rms_freq = rms_freq;
    sig_thresh.rms_val = rms_val;
    
    
    if plot_data==1
        % Raw Signal figure
        if fig1==true
            figure(1)
            plot(x_time,snipData') % Plot full data
            hold on;
            plot(x_time, avg_chanData,'-k'); %data average
            plot(s_x_time, s_avg, '-b'); % data with window mean
            hold off;
            axis tight
            if max(max(sig_thresh.rms_search_window))>70
                ylim([-70 (max(sig_thresh.rms_search_window)+20)]);
            else
                ylim([-70 70]);
            end
            drawnow
            set(gca,'ButtonDownFcn',@createnew_fig)
        end
        %RMS, mean, and avg line plot
        if fig2==true
            figure(2)
            %yyaxis right
            %plot(x_time, avg_chanData,'-k','LineWidth',0.8); %data average
            
            plot(s_x_time, s_avg, '-b', 'LineWidth', 0.8); % data with window mean
            hold on;
            plot(rms_x_time, rms_avg, '-g', 'LineWidth', 0.8); % data with rms window
            %hold on;
            %plot(x_time, s_avg, '-b');
            %thresh_low = linspace(sig_thresh.rms_thresh_low, sig_thresh.rms_thresh_low, size(snipData,2));
            %thresh_high = linspace(sig_thresh.rms_thresh_high, sig_thresh.rms_thresh_high, size(snipData,2));
            %thresh_low = linspace(sig_thresh.raw_thresh_low, sig_thresh.raw_thresh_low, size(snipData,2));
            %thresh_high = linspace(sig_thresh.raw_thresh_high, sig_thresh.raw_thresh_high, size(snipData,2));
            thresh_low = linspace(sig_thresh.s_thresh_low, sig_thresh.s_thresh_low, size(s_avg,2));
            thresh_high = linspace(sig_thresh.s_thresh_high, sig_thresh.s_thresh_high, size(s_avg,2));
            plot(s_x_time, thresh_low,'-r');
            plot(s_x_time, -thresh_low,'-r');
            %plot(x_time, thresh_low,'-r');
            %plot(x_time, -thresh_low,'-r');

            hold off;
            axis tight
            if max(sig_thresh.s_search_window)>4
                ylim([-1 (max(sig_thresh.s_search_window)+2)]);
            else
                ylim([-1 4]);
            end
            set(gca,'ButtonDownFcn',@createnew_fig)
        end
        
        % Search Window Plot
        hold on;
        color = [0, 0, 1];
        a=fill([search_offset*1000/fs, search_offset*1000/fs, (length(x_time)-search_offset)*1000/fs,...
            (length(x_time)-search_offset)*1000/fs], [-50, 50, 50, -50], color);
        a.FaceAlpha = 0.1;
        xlabel('time (msec)');
        ylabel('\muV');
        hold off;
        
        
        %rms_above_thresh = sum(rms_threshold.rms_search_window(rms_threshold.rms_search_window > rms_threshold.rms_thresh_low));
        %rms_above_thresh_count = length(rms_threshold.rms_search_window(rms_threshold.rms_search_window > rms_threshold.rms_thresh_low)); %SAME
        
        %Histogram Plot
        if fig3==true
            figure(3)
            hist(sig_thresh.raw_search_window,100)
            %hist(rms_threshold.rms_search_window,100)
            hold on;
            %plot(rms_threshold.rms_thresh_low, 1:250, 'r.'); %RMS line plot
            %plot(rms_threshold.rms_thresh_high, 1:250, 'r.');
            plot(sig_thresh.raw_thresh_low, 1:250, 'r.'); %RMS line plot
            plot(-sig_thresh.raw_thresh_low, 1:250, 'r.'); %RMS line plot
            %plot(0.1:0.001:rms_threshold.rms_thresh_high,60,'r.');
            %xlabel('Signal RMS');
            xlabel('Signal');
            ylabel('Frequency');
            hold off;
            drawnow
        end
        
        % average signal
        if fig4==true
            figure(4)
            plot(x_time, avg_chanData,'-k','LineWidth',0.8); %data average
            hold on
            thresh_low = linspace(sig_thresh.raw_thresh_low, sig_thresh.raw_thresh_low, size(avg_chanData,2));
            thresh_high = linspace(sig_thresh.raw_thresh_high, sig_thresh.raw_thresh_high, size(avg_chanData,2));
            plot(x_time, thresh_low,'-r');
            plot(x_time, -thresh_low,'-r');
            hold off
            axis tight
            ylim([-4 4]);
            xlabel('time (msec)');
            ylabel('\muV');
            % Search Window Plot
            hold on;
            color = [0, 0, 1];
            a=fill([search_offset*1000/fs, search_offset*1000/fs, (length(x_time)-search_offset)*1000/fs,...
                (length(x_time)-search_offset)*1000/fs], [-50, 50, 50, -50], color);
            a.FaceAlpha = 0.1;
            xlabel('time (msec)');
            ylabel('\muV');
            hold off;
            drawnow
        end
        
        % signal mean window
        if fig5==true
            figure(5)
            plot(s_x_time, s_avg, '-b', 'LineWidth', 0.8); % data with window mean
            hold on
            thresh_low = linspace(sig_thresh.s_thresh_low, sig_thresh.s_thresh_low, size(s_avg,2));
            thresh_high = linspace(sig_thresh.s_thresh_high, sig_thresh.s_thresh_high, size(s_avg,2));
            plot(s_x_time, thresh_low,'-r');
            plot(s_x_time, -thresh_low,'-r');
            hold off
            axis tight
            ylim([-4 4]);
            xlabel('time (msec)');
            ylabel('\muV');
            % Search Window Plot
            hold on;
            color = [0, 0, 1];
            a=fill([search_offset*1000/fs, search_offset*1000/fs, (length(x_time)-search_offset)*1000/fs,...
                (length(x_time)-search_offset)*1000/fs], [-50, 50, 50, -50], color);
            a.FaceAlpha = 0.1;
            xlabel('time (msec)');
            ylabel('\muV');
            hold off;
            drawnow
        end
        
        % rms signal window
        if fig6==true
            figure(6)
            plot(rms_x_time, rms_avg, '-g', 'LineWidth', 0.8); % data with rms window
            hold on
            thresh_low = linspace(sig_thresh.rms_thresh_low, sig_thresh.rms_thresh_low, size(rms_avg,2));
            thresh_high = linspace(sig_thresh.rms_thresh_high, sig_thresh.rms_thresh_high, size(rms_avg,2));
            plot(rms_x_time, thresh_low,'-r');
            plot(rms_x_time, -thresh_low,'-r');
            hold off
            axis tight
            ylim([-4 4]);
            xlabel('time (msec)');
            ylabel('\muV');
            % Search Window Plot
            hold on;
            color = [0, 0, 1];
            a=fill([search_offset*1000/fs, search_offset*1000/fs, (length(x_time)-search_offset)*1000/fs,...
                (length(x_time)-search_offset)*1000/fs], [-50, 50, 50, -50], color);
            a.FaceAlpha = 0.1;
            xlabel('time (msec)');
            ylabel('\muV');
            hold off;
            drawnow
        end
        
        % Overlay all together

        if fig7==true
            figure(7)
            plot(x_time, avg_chanData,'-k','LineWidth',0.8); %data average
            hold on
            plot(s_x_time, s_avg, '-b', 'LineWidth', 0.8); % data with window mean
            plot(rms_x_time, rms_avg, '-g', 'LineWidth', 0.8); % data with rms window
            hold off
            axis tight
            ylim([-4 4]);
            xlabel('time (msec)');
            ylabel('\muV');
            % Search Window Plot
            hold on;
            color = [0, 0, 1];
            a=fill([search_offset*1000/fs, search_offset*1000/fs, (length(x_time)-search_offset)*1000/fs,...
                (length(x_time)-search_offset)*1000/fs], [-50, 50, 50, -50], color);
            a.FaceAlpha = 0.1;
            xlabel('time (msec)');
            ylabel('\muV');
            hold off;
            drawnow
        end
        
        if fig1==true
            figure(1)
        end
        %rms_gauss_above_thresh = sum(rms_threshold.freq(rms_threshold.val > rms_threshold.rms_thresh_low));
        %rms_gauss_above_thresh_count = sum(rms_threshold.freq(rms_threshold.val > rms_threshold.rms_thresh_low)); %SAME
    end
end
end
% script to load and bin vitals data and save dataframe 
%% EKG/HR
allSubs = {'15-18', '13-18','14-18'};
trials = {[12,13],[12,13],[7,12]};
allHR = [];
allLabels = [];
subLabel = [];
trialLabel = [];
for i = 1:3
    subject = allSubs{i};
    for iTrial = trials{i}
        trialObj= mdf.load('subject',subject,'trial',iTrial,'mdf_type','trial');
        stimStart = trialObj.params.prestim_min;
        stimEnd = trialObj.params.prestim_min+trialObj.params.stim_min;
        
        vitalsChan = mdf.load('subject',subject,'trial',trialObj.trial,'mdf_type','vitals', 'chanLabel','A5');
        vitalWF = vitalsChan.wf;
        thresh = 3000; 
        [peaks, pkIdx] = findpeaks(vitalWF,'MinPeakHeight',thresh, 'MinPeakDistance',5000); % std Hr = 300 per minute
        pkTimes = pkIdx/30e3/60;

        binWidth = 0.25;            % binwidth*60 seconds
        bins = 0:binWidth:(length(vitalWF)/30e3/60);
        HR = histc(pkTimes,bins);
        
        G = ones(size(HR));
        G(bins<stimStart) = 0;
        allHR = [allHR; HR(bins<=stimEnd)*4];
        allLabels = [allLabels; G(bins<=stimEnd)];
        subLabel = [subLabel;i*ones(length(HR(bins<=stimEnd)),1)];
        trialLabel = [trialLabel; iTrial*ones(length(HR(bins<=stimEnd)),1)];

    end
end

%
subject = '14-18';
for iTrial = 5 
    trialObj= mdf.load('subject',subject,'trial',iTrial,'mdf_type','trial');
    vitalsChan = mdf.load('subject',subject,'trial',trialObj.trial,'mdf_type','vitals', 'chanLabel','A5');
    vitalWF = [vitalsChan(1).wf;vitalsChan(2).wf;vitalsChan(3).wf;vitalsChan(4).wf];
    
    thresh = 3000;
    [peaks, pkIdx] = findpeaks(vitalWF,'MinPeakHeight',thresh, 'MinPeakDistance',5000); % std Hr = 300 per minute
    pkTimes = pkIdx/30e3/60;

    binWidth = 0.25;            % binwidth*60 seconds
    bins = 0:binWidth:(length(vitalWF)/30e3/60);
    HR = histc(pkTimes,bins);

    allHR = [allHR; HR*4];
    allLabels = [allLabels; zeros(length(HR),1)];
    subLabel = [subLabel;3*ones(length(HR),1)];
    trialLabel = [trialLabel; 5*ones(length(HR),1)];

end

HRtable = array2table([allHR, allLabels, trialLabel,subLabel],'VariableNames',{'HR','label','trial','Subject'});
writetable(HRtable,'HRdata.csv')

%% Resp
allResp = [];
allLabelsResp = [];
subLabelResp = [];
trialLabelResp = [];
for i = 1:3
    subject = allSubs{i};
    for iTrial = trials{i} %[12, 13]
        trialObj= mdf.load('subject',subject,'trial',iTrial,'mdf_type','trial');
        stimStart = trialObj.params.prestim_min;
        stimEnd = trialObj.params.prestim_min+trialObj.params.stim_min;
        
        vitalsChan = mdf.load('subject',subject,'trial',trialObj.trial,'mdf_type','vitals', 'chanLabel','A7');
        vitalWF = -vitalsChan.wf;

        thresh = -80;
        [peaks, pkIdx] = findpeaks(vitalWF,'MinPeakHeight',thresh, 'MinPeakDistance',30e3); % std breath = 30 per minute
        pkTimes = pkIdx/30e3/60;

        binWidth = 0.25;            % binwidth*60 seconds
        bins = 0:binWidth:(length(vitalWF)/30e3/60);
        Resp = histc(pkTimes,bins);
        
        G = ones(size(Resp));
        G(bins<stimStart) = 0;
        allResp = [allResp; Resp(bins<=stimEnd)*4];
        allLabelsResp = [allLabelsResp; G(bins<=stimEnd)];
        subLabelResp = [subLabelResp;i*ones(length(Resp(bins<=stimEnd)),1)];
        trialLabelResp = [trialLabelResp; iTrial*ones(length(Resp(bins<=stimEnd)),1)];

    end
end

subject = '14-18';
for iTrial = 5
    trialObj= mdf.load('subject',subject,'trial',iTrial,'mdf_type','trial');
    
    vitalsChan = mdf.load('subject',subject,'trial',trialObj.trial,'mdf_type','vitals', 'chanLabel','A7');
    vitalWF = -[vitalsChan(1).wf;vitalsChan(2).wf;vitalsChan(3).wf;vitalsChan(4).wf];

    thresh = -80;
    [peaks, pkIdx] = findpeaks(vitalWF,'MinPeakHeight',thresh, 'MinPeakDistance',30e3); % std breath = 30 per minute
    pkTimes = pkIdx/30e3/60;
    plot(pkTimes,peaks,'.')
    ylim([-120, -50])

    binWidth = 0.25;            % binwidth*60 seconds
    bins = 0:binWidth:(length(vitalWF)/30e3/60);
    Resp = histc(pkTimes,bins);

    allResp = [allResp; Resp*4];
    allLabelsResp = [allLabelsResp; zeros(size(Resp))];
    subLabelResp = [subLabelResp; i*ones(size(Resp))];
    trialLabelResp = [trialLabelResp; 5*ones(size(Resp))];

end


Resptable = array2table([allResp, allLabelsResp, trialLabelResp, subLabelResp],'VariableNames',{'Resp','label','trial','Subject'});
writetable(Resptable,'Respdata.csv')

%% Temp
allSubs = {'15-18', '13-18','14-18'};
trials = {[12,13],[12,13],[7,12]};
allTemp = [];
allLabels = [];
subLabel = [];
trialLabel = [];
for i = 1:3
    subject = allSubs{i};
    x = 1;
    for iTrial = trials{i}
        trialObj= mdf.load('subject',subject,'trial',iTrial,'mdf_type','trial');
        stimStart = trialObj.params.prestim_min;
        stimEnd = trialObj.params.prestim_min+trialObj.params.stim_min;
        
        vitalsChan = mdf.load('subject',subject,'trial',trialObj.trial,'mdf_type','vitals', 'chanLabel','A6');
        vitalWF = vitalsChan.wf;
        vitalWF(vitalWF<-4000)=mean(vitalWF(1:10));
        vitalWF = (vitalWF-min(vitalWF))*(4)/(max(vitalWF)-min(vitalWF)) + 100+rand; %(vitalWF-min(vitalWF))/103;
        baseTemp = mean(vitalWF(time<=stimStart));
        stimTemp = mean(vitalWF((time>=stimStart) & (time<=stimEnd)));
        
      allTemp = [allTemp; baseTemp; stimTemp];
      allLabels = [allLabels; 0; 1];
      trialLabel = [trialLabel; x; x];
      subLabel = [subLabel; i; i];
      x=x+1;
    end
end
HRtable = array2table([allTemp, allLabels, trialLabel,subLabel],'VariableNames',{'Temp','label','trial','Subject'});
writetable(HRtable,'TemperatureData.csv')

%% BP

allSubs = {'15-18', '13-18','14-18'};
trials = {[12,13],[12,13],[7,12]};
allHR = [];
allLabels = [];
subLabel = [];
trialLabel = [];
for i = 1:3
    subject = allSubs{i};
        x= 1;
    
    for iTrial = trials{i}
        trialObj= mdf.load('subject',subject,'trial',iTrial,'mdf_type','trial');
        stimStart = trialObj.params.prestim_min;
        stimEnd = trialObj.params.prestim_min+trialObj.params.stim_min;
        
%         h = figure;
        vitalsChan = mdf.load('subject',subject,'trial',trialObj.trial,'mdf_type','vitals', 'chanLabel','A8');
        vitalWF = vitalsChan.wf;
        vitalWF(vitalWF<-4000)=mean(vitalWF(1:10));
        [b,a]=butter(2,[1,50]/(30e3/2),'bandpass');
        vitalWF = filtfilt(b,a,vitalWF);
        vitalWF = (vitalWF-min(vitalWF))*(54)/(max(vitalWF)-min(vitalWF)) + 110+rand*10;

        thresh = 150; %105; %input('set threshold ');
        [peaks, pkIdx] = findpeaks(vitalWF, 'MinPeakHeight',thresh,'MinPeakDistance',5000); % std Hr = 300 per minute
        pkTimes = pkIdx/30e3/60;

        allHR = [allHR; mean(peaks(pkTimes<=stimStart)); mean(peaks((pkTimes>=stimStart)&(pkTimes<=stimEnd)))];
        allLabels = [allLabels; 0;1];
        subLabel = [subLabel;i;i];
        trialLabel = [trialLabel; x;x];
        x=x+1;

    end
end


subject = '14-18';
for iTrial = 5 
    trialObj= mdf.load('subject',subject,'trial',iTrial,'mdf_type','trial');
    vitalsChan = mdf.load('subject',subject,'trial',trialObj.trial,'mdf_type','vitals', 'chanLabel','A8');
    vitalWF = [vitalsChan(1).wf;vitalsChan(2).wf;vitalsChan(3).wf;vitalsChan(4).wf];
        vitalWF(vitalWF<-4000)=mean(vitalWF(1:10));
        [b,a]=butter(2,[1,50]/(30e3/2),'bandpass');
        vitalWF = filtfilt(b,a,vitalWF);
        vitalWF = (vitalWF-min(vitalWF))*(54)/(max(vitalWF)-min(vitalWF)) + 110+rand*10;
        thresh =  150;
    [peaks, pkIdx] = findpeaks(vitalWF, 'MinPeakHeight', thresh,'MinPeakDistance',5000); % std Hr = 300 per minute
    pkTimes = pkIdx/30e3/60;
    
    allHR = [allHR; mean(peaks(pkTimes<=stimStart)); mean(peaks((pkTimes>=stimStart)&(pkTimes<=stimEnd)))];
    allLabels = [allLabels; 0;1];
    subLabel = [subLabel;i;i];
    trialLabel = [trialLabel; 2;2];

end

BPdata = array2table([allHR, allLabels, trialLabel,subLabel],'VariableNames',{'BP','label','trial','Subject'});
writetable(BPdata,'BPdata.csv')
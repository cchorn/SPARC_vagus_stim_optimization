% script to load GMA data and generate waterfall plots

loadData(subject, trial);
relevantWindow = [9,10];
freqWindow = [8,11];
baselineWindowMask = [res_baseline{3}.winFeat.winStart>=relevantWindow(1) & res_baseline{3}.winFeat.winStart<=relevantWindow(2)];
lastIdx = find(baselineWindowMask,1,'last');
figure; 
subplot(2,4,[3,4,7,8])
hold on
plotWaterfall([res_baseline{3}.winFeat.wf_spectroPower(:,baselineWindowMask), res12{3}.winFeat.wf_spectroPower(:,2:8)], ...
              [res_baseline{3}.winFeat.winStart(baselineWindowMask), res_baseline{3}.winFeat.winStart(lastIdx)+0.15 + res12{3}.winFeat.winStart(2:8)],...
              (6:0.1:15)/60)
fill3([6,15,15,6],[9,9,10,10],[0,0,0,0],[0,0,1], 'FaceAlpha',0.5)
fill3([6,15,15,6],[10,10,11.2,11.2],[0,0,0,0],[1,0,0], 'FaceAlpha',0.5)

normoBinIdx = [res_baseline{3}.winFeat.freqBins*60>freqWindow(1) & res_baseline{3}.winFeat.freqBins*60<freqWindow(2)];
normoPower_base = sum(sum(res_baseline{3}.winFeat.wf_spectroPower(normoBinIdx, baselineWindowMask)));
normoPower_stim = sum(sum(res12{3}.winFeat.wf_spectroPower(normoBinIdx, 2:8)));
subplot(2,4,[5,6])
bar([normoPower_base/(normoPower_base+normoPower_stim), normoPower_stim/(normoPower_base+normoPower_stim)])

%
relevantWindow = [4.5,6];
baselineWindowMask = [res_baseline{3}.winFeat.winStart>=relevantWindow(1) & res_baseline{3}.winFeat.winStart<=relevantWindow(2)];
lastIdx = find(baselineWindowMask,1,'last');
figure; 
subplot(2,4,[3,4,7,8])
hold on
plotWaterfall([res_baseline{3}.winFeat.wf_spectroPower(:,baselineWindowMask), res7{3}.winFeat.wf_spectroPower(:,2:8)*2], ...
              [res_baseline{3}.winFeat.winStart(baselineWindowMask), res_baseline{3}.winFeat.winStart(lastIdx)+0.15 + res7{3}.winFeat.winStart(2:8)],...
              (6:0.1:15)/60)
fill3([6,15,15,6],[4.5,4.5,6,6],[0,0,0,0],[0,0,1], 'FaceAlpha',0.5)
fill3([6,15,15,6],[6,6,7,7],[0,0,0,0],[1,0,0], 'FaceAlpha',0.5)

normoPower_base = sum(sum(res_baseline{3}.winFeat.wf_spectroPower(normoBinIdx, baselineWindowMask)))/1.5;
normoPower_stim = sum(sum(res7{3}.winFeat.wf_spectroPower(normoBinIdx, 1:8)));
subplot(2,4,[5,6])
bar([normoPower_base/(normoPower_base+normoPower_stim), normoPower_stim/(normoPower_base+normoPower_stim)])
ylim([0,1])


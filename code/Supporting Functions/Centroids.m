% Calculate the centroid of each set of active channels

%   Written by: Dylan Beam on 11/16/2020
%   Last edited by: Dylan Beam on 12/3/2020

%% Set up
%clc,clear

expmt        = 2; % choose animal (2 to 7)

PW_idx_1     = 1; % PW session for cuff pair 1:2 (choose from 1 to 3)
stim_idx_1   = 4; % stimulation amplitudes (select from PW session)

PW_idx_2     = 4; % PW session for cuff pair 3:4 (choose from 4 to 6)
stim_idx_2   = 4;


%% Order results by amplitude
% 
% %extract the stimulation history to temper with
% stim_hist = subject.stim_hist;
% 
% %loop through both cuff/pw pairs
% for trial = 1:6
%     
%     %extract the stimulation history 
%         
%     %loop through all amplitudes
%     for amp_trial = 1:10
%             
%         %identify the minimum amplitude for the run
%         [temp_min, index] = min(stim_hist(trial,:));
%         
%         %save the data from the trial in a new matrix
%         for channel = 1:32
%             ordered_stim(amp_trial, channel, trial) = subject.selectivityGrid(trial, index, channel);
%         end
%         
%         %rename the minimum in the stim_hist so a new minimum can be
%         %detected in the next pass
%         stim_hist(trial,index) = 5000;
%         
%     end      
%     
% end

%load expmt_list
subject = expmt_list{expmt,1};
ordered_stim = zeros(10,32,6);
%% Calculate the centroid for one cuff pair

%define what the layout of the channel looks like
geographic_grid = [ 1,5,9,13,17,21,25,29; 
                    2,6,10,14,18,22,26,30;
                    3,7,11,15,19,23,27,31;
                    4,8,12,16,20,24,28,32];
                
%create a reference variable of which channels are activated
%active_channels_1 = ordered_stim(stim_trial_1,:,amp_trial_1);
active_channels_1 = reshape(subject.selectivityGrid(PW_idx_1,stim_idx_1,:),1,32);
%active_channels_1 = subject.selectivityGrid(amp_trial_1, stim_trial_1, :);
row1 = [];
column1 = [];

%loop through all channels
for i = 1:32
    
    %skip pass if channel was inactive
    if active_channels_1(i) == 0
        continue
    end
    
    %find the locations for each active channel
    [row1(length(row1) + 1), column1(length(column1) + 1)] = gridLoc(geographic_grid,i);    
    
end

%average the row and column locations
rowAvg1 = mean(row1);
columnAvg1 = mean(column1);

%% Calculate the centroid for the second cuff pair

%create a reference variable of which channels are activated
%active_channels_2 = ordered_stim(stim_trial_2,:,amp_trial_2);
active_channels_2 = subject.selectivityGrid(PW_idx_2, stim_idx_2, :);

row2 = [];
column2 = [];

%loop through all channels
for i = 1:32
    
    %skip pass if channel was inactive
    if active_channels_2(i) == 0
        continue
    end
    
    %find the locations for each active channel
    [row2(length(row2) + 1), column2(length(column2) + 1)] = gridLoc(geographic_grid,i);    
    
end

%average the row and column locations
rowAvg2 = mean(row2);
columnAvg2 = mean(column2);

%fprintf("cuff 1:2 average position= (%0.2f, %0.2f)\n",rowAvg1,columnAvg1);
%fprintf("cuff 3:4 average position= (%0.2f, %0.2f)\n",rowAvg2,columnAvg2);
%% Extra Calculations

%convert to actual length (mm)
rowAvg1 = rowAvg1 * 0.4;
rowAvg2 = rowAvg2 * 0.4;
columnAvg1 = columnAvg1 * 0.4;
columnAvg2 = columnAvg2 * 0.4;
fprintf("cuff 1:2 average position = (%0.2f, %0.2f)mm\n",rowAvg1,columnAvg1);
fprintf("cuff 3:4 average position = (%0.2f, %0.2f)mm\n",rowAvg2,columnAvg2);

%calculate the distance between the centroids
dist = sqrt( (columnAvg1 - columnAvg2)^2 + (rowAvg1 - rowAvg2)^2);
fprintf("Distance between centroids: %0.2fmm\n",dist);
%% Plot results

figure
scatter([columnAvg1, columnAvg2], [rowAvg1, rowAvg2], 'filled');
title('Centroids (F21-19, PW = 1000us, 1500uA)')
ylabel('Electrode Location (y)')
xlabel('Electrode Location (x)')
% Limits to fit MEA pins pitch spacing in millimeters
xlim([0.4,3.2])
ylim([0.4,1.6])
% Limits for MEA channel spacing
%xlim([1,8])
%ylim([1,4])
% Evaluates the distance between active channels and the next nearest
% active channel

%   Written by: Dylan Beam on 10/20/2020
%   Last edited by: Dylan Beam on 12/4/2020

%% Set up
%clc,clear

%load expmt_list
subject = expmt_list{7,1};
PW_idx_1 = 2;
PW_idx_2 = 5;
stim_idx_1 = 3;
stim_idx_2 = 5;

%define what the layout of the channel looks like
geographic_grid = [ 1,5,9,13,17,21,25,29;
    2,6,10,14,18,22,26,30;
    3,7,11,15,19,23,27,31;
    4,8,12,16,20,24,28,32];

ordered_stim = zeros(size(subject.stim_hist,2),32,6);
diff_elec1(1:32,1:32) = 10;
diff_elec2(1:32,1:32) = 10;
same_elec1(1:32,1:32) = 10;
same_elec2(1:32,1:32) = 10;

fig_title = ['Distance Between Active Channels (F',subject.cohort,', PW = ',num2str(subject.pulseWidth(PW_idx_1)*1000),'us, ',...
    num2str(subject.stim_hist(PW_idx_1,stim_idx_1)),'uA and ', num2str(subject.stim_hist(PW_idx_2,stim_idx_2)),'uA)'];
%% Organize Data

%order results by stim amplitude for indexing

%extract the stimulation history to temper with
stim_hist = subject.stim_hist;

%loop through both cuff/pw pairs
for session = 1:6
    
    %     %extract the stimulation history
    %
    %     %loop through all amplitudes
    %     for amp_trial = 1:10
    %
    %         %identify the minimum amplitude for the run
    %         [temp_min, index] = min(stim_hist(session,:));
    %
    %         %save the data from the trial in a new matrix
    %         for channel = 1:32
    %             ordered_stim(amp_trial, channel, session) = subject.selectivityGrid(session, index, channel);
    %         end
    %
    %         %rename the minimum in the stim_hist so a new minimum can be
    %         %detected in the next pass
    %         stim_hist(session,index) = 5000;
    %
    %     end
    
    
    % Easiest way to know which index of PW session and stimulation trial is to
    % use the "max SI params" table.
    ordered_stim(:,:,session) = subject.selectivityGrid(session,:,:);
    
end

%eliminate overlapping channels

%loop through all channels
% for i = 1:32
%     
%     %if both cuffs have a response or no response, make equal to 0
%     if ordered_stim(stim_idx_1,i,PW_idx_1) == ordered_stim(stim_idx_2,i,PW_idx_2)
%         ordered_stim(stim_idx_1,i,PW_idx_1) = 0;
%         ordered_stim(stim_idx_2,i,PW_idx_2) = 0;
%     end
%     
% end

%% Identify the distances between one cuff pair and the next cuff pair at thresh and max stim

%loop through all channels to calculate distances
for i = 1:32
    
    if ordered_stim(stim_idx_1,i,PW_idx_1) == 0
        continue
    end
    
    %loop through all of the other channels
    for j = 1:32
        
        if ordered_stim(stim_idx_2,j,PW_idx_2) == 0
            continue
        elseif i == j
            continue
        end
        
        %identify locations for each electrode
        [cuffloc1x, cuffloc1y] = gridLoc(geographic_grid, j);
        [cuffloc2x, cuffloc2y] = gridLoc(geographic_grid, i);
        
        %calculate distances between electrodes
        diff_elec1(j,i) = sqrt( (cuffloc1x - cuffloc2x)^2 + (cuffloc1y - cuffloc2y)^2);
        diff_elec2(i,j) = sqrt( (cuffloc1x - cuffloc2x)^2 + (cuffloc1y - cuffloc2y)^2);
    end
    
end

%% Identify the distances between electrodes on one cuff pair, and other electrodes on the same cuff pair

%loop through all channels to calculate distances
for i = 1:32
    
    if ordered_stim(stim_idx_1,i,PW_idx_1) == 0
        continue
    end
    
    %loop through all of the other channels
    for j = 1:32
        
        if ordered_stim(stim_idx_1,j,PW_idx_1) == 0
            continue
        end
        
        if i == j
            continue
        end
        
        %identify locations for each electrode
        [cuffloc1x, cuffloc1y] = gridLoc(geographic_grid, j);
        [cuffloc2x, cuffloc2y] = gridLoc(geographic_grid, i);
        
        %calculate distances between electrodes
        same_elec1(j,i) = sqrt( (cuffloc1x - cuffloc2x)^2 + (cuffloc1y - cuffloc2y)^2);
    end
    
end

%loop through all channels to calculate distances
for i = 1:32
    
    if ordered_stim(stim_idx_2,i,PW_idx_2) == 0
        continue
    end
    
    %loop through all of the other channels
    for j = 1:32
        
        if ordered_stim(stim_idx_2,j,PW_idx_2) == 0
            continue
        end
        
        if i == j
            continue
        end
        
        %identify locations for each electrode
        [cuffloc1x, cuffloc1y] = gridLoc(geographic_grid, j);
        [cuffloc2x, cuffloc2y] = gridLoc(geographic_grid, i);
        
        %calculate distances between electrodes
        same_elec2(j,i) = sqrt( (cuffloc1x - cuffloc2x)^2 + (cuffloc1y - cuffloc2y)^2);
    end
    
end

%% Plot results

%histcounts for all data
edges = 0:1:8;

[N_diff1, edges] = histcounts(min(diff_elec1), edges);
[N_same1, edges] = histcounts(min(same_elec1), edges);
[N_diff2, edges] = histcounts(min(diff_elec2), edges);
[N_same2, edges] = histcounts(min(same_elec2), edges);

bar_plotting1(:,1) = N_diff1';
bar_plotting1(:,2) = N_same1';
bar_plotting2(:,1) = N_diff2';
bar_plotting2(:,2) = N_same2';

bar_plotting3 = bar_plotting1 + bar_plotting2;

% if false
%     figure(51)
%     plot_bar(bar_plotting1)
%     
%     figure(52)
%     plot_bar(bar_plotting2)
%     
%     figure(53)
%     plot_bar(bar_plotting3)
% end
% 
% %merge all active channel distances, cuff 1:2 & 3:4
% PW100_distances = ((F19_PW100_bar_plotting1+F19_PW100_bar_plotting2)+...
%                     +(F21_PW100_bar_plotting1+F21_PW100_bar_plotting2)+...
%                     +(F22_PW100_bar_plotting1+F22_PW100_bar_plotting2)+...
%                     +(F25_PW100_bar_plotting1+F25_PW100_bar_plotting2)+...
%                     +(F26_PW100_bar_plotting1+F26_PW100_bar_plotting2)+...
%                     +(F34_PW100_bar_plotting1+F34_PW100_bar_plotting2));
% 
% PW500_distances = ((F19_PW500_bar_plotting1+F19_PW500_bar_plotting2)+...
%                     +(F21_PW500_bar_plotting1+F21_PW500_bar_plotting2)+...
%                     +(F22_PW500_bar_plotting1+F22_PW500_bar_plotting2)+...
%                     +(F25_PW500_bar_plotting1+F25_PW500_bar_plotting2)+...
%                     +(F26_PW500_bar_plotting1+F26_PW500_bar_plotting2)+...
%                     +(F34_PW500_bar_plotting1+F34_PW500_bar_plotting2));
% 
% PW1000_distances = ((F19_PW1000_bar_plotting1+F19_PW1000_bar_plotting2)+...
%                     +(F21_PW1000_bar_plotting1+F21_PW1000_bar_plotting2)+...
%                     +(F22_PW1000_bar_plotting1+F22_PW1000_bar_plotting2)+...
%                     +(F25_PW1000_bar_plotting1+F25_PW1000_bar_plotting2)+...
%                     +(F26_PW1000_bar_plotting1+F26_PW1000_bar_plotting2)+...
%                     +(F34_PW1000_bar_plotting1+F34_PW1000_bar_plotting2));
% 
% %PW500_distances1 = (F19_PW500_bar_plotting1+F21_PW500_bar_plotting1+F22_PW500_bar_plotting1...
% %    +F25_PW500_bar_plotting1+F26_PW500_bar_plotting1+F34_PW500_bar_plotting1);
% 
% %PW1000_distances1 = (F19_PW1000_bar_plotting1+F21_PW1000_bar_plotting1+F22_PW1000_bar_plotting1...
% %    +F25_PW1000_bar_plotting1+F26_PW1000_bar_plotting1+F34_PW1000_bar_plotting1);
% 
% % figure(54)
% % fig_title = ['Distance Between Active Channels ,PW = 100us'];
% % plot_bar(PW100_distances, fig_title)
% % figure(55)
% % fig_title = ['Distance Between Active Channels ,PW = 500us'];
% % plot_bar(PW500_distances, fig_title)
% % figure(56)
% % fig_title = ['Distance Between Active Channels ,PW = 1000us'];
% % plot_bar(PW1000_distances, fig_title)
% 
% %merge all active channel distances, cuff 1:2
% %PW100_distances2 = (F19_PW100_bar_plotting2+F21_PW100_bar_plotting2+F22_PW100_bar_plotting2...
% %    +F25_PW100_bar_plotting2+F26_PW100_bar_plotting2+F34_PW100_bar_plotting2);
% %
% %PW500_distances2 = (F19_PW500_bar_plotting2+F21_PW500_bar_plotting2+F22_PW500_bar_plotting2...
% %    +F25_PW500_bar_plotting2+F26_PW500_bar_plotting2+F34_PW500_bar_plotting2);%
% %
% %PW1000_distances2 = (F19_PW1000_bar_plotting2+F21_PW1000_bar_plotting2+F22_PW1000_bar_plotting2...
% %    +F25_PW1000_bar_plotting2+F26_PW1000_bar_plotting2+F34_PW1000_bar_plotting2);
% 
% figure(57)
% fig_title = ['Distance Between Active Channels ,PW = 100us'];
% plot_bar(PW100_distances, fig_title)
% figure(58)
% fig_title = ['Distance Between Active Channels ,PW = 500us'];
% plot_bar(PW500_distances, fig_title)
% figure(59)
% fig_title = ['Distance Between Active Channels ,PW = 1000us'];
% plot_bar(PW1000_distances, fig_title)
% 
% function plot_bar(data, fig_title)
% 
% bar(data)
% title(fig_title)
% ylabel('Counts')
% xlabel('Distance between Active Channels (mm)')
% legend('Opposite Cuff Pair', 'Same Cuff Pair')
% ylim([0 80])
% xticks([1 2 3 4 5 6 7 8 10])
% xticklabels({'0-0.4','0.4-0.8','0.8-1.2', '1.2-1.6', '1.6-2.0', '2.0-2.4', '2.4-2.8', '2.8-3.2', 'NaN'})
% 
% 
% end
%% Identify the distances between old activation patterns and new activation patterns

% %define what the layout of the channel looks like
% geographic_grid = [ 1,5,9,13,17,21,25,29;
%                     2,6,10,14,18,22,26,30;
%                     3,7,11,15,19,23,27,31;
%                     4,8,12,16,20,24,28,32];
%
% %create a logical matrix to compare the amplitudes
% channel_comp = zeros(9,32,6);
%
% %create a matrix of the distances that are calculated
% distances = {};
%
% %create a matrix of nearest neighbor calculations
% near_neigh = {};
%
% %loop through all cuff/pulsewidths
% for pairs = 1:6
%
%     %loop through all stimulation amplitudes
%     for amp_comp = 1:9
%
%         channel_comp(amp_comp,:,pairs) = ordered_stim(amp_comp,:,pairs) == ordered_stim(amp_comp + 1, :, pairs);
%
%         %continue if there was no stimulation in the lower stim
%         if sum(ordered_stim(amp_comp,:,pairs)) == 0
%             continue
%         end
%
%         %If no changes in activation, continune to next pass
%         if sum(channel_comp(amp_comp,:,pairs)) == 32
%             continue
%         end
%
%         %idenfity the locations of activations from the lower amplitude
%         [initial_activation] = find(ordered_stim(amp_comp,:,pairs));
%
%         %identify where changes are made
%         [temp_diff] = find(~channel_comp(amp_comp,:,pairs));
%
%         %identify location where initial activation was present
%         for activation = 1:length(initial_activation)
%             [rows1(activation), columns1(activation)] = gridLoc(geographic_grid, initial_activation(activation));
%         end
%
%         %loop through all places changes have been found
%         for changes = 1:length(temp_diff)
%             [rows2(changes), columns2(changes)] = gridLoc(geographic_grid, temp_diff(changes));
%         end
%
%         %initiate variable for distances
%         temp_distance = [];
%
%         %calculate the euclidian distance between old and new points
%         for new = 1:length(temp_diff)
%             for old = 1:length(initial_activation)
%                 temp_distance(old, new) = sqrt( (rows1(old) - rows2(new))^2 + (columns1(old) - columns2(new))^2);
%             end
%         end
%
%         %continue to next array if no distances were calculated
%         if isempty(temp_distance)
%             continue
%         end
%
%         %save data to an external variable
%         distances{amp_comp,pairs} = temp_distance(:,:);
%
%         [y, x] = size(temp_distance);
%
%         %find the nearest neighbor for each term
%         for new_comps = 1:x
%
%             %don't count zero values
%             if min(temp_distance(:,new_comps)) == 0
%                 [temp_min, index] = min(temp_distance(:,new_comps));
%                 temp_distance(index, new_comps) = 5000;
%             end
%             if min(temp_distance(:,new_comps)) == 5000
%                 temp_neighbor(new_comps) = 0;
%                 continue
%             end
%
%             temp_neighbor(new_comps) = min(temp_distance(:,new_comps));
%         end
%
%         near_neigh{amp_comp,pairs} = temp_neighbor(:);
%
%         %clear variables for next pass
%         clear initial_activation temp_diff rows1 rows2 columns1 columns2 temp_neighbor temp_distance
%
%     end
%
% end


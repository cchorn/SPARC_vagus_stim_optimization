function [SI_combos_list, overlap_combos_list] = SI_search(expmt_list, expmt)
% For each pw config, compare each cuff 1-2 pair trial w/ each cuff 3-4 pair, calculate SI and Cost
% Separate trials into cuff 1-2 pair and cuff 3-4 pair
% Report max SI for 3each pw config (9 total)
%
%  Total Response (TR)        = # of electrode channels with response from A or B (0 to 32)
%  Total Response A & B       = # of electrode channels with response from each cuff pair
%  Total A/B Selective        = # of electrode channels with a response
%                                   from only one cuff pair
%  Total Shared Response (TS) = # of electrode channels with response from both stim
%       parameter configurations (0 to 32)
%  Total Inactivity (TI)      = # electrodes with response from neither
%       stim parameter configurations (0 to 32)
%
%   Efficiency        = (TR - TS)/(# total electrodes)
%   Balance        = (||# TA - #TB||)/TR
%
%  Selectivity Index = A*(Efficiency) + B*(Balance)
%
%  Jonathan Shulgach
%  Last updated: 5/6/2020
%
%

% Set plot_grid to true to view all combinations of channels responding to stim from both cuff pairs
plot_grid=false;
A = 0.5; %efficiency_term
B = 1-A; %balance_term

header = {'Efficiency','Balance','old_SI','SI','Cuff 1-2 Stim',...
    'Cuff 3-4 Stim','Cuff 1-2 PW','Cuff 3-4 PW','# Total Responses','# Total A Responses',...
    '# Total B Responses','# Total A Selective','# Total B Selective','# Total Inactive','# Total Shared Response',...
    'Shared Channels','Cuff 1-2 Grid Chan Response','Cuff 3-4 Grid Chan Response'};

for expmt=expmt % for each day
    SI_combos = [];
    SI_shared_chans_list = [];
    SI_response_1_2_grid_list = [];
    SI_response_3_4_grid_list = [];
    overlap_combos = [];
    exp_data = expmt_list{expmt};
    overlap_combo = 1;
    N_channels = length(exp_data.elec_list);
    
    
    % Get a summary of the data for each trial
    trial_data = get_trial_data(expmt_list, expmt);
    
    if size(exp_data.trial_list,1) > 2 % avoid data with less than 2 pw points
        
        
        for pw_1_2=1:3 % for all cuff 1-2 pw
            for pw_3_4=4:6 % for all cuff 3-4 pw
                
                N_1_2_trials = length(nonzeros(expmt_list{expmt}.stim_hist(pw_1_2, :)));
                N_3_4_trials = length(nonzeros(expmt_list{expmt}.stim_hist(pw_3_4, :)));
                selectGrid_1_2 = reshape(exp_data.selectivityGrid(pw_1_2,1:N_1_2_trials,1:N_channels), N_1_2_trials, N_channels);
                selectGrid_3_4 = reshape(exp_data.selectivityGrid(pw_3_4,1:N_3_4_trials,1:N_channels), N_3_4_trials, N_channels);
                new_combos = zeros((N_1_2_trials*N_3_4_trials),15);
                shared_chans_list = zeros((N_1_2_trials*N_3_4_trials),32);
                response_1_2_grid_list = zeros((N_1_2_trials*N_3_4_trials),32);
                response_3_4_grid_list = zeros((N_1_2_trials*N_3_4_trials),32);
                combo = 1;
                
                % "It is evident that no activity was recorded from many MEA channels, even at maximal stimulus intensities.
                % For discussions of selectivity, the denominator needs to be the maximal number of possible active channels,
                % and not 32." - Bill 1/2/21
                %
                % Get which channels had a possible response, could be all
                % 32
                %N_chans_1_2 = sum(reshape(expmt_list{expmt}.selectivityGrid(pw_1_2,:,:),size(expmt_list{expmt}.stim_hist,2),32),1);
                %N_chans_3_4 = sum(reshape(expmt_list{expmt}.selectivityGrid(pw_3_4,:,:),size(expmt_list{expmt}.stim_hist,2),32),1);
                
                %N_active_chans = length(find((N_chans_1_2+N_chans_3_4)~=0));
                
                % Go through each combination of trials to collect channel
                % response cases. Info collected from each trial combination
                % gets saved to main cell array
                for trial_1_2=1:N_1_2_trials  % for each trial comparison
                    for trial_3_4=1:N_3_4_trials
                        total_response = 0;
                        shared_chan_num = 0;
                        total_A_response = sum(selectGrid_1_2(trial_1_2,:));
                        total_B_response = sum(selectGrid_3_4(trial_3_4,:));
                        total_A_selective = 0;
                        total_B_selective = 0;
                        shared_chans = nan(1,32);
                        
                        stim_1_2 =  exp_data.stim_hist(pw_1_2, trial_1_2);
                        stim_3_4 = exp_data.stim_hist(pw_3_4, trial_3_4);
                        PW_1_2 = exp_data.pulseWidth(pw_1_2);
                        PW_3_4 = exp_data.pulseWidth(pw_3_4);
                        
                        for chan=1:32
                            % Collect the total count of activated channels between both cuff pairs
                            if (selectGrid_1_2(trial_1_2,chan) + selectGrid_3_4(trial_3_4,chan) > 0)
                                total_response = total_response + 1;
                            end
                            % Collect the total number of unique channels activated by cuff 1-2
                            if (selectGrid_1_2(trial_1_2,chan)==1 && selectGrid_3_4(trial_3_4,chan)~=1)
                                total_A_selective = total_A_selective + 1;
                            end
                            % Collect the total number of unique channels activated by cuff 3-4
                            if (selectGrid_1_2(trial_1_2,chan)~=1 && selectGrid_3_4(trial_3_4,chan)==1)
                                total_B_selective = total_B_selective + 1;
                            end
                            % Collect the total number of unique channels activated by both cuff 1-2 and cuff 3-4
                            if (selectGrid_1_2(trial_1_2,chan)==1 && selectGrid_3_4(trial_3_4,chan)==1)
                                shared_chan_num = shared_chan_num + 1; % # Shared response
                                shared_chans(shared_chan_num) = chan;
                            end
                        end
                        
                        if (shared_chan_num<=3 && total_B_response~=0 && total_A_response~=0 && stim_1_2~=0 && stim_1_2~=0)
                            %if (shared_response<=3 && shared_response>0) % 10% overlap of channels
                            
                            overlap_combos(overlap_combo).shared_chans = shared_chans;
                            overlap_combos(overlap_combo).stim_1_2 = stim_1_2;
                            overlap_combos(overlap_combo).stim_3_4 = stim_3_4;
                            overlap_combos(overlap_combo).pw_1_2 = PW_1_2;
                            overlap_combos(overlap_combo).pw_3_4 = PW_3_4;
                            overlap_combo = overlap_combo + 1;
                        end
                        
                        % Plot channel response grid
                        if plot_grid==true
                            fig_title = ['F',expmt_list{expmt}.cohort,' Cuff 1-2 PW', num2str(PW_1_2*1000),'us ',...
                                num2str(stim_1_2),'uA vs Cuff 3-4 PW', num2str(PW_3_4*1000), 'us ', num2str(stim_3_4),'uA'];
                            
                            if PW_1_2==0.5 && PW_3_4 == 0.5
                                if stim_1_2 == 400 && stim_3_4 == 400
                                    plot_response_grid(selectGrid_1_2(trial_1_2,:), selectGrid_3_4(trial_3_4,:),fig_title);
                                    a = input('continue?: ');
                                end
                            end
                        end
                        
                        % ========= SI calculation =====
                        efficiency = (total_response - shared_chan_num)/N_channels;
                        % efficiency = (total_response - shared_chan_num)/N_active_chans;
                        balance = (1- abs(total_A_response - total_B_response)/total_response);
                        SI_Cost =  A*efficiency + B*balance;
                        % =============================
                        
                        new_combos(combo, 1) = efficiency;
                        new_combos(combo, 2) = balance;
                        new_combos(combo, 3) = efficiency*balance;% original SI function
                        new_combos(combo, 4) = SI_Cost; % new SI function
                        new_combos(combo, 5) = exp_data.stim_hist(pw_1_2, trial_1_2);
                        new_combos(combo, 6) = exp_data.stim_hist(pw_3_4, trial_3_4);
                        new_combos(combo, 7) = exp_data.pulseWidth(pw_1_2);
                        new_combos(combo, 8) = exp_data.pulseWidth(pw_3_4);
                        new_combos(combo, 9) = total_response;
                        new_combos(combo, 10) = total_A_response;
                        new_combos(combo, 11) = total_B_response;
                        new_combos(combo, 12) = total_A_selective;
                        new_combos(combo, 13) = total_B_selective;
                        new_combos(combo, 14) = length(exp_data.elec_list) - total_response;% inactive
                        new_combos(combo, 15) = shared_chan_num;
                        
                        % Save shared channel responses for combo data
                        shared_chans_list(combo,:) = shared_chans;
                        
                        % Save channel response grids for combo data
                        response_1_2_grid_list(combo,:) = selectGrid_1_2(trial_1_2,:);
                        response_3_4_grid_list(combo,:) = selectGrid_3_4(trial_3_4,:);
                        
                        combo = combo + 1;
                    end
                end
                new_combos(isnan(new_combos))=0;
                SI_combos = [SI_combos; new_combos];
                
                shared_chans_list(isnan(shared_chans_list))=0;
                SI_shared_chans_list = [SI_shared_chans_list; shared_chans_list];
                SI_response_1_2_grid_list = [SI_response_1_2_grid_list; response_1_2_grid_list];
                SI_response_3_4_grid_list = [SI_response_3_4_grid_list; response_3_4_grid_list];
            end
        end
    end
    
    SI_combo_data = num2cell(SI_combos);
    add_on_data = num2cell(SI_shared_chans_list,2);
    add_on_data2 = num2cell(SI_response_1_2_grid_list,2);
    add_on_data3 = num2cell(SI_response_3_4_grid_list,2);
    
    temp = [SI_combo_data, add_on_data, add_on_data2, add_on_data3];
    SI_combos_list{expmt,1} = [header; temp];
    
    overlap_combos_list{expmt,1} = overlap_combos;
end

end

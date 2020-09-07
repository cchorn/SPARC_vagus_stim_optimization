function [avg_line_data, expmt_list] = get_selec_stim_bound(SI_combos_list, expmt_list)
% This function parses through the combo metadata, selects stim amplitudes
% with equal pulse widths between cuff contact pairs, and classifies data
% points between selective and non-selective.
%
%   INPUTS
%   ===================================================================
%   SI_combo_list :   (1x1 struct) experiment combo data
%   expmt         :   (numeric) experiment data header
%
%   OUTPUTS
%   ===================================================================
%   avg_line_data   :   (1x1 struct) struct holding selective metadata
%
%   Example
%   ===================================================================
% >> [avg_line_data, expmt_list] = get_selec_stim_bound(SI_combos_list, expmt_list)
%
%  Jonathan Shulgach
%  Last updated: 5/7/2020
%

% Changeable Parameters
% ========================================================================
stim_limit = 3000; % change stimulation amplitude limit to parse data points tested below
plot_fig = false; % Set plot_fig to true to view the shared PW selective boundary figures for each animal
plot_grid = false;% Set plot_grid to true to view all combinations of channels responding to stim from both cuff pairs

% ========================================================================

% Reference for struct arcitecture below...

%header = {'SI','Efficiency','Diff. stimA-stimB','Cost','Cuff 1-2 Stim',...
%    'Cuff 3-4 Stim','Cuff 1-2 PW','Cuff 3-4 PW','# Total Responses','# Total A Responses',...
%    '# Total B Responses','# Total A Selective','# Total B Selective','# Total Inactive','# Total Shared Response',...
%    'Shared Channels','Cuff 1-2 Grid Chan Response','Cuff 3-4 Grid Chan Response'};



pw_100_data = [];
pw_100_labels = [];
pw_100_animal = [];
pw_500_data = [];
pw_500_labels = [];
pw_500_animal = [];
pw_1000_data = [];
pw_1000_labels = [];
pw_1000_animal = [];
all_pw_100_min = [];
all_pw_500_min = [];
all_pw_1000_min = [];
all_pw100_sh_chan = [];
all_pw500_sh_chan = [];
all_pw1000_sh_chan = [];

avg_line_data = struct();
for expmt=2:length(SI_combos_list) % for each animal
    k=1;
    % Get a summary of the data for each trial
    trial_data = get_trial_data(expmt_list, expmt); % Can either look through in workspace or print in command window
    N_channels = size(expmt_list{expmt}.elec_list,2);
    
    % Get the minimum threshold for the animal and PW pairs
    expmt_list{expmt}.pw100_min = min([expmt_list{expmt}.minthresh(:,find(expmt_list{expmt}.pulseWidth==0.1))]);
    expmt_list{expmt}.pw400_min = min([expmt_list{expmt}.minthresh(:,find(expmt_list{expmt}.pulseWidth==0.4))]);
    expmt_list{expmt}.pw500_min = min([expmt_list{expmt}.minthresh(:,find(expmt_list{expmt}.pulseWidth==0.5))]);
    expmt_list{expmt}.pw1000_min = min([expmt_list{expmt}.minthresh(:,find(expmt_list{expmt}.pulseWidth==1))]);
    
    % Grab all the data
    SI_data = SI_combos_list{expmt};
    combo_data = cell2mat(SI_data(2:end,1:15));
    
    % Grab shared channel data
    shared_chan_data = cell2mat(SI_data(2:end,16));
    
    % Grab grid response data
    grid_1_2_chan_data = cell2mat(SI_data(2:end,17));
    grid_3_4_chan_data = cell2mat(SI_data(2:end,18));
    
    % Replace nan with zero
    combo_data(isnan(combo_data))=0;
    shared_chan_data(isnan(shared_chan_data))=0;
    
    % Grab all the data where SI output > 0
    %combo_data = combo_data(combo_data(:,1)~=0,:);
    %shared_chan_data = shared_chan_data(combo_data(:,1)~=0,:);
    %grid_1_2_chan_data = grid_1_2_chan_data(combo_data(:,1)~=0,:);
    %grid_3_4_chan_data = grid_3_4_chan_data(combo_data(:,1)~=0,:);
    
    % =====================================================================
    % =====================================================================
    % This portion uses all data points regardless of the PW used. This
    % amount of data may not provide meaningful classification results and
    % isnt used but is in here just in case...
    
    % Collect animal name and selective labels
    %all_TS = 100-combo_data(:,15)*100/32;
    %animal_num = linspace((expmt),(expmt),length(all_TS))';
    %selectivity_labels = classify_SI90_points(all_TS);
    %animal_names = classify_animal_names(animal_num);
    
    %store data points for large classification across animals
    %data = [combo_data(:,5),combo_data(:,6)];
    %all_data = [all_data; [combo_data(:,5),combo_data(:,6)]];
    %selectivity_labels = selectivity_labels;
    %all_labels = [all_labels; selectivity_labels];
    %all_animal_names = [all_animal_names; animal_names];
    
    % =====================================================================
    % =====================================================================
    
    % For each animal collect all the stim amp data with matching cuff pair PW
    for pw_1_2=reshape(unique(combo_data(:,7)),1,3) % for all cuff 1-2 PW values
        for pw_3_4=reshape(unique(combo_data(:,8)),1,3) % for all cuff 3-4 PW values
            
            
            pw_1_2_idx = find(expmt_list{expmt}.pulseWidth(1:3)==pw_1_2);
            pw_3_4_idx = find(expmt_list{expmt}.pulseWidth(4:6)==pw_3_4);
            pw_3_4_idx = 3 + pw_3_4_idx; % need offset
            
            % To plot the selective channel responses with the number of
            % responding channels from each cuff pair (recruitement
            % curves), collect channel responses from tested stim values
            chan_response_data = [];
            
            % Find the channels that responded from cuff 1-2 trials
            chan_response_data(:,1) = expmt_list{expmt}.stim_hist(pw_1_2_idx,:); % cuff 1-2 stim vals
            resp_1_2_grid = reshape(expmt_list{expmt}.selectivityGrid(pw_1_2_idx,:,:),size(expmt_list{expmt}.stim_hist,2),N_channels);
            for i=1:size(resp_1_2_grid,1)
                chan_response_data(i,2) = sum(resp_1_2_grid(i,:));
            end
            
            % Find the channels that responded from cuff 3-4 trials
            chan_response_data(:,3) = expmt_list{expmt}.stim_hist(pw_3_4_idx,:);
            resp_3_4_grid = reshape(expmt_list{expmt}.selectivityGrid(pw_3_4_idx,:,:),size(expmt_list{expmt}.stim_hist,2),N_channels);
            for i=1:size(resp_3_4_grid,1)
                chan_response_data(i,4) = sum(resp_3_4_grid(i,:));
            end
            chan_response_data = chan_response_data(chan_response_data(:,1)~=0,:);
            
            % Find the minimum stim amp where a response from stimulation
            % was recorded for the corresponding PW from both cuff pairs
            min_1_2_stim = min([expmt_list{expmt}.minthresh(:,pw_1_2_idx)]);
            min_3_4_stim = min([expmt_list{expmt}.minthresh(:,pw_3_4_idx)]);
            
            % Check whether data exists before collecting data and making plots
            if ~isempty(min_1_2_stim) && ~isempty(min_3_4_stim)
                
                % Grab all data with matching pw (minus header)
                grid_1_2_chans = []; grid_3_4_chans = []; temp = [];
                temp = combo_data((combo_data(:,7)==pw_1_2 & combo_data(:,8)==pw_3_4),:);
                shared_chans = shared_chan_data((combo_data(:,7)==pw_1_2 & combo_data(:,8)==pw_3_4),:);
                grid_1_2_chans = grid_1_2_chan_data((combo_data(:,7)==pw_1_2 & combo_data(:,8)==pw_3_4),:);
                grid_3_4_chans = grid_3_4_chan_data((combo_data(:,7)==pw_1_2 & combo_data(:,8)==pw_3_4),:);
                
                % Grab all data where # of channel responses is greater than zero from either cuff pair
                shared_chans = shared_chans(temp(:,15) > 0,:);
                grid_1_2_chans = grid_1_2_chans(temp(:,15) > 0,:);
                grid_3_4_chans = grid_3_4_chans(temp(:,15) > 0,:);
                temp = temp(temp(:,15) > 0,:);
                
                % Grab all data where # of channel responses from either cuff pair is 18 or less
                % (no longer reasonable because somatotopy isn't expressed)
                %shared_chans = shared_chans((temp(:,10) <18 & temp(:,11) <= 18),:);
                %grid_1_2_chans = grid_1_2_chans((temp(:,10) <18 & temp(:,11) <= 18),:);
                %grid_3_4_chans = grid_3_4_chans((temp(:,10) <18 & temp(:,11) <= 18),:);
                %temp = temp((temp(:,10) <18 & temp(:,11) <=18),:);
                
                % Grab all data where stim less than limit
                shared_chans = shared_chans(temp(:,5) < stim_limit,:);
                grid_1_2_chans = grid_1_2_chans(temp(:,5) < stim_limit,:);
                grid_3_4_chans = grid_3_4_chans(temp(:,5) < stim_limit,:);
                temp = temp(temp(:,5) < stim_limit,:);
                shared_chans = shared_chans(temp(:,6) < stim_limit,:);
                grid_1_2_chans = grid_1_2_chans(temp(:,6) < stim_limit,:);
                grid_3_4_chans = grid_3_4_chans(temp(:,6) < stim_limit,:);
                temp = temp(temp(:,6) < stim_limit,:);
                
                combo_temp = temp;
                
                % Collect stim amplitude data, and total number of channels responding to either cuff pair stim
                temp_stim_data = [combo_temp(:,5), combo_temp(:,10), combo_temp(:,6), combo_temp(:,11), combo_temp(:,15)];
                temp_stim_data_3chan = temp_stim_data(temp_stim_data(:,5)<=3,:);
                combo_temp_3chan = combo_temp(temp_stim_data(:,5)<=3,:);
                shared_chans_3chan = shared_chans(temp_stim_data(:,5)<=3,:);
                
                % Generate selectivity labels for stim amp data based on
                % 90% selectivity (10% chans w/ shared response)
                TS = 100-combo_temp(:,15)*100/32;
                selectivity_labels = classify_SI90_points(TS);
                
                % Generate list of animal names to accomodate data points
                animal_num = linspace(expmt,expmt,length(TS))';
                animal_names = classify_animal_names(animal_num);
                
                % Generate selectivity labels for 3 shared chan data and
                % animal names
                TS_3chan = 100-combo_temp_3chan(:,15)*100/32;
                selectivity_labels_3chan = classify_SI90_points(TS_3chan);
                animal_num_3chan = linspace(expmt,expmt,length(TS_3chan))';
                animal_names_3chan = classify_animal_names(animal_num_3chan);
                
                % Generate list of minimum stim amp for use with normalizing
                % stim amp data to threshold multiple
                minthresh_stim = [];
                minthresh_stim(1,:) = linspace(min_1_2_stim,min_1_2_stim,length(selectivity_labels));
                minthresh_stim(2,:) = linspace(min_3_4_stim,min_3_4_stim,length(selectivity_labels));
                minthresh_stim = minthresh_stim';
                
                % New list of minimum stim amplitudes for 3 shared chans
                minthresh_stim_3chan = [];
                minthresh_stim_3chan(1,:) = linspace(min_1_2_stim,min_1_2_stim,length(selectivity_labels_3chan));
                minthresh_stim_3chan(2,:) = linspace(min_3_4_stim,min_3_4_stim,length(selectivity_labels_3chan));
                minthresh_stim_3chan = minthresh_stim_3chan';
                
                %temp_idx = find(combo_temp(combo_temp(:,5)<3000 | combo_temp(:,6)<3000));
                %combo_temp = combo_temp(temp_idx,:);
                %shared_chans = shared_chans(temp_idx,:);
                
                % Have all shared channels lists converted to text
                for j=1:size(shared_chans,1)
                    % Display the number of shared channels with each data point
                    %chan_txt{j,1} = num2str(combo_temp(j,15));
                    chan_txt{j,1} = '';
                    % Display list of channel numbers shared (different from above)
                    %chan_txt{j,1} = erase(mat2str(shared_chans(j,:)),' 0');
                end
                
                % Collect PW-specific data across all animals for later use
                if pw_1_2 == 0.1 && pw_3_4 == 0.1
                    pw_100_data = [pw_100_data; temp_stim_data];
                    pw_100_labels = [pw_100_labels; selectivity_labels];
                    pw_100_animal = [pw_100_animal; animal_names];
                    all_pw_100_min = [all_pw_100_min; minthresh_stim];
                    all_pw100_sh_chan = [all_pw100_sh_chan; combo_temp(:,15)];
                    
                    %pw_100_data_3chan = [pw_100_data_3chan; temp_stim_data_3chan];
                    %pw_100_labels_3chan = [pw_100_labels_3chan; selectivity_labels_3chan];
                    %pw_100_animal_3chan = [pw_100_animal_3chan; animal_names_3chan];
                    %all_pw_100_min_3chan = [all_pw_100_min_3chan; minthresh_stim_3chan];
                    %all_pw100_sh_chan_3chan = [all_pw100_sh_chan_3chan; combo_temp_3chan(:,15)];
                    
                elseif pw_1_2 == 0.5 && pw_3_4 == 0.5
                    pw_500_data = [pw_500_data; temp_stim_data];
                    pw_500_labels = [pw_500_labels; selectivity_labels];
                    pw_500_animal = [pw_500_animal; animal_names];
                    all_pw_500_min = [all_pw_500_min; minthresh_stim];
                    all_pw500_sh_chan = [all_pw500_sh_chan; combo_temp(:,15)];
                elseif pw_1_2 == 1 && pw_3_4 == 1
                    pw_1000_data = [pw_1000_data; temp_stim_data];
                    pw_1000_labels = [pw_1000_labels; selectivity_labels];
                    pw_1000_animal = [pw_1000_animal; animal_names];
                    all_pw_1000_min = [all_pw_1000_min; minthresh_stim];
                    all_pw1000_sh_chan = [all_pw1000_sh_chan; combo_temp(:,15)];
                end
                
                % plot chan response data as grid
                if (plot_grid==true && ((pw_1_2==0.1 && pw_3_4==0.1)||(pw_1_2==0.5 && pw_3_4==0.5)||(pw_1_2==1 && pw_3_4==1)))
                    for i=1:size(temp_stim_data,1)
                        new_temp_stim_data = temp_stim_data;
                        new_temp_stim_data(:,3) = new_temp_stim_data(:,3)./minthresh_stim(:,2);
                        new_temp_stim_data(:,1) = new_temp_stim_data(:,1)./minthresh_stim(:,1);
                        fig_title = ['F',expmt_list{expmt}.cohort,' Cuff 1-2 PW', num2str(pw_1_2*1000),'us ',...
                            num2str(new_temp_stim_data(i,1)),'Thr vs Cuff 3-4 PW', num2str(pw_3_4*1000), 'us ',...
                            num2str(new_temp_stim_data(i,3)),'Thr'];
                        
                        plot_response_grid(grid_1_2_chans(i,:), grid_3_4_chans(i,:),fig_title);
                    end
                end
                % Plot channel responses with classifier or grouping
                if (plot_fig==true && ((pw_1_2==0.1 && pw_3_4==0.1)||(pw_1_2==0.5 && pw_3_4==0.5)||(pw_1_2==1 && pw_3_4==1)))
                    fig_title = [expmt_list{expmt}.cohort,' Cuff 1-2 PW ', num2str(pw_1_2*1000),...
                        'us vs Cuff 3-4 PW ', num2str(pw_3_4*1000), 'us',' Shared Response'];
                    try 
                    
                        selective_classifier(temp_stim_data, selectivity_labels, animal_names,...
                            fig_title, 'threshold', minthresh_stim, true, true, chan_txt);
                    
                        selective_classifier(temp_stim_data_3chan, selectivity_labels_3chan, animal_names_3chan,...
                            fig_title, 'threshold', minthresh_stim_3chan, true, true, chan_txt_3chan);
                    
                    
                        selective_plotter(temp_stim_data, combo_temp(:,15), fig_title, 'color', 'amp', minthresh_stim)
                        
                        selective_plotter(temp_stim_data_3chan, combo_temp_3chan(:,15), fig_title, 'color', 'amp', minthresh_stim)
                    catch ME
                        disp('Error Message:')
                        disp(ME.message)
                    end
                    
                    a=input('Continue? (Enter)');
                end
            end
        end
        k = k + 1;
    end
end

%--------------------------------------------------------------------
% Display all the PW data across all animals
%--------------------------------------------------------------------

% run classification of all data points from all animals with PW 500us
selective_classifier(pw_1000_data, pw_1000_labels, pw_1000_animal,...
    'PW 1000us Average Selectivity', 'amp', all_pw_1000_min, false);
a=input('Continue? (Enter)');

% run classification of all data points from all animals with PW 500us
selective_classifier(pw_100_data, pw_100_labels, pw_100_animal,...
    'PW 100us Average Selectivity', 'amp', all_pw_100_min, false);
a=input('Continue? (Enter)');

selective_classifier(pw_500_data, pw_500_labels, pw_500_animal,...
    'PW 500us Average Selectivity', 'amp', all_pw_500_min, false);
a=input('Continue? (Enter)');


selective_plotter(pw_500_data, all_pw500_sh_chan, 'PW 500us Average Channel Selectivity', 'color', 'amp', all_pw_500_min)
a=input('Continue? (Enter)');

selective_plotter(pw_1000_data, all_pw1000_sh_chan, 'PW 1000us Average Channel Selectivity', 'color', 'amp', all_pw_1000_min)
a=input('Continue? (Enter)');


%--------------------------------------------------------------------
% Display all the data with up to 3 shared channels
%--------------------------------------------------------------------
pw_100_data_3chan = pw_100_data(pw_100_data(:,5)<=3,:);
pw_500_data_3chan = pw_500_data(pw_500_data(:,5)<=3,:);
pw_1000_data_3chan = pw_1000_data(pw_1000_data(:,5)<=3,:);
pw_100_labels_3chan = pw_100_labels(pw_100_data(:,5)<=3,:);
pw_500_labels_3chan = pw_500_labels(pw_500_data(:,5)<=3,:);
pw_1000_labels_3chan = pw_1000_labels(pw_1000_data(:,5)<=3,:);
pw_100_animal_3chan = pw_100_animal(pw_100_data(:,5)<=3,:);
pw_500_animal_3chan = pw_500_animal(pw_500_data(:,5)<=3,:);
pw_1000_animal_3chan = pw_1000_animal(pw_1000_data(:,5)<=3,:);
all_pw_100_min_3chan = all_pw_100_min(pw_100_data(:,5)<=3,:);
all_pw_500_min_3chan = all_pw_500_min(pw_500_data(:,5)<=3,:);
all_pw_1000_min_3chan = all_pw_1000_min(pw_1000_data(:,5)<=3,:);
all_pw100_3shchan = all_pw100_sh_chan(pw_100_data(:,5)<=3,:);
all_pw500_3shchan = all_pw500_sh_chan(pw_500_data(:,5)<=3,:);
all_pw1000_3shchan = all_pw1000_sh_chan(pw_1000_data(:,5)<=3,:);


% run classification of all data points from all animals with PW 1000us
selective_classifier(pw_1000_data_3chan, pw_1000_labels_3chan, pw_1000_animal_3chan,...
    'PW 1000us Average Selectivity', 'amp', all_pw_1000_min_3chan, true);
a=input('Continue? (Enter)');

% run classification of all data points from all animals with PW 500us
selective_classifier(pw_500_data_3chan, pw_500_labels_3chan, pw_500_animal_3chan,...
    'PW 500us Average Selectivity', 'amp', all_pw_500_min_3chan, true);
a=input('Continue? (Enter)');

% run classification of all data points from all animals with PW 100us
selective_classifier(pw_100_data_3chan, pw_100_labels_3chan, pw_100_animal_3chan,...
    'PW 100us Average Selectivity', 'amp', all_pw_100_min_3chan, true);
a=input('Continue? (Enter)');


selective_plotter(pw_500_data_3chan, all_pw500_3shchan, 'PW 500us Average Channel Selectivity', 'color', 'amp', all_pw_500_min_3chan)
a=input('Continue? (Enter)');

selective_plotter(pw_1000_data_3chan, all_pw1000_3shchan, 'PW 1000us Average Channel Selectivity', 'color', 'amp', all_pw_1000_min_3chan)
a=input('Continue? (Enter)');



% run classification of all data points from all animals with PW 1000us
%selective_classifier(pw_1000_data, pw_1000_labels, pw_1000_animal,...
%    'PW 1000us Average Selectivity', rec_path,...
%    'Average Selectivity Figures', 'threshold', all_pw_1000_min, false, true);

%selective_plotter(pw_1000_data, all_pw1000_sh_chan, 'PW 1000us Average Channel Selectivity', 'color', 'amp', all_pw_1000_min)

% run classification of all data points from all animals
%coeff = selective_classifier(all_data, all_labels, [x_max y_max],...
%    'Total Average Boundary', rec_path,...
%    'Average Selectivity Figures', 'amp', nan, false);


% Next classify with units of charge
% run classification of all data points from all animals with PW 100us
%selective_classifier(pw_100_data, pw_100_labels, [x_max y_max],...
%    'PW 500us Total Average Boundary', rec_path,...
%    'Average Selectivity Figures', 'charge', 100, false);

% run classification of all data points from all animals with PW 500us
%selective_classifier(pw_500_data, pw_500_labels, [x_max y_max],...
%    'PW 500us Total Average Boundary', rec_path,...
%    'Average Selectivity Figures', 'charge', 500, false);

% run classification of all data points from all animals with PW 1000us
%selective_classifier(pw_1000_data, pw_1000_labels, [4000 4000],...
%    'PW 1000us Total Average Boundary', rec_path,...
%    'Average Selectivity Figures', 'charge', 1000, false);

% run classification of all data points from all animals
%coeff = selective_classifier(all_data, all_labels, [4000 4000],...
%    'Total Average Boundary', rec_path,...
%    'Average Selectivity Figures', 'charge', nan, false);


end

function animal_names = classify_animal_names(animal_num)

animal_names = {};
for i=1:length(animal_num)
    switch animal_num(i)
        case 2
            animal_names{i,1} = 'F21-19';
        case 3
            animal_names{i,1} = 'F22-19';
        case 4
            animal_names{i,1} = 'F19-19';
        case 5
            animal_names{i,1} = 'F25-19';
        case 6
            animal_names{i,1} = 'F26-19';
        case 7
            animal_names{i,1} = 'F34-19';
    end
end

end

function SI90 = classify_SI90_points(varargin)

TS_data = varargin{1};
if length(varargin)>1
    selective_percentage = varargin{2};
else
    selective_percentage = 90;
end

SI90 = {};
for i=1:length(TS_data)
    if TS_data(i)>=selective_percentage
        SI90{i,1} = 'Above 90';
    else
        SI90{i,1} = 'Below 90';
    end
end

end

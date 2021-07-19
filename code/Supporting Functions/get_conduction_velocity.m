function [cv_min_data, cv_max_data] = get_conduction_velocity(varargin)
% Generates conduction velocity figures from response bin data in header file

% Written by: Jonathan Shulgach
% Last updated: 5/6/21

%   INPUTS
%   ===================================================================
%   expmt_list      :   (struct) data file header for experiments
%   expmt           :   (numeric) experiment number, can be list
%   go_plot         :   (bool) (optional) enable/disable plotting figures
%
%   OUTPUTS
%   ===================================================================
%   cv_min_data         :   (1x1 struct) struct containing metadata from
%                          total channel responses and conduction velocity distributions
%
%   cv_max_data         :   (1x1 struct) struct containing metadata from
%                          total channel responses and conduction velocity distributions
%
%   EXAMPLE
%   ===================================================================
%   [cv_min_data, cv_max_data] = get_conduction_velocity(expmt_list, 2:7, true);

%custom enable/disable plotting, set to true to display certain data
plot_min_stim_hist = true;
plot_max_stim_hist = false;
plot_cuff_thresh_color_grid = false;
plot_all_chan_min_stim_hist = false;

% enable/disable displayed text
verbose = false;

analyze_min_stim = true;
analyze_max_stim = true;
N_channels  = 32;
expmt_list = varargin{1};
N_expmt = varargin{2};
if length(varargin)>2
    go_plot = varargin{3};
else
    go_plot=false;
end


stim_data = {'animal','cuff pair','PW','Channel','stim amp','C fiber','a-delta fiber','a-beta fiber','response bin sum','time bins'}; % initialize cell array and start first row with headers
cv_min_data = {'animal','cuff pair','PW','stim type','C fiber','a-delta fiber','response bin sum','time bins','cuff threshold stim','MEA Array','threshold stim trial'}; % initialize cell array and start first row with headers
cv_max_data = {'animal','cuff pair','PW','stim type','C fiber','a-delta fiber','response bin sum','time bins','cuff threshold stim','MEA Array','threshold stim trial'}; % initialize cell array and start first row with headers


for expmt=N_expmt
    cohort = expmt_list{expmt,1}.cohort;
    fprintf("F%s\n", cohort);
    skip_trials = expmt_list{expmt,1}.exclude_trials;
    rms_offset = floor(1000*expmt_list{expmt}.nerve_length/30); % 5/4/21 - rms_offset now set to time index at 30m/s per animal
    
    for session=1:6
        fprintf("session: %d ", session)
        cuff_pair = expmt_list{expmt,1}.cuff_list(session,:);
        
        if analyze_min_stim==true
            
            % --------------- All Channels threshold ----------------
            % Make a record of all channels which have responses and get threshold
            ch_thresh_stim = expmt_list{expmt}.minthresh(:,session);
            
            % Find trial associated with each stim value
            ch_min_trials = NaN(1,length(ch_thresh_stim));
            for i=1:length(ch_thresh_stim)
                if ~isnan(ch_thresh_stim(i))
                    ses_trial = find(expmt_list{expmt}.stim_hist(session,:)==ch_thresh_stim(i));
                    trial = expmt_list{expmt}.trial_list(session,1) + ses_trial - 1;
                    
                    %Might have multiple trials with same stim amplitude, so be sure to
                    %select the one that's not excluded
                    if length(trial)>1
                        trial = trial(~ismember(trial,skip_trials));
                    end
                    ch_min_trials(i) = trial(1);
                end
            end
            
            % --------------- Cuff threshold ----------------
            % Make a record of which channels have responses
            
            % First find lowest stim with a response detected, and get corresponding channels
            cuff_thresh_stim = min(expmt_list{expmt}.minthresh(:,session));
            N_valid_channels = find(expmt_list{expmt}.minthresh(:,session)==cuff_thresh_stim);
            
            % Find trial associated with lowest stim value
            ses_trial = find(expmt_list{expmt}.stim_hist(session,:)==cuff_thresh_stim);
            trial = expmt_list{expmt}.trial_list(session,1) + ses_trial - 1;
            
            %Might have multiple trials with same stim amplitude, so be sure to
            %select the one that's not excluded
            if length(trial)>1
                trial = trial(~ismember(trial,skip_trials));
            end
            
            PW = find_pw(expmt_list{expmt,1}, trial);
            fprintf("PW %d\n", PW)
            A_delta_fiber_count = 0;
            A_beta_fiber_count = 0;
            C_fiber_count = 0;
            MEA_chan_data = {};
            min_bin_sum = zeros(1,size(expmt_list{expmt}.cv(10,trial).resp_bins,2));
            for chan=1:N_channels
                %fprintf("chan %d\n", chan);
                
                if ismember(chan,N_valid_channels)
                    MEA_chan_data{chan,1} = expmt_list{expmt}.cv(chan,trial).resp_bins;
                    %MEA_chan_data{chan,1} = expmt_list{expmt}.cv(chan,ch_min_trials(chan)).resp_bins;
                    resp_bins = expmt_list{expmt}.cv(chan,trial).resp_bins;
                    time_bins = expmt_list{expmt}.cv(chan,trial).time_bins;
                    min_bin_sum = min_bin_sum + resp_bins;
                    %fprintf("chan %d: %d\n",chan,sum(resp_bins))
                    [C_fiber_resp, A_delta_fiber_resp, A_beta_fiber_resp] = check_fiber_type(resp_bins, time_bins, verbose);
                        
                    if C_fiber_resp==0 && A_delta_fiber_resp==0
                        a=1;
                    end
                    C_fiber_count = C_fiber_count + C_fiber_resp;
                    A_delta_fiber_count =  A_delta_fiber_count +  A_delta_fiber_resp;
                    A_beta_fiber_count = A_beta_fiber_count + A_beta_fiber_resp;
                else
                    C_fiber_resp = false;
                    A_delta_fiber_resp = false;
                    A_beta_fiber_resp = false;
                end
                if plot_cuff_thresh_color_grid==true
                    fig_title = ['F', expmt_list{expmt}.cohort,' ', cuff_pair, ' PW', num2str(PW*1000), 'us Electrode CV at Cuff Threshold'];
                    color_cv_response_grid(chan, fig_title, 11, C_fiber_resp, A_delta_fiber_resp, A_beta_fiber_resp)
                end
            end
            % Update master data cell table
            temp_data = {['F',cohort], cuff_pair, PW, 'min', C_fiber_count, A_delta_fiber_count, min_bin_sum, time_bins, cuff_thresh_stim, MEA_chan_data, trial};
            cv_min_data = [cv_min_data; temp_data];
        end
        
        
        if analyze_max_stim==true
            % ------------------- Max stim responses ---------------------------
            % Make a record of which channels have responses
            % First find highest stim with a response detected, and get corresponding channels
            temp_chan_list = 1:32;
            N_valid_channels = temp_chan_list(~ismember(temp_chan_list,expmt_list{expmt}.exclude_channels));
            
            % Find trial associated with highest stim value, make sure its
            % valid
            
            % -- Set limit var to true below to implement max stim limit
            use_max_stim_lim = true;
            if use_max_stim_lim==true
                max_stim = expmt_list{expmt}.stim_hist(session, 1);
                if expmt_list{expmt}.stim_hist(session, 1) > 3000
                    trial = expmt_list{expmt}.trial_list(session,2);
                else
                    trial = expmt_list{expmt}.trial_list(session,1);
                end
            else
                trial = expmt_list{expmt}.trial_list(session,1);
            end
            
            PW = find_pw(expmt_list{expmt}, trial);
            A_delta_fiber_count = 0;
            C_fiber_count = 0;
            MEA_chan_data = {};
            max_bin_sum = zeros(1,size(expmt_list{expmt}.cv(10,trial).resp_bins,2));
            for chan=1:N_channels
                A_beta_fiber_resp = false;
                if ismember(chan,N_valid_channels)
                    MEA_chan_data{chan,1} = expmt_list{expmt}.cv(chan,trial).resp_bins;
                    time_bins = expmt_list{expmt}.cv(chan,trial).time_bins;
                    resp_bins = expmt_list{expmt}.cv(chan,trial).resp_bins;
                    max_bin_sum = max_bin_sum + resp_bins;
                    %fprintf("chan %d: %d\n",chan,sum(resp_bins))
                    [C_fiber_resp, A_delta_fiber_resp, A_beta_fiber_resp] = check_fiber_type(resp_bins, time_bins, verbose);
                    C_fiber_count = C_fiber_count + C_fiber_resp;
                    A_delta_fiber_count =  A_delta_fiber_count +  A_delta_fiber_resp;
                else
                    C_fiber_resp = false;
                    A_delta_fiber_resp = false;
                end
                if plot_cuff_thresh_color_grid==true
                    fig_title = ['F', expmt_list{expmt}.cohort,' ', cuff_pair, ' PW', num2str(PW*1000), 'us Electrode CV at Max Stim'];
                    color_cv_response_grid(chan, fig_title, 11, C_fiber_resp, A_delta_fiber_resp, A_beta_fiber_resp)
                end
            end
            % Update master data cell table
            temp_data = {['F',cohort], cuff_pair, PW, 'max', C_fiber_count, A_delta_fiber_count, max_bin_sum, time_bins, max_stim, MEA_chan_data, trial};
            cv_max_data = [cv_max_data; temp_data];
            
        end
        
        if plot_cuff_thresh_color_grid==true
            a=input('Next PW? (Enter): ');
        end
        
    end
end

if plot_all_chan_min_stim_hist==true && go_plot==true
    for pw = unique([stim_data{2:end,3}])
        fig_title = ['PW',num2str(pw*1000),'us Conduction Velocity Responses at Threshold'];
        
        % Histograms for stim response classified by fiber type response per chan
        stim_amp_histogram(stim_data, pw, 45, fig_title);
        
        a=input('Next PW? (Enter): ');
    end
    
end

%cuff_types = ['1:2';'3:4'];
if plot_min_stim_hist==true && go_plot==true
    %for cuff=1:size(cuff_types,1)
    for pw = unique([cv_min_data{2:end,3}])
        fig_title = ['PW',num2str(pw*1000),'us Conduction Velocity Responses at Threshold'];
        
        % Histogram for CVs detected from 0 - 25m/s
        cv_histogram(cv_min_data, 'min', pw, 45, fig_title, false);
        
        % Histogram for fiber type responses at cuff threshold
        fiber_type_histogram(cv_min_data,'min', pw, 60, fig_title);
        
        % Selective Fiber Type Histogram, shows 3 columns {"C", "a-delta", "both fibers"}
        fig_title = ['PW',num2str(pw*1000),'us Selective Fiber Type CV Responses at Threshold'];
        selective_fiber_type_histogram(cv_min_data,'min', pw, 60, fig_title);
        
        % Histogram for stim amplitude from responses at threshold
        %min_stim_amp_histogram(cv_min_data, pw, 60, fig_title, false);
        
        a=input('Next PW? (Enter): ');
    end
    %end
end
if plot_max_stim_hist==true && go_plot==true
    %for cuff=1:size(cuff_types,1)
    for pw = unique([cv_max_data{2:end,3}])
        %fig_title = [cuff_types(cuff,:), ' PW',num2str(pw*1000),'us Conduction Velocity Responses at Max Stim'];
        fig_title = ['PW',num2str(pw*1000),'us Conduction Velocity Responses at Max Stim'];
        
        % Histogram for CVs detected from 0 - 25m/s
        cv_histogram(cv_max_data, 'max', pw, 180, fig_title, false);
        
        % Histogram for fiber type responses at cuff threshold
        fiber_type_histogram(cv_max_data,'max', pw, 300, fig_title);
        
        % Selective Fiber Type Histogram, shows 3 columns {"C", "a-delta", "both fibers"}
        fig_title = ['PW',num2str(pw*1000),'us Selective Fiber Type CV Responses at Max Stim'];
        selective_fiber_type_histogram(cv_max_data,'max', pw, 300, fig_title);
        
        a=input('Next PW? (Enter): ');
    end
    %end
end

end


function color_cv_response_grid(chan, fig_title, fig_id, C_fiber_resp, A_delta_fiber_resp, A_beta_fiber_resp)

% Custom array for displaying channel outputs in appropriate subplots to
% match Utah MEA chnanel configuration on screen
elec_chan_map = [...
    1,  9, 17, 25,...
    2, 10, 18, 26,...
    3, 11, 19, 27,...
    4, 12, 20, 28,...
    5, 13, 21, 29,...
    6, 14, 22, 30,...
    7, 15, 23, 31,...
    8, 16, 24, 32];

N_channels  = 32;
sub_size = numSubplots(N_channels);

figure(fig_id)
grid_fig = subplot(sub_size(1),sub_size(2), elec_chan_map(chan));
sgtitle(fig_title);
if C_fiber_resp==1 && A_delta_fiber_resp==1 && A_beta_fiber_resp==1
    rectangle('Position',[0,1,2,1],'FaceColor','r','EdgeColor','k','LineWidth',1)
    rectangle('Position',[0,0,1,1],'FaceColor','b','EdgeColor','k','LineWidth',1)
    rectangle('Position',[1,0,1,1],'FaceColor','g','EdgeColor','k','LineWidth',1)
elseif C_fiber_resp==0 && A_delta_fiber_resp==0 && A_beta_fiber_resp==0
    rectangle('Position',[0,0,2,2],'FaceColor','w','EdgeColor','k','LineWidth',1)
elseif C_fiber_resp==1 && A_delta_fiber_resp==1
    patch([0,0,2],[0,2,2],'red');
    patch([0,2,2],[0,0,2],'blue');
elseif C_fiber_resp==1 && A_beta_fiber_resp==1
    patch([0,0,2],[0,2,2],'red');
    patch([0,2,2],[0,0,2],'green');
elseif A_beta_fiber_resp==1 && A_delta_fiber_resp==1
    patch([0,0,2],[0,2,2],'blue');
    patch([0,2,2],[0,0,2],'green');
elseif C_fiber_resp==1
    rectangle('Position',[0,0,2,2],'FaceColor','r','EdgeColor','k','LineWidth',1)
elseif A_delta_fiber_resp==1
    rectangle('Position',[0,0,2,2],'FaceColor','b','EdgeColor','k','LineWidth',1)
elseif A_beta_fiber_resp==1
    rectangle('Position',[0,0,2,2],'FaceColor','.5 0 .5','EdgeColor','k','LineWidth',1)
else
    rectangle('Position',[0,0,2,2],'FaceColor','w','EdgeColor','k','LineWidth',1)
end
set(gca,'XTick',[])
set(gca,'YTick',[])
end

function stim_amp_histogram(data, pw, y_lim, fig_title)
% Generate separate histograms, overlayed on each other
data = data(2:end,:);
temp_data = data([data{:,3}]==pw,:);% grab data for specific pw

%remove data points at max stim
stim_lim = 3000;
temp_data = temp_data([temp_data{:,5}] < stim_lim,:);

% find unique animal names and rearrange for consistent order across figures
temp_names = temp_data(:,1)';
[names,i] = unique(temp_names);
animal_names = temp_names(sort(i));

% Group stim amps into respective fiber type response datasets
c_fiber_data = temp_data([temp_data{:,6}]==true,:);
a_delta_fiber_data = temp_data([temp_data{:,7}]==true,:);
a_beta_fiber_data = temp_data([temp_data{:,8}]==true,:);

Y_1 = [c_fiber_data{:,5}];
%Y_1 = [c_fiber_data{:,5}]/norm([c_fiber_data{:,5}]);
Y_2 = [a_delta_fiber_data{:,5}];
%Y_2 = [a_delta_fiber_data{:,5}]/norm([a_delta_fiber_data{:,5}]);

h1 = histogram(Y_1);
hold on
h2 = histogram(Y_2);
h1.Normalization = 'probability';
h1.BinWidth = 100;
h2.Normalization = 'probability';
h2.BinWidth = 100;
hold off;

%cvFig = gcf;
title(fig_title)
legend('C Fiber','a-delta Fiber');
%assign_animal_colors(cvFig, animal_names);
ylabel('% of Responses');
xlabel('Stim (uA)')


end

function min_stim_amp_histogram(data, pw, y_lim, fig_title,shrink)

% Generate separate histograms, overlayed on each other
data = data(2:end,:);
temp_data = data([data{:,3}]==pw,:);% grab data for specific pw

%remove data points at max stim
%stim_lim = 3000;
%temp_data = temp_data([temp_data{:,5}] < stim_lim,:);

% find unique animal names and rearrange for consistent order across figures
temp_names = temp_data(:,1)';
[names,i] = unique(temp_names);
animal_names = temp_names(sort(i));

% Get percentage of fiber-selective threshold responses
only_c_fiber = temp_data([temp_data{:,6}]==0,:);
only_a_delta_fiber = temp_data([temp_data{:,5}]==0,:);
[N_response,~] = size(temp_data);

%fprintf('%% distribution of Selective Fiber Types at Threshold: \nC Fiber: %d\na-delta Fiber: %d\n',
display([100*size(only_c_fiber,1)/N_response,"%"])
display([100*size(only_a_delta_fiber,1)/N_response,"%"])

% Group stim amps into respective fiber type response datasets
c_fiber_data = temp_data([temp_data{:,5}]>0,:);
a_delta_fiber_data = temp_data([temp_data{:,6}]>0,:);
%a_beta_fiber_data = temp_data([temp_data{:,7}]==true,:);

Y_1 = [c_fiber_data{:,9}];
fprintf('Mean C Fiber threhsold stim amplitude: %0.2d\n',mean(Y_1));
%Y_1 = [c_fiber_data{:,5}]/norm([c_fiber_data{:,5}]);
Y_2 = [a_delta_fiber_data{:,9}];
fprintf('Mean a-delta threshold stim amplitude: %0.2d\n',mean(Y_2));
%Y_2 = [a_delta_fiber_data{:,5}]/norm([a_delta_fiber_data{:,5}]);

figure(33);
h1 = histogram(Y_1);
hold on
h2 = histogram(Y_2);
h1.Normalization = 'probability';
h1.BinWidth = 50;
h2.Normalization = 'probability';
h2.BinWidth = 50;
%h1.FaceColor = [43/255 140/255 190/255];
h1.FaceAlpha=0.8;
%h2.FaceColor = [241/255 238/255 246/255];
h2.FaceAlpha = 0.8;
hold off;

cvFig = gcf;
%yticks([1 4 6 10])
YTickLabel = cvFig.CurrentAxes.YTick;% = cvFig.CurrentAxes.YTick*100;
yticklabels(string(YTickLabel*100));
title(fig_title)
legend('C Fiber','a-delta Fiber');
%assign_animal_colors(cvFig, animal_names);
ylabel('% of Selective Fiber Type Response');
xlabel('Stim (uA)')

end

function fiber_type_histogram(data, stim_type, pw, y_lim, fig_title)
% Generate histogram with the 2 fiber types
% grab min stim data
temp_data = data(strcmp(data(:,4),stim_type),:); %choose min stim
%temp_data = temp_data(strcmp(temp_data(:,2),cuff_pair),:); % choose cuff pair
temp_data = temp_data([temp_data{:,3}]==pw,:);% grab data for specific pw

% find unique animal names and rearrange for consistent order across figures
temp_names = temp_data(:,1)';
[~,i] = unique(temp_names);
animal_names = temp_names(sort(i));

% combine response data between cuff pairs so only one solid color for both
for i=1:length(animal_names)
    row_idx = find(strcmp(temp_data(:,1),animal_names(1,i)));% Get data for this animal
    resp_data{i,1} = sum(reshape([temp_data{row_idx,7}], length(temp_data{row_idx(1),7}), length(row_idx)),2)';
end
%resp_data = temp_data(:,7);

% response bins may be different lengths, so we need to pad the shortest vectors with zeros
maxSize = max(cellfun(@numel, resp_data));               % Get the maximum vector size
fcn = @(x) [x zeros(1, maxSize-numel(x))];             % Create an anonymous function
rmat = cellfun(fcn, resp_data, 'UniformOutput', false);  % Pad each cell with NaNs
temp_Y = vertcat(rmat{:});

Y = [sum(temp_Y(:,1:6),2)'; sum(temp_Y(:,7:end),2)'];


%Y = [[resp_data{:,1}];  [resp_data{:,2}]]';
figure(30);
set(gcf,'Position',[100 100 500 500])
X = categorical({'C fibers','A\delta fibers'});
X = reordercats(X,{'C fibers','A\delta fibers'});
bar(X,Y,0.4,'stacked')
cvFig = gcf;
cvFig.Children.YLim = [0 y_lim];
assign_animal_colors(cvFig, animal_names);
legend(animal_names);
title(fig_title)
ylabel('# of Responses');
xlabel('Stim (uA)')

end

function cv_histogram(data, stim_type, pw, y_lim, fig_title, shrink)

% grab min stim data
temp_data = data(strcmp(data(:,4),stim_type),:); %choose min stim
%temp_data = temp_data(strcmp(temp_data(:,2),cuff_pair),:); % choose cuff pair
temp_data = temp_data([temp_data{:,3}]==pw,:);% grab data for specific pw

% find unique animal names and rearrange for consistent order across figures
temp_names = temp_data(:,1)';
[names,i] = unique(temp_names);
animal_names = temp_names(sort(i));

% combine response data between cuff pairs so only one solid color for both
for i=1:length(animal_names)
    row_idx = find(strcmp(temp_data(:,1),animal_names(1,i)));% Get data for this animal
    resp_data{i,1} = sum(reshape([temp_data{row_idx,7}], length(temp_data{row_idx(1),7}), length(row_idx)),2)';
end
%resp_data = temp_data(:,7);

% response bins may be different lengths, so we need to pad the shortest vectors with zeros
maxSize = max(cellfun(@numel, resp_data));               % Get the maximum vector size
fcn = @(x) [x zeros(1, maxSize-numel(x))];             % Create an anonymous function
rmat = cellfun(fcn, resp_data, 'UniformOutput', false);  % Pad each cell with NaNs
temp_Y = vertcat(rmat{:});

%temp_Y = reshape([temp_data{:,7}],50,size(temp_data,1))';
if shrink==true
    D_num = 2;
    for i=1:size(temp_Y,1)
        X = linspace(0,25,50/D_num);
        Y(i,:) = sum(reshape(temp_Y(i,:),D_num,[]),1);
    end
else
    X = linspace(0.5,size(temp_Y,2)/2,size(temp_Y,2));
    Y = temp_Y;
end

% also need to sort data order if legen order is set
figure(31);
bar(X,Y,'stacked');
cvFig = gcf;
cvFig.Children.YLim = [0 y_lim];
title(fig_title)
legend(animal_names);
assign_animal_colors(cvFig, animal_names);
ylabel('# of Responses');
xlabel('CV (m/sec)')
xlim([0,30]);

end

function selective_fiber_type_histogram(data, stim_type, pw, y_lim, fig_title)
% Generate histogram with 3 columns: C fibers, a-delta fibers, and both fiber types

% grab min stim data
temp_data = data(strcmp(data(:,4),stim_type),:); %choose min stim
%temp_data = temp_data(strcmp(temp_data(:,2),cuff_pair),:); % choose cuff pair
temp_data = temp_data([temp_data{:,3}]==pw,:);% grab data for specific pw

% find unique animal names and rearrange for consistent order across figures
temp_names = temp_data(:,1)';
[~,i] = unique(temp_names);
animal_names = temp_names(sort(i));

% combine response data between cuff pairs so only one solid color for both
c_fiber_sum = [];
a_delta_fiber_sum = [];
both_fiber_sum = [];
for i=1:length(animal_names)
    row_idx = find(strcmp(temp_data(:,1),animal_names(1,i)));% Get data for this animal
    resp_data{i,1} = sum(reshape([temp_data{row_idx,7}], length(temp_data{row_idx(1),7}), length(row_idx)),2)';
    
    % Collect data from both cuff pairs for each animal
    MEA_resp_data = temp_data(row_idx, 10);
    for j=1:size(MEA_resp_data,1)

        c_fiber_ct = 0;
        a_delta_fiber_ct = 0;
        a_beta_fiber_ct = 0;
        both_fiber_ct = 0;            

        % Find Selective CV response behavior for each channel
        for chan=1:size(MEA_resp_data{j,1},1)

            temp_MEA_data = MEA_resp_data{j,1}{chan,1};
            %fprintf("chan%d resp: %d\n", chan, sum(temp_MEA_data));
            if ~isempty(temp_MEA_data)
                
                %Check that there is channel data
                if sum(temp_MEA_data(1:7))>0 && sum(temp_MEA_data(8:end))==0
                    c_fiber_ct = c_fiber_ct + 1;
                end
                if sum(temp_MEA_data(8:end))>0 && sum(temp_MEA_data(1:7))==0
                    a_delta_fiber_ct = a_delta_fiber_ct + 1;
                end
                if sum(temp_MEA_data(1:7))>0 && sum(temp_MEA_data(8:end))>0
                    both_fiber_ct = both_fiber_ct + 1;
                end
            end
        end 
        c_fiber(j) = c_fiber_ct;
        a_delta_fiber(j) = a_delta_fiber_ct;
        both_fiber(j) = both_fiber_ct;

    end
    c_fiber_sum(i) = sum(c_fiber);
    a_delta_fiber_sum(i) = sum(a_delta_fiber);
    both_fiber_sum(i) = sum(both_fiber);
    
end

% response bins may be different lengths, so we need to pad the shortest vectors with zeros
%maxSize = max(cellfun(@numel, resp_data));               % Get the maximum vector size
%fcn = @(x) [x zeros(1, maxSize-numel(x))];             % Create an anonymous function
%rmat = cellfun(fcn, resp_data, 'UniformOutput', false);  % Pad each cell with NaNs
%temp_Y = vertcat(rmat{:});

%
% old_Y = [sum(temp_Y(:,1:6),2)'; sum(temp_Y(:,7:end),2)'];
Y = [c_fiber_sum; a_delta_fiber_sum; both_fiber_sum];
%Y = [[resp_data{:,1}];  [resp_data{:,2}]]';
figure(36);
set(gcf,'Position',[100 100 500 500])
X = categorical({'C fibers','A\delta fibers', 'Both Fibers'});
X = reordercats(X,{'C fibers','A\delta fibers','Both Fibers'});
bar(X,Y,0.4,'stacked')
cvFig = gcf;
cvFig.Children.YLim = [0 y_lim];
assign_animal_colors(cvFig, animal_names);
legend(animal_names);
title(fig_title)
ylabel('# of Responses');
xlabel('Stim (uA)')

end

function [C_fiber_resp, A_delta_fiber_resp, A_beta_fiber_resp] = check_fiber_type(resp_bins, time_bins, verbose)

if verbose==true
    fprintf("\nChan %d - Total bins with response: %d\n", chan, sum(resp_bins));
end

% Assuming bin step size is 0.5m/s, first element is time 0, second element
% is 0.5m/s, third is 1m/s, etc. 7th element is 3m/s, which is the max CV
% for vagus C-Fibers (Lee Fisher - 6/23/20)
if sum(resp_bins(1:7))>0
    if verbose==true
        fprintf("C-fiber: YES");
    end
    C_fiber_resp = true;
else
    if verbose==true
        fprintf("C-fiber: NO");
    end
    C_fiber_resp = false;
end
if sum(resp_bins(7:end))>0
    if verbose==true
        fprintf("\t aDelta-fiber: YES\n");
    end
    A_delta_fiber_resp = true;
else
    if verbose==true
        fprintf("\t aDelta-fiber:  NO\n");
    end
    A_delta_fiber_resp = false;
end
A_beta_fiber_resp = false;

end

function [FaceColor] = assign_animal_colors(cvFig, animal_names)
fig_h = cvFig.Children;
if length(fig_h)>1
    bar_h = fig_h(end).Children;
else
    bar_h = fig_h.Children;
end
N_animals = size(animal_names,2);
for i=1:N_animals
    animal = animal_names{i};
    switch(animal)
        % Weird
        case 'F21-19'
            bar_h(N_animals-i+1).FaceColor = [4/255 90/255 141/255];
        case 'F22-19'
            bar_h(N_animals-i+1).FaceColor = [43/255 140/255 190/255];
        case 'F29-19'
            bar_h(N_animals-i+1).FaceColor = [116/255 169/255 207/255];
        case 'F25-19'
            bar_h(N_animals-i+1).FaceColor = [166/255 189/255 219/255];
        case 'F26-19'
            bar_h(N_animals-i+1).FaceColor = [208/255 209/255 230/255];
        case 'F34-19'
            bar_h(N_animals-i+1).FaceColor = [241/255 238/255 246/255];
    end
end

end

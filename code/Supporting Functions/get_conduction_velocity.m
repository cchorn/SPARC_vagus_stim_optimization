function get_conduction_velocity(varargin)
% Generates conduction velocity figures from response bin data in header file

% Written by: Jonathan Shulgach
% Last updated: 8/13/20

plot_min_stim_hist = true;
plot_max_stim_hist = true;
plot_cuff_thresh_color_grid = false;
verbose = false;
analyze_min_stim = true;
analyze_max_stim = true;
N_channels  = 32;
rms_offset = 15;%milliseconds
expmt_list = varargin{1};
N_expmt = varargin{2};

cv_data = {'animal','cuff pair','PW','stim type','C fiber','a-delta fiber','response bin sum','time bins'}; % initialize cell array and start first row with headers

for expmt=N_expmt
    cohort = expmt_list{expmt,1}.cohort;
    fprintf("F%s\n", cohort);
    skip_trials = expmt_list{expmt,1}.exclude_trials;
    
    for session=1:6
        fprintf("session: %d ", session)
        cuff_pair = expmt_list{expmt,1}.cuff_list(session,:);
        
        if analyze_min_stim==true
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
            C_fiber_count = 0;
            min_bin_sum = zeros(1,size(expmt_list{expmt}.cv(10,trial).resp_bins,2));
            for chan=1:N_channels
                %fprintf("chan %d\n", chan);
                A_beta_fiber_resp = false;
                if ismember(chan,N_valid_channels)
                    resp_bins = expmt_list{expmt}.cv(chan,trial).resp_bins;
                    time_bins = expmt_list{expmt}.cv(chan,trial).time_bins;
                    min_bin_sum = min_bin_sum + resp_bins;
                    fprintf("chan %d: %d\n",chan,sum(resp_bins))
                    [C_fiber_resp, A_delta_fiber_resp, A_beta_fiber_resp] = check_fiber_type(resp_bins, time_bins, verbose);
                    C_fiber_count = C_fiber_count + C_fiber_resp;
                    A_delta_fiber_count =  A_delta_fiber_count +  A_delta_fiber_resp;
                else
                    C_fiber_resp = false;
                    A_delta_fiber_resp = false;
                end
                if plot_cuff_thresh_color_grid==true
                    fig_title = ['F', expmt_list{expmt}.cohort,' ', cuff_pair, ' PW', num2str(PW*1000), 'us Electrode CV at Cuff Threshold'];
                    color_cv_response_grid(chan, fig_title, 11, C_fiber_resp, A_delta_fiber_resp, A_beta_fiber_resp)
                end
            end
            % Update master data cell table
            temp_data = {['F',cohort], cuff_pair, PW, 'min', C_fiber_count, A_delta_fiber_count, min_bin_sum, time_bins};
            cv_data = [cv_data; temp_data];
        end
        
        
        if analyze_max_stim==true
            % ------------------- Max stim responses ---------------------------
            % Make a record of which channels have responses
            % First find highest stim with a response detected, and get corresponding channels
            temp_chan_list = 1:32;
            N_valid_channels = temp_chan_list(~ismember(temp_chan_list,expmt_list{expmt}.exclude_channels));
            
            % Find trial associated with highest stim value
            trial = expmt_list{expmt}.trial_list(session,1);
            
            PW = find_pw(expmt_list{expmt}, trial);
            A_delta_fiber_count = 0;
            C_fiber_count = 0;
            max_bin_sum = zeros(1,size(expmt_list{expmt}.cv(10,trial).resp_bins,2));
            for chan=1:N_channels
                A_beta_fiber_resp = false;
                if ismember(chan,N_valid_channels)
                    
                    time_bins = expmt_list{expmt}.cv(chan,trial).time_bins;
                    resp_bins = expmt_list{expmt}.cv(chan,trial).resp_bins;
                    max_bin_sum = max_bin_sum + resp_bins;
                    fprintf("chan %d: %d\n",chan,sum(resp_bins))
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
            temp_data = {['F',cohort], cuff_pair, PW, 'max', C_fiber_count, A_delta_fiber_count, max_bin_sum, time_bins};
            cv_data = [cv_data; temp_data];
            
        end
        
        if plot_cuff_thresh_color_grid==true
            a=input('Next PW? (Enter): ');
        end
        
    end
end

%cuff_types = ['1:2';'3:4'];
if plot_min_stim_hist==true
    %for cuff=1:size(cuff_types,1)
        for pw = unique([cv_data{2:end,3}])
            fig_title = ['PW',num2str(pw*1000),'us Conduction Velocity Responses at Threshold'];
            
            % Histogram for CVs detected from 0 - 25m/s
            cv_histogram(cv_data, 'min', pw, fig_title, false);
            
            % Histogram for fiber type responses at cuff threshold
            fiber_type_histogram(cv_data,'min', pw, fig_title);
            
            a=input('Next PW? (Enter): ');
        end
    %end
end
if plot_max_stim_hist==true
    %for cuff=1:size(cuff_types,1)
        for pw = unique([cv_data{2:end,3}])
            %fig_title = [cuff_types(cuff,:), ' PW',num2str(pw*1000),'us Conduction Velocity Responses at Max Stim'];
            fig_title = ['PW',num2str(pw*1000),'us Conduction Velocity Responses at Max Stim'];
            
            % Histogram for CVs detected from 0 - 25m/s
            cv_histogram(cv_data, 'max', pw, fig_title, false);
            
            % Histogram for fiber type responses at cuff threshold
            fiber_type_histogram(cv_data,'max', pw, fig_title);
            
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

function fiber_type_histogram(data, stim_type, pw, fig_title)

% grab min stim data
temp_data = data(strcmp(data(:,4),stim_type),:); %choose min stim
%temp_data = temp_data(strcmp(temp_data(:,2),cuff_pair),:); % choose cuff pair

temp_data = temp_data([temp_data{:,3}]==pw,:);% grab data for specific pw
temp_names = temp_data(:,1)';
[names,i] = unique(temp_names);
animal_names = temp_names(sort(i));

% combine response data between cuff pairs so only one solid color for both
for i=1:length(animal_names)
    row_idx = find(strcmp(temp_data(:,1),animal_names(1,i)));
    resp_data{i,1} = sum([temp_data{row_idx,5}]);
    resp_data{i,2} = sum([temp_data{row_idx,6}]);
end

Y = [[resp_data{:,1}];  [resp_data{:,2}]]';
figure(30);
set(gcf,'Position',[100 100 500 500])
X = categorical({'C fibers','A\delta fibers'});
X = reordercats(X,{'C fibers','A\delta fibers'});
bar(X,Y',0.4,'stacked')
legend(animal_names);
title(fig_title)
ylabel('# of Responses');

end

function cv_histogram(data, stim_type, pw, fig_title, shrink)
D_num = 2;
% grab min stim data
temp_data = data(strcmp(data(:,4),stim_type),:); %choose min stim
%temp_data = temp_data(strcmp(temp_data(:,2),cuff_pair),:); % choose cuff pair

temp_data = temp_data([temp_data{:,3}]==pw,:);% grab data for specific pw
temp_names = temp_data(:,1)';
[names,i] = unique(temp_names);
animal_names = temp_names(sort(i));

% combine response data between cuff pairs so only one solid color for both
for i=1:length(animal_names)
    row_idx = find(strcmp(temp_data(:,1),animal_names(1,i)));
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
    for i=1:size(temp_Y,1)
        X = linspace(0,25,50/D_num);
        Y(i,:) = sum(reshape(temp_Y(i,:),D_num,[]),1);
    end
else
    X = linspace(0,size(temp_Y,2)/2,size(temp_Y,2));
    Y = temp_Y;
end
figure(31);
bar(X,Y','stacked')
title(fig_title)
legend(animal_names);
ylabel('# of Responses');
xlabel('CV (m/sec)')

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

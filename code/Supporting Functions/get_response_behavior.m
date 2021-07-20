function [resp_data] = get_response_behavior(varargin)
% Generates figures from response bin data in header file
%
% Inputs:
%   expmt_list : (struct) header file
%   N_expmt : (int) experiment number to analyze
%
% Written by: Jonathan Shulgach
% Last updated: 4/22/21

verbose = false;
analyze_min_stim = true;
analyze_max_stim = false;
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
            % --------------- All Channels threshold ----------------
            % Make a record of all channels which have responses and get threshold
            cuff_thresh_stim = expmt_list{expmt}.minthresh(:,session);
            
            % Find trial associated with each stim value
            trials = NaN(1,length(cuff_thresh_stim));
            for i=1:length(cuff_thresh_stim)
                if isnan(cuff_thresh_stim(i))
                    trials(i) = NaN;
                else
                    ses_trial = find(expmt_list{expmt}.stim_hist(session,:)==cuff_thresh_stim(i));
                    trial = expmt_list{expmt}.trial_list(session,1) + ses_trial - 1;
                    
                    %Might have multiple trials with same stim amplitude, so be sure to
                    %select the one that's not excluded
                    if length(trial)>1
                        trial = trial(~ismember(trial,skip_trials));
                    end
                    
                    trials(i) = trial(1);
                end
            end
            PW = find_pw(expmt_list{expmt,1}, trial);
            fprintf("PW %d\n", PW)
            A_delta_fiber_count = 0;
            C_fiber_count = 0;
            min_bin_sum = zeros(1,size(expmt_list{expmt}.cv(10,trial).resp_bins,2));
            
            for chan=1:N_channels
                %fprintf("chan %d\n", chan);
                A_beta_fiber_resp = false;
                if ~isnan(trials(chan))
                    resp_bins = expmt_list{expmt}.cv(chan,trials(chan)).resp_bins;
                    time_bins = expmt_list{expmt}.cv(chan,trials(chan)).time_bins;
                    min_bin_sum = min_bin_sum + resp_bins;
                    fprintf("chan %d: %d\n",chan,sum(resp_bins))
                    [C_fiber_resp, A_delta_fiber_resp, A_beta_fiber_resp] = check_fiber_type(resp_bins, time_bins, verbose);
                    C_fiber_count = C_fiber_count + C_fiber_resp;
                    A_delta_fiber_count =  A_delta_fiber_count +  A_delta_fiber_resp;
                else
                    C_fiber_resp = false;
                    A_delta_fiber_resp = false;
                end
         %       if plot_cuff_thresh_color_grid==true
         %           fig_title = ['F', expmt_list{expmt}.cohort,' ', cuff_pair, ' PW', num2str(PW*1000), 'us Electrode CV at Cuff Threshold'];
         %           color_cv_response_grid(chan, fig_title, 11, C_fiber_resp, A_delta_fiber_resp, A_beta_fiber_resp)
         %       end
            end
         %   % Update master data cell table
         %   temp_data = {['F',cohort], cuff_pair, PW, 'min', C_fiber_count, A_delta_fiber_count, min_bin_sum, time_bins, cuff_thresh_stim};
         %   cv_data = [cv_data; temp_data];
            
        end
        
        
        if analyze_max_stim==true
            
        end
        
    end
    
end
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

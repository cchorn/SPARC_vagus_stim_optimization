function opt_params = get_opt_params(SI_combos_list, expmt_list, N_expmt, varargin)
% 1) Obtain and plot SI figure with associated stim params, for each PW
% combination (mesh or surface plot)
%
% 2)  Create a bar plot with total number of selectively activated
% electrode channels fro highest SI outputs over all PW, for each animal
%
%   INPUTS
%   ===================================================================
% SI_combos_list     : (Nx1 cell) cell array with response combination data
% expmt_list         : (Nx1 cell) header file
% N_expmt            : (1xN num)  vector of experiments to analyze
% selective_chan_lim : (num) (Optional) number of channels allowed for
%                          response overlap
% generate_plot      : (bool) (optional) enable/disable plots, default enable
% 
%   OUTPUTS
%   ===================================================================
%   opt_params       : (struct) with fields for each PW, containing SI optimizations    
%
%   EXAMPLE
%   ===================================================================
%   opt_params_list = get_opt_params(SI_combos_list, expmt_list, 2:7, 3);

% ---- changeable params ----

plot_SI_color   = false;
plot_PW_bar  = false;
plot_PW_bar_selec_chan = true;
plot_PW_bar_bal = false;
plot_PW_bar_bal_selec_chan = false;
non_zero_resp = false;
plot_max_stim = false;
plot_min_stim = false;

global selective_chan_lim
if ~isempty(varargin)
    selective_chan_lim = varargin{1};
else
    selective_chan_lim = 3;
end
if length(varargin)>1
    go_plot = varargin{2};
else
    go_plot = true;
end

selec_percent = 100*(round((1-selective_chan_lim/32)*20)/20); % 5 percent accuraccy

% ---------------------------
opt_params = struct();
bar_animals  = cell(0);
%rec_path     = 'C:\Users\shulgac\Documents\MATLAB\Data';

for anim_expmt=N_expmt
    
    data = SI_combos_list{anim_expmt};
    combo_data = cell2mat(data(2:end,:));
    combo_data(isnan(combo_data))=0;
    expmt = anim_expmt - 1;
    bar_animals{1,(expmt)} = expmt_list{anim_expmt,1}.cohort;
    
    header_data = data(1,:);
    % {'Response Power','Efficiency','Diff. stimA-stimB','SI Score','Cuff 1-2 Stim',
    % 'Cuff 3-4 Stim','Cuff 1-2 PW','Cuff 3-4 PW','# Total resps','# Total A resps'
    % '# Total B resps','# Total A Selective','# Total B Selective','# Total Inactive','# Total Shared resp'};
    
    % Collect SI data for each PW pair combo
    for pw_1_2=reshape(unique(combo_data(:,7)),1,3)         % for all cuff 1-2 pw
        for pw_3_4=reshape(unique(combo_data(:,8)),1,3)     % for all cuff 3-4 pw
            
            % Only for shared PW values
            if pw_1_2==pw_3_4
                %file_directory = [rec_path, '\F', expmt_list{expmt,1}.cohort];
                temp = combo_data((combo_data(:,7)==pw_1_2 & combo_data(:,8)==pw_3_4),:); % All data with matching pw (minus header)
                
                %temp_data = temp; %comment this out if setting limit
                
                % Remove highest stim value as optimal stim param
                %------------------------------------------------------------
                cuff_1_2_lim = 3000; cuff_3_4_lim = 3000;
                temp_data = temp((temp(:,5)<=cuff_1_2_lim & temp(:,6)<=cuff_3_4_lim),:);
                %temp_data = temp((temp(:,5)<max(temp(:,5))&temp(:,6)<max(temp(:,6))),:);
                %--------------------------------------------------------------
                temp_SI = temp_data(:,4);
                
                if non_zero_resp==true
                    % enforce both sides having to respond for valid max SI
                    for i=1:length(temp_SI)
                        if temp_data(i,10)==0 || temp_data(i,11)==0 || temp_data(i,12)==0 || temp_data(i,13)==0
                            temp_data(:,4)=0;
                        end
                    end
                end
                
                
                % Data collection for Max SI score
                title_text = ['F',expmt_list{anim_expmt,1}.cohort,' PW ', num2str(pw_1_2*1000),'\mus SI Score'];
                temp = optimize(temp_data, 'Optimize', 'SI', 'Parameter', 'Selective Response', 'max');
                opt_params.(['PW_',num2str(pw_1_2*1000)]).max_SI.max_SI(1,expmt) = temp;
                
                temp = optimize(temp_data, 'Optimize', 'Max Stim', 'Parameter', 'Selective Response');
                opt_params.(['PW_',num2str(pw_1_2*1000)]).max_stim.max_SI(1,expmt) = temp;
                %opt_params.(['PW_',num2str(pw_1_2*1000)]).threshold.max_SI = optimize(temp_data, expmt, 'Plot SI', title_text, 'Optimize', 'SI', 'Parameter', 'Selective Response', 'max');
                
                temp = optimize(temp_data, 'Optimize', 'Threshold', 'Parameter', 'Selective Response');
                opt_params.(['PW_',num2str(pw_1_2*1000)]).threshold.max_SI(1,expmt) = temp;
                
                % Data collection for max SI score given selective channel overlap
                temp_data_selec_chan = temp_data(temp_data(:,15)<=selective_chan_lim,1:15);
                opt_params.(['PW_',num2str(pw_1_2*1000)]).max_SI.max_SI_selec_chan(1,expmt) = optimize(temp_data_selec_chan, 'Optimize', 'SI', 'Parameter', 'Selective Response', 'max');
                
                % Data collection for max SI score given selective channel overlap
                temp_data_selec_chan = temp_data(temp_data(:,15)<=selective_chan_lim,1:15);
                opt_params.(['PW_',num2str(pw_1_2*1000)]).max_SI.max_SI_selec_chan(1,expmt) = optimize(temp_data_selec_chan, 'Optimize', 'SI', 'Parameter', 'Selective Response', 'max');
                
                % Data collection for best balanced response
                opt_params.(['PW_',num2str(pw_1_2*1000)]).threshold.max_bal(1, expmt) = optimize(temp_data, 'Optimize', 'Balance', 'Parameter', 'Selective Response', 'max');
                
                % Data collection for best balanced response and no more than 3 channel overlap
                opt_params.(['PW_',num2str(pw_1_2*1000)]).threshold.max_bal_selec_chan(1,expmt) = optimize(temp_data_selec_chan, 'Optimize', 'Balance', 'Parameter', 'Selective Response', 'max');
                
                % ----------------------------------------------------------------------------
                
                if plot_SI_color==true
                    x_3chan = temp_data_selec_chan(:,5);
                    y_3chan = temp_data_selec_chan(:,6);
                    z_3chan = temp_data_selec_chan(:,4);
                    %SI_surface_plot(temp_data(:,5), temp_data(:,6), temp_data(:,4), fig_title);
                    %SI_surface_plot(temp_data(:,5), temp_data(:,6), temp_data(:,4), title_text, x_3chan, y_3chan, z_3chan);
                    SI_surface_plot2(temp_data, title_text, x_3chan, y_3chan, z_3chan);
                    a = input("\nContinue? (Enter): ");
                end
                
                % =======================================================
            end
        end
    end
    
    %fprintf('===================================================\n')
end

%fprintf('\nPW100us best cost: \n')
%fprintf('%0.1f |', opt_params.PW_100.max_SI.SI)
%fprintf('\nPW500us best cost: \n')
%fprintf('%0.1f |', opt_params.PW_500.max_SI.SI)
%fprintf('\nPW1000us best cost: \n')
%fprintf('%0.1f |', opt_params.PW_1000.max_SI.SI)
%fprintf('\n');



pw_list=reshape(unique(combo_data(:,7)),1,3);

for pw=1:length(pw_list) % dataset collection for each animal
    
    %max stim data
    param_data = filter_data(opt_params.(['PW_',num2str(pw_list(pw)*1000)]).max_stim.max_SI);
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).max_stim_data          = ([param_data.A_selective; param_data.B_selective])';
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).max_stim_overlap_data  = ([param_data.shared_chans])';
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).max_stim_stim_data     = [param_data.stim_data_A; param_data.stim_data_B];
    
    param_data = filter_data(opt_params.(['PW_',num2str(pw_list(pw)*1000)]).threshold.max_SI);
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).threshold_data          = ([param_data.A_selective; param_data.B_selective])';
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).threshold_overlap_data  = ([param_data.shared_chans])';
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).threshold_stim_data     = [param_data.stim_data_A; param_data.stim_data_B];
    
    param_data = filter_data(opt_params.(['PW_',num2str(pw_list(pw)*1000)]).max_SI.max_SI);
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).data          = ([param_data.A_selective; param_data.B_selective])';
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).overlap_data  = ([param_data.shared_chans])';
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).stim_data     = [param_data.stim_data_A; param_data.stim_data_B];
    
    param_data = filter_data(opt_params.(['PW_',num2str(pw_list(pw)*1000)]).max_SI.max_SI_selec_chan);
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).data_selec_chan    = ([param_data.A_selective; param_data.B_selective])';
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).overlap_data_selec_chan = ([param_data.shared_chans])';
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).stim_data_selec_chan   = [param_data.stim_data_A; param_data.stim_data_B];
    
    param_data = filter_data(opt_params.(['PW_',num2str(pw_list(pw)*1000)]).threshold.max_bal);
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).data_bal    = ([param_data.A_selective; param_data.B_selective])';
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).overlap_data_bal = ([param_data.shared_chans])';
    
    param_data = filter_data(opt_params.(['PW_',num2str(pw_list(pw)*1000)]).threshold.max_bal_selec_chan);
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).data_bal_selec_chan    = ([param_data.A_selective; param_data.B_selective])';
    bar_data.(['PW_',num2str(pw_list(pw)*1000)]).overlap_data_bal_selec_chan = ([param_data.shared_chans])';
    
end
if plot_min_stim==true && go_plot==true
    figure(35)
    stacked_bar_plot(bar_data.PW_100.threshold_data, bar_data.PW_100.threshold_overlap_data, bar_animals, 'PW100us Min Stim Total # Channels', true)
    figure(36)
    stacked_bar_plot(bar_data.PW_500.threshold_data, bar_data.PW_500.threshold_overlap_data, bar_animals, 'PW500us Min Stim Total # Channels', true)
    figure(37)
    stacked_bar_plot(bar_data.PW_1000.threshold_data, bar_data.PW_1000.threshold_overlap_data, bar_animals, 'PW1000us Min Stim Total # Channels', true)    
end
if plot_max_stim==true && go_plot==true
    figure(38)
    stacked_bar_plot(bar_data.PW_100.max_stim_data, bar_data.PW_100.max_stim_overlap_data, bar_animals, 'PW100us Max Stim Max Chan Response', true)
    figure(39)
    stacked_bar_plot(bar_data.PW_500.max_stim_data, bar_data.PW_500.max_stim_overlap_data, bar_animals, 'PW500us Max Stim Max Chan Response', true)
    figure(40)
    stacked_bar_plot(bar_data.PW_1000.max_stim_data, bar_data.PW_1000.max_stim_overlap_data, bar_animals, 'PW1000us Max Stim Max Chan Response', true)
end
if plot_PW_bar==true && go_plot==true
    figure(41)
    %bar_plot(pw100_bar_data, bar_animals, 'PW100us Max SI Total # Channels', true)
    stacked_bar_plot(bar_data.PW_100.data, bar_data.PW_100.overlap_data, bar_animals, 'PW100us Max SI Max Chan Response', true)
    figure(42)
    %bar_plot(pw500_bar_data, bar_animals, 'PW500us Max SI Total # Channels', true)
    stacked_bar_plot(bar_data.PW_500.data, bar_data.PW_500.overlap_data, bar_animals, 'PW500us Max SI Max Chan Response', true)
    figure(43)
    %bar_plot(pw1000_bar_data, bar_animals, 'PW1000us Max SI Total # Channels', true)
    stacked_bar_plot(bar_data.PW_1000.data, bar_data.PW_1000.overlap_data, bar_animals, 'PW1000us Max SI Max Chan Response', true)
end
if plot_PW_bar_selec_chan==true && go_plot==true
    figure(44);
    stacked_bar_plot(bar_data.PW_100.data_selec_chan, bar_data.PW_100.overlap_data_selec_chan, bar_animals,...
        ['PW100us Max SI Max Chan Response with ',num2str(selec_percent),'% Selectivity'], true);
    figure(45);
    stacked_bar_plot(bar_data.PW_500.data_selec_chan, bar_data.PW_500.overlap_data_selec_chan, bar_animals,...
        ['PW500us Max SI Max Chan Response with ',num2str(selec_percent),'% Selectivity'], true);
    figure(46);
    stacked_bar_plot(bar_data.PW_1000.data_selec_chan, bar_data.PW_1000.overlap_data_selec_chan, bar_animals,...
        ['PW1000us Max SI Max Chan Response with ',num2str(selec_percent),'% Selectivity'], true);
end
if plot_PW_bar_bal==true && go_plot==true
    figure(47)
    stacked_bar_plot(bar_data.PW_100.data_bal, bar_data.PW_100.overlap_data_bal, bar_animals, 'PW100us Total # Responding Channels at Optimized Balance', true)
    figure(48)
    stacked_bar_plot(bar_data.PW_500.data_bal, bar_data.PW_500.overlap_data_bal, bar_animals, 'PW500us Total # Responding Channels at Optimized Balance', true)
    figure(49)
    stacked_bar_plot(bar_data.PW_1000.data_bal, bar_data.PW_1000.overlap_data_bal, bar_animals, 'PW1000us Total # Responding Channels at Optimized Balance', true)
    
end
if plot_PW_bar_bal_selec_chan==true && go_plot==true
    figure(50)
    stacked_bar_plot(bar_data.PW_100.data_bal_selec_chan, bar_data.PW_100.overlap_data_bal_selec_chan, bar_animals, 'PW100us Total # Responding Channels at Optimized Balance', true)
    figure(51)
    stacked_bar_plot(bar_data.PW_500.data_bal_selec_chan, bar_data.PW_500.overlap_data_bal_selec_chan, bar_animals, 'PW500us Total # Responding Channels at Optimized Balance', true)
    figure(52)
    stacked_bar_plot(bar_data.PW_1000.data_bal_selec_chan, bar_data.PW_1000.overlap_data_bal_selec_chan, bar_animals, 'PW1000us Total # Responding Channels at Optimized Balance', true)    
end

end

function opt_params = filter_data(opt_params)

for i=1:size(opt_params,2)
    if isempty(opt_params(i).SI)
        opt_params(i).SI=0;
    end
    if isempty(opt_params(i).shared_chans)
        opt_params(i).shared_chans=0;
    end
    if isempty(opt_params(i).A_selective)
        opt_params(i).A_selective=0;
    end
    if isempty(opt_params(i).B_selective)
        opt_params(i).B_selective=0;
    end
    if isempty(opt_params(i).A_resp)
        opt_params(i).A_resp=0;
    end
    if isempty(opt_params(i).B_resp)
        opt_params(i).B_resp=0;
    end
    if isempty(opt_params(i).stim_data_A)
        opt_params(i).stim_data_A=0;
    end
    if isempty(opt_params(i).stim_data_B)
        opt_params(i).stim_data_B=0;
    end
end

end

function stacked_bar_plot(bar_data, bar_overlap_data, bar_animals, figure_title, dual)
global selective_chan_lim

%figure; % can create figure outside function beforehand
% Merge both columns into one column to zip the selective response data
temp_zip = [bar_data(:,1), bar_data(:,2)].';
zip_data = temp_zip(:);

% Fill out shared channel response data, duplicating every element
shared_resp = repelem(bar_overlap_data,2);

% Actual binned data
Y = [zip_data, shared_resp];

temp_Y = reshape([reshape(Y,2,[]); zeros(1,numel(Y)/2)],[],2);
if dual==true
    resp_after_selec_chan = zeros(size(temp_Y,1),1);
    for i=1:length(resp_after_selec_chan)
        if temp_Y(i,2) > selective_chan_lim
            resp_after_selec_chan(i,1) = temp_Y(i,2)-selective_chan_lim;
            temp_Y(i,2) = selective_chan_lim;
        end
    end
    
    temp_Y = [temp_Y, resp_after_selec_chan];
end

b_fig = bar(temp_Y,'stacked');
if length(b_fig)==3
    b_fig(3).FaceColor = [0.7 0.7 0.7];
    legend('Selective Response',[num2str(selective_chan_lim),' Shared Channels'],['More than ',num2str(selective_chan_lim),' Shared Channels']);
else
    legend('Selective Response',[num2str(selective_chan_lim),' Shared Channels']);%'More than 3 Shared Channels')    
end
set(gca,'XTick',1.5:3:19,'XTickLabel',bar_animals);  %# Modify axes
%hB = bar(bar_data);
%hAx = gca;
%set(gca, 'XTickLabel',bar_animals, 'XTick',1:numel(bar_animals))
title(figure_title);
ylim([0 32])
xlabel('Animals');
ylabel('# of Selectively Channels Activated');

end

function bar_plot(bar_data, bar_animals, figure_title, dual)
hB = bar(bar_data);
hAx = gca;
set(gca, 'XTickLabel',bar_animals, 'XTick',1:numel(bar_animals))
title(figure_title);
if dual==true
    %labels = {'Cuff 1-2','Cuff 3-4';'Cuff 1-2','Cuff 3-4';'Cuff 1-2','Cuff 3-4';...
    %'Cuff 1-2','Cuff 3-4';'Cuff 1-2','Cuff 3-4';'Cuff 1-2','Cuff 3-4'};
    legend('Cuff 1-2','Cuff 3-4');
end
ylim([0 32])
xlabel('Animals');
ylabel('# of Selectively Channels Activated');

end

function SI_surface_plot(x, y, z, fig_title, varargin)
% param: x data column (cuff stim)
% param: y data column (cuff stim)
% param: z data column (SI score)

figure(2)
%surf(x,y,v); % Add surf
dt = delaunayTriangulation(x,y) ;
tri = dt.ConnectivityList ;
xi = dt.Points(:,1) ;
yi = dt.Points(:,2) ;
F = scatteredInterpolant(x,y,z);
zi = F(xi,yi) ;
trisurf(tri,xi,yi,zi)
view(2) % top view
%view(3) % angle view
shading interp %flat
drawnow
hold on;

%All data points
plot3(x, y, z, 'ob','MarkerSize', 1, 'LineWidth', 3);

%Point with max SI
[maxcostval, maxcostidx] = get_opt_SI(x, y, z, 'min');
plot3(x(maxcostidx), y(maxcostidx), maxcostval, 'or', 'MarkerSize', 5, 'LineWidth', 2.0);

%Point with max SI given 3 chan overlap
if isempty(varargin)==false
    x_3chan = varargin{1}; y_3chan = varargin{2}; z_3chan = varargin{3};
    [maxcostval_3chan, maxcostidx_3chan] = get_opt_SI(x_3chan, y_3chan, z_3chan, 'min');
    plot3(x_3chan(maxcostidx_3chan), y_3chan(maxcostidx_3chan), maxcostval_3chan, 'om', 'MarkerSize', 7, 'LineWidth', 2.0);
end

fprintf("max position| stim 1-2: %d, stim 3-4: %d",x(maxcostidx), y(maxcostidx))
hold off;
ylabel('\muA');
xlabel('\muA');
title(fig_title);

end

function SI_surface_plot2(temp_data, fig_title, varargin)
% param: x data column (cuff stim)
% param: y data column (cuff stim)
% param: z data column (SI score)

figure(2)

x = temp_data(:,5);
y = temp_data(:,6);
z = temp_data(:,4);

%surf(x,y,v); % Add surf
dt = delaunayTriangulation(x,y) ;
tri = dt.ConnectivityList ;
xi = dt.Points(:,1) ;
yi = dt.Points(:,2) ;
F = scatteredInterpolant(x,y,z);
zi = F(xi,yi) ;
trisurf(tri,xi,yi,zi)
view(2) % top view
%view(3) % angle view
shading interp %flat
drawnow
hold on;

% collect non selective data
cover = temp_data((temp_data(:,15)>3),:);
% Get outermost points from data
if size(cover,1)>=3
    k = convhull(cover(:,5),cover(:,6));
    % shade out the region of >3 shared channels
    patch(cover(k,5),cover(k,6),repmat(max(z),length(k),1),[0.5 0.5 0.5])
else
    plot3(cover(:,5),cover(:,6),repmat(max(z),size(cover,1),1),...
        'Color',[0.5 0.5 0.5],'LineWidth',5);
end

%All data points
plot3(x, y, z, 'ob','MarkerSize', 1, 'LineWidth', 3);

%Point with max SI
[maxcostval, maxcostidx] = get_opt_SI(x, y, z, 'min');
plot3(x(maxcostidx), y(maxcostidx), maxcostval, 'or', 'MarkerSize', 5, 'LineWidth', 2.0);

%Point with max SI given 3 chan overlap
if isempty(varargin)==false
    x_3chan = varargin{1}; y_3chan = varargin{2}; z_3chan = varargin{3};
    [maxcostval_3chan, maxcostidx_3chan] = get_opt_SI(x_3chan, y_3chan, z_3chan, 'min');
    plot3(x_3chan(maxcostidx_3chan), y_3chan(maxcostidx_3chan), maxcostval_3chan, 'om', 'MarkerSize', 7, 'LineWidth', 2.0);
end

% adjust axis range limits
if max(x)>1000
    %xlim([0 1000]);
end
if max(y)>1000
    %ylim([0 1000]);
end

fprintf("max position| stim 1-2: %d, stim 3-4: %d",x(maxcostidx), y(maxcostidx))
hold off;
ylabel('\muA 3-4');
xlabel('\muA 1-2');
title(fig_title);

end

function [opt_params] = optimize(data, varargin)

if isempty(varargin)
    error("No command entered!");
    return
end

temp_cmd_idx = strcmp('Optimize', varargin);
cmd_idx = find(temp_cmd_idx==1) + 1;
if length(varargin)<=cmd_idx
    error("Missing parameter input for Command");
else
    command = varargin{cmd_idx};
end
temp_param_idx = strcmp('Parameter', varargin);
param_idx = find(temp_param_idx==1) + 1;
param_cmd_idx = param_idx + 1;

switch(varargin{param_idx})
    case 'Selective Response'
        x = data(:,12); y =  data(:,13);
    otherwise
        disp('Only Selective Response for now...')
end

if length(varargin)<=param_cmd_idx
    param_cmd = 'max';
else
    param_cmd = varargin{param_cmd_idx};
end

if cmd_idx==0
    disp("No 'Optimize' command entered");
else
    switch(command)
        case 'SI'
            z = data(:,4);
            [~, idx] = get_opt_SI(x, y, z, param_cmd);
            if sum(ismember(varargin, 'Plot SI'))>0
                temp_plot_idx = strcmp('Plot SI', varargin);
                plot_idx = double(find(temp_plot_idx==1) + 1);
                if isempty(plot_idx)
                    fig_title = '';
                else
                    fig_title = varargin{plot_idx};
                end
                SI_surface_plot(data(:,5), data(:,6), data(:,4), fig_title);
                a = input("\nContinue? (Enter): ");
            end
        case 'Balance'
            [idx] = get_balanced_response(x, y, param_cmd);
            
        case 'Max Stim'
            % get params from max stim trial
            x = data(:,5); y =  data(:,6);          
            [idx] = get_max_stim(x, y, param_cmd);
        case 'Threshold'
            x = data(:,5); y =  data(:,6);
            z = data(:,9);
            [idx, x_stim, y_stim] = get_min_stim(x, y, z);            
        otherwise
            beep;
            disp('Invalid command');
    end
    
    opt_params.SI = data(idx,4);
    opt_params.shared_chans = data(idx,15);
    opt_params.A_selective = data(idx,12);
    opt_params.B_selective = data(idx,13);
    opt_params.A_resp = data(idx,10);
    opt_params.B_resp = data(idx,11);
    opt_params.stim_data_A = data(idx, 5);
    opt_params.stim_data_B = data(idx, 6);
    
end

end

function [idx] = get_max_stim(x,y,varargin)
 cuff_1_2_max = x==max(x);
 cuff_3_4_max = y==max(y);
 [~,idx] = max(cuff_1_2_max + cuff_3_4_max);
end

function [maxcostidx,x_stim,y_stim] = get_min_stim(x,y,z)
 %non zero values
 z(z==0) = nan;
 
 if sum(z==min(z))>1
        logic = double(z == min(z)); %convert logic array of max values to double
        % Find max SI with highest parameter values
        select_chan_sum = (x + y).*logic; % get sum of selective channels
        [~,maxcostidx] = max(select_chan_sum); %find new index of highest val where total channels responding is higher
        maxcostval = max(z(maxcostidx)); % then get max cost val at new index
 else
     [maxcostval, maxcostidx] = max(z);
 end
 x_stim = x(maxcostidx);
 y_stim = y(maxcostidx);
end

function [idx] = get_balanced_response(x, y, varargin)

diff_resp = abs(x - y);
if sum(diff_resp==min(diff_resp))>1
    logic_idx = diff_resp==min(min(diff_resp));
    if strcmp(varargin{1}, 'min')
        % If multiple indices have a response with the same number
        stim_sq = sqrt((x.*logic_idx).^2 + (y.*logic_idx).^2);% get sum of stim amplitudes
        stim_sq(stim_sq==0) = NaN; % replace zeros with NaN
        [~,idx] = min(stim_sq); % Goal is to find SI with smallest stim amps
    else
        select_chan_sum = (x + y).*logic_idx; % get sum of selective channels
        [~,idx] = max(select_chan_sum); %find new index of highest val where total channels responding is higher
    end
else
    [~, idx] = min(diff_resp);
end

end

function [maxcostval, maxcostidx] = get_opt_SI(x, y, z, varargin)

% If multiple SI scores share the same maximum value, this finds the row index of the
% highest SI value and saves the parameters from that run
if sum(z==max(z))>1
    logic = double(z == max(z)); %convert logic array of max values to double
    if strcmp(varargin{1}, 'min')
        % See if multiple points have same SI value and if so, choose the one with
        % the lowest parameter values
        stim_sq = sqrt((x.*logic).^2 + (y.*logic).^2);% get sum of stim amplitudes
        stim_sq(stim_sq==0) = NaN; % replace zeros with NaN
        [~,maxcostidx] = min(stim_sq); % Goal is to find SI with smallest stim amps
        maxcostval = z(maxcostidx);
    else
        % Find max SI with highest parameter values
        select_chan_sum = (x + y).*logic; % get sum of selective channels
        [~,maxcostidx] = max(select_chan_sum); %find new index of highest val where total channels responding is higher
        maxcostval = max(z(maxcostidx)); % then get max cost val at new index
    end
else
    [maxcostval, maxcostidx] = max(z);
end

end
function opt_params_list = get_opt_params(SI_combos_list, expmt_list, N_expmt)
% 1) Obtain and plot SI figure with associated stim params, for each PW
% combination (mesh or surface plot)
%
% 2)  Create a bar plot with total number of selectively activated
% electrode channels fro highest SI outputs over all PW, for each animal
% Example:
% >> opt_params_list = get_opt_params(SI_combos_list, expmt_list, 2:7)

% ---- changeable params ----

plot_mesh    = false;
plot_color   = true;
plot_all_bar = false;
plot_PW_bar  = false;
non_zero_resp = true;
select_3chan = true;
plot_3chan_overlap = false;
%save_fig2    = false;
%save_fig3    = false;

opt_params_list = struct();
bar_animals  = cell(0);
bar_data     = zeros(length(SI_combos_list)-1, 2);
%rec_path     = 'C:\Users\shulgac\Documents\MATLAB\Data';
%grid_space   = 40;
for expmt=N_expmt % for each day
    data = SI_combos_list{expmt};
    
    combo_data = cell2mat(data(2:end,:));
    combo_data(isnan(combo_data))=0;
    j = 1;
    stim_1_2_list = [];
    stim_3_4_list = [];
    pw_1_2_list = [];
    pw_3_4_list = [];
    si_list = [];
    
    bar_animals{1,(expmt-1)} = expmt_list{expmt,1}.cohort;
    
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
                
                % Remove highest stim value as optimal stim param
                %cuff_1_2_lim = 1500; cuff_3_4_lim = 1500;
                %temp_data = temp((temp(:,5)<cuff_1_2_lim & temp(:,6)<cuff_3_4_lim),:);
                %temp_data = temp((temp(:,5)<max(temp(:,5))&temp(:,6)<max(temp(:,6))),:);
                %stim_temp_logic = double((temp_data(:,5)<max(temp_data(:,5))&temp_data(:,6)<max(temp_data(:,6))));
                %temp_data = temp_data(:,4)*stim_temp_logic;
                temp_data = temp;
                %temp_cost = temp_data(:,4);
                
                % ---------- 3 chan overlap only ---------------
                % Get data where shared channel resps <= 3
                %temp_data = temp_data(temp_data(:,15)<=3,1:15);
                
                temp_SI = temp_data(:,4);
                temp_SI_3chan = temp_data(:,4);
                
                % Find the responses corresponding to the highest SI score
                if non_zero_resp==true
                    % enforce both sides having to respong for valid max SI
                    for i=1:length(temp_SI)
                        if temp_data(i,10)==0 || temp_data(i,11)==0 || temp_data(i,12)==0 || temp_data(i,13)==0
                            %temp_data(i,4)=0;
                            temp_SI(i)=0;
                            temp_SI_3chan(i)=0;
                        end
                    end
                end
                [~, temp_SI_max] = max(temp_SI);
                
                if select_3chan==true
                    for i=1:length(temp_SI_3chan)
                        if temp_data(i,15)>3
                            %temp_data(i,4)=0;
                            temp_SI_3chan(i)=0;
                        end
                    end
                end
                % In case temp_SI_chan doesn't have any values greater than
                % zero, or there aren't any selective channel responses
                [~, temp_SI_3chan_max] = max(temp_SI_3chan);
                if max(temp_SI_3chan) == 0
                    cuff_1_2_chan_resp = 0;
                    cuff_3_4_chan_resp = 0;
                else
                    cuff_1_2_chan_resp = temp_data(temp_SI_3chan_max,10);
                    cuff_3_4_chan_resp = temp_data(temp_SI_3chan_max,11);
                end
                % If multiple SI scores share the same maximum value, this finds the row index of the
                % highest SI value and saves the parameters from that run
                if sum(temp_SI==max(temp_SI)) > 1
                    logic = double(temp_SI == max(temp_SI)); %convert logic array of max values to double
                    select_chan_sum = (temp_data(:,12) + temp_data(:,13)).*logic; % get sum of selective channels
                    [~,maxcostidx] = max(select_chan_sum); %find new index of highest val where total channels responding is higher
                    maxcostval = max(temp_SI(maxcostidx)); % then get max cost val at new index
                else
                    [maxcostval, maxcostidx] = max(temp_SI);
                end
                
                % Store metadata pertaining to max SI value for shared PW
                si_list = [si_list, temp_data(maxcostidx,4)];
                pw_1_2_list = [pw_1_2_list, pw_1_2];
                pw_3_4_list = [pw_3_4_list, pw_3_4];
                stim_1_2_list = [stim_1_2_list, temp_data(maxcostidx,5)];
                stim_3_4_list = [stim_3_4_list, temp_data(maxcostidx,6)];
                shared_chans = temp_data(maxcostidx,15);
                
                
                j = j + 1;
                
                opt_params_list.(['PW_',num2str(pw_1_2*1000)]).best_cost(1,expmt-1) = maxcostval;
                opt_params_list.(['PW_',num2str(pw_1_2*1000)]).best_cost_shared_chans(1,expmt-1) = shared_chans;
                opt_params_list.(['PW_',num2str(pw_1_2*1000)]).best_cost_A_selective(1,expmt-1) = temp_data(maxcostidx,12);
                opt_params_list.(['PW_',num2str(pw_1_2*1000)]).best_cost_B_selective(1,expmt-1) = temp_data(maxcostidx,13);
                opt_params_list.(['PW_',num2str(pw_1_2*1000)]).best_cost_A_resp(1,expmt-1) = temp_data(maxcostidx,10);
                opt_params_list.(['PW_',num2str(pw_1_2*1000)]).best_cost_B_resp(1,expmt-1) = temp_data(maxcostidx,11);
                opt_params_list.(['PW_',num2str(pw_1_2*1000)]).best_A_resp_3chan_overlap(1,expmt-1) = cuff_1_2_chan_resp;
                opt_params_list.(['PW_',num2str(pw_1_2*1000)]).best_B_resp_3chan_overlap(1,expmt-1) = cuff_3_4_chan_resp;
                
                % ----------------------------------------------------------------------------
                          
                % Create 3D mesh plots of stim amplitudes between cuff pairs
                % 1-2 and 3-4, plot the corresponding optimized cost value
                if plot_mesh == true
                    figure(1)
                    mesh(xq,yq,vq); % Add mesh
                    grid on
                    hold on
                    plot3(temp_data(:,5), temp_data(:,6), temp_data(:,4), 'o');
                    hold off
                    xlabel('Cuff 1:2 stim (\muA)');
                    ylabel('Cuff 3:4 stim (\muA)');
                    zlabel('Cost');
                    title_text = [expmt_list{expmt,1}.cohort,' Cuff 1-2 PW ', num2str(pw_1_2*1000),'\us vs Cuff 3-4 PW ', num2str(pw_3_4*1000), '\us'];
                    title(title_text);
                    view(2)
                    drawnow
                    %saveas(gcf, [file_directory, '\', title_text, '.fig'], 'fig');
                end
                if plot_color == true
                    title_text = ['F',expmt_list{expmt,1}.cohort,' PW ', num2str(pw_1_2*1000),'\mus SI Score'];
                    x = temp_data(:,5);
                    y = temp_data(:,6);
                    %z = temp_data(:,4);
                    z = temp_SI;
                    SI_surface_plot(x, y, z, title_text);
                    a = input("Continue? (Enter): ");
                end
                a=1;
            end
        end
    end
    
    j = j - 1;
    %bar_data_selective((expmt-1),:) = [best_cost_A_selective(expmt-1),best_cost_B_selective(expmt-1)];
    %bar_data((expmt-1),:) = [best_cost_A_resp(expmt-1),best_cost_B_resp(expmt-1)];
    
    opt_params_list.(['animal_',num2str(expmt)]).selective.cohort = expmt_list{expmt,1}.cohort;
    opt_params_list.(['animal_',num2str(expmt)]).selective.SI = si_list';
    opt_params_list.(['animal_',num2str(expmt)]).selective.pw_1_2 = pw_1_2_list';
    opt_params_list.(['animal_',num2str(expmt)]).selective.pw_3_4 = pw_3_4_list';
    opt_params_list.(['animal_',num2str(expmt)]).selective.stim_1_2 = stim_1_2_list';
    opt_params_list.(['animal_',num2str(expmt)]).selective.stim_3_4 = stim_3_4_list';
    opt_params_list.(['animal_',num2str(expmt)]).selective.shared_chans = shared_chans';
    %opt_params_list.(['animal_',num2str(expmt)]).selective.Chan_resp_A = 
    %opt_params_list.(['animal_',num2str(expmt)]).selective.Chan_resp_A_selec = 
    %opt_params_list.(['animal_',num2str(expmt)]).selective.Chan_resp_B = 
    %opt_params_list.(['animal_',num2str(expmt)]).selective.Chan_resp_B_select = 
    
    
    
    %Save overall best position data to text file
    %fileID = fopen(fullfile(file_directory, 'Best Parameters.txt'), 'wt');
    %fprintf(fileID, 'Best Cost: %.4f\n', avg_bestcost/j);
    %fprintf(fileID, 'Best stim 1-2: %4.0fuA\n', stim_1_2_pos/j);
    %fprintf(fileID, 'Best stim 3-4: %4.0fuA\n', stim_3_4_pos/j);
    %fclose(fileID);
    fprintf('===================================================\n')
end

fprintf('\nPW100us best cost: \n')
fprintf('%0.1f ', opt_params_list.PW_100.best_cost)
fprintf('\nPW500us best cost: \n')
fprintf('%0.1f ', opt_params_list.PW_500.best_cost)
fprintf('\nPW1000us best cost: \n')
fprintf('%0.1f ', opt_params_list.PW_1000.best_cost)

if plot_3chan_overlap==true
    for i=1:6
        bar_data_pw100_3chan(i,:) = [opt_params_list.PW_100.best_A_resp_3chan_overlap(1,i), opt_params_list.PW_100.best_B_resp_3chan_overlap(1,i)];
        bar_data_pw500_3chan(i,:) = [opt_params_list.PW_500.best_A_resp_3chan_overlap(1,i), opt_params_list.PW_500.best_B_resp_3chan_overlap(1,i)];
        bar_data_pw1000_3chan(i,:) = [opt_params_list.PW_1000.best_A_resp_3chan_overlap(1,i), opt_params_list.PW_1000.best_B_resp_3chan_overlap(1,i)];
    end
    figure(33)
    bar_plot(bar_data_pw100_3chan, bar_animals, 'PW100us Total # Selective Channels with Maximized SI', true)
    figure(34)
    bar_plot(bar_data_pw500_3chan, bar_animals, 'PW500us Total # Selective Channels with Maximized SI', true)
    figure(35)
    bar_plot(bar_data_pw1000_3chan, bar_animals, 'PW1000us Total # Selective Channels with Maximized SI', true)
    
end

if plot_PW_bar==true
    for i=1:6
        pw100_bar_data = ([opt_params_list.PW_100.best_cost_A_resp; opt_params_list.PW_100.best_cost_B_resp])';
        pw500_bar_data = ([opt_params_list.PW_500.best_cost_A_resp; opt_params_list.PW_500.best_cost_B_resp])';
        pw1000_bar_data = ([opt_params_list.PW_1000.best_cost_A_resp; opt_params_list.PW_1000.best_cost_B_resp])';
    end
    figure(30)
    bar_plot(pw100_bar_data, bar_animals, 'PW100us Total # of Channels Responding with Maximized SI', true)
    figure(31)
    bar_plot(pw500_bar_data, bar_animals, 'PW500us Total # of Channels Responding with Maximized SI', true)
    figure(32)
    bar_plot(pw1000_bar_data, bar_animals, 'PW1000us Total # of Channels Responding with Maximized SI', true)
end

%Bar graph of best cost and selectivity per animal regardless of PW
if plot_all_bar == true
    
    figure(3)
    hB = bar(bar_data); %for percentage
    %hB = bar(bar_data);
    hAx = gca;
    set(gca, 'XTickLabel',bar_animals, 'XTick',1:numel(bar_animals))
    labels = {'Cuff 1-2','Cuff 3-4';'Cuff 1-2','Cuff 3-4';'Cuff 1-2','Cuff 3-4';...
        'Cuff 1-2','Cuff 3-4';'Cuff 1-2','Cuff 3-4';'Cuff 1-2','Cuff 3-4'};
    title('Total # Selectively Activated Channels Between Cuff Pairs');
    legend('Cuff 1-2','Cuff 3-4');
    %title('# Selectivity Plot between Cuff pairs');
    ylim([0 32])
    xlabel('Animals');
    ylabel('# of Channels Responding with Maximized SI');
end

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

function SI_surface_plot(x, y, z, fig_title)
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
shading interp
drawnow
hold on;

%Mark the largest SI score on the surface plot
[maxcostval, maxcostidx] = max(z);
plot3(x(maxcostidx), y(maxcostidx), maxcostval, 'or', 'MarkerSize', 5, 'LineWidth', 1.5);
%All data points
plot3(x, y, z, 'ob','MarkerSize', 1, 'LineWidth', 3);
fprintf("max position| stim 1-2: %d, stim 3-4: %d",x(maxcostidx), y(maxcostidx))
hold off;
ylabel('\muA');
xlabel('\muA');
title(fig_title);

end

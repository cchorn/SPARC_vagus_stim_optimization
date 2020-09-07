function selective_plotter(data, arg1, fig_title, plot_type, unit_type, arg2)

% generates shared channel figure with colored, heatmap, surface, or
% contour plot of shared channels

% Inputs
% =======================================================================
% data      : (Nx4) numerical array containing stim data and channel
%                   responses
% arg1      : (Nx1) parameter for scaling data (minimum threshold stim or charge)
% fig_title : (char) figure title
% plot_type : (char) plot type selector
% unit_type : (char) unit selector
% arg2      : (Nx1) additional scaling param


stim_data = data(:,[1,3]);
switch unit_type
    case 'charge'
        stim_data = (arg1*stim_data)/1000000;
        units = '(C)';
        x_label = ['Cuff 1-2 Stim', units];
        y_label = ['Cuff 3-4 Stim', units];
    case 'amp'
        grid_space = 40;
        units = '(\muA)';
        x_label = ['Cuff 1-2 Stim', units];
        y_label = ['Cuff 3-4 Stim', units];
    case 'threshold'
        % arg1 must be a list of min values to divide data by, since each
        % animal and pw has its own min thresh value
        stim_data = stim_data./arg2;
        grid_space = 0.1;
        units = 'Multiples of Threshold';
        x_label = ['Cuff 1-2 Threshold'];
        y_label = ['Cuff 3-4 Threshold'];

end
xmax = max(stim_data(:,1)); ymax = max(stim_data(:,2));
N_data = length(stim_data(:,1));

switch plot_type
    case 'color'
        cont_fig = figure(1);
        ax_h = gca;
        [xq,yq] = meshgrid(0:grid_space:max(stim_data(:,1)), 0:grid_space:max(stim_data(:,2))); %create grid array of stim values with 40uA increment
        vq = round(griddata(stim_data(:,1),stim_data(:,2),arg1,xq,yq)); % plot stim data
        [M,c] = contourf(ax_h, xq,yq,vq,'ShowText','on');
        if max(max(vq))<5
            c.LevelStep = 1;
        end
        cbh = colorbar('v');
        oldcmap = colormap(jet(2*N_data));
        colormap(flipud(oldcmap(1/3*N_data:5/3*N_data,:)));

    case 'mesh'
        figure(1)
        grid_space = 40;
        [xq,yq] = meshgrid(0:grid_space:max(stim_data(:,1)), 0:grid_space:max(stim_data(:,2))); %create grid array of stim values with 40uA increment
        vq = griddata(stim_data(:,1),stim_data(:,2),arg1,xq,yq); % plot stim data
        
        mesh(xq,yq,vq); % Add mesh
        grid on
        hold on
        plot3(stim_data(:,1), stim_data(:,2), arg1, 'o');
        hold off
        xlabel('Cuff 1:2 stim (\muA)');
        ylabel('Cuff 3:4 stim (\muA)');
        zlabel('Cost');
        title_text = [expmt_list{expmt,1}.cohort,' Cuff 1-2 PW ', num2str(pw_1_2*1000),'\us vs Cuff 3-4 PW ', num2str(pw_3_4*1000), '\us'];
        title(title_text);
        view(2)
        drawnow
    case 'surface'
        %imagesc(temp(:,1), temp(:,2), arg1);
        figure(2)
        %surf(xq,yq,vq); % Add surf
        x = stim_data(:,1);
        y = stim_data(:,2);
        z = arg1;
        dt = delaunayTriangulation(x,y) ;
        tri = dt.ConnectivityList ;
        xi = dt.Points(:,1) ;
        yi = dt.Points(:,2) ;
        F = scatteredInterpolant(x,y,z);
        zi = F(xi,yi) ;
        trisurf(tri,xi,yi,zi)
        view(2)
        shading interp
        drawnow
        hold on;
        %Biggest value
        plot3(stim_data(maxcostidx,5), stim_data(maxcostidx,6), maxcostval, 'or', 'MarkerSize', 5, 'LineWidth', 1.5);
        %All data points
        plot3(stim_data(:,5), stim_data(:,6), stim_data(:,4), 'ob','MarkerSize', 1, 'LineWidth', 3);
        fprintf("max position| stim 1-2: %d, stim 3-4: %d",stim_data(maxcostidx,5), stim_data(maxcostidx,6))
        hold off;
        ylabel('\muA');
        xlabel('\muA');
        title_text = [expmt_list{expmt,1}.cohort,' Cuff 1_2 PW ', num2str(pw_1_2*1000),'uA vs Cuff 3_4 PW ', num2str(pw_3_4*1000), 'uA'];
        title([expmt_list{expmt,1}.cohort,' Cuff 1-2 PW ', num2str(pw_1_2*1000),'\muA vs Cuff 3-4 PW ', num2str(pw_3_4*1000), '\muA']);
        
end

title(fig_title);
xlim([min(stim_data(:,1)),max(stim_data(:,1))])
ylim([min(stim_data(:,2)),max(stim_data(:,2))])
xlabel(x_label);
ylabel(y_label);
hold off;
drawnow;

end


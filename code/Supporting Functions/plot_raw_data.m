function plot_raw_data(varargin)

snipData = varargin{1};
STA_t_min = varargin{2};
STA_t_max = varargin{3};
chan = varargin{4};
search_window = varargin{5};

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

try
    figure(1);
    
    % create a list of numbers between two values, with a specified number of values in between
    x_time = linspace(STA_t_min,STA_t_max,size(snipData,2))*1000/fs;
    subplot(sub_size(1),sub_size(2), elec_chan_map(chan));
    hold on;
    plot(x_time, snipData') % Plot full data
    plot(x_time, mean(snipData),'-k'); %data average
    hold off;
    title(['STA Elec Chan  ',num2str(chan)]);
    axis tight
    if max(max(search_window))>70
        ylim([-(max(search_window)+20) (max(search_window)+20)]);
    else
        ylim([-70 70]);
    end
    drawnow
    set(gca,'ButtonDownFcn',@createnew_fig)
catch ME
    sprintf("%s",ME.message)
end

end

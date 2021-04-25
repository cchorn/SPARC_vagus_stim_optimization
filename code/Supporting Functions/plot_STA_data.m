function plot_STA_data(varargin)

snipData = varargin{1};
STA_t_min = varargin{2};
STA_t_max = varargin{3};
chan = varargin{4};
search_window = varargin{5};
thresh = varargin{6};

if ~isempty(varargin{7})
    peak_time = varargin{7};
    max_val = varargin{8};
end

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
    figure(2);
    
    % create a list of numbers between two values, with a specified number of values in between
    x_time = linspace(STA_t_min,STA_t_max,size(snipData,2))*1000/fs;
    subplot(sub_size(1),sub_size(2), elec_chan_map(chan));
    hold on;
    plot(x_time, mean(snipData), '-k', 'LineWidth', 0.7); % data with window mean
    thresh = linspace(thresh, thresh, size(mean(snipData),2));
    plot(x_time, thresh,'-r','LineWidth', 1);
    plot(x_time, -thresh,'-r','LineWidth', 1);
    plot(peak_time, max_val,'or', 'MarkerSize',8,'MarkerFaceColor','r');
    
    hold off;
    axis tight
    if max(max(search_window))>4
        ylim([-1 (max(search_window)+1)]);
    else
        ylim([-1 3]);
    end
    drawnow
    
    % Allow creation of new figure if clicked on
    set(gca,'ButtonDownFcn',@createnew_fig)
catch ME
    sprintf("%s", ME.message);
end

end


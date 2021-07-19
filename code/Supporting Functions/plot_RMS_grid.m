function plot_RMS_grid(varargin)

Vrms_list = varargin{1};
Vrms_time = varargin{2};
fig_title = varargin{3};
chan = varargin{4};
fs=30e3;
N_channels  = 32;
sub_size = numSubplots(N_channels);
plot_cv_regions = true;
no_axes = false; % Removes borders and minimizes empty space, no axes or numbers

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

figure(27)
subplot(sub_size(1),sub_size(2), elec_chan_map(chan));
plot(Vrms_time, Vrms_list)


if length(varargin)>4
    threshold = varargin{5};
    rms_offset = varargin{6};
    hold on
    thresh_line = linspace(threshold,threshold,length(Vrms_time));
    plot(Vrms_time,thresh_line,'-r','LineWidth', 1);
    %[max_val, max_idx] = max(Vrms_list(rms_offset:end-rms_offset));
    if length(varargin)>6
        bin_range = varargin{7};
    else
        bin_range = 50;
    end
    if length(varargin)>7 && varargin{8}==true
        % extra offset towards end of signal in case stim artifacts exist (needed for F21-19 Trial49)
        rms_end_offset = round(varargin{9}*10);
    else
        rms_end_offset = 0;
    end
    
    [max_val, max_idx] = find_peaks(Vrms_list(rms_offset:end-rms_offset-rms_end_offset),threshold, fs, bin_range);
    if ~isempty(max_idx)
        if max_val > threshold
            for i=1:length(max_val)
                if max_val(i)>3
                    max_val(i)=3;
                end
                plot((rms_offset/10+Vrms_time(max_idx(i))),max_val(i),'or')
            end
        end
    end
    hold off
end

if no_axes == true
    
    %drawing squares like this...
    x_start = floor((chan-0.1)/4)/8;
    %x_start = ceil(chan/4)-1;
    x_end = 1/8;
    if mod(chan,4)==1
        j=3;
    elseif mod(chan,4)==2
        j=2;
    elseif mod(chan,4)==3
        j=1;
    elseif mod(chan,4)==0
        j=0;
    end
    y_start = j/4;
    %y_start = 4 - mod(chan-1,4);
    y_end = 1/4;%5 - mod(chan,4);
    
    ax = gca;
    set(ax,'XTick',[], 'YTick', []);
    ax.Position = [x_start, y_start, x_end, y_end];
    %ax.Position = [0 0.75 0.125 0.25];
    drawnow;
else
    sgtitle(fig_title)
    
    ylim([0 3]);
    if ~mod(chan,4)==0
        set(gca,'XTick',[])
    else
        xlabel('time (msec)');
    end
    if chan>4
        set(gca,'YTick',[])
    else
        ylabel('\muV')
    end
end

end

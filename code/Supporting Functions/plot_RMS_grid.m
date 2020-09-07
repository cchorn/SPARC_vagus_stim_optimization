function plot_RMS_grid(varargin)

Vrms_list = varargin{1};
Vrms_time = varargin{2};
fig_title = varargin{3};
chan = varargin{4};
Fs=30e3;
N_channels  = 32;
sub_size = numSubplots(N_channels);

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
sgtitle(fig_title)
subplot(sub_size(1),sub_size(2), elec_chan_map(chan));
plot(Vrms_time, Vrms_list)

if length(varargin)>4
    thresh = varargin{5};
    rms_offset = varargin{6};
    hold on
    thresh_line = linspace(thresh,thresh,length(Vrms_time));
    plot(Vrms_time,thresh_line,'-r','LineWidth', 1);
    %[max_val, max_idx] = max(Vrms_list(rms_offset:end-rms_offset));
    if length(varargin)>6
        bin_range = varargin{7};
    else
        bin_range = 0.05;
    end
    [max_val, max_idx] = find_peaks(Vrms_list(rms_offset:end-rms_offset),thresh, Fs, bin_range);
    if ~isempty(max_idx)
        if max_val > thresh
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

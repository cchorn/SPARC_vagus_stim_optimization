function plot_RMS_figure(varargin)

Vrms_list = varargin{1};
Vrms_time = varargin{2};
fig_title = varargin{3};
Fs = 30e3;

figure(23);
plot(Vrms_time, Vrms_list)
xlabel('time (msec)');
ylabel('\muV')
title(fig_title);

if length(varargin)>3
    threshold = varargin{4};
    rms_offset = varargin{5};
    hold on
    thresh_line = linspace(threshold,threshold,length(Vrms_time));
    plot(Vrms_time,thresh_line,'-r','LineWidth', 1);
    if length(varargin)>7
        bin_range = varargin{8};
    else
        bin_range = 0.05;
    end
    [max_val, max_idx] = find_peaks(Vrms_list(rms_offset:end-rms_offset),threshold, Fs, bin_range);
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

axis tight
ylim([0 3]);
drawnow

% Plot CV regions
if length(varargin)>5
    resp_bins = varargin{6};
    time_bins = varargin{7};
    plot_CV_regions(Vrms_time, rms_offset, resp_bins, time_bins);
end

end


function plot_chan_grid(Vrms_list, Vrms_time, fig_title, chan, trial, threshold)

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
    titleName = sprintf('Trial %d RMS Window %d us, step %.3g us', trial, windowSize/0.03, stepSize/0.03);
    sgtitle(titleName)
    subplot(sub_size(1),sub_size(2), elec_chan_map(chan));
    plot(Vrms_time, Vrms_list)
    if ismember(chan,skip_channels)==0
        hold on
        plot(Vrms_time,threshold,'-r','LineWidth', 1);
        if max_idx~=0
            plot((rms_offset/10+Vrms_time(max_idx)),threshold(1),'or')
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

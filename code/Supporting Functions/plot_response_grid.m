function plot_response_grid(resp_grid_1, resp_grid_2, fig_title)

% Generates a 4x32 square grid with each channel represented as red (cuff
% 1-2 response), blue (cuff 3-4 response), and purple (response from both cuff pairs)

% resp_grid_1  : (1x32 num) array of channel responses from cuff 1-2
% resp_grid_2  : (1x32 num) array of channel responses from cuff 3-4
% fig_title    : (char) figure title


N_electrodes  = 32;

sub_size = numSubplots(N_electrodes);

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

for chan=1:N_electrodes
    figure(10)
    grid_fig = subplot(sub_size(1),sub_size(2), elec_chan_map(chan));
    sgtitle(fig_title);
    if resp_grid_1(chan)==1 && resp_grid_2(chan)==0
        rectangle('Position',[0,0,2,2],'FaceColor','r','EdgeColor','k','LineWidth',3)
    elseif resp_grid_1(chan)==0 && resp_grid_2(chan)==1
        rectangle('Position',[0,0,2,2],'FaceColor','b','EdgeColor','k','LineWidth',3)
    elseif resp_grid_1(chan)==1 && resp_grid_2(chan)==1
        rectangle('Position',[0,0,2,2],'FaceColor','.5 0 .5','EdgeColor','k','LineWidth',3)
    else
        rectangle('Position',[0,0,2,2],'FaceColor','w','EdgeColor','k','LineWidth',3)
    end
    set(gca,'XTick',[])
    set(gca,'YTick',[])

end
drawnow;
%a=input('next combo: ');
end

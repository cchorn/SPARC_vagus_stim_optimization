function plot_raw_grid(varargin)

data = varargin{1};
x_time = varargin{2};
fig_title = varargin{3};
chan = varargin{4};

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

figure(22)
sgtitle(fig_title)
subplot(sub_size(1),sub_size(2), elec_chan_map(chan));
plot(x_time, data)
axis tight
ylim([-250 250]);
drawnow

%ylim([0 3]);
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


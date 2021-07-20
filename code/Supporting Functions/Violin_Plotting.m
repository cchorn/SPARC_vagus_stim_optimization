%load PW100 
%load PW500 
%load PW1000

% From "Max SI params.xlsx" table, PW arrays are bellow
%PW100 = [NaN, NaN;
%         0.40,0.05;
%         NaN, NaN;
%         0.68,0.21;
%         1.21,0.16;
%         0.22,0.26]
%PW500 = [NaN,0.09;
%         0.53,0.09;
%         1.26,0.07;
%         0.41,0.15;
%         0.27,0.03;
%         0.57,0.00]
%PW1000 = [0.00,0.06;
%          0.57,0.07;
%          1.23,0.04;
%          NaN,0.08;
%          0.6,0.1;
%          0.56,0.19]

%set universal variables
ylimits = [-0.5,2];
ylab = 'Distance between Centroids (mm)';
xtick = [1, 2];
xticklab = {'Max SI', 'Max Stim'};
PW100xSI = [1,1,1,1];
PW100xStim = [2,2,2,2];
xSI = [1,1,1,1,1,1];
xStim = [2,2,2,2,2,2];


%% PW 100 Plotting
violin(PW100, 'mc', '', 'plotlegend', 0)
hold on
scatter(PW100xSI, PW100(:,1), 'filled', 'black');
scatter(PW100xStim, PW100(:,2), 'filled', 'black');
hold off
title('Distance between Centroids; PW = 100us')
ylim(ylimits)
ylabel(ylab)
xticks(xtick)
xticklabels(xticklab)

%% PW 500 Plotting
figure
violin(PW500, 'mc', '', 'plotlegend', 0)
hold on
scatter(xSI, PW500(:,1), 'filled', 'black');
scatter(xStim, PW500(:,2), 'filled', 'black');
hold off
title('Distance between Centroids; PW = 500us')
ylim(ylimits)
ylabel(ylab)
xticks(xtick)
xticklabels(xticklab)

%% PW 1000 Plotting
figure
violin(PW1000, 'mc', '', 'plotlegend', 0)
hold on
scatter(xSI, PW1000(:,1), 'filled', 'black');
scatter(xStim, PW1000(:,2), 'filled', 'black');
hold off
title('Distance between Centroids; PW = 1000us')
ylim(ylimits)
ylabel(ylab)
xticks(xtick)
xticklabels(xticklab)

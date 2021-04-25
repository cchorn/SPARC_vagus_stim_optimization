function min_thresh_plot(expmt_list)

rec_path = 'C:\Users\shulgac\Documents\MATLAB\Data';
close all;
save_fig2 = false;
x_offset = 20;
y_offset = 20;
rheo1_list = [];
rheo2_list = [];
chronax1_list = [];
chronax2_list = [];
plot_min_thresh = false;

for expmt=1:length(expmt_list) % for each day
    disp(expmt)
    exp_data = expmt_list{expmt};
    if plot_min_thresh==true
        fig(expmt) = figure;
    end
    file_directory = [rec_path, '\F', expmt_list{expmt,1}.cohort];
    
    
    if size(exp_data.minthresh,2) < 3
        for i=1:2 % Plot
            min_val(i) = min(exp_data.minthresh(:,i));
            min_pw(i) = exp_data.pulseWidth(i)*1000;
        end
        if plot_min_thresh==true
            plot(min_pw(2:end), min_val(2:end), 'o','MarkerSize',(2 + i*1.5));
            plottext=strcat('(',num2str(min_pw(2:end)'),',',num2str(min_val(2:end)'),')');
            text(min_pw(2:end)+x_offset, min_val(2:end)+y_offset, plottext,'horiz','left','vert','bottom')
            legend('Cuff 1-2');
        end
    else
        %Plot 1st 3 data points
        for i=1:3 % Plot for cuff 1-2
            min_val(i) = min(exp_data.minthresh(:,i));
            min_pw(i) = exp_data.pulseWidth(i)*1000;
        end
        
        % fit line conditions
        x = 0:1:3000;
        lb = [0];
        ub = [];
        modelfun = @(b,xdata) b(1) + b(2)./xdata;
        x0 = [500,200];
        if plot_min_thresh==true
            plot(min_pw, min_val, 'ro','MarkerSize',(3 + i*1.5));
            hold on;
        end
        coeff1 = lsqcurvefit(modelfun,x0,min_pw,min_val,lb,ub);
        %[coeff1,normresm] = fmincon(@(b) norm(min_val - modelfun(b,min_pw')), x0, [],[],[],[],[], 0,[]);
        %SSE = sum((min_val' - modelfun(coeff1,min_pw')).^2);
        yFit1 = coeff1(1)+coeff1(2)./x;
        if plot_min_thresh==true
            p1 = plot(x, yFit1, '-r');
            p1.LineWidth = 2.5;
        end
        % Find Rheobase and Chronaxie locations
        yFit1 = round(yFit1);
        rheo1 = round(coeff1(1)+coeff1(2)/1500);
        chronax1 = find(yFit1>=2*rheo1, 1, 'last');
        
        rheo1_list = [rheo1_list, rheo1];
        chronax1_list = [chronax1_list, chronax1];
        
        if plot_min_thresh==true
            plot(0:1500, linspace(rheo1,rheo1,size(0:1500,2)), 'color', [0 0.5 0], 'linestyle','--');
            plot(0:chronax1, linspace(yFit1(chronax1),yFit1(chronax1),size(0:chronax1,2)),'--k');
            plot(linspace(chronax1,chronax1,size(0:yFit1(chronax1),2)), 0:yFit1(chronax1),'--k');
            
            plottext=strcat('(',num2str(min_pw'),',',num2str(min_val'),')');
            text(min_pw+x_offset, min_val+y_offset, plottext,'horiz','left','vert','bottom')
        end
        % Plot last 3 data points, cuff 3-4
        for i=1:3
            min_val(i) = min(exp_data.minthresh(:,i+3));
            min_pw(i) = exp_data.pulseWidth(i+3)*1000;
        end
        
        if plot_min_thresh==true
            plot(min_pw, min_val, 'bo','MarkerSize',(3 + 3));
        end
        coeff2 = lsqcurvefit(modelfun,x0,min_pw,min_val,lb,ub);
        %[coeff2,normresm] = fmincon(@(b) norm(min_val - modelfun(b,min_pw')), beta0, [],[],[],[],[],0,[]);
        %SSE = sum((min_val' - modelfun(coeff1,min_pw')).^2);
        yFit2 = coeff2(1)+coeff2(2)./x;
        if plot_min_thresh==true
            p2 = plot(x, yFit2, '-b');
            p2.LineWidth = 1.5;
        end
        
        % Find Rheobase and Chronaxie locations
        yFit2 = round(yFit2);
        rheo2 = round(coeff2(1)+coeff2(2)/1500);
        chronax2 = find(yFit2>=2*rheo2, 1, 'last'); % Get the last value at 2xRheo value
        
        rheo2_list = [rheo2_list, rheo2];
        chronax2_list = [chronax2_list, chronax2];
        
        if plot_min_thresh==true
            plot(0:1500, linspace(rheo2,rheo2,size(0:1500,2)), 'color', [0 0.5 0], 'linestyle','--');
            plot(0:chronax2, linspace(yFit2(chronax2),yFit2(chronax2),size(0:chronax2,2)),'--k');
            plot(linspace(chronax2,chronax2,size(0:yFit2(chronax2),2)), 0:yFit2(chronax2),'--k');
            
            plottext=strcat('(',num2str(min_pw'),',',num2str(min_val'),')');
            text(min_pw+x_offset, min_val+y_offset, plottext,'horiz','left','vert','bottom')
            
            hold off;
            drawnow;
        end
        
        %eq1 = [num2str(round(coeff1(1))),' + ', num2str(round(coeff1(2))), 'e^{', num2str(double(coeff1(3))),'x}'];
        %eq2 = [num2str(round(coeff2(1))),' + ', num2str(round(coeff2(2))), 'e^{', num2str(double(coeff2(3))),'x}'];
        %legend('Cuff 1-2',eq1,'Cuff 3-4', eq2);
        
    end
    
    if plot_min_thresh==true
        title(strcat(['F',exp_data.cohort]));
        ylabel('\muA');
        xlabel('Pulse Width (\mus)');
        xlim([0 1500]);
        ylim([0 2000]);
        if max(min_val) > 2000
            ylim([0 (max(min_val) + 200)]);
        end
        title_text = [expmt_list{expmt}.cohort,' Minimum Threshold Plot'];
    end
    new_folder = 'Minimum Threshold Figures';
    %cd(['F',expmt_list{expmt,1}.cohort])
    %if ~exist(new_folder, 'dir')
    %    mkdir(new_folder)
    %end
    if save_fig2==true
        saveas(gcf, [file_directory, '\', new_folder, '\', title_text, '.fig'], 'fig');
        saveas(gcf, [file_directory, '\', new_folder, '\', title_text, '.png'], 'png');
        
    end
    %cd ..
end

rheo_list = floor(([rheo1_list, rheo2_list])./10)*10;
chronax_list = floor(([chronax1_list, chronax2_list])./10)*10;

figure(8)
hist(rheo_list, 8)
title('Rheobase Histogram across animals')
xlim([0 800])
xlabel('Pulse Width (\muS)')
ylabel('\muA')

figure(9)
hist(chronax_list, 16)
title('Chronaxie Histogram across animals')
xlim([0 800])
xlabel('Pulse Width (\muS)')
ylabel('Count')

a=1;

%close all;

end

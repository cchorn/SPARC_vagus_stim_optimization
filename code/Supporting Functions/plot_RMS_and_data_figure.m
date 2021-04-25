function plot_RMS_and_data_figure();

    figure(6)
%                     sgtitle(sprintf('Cohort %s | Trial %d | Channel %d |          ',strrep(cohort,'_','-'),trial,chan))
%                     subplot(1,2,1)
%                     plot(x_time, STA_data, '-k','LineWidth',0.1);
%                     title('Filtered signal')
%                     ylim([-20,20]);
%                     xlim([STA_t_min,STA_t_max]);
%                     drawnow
%                     titleName = sprintf('RMS Window %d us, step %.3g us', windowSize/0.03, stepSize/0.03);
%                     subplot(1,2,2)
%                     plot(Vrms_time, Vrms_list)
%                     hold on
%                     plot(Vrms_time,x_threshold,'-r','LineWidth', 1);
%                     if max_idx~=0
%                         plot((rms_offset+max_idx)/10,x_threshold(1),'or')
%                     end
%                     hold off
%                     title(titleName)
%                     axis tight
%                     ylim([-1,5])
%                     drawnow
%                     
%                     if plot_cv_shades==true
%                         %C fibers
%                         hold on;
%                         t_start = time_bins(4);
%                         t_end = Vrms_time(end-rms_offset);
%                         x_points = [t_start,t_start,t_end,t_end];
%                         y_points = [0, 3,3,0];
%                         c_region = fill(x_points,y_points,[1, 0, 0]);
%                         c_region.FaceAlpha = 0.05;
%                         
%                         %a-delta
%                         t_start = time_bins(60);
%                         if t_start < rms_offset/10
%                             t_start = rms_offset/10;
%                         end
%                         t_end = time_bins(4);
%                         x_points = [t_start,t_start,t_end,t_end];
%                         y_points = [0, 3,3,0];
%                         a_delta_region = fill(x_points,y_points,[0, 0, 1]);
%                         a_delta_region.FaceAlpha = 0.05;
%                         
%                         % shade all the time bin regions with a stim response
%                         for k=1:length(resp_bins)-1
%                             if resp_bins(k)==1 && time_bins(k)>rms_offset/10
%                                 t_start = time_bins(k);
%                                 t_end = time_bins(k+1);
%                                 x_points = [t_start,t_start,t_end,t_end];
%                                 y_points = [0, 3,3,0];
%                                 a_delta_region = fill(x_points,y_points,[0, 0.5, 0.5]);
%                                 a_delta_region.FaceAlpha = 0.08;
%                             elseif resp_bins(k)==0 && time_bins(k)>rms_offset/10
%                                 t_start = time_bins(k);
%                                 t_end = time_bins(k+1);
%                                 x_points = [t_start,t_start,t_end,t_end];
%                                 y_points = [0, 3,3,0];
%                                 a_delta_region = fill(x_points,y_points,[0.9, 0.4, 0.4]);
%                                 a_delta_region.FaceAlpha = 0.08;
%                             end
%                         end
%                         hold off;

end
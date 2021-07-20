function plot_CV_regions(Vrms_time, rms_offset, resp_bins, time_bins)

%C fibers
hold on;
t_start = floor(10*time_bins(4))/10;
t_end = floor(10*Vrms_time(end-rms_offset))/10;
t_start_idx = find(floor(10*Vrms_time)/10 == t_start,1); % get index of time bin mark WRT Vrms time
t_end_idx = find(floor(10*Vrms_time)/10 == t_end,1);
x_points = [t_start,t_start,t_end,t_end];
y_points = [0, 3,3,0];
c_region = fill(x_points,y_points,[1, 0, 0]);
c_region.FaceAlpha = 0.05;

%a-delta
t_start = floor(10*time_bins(end))/10;
if t_start < floor(10*Vrms_time(rms_offset))/10
    t_start = floor(10*Vrms_time(rms_offset))/10;
end
t_end = floor(10*time_bins(4))/10;
x_points = [t_start, t_start, t_end, t_end];
y_points = [0, 3,3,0];
a_delta_region = fill(x_points,y_points,[0, 0, 1]);
a_delta_region.FaceAlpha = 0.05;

% shade all the time bin regions with a stim response
for k=1:length(resp_bins)-1
    t_start = time_bins(k);
    t_end = time_bins(k+1);
    x_points = [t_start,t_start,t_end,t_end];
    y_points = [0, 3,3,0];
    if resp_bins(k)==1% && time_bins(k)>rms_offset/10
        a_delta_region = fill(x_points,y_points,[0, 0.5, 0.5]);
    elseif resp_bins(k)==0% && time_bins(k)>rms_offset/10
        a_delta_region = fill(x_points,y_points,[0.9, 0.4, 0.4]);
    end
    a_delta_region.FaceAlpha = 0.08;
end
hold off;

end


function plot_response_distribution(expmt_list, N_expmt, opt_params)


for anim_expmt=N_expmt
    
    expmt = anim_expmt - 1;
    bar_animals{1,(expmt)} = expmt_list{anim_expmt,1}.cohort;

end

pw_list=[100, 500, 1000];
for pw=1:length(pw_list) % dataset collection for each animal
    % Threshold data
    param_data = filter_data(opt_params.(['PW_',num2str(pw_list(pw))]).threshold.max_SI);
    min_data.(['PW_',num2str(pw_list(pw))]).resp_data          = ([param_data.A_selective; param_data.B_selective])';
    min_data.(['PW_',num2str(pw_list(pw))]).overlap_data  = ([param_data.shared_chans])';
    min_data.(['PW_',num2str(pw_list(pw))]).stim_data     = [param_data.stim_data_A; param_data.stim_data_B];
        
    % max stim data
    param_data = filter_data(opt_params.(['PW_',num2str(pw_list(pw))]).max_stim.max_SI);
    max_data.(['PW_',num2str(pw_list(pw))]).resp_data          = ([param_data.A_selective; param_data.B_selective])';
    max_data.(['PW_',num2str(pw_list(pw))]).overlap_data  = ([param_data.shared_chans])';
    max_data.(['PW_',num2str(pw_list(pw))]).stim_data     = [param_data.stim_data_A; param_data.stim_data_B];
    
end

% histogram 
figure(47);
overlap_histogram(min_data.PW_100, max_data.PW_100, bar_animals,'PW 100 # Responding Channels');
figure(48);
overlap_histogram(min_data.PW_500, max_data.PW_500, bar_animals,'PW 500 # Responding Channels');
figure(49);
overlap_histogram(min_data.PW_1000, max_data.PW_1000, bar_animals,'PW 1000 # Responding Channels');

end

function overlap_histogram(min_data, max_data, bar_animals, figure_title)

global selective_chan_lim

min_temp = min_data.resp_data + min_data.overlap_data;
min_bin = reshape(min_temp, 1,[]);
min_bin = min_bin(~min_bin==0);

max_temp = max_data.resp_data + max_data.overlap_data;
max_bin = reshape(max_temp, 1,[]);
max_bin = max_bin(~max_bin==0);

%h1 = histogram(min_bin);
h1 = histfit(min_bin,length(max_bin),'poisson');
h1(1).FaceAlpha = 0.5;
h1(2).Color = [0 128/255 1];
h1(2).Visible = 'off';
xlim([0 32]);
ylim([0 10]);
hold on;

%h2 = histogram(max_bin);
h2 = histfit(max_bin,length(max_bin),'poisson');
h2(1).FaceAlpha = 0.5;
h2(2).Color = [1 191/255 0];
h2(2).Visible = 'off';
xlim([0 32]);
ylim([0 10]);
hold off;

%h1.Normalization = 'probability';
%h1.BinWidth = 0.25;
%h2.Normalization = 'probability';
%h2.BinWidth = 0.25;


legend([h1(1), h2(1)],'Threshold','Max Stim');%'More than 3 Shared Channels')

title(figure_title);
ylabel('# Cuff Channels');
xlabel('# responding MEA Channels');

end

function opt_params = filter_data(opt_params)

for i=1:size(opt_params,2)
    if isempty(opt_params(i).SI)
        opt_params(i).SI=0;
    end
    if isempty(opt_params(i).shared_chans)
        opt_params(i).shared_chans=0;
    end
    if isempty(opt_params(i).A_selective)
        opt_params(i).A_selective=0;
    end
    if isempty(opt_params(i).B_selective)
        opt_params(i).B_selective=0;
    end
    if isempty(opt_params(i).A_resp)
        opt_params(i).A_resp=0;
    end
    if isempty(opt_params(i).B_resp)
        opt_params(i).B_resp=0;
    end
    if isempty(opt_params(i).stim_data_A)
        opt_params(i).stim_data_A=0;
    end
    if isempty(opt_params(i).stim_data_B)
        opt_params(i).stim_data_B=0;
    end
end

end

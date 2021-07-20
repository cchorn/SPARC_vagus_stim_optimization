function analogData = getAnalogData(expmt_list, expmt, rec_filePath, N_channels, use_cloud_path, use_trellis_file)

% Grab analogData
trial = str2double(rec_filePath(end-5:end-4)); % Usually the datafile number is the trial number
session = find_session(expmt_list, expmt, trial);
entity = expmt_list{expmt}.entity(session,:);  % necessary for reading trellis data files, specified as "raw", or "hi-res"

cohort = expmt_list{expmt}.cohort;
cohort_path = strrep(expmt_list{expmt}.cohort,'-','_');

if use_cloud_path == true
    if use_trellis_file == true
        % Neuroshare function to read all channels on continuous recording file. rec_filePath contains directory
        % to file starting with "C:/...", entity type for reading data based om sampling rate (ex: 'hi-res' = 2kHz,
        % 'raw' = 30Kz), elec_list is number array (1:32)
        analogData = read_continuousData(rec_filePath, char(entity), N_channels);
    else
        
        username = char(java.lang.System.getProperty('user.name'));
        box_path    = ['C:\Users\',username,'\Box\Data_for_Dylan'];
        mat_file_name = sprintf('F%s_trial_%d_data',cohort_path,trial);
        matfilePath = fullfile([box_path,'\',cohort,'\.mat data files\', mat_file_name,'.mat']);
        temp = load(matfilePath);
        analogData = temp.(mat_file_name);
    end
    
else
    % Neuroshare function to read all channels on continuous recording file. rec_filePath contains directory
    % to file starting with "C:/...", entity type for reading data based om sampling rate (ex: 'hi-res' = 2kHz,
    % 'raw' = 30Kz), elec_list is number array (1:32)
    analogData = read_continuousData(rec_filePath, char(entity), N_channels);
end

end


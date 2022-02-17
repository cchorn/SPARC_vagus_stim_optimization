function [json_obj] = convert_exp_to_json_metadata(expmt_list, save_path)
% Function that parses the struct containing ferret stim data into JSON
% readable format for SPARC data server

% Syntax:
%    json_data = convert_exp_to_json_metadata(expmt_list)
%
% Inputs:
%    expmt_list -- data struct containing information about animal runs
%    save_path -- boolean to enable saving animal data as json files for
%                 each animal
% Outputs:
%    json_data -- cell array of json data struct variables for each animal
%

%       "absolute_time_started": "2018-01-23T16:35:10,000000-05:00",
%     "subject_id": "14-18",
%     "paddles_placement": [
%       {
%         "paddle_id": 1,
%         "location": "S4",
%         "location_full_name": "Stomach Segment 4"
%       },
%     "body_mass_unit": "kg",
%     "recording_channels": [
%         {
%         "type_full_name": "Gastric Myoelectric Activity",
%         "channel_id": 1,
%         "type": "GMA",
%         "source_channel_id": 161,
%         "paddle_id": 1
%       },


json_obj ={};
for i=1:size(expmt_list,1)
    exp_data = expmt_list{i,1};
    
    s = struct();
    file_id = [];
    
    % append the actual trial numbers to the list of usable run names
    for session=1:size(exp_data.trial_list,1)
        
        session_start = exp_data.trial_list(session,1);
        session_end = exp_data.trial_list(session,2);
        
        for trial=session_start:1:session_end
            file_id = [file_id, trial];
        end
    end
    
    for j=1:1:length(file_id)
        
        % Convert file_ids to field names starting with 'run'
        run_id = horzcat('run_',sprintf('%04d',file_id(j)));
        s.runs.(run_id).raw_data_types = ["uint32","int16"];
        year = exp_data.date_dash(5:8);
        month = exp_data.date_dash(1:2);
        day = exp_data.date_dash(3:4);
        s.runs.(run_id).absolute_time_started = horzcat(year,'-',month,'-',day);
        s.runs.(run_id).run_id = file_id(j);
        s.runs.(run_id).source_data_types = ["matlab double"];
        s.runs.(run_id).recording_frequency =  30000;
        
        [stim_amp, pw, chan, session] = get_trial_data(exp_data, file_id(j));
        s.runs.(run_id).params.poststim_min = 0;
        s.runs.(run_id).params.amplitude = stim_amp;
        s.runs.(run_id).params.chan = chan;
        s.runs.(run_id).params.frequency = "2";
        s.runs.(run_id).params.prestim_min = 0;
        s.runs.(run_id).params.stim_min = 0;
        s.runs.(run_id).params.PW = pw;
        s.runs.(run_id).raw_data_file_id = file_id(j);
        s.runs.(run_id).recording_frequency_unit = "Hz";
        s.runs.(run_id).trial_type = "stim";
        
        % Build metadata 
        %s.metadata.organism_RRID = "FILL IN (ex: RRID:NCBITaxon_9669)"; % not using RRID
        s.metadata.sex = "M";
        s.metadata.strain = "n/a";
        s.metadata.body_mass = "FILL IN";
        s.metadata.age_range_max =  "FILL IN NUMBER";
        s.metadata.species = "Mustela putorius furo";
        s.metadata.group = "control";
        s.metadata.genotype = "n/a";
        s.metadata.subject_id = string(exp_data.cohort);
        s.metadata.body_mass_unit = "kg";
        s.metadata.experimental_log_file = "n/a";
        s.metadata.protocols.protocol_title = "SPARC - Acute surgery and experimentation of the gastrointestinal tract and vagus nerve in the ferret";
        s.metadata.protocols.protocol_url_or_doi = "https://www.protocols.io/view/sparc-acute-surgery-and-experimentation-of-the-gas-6a7hahn";
        s.metadata.age_category = "prime adult stage";
        s.metadata.number_of_runs = length(file_id);
        s.metadata.recording_channels = linspace(1,32,32); % same channel names
        s.metadata.cuff_placement.distance_to_nodose = "FILL IN";
        s.metadata.cuff_placement.distance_to_nodose_unit = "mm";
        s.metadata.cuff_placement.location = [1, 2, 3, 4];
        s.metadata.cuff_placement.cuff_contact_id = [1, 2, 3, 4];
        s.metadata.age_range_min = "FILL IN";
        s.metadata.age = "FILL IN";
        s.metadata.age_unit = "days";
        s.metadata.experiment_type = "acute";
        s.metadata.software.software_version = "1.8.3.294";
        s.metadata.software.software_URL = "https://rippleneuro.com/support/software-downloads-updates/";
        s.metadata.software.software_vendor = "Ripple Neuro";
        s.metadata.software.software_RRID = "";
        s.metadata.software.software = "Trellis";
        
    end
    
    json_obj{i,1} = s;
    
end

if nargin>1
    if islogical(save_path)
        if save_path == true
            folder_path = uigetdir('C:\','Select path to save json files');
            mkdir(fullfile(folder_path,'Animal JSON Data'));

            % Create new JSON file for each animal
            for k=1:1:size(json_obj,1)
                animal_id = expmt_list{k,1}.cohort;
                mkdir(fullfile(folder_path,'Animal JSON Data',animal_id));
                file = fullfile(folder_path, 'Animal JSON Data', animal_id, 'metadata.json');
                fid = fopen(file,'w');
                encodedJSON = jsonencode(json_obj{k,1});
                %encodedJSON = jsonencode(s,'PrettyPrint',true);
                fprintf(fid, encodedJSON); 
                fclose('all'); 
            end
        end
    end
end

end


function [stim_amp, pw, chan, session] = get_trial_data(exp_data, trial)
% Helper function to collect stimulation amplitude, pulse width, and
% channels

session  = [];
session_trial = [];
chan = "";
for i=1:size(exp_data.trial_list,1)
    t = exp_data.trial_list(i,:);
    if trial>=t(1) && trial<=t(2)
        trial_list = t(1):t(2);
        session_trial = find(trial_list==trial);
        session = i;
        chan = string(strrep(exp_data.cuff_list(session,1:end),':','-'));
        break
    end
end
pw = string(exp_data.pulseWidth(session));
stim_amp = exp_data.stim_hist(session, session_trial)/1000;

end

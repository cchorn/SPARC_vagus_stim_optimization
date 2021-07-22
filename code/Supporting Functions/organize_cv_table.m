function cv_table = organize_cv_table(expmt_list, data)
% Calculate and organize data output from "get_conduction_velocity" into a
% organized struct containing the cv speeds for each channel and every response
%
% Inputs
% =============================
% expmt_list: (struct) header file with animal metadata
% cv_min_data: (struct) data struct containing cv metadata from get_conduction_velocity
%
% Outputs
% =============================
% cv_table: (struct) reorganized cv data

header_info = data(1,:);
data = data(2:end,:);

animal_list = unique(data(:,1));
for i=1:size(animal_list,1) % each animal
    animal_name = animal_list{i}(2:end);
    idx = find(contains(data(:,1), animal_name)); % matching rows with animal string
    animal_data = data(idx,:);
    
    % access exp data fron index of animal name
    animal_names = cellfun(@(x) x.cohort, expmt_list, 'UniformOutput',0);
    animal_idx = strcmp(animal_names, animal_name); % iError occuring here! Find out why F25-19 matches with F21-19 string????
    exp_data = expmt_list{animal_idx};
    
    PW=unique([animal_data{:,3}]);
    for j=1:length(PW) % each PW
        animal_pw_data = animal_data([animal_data{:,3}]==PW(j),:); % get rows of data where PW matches
        
        cuff_pair = unique(animal_pw_data(:,2));
        for k=1:size(cuff_pair) % each cuff pair
            idx = find(contains(animal_pw_data(:,2), cuff_pair{k}));
            animal_pw_cuff_data = animal_pw_data(idx,:);
            
            % Separate all channels with responses and ignore empty cells
            chan = find(~cellfun(@isempty,animal_pw_cuff_data{1,10}));
            
            for m=1:size(chan)
                % Collect time indices where responses occured for that channel
                temp_data = animal_pw_cuff_data{:,10};
                resp_idx = temp_data{chan(m)};
                t = resp_idx*data{1,8}';
                
                % Need to get the right cuff pair name for field
                cuff_pair_field = ['cuff_', strrep(cuff_pair{k},':','_')];
                
                % animal name
                temp_name = exp_data.cohort;
                animal_field_name = ['F', strrep(temp_name, '-', '_')];
                
                % and pw name
                pw_field = ['pw_', num2str(1000*PW(j)),'us'];
                
                % Calculate and store cv data for set channel
                L = exp_data.nerve_length;
                acc = 0.5; % cv step accuracy adjustment for rounding values
                cv_table.(animal_field_name).(pw_field).(cuff_pair_field).(['Channel_', num2str(chan(m))]) = round(1000*(L/t)*(1/acc))*acc; % need to convert time from ms to sec
            end
        end
    end
end

end


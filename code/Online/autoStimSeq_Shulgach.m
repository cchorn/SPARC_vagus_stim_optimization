function autoStimSeq_Shulgach(cathElec, anodElec, stimAmp, preQuietT_s, stimT_s, postQuietT_s, freq, phaseDur_ms, elec_list, stim_port, rec_chan, entity, STA_bin, thresh_bin, enable_stim)
%   Neural Rehabilitation Engineering Lab
%   Jonathan Shulgach
%   Last Updated: 5/2/19
%
%   DESCRIPTION
%   ===================================================================
%   Automated optimization bisectional algorithm that:
%   1) performs stimulation trials using Trellis API
%   2) performs STA or spike threshold analysis on recording data
%   3) Adjust stimulation parameter until difference of 60uA
%
%   INPUTS
%   ===================================================================
%   cathElec        :   (int) cathode electrode 1 to 8
%   anodElec        :   (int) anode electrode 1 to 8
%   cathAmp         :   (1x2 numeric) stimulation values for stim trials
%                           where first value is current stim, second is
%                           previous stim value to define search window
%   preQuietT_s     :   (numeric) pre stim duration in seconds
%   stimT_s         :   (numeric) stim duration in seconds
%   postQuietT_s    :   (numeric) post stim duration in seconds
%   freq            :   (int) stimulatio nfrequency
%   phaseDur_ms     :   (int) phase duration of stim pulse
%                             (see assumption 4)
%   elec_list       :   (numeric) list of electrode channels to record
%                           [1:32]
%   stim_port       :   (char) grapevine port that stimulation headstage is
%                           plugged into
%   stim_chan       :   (int) channel
%   entity          :   (char vector) type of recording data to load
%   STA_bin         :   (1x2 numeric) time window for STA analysis
%   thresh_bin      :   (1x2 numeric) time window for threshold detection
%   enable_stim     :   (log) enabele or skip stim event
%
%   OUTPUTS
%   ===================================================================
%
%   EXAMPLE
%   ===================================================================
%   Ex 1: (using stim channels 1-2, 3mA 60sec stim at 2Hz, 1ms pw, electrodes
%           1-32, grapevine stim port B, rec port 1, raw analog data, STA and 
%           spike threshold time window 2ms prior, 400ms post, stim enabled)
%   autoStimSeq_Shulgach(1, 2, [3000], 0, 60, 0, 2, 1, [1:32], 'b', 1, 1, 'raw', [2 400], [2 400], true)
%
%   Ex 2: (using stim channels 3-4, 760uA stim bound with 100uA, 60sec stim at 2Hz, 500us pw, electrodes
%           1-32, grapevine stim port B, rec port 1, raw analog data, STA and 
%           spike threshold time window 2ms prior, 400ms post, stim disabled)
%   autoStimSeq_Shulgach(3, 4, [760 100], 0, 60, 0, 2, 0.5, [1:32], 'b', 1, 'raw', [2 400], [2 400], false)
%
%   NOTES
%   ===================================================================
%   1) Figure can be saved to data directory by uncommenting save section 
%


if ~isdir('C:\GIT')
    addpath(genpath('C:\Users\Jonathan\Desktop\Jonathan\CMU\RNEL'));
    addpath(genpath('C:\GIT'));
end
%Change to corresponding computer (use if on lab machine)
%code_dir = 'C:\Users\shulgac\Documents\MATLAB\Experiment Code';
code_dir = 'C:\Users\nanivadekarac\Documents\ferret_repo\Scripts_and_Programming_code\matlab_toolbox\experiment_design\Experiment Code';
if ~isdir(code_dir)
    addpath(code_dir);
end

% Configure the stim channel according to the cuff electrode channels
if cathElec==3
    stim_chan = 9;
else
    stim_chan = 1;
end

%Check for lower bound
if length(stimAmp)==1
    prev_cathAmp = 0;
else
    prev_cathAmp = stimAmp(2);
end
cathAmp = stimAmp(1);

%Recording channel configuration (In case wire is dual port)
if rec_chan==2
    elec_list = elec_list + 32;
end

%Entity channel configuration (future update: check entity type exists)
if strcmpi(entity,'raw')
    %To-do
elseif strcmpi(entity,'analog')
    %To-do
elseif strcmpi(entity,'hi-res')
    % To-do
end

first_file = 1;
trial_count = 1;
threshold = 30;
stim_diff = 30; % uA

%--------------------------------------------------------------------------
% Perform quiet trial for baseline data - TO-DO for next version
%[base_file, base_path] = uigetfile('*.nev','Please select baseline data');
%if isequal(base_file,0)
%    disp('User selected Cancel');
base_filepath = '';
%return
%else
%base_path = 'C:\Users\Jonathan\Desktop\Jonathan\CMU\RNEL\Data';
%base_path = 'C:\Users\shulgac\Documents\MATLAB\Data'; % use for testing
%base_file = 'datafile0009.nev';% use for testing
%base_filepath = fullfile(base_path,base_file);
%end
%--------------------------------------------------------------------------

fprintf("============================================================\n")
fprintf("Experiment Begin: %s, at %s\n", date, datestr(now,'HH:MM:SS'))
fprintf("Cuff Channels: %d-%d | Starting Stim: %duA | Frequency: %dHz | Phase duration: %0.3fms\n",cathElec, anodElec, stimAmp(1), freq, phaseDur_ms)
fprintf("Electrodes: %s\n", num2str(elec_list))
fprintf("Stim port: %s | Stim channel: %d | recording channel: %d | entity: %s\n", stim_port, stim_chan, rec_chan, entity)
fprintf("STA time bin: %dms to %dms | Threshold time bin: %dms to %dms\n",-STA_bin(1), STA_bin(2), -thresh_bin(1), thresh_bin(2))
LogicalStr = {'false', 'true'};
fprintf("Stimulation: %s\n", LogicalStr{enable_stim + 1})
text = 'Continue (y or n)?: ';
prompt = input(text, 's');
if strcmpi(prompt,'y')
    % Begin bisection algorithm
    while abs((sum(prev_cathAmp) - sum(cathAmp))) > stim_diff
        cathAmp = round(cathAmp/20)*20; % Force input to be multiple of 20
        if mod(cathAmp,2)==1 
            fprintf("Error with stimulus, not increment of 20uA. Recommend restarting at \n", round(cathAmp/20)*20)
            break;
        end
        fprintf("------------------------------------------------------------\n")
        fprintf("Performing stimulation trial %d | Trial stim: %duA\n", trial_count, cathAmp)
        cathAmp_step = getStimSteps(cathAmp); % Get new 1x4 step stim values
        if enable_stim
            stimseqDesign(cathElec, anodElec, cathAmp_step, preQuietT_s, stimT_s, postQuietT_s, freq, phaseDur_ms);
        end
        fprintf("Stimulation trial done\n")
        fprintf("------------------------------------------------------------\n")
        if first_file==1        %Only ask for filepath once, most recent file
            [rec_file, rec_path] = uigetfile('*.nev','Please select stim recording data');
            if isequal(rec_file,0)
                disp('User selected Cancel');
                break;
            else
                %    rec_path = 'C:\Users\Jonathan\Desktop\Jonathan\CMU\RNEL\Data';
                %    rec_path = 'C:\DataTanks\2019\20-19\Trellis'; % use for testing
                %    rec_file = 'datafile0009.nev';% use for testing
                rec_filepath = fullfile(rec_path,rec_file);
                first_file = 0;
            end
        else
            new_file_num = str2num(rec_file(end-7:end-4))+1;
            rec_file = [rec_file(1:8),num2str(new_file_num,'%04d'),'.nev'];
            rec_filepath = fullfile(rec_path, rec_file);
        end
        if exist(rec_filepath, 'file') ==0
            error('could not find file %s', rec_file)
            break;
        else
            disp(rec_file)
            [stim_decision, AP_list, spike_figure] = spike_threshold_detection_Shulgach(rec_filepath, base_filepath, elec_list, threshold, stim_port, stim_chan, rec_chan, entity, STA_bin, thresh_bin);
            fprintf("----------------------------------------------------\n")
            fprintf("Total spikes detected: %d\n", sum(AP_list))
            if sum(AP_list) >= 30
                fprintf("WARNING: Abnormally large amount of spikes detected....Recommend Restart\n")
            end
            if stim_decision==1 % spike seen, lower stim
                new_cathAmp = round(round(cathAmp - abs(prev_cathAmp - cathAmp)/2)/20)*20;
                fprintf("Decreasing Stimulation - New Stimulation: %duA", round(new_cathAmp/20)*20);
            else
                new_cathAmp = round(round(cathAmp + abs(prev_cathAmp - cathAmp)/2)/20)*20;
                fprintf("Increasing Stimulation - New Stimulation: %duA", round(new_cathAmp/20)*20);
            end
        end
        trial_count = trial_count + 1;
        fprintf(" | Stim difference: %duA\n", abs(cathAmp - round(new_cathAmp/20)*20))
        
        %Save figure
        %if ~isdir(fullfile(rec_path,'Figures'))
        %    folder_dir = fullfile(rec_path,'Figures');
        %    mkdir(folder_dir);
        %else
        %    folder_dir = [rec_path,'\Figures'];
        %end
        %current_dir = pwd;
        %cd(folder_dir);
        %saveas(spike_figure, ['STA Grid Plot ', rec_file(end-7:end-4)],'png');
        %cd(current_dir);
        
        beep
        % Manual Check (can disable)
        text = 'Continue (y or n), manual override (m), or restart (r)?: ';
        prompt = input(text, 's');
        if strcmpi(prompt,'n')
            break;
        elseif strcmpi(prompt,'r')
            %restart loop with same stim params
        elseif strcmpi(prompt,'y')
            prev_cathAmp = cathAmp;
            cathAmp = round(new_cathAmp/20)*20;
        elseif strcmpi(prompt,'m')
            text = 'higher or lower (h or l)?: ';
            manual_prompt = input(text, 's');
            if strcmpi(manual_prompt,'h')
                new_cathAmp = round(round(cathAmp + abs(prev_cathAmp - cathAmp)/2)/20)*20;
                fprintf("Increasing Stimulation - New Stimulation: %duA", new_cathAmp);
            elseif strcmpi(manual_prompt,'l')
                new_cathAmp = round(round(cathAmp - abs(prev_cathAmp - cathAmp)/2)/20)*20;
                fprintf("Decreasing Stimulation - New Stimulation: %duA", new_cathAmp);
            else
                fprintf("Unknown input, continuing with recommended stimulus")
            end
            prev_cathAmp = cathAmp;
            cathAmp = round(new_cathAmp/20)*20;
            fprintf(" | Stim difference: %duA\n", abs(prev_cathAmp - cathAmp))
        end
    end
end
fprintf("============================================================\n")
fprintf("Final stimulatim: %duA\n", round(prev_cathAmp/20)*20)
fprintf("done\n")


function stimTimes = getStimEvents(expmt_list, expmt, rec_filePath)
% The sole purpose of this function is to obtain the stimevents from a file


trial = str2double(rec_filePath(end-5:end-4)); % Usually the datafile number is the trial number
session = find_session(expmt_list, expmt, trial);
stim_port   = expmt_list{expmt}.stim_port;


% Channel 1 on headstage was used to stim nerve cuff electrode contacts 1-2, channel 9 on headstage for contacts 3-4
if str2double(expmt_list{expmt}.cuff_list(session,1))==3
    stim_chan = 9;
else
    stim_chan = 1;
end
if strcmpi(stim_port,'a')
    stimChan = stim_chan;
elseif strcmpi(stim_port,'b')
    stimChan = 128 + stim_chan;
elseif strcmpi(stim_port,'c')
    stimChan = 256 + stim_chan;
elseif strcmpi(stim_port,'d')
    stimChan = 384 + stim_chan;
end

% Neuroshare function to read stimulation events file, rec_filePath contains directory to file starting
    %with "C:/...", stimChan has channel ID to read from
stimTimes = read_stimEvents(rec_filePath, stimChan);

end


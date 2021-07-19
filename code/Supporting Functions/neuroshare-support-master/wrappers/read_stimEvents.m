function [stimEvts, channels, stimWF] = read_stimEvents(dataPath, channels)
    
    %   DESCRIPTION
    %   ===================================================================
    %   Reads stim events for input channels from specified NEV file. 
    %   REMEMBER the stim channel numbering changes based on the FE port 
    %   where the nanostim is connected. Stim channels are numbered as 
    %   follows: (128*(port-1)+chan) where port A=1, B=2,...
    %
    %   INPUTS
    %   ===================================================================
    %   dataPath        : (char) path to data file
    %   channels        : (int) 1xn vector of stim channels
    %
    %   OUTPUT
    %   ===================================================================
    %   stimEvts      : (1xn) cell array of stim times
    %   channels      : (1xn) vector of stim channels read
    %   stimWF        : (1xn) cell array of stim voltage waveforms
    %
    %   ACN created 11/16 
    %   ACN modified 2/17
    
    if exist(dataPath)
        [~,hFile] = ns_OpenFile(dataPath);
    else
        error('File not found')
    end
    
    if isempty(channels)
        channels = [hFile.Entity(cellfun(@(x) ~isempty(strfind(x, 'stim')), {hFile.Entity.Label})).ElectrodeID]-5120;
    end
    if ~isempty(find(channels==0))
        channels(find(channels==0)) = []; 
    end
    
    numChannels = length(channels);
    stimEvts = cell(1,numChannels);
    stimWF = cell(1,numChannels);
    
    for iChan = 1:numChannels
        fprintf('Reading nev file, channel: %02d\n', channels(iChan))
        chanEntityIdx = find([hFile.Entity(:).ElectrodeID] == channels(iChan) + 5120);
        if isempty(chanEntityIdx)
            disp({hFile.Entity.Label})
            error('Could not find entity %d. Entities found in file are listed above', channels(iChan))
        else
            numStimEvts = hFile.Entity(chanEntityIdx).Count;
            stimEvts{iChan} = zeros(1,numStimEvts);
            for i = 1:numStimEvts
                [~,stimEvts{iChan}(i),stimWF{iChan}(i,:),~] = ns_GetSegmentData(hFile, chanEntityIdx, i);
            end
        end
    end
    ns_CloseFile(hFile);
end
function [digitalData] = read_digitalEvents(dataPath, channels)
    
    %   DESCRIPTION
    %   ===================================================================
    %   Reads digital events for input channels from specified NEV file. 
    %   
    %   INPUTS
    %   ===================================================================
    %   dataPath        : (char) path to data file
    %   channels        : (int) 1xn vector of SMA channels
    %
    %   OUTPUT
    %   ===================================================================
    %   digitalData      : (1xn) cell array of stim times
    %
    %   ACN created 11/16 
    %   ACN modified 2/17
    
    if exist(dataPath)
        [~,hFile] = ns_OpenFile(dataPath);
    else
        error('File not found')
    end
    
    numChannels = length(channels);
    digitalData = struct;
    
    for iSMA = 1:numChannels
        chanLabel = sprintf('SMA %d', channels(iSMA));
        fprintf('Reading nev file, channel: %02d\n', channels(iSMA))
        chanEntityIdx = find(cellfun(@(x) strcmp(x, chanLabel), {hFile.Entity.Reason}));

        if isempty(chanEntityIdx)
            disp({hFile.Entity.Label})
            error('Could not find entity %d. Entities found in file are listed above', channels(iChan))
        else
            [~, entityInfo] = ns_GetEntityInfo(hFile, chanEntityIdx);

            numCount = entityInfo.ItemCount;
            digitalData.data{iSMA} = NaN(1, numCount);
            digitalData.timeStamp{iSMA} = NaN(1, numCount); 
            digitalData.dataSize{iSMA} = NaN(1, numCount);
            for i = 1:numCount
                [~, digitalData.timeStamp{iSMA}(i), digitalData.data{iSMA}(i), digitalData.dataSize{iSMA}(i)] = ns_GetEventData(hFile, chanEntityIdx, i);
            end 
        end
    end
    ns_CloseFile(hFile);
        
end
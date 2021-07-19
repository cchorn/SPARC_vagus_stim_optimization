function [analogData, timeVec] = read_continuousData(dataPath, datastream, channels)
    
    %   DESCRIPTION
    %   ===================================================================
    %   Reads continuous data for input channels from specified NS5/NF3 file 
    %
    %   INPUTS
    %   ===================================================================
    %   dataPath        : (char) path to data file
    %   dataStream      : (char) 'hi-res'/'raw'/'analog'
    %   channels        : (int) 1xn vector of channels to be read
    %
    %   OUTPUT
    %   ===================================================================
    %   analogData      : (nxm) array of data points for n channels
    %   timeVec         : (1xm) time vector
    %
    %   ACN created 11/16 
    %   ACN modified 2/17
    
    if exist(dataPath)
        [~,hFile] = ns_OpenFile(dataPath);
    else
        error('File not found')
    end
    
    startIdx = 1;
    numChannels = length(channels);
    if strcmp(datastream, 'hi-res')
        fs = 2e3;
        type_mask = ismember({hFile.FileInfo.Type}, 'nf3');
        endIdx = hFile.FileInfo(type_mask).TimeStamps(2);
        timeVec = linspace(0,endIdx/fs,endIdx);
    elseif strcmp(datastream, 'raw')
        fs = 30e3;
        type_mask = ismember({hFile.FileInfo.Type}, 'ns5');
        endIdx = hFile.FileInfo(type_mask).TimeStamps(2);
        timeVec = linspace(0,endIdx/fs,endIdx);
    elseif ismember(datastream, 'analog')                       % analog defaults to 30k. 1K is extracted only if 30k does not exist
        type_mask = ismember({hFile.FileInfo.Type}, 'ns5');
        if any(type_mask)
            fs = 30e3;
            endIdx = hFile.FileInfo(type_mask).TimeStamps(2);
        else
            warning('30K analog channel not found. Checking for 1K channel')
            type_mask = ismember({hFile.FileInfo.Type}, 'ns2');
            if any(type_mask)
                fs = 1e3;
                endIdx = hFile.FileInfo(type_mask).TimeStamps(2);
            else
                error('No 1K data found for this channel')
            end
        end
        timeVec = linspace(0,endIdx/fs,endIdx);
    else
        error('Invalid data stream')
    end
    analogData = zeros(numChannels, endIdx);
    
    for iChan = 1:numChannels
        fprintf('Reading %s file, channel: %02d\n', hFile.FileInfo(type_mask).Type, channels(iChan))
        chanEntityIdx = find(strcmp({hFile.Entity.Label}, sprintf('%s %d',datastream, channels(iChan))));
        if isempty(chanEntityIdx)
            disp({hFile.Entity.Label})
            error('Could not find entity %s %d. Entities found in file are listed above',datastream, channels(iChan))
        else
            [~, ~, analogData(iChan, :)] = ns_GetAnalogData(hFile, chanEntityIdx, startIdx, endIdx);
        end
    end
    ns_CloseFile(hFile);
    
end
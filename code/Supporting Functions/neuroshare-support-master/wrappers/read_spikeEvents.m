function [wfSnip, Tstamp, channels] = read_spikeEvents(dataPath, FEport)

    %   DESCRIPTION
    %   ===================================================================
    %   Reads spike events for input channels from specified NEV file. 
    %
    %   INPUTS
    %   ===================================================================
    %   dataPath        : (char) path to data file
    %   FEport          : (string) front end port for spike recording
    %
    %   OUTPUT
    %   ===================================================================
    %   wfSnip          : (1xn) cell array of spike waveforms
    %   Tstamp          : (1xn) cell array of stim times
    %   channels        : (1xn) vector of stim channels read
    %
    %   ACN created 3/19
    
    [~,hFileEV] = ns_OpenFile(dataPath);
    channels = 1:32;
    numChans = 32;
    Tstamp = cell(1,numChans);
    wfSnip = cell(1,numChans);
            
    for iChan = 1:numChans
        fprintf('Reading Nodose snippet %02d\n', channels(iChan))
        chanLabel = sprintf('%s-%02d spk', FEport, channels(iChan));          
        chanIdx = find(strcmp({hFileEV.Entity.Label}, chanLabel));
        if ~isempty(chanIdx)
            numSnip = hFileEV.Entity(chanIdx).Count;
            Tstamp{iChan} = zeros(1, numSnip);
            wfSnip{iChan} = zeros(52, numSnip);
            for iSnip = 1:numSnip
                [~,Tstamp{iChan}(iSnip), wfSnip{iChan}(:,iSnip), ~, ~] = ns_GetSegmentData(hFileEV, chanIdx,iSnip);
            end
        end
    end
    
    ns_CloseFile(hFileEV);
end

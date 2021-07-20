function [commandedSteps, commandedRes]  = getStimSteps(commandedAmp_uA)

    stimLookup12510 = [[zeros(3,127); 1:127]';
                [zeros(2,127); 1:127; 127*ones(1,127)]';
                [zeros(1,127); 1:127; 127*ones(1,127); 127*ones(1,127)]';
                [1:127; 127*ones(1,127); 127*ones(1,127); 127*ones(1,127)]'];
            
    stimLookup20 = [[zeros(3,75); 1:75]';
                [zeros(2,75); 1:75; 75*ones(1,75)]';
                [zeros(1,75); 1:75; 75*ones(1,75); 75*ones(1,75)]';
                [1:75; 75*ones(1,75); 75*ones(1,75); 75*ones(1,75)]'];

    stimAmp1 = sum(stimLookup12510,2)*1;
    stimAmp2 = sum(stimLookup12510,2)*2;
    stimAmp5 = sum(stimLookup12510,2)*5;
    stimAmp10 = sum(stimLookup12510,2)*10;
    stimAmp20 = sum(stimLookup20,2)*20;

    allPossibleAmps = unique([stimAmp1; stimAmp2; stimAmp5; stimAmp10; stimAmp20]);
    if ~ismember(commandedAmp_uA, allPossibleAmps)
        [~, ix] = min(abs(allPossibleAmps - commandedAmp_uA));
        errordlg(sprintf('invalid stim amp. closest allowable amps: %d or %d', allPossibleAmps([ix, ix+1])))
        error('invalid stim amp. closest allowable amps: %d or %d', allPossibleAmps([ix, ix+1]))
    end
    
    numAllAmps = length(allPossibleAmps);
    stimSteps = cell(numAllAmps, 5);
    for i = 1:numAllAmps
        ampVal = allPossibleAmps(i);

        if any(stimAmp1 == ampVal)
            stimSteps{i,1} = stimLookup12510(stimAmp1 == ampVal,:);
        end

        if any(stimAmp2 == ampVal)
            stimSteps{i,2} = stimLookup12510(stimAmp2 == ampVal,:);
        end
        if any(stimAmp5 == ampVal)
            stimSteps{i,3} = stimLookup12510(stimAmp5 == ampVal,:);
        end

        if any(stimAmp10 == ampVal)
            stimSteps{i,4} = stimLookup12510(stimAmp10 == ampVal,:);
        end

        if any(stimAmp20 == ampVal)
            stimSteps{i,5} = stimLookup20(stimAmp20 == ampVal,:);
        end
    end
    
    possibleSteps = stimSteps(allPossibleAmps == commandedAmp_uA,:);
    commandedRes = max(find(cellfun(@(x) ~isempty(x), possibleSteps)));
    commandedSteps = possibleSteps{commandedRes};
end
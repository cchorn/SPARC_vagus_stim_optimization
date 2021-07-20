function [] = stimseqDesign(cathElec, anodElec, cathAmp_uA, preQuietT_s, stimT_s, postQuietT_s, freq, phaseDur_ms, varargin)
    %   DESCRIPTION
    %   ===================================================================
    %   File save start and stop and stimulation design for ferret stim 
    %   trials.
    %  
    %   INPUTS
    %   ===================================================================
    %   cathElec        :   (int) cathode electrode 1 to 8
    %   anodElec        :   (int) anode electrode 1 to 8
    %   cathAmp_uA      :   (int) stimulation amplitude in uA
    %   preQuietT_s     :   (numeric) pre stim duration in seconds
    %   stimT_s         :   (numeric) stim duration in seconds
    %   postQuietT_s    :   (numeric) post stim duration in seconds
    %   freq            :   (int) stimulation frequency
    %   phaseDur_ms     :   (int) phase duration of stim pulse 
    %                             (see assumption 4)
    %
    %   OPTIONAL INPUTS
    %   ===================================================================
    %   digTrigger_enable : (boolean) trigger rec start/stop using SMA I/O (default: false)
    %   digTrigger_chan   : (numeric) SMA channel to trigger rec start/stop (default: 1)
    %   stimFE            : (string) A/B/C/D (default: B)
    %
    %   ASSUMPTIONS
    %   ===================================================================
    %   1) nano2+stim is connected to the NIP FE port B i.e. stim channels 
    %      are numbered 129,130,131...
    %   3) 4 channels are shorted together such that nano2 --> elec is:
    %      1:4 --> 1, 5:8 -->2, ...
    %   4) stim pulses are biphasic symmetric
    %
    %   EXAMPLE
    %   ===================================================================
    %   stimDesign(1, 2, 5000, 60, 5*60, 5*60, 15, 0.1)
    %
    %   this will set contact 1 as a cathode, 2 as anode, 5000uA 
    %   stimulation on nano 2 channels 1,2,3,4. quiet recording for 1 min, 
    %   stim for 5 mins, post stim quiet for 5 mins. stim frequency is 
    %   15 Hz and stim pulse width is 100us
    %
    %   NOTE
    %   ===================================================================
    %   Both xippmex(stim) and xippmex(stimseq) do not allow programming
    %   the stimulator for more than 3 and 4 minutes respectively. This
    %   function constructs a stim pulse for one minute and runs a for loop
    %   to deliver multiple minutes of stimulation
    %
    %   U18 STIM PROTOCOL
    %   ===================================================================
    %   % ch1-2
    %   stimseqDesign(1, 2, 100, 10, 2*60, 5*60, 15, 0.1);
    %   stimseqDesign(1, 2, 1000, 10, 2*60, 5*60, 15, 0.1);
    %   stimseqDesign(1, 2, 5000, 10, 2*60, 5*60, 15, 0.1);
    %   stimseqDesign(1, 2, 1000, 10, 2*60, 5*60, 20, 0.1);
    %   stimseqDesign(1, 2, 1000, 10, 2*60, 5*60, 30, 0.1);
    
    %DEFINE_CONSTANTS
    digTrigger_enable = false;
    digTrigger_chan   = 1;
    stimFE = 'B';
    %END_DEFINE_CONSTANTS
    
    [cathAmp_step, stimRes] = getStimSteps(cathAmp_uA);
    
    if strcmpi(stimFE,'A')
        elecShift = 0;
    elseif strcmpi(stimFE,'B')
        elecShift = 128;
    elseif strcmpi(stimFE,'C')
        elecShift = 256;
    elseif strcmpi(stimFE,'D')
        elecShift = 384;
    end
    
    status=1;
    %status = xippmex;
    if status == 0
        error('We have a problem...You dont say?')
    end
    %xippmex('stim', 'enable',0)
    %xippmex('digout',4,0)
    
    if digTrigger_enable
        xippmex('digout', digTrigger_chan, 0);
        pause(1);
    end
    
    shortedChans = 4;                                                           % shorted channels
    nano_cath = elecShift+[((cathElec-1)*shortedChans+1):cathElec*shortedChans];
    nano_anod = elecShift+[((anodElec-1)*shortedChans+1):anodElec*shortedChans];
    
    stimElectrodes = [nano_cath, nano_anod];
    %if xippmex('stim','res',stimElectrodes) ~= stimRes
    %    xippmex('stim','res',stimElectrodes, stimRes);
    %    pause(1)
    %end
    %xippmex('stim', 'enable',1)
    
    time_period_samples = 1/freq*30e3;
    phase_step = phaseDur_ms*1e-3*30e3;
    numRepeats = freq*stimT_s;
    if numRepeats > 4094
        numLoops = ceil(numRepeats/4094);
        repeatsMat = [4094*ones(1,numLoops-1), rem(numRepeats,4094)];
    else
        numLoops = 1;
        repeatsMat = numRepeats;
    end
    
    pulseStruct_cat = [struct('length', phase_step, 'ampl', num2cell(cathAmp_step), 'pol', 0, ...       % negative first
                            'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);...
                        struct('length', phase_step, 'ampl', num2cell(cathAmp_step), 'pol', 1, ...
                            'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1)]';
    
    pulseStruct_anod = [struct('length', phase_step, 'ampl', num2cell(cathAmp_step), 'pol', 1, ...      % positive first
                            'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1);...
                         struct('length', phase_step, 'ampl', num2cell(cathAmp_step), 'pol', 0, ...
                            'fs', 0, 'enable', 1, 'delay', 0, 'ampSelect', 1)]';
                        
    cmd1 = struct('elec', num2cell(nano_cath), 'period', time_period_samples, 'repeats', [], 'seq', []);
    cmd2 = struct('elec', num2cell(nano_anod), 'period', time_period_samples, 'repeats', [], 'seq', []);
    cmd = [cmd1 cmd2];
    
    for iChan = 1:4
        cmd(iChan).seq = pulseStruct_cat(iChan,:);
        cmd(iChan+4).seq = pulseStruct_anod(iChan,:);
    end

    if digTrigger_enable
        xippmex('digout', digTrigger_chan, 1);
        disp('Recording started...')
        pause(1)
    else
        disp('Starting trial. Verify that remote control is enabled')
    %    xippmex('trial', 'recording');
    end
    
    disp('Pre-stim quiet recording')
    pause(preQuietT_s)
    
    disp('starting stim')
    %xippmex('digout',4,1)
    for iLoop = 1:numLoops
        for iCmd = 1:length(cmd)
            cmd(iCmd).repeats = repeatsMat(iLoop);
        end
        xippmex('stimseq', cmd);
        pause(repeatsMat(iLoop)*time_period_samples/30e3)
    end
    disp('stim ended')
    xippmex('digout',4,0)
    
    disp('Post-stim quiet recording')
    pause(postQuietT_s-1)
    
    if digTrigger_enable
        xippmex('digout', digTrigger_chan, 0);
        disp('Recording stopped, end of trial') 
        pause(1)
    else
        disp('Recording stopped')
        xippmex('trial', 'stopped');
    end
end

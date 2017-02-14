% Copyright 2017 Chih-Yu (John) Yang
% This script read in abf files and do some simple analyses as an example

% Set clampexFolderPath to the folder path where you store your recordings
% (abf files). Usually if you set the Clampex correct, there will be
% sequential numbers at the suffix. Select the number that you want to
% analyse and assign it to patchNum.
    
    % clampexFolderPath = 'your folder path'
    % patchNum = a number that is your abf file recording number that you
    % want to analyse
    
    filePath = ClampexData.findFile(clampexFolderPath,patchNum);
    
    assignNo = [];
    
    % Clamp will be ClampexData object, which will have the attributes and
    % functions specified in ClamexData.mat
    clamp = patchData(filePath,assignNo,'e');
    
    % Here are some of the recording protocols that I have used. You should
    % see if they match your protocols and run them.
    analyseMethod = 'VC input'; % 'IC AP family'
    
    % Usually I write in Excel file the patchNum and the protocol, then I
    % can do an automatic check of what the recording protocol is.
    % Here, just to show how data is 
    if strcmp(analyseMethod, 'VC input')
        % VC input recording that is aimed to get the membrane resistance
        patch = clamp.patch;
        base = median(patch);
        
        cPatch = patch - base; % correction of VC input potential to 0
        nPeak = -max(abs(cPatch));
        pPeak = max(cPatch);
        transientV = (pPeak - nPeak)/2;
        VCInputTransient = transientV;
        
    elseif ~isempty(strfind(protocolCheck{1}, 'IC AP family'))
        %% Current clamp to record cell responses to intracellular injection
        stimTime = [];
        clamp.getStimTime(stimTime);
        % patch = clamp.patch;
        clamp.getIntSpike(30,0.001); % first input = relative spike threshold from baseline (mV), second is the minimum time before detecting the next spike
        
        spike = clamp.spikeData; % A cell format storing the spike time data point, spike amplitude and spike time in time
        
    elseif ~isempty(strfind(protocolCheck{1}, 'LED'))
        %% Current clamp to record cell responses to light stimulation. I think only people involved in visual neuroscience will need this
        stimTime = [];
        clamp.getStimTime(stimTime);
        clamp.getLightSpike(30,0.001); % first input = relative spike threshold from baseline (mV), second is the minimum time before detecting the next spike
       
        spike = clamp.spikeData; % A cell format storing the spike time data point, spike amplitude and spike time in time
        
        
    elseif ~isempty(strfind(protocolCheck{1}, 'SPS'))
        %% Current clamp to record cell responses to single pulse stimulation
        clamp.getSPSSpike(30,0.001); % first input = relative spike threshold from baseline (mV), second is the minimum time before detecting the next spike
        clamp.getSpikeRate;
        
        spike = clamp.spikeData; 
        
    end
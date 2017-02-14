% Copyright 2017 Chih-Yu (John) Yang
% This Class inherits from ClampexData. It is mainly for the analysis of
% retinal ganglion cell responses to electrical stimulation. There are also
% several plot functions that I have written to facilitate the
% visualisation of data.


classdef patchData < ClampexData
    
    % patchData is a subclass of ClampexData
    
    
    properties
        
        spikeData = {};
        % Spike data contains the amplitude and timing of the spikes. It is
        % organised into 3 rows of cells. The 1st row of cells contains the
        % timing and column number represents the sweep number. 2nd row
        % contains the amplitude. The 3rd row contains the data point
        % number. At the moment I set the cut-off threshold for spike to be
        % Vmrest + 40mV for episodic stimulation. For gap-free recording, I
        % set the cut-off threshold for spikes as Vmrest + 40mV as well
        spikeletData = {};
        % For spikes that have smaller amplitude than the usual spikes. It
        % is organised similarly to spikeData.At the moment I set the
        % cut-off threshold for spike to be Vmrest + 10mV.For gap-free
        % recording, I set the cut-off threshold for spikelets as
        % Vmrest + 10mV as well
        
        baseSpike = {};
        % Further sorting of spikes into baseline spikes, stim spikes and
        % after stimulation spikes
        baseSpikelet = {};
        
        stimSpike = {};
        stimSpikelet = {};
        
        AStimSpike = {};
        AStimSpikelet = {};
        
        % gdSpikeData = {};
        
        imSig = '';
        % the object image signals from imSignal (a subclass of CaImage)
        
        
        framesPerSweep = [];
        % The time of each integer frame number
        
        
        F0 = [];
        DF = [];
        DFF = [];
        % data that go through the DF/F0 process. F0 is the baseline
        % fluorescence before the stimulation. DF = F - F0
        
        sweepFrame = [];
        % The frame number in each sweep
        sweepFrameTime = [];
        
        %% Intracellular stimulation protocols
        stimType = '';
        % either LED, intStim or extStim
        stimInitAmp = '';
        stimDAmp = '';
        stimAmp = [];
        % stimulation amplitude is deduced from stimInitAmp and stimDAmp
        intStimS = [];
        intStimE = [];
        
        intStimSTime = [];
        intStimETime = [];
        
        stimDataPoint = [];
        
        intStimSpikeNumber;
        %
        %         framesPerStim = '';
        %         % Number of frames during stimulatin
        %         stimFrame = [];
        %         % The frame numbers when stimulation is given. It is organised in
        %         % nSweeps x framesPerStim
        %         stimFrameTime = [];
        %         % The time of the stimulation frames
        %         baselineFrame = [];
        %         % The baseline frames before the start of the stimulation
        %         baselineFrameTime = [];
        
        
        
        %% Baseline times and properties
        baselineTime = {};
        % The initial baseline time obtained from timing the data points
        % by the data point acquisition time
        baselineDataPoint = [];
        % The initial baseline data points recorded before stimulation. It
        % is organised into nSweeps x nBaselineDataPoints
        baselineMax = [];
        restVm = [];
        % resting membrane voltage at each sweep. It is found by using min
        % function during for Vm within the baseline time
        
        %% Stimulation times and properties
        
        stimTime = {};
        stimStartTime = [];
        stimEndTime = [];
        stimStartDataPoint = [];
        stimEndDataPoint = [];
        
        stimSpikeNumber = [];
        
        intStimAmp = [];
        % For recording intracellular stimulation amplitudes
        
        intStimDur = [];
        % For recording intracellular stimulation durations
        
        stimMax = [];
        % The maximum membrane voltage reached during stimulation. This is
        % mainly for depolarising stimulation
        stimMin =[];
        % The minimum membrane potential during stimulation. It should find
        % the raised and sustained potential during depolarising
        % stimulation. For hyperpolarising stimulation, median is used to
        % find the hyperpolarising potential.
        
        LEDONTime = {};
        % The time when LED is turned on
        LEDONDataPoint = {};
        % The data points when LED is turned on
        LEDOFFTime = {};
        % The time points when LED is turned off (in the LED sequence phase)
        LEDOFFDataPoint = {};
        % The data points when LED is turned off (in the LED sequence phase)
        
        LEDONSpike = [];
        LEDONSpikelet = [];
        LEDOFFSpike = [];
        LEDOFFSpikelet = [];
        

        
        nStim = 0;
        % Note that patch nStim is different from CaImage nStim. This
        % problem occurs because sometimes imaging stops before stimulation
        % stops, and my definition of nStim in CaImage is the number of
        % sitmulation that occured during imaging. Any stimulation outside
        % of imaging time is not included in nStim.
        
        %% AStimulation times and properties
        
        AStimTime = {};
        AStimDataPoint = [];
        AStimMax = [];
        % The maximum membrane voltage reached after stimulation
        AStimMin = [];
        % The minimum membrane voltage after stimulation, this is mainly to
        % find the resting membrane potential after stimulation
        
        folderPath = '';
        % the folder path name for batch analysing files within the folder
        nFile = '';
        % number of files in the folder for batch analysis
        
        spikeCutOff = [];
        % organise into nSweeps x 3. The 1st column value is used for
        % baseline cut-off, 2nd for stim cut-off and 3rd for
        % afterhyperpolarisation
        
        spikeletCutOff = [];
        
        %% Stimlation time peaks (both artefacts and possible spikes)
        stimPosPeak = {};
        
        stimNegPeak = {};
        
        stimBaseline = [];
        
        
        % For spike rate plotting
        
        % Individual spike rate = time different between current AP and
        % previous AP (Ind).
        
        % General spike rate = numbers of APs in each section then divided
        % by time (Gene).
        
        % Common parameters for all different stimulation protocols
        baselineIndSpikeRate = [];
        AStimIndSpikeRate = [];
        
        baselineGeneSpikeRate
        AStimGeneSpikeRate
        
        % Specific parameters for each different stimulation protocol
        
        %% For intracellular stimulation
        stimIndSpikeRate = [];
        
        stimGeneSpikeRate
        
        %% For LED stimulation
        LEDONIndSpikeRate = [];
        LEDOFFIndSpikeRate = [];
        
        LEDONGeneSpikeRate;
        LEDOFFGeneSpikeRate;
        
        LEDONSpikeNumber;
        LEDOFFSpikeNumber;
        
        %% For HFS stimulation
        
        HFSIndSpikeRate = [];
        HFSGeneSpikeRate;
        
        HFSMedianMV;
        
        spikeRate = {};
        spikeRateTime = {};
        
        nSpikeRatePt = 0;
        spikeRateWindowSize = 0;
        
        % ExtStim properties
        patchExtStimCorrected = [];
        extStimSpikeDataPt = [];
        
        spikeLatency = {}; % The spike latency after stimulation
        
        HFSfilteredSig = [];
        FFTFreq = 0;
        
        % HFS Stimulation
        HFSDataPt = [];
        HFSTime = [];
        indSpikeRate = [];
        generalSpikeRate = [];
        
        HFSGapDataPoint;
        HFSGapTime;
        
        HFSGapIndSpikeRate;
        HFSGapGeneSpikeRate;
        
        HFSSpikeNumber;
        
        HFSBaseSponSpikeN;
        % HFS spontaneous spike number during baseline
        
        HFSAStimSponSpikeN;
        % HFS spontaneous spike number after stimulation
        
        HFSGapSpikelet;
        
        HFSSpike = {};
        % HFS spikes in original signal
        HFSFiltSpike = {};
        % HFS spikes in filtered signal
        
        HFSGapSpike = {};
        % The interval between each HFS 
        
        HFSSweepMedianMV = [];
        
        HFSGapSpikeNumber;
        
        HFSSpikelet = {};
        
        HFSPotential;
        % the median value of membrane potential during HFS
        
        HFSPotentialPerAmp;
        
        HFSDiffVm;
        % Difference between resting membrnae potential and mbm potential
        % during HFS
        
        
        protocol = '';
        
        
        %% For SPS stimulation
        SPSDataPt;
        SPSIndSpikeRate;
        SPSGeneSpikeRate;
        SPSSpikeNumber;
        
        SPSSpikelet;
        SPSSpike;
        
        SPSTime;
        
        SPSGapDataPoint;
        SPSGapIndSpikeRate;
        SPSGapGeneSpikeRate;
        
        stimulatorPos;
        % Stimulator position for mapping 
        
        noExtStimSweepN

        perfusion;
        % Perfusion type. 0 = normal, 1 = others to be specified
        
        JPCVal = 0;
        % junction potential correction; default value is 0
        
        
        sweepMedianMV = [];
        % Median of membrane voltage of every sweep
        
        avgMedianMV;
        % Average of all the medians of membrane voltage of recording sweeps.
        
        stdMedianMV;
        % Standard deviation of membrane voltage of all the recording
        % sweeps
        
        SEConfig;
        % Stimulating electrode configuration used
        
    end
    
    
    properties(Dependent, SetAccess=private)
        %         LEDONSpike
        %         LEDONSpikelet
        %         LEDOFFSpike
        %         LEDOFFSpikelet
    end
    
    % Data -> getStimTime for all intracellular, extracellular and LED
    % stimulation -> getSpike and sort the spikes into spikes that occured
    % in baseline, stimulation, or after stimulation time
    
    methods(Static)
        
        function obj = loadobj(s)
            if isstruct(s)
                obj = patchData;
                obj = reload(obj, s);
                obj.stimType = s.stimType;
                obj.spikeData = s.spikeData;
                obj.spikeletData = s.spikeletData;
                
                obj.spikeCutOff = s.spikeCutOff;
                obj.spikeletCutOff = s.spikeletCutOff;
                
                
                obj.stimStartTime = s.stimStartTime;
                obj.stimEndTime = s.stimEndTime;
                obj.stimStartDataPoint = s.stimStartDataPoint;
                obj.stimEndDataPoint = s.stimEndDataPoint;
                
                obj.baseSpike = s.baseSpike;
                obj.baseSpikelet = s.baseSpikelet;
                
                obj.stimSpike = s.stimSpike;
                obj.stimSpikelet = s.stimSpikelet;
                
                obj.AStimSpike = s.AStimSpike;
                obj.AStimSpikelet = s.AStimSpikelet;
                
                obj.baselineTime = s.baselineTime;
                obj.baselineDataPoint = s.baselineDataPoint;
                
                
                obj.stimTime = s.stimTime;
                obj.stimDataPoint = s.stimDataPoint;
                
                obj.AStimTime = s.AStimTime;
                obj.AStimDataPoint = s.AStimDataPoint;
                
                
                obj.restVm = s.restVm;
                
                obj.LEDONTime = s.LEDONTime;
                obj.LEDONDataPoint = s.LEDONDataPoint;
                obj.LEDOFFTime = s.LEDOFFTime;
                obj.LEDOFFDataPoint = s.LEDOFFDataPoint;
                obj.nStim = s.nStim;
                
                obj.spikeRate = s.spikeRate;
                obj.spikeRateTime = s.spikeRateTime;
                
                obj.nSpikeRatePt = s.nSpikeRatePt;
                obj.spikeRateWindowSize = s.spikeRateWindowSize;
                obj.patch = s.patch;
                obj.wholeData = s.wholeData;
                obj.dataTime = s.dataTime;
                
            else
                obj = s;
            end
            
        end
        
    end
    
    
    methods
        function obj = patchData(path, assignNo, abfEnd, varargin)
            
            
            obj@ClampexData(path,assignNo, abfEnd, varargin{:});
            
            %             if obj.stimFlag == 1
            %                obj.getStimTime;
            %             end
            %             obj.imSig = imSig;
        end
        function JPCorr(obj)
            obj.patch = obj.patch - obj.JPCVal;
            
        end
        
        
        %         function ONSpike = get.LEDONSpike(obj)
        %             ONSpike = obj.stimSpike;
        %         end
        %         function ONSpikelet = get.LEDONSpikelet(obj)
        %             ONSpikelet = obj.stimSpikelet;
        %         end
        %         function OFFSpike = get.LEDOFFSpike(obj)
        %             OFFSpike = obj.AStimSpike;
        %         end
        %         function OFFSpikelet = get.LEDOFFSpikelet(obj)
        %              OFFSpikelet = obj.AStimSpikelet;
        %         end
        
        function save(obj, saveName)
            s.stimType = obj.stimType;
            
            %             s.spikeData = obj.spikeData;
            %             s.spikeletData = obj.spikeletData;
            
            s.spikeCutOff = obj.spikeCutOff;
            s.spikeletCutOff = obj.spikeletCutOff;
            
            
            s.stimStartTime = obj.stimStartTime;
            s.stimEndTime = obj.stimEndTime;
            s.stimStartDataPoint = obj.stimStartDataPoint;
            s.stimEndDataPoint = obj.stimEndDataPoint;
            
            s.baseSpike = obj.baseSpike;
            s.baseSpikelet = obj.baseSpikelet;
            
            s.stimSpike = obj.stimSpike;
            s.stimSpikelet = obj.stimSpikelet;
            
            s.AStimSpike = obj.AStimSpike;
            s.AStimSpikelet = obj.AStimSpikelet;
            
            s.baselineTime = obj.baselineTime;
            s.baselineDataPoint = obj.baselineDataPoint;
            
            
            s.stimTime = obj.stimTime;
            s.stimDataPoint = obj.stimDataPoint;
            
            s.AStimTime = obj.AStimTime;
            s.AStimDataPoint = obj.AStimDataPoint;
            
            
            s.restVm = obj.restVm;
            
            s.LEDONTime = obj.LEDONTime;
            s.LEDONDataPoint = obj.LEDONDataPoint;
            s.LEDOFFTime = obj.LEDOFFTime;
            s.LEDOFFDataPoint = obj.LEDOFFDataPoint;
            s.nStim = obj.nStim;
            
            s.spikeRate = obj.spikeRate;
            s.spikeRateTime = obj.spikeRateTime;
            
            s.nSpikeRatePt = obj.nSpikeRatePt;
            s.spikeRateWindowSize = obj.spikeRateWindowSize;
            s.patch = obj.patch;
            s.extStim = obj.extStim;
            s.intStim = obj.intStim;
            s.LED = obj.LED;
            %             s.wholeData = obj.wholeData;
            %             s.dataTime = obj.dataTime;
            %             s.recordTime = obj.recordTime;
            s.spikeLatency = obj.spikeLatency;
            
            s.SR = obj.SR;
            
            s.HFSfilteredSig = obj.HFSfilteredSig;
            s.FFTFreq = obj.FFTFreq;
            
            save(saveName, 's');
            
        end
        
        function obj = reload(obj, s)
            obj = reload@ClampexData(obj,s);
        end
        
        function getStimTime(obj,varargin)
            % extract the stimulation protocols if they are recorded
            
            
            stimTime = [];
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'stimTime')
                    stimTime = varargin{i + 1};
                end
            end
            
            %% Intracellular stimulation
            if ~isempty(obj.intStim)
                nSweep = 0;
                if ndims(obj.wholeData) == 3
                    nSweep = size(obj.wholeData,3);
                else
                    nSweep = 1;
                end
                for i = 1:nSweep
                    intStimTemp = obj.wholeData(:,2,i);
                    intStimBaseValue = median(intStimTemp);
                    intStimTemp = intStimTemp - intStimBaseValue;
                    % sometimes the recorded injected current base is not zero,
                    % have to correct this.
                    
                    %                 stimPeak = findpeaks(intStimTemp,'MINPEAKHEIGHT',5);
                    
                    diffIntStim = diff(abs(intStimTemp));
                    stimSTemp = find(diffIntStim > 5);
                    diffstimSTemp = diff([0, stimSTemp']);
                    stimS = stimSTemp(diffstimSTemp  > 1);
                    
                    [~, maxI] = max(intStimTemp(stimS + 5));
                    stimS = stimS(maxI);
                    
                    stimETemp = find(diffIntStim < -5);
                    diffstimETemp = diff([stimETemp', length(intStimTemp)]);
                    stimE = stimETemp(diffstimETemp  > 1);
                    
                    [~, maxI] = max(intStimTemp(stimE - 5));
                    stimE = stimE(maxI);
                    
                    if ~isempty(stimS)
                        obj.intStimAmp(i) = intStimTemp(stimS + 5);
                        obj.stimStartDataPoint(i) = stimS;
                        obj.stimEndDataPoint(i) = stimE;
                        
                        if nSweep == 1
                            obj.stimStartTime(i) = stimS/obj.SR;
                            obj.stimEndTime(i) = stimE/obj.SR;
                        else
                            obj.stimStartTime(i) = obj.gapFreeTime(stimS);
                            obj.stimEndTime(i) = obj.gapFreeTime(stimE);
                        end
                    elseif ~isempty(obj.stimStartDataPoint)
                        obj.intStimAmp(i) = 0;
                        stimS = obj.stimStartDataPoint(i - 1);
                        stimE = obj.stimEndDataPoint(i - 1);
                        obj.stimStartDataPoint(i) = stimS;
                        obj.stimEndDataPoint(i) = stimE;
                        
                        obj.stimStartTime(i) = obj.gapFreeTime(stimS);
                        obj.stimEndTime(i) = obj.gapFreeTime(stimE);
                        
                    else
                        stimS = 25000;
                        stimE = 50000;
                        obj.stimStartDataPoint(i) = stimS;
                        obj.stimEndDataPoint(i) = stimE;
                        
                        obj.stimStartTime(i) = obj.gapFreeTime(stimS);
                        obj.stimEndTime(i) = obj.gapFreeTime(stimE);

                    end
                    % several data point addition to ensure the intracellular
                    % amplitude obtained is on the flat stim onset time
                    %                     if (~isempty(obj.intStimAmp > 0))&(~isempty(obj.intStimAmp < 0))
                    %                         negI = find(obj.intStimAmp < 0);
                    %                         obj.intStimAmp(i) =...
                    %                             [obj.intStimAmp(1:negI(end)), 0, obj.intStimAmp((negI(end) + 1):end)];
                    %                     end
                    
                    
                    
                    
                end
                
                for i = 1:length(obj.stimStartTime)
                    obj.intStimDur(i) = obj.stimEndTime(i) - obj.stimStartTime(i);
                end
                obj.stimType = 'intStim';
            end
            
            
            %% Extracellular stimulation
            if ~isempty(obj.extStim)
                if ndims(obj.wholeData) == 3
                    nSweep = size(obj.wholeData,3);
                else
                    nSweep = 1;
                end
                % The detection level is quite sensitive to the noise, so
                % might need to constantly change it. Also, I only assign
                % 10 stimulation at the moment, so anything beyond 10 is
                % artefact.
                for i = 1:nSweep
                    if size(obj.wholeData, 2) == 4 % if 4 channels are recorded, that means one channel is used for recording the trigger.
                        extStimTemp = obj.wholeData(:,3,i);
                    else
                        extStimTemp = obj.wholeData(:,2,i);
                    end
                    detectLevel = max(extStimTemp); % 20150204
                    %                 diffExtStim = diff(obj.extStim);
                    norExtStim = extStimTemp - median(extStimTemp);
                    allStimPoint = find(abs(norExtStim) > 0.7*detectLevel);
                    
                    [pStimPeaks, pStimLocs] =...
                        findpeaks(norExtStim, 'MINPEAKHEIGHT', 0.7*detectLevel,'MINPEAKDISTANCE',round(0.01*obj.SR));
                    [nStimPeaks, nStimLocs] =...
                        findpeaks(-norExtStim,'MINPEAKHEIGHT', 0.7*detectLevel, 'MINPEAKDISTANCE', round(0.01*obj.SR));
                    if isempty(stimTime)
                        if size(obj.wholeData, 2) == 4
                            % If using trigger, note that there are two
                            % triggers, one for start recording, and the
                            % other for  time stamping stimulation, so take
                            % the second trigger as the stimulation time.
                            obj.stimStartDataPoint(i) = pStimLocs(end);
                            obj.stimEndDataPoint(i) = pStimLocs(end) + 10;
                            
                        else
                            obj.stimStartDataPoint(i) = nStimLocs;
                            obj.stimEndDataPoint(i) = pStimLocs;
                        end
                        
                    else
                        obj.stimStartDataPoint(i) = round(stimTime*obj.SR);
                        obj.stimEndDataPoint(i) = round(stimTime*obj.SR) + 5;
                        % manually assign stimulation time, this is mainly
                        % for the control recording with 0 amplitude
                        % stimulation
                    end
                    obj.stimStartTime(i) = obj.stimStartDataPoint(i)/obj.SR;
                    obj.stimEndTime(i) = obj.stimEndDataPoint(i)/obj.SR;
                end
                
                obj.stimType = 'extStim';
                
            end
            
            
            
            
            %% LED (light stimulation), may have multiple stimulation in each sweep, so use cell structure for stimStartDataPoint
            if ~isempty(obj.LED)
                
                if ndims(obj.wholeData) == 2
                    diffLEDStim = diff(obj.LED);
                    allStimPoint = find(diffLEDStim > 1);
                    diffStimTime = diff([0, obj.dataTime(allStimPoint)]);
                    %                 findPulse = find(diffStimTime > 0.1);
                    % Time between each stimulation has to be at least 0.5ms long
                    obj.stimStartDataPoint = allStimPoint(diffStimTime > 0.1);
                    obj.stimStartTime = obj.dataTime(obj.stimStartDataPoint);
                    
                    allStimPoint = find(diffLEDStim < -1);
                    diffStimTime = diff([obj.dataTime(allStimPoint), obj.dataTime(end)]);
                    %                 findPulse = find(diffStimTime > 0.1);
                    
                    %             if findPulse(end) ~= length(diffStimTime)
                    obj.stimEndDataPoint = allStimPoint(diffStimTime > 0.1);
                    %             else
                    %                 pulseEndPoint = allStimPoint(findPulse(1:(end-1))- 1);
                    %                 pulseEndPoint(end) = allStimPoint(findPulse(end));
                    %             end
                    obj.stimEndTime = obj.dataTime(obj.stimEndDataPoint);
                    
                else
                    nSweep = size(obj.wholeData,3);
                    
                    for i = 1:nSweep
                        diffLEDStim = diff(obj.LED(i,:));
                        allStimPoint = find(diffLEDStim > 1);
                        diffStimTime = diff([0, obj.dataTime(allStimPoint)]);
                        %                 findPulse = find(diffStimTime > 0.1);
                        % Time between each stimulation has to be at least 0.5ms long
                        obj.stimStartDataPoint{i} = allStimPoint(diffStimTime > 0.1);
                        obj.stimStartTime{i} = obj.dataTime(obj.stimStartDataPoint{i});
                        
                        allStimPoint = find(diffLEDStim < -1);
                        diffStimTime = diff([obj.dataTime(allStimPoint), obj.dataTime(end)]);
                        %                 findPulse = find(diffStimTime > 0.1);
                        
                        %             if findPulse(end) ~= length(diffStimTime)
                        obj.stimEndDataPoint{i} = allStimPoint(diffStimTime > 0.1);
                        %             else
                        %                 pulseEndPoint = allStimPoint(findPulse(1:(end-1))- 1);
                        %                 pulseEndPoint(end) = allStimPoint(findPulse(end));
                        %             end
                        obj.stimEndTime{i} = obj.dataTime(obj.stimEndDataPoint{i});
                        
                    end
                end
                
                obj.stimType = 'LED';
            end
            
            
            obj.stimStartEndCorrection;
            
            obj.nStim = length(obj.stimStartTime);
            
            if strcmp(obj.stimType, 'extStim')|strcmp(obj.stimType, 'intStim')
                obj.extStimDissection;
                % dissection of an episodic mode sweep into: baseline, stimulation
                % and after-stimulation periods.
            end
            
            if strcmp(obj.stimType, 'LED')
                obj.lightStimDissection
                % dissection of light stimulation
            end
            
        end
        
        
        function extStimDissection(obj)
                %% For intracellular recording and episodic mode extStim (SPS)
                obj.stimDataPoint = cell(1,obj.ndim);
                obj.AStimDataPoint = cell(1,obj.ndim);
                obj.baselineDataPoint = cell(1,obj.ndim);
                for i = 1:obj.nSweeps
                    % constructing stimulation time and the max/min membrane
                    % potential
                    %                         obj.intStimAmp(i) = obj.stimInitAmp(i) + (i - 1)*obj.stimDAmp;
                    
                    %                     obj.stimTime(i,:) = [(obj.STST*(i - 1) + obj.stimStartTime(1)), (obj.STST*(i - 1) + obj.stimEndTime(1))];
                    %                     obj.stimDataPoint(i,:) = round(obj.stimTime(i,1)*obj.SR):round(obj.stimTime(i,2)*obj.SR - 1);
                    
                    
                    obj.stimDataPoint{i} = ...
                        obj.stimStartDataPoint(i):obj.stimEndDataPoint(i);
                    
                    if obj.nSweeps == 1
                        obj.stimTime{i} = obj.stimDataPoint{i}/obj.SR;
                    else
                        obj.stimTime{i} = obj.STST*(i - 1) + obj.stimDataPoint{i}/obj.SR;
                    end
                    
                    
                    
                    obj.stimMax(i) = max(obj.patch(i,obj.stimDataPoint{i}));
                    
                    %                     if obj.intStimAmp(i) > 0
                    %                         obj.stimMin(i) = median(obj.patch(i,obj.stimDataPoint{i}));
                    %                     else
                    %                         obj.stimMin(i) = mean(obj.patch(i,obj.stimDataPoint{i}));
                    %                     end
                    
                    % constructing after stimulation time and the max/min membrane
                    % potential
                    %                     obj.AStimTime(i,:) = [(obj.STST*(i - 1) + obj.stimEndTime(1)), obj.recordTime(i,end)];
                    %                     obj.AStimDataPoint(i,:) = round(obj.AStimTime(i,1)*obj.SR):round(obj.AStimTime(i,2)*obj.SR);
                    if obj.nSweeps == 1
                        obj.AStimDataPoint{i} = ...
                            obj.stimEndDataPoint(i):length(obj.patch);
                    else
                        obj.AStimDataPoint{i} = ...
                            obj.stimEndDataPoint(i):size(obj.recordTime,2);
                    end
                    
                    if obj.nSweeps == 1
                        obj.AStimTime{i} = obj.AStimDataPoint{i}/obj.SR;
                    else
                        obj.AStimTime{i} = obj.STST*(i - 1) + obj.AStimDataPoint{i}/obj.SR;
                    end
                    
                    
                    %                     obj.stimTime{i} = obj.STST*(i - 1) + obj.stimDataPoint(i,:)/obj.SR;
                    
                    
                    obj.AStimMax(i) = max(obj.patch(i,obj.AStimDataPoint{i}));
                    
                    obj.AStimMin(i) = min(obj.patch(i,obj.AStimDataPoint{i}));
                    
                    % constructing baseline time and the max/min membrane
                    % potential
                    %                     obj.baselineTime(i,:) = [obj.STST*(i - 1), obj.STST*(i - 1) + obj.stimStartTime(1)];
                    %                     obj.baselineDataPoint(i,:) = round(obj.baselineTime(i,1)*obj.SR + 1):round(obj.baselineTime(i,2)*obj.SR - 1);
                    
                    obj.baselineDataPoint{i} = ...
                        1:obj.stimStartDataPoint(i);
                    
                    if obj.nSweeps == 1
                        obj.baselineTime{i} = obj.baselineDataPoint{i}/obj.SR;
                    else
                        obj.baselineTime{i} = obj.STST*(i - 1) + obj.baselineDataPoint{i}/obj.SR;
                    end
                    
                    obj.baselineMax(i) = max(obj.patch(i,obj.baselineDataPoint{i}));
                    obj.restVm(i) = min(obj.patch(i,obj.baselineDataPoint{i}));
                    % remember that for each sweep the patch data point starts
                    % from the beginning, so just need to use the 1st
                    % baselineDataPoint to get the baseline value for each
                    % sweep
                end
        end
        
        function lightStimDissection(obj)
            
            if ndims(obj.wholeData) == 2
                
                obj.baselineDataPoint = 1:obj.stimStartDataPoint(1);
                obj.baselineTime = obj.baselineDataPoint/obj.SR;
                for i = 1:obj.nSweeps
                    % obj.baselineDataPoint = round(obj.baselineTime(i,1)*obj.SR + 1):round(obj.baselineTime(i,2)*obj.SR - 1);
                    % obj.baselineTime = [0, obj.stimStartTime(1)];
                    
                    % obj.LEDONTime(i,:) = [obj.stimEndTime(1)];
                    % obj.LEDONDataPoint(i,:) = round(obj.stimTime(i,1)*obj.SR):round(obj.stimTime(i,2)*obj.SR - 1);
                    if obj.nSweeps ~= 1
                        obj.LEDONDataPoint{i,:} = obj.stimStartDataPoint(i,:):obj.stimEndDataPoint(i,:);
                        obj.LEDONTime{i,:} = obj.LEDONDataPoint{i,:}/obj.SR;
                    
                    else
                        obj.LEDONDataPoint{i,1} = obj.stimStartDataPoint(i,:):obj.stimEndDataPoint(i,:);
                        obj.LEDONTime{i,1} = obj.LEDONDataPoint{i,:}/obj.SR;
                    end
                    
                    
                    %                         obj.LEDOFFTime(i,:) = [(obj.STST*(i - 1) + obj.stimEndTime(1)), obj.recordTime(i,end)];
                    %                         obj.LEDOFFDataPoint(i,:) = round(obj.AStimTime(i,1)*obj.SR):round(obj.AStimTime(i,2)*obj.SR);
                    if obj.nSweeps ~= 1
                    for j = 1:length(obj.LEDONDataPoint{i,:})
                        if i ~= length(obj.stimStartTime)
                            obj.LEDOFFDataPoint{i,j} = obj.stimEndDataPoint(i,j):obj.stimStartDataPoint(i,j + 1);
                        else
                            obj.LEDOFFDataPoint{i,j} = obj.stimEndDataPoint(i,j):(obj.stimEndDataPoint(i,j) + 0.5*obj.SR);
                        end
                        obj.LEDOFFTime{i,j} = obj.LEDOFFDataPoint{i,j}/obj.SR;
                        
                    end
                    else
                        if i ~= length(obj.stimStartTime)
                            obj.LEDOFFDataPoint{i,1} = obj.stimEndDataPoint(i,1):obj.stimStartDataPoint(i,1 + 1);
                        else
                            obj.LEDOFFDataPoint{i,1} = obj.stimEndDataPoint(i,1):(obj.stimEndDataPoint(i,1) + 0.5*obj.SR);
                        end
                    end
                end
                
            elseif ndims(obj.wholeData) == 3
                for i = 1:size(obj.stimStartTime,2)
                    %                         obj.baselineDataPoint = round(obj.baselineTime(i,1)*obj.SR + 1):round(obj.baselineTime(i,2)*obj.SR - 1);
                    %                         obj.baselineTime = [0, obj.stimStartTime(1)];
                    
                    % obj.LEDONTime(i,:) = [obj.stimEndTime(1)];
                    % obj.LEDONDataPoint(i,:) = round(obj.stimTime(i,1)*obj.SR):round(obj.stimTime(i,2)*obj.SR - 1);
                    nLED = length(obj.stimStartTime{i});
                    for j = 1:nLED
                        obj.LEDONDataPoint{i,j} = obj.stimStartDataPoint{i}(j):obj.stimEndDataPoint{i}(j);
                        obj.LEDONTime{i,j} = obj.LEDONDataPoint{i,j}/obj.SR;
                        
                        %                         obj.LEDOFFTime(i,:) = [(obj.STST*(i - 1) + obj.stimEndTime(1)), obj.recordTime(i,end)];
                        %                         obj.LEDOFFDataPoint(i,:) = round(obj.AStimTime(i,1)*obj.SR):round(obj.AStimTime(i,2)*obj.SR);
                        if j ~= length(obj.stimStartTime{i})
                            obj.LEDOFFDataPoint{i,j} = obj.stimEndDataPoint{i}(j):obj.stimStartDataPoint{i}(j + 1);
                        else
                            obj.LEDOFFDataPoint{i,j} = obj.stimEndDataPoint{i}(j):(obj.stimEndDataPoint{i}(j) + 0.5*obj.SR);
                        end
                        obj.LEDOFFTime{i,j} = obj.LEDOFFDataPoint{i,j}/obj.SR;
                        
                        
                        
                        
                    end
                    obj.baselineDataPoint{i} = 1:obj.stimStartDataPoint{i}(1);
                    obj.baselineTime{i} = obj.recordTime(1,obj.baselineDataPoint{i});
                    obj.AStimDataPoint{i} = obj.LEDOFFDataPoint{i,nLED}(end):obj.nSamples;
                    obj.AStimTime{i} = obj.recordTime(1,obj.AStimDataPoint{i});
                    
                end
            end
        end
        
        
        
        
        
        
        function stimStartEndCorrection(obj)
            % To correct the problem
            % 1. where the stimulation started before recording started
            % (where the trigger will be on high before recording).
            % This is done by setting the stimulation start time at dataPoint 1 as the starting point
            % 2. Where the stimulation ended beyond the end of recording
            % (the trigger will remain on high till the end and no
            % stimulation end will be recorded). This is corrected by
            % setting the recording end point as stimulation end time.
            
            
            if (~isempty(obj.stimEndTime))&(~isempty(obj.stimStartTime))
                
                if iscell(obj.stimEndTime)
                    % Specifically for LED stimulation with episodic mode
                    for i = 1:size(obj.stimEndTime,2)
                        
                        if isempty(obj.stimStartTime{i}(1))
                            obj.stimStartDataPoint{i} = 0;
                            obj.stimStartTime{i} = 0;
                            continue
                        end
                        if obj.stimEndTime{i}(1) < obj.stimStartTime{i}(1)
                            obj.stimStartDataPoint{i} = [1, obj.stimStartDataPoint{i}];
                            obj.stimStartTime{i} = [obj.dataTime(1), obj.stimStartTime{i}];
                        end
                        
                        if obj.stimEndTime{i}(end) < obj.stimStartTime{i}(end)
                            obj.stimEndDataPoint{i} = [obj.stimEndDataPoint{i}, obj.nSamples];
                            obj.stimEndTime{i} = [obj.stimEndTime{i}, obj.dataTime(end)];
                        end
                        
                    end
                else
                    % LED stimulation with gap-free mode or either
                    % modes for extStim
                    
                    if obj.stimEndTime(1) < obj.stimStartTime(1)
                        obj.stimStartDataPoint = [1, obj.stimStartDataPoint];
                        obj.stimStartTime = [obj.dataTime(1), obj.stimStartTime];
                    end
                    
                    if obj.stimEndTime(end) < obj.stimStartTime(end)
                        obj.stimEndDataPoint = [obj.stimEndDataPoint, obj.nSamples];
                        obj.stimEndTime = [obj.stimEndTime, obj.dataTime(end)];
                    end
                end
            end
        end
        
        function corrExtStim(obj)
            % correct the current clamp voltage of the cell by subtracting
            % the extracellular stimulation voltage recorded by the other
            % digitiser input
            
            % Method
            % 1. make sections of patch and extStim both go to baseline
            % (minus 1st entry value)
            % 2. Find patch peaks (because it is not exactly the same as
            % the extStim peak)
            % 3. adjust the maximal point time shift of cell patch voltage
            % 4. compare the peak values between the patch and extStim,
            % times a multiplication factor to the extStim so that the peak
            % values match the patch clamp voltage ones
            % 5. Subtract the adjusted extStim from the patch voltage to
            % get the corrected voltage
            
            obj.patchExtStimCorrected = obj.patch;
            extentTime = 0.001;
            for i = 1:length(obj.stimStartDataPoint)
                dataPtForSub =...
                    (obj.stimStartDataPoint(i) - extentTime*obj.SR):(obj.stimEndDataPoint(i) + extentTime*obj.SR);
                dataPtForSub(dataPtForSub < 0) = [];
                dataPtForSub(dataPtForSub > length(obj.dataTime)) = [];
                dataPtForSub = round(dataPtForSub);
                
                tempPatch = obj.patch;
                baselineCorrVal = tempPatch(dataPtForSub(1));
                tempPatch = tempPatch - baselineCorrVal;
                % norming the patch data with the first data point for
                % adjustment
                
                tempExtStim = obj.extStim;
                tempExtStim = tempExtStim - tempExtStim(dataPtForSub(1));
                %                 baselineDiff = tempPatch(dataPtForSub(1)) - tempExtStim(dataPtForSub(1));
                
                
                %                 tempExtStim = tempExtStim + baselineDiff;
                
                [nPatchPeaks, nPatchLocs] = findpeaks(-tempPatch(dataPtForSub),'MINPEAKDISTANCE',round(0.001*obj.SR));
                [pPatchPeaks, pPatchLocs] = findpeaks(tempPatch(dataPtForSub),'MINPEAKDISTANCE',round(0.001*obj.SR));
                
                multiF =...
                    tempPatch((dataPtForSub(1) - 1) + nPatchLocs(1))/tempExtStim(obj.stimStartDataPoint(i));
                
                adjExtStim = tempExtStim(dataPtForSub)*multiF;
                
                maxShift = (dataPtForSub(1) - 1) + nPatchLocs(1) - obj.stimStartDataPoint(i);
                % There is a time delay of in the patch clamp voltage
                % maximum compared to extStim max, so have to make the
                % shift.
                
                obj.patchExtStimCorrected(dataPtForSub + maxShift) =...
                    tempPatch(dataPtForSub + maxShift) - adjExtStim;
                obj.patchExtStimCorrected(dataPtForSub + maxShift) =...
                    obj.patchExtStimCorrected(dataPtForSub + maxShift) + baselineCorrVal;
            end
        end
        
        function getExtStimSpike(obj)
            if  isempty(obj.patchExtStimCorrected)
                obj.corrExtStim;
            end
            
            tempPatch = obj.patchExtStimCorrected;
            spikeThresh = -20;
            
            extentTime = 0.2;
            for i = 1:length(obj.stimStartDataPoint)
                dataPtForSub =...
                    obj.stimStartDataPoint(i):(obj.stimEndDataPoint(i) + extentTime*obj.SR);
                dataPtForSub(dataPtForSub < 0) = [];
                dataPtForSub(dataPtForSub > length(obj.dataTime)) = [];
                dataPtForSub = round(dataPtForSub);
                
                [spikePeaks, spikeLocs] = findpeaks(tempPatch(dataPtForSub),'MINPEAKHEIGHT', spikeThresh, 'MINPEAKDISTANCE',round(0.001*obj.SR));
                %                spikeI = find(tempPatch(dataPtForSub) > spikeThresh) + dataPtForSub(1) - 1;
                if ~isempty(spikeLocs)
                    obj.extStimSpikeDataPt((end + 1):(end + length(spikeLocs))) = dataPtForSub(1) - 1 + spikeLocs;
                end
            end
            
        end
        
        
%         function getSPSSpike(obj, cutOffThreshPlus, minGap)
%             obj.spikeletCutOff = -50;
%             
%             if nargin > 1
%                 obj.spikeCutOff = cutOffThreshPlus;
%                 % how much above baseline
%             else
%                 obj.spikeCutOff = 20;
%             end
%             
%             if ~isempty(obj.patch)
%                 
%             
%             
%         end
        
%         function getHFSSpike(obj, cutOffThreshPlus, shortMinGap, longMinGap, HFSThreshPercent)
        function getHFSSpike(obj, varargin)
            % for HFS spike detection
            
%             obj.spikeletCutOff = -50;
%             spikeletCutOff = 10;
%             
%             highSpikeCutOff = 40;

            spikeCutOff = 10;
            
            sponSpikeCutoff = [];
            
            minGap = 0.001;
            
            APThreshFac = 0.2;
            
            smoothN = 0;
            
            analyseStimN = 1:obj.nSweeps;
            
            extendStimDur = 0;
            
            % 4 parameters AP finding system
            % high spike cut off is associated with short minimum gap
            % between peaks. This is to dreduce the chance of finding false
            % peaks around a true peak: by making sure that another peak
            % near a true peak has to have high amplitde.
            
            % low spike cut off is to be used with long minimum peak cut so
            % to find APs with smaller amplitudes (mainly with high
            % frequency down curve effect).
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'spikeCutOff')
                    spikeCutOff = varargin{i + 1};
                end
                
                if strcmp(varargin{i}, 'sponSpikeCutOff')
                    sponSpikeCutoff = varargin{i + 1};
                end
                
                
                if strcmp(varargin{i}, 'minGap')
                    minGap = varargin{i + 1};
                end

                if strcmp(varargin{i}, 'APThreshFac')
                    APThreshFac = varargin{i + 1};
                end
                              
                if strcmp(varargin{i}, 'stimN')
                    analyseStimN = varargin{i + 1};
                end
                
                if strcmp(varargin{i}, 'smoothN')
                    smoothN = varargin{i + 1};
                    % To smooth the filtered data (get better result)
                end
                
                if strcmp(varargin{i}, 'extendStimDur')
                   % This is to extend stimulation period for analysis for low frequency stimulation, for example, 50 Hz
                   extendStimDur = varargin{i + 1}*obj.SR; % time (is s) times the sampling rate to get data point;
                end
                
            end
            
            if nargin > 1
                obj.spikeCutOff = spikeCutOff;
                % how much above baseline
            else
                obj.spikeCutOff = 20;
            end
            
            
            noExtStimSweep = [];
            
            if ~isempty(obj.patch)
                %                 totalSpikeTemp = [];
                if (~isempty(obj.intStim))|(~isempty(obj.extStim))
                    %% Gap-free mode
                    if obj.ndim == 2
                        
                        HFSStimDiff = diff(obj.extStimStartPt);
                        % This is for gap free mode HFS. Using the diff function to
                        % find times at which HFS stim started
                        findHFSStimDiff = find(HFSStimDiff > 0.5*obj.SR);
                        % if the data point difference is bigger than 0.5s
                        % (0.5*obj.SR), then it is considered as a different
                        % stimulation
                        
                        % need to write code for episodic mode
                        
                        HFSStartStim = [obj.extStimStartPt(1), obj.extStimStartPt(findHFSStimDiff + 1)];
                        HFSEndStim = [obj.extStimStartPt(findHFSStimDiff), obj.extStimStartPt(end)];
                        
                        
                        if isempty(obj.HFSfilteredSig)
                            tempPatch = obj.patch(:);
                        else
                            tempPatch = obj.HFSfilteredSig(:);
                        end
                        
                        
                        
                        
                        %% Baseline spikes
                        obj.baselineDataPoint{1} = 1:HFSStartStim(1);
                        
                        baselineBase = median(tempPatch(1:HFSStartStim(1)));
                        % For 20151002 experiment where HFS was
                        % recorded with gap free mode
                        baselineCutOff = baselineBase + obj.spikeCutOff;
                        
                        if isempty(sponSpikeCutoff)
                            if max(tempPatch) > 0
                                sponSpikeCutoff = 40;
                            else
                                sponSpikeCutoff = 0.8*(max(tempPatch) - baselineBase);
                            end
                        end
                        
                        [baseSpikeAmp, baseSpikeDataPoint] = ...
                            findpeaks(tempPatch(1:HFSStartStim(1)), 'MINPEAKHEIGHT', sponSpikeCutoff + baselineBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                        % use the 1st column of baselineDataPoint because
                        % patch is in sweep format and the datapoint starts
                        % from 1 for each sweep
                       
%                         adjustBaseDataPoint = baseSpikeDataPoint;
%                         
%                         [baseSpikeletAmp, baseSpikeletDataPoint] = obj.sortSpike(tempPatch, baseSpikeAmp, adjustBaseDataPoint,obj.spikeCutOff);
%                         
%                         
%                         obj.baseSpikelet{1,1} = obj.dataTime(baseSpikeletDataPoint);
%                         obj.baseSpikelet{2,1} = baseSpikeletAmp;
%                         obj.baseSpikelet{3,1} = baseSpikeletDataPoint;
%                         
%                         baseSpikeDataPoint(baseSpikeAmp < baselineCutOff) = [];
%                         baseSpikeAmp(baseSpikeAmp < baselineCutOff) = [];
                        
                        blArtefactEntry = find(baseSpikeAmp > 10);
                        if ~isempty(blArtefactEntry)
                            baseSpikeDataPoint(blArtefactEntry) = [];
                            baseSpikeAmp(blArtefactEntry) = [];
                        end
                        
                        
                        obj.baseSpike{1,1} = baseSpikeDataPoint/obj.SR;
                        obj.baseSpike{2,1} = baseSpikeAmp;
                        obj.baseSpike{3,1} = baseSpikeDataPoint;
                        
                        
                        %% After stimulation  spikes
                        obj.AStimDataPoint{1} = HFSEndStim(end):obj.nSamples;
                        
                        AStimBase = median(tempPatch(HFSEndStim(end):end));
                        AStimCutOff = AStimBase + obj.spikeCutOff;
                        
                        [AStimSpikeAmp, AStimDataPoint] = ...
                            findpeaks(tempPatch(HFSEndStim(end):end), 'MINPEAKHEIGHT', sponSpikeCutoff +  AStimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                        
%                         adjustAStimDataPoint = HFSEndStim(end) + AStimDataPoint;
%                         [AStimSpikeletAmp, AStimSpikeletDataPoint] = obj.sortSpike(tempPatch, AStimSpikeAmp, adjustAStimDataPoint,obj.spikeCutOff);
%                         
%                         obj.AStimSpikelet{1,1} = obj.dataTime(AStimSpikeletDataPoint);
%                         obj.AStimSpikelet{2,1} = AStimSpikeletAmp;
%                         obj.AStimSpikelet{3,1} = AStimSpikeletDataPoint;
%                         
%                         AStimDataPoint(AStimSpikeAmp < AStimCutOff) = [];
%                         AStimSpikeAmp(AStimSpikeAmp < AStimCutOff) = [];
                        
                        artefactEntry = find(AStimSpikeAmp > 10);
                        if ~isempty(artefactEntry)
                            AStimDataPoint(artefactEntry) = [];
                            AStimSpikeAmp(artefactEntry) = [];
                        end
                        
                        
                        obj.AStimSpike{1,1} = (HFSEndStim(end) + AStimDataPoint - 1)/obj.SR;
                        obj.AStimSpike{2,1} = AStimSpikeAmp;
                        obj.AStimSpike{3,1} = HFSEndStim(end) + AStimDataPoint - 1;
                        
                        for HFSStimI = 1:length(HFSStartStim)
                            
                            obj.HFSDataPt{HFSStimI} = HFSStartStim(HFSStimI):HFSEndStim(HFSStimI);
                            obj.HFSTime{HFSStimI} = obj.HFSDataPt{HFSStimI}/obj.SR;
                            
                            
                            
                            
                            
                            
                            if ~isempty(obj.extStimStartPt)
                                %                             if length(HFSStimN) > round(0.001*obj.SR)
                                
                                stimBase = median(tempPatch(obj.HFSDataPt{HFSStimI}));
                                stimCutOff = stimBase + obj.spikeCutOff;
                                
                                timeDelay = 0.005;
                                timeDelayDataPt = timeDelay*obj.SR;
                                
                                adpatedHFSDataPt = (obj.HFSDataPt{HFSStimI}(1) + timeDelayDataPt):obj.HFSDataPt{HFSStimI}(end);
                                
                                maxAdaptedAPAmp = max(tempPatch(obj.HFSDataPt{HFSStimI}));
                                
                                if stimBase < -50
                                    % When the baseline is not too far from
                                    % -60, it means the HFS amplitude is
                                    % not too high yet
                                    
                                
                                [adaptedAPAmp, ~] = ...
                                    findpeaks(tempPatch(adpatedHFSDataPt), 'MINPEAKHEIGHT', spikeCutOff + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                                else
                                    % When the baseline is far away from
                                    % -60, it means the HFS amplitude is
                                    % strong and AP should be low
                                    [adaptedAPAmp, ~] = ...
                                    findpeaks(tempPatch(adpatedHFSDataPt), 'MINPEAKHEIGHT', lowSpikeCutOff + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                                end
                                % To find the 1st APs that are not in the
                                % initial stage of HFS.
                                
                                [stimSpikeAmp, stimSpikeDataPoint] = ...
                                    findpeaks(tempPatch(obj.HFSDataPt{HFSStimI}), 'MINPEAKHEIGHT', spikeletCutOff + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                                
                                
                                
                                
                                
                                adjustStimDataPoint = obj.HFSDataPt{HFSStimI}(1) + stimSpikeDataPoint;
                                
                                [stimSpikeletAmp, stimSpikeletDataPoint] = obj.sortSpike(tempPatch, stimSpikeAmp, adjustStimDataPoint,obj.spikeCutOff);
                                
                                obj.HFSSpikelet{1,HFSStimI} = obj.dataTime(stimSpikeletDataPoint);
                                % Gap-free mode = dataTime
                                % episodic mode: both dataTime and recordTime
                                % are available
                                obj.HFSSpikelet{2,HFSStimI} = stimSpikeletAmp;
                                obj.HFSSpikelet{3,HFSStimI} = stimSpikeletDataPoint;
                                
                                if ~isempty(adaptedAPAmp)
                                    maxAdaptedAPAmp = max(adaptedAPAmp);
                                    if (maxAdaptedAPAmp - stimBase) > 60
                                        updateStimCutOff = -10;
                                        [stimSpikeAmp, stimSpikeDataPoint] = ...
                                            findpeaks(tempPatch(obj.HFSDataPt{HFSStimI}),...
                                            'MINPEAKHEIGHT', (updateStimCutOff - stimBase)*highAPThreshFac + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
%                                         longThreshFac = 0.7;
                                    elseif (maxAdaptedAPAmp - stimBase) > 40
                                        updateStimCutOff = maxAdaptedAPAmp;
                                        [stimSpikeAmp, stimSpikeDataPoint] = ...
                                            findpeaks(tempPatch(obj.HFSDataPt{HFSStimI}),...
                                            'MINPEAKHEIGHT', (updateStimCutOff - stimBase)*highAPThreshFac + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
%                                         longThreshFac = 0.7;
                                    else
                                        updateStimCutOff = maxAdaptedAPAmp;
                                        [stimSpikeAmp, stimSpikeDataPoint] = ...
                                            findpeaks(tempPatch(obj.HFSDataPt{HFSStimI}),...
                                            'MINPEAKHEIGHT', (updateStimCutOff - stimBase)*lowAPThreshFac + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
%                                         longThreshFac = 0.2;
                                    end
                                else
                                    updateStimCutOff = -10;
                                    [stimSpikeAmp, stimSpikeDataPoint] = ...
                                            findpeaks(tempPatch(obj.HFSDataPt{HFSStimI}),...
                                            'MINPEAKHEIGHT', (updateStimCutOff - stimBase)*highAPThreshFac + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                                    
                                end
                                


                                  stimSpikeAmp = stimSpikeAmp';
                                  stimSpikeDataPoint = stimSpikeDataPoint';
                                

                                
                                
                                obj.HFSSpike{1,HFSStimI} = (obj.HFSDataPt{HFSStimI}(1) + stimSpikeDataPoint)/obj.SR;
                                % Time point
                                obj.HFSSpike{2,HFSStimI} = stimSpikeAmp;
                                % spike amplitude
                                obj.HFSSpike{3,HFSStimI} = obj.HFSDataPt{HFSStimI}(1) + stimSpikeDataPoint;
                                % data point
                                
                                %                                 obj.HFSSpikeNumber(HFSStimN) = length(stimSpikeDataPoint);
                                
                                %                                 obj.generalSpikeRate(HFSStimN) = length(stimSpikeDataPoint)/(HFSTime(end) - HFSTime(1));
                                % generalSpikeRate is simply the number of
                                % spikes divided by the stimulation time
                                
                                %                                 APTimeDiff = diff(obj.stimSpike{1,HFSStimN});
                                %
                                %                                 obj.indSpikeRate{HFSStimN} = 1./APTimeDiff;
                                %                                 % spike rate obtained by looking at the time
                                %                                 % gap between each AP. Divide 1 by this time to
                                %                                 % get the rate.
                                
                                
                                if HFSStimI ~= length(HFSStartStim)
                                    
                                    obj.HFSGapDataPoint{HFSStimI} = HFSEndStim(HFSStimI):HFSStartStim(HFSStimI + 1);
                                    obj.HFSGapTime{HFSStimI} = obj.HFSGapDataPoint{HFSStimI}/obj.SR;
                                    
                                    stimGapBase = median(tempPatch(obj.HFSGapDataPoint{HFSStimI}));
                                    stimGapCutOff = stimGapBase + obj.spikeCutOff;
                                    
                                    [stimGapSpikeAmp, stimGapSpikeDataPoint] = ...
                                        findpeaks(tempPatch(obj.HFSGapDataPoint{HFSStimI}), 'MINPEAKHEIGHT', spikeletCutOff + stimGapBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                                                                        
                                    
                                    
                                    adjustStimGapDataPoint = obj.HFSGapDataPoint{HFSStimI}(1) + stimGapSpikeDataPoint;
                                    
                                    [stimGapSpikeletAmp, stimGapSpikeletDataPoint] = obj.sortSpike(tempPatch, stimGapSpikeAmp, adjustStimGapDataPoint,obj.spikeCutOff);
                                    
                                    obj.HFSGapSpikelet{1,HFSStimI} = obj.dataTime(stimGapSpikeletDataPoint);
                                    % Gap-free mode = dataTime
                                    % episodic mode: both dataTime and recordTime
                                    % are available
                                    obj.HFSGapSpikelet{2,HFSStimI} = stimGapSpikeletAmp;
                                    obj.HFSGapSpikelet{3,HFSStimI} = stimGapSpikeletDataPoint;
                                    
                                    stimGapSpikeDataPoint(stimGapSpikeAmp < stimCutOff) = [];
                                    stimGapSpikeAmp(stimGapSpikeAmp < stimCutOff) = [];
                                    
                                    
                                    obj.HFSGapSpike{1,HFSStimI} = (obj.HFSGapDataPoint{HFSStimI}(1) + stimGapSpikeDataPoint)/obj.SR;
                                    % Time point
                                    obj.HFSGapSpike{2,HFSStimI} = stimGapSpikeAmp;
                                    % spike amplitude
                                    obj.HFSGapSpike{3,HFSStimI} = obj.HFSGapDataPoint{HFSStimI}(1) + stimGapSpikeDataPoint;
                                    %
                                    %
                                end
                                
                            end
                            
                        end

                        
                    else % episodic recording mode
                        for HFSStimI = 1:length(analyseStimN)
                            
                            if obj.extStimStartPt(analyseStimN(HFSStimI)) == 0
                                obj.extStimStartPt(analyseStimN(HFSStimI)) = obj.extStimStartPt(2);  
                                obj.extStimEndPt(analyseStimN(HFSStimI)) = obj.extStimEndPt(2);
                                noExtStimSweep = [noExtStimSweep, analyseStimN(HFSStimI)];
%                                 
                            end
                            
                            tempPatchOrg = obj.patch(analyseStimN(HFSStimI),:);
                            % original data
                            
                            if isempty(obj.HFSfilteredSig)
                                tempPatch = obj.patch(analyseStimN(HFSStimI),:);
                            else
                                tempPatch = obj.HFSfilteredSig(analyseStimN(HFSStimI),:);
                            end
                            %% Baseline spikes
                            
                            obj.baselineDataPoint{analyseStimN(HFSStimI)} = 1:obj.extStimStartPt(analyseStimN(HFSStimI));
                            
                            baselineBase = median(tempPatch(obj.baselineDataPoint{analyseStimN(HFSStimI)}));
                            maxBaselineDiff = max(tempPatch(obj.baselineDataPoint{analyseStimN(HFSStimI)})) - baselineBase;
                            
                            artefactAmp = 10;
                            
                            if isempty(sponSpikeCutoff)
                                if max(tempPatch(obj.baselineDataPoint{analyseStimN(HFSStimI)})) > 0
                                    artefactAmp = 30;
                                    if max(tempPatch(obj.baselineDataPoint{analyseStimN(HFSStimI)})) > 20 % in case of artefact
                                        sponBaseSpikeCutoff = 20;
                                    else
                                        sponBaseSpikeCutoff = 0.5*maxBaselineDiff;
                                    end
                                elseif maxBaselineDiff < 20
                                    sponBaseSpikeCutoff = 20;
                                else
                                    sponBaseSpikeCutoff = 0.5*maxBaselineDiff;
                                end
                            end
                            
                            baselineCutOff = baselineBase + sponBaseSpikeCutoff;
                            

                            
                            [baseSpikeAmp, baseSpikeDataPoint] = ...
                                findpeaks(tempPatch(obj.baselineDataPoint{analyseStimN(HFSStimI)}), 'MINPEAKHEIGHT', baselineCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                            % use the 1st column of baselineDataPoint because
                            % patch is in sweep format and the datapoint starts
                            % from 1 for each sweep
                            
%                             adjustBaseDataPoint = obj.baselineDataPoint{analyseStimN(HFSStimI)}(1) + baseSpikeDataPoint;
%                             
%                             [baseSpikeletAmp, baseSpikeletDataPoint] = obj.sortSpike(tempPatch, baseSpikeAmp, adjustBaseDataPoint,obj.spikeCutOff);
%                             
%                             
%                             obj.baseSpikelet{1,analyseStimN(HFSStimI)} = obj.dataTime(baseSpikeletDataPoint);
%                             obj.baseSpikelet{2,analyseStimN(HFSStimI)} = baseSpikeletAmp;
%                             obj.baseSpikelet{3,analyseStimN(HFSStimI)} = baseSpikeletDataPoint;
%                             
%                             baseSpikeDataPoint(baseSpikeAmp < baselineCutOff) = [];
%                             baseSpikeAmp(baseSpikeAmp < baselineCutOff) = [];
                            
                            blArtefactEntry = find(baseSpikeAmp > artefactAmp);
                            if ~isempty(blArtefactEntry)
                                baseSpikeDataPoint(blArtefactEntry) = [];
                                baseSpikeAmp(blArtefactEntry) = [];
                            end
                            
                            
                            obj.baseSpike{1,analyseStimN(HFSStimI)} = (obj.baselineDataPoint{analyseStimN(HFSStimI)}(1) + baseSpikeDataPoint - 1)/obj.SR;
                            obj.baseSpike{2,analyseStimN(HFSStimI)} = baseSpikeAmp;
                            obj.baseSpike{3,analyseStimN(HFSStimI)} = obj.baselineDataPoint{analyseStimN(HFSStimI)}(1) + baseSpikeDataPoint - 1;

                            %% After stimulation  spikes
                            
                            obj.AStimDataPoint{analyseStimN(HFSStimI)} = obj.extStimEndPt(analyseStimN(HFSStimI)):obj.nSamples;
                            
                            AStimBase = median(tempPatch(obj.AStimDataPoint{analyseStimN(HFSStimI)}));
                            maxAstimDiff = max(tempPatch(obj.AStimDataPoint{analyseStimN(HFSStimI)})) - AStimBase;
                            
                            artefactAmp = 10;
                            
                            if isempty(sponSpikeCutoff)
                                if max(tempPatch(obj.AStimDataPoint{analyseStimN(HFSStimI)})) > 0
                                    artefactAmp = 30;
                                    if max(tempPatch(obj.AStimDataPoint{analyseStimN(HFSStimI)})) > 20 % in case of artefact
                                        sponAStimSpikeCutoff = 20;
                                    else
                                        sponAStimSpikeCutoff = 0.5*maxAstimDiff;
                                    end
                                elseif maxAstimDiff < 30
                                    sponAStimSpikeCutoff = 20;
                                else
                                    sponAStimSpikeCutoff = 0.5*maxAstimDiff;
                                end
                            end
                            
                            AStimCutOff = AStimBase + sponAStimSpikeCutoff;
                            
                            AStimDataPoint = obj.AStimDataPoint{analyseStimN(HFSStimI)};
%                             AStimDataPoint(1:250) = [];
                            
                            
                            [AStimSpikeAmp, AStimSpikeDataPoint] = ...
                                findpeaks(tempPatch(AStimDataPoint), 'MINPEAKHEIGHT', AStimCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                            
%                             AStimBase = median(tempPatch(obj.AStimDataPoint{analyseStimN(HFSStimI)}));
%                             AStimCutOff = AStimBase + obj.spikeCutOff;                            
%                             
%                             adjustAStimDataPoint = obj.AStimDataPoint{analyseStimN(HFSStimI)}(1) + AStimDataPoint;
%                             [AStimSpikeletAmp, AStimSpikeletDataPoint] = obj.sortSpike(tempPatch, AStimSpikeAmp, adjustAStimDataPoint,obj.spikeCutOff);
%                             
%                             obj.AStimSpikelet{1,analyseStimN(HFSStimI)} = obj.recordTime(analyseStimN(HFSStimI),AStimSpikeletDataPoint);
%                             obj.AStimSpikelet{2,analyseStimN(HFSStimI)} = AStimSpikeletAmp;
%                             obj.AStimSpikelet{3,analyseStimN(HFSStimI)} = AStimSpikeletDataPoint;
%                             
%                             AStimDataPoint(AStimSpikeAmp < AStimCutOff) = [];
%                             AStimSpikeAmp(AStimSpikeAmp < AStimCutOff) = [];

                            AStimSpikeAmp(AStimSpikeDataPoint < 250) = [];
                            AStimSpikeDataPoint(AStimSpikeDataPoint < 250) = [];
                            
                            % To rid of spikes within the period just after
                            % stimulation since the stim box has some
                            % residual noise straight after stimulation
                            
                            
                            artefactEntry = find(AStimSpikeAmp > artefactAmp);
                            if ~isempty(artefactEntry)
                                AStimSpikeDataPoint(artefactEntry) = [];
                                AStimSpikeAmp(artefactEntry) = [];
                            end
                            
                            
                            obj.AStimSpike{1,analyseStimN(HFSStimI)} = (obj.AStimDataPoint{analyseStimN(HFSStimI)}(1) + AStimSpikeDataPoint - 1)/obj.SR;
                            obj.AStimSpike{2,analyseStimN(HFSStimI)} = AStimSpikeAmp;
                            obj.AStimSpike{3,analyseStimN(HFSStimI)} = obj.AStimDataPoint{analyseStimN(HFSStimI)}(1) + AStimSpikeDataPoint - 1;
                                
                            
                            %% HFS spikes
                            
                            if ~isempty(obj.extStimStartPt)
                                obj.HFSDataPt{analyseStimN(HFSStimI)} = obj.extStimStartPt(analyseStimN(HFSStimI)):obj.extStimEndPt(analyseStimN(HFSStimI)) + extendStimDur;
                                obj.HFSTime{analyseStimN(HFSStimI)} = obj.HFSDataPt{analyseStimN(HFSStimI)}/obj.SR;
                                %                             if length(HFSStimN) > round(0.001*obj.SR)
                                
                                stimBase = median(tempPatch(obj.HFSDataPt{analyseStimN(HFSStimI)}));
                                stimCutOff = stimBase + obj.spikeCutOff;
                                
                                timeDelay = 0.005;
                                % Add time delay to find an AP that is away
                                % from the initial part of the stimulation
                                timeDelayDataPt = timeDelay*obj.SR;
                                
                                adpatedHFSDataPt = (obj.HFSDataPt{analyseStimN(HFSStimI)}(1) + timeDelayDataPt):obj.HFSDataPt{analyseStimN(HFSStimI)}(end);
                                
                                maxAdaptedAPAmp = max(tempPatch(obj.HFSDataPt{analyseStimN(HFSStimI)}));
                                
%                                 if stimBase < -50
                                    % When the baseline is not too far from
                                    % -60, it means the HFS amplitude is
                                    % not too high yet
                                    
                                
                                [~, ~] = ...
                                    findpeaks(tempPatch(adpatedHFSDataPt), 'MINPEAKHEIGHT', spikeCutOff + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
%                                 else
%                                     % When the baseline is far away from
%                                     % -60, it means the HFS amplitude is
%                                     % strong and AP should be low
%                                     [adaptedAPAmp, ~] = ...
%                                     findpeaks(tempPatch(adpatedHFSDataPt), 'MINPEAKHEIGHT', lowSpikeCutOff + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
%                                 end
                                % To find the 1st APs that are not in the
                                % initial stage of HFS.
                                
                                [stimSpikeAmp, stimSpikeDataPoint] = ...
                                    findpeaks(tempPatch(obj.HFSDataPt{analyseStimN(HFSStimI)}), 'MINPEAKHEIGHT', spikeCutOff + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                                
                                adaptedSpikeAmp = median(stimSpikeAmp);
                                
%                                 adjustStimDataPoint = obj.HFSDataPt{analyseStimN(HFSStimI)}(1) + stimSpikeDataPoint;
%                                 
%                                 [stimSpikeletAmp, stimSpikeletDataPoint] = obj.sortSpike(tempPatch, stimSpikeAmp, adjustStimDataPoint,obj.spikeCutOff);
%                                 
%                                 obj.HFSSpikelet{1,analyseStimN(HFSStimI)} = obj.dataTime(stimSpikeletDataPoint);
%                                 % Gap-free mode = dataTime
%                                 % episodic mode: both dataTime and recordTime
%                                 % are available
%                                 obj.HFSSpikelet{2,analyseStimN(HFSStimI)} = stimSpikeletAmp;
%                                 obj.HFSSpikelet{3,analyseStimN(HFSStimI)} = stimSpikeletDataPoint;
                                
%                                 if ~isempty(adaptedAPAmp)
%                                     maxAdaptedAPAmp = max(adaptedAPAmp);
%                                     if (maxAdaptedAPAmp - stimBase) > 60
%                                         updateStimCutOff = -10;
%                                         [stimHighSpikeAmp, stimHighSpikeDataPoint] = ...
%                                             findpeaks(tempPatch(obj.HFSDataPt{analyseStimN(HFSStimI)}),...
%                                             'MINPEAKHEIGHT', (updateStimCutOff - stimBase)*highAPThreshFac + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
% %                                         longThreshFac = 0.7;
%                                     elseif (maxAdaptedAPAmp - stimBase) > 40
%                                         updateStimCutOff = maxAdaptedAPAmp;
%                                         [stimHighSpikeAmp, stimHighSpikeDataPoint] = ...
%                                             findpeaks(tempPatch(obj.HFSDataPt{analyseStimN(HFSStimI)}),...
%                                             'MINPEAKHEIGHT', (updateStimCutOff - stimBase)*highAPThreshFac + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
% %                                         longThreshFac = 0.7;
%                                     else
%                                         updateStimCutOff = maxAdaptedAPAmp;
%                                         [stimHighSpikeAmp, stimHighSpikeDataPoint] = ...
%                                             findpeaks(tempPatch(obj.HFSDataPt{analyseStimN(HFSStimI)}),...
%                                             'MINPEAKHEIGHT', (updateStimCutOff - stimBase)*lowAPThreshFac + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
% %                                         longThreshFac = 0.2;
%                                     end
%                                 else
%                                     updateStimCutOff = -10;
%                                     [stimHighSpikeAmp, stimHighSpikeDataPoint] = ...
%                                             findpeaks(tempPatch(obj.HFSDataPt{analyseStimN(HFSStimI)}),...
%                                             'MINPEAKHEIGHT', (updateStimCutOff - stimBase)*highAPThreshFac + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));
%                                     
%                                 end

                                if smoothN ~= 0
                                    tempPatch = smooth(tempPatch,smoothN);
                                    
                                    
                                end


                                [stimSpikeAmp, stimSpikeDataPoint] = ...
                                            findpeaks(tempPatch(obj.HFSDataPt{analyseStimN(HFSStimI)}),...
                                            'MINPEAKHEIGHT', (adaptedSpikeAmp - stimBase)*APThreshFac + stimBase, 'MINPEAKDISTANCE', round(minGap*obj.SR));

%                                 orgPeakDataPt = tempPatchOrg(stimSpikeDataPoint);
                                %% HFS filtered spikes
%                                 stimArtefactEntry = find(stimSpikeAmp > 10);
%                                 
%                                 if (~isempty(stimArtefactEntry))
%                                     stimSpikeDataPoint(stimArtefactEntry) = [];
%                                     stimSpikeAmp(stimArtefactEntry) = [];
%                                 end
                                
                                obj.HFSFiltSpike{1,analyseStimN(HFSStimI)} = (obj.HFSDataPt{analyseStimN(HFSStimI)}(1) + stimSpikeDataPoint)/obj.SR;
                                % Time point
                                obj.HFSFiltSpike{2,analyseStimN(HFSStimI)} = stimSpikeAmp;
                                % spike amplitude
                                obj.HFSFiltSpike{3,analyseStimN(HFSStimI)} = obj.HFSDataPt{analyseStimN(HFSStimI)}(1) + stimSpikeDataPoint;
                                % data point

                                stimSpikeTempAmp = [];
                                stimSpikeTempDataPoint = [];
                                
                                if ~isempty(stimSpikeDataPoint)
                                    % Get the original signal peaks by
                                    % checking max around filtered signal
                                    % peak
                                    for i = 1:length(stimSpikeDataPoint)
                                        tempOrgPeakTest =...
                                            (stimSpikeDataPoint(i) - round(minGap*obj.SR*0.4)):(stimSpikeDataPoint(i) + round(minGap*obj.SR*0.4));
                                        % extending points around filtered
                                        % peak (slightly less than
                                        % minGap/2)
                                        testHFSPt = obj.HFSDataPt{analyseStimN(HFSStimI)}(1) + tempOrgPeakTest;
                                        
                                        [tempOrgPeak, tempOrgPeakPt] = max(tempPatchOrg(testHFSPt));
                                        
                                        testMaxHFSPt = testHFSPt(1) + tempOrgPeakPt - 1; 
                                        
                                        testPeak = tempPatch(testMaxHFSPt);
                                        % Get the peak position found in
                                        % original data and use that
                                        % position to find the
                                        % corresponding filtered peak with
                                        % this position to test for the if
                                        % statement below 
                                        
                                        if (testPeak - stimBase)/(stimSpikeAmp(i) - stimBase) > 1
                                            % If the new peak is much
                                            % higher than the old one, it
                                            % means the old one is a false
                                            % peak
                                            stimSpikeTempAmp(i) = NaN;
                                            stimSpikeTempDataPoint(i) = NaN;
                                            
                                        else
%                                             tempOrgPeakPt = tempOrgPeakPt + obj.HFSDataPt{analyseStimN(HFSStimI)}(1);
                                            tempOrgPeakPt = tempOrgPeakPt + tempOrgPeakTest(1);
                                            stimSpikeTempAmp(i) = tempOrgPeak;
                                            stimSpikeTempDataPoint(i) = tempOrgPeakPt;
                                            
                                        end
                                    end
                                    
                                    stimSpikeAmp = stimSpikeTempAmp;
                                    stimSpikeDataPoint = stimSpikeTempDataPoint;
                                end
                                
                                
                                obj.HFSSpike{1,analyseStimN(HFSStimI)} = (obj.HFSDataPt{analyseStimN(HFSStimI)}(1) + stimSpikeDataPoint)/obj.SR;
                                % Time point
                                obj.HFSSpike{2,analyseStimN(HFSStimI)} = stimSpikeAmp;
                                % spike amplitude
                                obj.HFSSpike{3,analyseStimN(HFSStimI)} = obj.HFSDataPt{analyseStimN(HFSStimI)}(1) + stimSpikeDataPoint;
                                % data point
                            end
%                             end

                            
                        end
                    end
                    obj.noExtStimSweepN = noExtStimSweep;
                end
            end
        end
        
        function getHFSPotential(obj)
            count = 1;
            
            for HFSSweepI = 1:obj.nSweeps
                tempHFSDataPt = obj.extStimStartPt(HFSSweepI):obj.extStimEndPt(HFSSweepI);
                tempPatch = obj.HFSfilteredSig(HFSSweepI,tempHFSDataPt);
            
                obj.HFSPotential(HFSSweepI) = median(tempPatch);
                
                if HFSSweepI > 1
                    if mod(HFSSweepI, 3) == 1
                        sumHFSPotential = obj.HFSPotential(HFSSweepI) + ...
                            obj.HFSPotential(HFSSweepI - 1) + ...
                            obj.HFSPotential(HFSSweepI - 2);
                        obj.HFSPotentialPerAmp(count) = sumHFSPotential/3;
                        count = count + 1;
                    end
                end
                
            end
            
            
        end
        
        function getHFSSponSpikeN(obj)
            
            for stimI = 1:size(obj.baseSpike,2)
                obj.HFSBaseSponSpikeN(stimI) = length(obj.baseSpike{2,stimI});
                obj.HFSAStimSponSpikeN(stimI) = length(obj.AStimSpike{2,stimI});
                
%                 obj.HFSSponSpikeN(stimI) = baseSponSpikeN(stimI) +...
%                     AStimSponSpikeN(stimI);
                
            end
            
            
            
        end
        
        function getLightSpike(obj, cutOffThreshPlus, minGap)
            if (~isempty(obj.LED))&(~isempty(obj.stimStartDataPoint))
                
                
                obj.spikeletCutOff = -40;
                
                if nargin > 1
                    obj.spikeCutOff = cutOffThreshPlus;
                else
                    obj.spikeCutOff = -20;
                end
                
                
                if obj.ndim == 2
                    tempPatch = obj.patch;
                    %% Baseline spikes
                    baselineBase = median(tempPatch(obj.baselineDataPoint));
                    baselineCutOff = baselineBase + obj.spikeCutOff;
                    
                    [baseSpikeAmp, baseSpikeDataPoint] = ...
                        findpeaks(tempPatch(obj.baselineDataPoint), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                    % use the 1st column of baselineDataPoint because
                    % patch is in sweep format and the datapoint starts
                    % from 1 for each sweep
                    
                    adjustBaseDataPoint = obj.baselineDataPoint(1) + baseSpikeDataPoint;
                    
                    [baseSpikeletAmp, baseSpikeletDataPoint] = obj.sortSpike(tempPatch, baseSpikeAmp, adjustBaseDataPoint,obj.spikeCutOff);
                    
                    
                    obj.baseSpikelet{1,1} = obj.dataTime(baseSpikeletDataPoint);
                    obj.baseSpikelet{2,1} = baseSpikeletAmp;
                    obj.baseSpikelet{3,1} = baseSpikeletDataPoint;
                    
                    baseSpikeDataPoint(baseSpikeAmp < baselineCutOff) = [];
                    baseSpikeAmp(baseSpikeAmp < baselineCutOff) = [];
                    
                    
                    obj.baseSpike{1,1} = obj.dataTime(obj.baselineDataPoint(1) + baseSpikeDataPoint - 1);
                    obj.baseSpike{2,1} = baseSpikeAmp;
                    obj.baseSpike{3,1} = obj.baselineDataPoint(1) + baseSpikeDataPoint;
                    
                    for i = 1:length(obj.stimStartTime)
                        
                        
                        %% Light ON spike
                        [LEDONSpikeAmp, LEONDataPoint] = ...
                            findpeaks(tempPatch(obj.LEDONDataPoint{i}), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                        
                        LEDONBase = median(tempPatch(obj.LEDONDataPoint{i}));
                        LEDONCutOff = LEDONBase + obj.spikeCutOff;
                        
                        adjustLEDONDataPoint = obj.LEDONDataPoint{i}(1) + LEONDataPoint;
                        
                        [LEDONSpikeletAmp, LEDONSpikeletDataPoint] = obj.sortSpike(tempPatch, LEDONSpikeAmp, adjustLEDONDataPoint,obj.spikeCutOff);
                        
                        obj.LEDONSpikelet{1,i} = obj.dataTime(LEDONSpikeletDataPoint);
                        obj.LEDONSpikelet{2,i} = LEDONSpikeletAmp;
                        obj.LEDONSpikelet{3,i} = LEDONSpikeletDataPoint;
                        
                        LEONDataPoint(LEDONSpikeAmp < LEDONCutOff) = [];
                        LEDONSpikeAmp(LEDONSpikeAmp < LEDONCutOff) = [];
                        
                        
                        obj.LEDONSpike{1,i} = obj.dataTime(obj.LEDONDataPoint{i}(1) + LEONDataPoint - 1);
                        obj.LEDONSpike{2,i} = LEDONSpikeAmp;
                        obj.LEDONSpike{3,i} = obj.LEDONDataPoint{i}(1) + LEONDataPoint;
                        
                        
                        
                        
                        
                        %% Light OFF spike
                        % recording spontaneous or LED responses, assuming
                        % spike threshold to be -20mV
                        
                        [LEDOFFSpikeAmp, LEDOFFSpikeDataPoint] = ...
                            findpeaks(tempPatch(obj.LEDOFFDataPoint{i}), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                        
                        LEDOFFBase = median(tempPatch(obj.LEDOFFDataPoint{i}));
                        LEDOFFCutOff = LEDOFFBase + obj.spikeCutOff;
                        
                        adjustLEDOFFDataPoint = obj.LEDOFFDataPoint{i}(1) + LEDOFFSpikeDataPoint;
                        [LEDOFFSpikeletAmp, LEDOFFSpikeletDataPoint] = obj.sortSpike(tempPatch, LEDOFFSpikeAmp, adjustLEDOFFDataPoint,obj.spikeCutOff);
                        
                        obj.LEDOFFSpikelet{1,i} = obj.dataTime(LEDOFFSpikeletDataPoint);
                        obj.LEDOFFSpikelet{2,i} = LEDOFFSpikeletAmp;
                        obj.LEDOFFSpikelet{3,i} = LEDOFFSpikeletDataPoint;
                        
                        LEDOFFSpikeDataPoint(LEDOFFSpikeAmp < LEDOFFCutOff) = [];
                        LEDOFFSpikeAmp(LEDOFFSpikeAmp < LEDOFFCutOff) = [];
                        
                        
                        obj.LEDOFFSpike{1,i} = obj.dataTime(obj.LEDOFFDataPoint{i}(1) + LEDOFFSpikeDataPoint - 1);
                        obj.LEDOFFSpike{2,i} = LEDOFFSpikeAmp;
                        obj.LEDOFFSpike{3,i} = obj.LEDOFFDataPoint{i}(1) + LEDOFFSpikeDataPoint;
                        
                        %                         obj.spikeData{1,i} = [obj.stimSpike{1,i}, obj.AStimSpike{1,i}];
                        %                         obj.spikeData{2,i} = [obj.stimSpike{2,i}, obj.AStimSpike{2,i}];
                        %                         obj.spikeData{3,i} = [obj.stimSpike{3,i}, obj.AStimSpike{3,i}];
                        %
                        %                         obj.spikeletData{1,i} = [obj.stimSpikelet{1,i}, obj.AStimSpikelet{1,i}];
                        %                         obj.spikeletData{2,i} = [obj.stimSpikelet{2,i}, obj.AStimSpikelet{2,i}];
                        %                         obj.spikeletData{3,i} = [obj.stimSpikelet{3,i}, obj.AStimSpikelet{3,i}];
                        
                    end
                    
                    
                elseif obj.ndim == 3
                    for episodicI = 1:obj.nSweeps
                        tempPatch = obj.patch(episodicI,:);
                        %% Baseline spikes
                        if ~isempty(obj.baselineDataPoint)
                            baselineBase = median(tempPatch(obj.baselineDataPoint{episodicI}));
                            baselineCutOff = baselineBase + obj.spikeCutOff;
                            
                            [baseSpikeAmp, baseSpikeDataPoint] = ...
                                findpeaks(tempPatch(obj.baselineDataPoint{episodicI}), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                            [msg id] = lastwarn;
                            warning('off',id);
                            % use the 1st column of baselineDataPoint because
                            % patch is in sweep format and the datapoint starts
                            % from 1 for each sweep
                            
                            adjustBaseDataPoint = obj.baselineDataPoint{episodicI}(1) + baseSpikeDataPoint;
                            
                            [baseSpikeletAmp, baseSpikeletDataPoint] = obj.sortSpike(tempPatch, baseSpikeAmp, adjustBaseDataPoint,obj.spikeCutOff);
                            
                            
                            obj.baseSpikelet{1,episodicI} = obj.dataTime(baseSpikeletDataPoint);
                            obj.baseSpikelet{2,episodicI} = baseSpikeletAmp;
                            obj.baseSpikelet{3,episodicI} = baseSpikeletDataPoint;
                            
                            baseSpikeDataPoint(baseSpikeAmp < baselineCutOff) = [];
                            baseSpikeAmp(baseSpikeAmp < baselineCutOff) = [];
                            
                            
                            obj.baseSpike{1,episodicI} = (obj.baselineDataPoint{episodicI}(1) + baseSpikeDataPoint - 1)/obj.SR;
                            obj.baseSpike{2,episodicI} = baseSpikeAmp;
                            obj.baseSpike{3,episodicI} = obj.baselineDataPoint{episodicI}(1) + baseSpikeDataPoint - 1;
                        end
                        
                        %% After stimulation  spikes
                        if ~isempty(obj.AStimDataPoint)
                            
                            [AStimSpikeAmp, AStimDataPoint] = ...
                                findpeaks(tempPatch(obj.AStimDataPoint{episodicI}), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                            
                            [msg id] = lastwarn;
                            warning('off',id);
                            
                            AStimBase = median(tempPatch(obj.AStimDataPoint{episodicI}));
                            AStimCutOff = AStimBase + obj.spikeCutOff;
                            
                            
                            adjustAStimDataPoint = obj.AStimDataPoint{episodicI}(1) + AStimDataPoint;
                            [AStimSpikeletAmp, AStimSpikeletDataPoint] = obj.sortSpike(tempPatch, AStimSpikeAmp, adjustAStimDataPoint,obj.spikeCutOff);
                            
                            obj.AStimSpikelet{1,episodicI} = obj.recordTime(episodicI,AStimSpikeletDataPoint);
                            obj.AStimSpikelet{2,episodicI} = AStimSpikeletAmp;
                            obj.AStimSpikelet{3,episodicI} = AStimSpikeletDataPoint;
                            
                            AStimDataPoint(AStimSpikeAmp < AStimCutOff) = [];
                            AStimSpikeAmp(AStimSpikeAmp < AStimCutOff) = [];
                            
                            
                            obj.AStimSpike{1,episodicI} = (obj.AStimDataPoint{episodicI}(1) + AStimDataPoint - 1)/obj.SR;
                            obj.AStimSpike{2,episodicI} = AStimSpikeAmp;
                            obj.AStimSpike{3,episodicI} = obj.AStimDataPoint{episodicI}(1) + AStimDataPoint - 1;
                            
                            
                        end
                        %                         for i = 1:length(obj.stimStartTime)
                        
                        for LEDNum = 1:size(obj.LEDONDataPoint,2)
                            %% Light ON spike
                            [LEDONSpikeAmp, LEDONSpikeDataPoint] = ...
                                findpeaks(tempPatch(obj.LEDONDataPoint{episodicI,LEDNum}), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                            
                            [msg id] = lastwarn;
                            warning('off',id);
                            
                            LEDONBase = median(tempPatch(obj.LEDONDataPoint{episodicI,LEDNum}));
                            LEDONCutOff = LEDONBase + obj.spikeCutOff;
                            
                            adjustLEDONDataPoint = obj.LEDONDataPoint{episodicI,LEDNum}(1) + LEDONSpikeDataPoint;
                            
                            [LEDONSpikeletAmp, LEDONSpikeletDataPoint] = obj.sortSpike(tempPatch, LEDONSpikeAmp, adjustLEDONDataPoint,obj.spikeCutOff);
                            
                            obj.LEDONSpikelet{1, episodicI, LEDNum} = obj.dataTime(LEDONSpikeletDataPoint);
                            obj.LEDONSpikelet{2, episodicI, LEDNum} = LEDONSpikeletAmp;
                            obj.LEDONSpikelet{3, episodicI, LEDNum} = LEDONSpikeletDataPoint;
                            
                            LEDONSpikeDataPoint(LEDONSpikeAmp < LEDONCutOff) = [];
                            LEDONSpikeAmp(LEDONSpikeAmp < LEDONCutOff) = [];
                            
                            
                            obj.LEDONSpike{1, episodicI, LEDNum} = (obj.LEDONDataPoint{episodicI,LEDNum}(1) + LEDONSpikeDataPoint - 1)/obj.SR;
                            obj.LEDONSpike{2, episodicI, LEDNum} = LEDONSpikeAmp;
                            obj.LEDONSpike{3, episodicI, LEDNum} = obj.LEDONDataPoint{episodicI,LEDNum}(1) + LEDONSpikeDataPoint - 1;
                            
                            obj.LEDONIndSpikeRate{episodicI, LEDNum} = 1./diff(obj.LEDONSpike{1, episodicI, LEDNum});
                            
                            
                            
                            %% Light OFF spike
                            % recording spontaneous or LED responses, assuming
                            % spike threshold to be -20mV
                            
                            [LEDOFFSpikeAmp, LEDOFFSpikeDataPoint] = ...
                                findpeaks(tempPatch(obj.LEDOFFDataPoint{episodicI,LEDNum}), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                            
                            [msg id] = lastwarn;
                            warning('off',id);
                            
                            LEDOFFBase = median(tempPatch(obj.LEDOFFDataPoint{episodicI,LEDNum}));
                            LEDOFFCutOff = LEDOFFBase + obj.spikeCutOff;
                            
                            adjustLEDOFDataPoint = obj.LEDOFFDataPoint{episodicI,LEDNum}(1) + LEDOFFSpikeDataPoint;
                            [LEDOFFSpikeletAmp, LEDOFFSpikeletDataPoint] = obj.sortSpike(tempPatch, LEDOFFSpikeAmp, adjustLEDOFDataPoint,obj.spikeCutOff);
                            
                            obj.LEDOFFSpikelet{1, episodicI, LEDNum} = obj.dataTime(LEDOFFSpikeletDataPoint);
                            obj.LEDOFFSpikelet{2, episodicI, LEDNum} = LEDOFFSpikeletAmp;
                            obj.LEDOFFSpikelet{3, episodicI, LEDNum} = LEDOFFSpikeletDataPoint;
                            
                            LEDOFFSpikeDataPoint(LEDOFFSpikeAmp < LEDOFFCutOff) = [];
                            LEDOFFSpikeAmp(LEDOFFSpikeAmp < LEDOFFCutOff) = [];
                            
                            
                            obj.LEDOFFSpike{1, episodicI, LEDNum} = (obj.LEDOFFDataPoint{episodicI,LEDNum}(1) + LEDOFFSpikeDataPoint - 1)/obj.SR;
                            obj.LEDOFFSpike{2, episodicI, LEDNum} = LEDOFFSpikeAmp;
                            obj.LEDOFFSpike{3, episodicI, LEDNum} = obj.LEDOFFDataPoint{episodicI,LEDNum}(1) + LEDOFFSpikeDataPoint - 1;
                            
                            obj.LEDOFFIndSpikeRate{episodicI, LEDNum} = 1./diff(obj.LEDOFFSpike{1, episodicI, LEDNum});
                        end
                    end
                end
            end
        end
        
        function getSpikeRate(obj, varargin)
            
            % Common parameters for all different stimulation protocols
            
            if obj.ndim == 2
                for stimN = 1:size(obj.baseSpike, 2)
                    %                 baseSpike = obj.baseSpike{1,baseN};
                    baseSpikeTime = [0, obj.baseSpike{1,stimN}'];
                    
                    obj.baselineIndSpikeRate{stimN} = 1./diff(baseSpikeTime);
                    
                    obj.baselineGeneSpikeRate(stimN) = ...
                        numel(baseSpikeTime)/((obj.baselineDataPoint{stimN}(end) - obj.baselineDataPoint{stimN}(1))/obj.SR);
                    
                    AStimSpikeTime = [obj.AStimDataPoint{stimN}(1)/obj.SR, obj.AStimSpike{1,stimN}'];
                    
                    obj.AStimIndSpikeRate{stimN} = 1./diff(AStimSpikeTime);
                    
                    obj.AStimGeneSpikeRate(stimN) = ...
                        numel(AStimSpikeTime)/((obj.AStimDataPoint{stimN}(end) - obj.AStimDataPoint{stimN}(1))/obj.SR);
                    
                end
              
            else
                
                for sweepN = 1:size(obj.baseSpike, 2)
                    %                 baseSpike = obj.baseSpike{1,baseN};
                    
                    if ~isempty(obj.baseSpike{1,sweepN})
                        baseSpikeTime = [0, obj.baseSpike{1,sweepN}];
                        
                        obj.baselineIndSpikeRate{sweepN} = 1./diff(baseSpikeTime);
                        
                        obj.baselineGeneSpikeRate(sweepN) = ...
                            numel(baseSpikeTime)/((obj.baselineDataPoint{sweepN}(end) - obj.baselineDataPoint{sweepN}(1))/obj.SR);
                    else
                        obj.baselineGeneSpikeRate(sweepN) = 0;
                    end
                    
                    if ~isempty(obj.AStimSpike{1,sweepN})
                        AStimSpikeTime = [obj.AStimDataPoint{sweepN}(1)/obj.SR,obj.AStimSpike{1,sweepN}];
                        
                        obj.AStimIndSpikeRate{sweepN} = 1./diff(AStimSpikeTime);
                        
                        obj.AStimGeneSpikeRate(sweepN) = ...
                            numel(AStimSpikeTime)/((obj.AStimDataPoint{sweepN}(end) - obj.AStimDataPoint{sweepN}(1))/obj.SR);
                    else
                        
                        obj.AStimGeneSpikeRate(sweepN) = 0;
                    end
                end
            end
            
            
            
            % Specific parameters for each different stimulation protocol
                
            if ~isempty(strfind(obj.protocol, 'IC AP family'))
                obj.getIntSpikeRate(varargin);
                
            elseif ~isempty(strfind(obj.protocol, 'LED'))
                obj.getLightSpikeRate(varargin);
                
            elseif ~isempty(strfind(obj.protocol, 'SPS'))
                obj.getSPSSpikeRate(varargin);
                
            elseif ~isempty(strfind(obj.protocol, 'HFS'))
                obj.getHFSSpikeRate(varargin);
                
            end            
        end
        
        function getIntSpikeRate(obj, varargin)
            for sweepN = 1:obj.nSweeps
                
                stimSpikeTime = [obj.stimDataPoint{sweepN}(1)/obj.SR, obj.stimSpike{1,sweepN}];
                
                obj.stimIndSpikeRate{sweepN} = 1./diff(stimSpikeTime);
                
                obj.stimGeneSpikeRate(sweepN) = ...
                    numel(stimSpikeTime)/((obj.stimDataPoint{sweepN}(end) - obj.stimDataPoint{sweepN}(1))/obj.SR);
                obj.stimSpikeNumber(sweepN) = numel(stimSpikeTime);
                
            end
        end
        
        function getLightSpikeRate(obj, varargin)
            
            if obj.ndim == 2
                for LEDN = 1:size(obj.LEDONSpike, 2)
                    LEDONSpikeTime = [obj.LEDONDataPoint{LEDN}(1)/obj.SR, obj.LEDONSpike{1,LEDN}];
                    
                    obj.LEDONIndSpikeRate{LEDN} = 1./diff(LEDONSpikeTime);
                    
                    obj.LEDONGeneSpikeRate(LEDN) = ...
                        numel(LEDONSpikeTime)/((obj.LEDONDataPoint{LEDN}(end) - obj.LEDONDataPoint{LEDN}(1))/obj.SR);
                    
                    obj.LEDONSpikeNumber(LEDN) = numel(LEDONSpikeTime);
                    
                    LEDOFFSpikeTime = [obj.LEDOFFDataPoint{LEDN}(1)/obj.SR, obj.LEDOFFSpike{1,LEDN}];
                    
                    obj.LEDOFFIndSpikeRate{LEDN} = 1./diff(LEDOFFSpikeTime);
                    
                    obj.LEDOFFGeneSpikeRate(LEDN) = ...
                        numel(LEDOFFSpikeTime)/((obj.LEDOFFDataPoint{LEDN}(end) - obj.LEDOFFDataPoint{LEDN}(1))/obj.SR);
                    
                    obj.LEDOFFSpikeNumber(LEDN) = numel(LEDOFFSpikeTime);
                    
                    
                    
                end
            else
                for sweepN = 1:size(obj.LEDONSpike, 2)
                    for LEDN = 1:size(obj.LEDONSpike, 3)
                        LEDONSpikeTime = [obj.LEDONDataPoint{LEDN}(1)/obj.SR, obj.LEDONSpike{1,sweepN, LEDN}];
                        
                        obj.LEDONIndSpikeRate{sweepN, LEDN} = 1./diff(LEDONSpikeTime);
                        
                        obj.LEDONGeneSpikeRate(sweepN, LEDN) = ...
                            numel(LEDONSpikeTime)/((obj.LEDONDataPoint{sweepN, LEDN}(end) - obj.LEDONDataPoint{sweepN, LEDN}(1))/obj.SR);
                        
                        obj.LEDONSpikeNumber(sweepN, LEDN) = numel(LEDONSpikeTime);
                        
                        LEDOFFSpikeTime = [obj.LEDOFFDataPoint{LEDN}(1)/obj.SR, obj.LEDOFFSpike{1,sweepN, LEDN}];
                        
                        obj.LEDOFFIndSpikeRate{sweepN, LEDN} = 1./diff(LEDOFFSpikeTime);
                        
                        obj.LEDOFFGeneSpikeRate(sweepN, LEDN) = ...
                            numel(LEDOFFSpikeTime)/((obj.LEDOFFDataPoint{sweepN, LEDN}(end) - obj.LEDOFFDataPoint{sweepN, LEDN}(1))/obj.SR);
                        
                        obj.LEDOFFSpikeNumber(sweepN, LEDN) = numel(LEDOFFSpikeTime);
                        
                    end
                end
                
            end

        end
        
        function getMedianMbmV(obj)
            % Get the median membrane voltage of each recording
%             sumVal = 0;
            
            for sweepI = 1:obj.nSweeps
                
                obj.sweepMedianMV(sweepI) = median(obj.patch(sweepI,:));
                % Median of membrane voltage of every sweep
        
%                 sumVal = sumVal + obj.sweepMedianMV(sweepI);
            end
            
            obj.avgMedianMV = mean(obj.sweepMedianMV);
            
            obj.stdMedianMV = std(obj.sweepMedianMV);
            
            
        end
        
        function getDiffVm(obj)
        % To get the difference in membrane voltage between resting membrane potential and HFSstimulation period.   
        count = 1;
        obj.HFSDiffVm = [];
            for sweepI = 1:obj.nSweeps
                
                obj.HFSSweepMedianMV(sweepI) = median(obj.patch(sweepI,obj.baselineDataPoint{sweepI}));
                % Median of membrane voltage of every sweep
                obj.HFSMedianMV(sweepI) = median(obj.patch(sweepI,obj.HFSDataPt{sweepI}));
%                 sumVal = sumVal + obj.sweepMedianMV(sweepI);
                temp(sweepI) = obj.HFSSweepMedianMV(sweepI) - obj.HFSMedianMV(sweepI);
                if (mod(sweepI, 3) == 1)&(sweepI ~= 1)
                    obj.HFSDiffVm(count) = (temp(sweepI - 2) + temp(sweepI - 1) + temp(sweepI))/3;
                    count = count + 1;
                end
            end
            
        
        
        end
        
        
        
        function getSPSSpikeRate(obj, varargin)
            if obj.ndim == 2
                for SPSN = 1:size(obj.HFSSpike, 2)
                    SPSSpikeTime = [obj.SPSDataPt{SPSN}(1)/obj.SR, obj.SPSSpike{1,SPSN}'];
                    
                    obj.SPSIndSpikeRate{SPSN} = 1./diff(SPSSpikeTime);
                    
                    obj.SPSGeneSpikeRate(SPSN) = ...
                        numel(SPSSpikeTime)/((obj.SPSDataPt{SPSN}(end) - obj.SPSDataPt{SPSN}(1))/obj.SR);
                    
                    obj.SPSSpikeNumber(SPSN) = numel(SPSSpikeTime);
                    
                    if SPSN ~= size(obj.SPSSpike, 2)
                        SPSGapSpikeTime = [obj.SPSGapDataPoint{SPSN}(1)/obj.SR, obj.SPSGapSpike{1,SPSN}'];
                        
                        obj.SPSGapIndSpikeRate{SPSN} = 1./diff(SPSGapSpikeTime);
                        
                        obj.SPSGapGeneSpikeRate(SPSN) = ...
                            numel(SPSGapSpikeTime)/((obj.SPSGapDataPoint{SPSN}(end) - obj.SPSGapDataPoint{SPSN}(1))/obj.SR);
                        
                        obj.SPSGapSpikeNumber(SPSN) = numel(SPSGapSpikeTime);
                    end
                end

                
                
            else
                if ~isempty(obj.SPSDataPt)
                    for sweepN = 1:obj.nSweeps
                        if isempty(obj.SPSDataPt{sweepN})
                            continue;
                        end
                        
                        SPSSpikeTime = [obj.SPSDataPt{sweepN}(1)/obj.SR, obj.SPSSpike{1,sweepN}];
                        
                        obj.SPSIndSpikeRate{sweepN} = 1./diff(SPSSpikeTime);
                        
                        obj.SPSGeneSpikeRate(sweepN) = ...
                            numel(SPSSpikeTime)/((obj.SPSDataPt{sweepN}(end) - obj.SPSDataPt{sweepN}(1))/obj.SR);
                        
                        obj.SPSSpikeNumber(sweepN) = numel(SPSSpikeTime);
                        
                    end
                end
                
            end
            
            
            
            
            
            SPSIndSpikeRate = [];
            SPSGeneSpikeRate = '';
        end
        
        
        function getHFSSpikeRate(obj, varargin)
            
            if obj.ndim == 2
                for HFSN = 1:size(obj.HFSSpike, 2)
                    HFSSpikeTime = [obj.HFSDataPt{HFSN}(1)/obj.SR, obj.HFSSpike{1,HFSN}];
                    
                    obj.HFSIndSpikeRate{HFSN} = 1./diff(HFSSpikeTime);
                    
                    obj.HFSGeneSpikeRate(HFSN) = ...
                        numel(HFSSpikeTime)/((obj.HFSDataPt{HFSN}(end) - obj.HFSDataPt{HFSN}(1))/obj.SR);
                    
                    obj.HFSSpikeNumber(HFSN) = numel(HFSSpikeTime);
                    
                    if HFSN ~= size(obj.HFSSpike, 2)
                        HFSGapSpikeTime = [obj.HFSGapDataPoint{HFSN}(1)/obj.SR, obj.HFSGapSpike{1,HFSN}'];
                        
                        obj.HFSGapIndSpikeRate{HFSN} = 1./diff(HFSGapSpikeTime);
                        
                        obj.HFSGapGeneSpikeRate(HFSN) = ...
                            numel(HFSGapSpikeTime)/((obj.HFSGapDataPoint{HFSN}(end) - obj.HFSGapDataPoint{HFSN}(1))/obj.SR);
                        
                        obj.HFSGapSpikeNumber(HFSN) = numel(HFSGapSpikeTime);
                    end
                end
                
            else
                for sweepN = 1:obj.nSweeps
                    
                    if isempty(obj.HFSDataPt{sweepN})
                        continue;
                    end
                    
                    HFSSpikeTime = [obj.HFSDataPt{sweepN}(1)/obj.SR, obj.HFSSpike{1,sweepN}];
                    
                    obj.HFSIndSpikeRate{sweepN} = 1./diff(HFSSpikeTime);
                    
                    obj.HFSGeneSpikeRate(sweepN) = ...
                        numel(HFSSpikeTime)/((obj.HFSDataPt{sweepN}(end) - obj.HFSDataPt{sweepN}(1))/obj.SR);
                    
                    obj.HFSSpikeNumber(sweepN) = numel(HFSSpikeTime) - 1;
                    
                end
                       
            end
            
            
            
            
            
            HFSIndSpikeRate = [];
            HFSGeneSpikeRate = '';
        end
        
        function getIntSpike(obj, cutOffThreshPlus, minGap)
            % This is for single pulse stimulation (SPS)
            
            % rewritten get spikes on 20150608, mainly changing the spike
            % threshold to fixed values. Make it a general method for both
            % intracellular stim and extraceullar stim since both will be
            % in episodic mode now.
            % minGap (in s) is the minimum distance between spike for findpeak
            % function
            if ~isempty(obj.patch)
                %                 totalSpikeTemp = [];
                if (~isempty(obj.intStim))|(~isempty(obj.extStim))
                    for i = 1:obj.nSweeps
                        
                        obj.spikeletCutOff = -50;
                        
                        if nargin > 1
                            obj.spikeCutOff = cutOffThreshPlus;
                            % how much above baseline
                        else
                            obj.spikeCutOff = 20;
                        end

                        tempPatch = obj.patch(i,:);

                        %% Baseline spikes
                        %                         adjustBaselineDataPoint =...
                        %                             obj.baselineDataPoint(i,:) - obj.baselineDataPoint(i,1) + 1;
                        if ~isempty(obj.baselineDataPoint)
                            [baseSpikeAmp, baseSpikeDataPoint] = ...
                                findpeaks(tempPatch(obj.baselineDataPoint{i}), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                            % use the 1st column of baselineDataPoint because
                            % patch is in sweep format and the datapoint starts
                            % from 1 for each sweep
                            
                            baselineBase = median(tempPatch(obj.baselineDataPoint{i}));
                            baselineCutOff = baselineBase + obj.spikeCutOff;
                            
                            adjustBaseDataPoint = obj.baselineDataPoint{i}(1) + baseSpikeDataPoint;
                            
                            [baseSpikeletAmp, baseSpikeletDataPoint] = obj.sortSpike(tempPatch, baseSpikeAmp, adjustBaseDataPoint,obj.spikeCutOff);
                            
                            
                            obj.baseSpikelet{1,i} = obj.recordTime(i,baseSpikeletDataPoint);
                            obj.baseSpikelet{2,i} = baseSpikeletAmp;
                            obj.baseSpikelet{3,i} = baseSpikeletDataPoint;
                            
                            baseSpikeDataPoint(baseSpikeAmp < baselineCutOff) = [];
                            baseSpikeAmp(baseSpikeAmp < baselineCutOff) = [];
                            
                            
                            obj.baseSpike{1,i} = (obj.baselineDataPoint{i}(1) + baseSpikeDataPoint - 1)/obj.SR;
                            obj.baseSpike{2,i} = baseSpikeAmp;
                            obj.baseSpike{3,i} = obj.baselineDataPoint{i}(1) + baseSpikeDataPoint - 1;
                        end
                        
                        
                        %% Stim spikes
                        %                         adjustStimDataPoint =...
                        %                             obj.stimDataPoint(i,:) - obj.baselineDataPoint(i,1) + 1;
                        
                        if ~isempty(obj.stimDataPoint)
                            if length(obj.stimDataPoint{i}) > round(0.001*obj.SR)
                                [stimSpikeAmp, stimSpikeDataPoint] = ...
                                    findpeaks(tempPatch(obj.stimDataPoint{i}), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                                
                            else
                                [stimSpikeAmp, stimSpikeDataPoint] = ...
                                    findpeaks(tempPatch(obj.stimDataPoint{i}), 'MINPEAKHEIGHT', obj.spikeletCutOff);
                            end
                            
                            stimBase = median(tempPatch(obj.stimDataPoint{i}));
                            stimCutOff = stimBase + obj.spikeCutOff;
                            
                            adjustStimDataPoint = obj.stimDataPoint{i}(1) + stimSpikeDataPoint;
                            
                            [stimSpikeletAmp, stimSpikeletDataPoint] = obj.sortSpike(tempPatch, stimSpikeAmp, adjustStimDataPoint,obj.spikeCutOff);
                            
                            obj.stimSpikelet{1,i} = obj.recordTime(i,stimSpikeletDataPoint);
                            obj.stimSpikelet{2,i} = stimSpikeletAmp;
                            obj.stimSpikelet{3,i} = stimSpikeletDataPoint;
                            
                            stimSpikeDataPoint(stimSpikeAmp < stimCutOff) = [];
                            stimSpikeAmp(stimSpikeAmp < stimCutOff) = [];
                            
                            
                            obj.stimSpike{1,i} = (obj.stimDataPoint{i}(1) + stimSpikeDataPoint - 1)/obj.SR;
                            obj.stimSpike{2,i} = stimSpikeAmp;
                            obj.stimSpike{3,i} = obj.stimDataPoint{i}(1) + stimSpikeDataPoint - 1;
                            
                        end
                        
                        %% After stim spikes
                        
                        %                         adjustAStimDataPoint =...
                        %                             obj.AStimDataPoint(i,:) - obj.baselineDataPoint(i,1) + 1;
                        
                        if ~isempty(obj.AStimDataPoint)
                            
                            [AStimSpikeAmp, AStimSpikeDataPoint] = ...
                                findpeaks(tempPatch(obj.AStimDataPoint{i}), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                            
                            sponBase = median(tempPatch(obj.AStimDataPoint{i}));
                            sponCutOff = sponBase + obj.spikeCutOff;
                            
                            
                            adjustAStimDataPoint = obj.AStimDataPoint{i}(1) + AStimSpikeDataPoint;
                            [AStimSpikeletAmp, AStimSpikeletDataPoint] = obj.sortSpike(tempPatch, AStimSpikeAmp, adjustAStimDataPoint,obj.spikeCutOff);
                            
                            obj.AStimSpikelet{1,i} = obj.recordTime(i,AStimSpikeletDataPoint);
                            obj.AStimSpikelet{2,i} = AStimSpikeletAmp;
                            obj.AStimSpikelet{3,i} = AStimSpikeletDataPoint;
                            
                            AStimSpikeDataPoint(AStimSpikeAmp < sponCutOff) = [];
                            AStimSpikeAmp(AStimSpikeAmp < sponCutOff) = [];
                            
                            
                            obj.AStimSpike{1,i} = (obj.AStimDataPoint{i}(1) + AStimSpikeDataPoint - 1)/obj.SR;
                            obj.AStimSpike{2,i} = AStimSpikeAmp;
                            obj.AStimSpike{3,i} = obj.AStimDataPoint{i}(1) + AStimSpikeDataPoint - 1;
                            
                            obj.spikeData{1,i} = [obj.baseSpike{1,i}, obj.stimSpike{1,i}, obj.AStimSpike{1,i}];
                            obj.spikeData{2,i} = [obj.baseSpike{2,i}, obj.stimSpike{2,i}, obj.AStimSpike{2,i}];
                            obj.spikeData{3,i} = [obj.baseSpike{3,i}, obj.stimSpike{3,i}, obj.AStimSpike{3,i}];
                            
                            obj.spikeletData{1,i} = [obj.baseSpikelet{1,i}, obj.stimSpikelet{1,i}, obj.AStimSpikelet{1,i}];
                            obj.spikeletData{2,i} = [obj.baseSpikelet{2,i}, obj.stimSpikelet{2,i}, obj.AStimSpikelet{2,i}];
                            obj.spikeletData{3,i} = [obj.baseSpikelet{3,i}, obj.stimSpikelet{3,i}, obj.AStimSpikelet{3,i}];
                            
                        end
                        
                    end
                end
                
            else
                obj.spikeletCutOff = -40;
                
                if nargin > 1
                    obj.spikeCutOff = cutOffThreshPlus;
                else
                    obj.spikeCutOff = -20;
                end
                tempPatch = obj.patch;
                [spikeAmp, spikeDataPoint] = ...
                    findpeaks(tempPatch, 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                %                 gdSpikeTemp = selectspk(totalSpikeTemp);
                
                sponBase = median(tempPatch);
                sponCutOff = sponBase + obj.spikeCutOff;
                
                [spikeletAmp, spikeletDataPoint] = obj.sortSpike(tempPatch, spikeAmp, spikeDataPoint,obj.spikeCutOff);
                
                obj.spikeletData{1,1} = obj.dataTime(spikeletDataPoint);
                obj.spikeletData{2,1} = spikeletAmp;
                obj.spikeletData{3,1} = spikeletDataPoint;
                
                spikeDataPoint(spikeAmp < sponCutOff) = [];
                spikeAmp(spikeAmp < sponCutOff) = [];
                
                
                obj.spikeData{1} = obj.dataTime(spikeDataPoint);
                obj.spikeData{2} = spikeAmp;
                obj.spikeData{3} = spikeDataPoint;
                
            end
        end
        
        function getSPSSpike(obj, cutOffThreshPlus, minGap,varargin)
            % This is for single pulse stimulation (SPS)
            
            % rewritten get spikes on 20150608, mainly changing the spike
            % threshold to fixed values. Make it a general method for both
            % intracellular stim and extraceullar stim since both will be
            % in episodic mode now.
            % minGap (in s) is the minimum distance between spike for findpeak
            % function
            
            stimExtend = 20;
            % This is for extending the stimulation time for detecting peak
            % (mainly for inhibition investigation)
            
            coupledTimeFrame = 0.002; % 1ms time frame for a spike to be considered as an AP.
            
            
            artSpikeAmp = 30;
            
            if ~isempty(varargin)
                % Manually type in HFS Amp
                for i = 1:length(varargin)
                    
                    if strcmp(varargin{i}, 'HFSN')
                        HFSN = varargin{i + 1};
                    end
                    
                    if strcmp(varargin{i}, 'patchN')
                        patchN = varargin{i + 1};
                    end
                    
                    if strcmp(varargin{i}, 'stimExtend')
                        stimExtend = varargin{i + 1};
                    end
                    
                    if strcmp(varargin{i}, 'coupledTimeFrame')
                        coupledTimeFrame = varargin{i + 1};
                    end
                    
                    if strcmp(varargin{i}, 'artSpikeAmp')
                        artSpikeAmp = varargin{i + 1};
                        % artefact spiking amplitude, sometimes there is
                        % spike that is due to unknown reason
                    end
                    
                    if strcmp(varargin{i}, 'extSkipPt')
                        extSkipPt = varargin{i + 1};
                        % number of points after the start of stimulationso to
                        % avoid recognising artefacts as spike. Default is zero
                        % because it rarely occurs.
                    end
                    
                    
                    
                end
            end
            
            if ~isempty(obj.patch)
                %                 totalSpikeTemp = [];
                if (~isempty(obj.intStim))|(~isempty(obj.extStim))
                    for SPSStimN = 1:obj.nSweeps
                        
                        if obj.extStimStartPt(SPSStimN) == 0
                            continue;
                        end
                        
                        obj.spikeletCutOff = -50;
                        
                        if nargin > 1
                            obj.spikeCutOff = cutOffThreshPlus;
                            % how much above baseline
                        else
                            obj.spikeCutOff = 20;
                        end

                        tempPatch = obj.patch(SPSStimN,:);

                        %% Baseline spikes
                        %                         adjustBaselineDataPoint =...
                        %                             obj.baselineDataPoint(i,:) - obj.baselineDataPoint(i,1) + 1;
                        
                        obj.baselineDataPoint{SPSStimN} = 1:obj.extStimStartPt(SPSStimN);
                        
                        baselineBase = median(tempPatch(obj.baselineDataPoint{SPSStimN}));
                        baselineCutOff = baselineBase + obj.spikeCutOff;
%                         SPSStimN
%                         if SPSStimN == 100
%                             SPSStimN
%                         end
                        if obj.baselineDataPoint{SPSStimN} < round(minGap*obj.SR)
                            continue;
                        end
                        [baseSpikeAmp, baseSpikeDataPoint] = ...
                            findpeaks(tempPatch(obj.baselineDataPoint{SPSStimN}), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                        % use the 1st column of baselineDataPoint because
                        % patch is in sweep format and the datapoint starts
                        % from 1 for each sweep
                        
                        [msg id] = lastwarn;
                        if ~isempty(id)
                            warning('off',id);
                        end
                        
                        adjustBaseDataPoint = obj.baselineDataPoint{SPSStimN}(1) + baseSpikeDataPoint;
                        
                        [baseSpikeletAmp, baseSpikeletDataPoint] = obj.sortSpike(tempPatch, baseSpikeAmp, adjustBaseDataPoint,obj.spikeCutOff);
                        
                        
                        obj.baseSpikelet{1,SPSStimN} = obj.recordTime(SPSStimN,baseSpikeletDataPoint);
                        obj.baseSpikelet{2,SPSStimN} = baseSpikeletAmp;
                        obj.baseSpikelet{3,SPSStimN} = baseSpikeletDataPoint;
                        
                        baseSpikeDataPoint(baseSpikeAmp < baselineCutOff) = [];
                        baseSpikeAmp(baseSpikeAmp < baselineCutOff) = [];
                        
                        blArtefactEntry = find(baseSpikeAmp > artSpikeAmp);
                        if ~isempty(blArtefactEntry)
                            baseSpikeDataPoint(blArtefactEntry) = [];
                            baseSpikeAmp(blArtefactEntry) = [];
                        end
                        
                        
                        obj.baseSpike{1,SPSStimN} = (obj.baselineDataPoint{SPSStimN}(1) + baseSpikeDataPoint - 1)/obj.SR;
                        obj.baseSpike{2,SPSStimN} = baseSpikeAmp;
                        obj.baseSpike{3,SPSStimN} = obj.baselineDataPoint{SPSStimN}(1) + baseSpikeDataPoint - 1;

                        
                        
%                         %% Stim spikes
%                         %                         adjustStimDataPoint =...
%                         %                             obj.stimDataPoint(SPSStimN,:) - obj.baselineDataPoint(SPSStimN,1) + 1;
%                         obj.SPSDataPt{SPSStimN} = obj.extStimStartPt(SPSStimN):obj.extStimEndPt(SPSStimN);
%                         obj.SPSTime{SPSStimN} = obj.SPSDataPt{SPSStimN}/obj.SR;
%                         
%                         
%                         
%                         if length(obj.SPSDataPt{SPSStimN}) > round(0.001*obj.SR)
%                             [stimSpikeAmp, stimSpikeDataPoint] = ...
%                                 findpeaks(tempPatch(obj.SPSDataPt{SPSStimN}), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
%                         elseif length(obj.SPSDataPt{SPSStimN}) > 3 
%                           [stimSpikeAmp, stimSpikeDataPoint] = ...
%                                 findpeaks(tempPatch(obj.SPSDataPt{SPSStimN}), 'MINPEAKHEIGHT', obj.spikeletCutOff);
%                             
%                             
%                             stimBase = median(tempPatch(obj.SPSDataPt{SPSStimN}));
%                             stimCutOff = stimBase + obj.spikeCutOff;
%                             
%                             adjustStimDataPoint = obj.SPSDataPt{SPSStimN}(1) + stimSpikeDataPoint;
%                             
%                             [stimSpikeletAmp, stimSpikeletDataPoint] = obj.sortSpike(tempPatch, stimSpikeAmp, adjustStimDataPoint,obj.spikeCutOff);
%                             
%                             obj.SPSSpikelet{1,SPSStimN} = obj.recordTime(SPSStimN,stimSpikeletDataPoint);
%                             obj.SPSSpikelet{2,SPSStimN} = stimSpikeletAmp;
%                             obj.SPSSpikelet{3,SPSStimN} = stimSpikeletDataPoint;
%                             
%                             stimSpikeDataPoint(stimSpikeAmp < stimCutOff) = [];
%                             stimSpikeAmp(stimSpikeAmp < stimCutOff) = [];
%                             
%                             
%                             obj.SPSSpike{1,SPSStimN} = (obj.SPSDataPt{SPSStimN}(1) + stimSpikeDataPoint - 1)/obj.SR;
%                             obj.SPSSpike{2,SPSStimN} = stimSpikeAmp;
%                             obj.SPSSpike{3,SPSStimN} = obj.SPSDataPt{SPSStimN}(1) + stimSpikeDataPoint - 1;
%                         end
                        
                        
                        %% After stim spikes
                        
                        %                         adjustAStimDataPoint =...
                        %                             obj.AStimDataPoint(SPSStimN,:) - obj.baselineDataPoint(SPSStimN,1) + 1;
                        
%                         obj.AStimDataPoint{SPSStimN} = obj.extStimEndPt(SPSStimN):obj.nSamples;
                        obj.AStimDataPoint{SPSStimN} = (obj.extStimEndPt(SPSStimN) + extSkipPt):obj.nSamples;
                        
                        if SPSStimN == 40
                            SPSStimN
                        end
                        
                        [AStimSpikeAmp, AStimSpikeDataPoint] = ...
                            findpeaks(tempPatch(obj.AStimDataPoint{SPSStimN}), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(minGap*obj.SR));
                        
                        AStimBase = median(tempPatch(obj.AStimDataPoint{SPSStimN}));
                        AStimCutOff = AStimBase + obj.spikeCutOff;
                        
                        adjustAStimDataPoint = obj.AStimDataPoint{SPSStimN}(1) + AStimSpikeDataPoint;
                        [AStimSpikeletAmp, AStimSpikeletDataPoint] = obj.sortSpike(tempPatch, AStimSpikeAmp, adjustAStimDataPoint,obj.spikeCutOff);
                        
                        obj.AStimSpikelet{1,SPSStimN} = obj.recordTime(SPSStimN,AStimSpikeletDataPoint);
                        obj.AStimSpikelet{2,SPSStimN} = AStimSpikeletAmp;
                        obj.AStimSpikelet{3,SPSStimN} = AStimSpikeletDataPoint;
                        
                        AStimSpikeDataPoint(AStimSpikeAmp < AStimCutOff) = [];
                        AStimSpikeAmp(AStimSpikeAmp < AStimCutOff) = [];
                        
                        artefactEntry = find(AStimSpikeAmp > artSpikeAmp);
                        if ~isempty(artefactEntry)
                            AStimSpikeDataPoint(artefactEntry) = [];
                            AStimSpikeAmp(artefactEntry) = [];
                        end
                        
                        obj.AStimSpike{1,SPSStimN} = (obj.AStimDataPoint{SPSStimN}(1) + AStimSpikeDataPoint - 1)/obj.SR;
                        obj.AStimSpike{2,SPSStimN} = AStimSpikeAmp;
                        obj.AStimSpike{3,SPSStimN} = obj.AStimDataPoint{SPSStimN}(1) + AStimSpikeDataPoint - 1;
                        
                        % If AStimSpike is within 10ms after SPS time, then
                        % it is considered as SPS coupled AP. Write the
                        % code
                        
                        if SPSStimN == 17
                            SPSStimN
                        end
                        
                        
                        afterStimFrame = 0.00001;
%                         AStimSpikeDataPoint
                        spikeWithinCoupledRange = find(AStimSpikeDataPoint < coupledTimeFrame*obj.SR);
%                         coupledSpikeAfterStim = find(AStimSpikeDataPoint > afterStimFrame*obj.SR);
                        coupledSpikeAfterStim = spikeWithinCoupledRange;
                        
                        coupledSpikeCheck = ismember(spikeWithinCoupledRange, coupledSpikeAfterStim);
                        coupledSpike = spikeWithinCoupledRange(coupledSpikeCheck);
                        
                        if ~isempty(coupledSpike)
                            obj.SPSSpike{1,SPSStimN} = (obj.AStimDataPoint{SPSStimN}(1) + AStimSpikeDataPoint(coupledSpike) - 1)/obj.SR;
                            obj.SPSSpike{2,SPSStimN} = AStimSpikeAmp(coupledSpike);
                            obj.SPSSpike{3,SPSStimN} = obj.AStimDataPoint{SPSStimN}(1) + AStimSpikeDataPoint(coupledSpike) - 1;
                        end
                        
%                         obj.spikeData{1,SPSStimN} = [obj.baseSpike{1,SPSStimN}, obj.stimSpike{1,SPSStimN}, obj.AStimSpike{1,SPSStimN}];
%                         obj.spikeData{2,SPSStimN} = [obj.baseSpike{2,SPSStimN}, obj.stimSpike{2,SPSStimN}, obj.AStimSpike{2,SPSStimN}];
%                         obj.spikeData{3,SPSStimN} = [obj.baseSpike{3,SPSStimN}, obj.stimSpike{3,SPSStimN}, obj.AStimSpike{3,SPSStimN}];
%                         
%                         obj.spikeletData{1,SPSStimN} = [obj.baseSpikelet{1,SPSStimN}, obj.stimSpikelet{1,SPSStimN}, obj.AStimSpikelet{1,SPSStimN}];
%                         obj.spikeletData{2,SPSStimN} = [obj.baseSpikelet{2,SPSStimN}, obj.stimSpikelet{2,SPSStimN}, obj.AStimSpikelet{2,SPSStimN}];
%                         obj.spikeletData{3,SPSStimN} = [obj.baseSpikelet{3,SPSStimN}, obj.stimSpikelet{3,SPSStimN}, obj.AStimSpikelet{3,SPSStimN}];
                        
                        % 20160206 update for long pulse duration analysis
                        
                        stimMinGap = 1/obj.SR;
                        stimNegDataPt = (obj.extStimStartPt(SPSStimN) + 3):(obj.extStimEndPt(SPSStimN) + 10);
                        stimPosDataPt = (obj.extStimStartPt(SPSStimN) + 3):(obj.extStimEndPt(SPSStimN) + stimExtend);
                        % Extend the last trigger a bit in order to find
                        % the peak
                        
%                         [stimNegSpikeAmp, stimNegSpikeDataPoint] = ...
%                             findpeaks(-tempPatch(stimDataPt), 'MINPEAKHEIGHT', obj.spikeletCutOff, 'MINPEAKDISTANCE', round(stimMinGap*obj.SR));
                        
%                         [stimPosSpikeAmp, stimPosSpikeDataPoint] = ...
%                             findpeaks(tempPatch(stimDataPt));
                        
                        diffStimPosTempPatch = diff(tempPatch(stimPosDataPt));
                        diffStimNegTempPatch = diff(tempPatch(stimNegDataPt));
                        
                        
                        stimPosSpikeDataPoint = [];
                        stimPosSpikeAmp = [];
                        
                        stimNegSpikeDataPoint = [];
                        stimNegSpikeAmp = [];
                        
                        extStimStartPt = obj.extStimStartPt(SPSStimN);
                        
                        obj.stimBaseline(SPSStimN) =...
                            mean(tempPatch((extStimStartPt - 1000):extStimStartPt));
                        
                        for i = 1:(length(diffStimPosTempPatch) - 1)
                            if (diffStimPosTempPatch(i) > 0) &&(diffStimPosTempPatch(i + 1) <= 0)
                                stimPosSpikeDataPoint(end + 1) = stimPosDataPt(i + 1);
                                stimPosSpikeAmp(end + 1) = tempPatch(stimPosDataPt(i + 1));
                            end
                        end
                        
                        for i = 1:(length(diffStimNegTempPatch) - 1)
                            if (diffStimNegTempPatch(i) < 0) &&(diffStimNegTempPatch(i + 1) >= 0)
                                % small bump finding
                                stimNegSpikeDataPoint(end + 1) = stimNegDataPt(i + 1);
                                stimNegSpikeAmp(end + 1) = tempPatch(stimNegDataPt(i + 1));
                            end
                            
                            
                            
                            if (i > 1)
                                lessThanPrev = diffStimNegTempPatch(i - 1) > diffStimNegTempPatch(i);
                                lessThanNext = diffStimNegTempPatch(i + 1) > diffStimNegTempPatch(i);
                                diffBetPrev = diffStimNegTempPatch(i - 1) - diffStimNegTempPatch(i);
                                
                                bumpIncrease = diffStimNegTempPatch(i + 1) > 2*diffStimNegTempPatch(i);
                                
                                
                                if (lessThanPrev)&&(lessThanNext)&&(bumpIncrease)&&(diffStimNegTempPatch(i) > 0)
                                    % small bump finding
                                    stimNegSpikeDataPoint(end + 1) = stimNegDataPt(i + 1);
                                    stimNegSpikeAmp(end + 1) = tempPatch(stimNegDataPt(i + 1));
                                end
                            end
                            
                        end
                        
                        
%                         if ~isempty(stimPosSpikeAmp)
                            obj.stimPosPeak{1,SPSStimN} = (stimPosSpikeDataPoint)/obj.SR;
                            obj.stimPosPeak{2,SPSStimN} = stimPosSpikeAmp;
                            obj.stimPosPeak{3,SPSStimN} = stimPosSpikeDataPoint;
%                         end
% 
%                         [stimNegSpikeAmp, stimNegSpikeDataPoint] = ...
%                             findpeaks(-tempPatch(stimDataPt));
%                         
%                         stimNegSpikeAmp = -stimNegSpikeAmp;
%                         
%                         if ~isempty(stimNegSpikeAmp)
                            obj.stimNegPeak{1,SPSStimN} = (stimNegSpikeDataPoint)/obj.SR;
                            obj.stimNegPeak{2,SPSStimN} = stimNegSpikeAmp;
                            obj.stimNegPeak{3,SPSStimN} = stimNegSpikeDataPoint;
%                         end                        
                    end
                end               
            end
        end
        
        function getSpikeLatency(obj)
            if obj.ndim == 3
                for i = 1:obj.nSweeps
                    obj.spikeLatency{i} = obj.AStimSpike{1,i} - obj.stimStartTime(i);
                end
            else
                
            end
            
        end
        
        function HFSfiltering(obj,samplingRate, stimFreq)
            nData = size(obj.patch,2);
            obj.HFSfilteredSig = zeros(obj.nSweeps, nData);
            for i = 1:size(obj.patch,1)
                
                y = obj.patch(i,:);
                obj.HFSfilteredSig(i,:) = y;
                
                %                 samplingRate = 50000; % Sampling frequency of Clampex data
                %                 T = 1/samplingRate; % Sample time
                %                 nData = 65000; % Length of signal
                %                 t = (0:nData-1)*T; % Time vector
                
                
                
                
                NFFT = 2^nextpow2(nData); % Next power of 2 from length of y
                %                 Y = fft(y,NFFT)/nData;
                obj.FFTFreq = samplingRate/2*linspace(0,1,NFFT/2+1);
                
                for j = 1:5
                    % delete 4 HFS components, do not delete 10000Hz
                    % component because that is APs frequency
                    if (j*stimFreq < 9500)|(j*stimFreq > 10500)
                        [b,a] = ellip (3, .01, 100, [(j*stimFreq - 50), (j*stimFreq + 50)]/(samplingRate/2));
                        N = filter (b, a, y);
                        obj.HFSfilteredSig(i,:) = obj.HFSfilteredSig(i,:) - N;
                    end
                end
                
            end
            
        end
        
        function getWindowSpikeRate(obj,spikeRateWindowSize)
            % windowSize is how long of a duration for counting spiking
            % rate. The smaller the windowSize it is, the finer the spiking
            % rate is. Make windowSize in time, not data point.
            
            obj.nSpikeRatePt = 1000; % make 100 data points for spike rate
            if isempty(obj.stimSpike)
                obj.getStimSpike;
            end
            
            
            
            obj.spikeRateWindowSize = spikeRateWindowSize;
            
            moveStep = round(size(obj.recordTime,2)/obj.nSpikeRatePt);
            % moveStep is moving in data point
            if obj.nSweeps > 1
                for i = 1:obj.nSweeps
                    obj.spikeRate{i} = zeros(1,obj.nSpikeRatePt);
                    
                    for j = 1:(obj.nSpikeRatePt + 1)
                        currLowerLoc = (j - 1)*moveStep - round(spikeRateWindowSize/2) + 1;
                        currHighLoc = (j - 1)*moveStep + round(spikeRateWindowSize/2);
                        % Need to make the window go before and after the
                        % mid point
                        adjustWindowSize = spikeRateWindowSize;
                        if currLowerLoc < 1
                            currLowerLoc = 1;
                            adjustWindowSize = currHighLoc - currLowerLoc;
                        end
                        
                        if currHighLoc > size(obj.recordTime,2)
                            currHighLoc = size(obj.recordTime,2);
                            adjustWindowSize = currHighLoc - currLowerLoc;
                        end
                        stimStartDataPoint = obj.stimDataPoint{1,i}(1);
                        if ((j - 1)*moveStep < stimStartDataPoint)&(currHighLoc > stimStartDataPoint)
                            % Separate the modulation phase from baseline
                            % calculation to eliminate boundary effect
                            currHighLoc = stimStartDataPoint;
                            adjustWindowSize = currHighLoc - currLowerLoc;
                        end
                        if ((j - 1)*moveStep > stimStartDataPoint)&(currLowerLoc < stimStartDataPoint)
                            currLowerLoc = stimStartDataPoint;
                            adjustWindowSize = currHighLoc - currLowerLoc;
                        end
                        
                        % 300ms for modulation phase, so the after
                        % stimulation time starts from 300ms*SR
                        AStimeDataPoint = stimStartDataPoint + 0.3*obj.SR;
                        if ((j - 1)*moveStep > AStimeDataPoint)&(currLowerLoc < AStimeDataPoint)
                            currLowerLoc = AStimeDataPoint;
                            adjustWindowSize = currHighLoc - currLowerLoc;
                        end
                        if ((j - 1)*moveStep < AStimeDataPoint)&(currHighLoc > AStimeDataPoint)
                            currHighLoc = AStimeDataPoint;
                            adjustWindowSize = currHighLoc - currLowerLoc;
                        end
                        
                        windowSpikes = find((obj.spikeData{3,i} >= currLowerLoc)&(obj.spikeData{3,i} <= currHighLoc));
                        obj.spikeRate{i}(j) = length(windowSpikes)/spikeRateWindowSize*obj.SR;
                        
                        obj.spikeRateTime{i}(j) = obj.recordTime(i,currLowerLoc);
                        
                        
                    end
                    
                end
                
            else
                for i = 1:size(obj.spikeData,2)
                    for j = 1:obj.nSpikeRatePt
                        currLowerLoc = (j - 1)*moveStep + 1;
                        currHighLoc = (j - 1)*moveStep + spikeRateWindowSize;
                        windowSpikes = find((obj.spikeData{3,i} >= currLowerLoc)&(obj.spikeData{3,i} <= currHighLoc));
                        obj.spikeRate{i}(j) = length(windowSpikes)/spikeRateWindowSize*obj.SR;
                        
                        obj.spikeRateTime{i}(j) = obj.dataTime(currLowerLoc);
                    end
                end
            end
            
            
        end
        

        
        function plotSpikeRate(obj, varargin)
            
            obj.getSpikeRate;
            
            if ~isempty(strfind(obj.protocol, 'IC AP family'))
                obj.plotIntSpikeRate(varargin);                
            elseif ~isempty(strfind(obj.protocol, 'LED'))
                obj.plotLightSpikeRate(varargin);                
            elseif ~isempty(strfind(obj.protocol, 'SPS'))
                obj.plotSPSSpikeRate(varargin);                
            elseif ~isempty(strfind(obj.protocol, 'HFS'))
                obj.plotHFSSpikeRate(varargin);                
            end
            
            
        end
        
        function plotIntSpikeRate(obj, varargin)
            % Simply plot the spike rate of each sweep at different phases
            % (baseline, stimulation and after stimulation).
            
            plotSweep = 0;
            
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    if strcmp(varargin{1}{i}, 'sweepN')
                        plotSweep = varargin{1}{i + 1};
                    end
                end
            end
            
            obj.getSpikeRate;
            
            if plotSweep == 0
                figure;
                for sweepN = 1:obj.nSweeps
                    subplot(3,ceil(obj.nSweeps/3),sweepN);
%                     plot(obj.recordTime(sweepN,:),obj.patch(sweepN,:),'.k');
                    hold on
                    if ~isempty(obj.baseSpike{1,sweepN})
                        
                        plot(obj.recordTime(sweepN,1) - 1/obj.SR + obj.baseSpike{1,sweepN},obj.baselineIndSpikeRate{sweepN},'.b');
                        
                        baselineTimeStart = obj.recordTime(sweepN,1) + (obj.baselineDataPoint{sweepN}(1) - 1)/obj.SR;
                        baselineTimeEnd = obj.recordTime(sweepN,1) + (obj.baselineDataPoint{sweepN}(end) - 1)/obj.SR;
                        
                        line([baselineTimeStart, baselineTimeEnd],...
                            [obj.baselineGeneSpikeRate(sweepN), obj.baselineGeneSpikeRate(sweepN)], 'color', 'b')

                    end
                    
                    if ~isempty(obj.stimSpike{1,sweepN})
                        
                        plot(obj.recordTime(sweepN,1) - 1/obj.SR + obj.stimSpike{1,sweepN},obj.stimIndSpikeRate{sweepN},'.r');
                        
                        stimTimeStart = obj.recordTime(sweepN,1) + (obj.stimDataPoint{sweepN}(1) - 1)/obj.SR;
                        stimTimeEnd = obj.recordTime(sweepN,1) + (obj.stimDataPoint{sweepN}(end) - 1)/obj.SR;
                        
                        line([stimTimeStart, stimTimeEnd],...
                            [obj.stimGeneSpikeRate(sweepN), obj.stimGeneSpikeRate(sweepN)], 'color', 'r')
                    end
                    
                    if ~isempty(obj.AStimSpike{1,sweepN})
                        plot(obj.recordTime(sweepN,1) - 1/obj.SR + obj.AStimSpike{1,sweepN},obj.AStimIndSpikeRate{sweepN},'.k');
                        
                        AStimTimeStart = obj.recordTime(sweepN,1) + (obj.AStimDataPoint{sweepN}(1) - 1)/obj.SR;
                        AStimTimeEnd = obj.recordTime(sweepN,1) + (obj.AStimDataPoint{sweepN}(end) - 1)/obj.SR;
                        
                        line([AStimTimeStart, AStimTimeEnd],...
                            [obj.AStimGeneSpikeRate(sweepN), obj.AStimGeneSpikeRate(sweepN)], 'color', 'k')
                    end
                    
                    
                    
                    %                     if ~isempty(obj.spikeletData{1,i})
                    %                         hold on
                    %                         plot(obj.spikeletData{1,i},obj.spikeletData{2,i},'.g');
                    %                     end
                    
                    if ~isempty(obj.intStimAmp)
                        titleString = sprintf('Sweep #%d, current = %3.1fpA', sweepN, obj.intStimAmp(sweepN));
                        title(titleString, 'FontSize', 8);
                    end
                    xlabel('Time (s)','FontSize', 8);
                    ylabel('Spike rate (Hz)','FontSize', 8);
                    
                    maxGeneSpikeRate = max([obj.baselineGeneSpikeRate(sweepN),...
                        obj.stimGeneSpikeRate(sweepN),...
                        obj.AStimGeneSpikeRate(sweepN)]);
                    
                    xlim([obj.recordTime(sweepN,1), obj.recordTime(sweepN,end)]);
                    ylim([0, maxGeneSpikeRate*2]);
                    
                end
            else
                figure
                                    hold on
                    if ~isempty(obj.baseSpike{1,plotSweep})
                        
                        plot(obj.recordTime(plotSweep,1) - 1/obj.SR + obj.baseSpike{1,plotSweep},obj.baselineIndSpikeRate{plotSweep},'.b');
                        
                        baselineTimeStart = obj.recordTime(plotSweep,1) + (obj.baselineDataPoint{plotSweep}(1) - 1)/obj.SR;
                        baselineTimeEnd = obj.recordTime(plotSweep,1) + (obj.baselineDataPoint{plotSweep}(end) - 1)/obj.SR;
                        
                        line([baselineTimeStart, baselineTimeEnd],...
                            [obj.baselineGeneSpikeRate(plotSweep), obj.baselineGeneSpikeRate(plotSweep)], 'color', 'b')

                    end
                    
                    if ~isempty(obj.stimSpike{1,plotSweep})
                        
                        plot(obj.recordTime(plotSweep,1) - 1/obj.SR + obj.stimSpike{1,plotSweep},obj.stimIndSpikeRate{plotSweep},'.r');
                        
                        stimTimeStart = obj.recordTime(plotSweep,1) + (obj.stimDataPoint{plotSweep}(1) - 1)/obj.SR;
                        stimTimeEnd = obj.recordTime(plotSweep,1) + (obj.stimDataPoint{plotSweep}(end) - 1)/obj.SR;
                        
                        line([stimTimeStart, stimTimeEnd],...
                            [obj.stimGeneSpikeRate(plotSweep), obj.stimGeneSpikeRate(plotSweep)], 'color', 'r')
                    end
                    
                    if ~isempty(obj.AStimSpike{1,plotSweep})
                        plot(obj.recordTime(plotSweep,1) - 1/obj.SR + obj.AStimSpike{1,plotSweep},obj.AStimIndSpikeRate{plotSweep},'.k');
                        
                        AStimTimeStart = obj.recordTime(plotSweep,1) + (obj.AStimDataPoint{plotSweep}(1) - 1)/obj.SR;
                        AStimTimeEnd = obj.recordTime(plotSweep,1) + (obj.AStimDataPoint{plotSweep}(end) - 1)/obj.SR;
                        
                        line([AStimTimeStart, AStimTimeEnd],...
                            [obj.AStimGeneSpikeRate(plotSweep), obj.AStimGeneSpikeRate(plotSweep)], 'color', 'k')
                    end
                    
                    maxGeneSpikeRate = max([obj.baselineGeneSpikeRate(plotSweep),...
                        obj.stimGeneSpikeRate(plotSweep),...
                        obj.AStimGeneSpikeRate(plotSweep)]);
                    
                    %                     if ~isempty(obj.spikeletData{1,i})
                    %                         hold on
                    %                         plot(obj.spikeletData{1,i},obj.spikeletData{2,i},'.g');
                    %                     end
                    
                    if ~isempty(obj.intStimAmp)
                        titleString = sprintf('Sweep #%d, current = %3.1fpA', plotSweep, obj.intStimAmp(plotSweep));
                        title(titleString, 'FontSize', 8);
                    end
                    xlabel('Time (s)','FontSize', 8);
                    ylabel('Spike rate (Hz)','FontSize', 8);    
                    
                    xlim([obj.recordTime(plotSweep,1), obj.recordTime(plotSweep,end)]);
                    ylim([0, maxGeneSpikeRate]);
                
            end
        end
        
        function plotIntSpikeRateOverview(obj, varargin)
            % Plot the intracellular spike rates (number of spikes in 0.5s) against
            % 1. Intracellular current injected
            % 2. The membrane voltage reached during current injection
            % 3. The membrane difference between 2. and resting Vm
            % Each sweep will be a point, and this function will plot all
            % the sweeps to characterise the relationship between current
            % injection amount and spiking rate.
            
            plotPhase = 'stim';
            % can plot for baseline, stimulation or after stimulation, the
            % default is to plot the stimulation phase
            
            plotAgainst = 'IC';
            
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    if strcmp(varargin{i}, 'phase')
                        plotPhase = varargin{1}{i + 1};
                        % 'baseline'
                        % 'stim'
                        % 'AStim'
                    end
                    if strcmp(varargin{i}, 'against')
                        plotAgainst = varargin{i + 1};
                        % 'IC' = intracellular current
                        % 'Vm' = membrane potential reached
                        % 'VmD' = different in membrane voltage and resting membrane voltage
                    end
                end
            end
            
            obj.getSpikeRate;
            
            spikeNum = zeros(1,obj.nSweeps);
            spikeRate = zeros(1,obj.nSweeps);
            
            xVal = zeros(1,obj.nSweeps);
            
            for sweepI = 1:obj.nSweeps
                % Get the spike rate at desired phase
                if strcmp(plotPhase, 'baseline')
                    spikeNum = length(obj.baseSpike{1, sweepI});
                    baselineTimeStart = (obj.baselineDataPoint{sweepI}(1) - 1)/obj.SR;
                    baselineTimeEnd = (obj.baselineDataPoint{sweepI}(end) - 1)/obj.SR;
                    
                    baseTime = baselineTimeEnd - baselineTimeStart;
                    
                    spikeRate(sweepI) = spikeNum/baseTime;
                    
                    
                elseif strcmp(plotPhase, 'stim')
                    spikeNum = length(obj.stimSpike{1, sweepI});
                    
                    stimTimeStart = (obj.stimDataPoint{sweepI}(1) - 1)/obj.SR;
                    stimTimeEnd = (obj.stimDataPoint{sweepI}(end) - 1)/obj.SR;
                    
                    stimTime = stimTimeEnd - stimTimeStart;
                    
                    spikeRate(sweepI) = spikeNum/stimTime;
                    
                elseif strcmp(plotPhase, 'AStim')
                    spikeNum = length(obj.AStimSpike{1, sweepI});
                    
                    AStimTimeStart = (obj.AStimDataPoint{sweepI}(1) - 1)/obj.SR;
                    AStimTimeEnd = (obj.AStimDataPoint{sweepI}(end) - 1)/obj.SR;
                    
                    AStimTime = AStimTimeEnd - AStimTimeStart;
                    
                    spikeRate(sweepI) = spikeNum/AStimTime;
                end
                
                baseStartDataPt = obj.baselineDataPoint{sweepI}(1);
                baseEndDataPt = obj.baselineDataPoint{sweepI}(end);
                
                baseVmMedian = median(obj.patch(sweepI, baseStartDataPt:baseEndDataPt));
                
                stimStartDataPt = obj.stimDataPoint{sweepI}(1);
                stimEndDataPt = obj.stimDataPoint{sweepI}(end);
                
                stimVmMedian = median(obj.patch(sweepI,stimStartDataPt:stimEndDataPt));
                
                if strcmp(plotAgainst, 'IC')
                    xVal(sweepI) = obj.intStimAmp(sweepI);
                    if sweepI == 1
                        baseVmMedian
                    end
                    
                elseif strcmp(plotAgainst, 'Vm')
                    xVal(sweepI) = stimVmMedian;
                    
                elseif strcmp(plotAgainst, 'VmD')
                    xVal(sweepI) = stimVmMedian - baseVmMedian;
                    
                    
                end
                
                % Get the properties to be plotted against
            end
            
            if strcmp(plotPhase, 'baseline')
                yString = 'Spike rate (Hz) during baseline';
            elseif strcmp(plotPhase, 'stim')
                yString = 'Spike rate (Hz) during stim';
            elseif strcmp(plotPhase, 'AStim')
                yString = 'Spike rate (Hz) during AStim';
            end
            
            if strcmp(plotAgainst, 'IC')
                xString = 'Intracellular injected current (\muA)';
            elseif strcmp(plotAgainst, 'Vm')
                xString = 'Membrane potential reached during current injection (mV)';
            elseif strcmp(plotAgainst, 'VmD')
                xString = 'Membrane potential difference to resting membrane potential (mV)';
            end
            
            
            plot(xVal, spikeRate, '.r', 'MarkerSize', 30);
             
            xlabel(xString, 'FontSize', 20);
            ylabel(yString, 'FontSize', 20);
            
               
        end

       
    
    
    function plotLightSpikeRate(obj, varargin)
            plotSweep = 0;
            
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    if strcmp(varargin{1}{i}, 'sweepN')
                        plotSweep = varargin{1}{i + 1};
                    end
                end
            end
            
            obj.getSpikeRate;
            
            figure
            
            %             if plotSweep == 0
            for sweepN = 1:size(obj.LEDONSpike, 2)
                if plotSweep ~= 0
                    sweepN = plotSweep;
                end
                
                
                
                if (obj.ndim == 3)&(plotSweep == 0)
                    subplot(3,ceil(obj.nSweeps/3),sweepN);
                end
                %                     plot(obj.recordTime(sweepN,:),obj.patch(sweepN,:),'.k');
                hold on
                if ~isempty(obj.baseSpike{1,sweepN})
                    
                    plot(obj.recordTime(sweepN,1) - 1/obj.SR + obj.baseSpike{1,sweepN},obj.baselineIndSpikeRate{sweepN},'.b');
                    
                    baselineTimeStart = obj.recordTime(sweepN,1) + (obj.baselineDataPoint{sweepN}(1) - 1)/obj.SR;
                    baselineTimeEnd = obj.recordTime(sweepN,1) + (obj.baselineDataPoint{sweepN}(end) - 1)/obj.SR;
                    
                    line([baselineTimeStart, baselineTimeEnd],...
                        [obj.baselineGeneSpikeRate(sweepN), obj.baselineGeneSpikeRate(sweepN)], 'color', 'b')
                    
                end
                
                if ~isempty(obj.LEDONSpike{1,sweepN})
                    for LEDN = 1:size(obj.LEDONSpike,3)
                        plot(obj.recordTime(sweepN,1) - 1/obj.SR + obj.LEDONSpike{1,sweepN, LEDN},obj.LEDONIndSpikeRate{sweepN, LEDN},'.r');
                        
                        LEDONTimeStart = obj.recordTime(sweepN,1) + (obj.LEDONDataPoint{sweepN, LEDN}(1) - 1)/obj.SR;
                        LEDONTimeEnd = obj.recordTime(sweepN,1) + (obj.LEDONDataPoint{sweepN, LEDN}(end) - 1)/obj.SR;
                        
                        line([LEDONTimeStart, LEDONTimeEnd],...
                            [obj.LEDONGeneSpikeRate(sweepN, LEDN), obj.LEDONGeneSpikeRate(sweepN, LEDN)], 'color', 'r')
                    end
                end
                
                if ~isempty(obj.LEDOFFSpike{1,sweepN})
                    for LEDN = 1:size(obj.LEDONSpike,3)
                        plot(obj.recordTime(sweepN,1) - 1/obj.SR + obj.LEDOFFSpike{1,sweepN, LEDN},obj.LEDOFFIndSpikeRate{sweepN, LEDN},'.k');
                        
                        LEDOFFTimeStart = obj.recordTime(sweepN,1) + (obj.LEDOFFDataPoint{sweepN, LEDN}(1) - 1)/obj.SR;
                        LEDOFFTimeEnd = obj.recordTime(sweepN,1) + (obj.LEDOFFDataPoint{sweepN, LEDN}(end) - 1)/obj.SR;
                        
                        line([LEDOFFTimeStart, LEDOFFTimeEnd],...
                            [obj.LEDOFFGeneSpikeRate(sweepN, LEDN), obj.LEDOFFGeneSpikeRate(sweepN, LEDN)], 'color', 'k')
                    end
                end
                
                if ~isempty(obj.AStimSpike{1,sweepN})
                    plot(obj.recordTime(sweepN,1) - 1/obj.SR + obj.AStimSpike{1,sweepN},obj.AStimIndSpikeRate{sweepN},'.m');
                    
                    AStimTimeStart = obj.recordTime(sweepN,1) + (obj.AStimDataPoint{sweepN}(1) - 1)/obj.SR;
                    AStimTimeEnd = obj.recordTime(sweepN,1) + (obj.AStimDataPoint{sweepN}(end) - 1)/obj.SR;
                    
                    line([AStimTimeStart, AStimTimeEnd],...
                        [obj.AStimGeneSpikeRate(sweepN), obj.AStimGeneSpikeRate(sweepN)], 'color', 'm')
                end
                
                
                
                %                     if ~isempty(obj.spikeletData{1,i})
                %                         hold on
                %                         plot(obj.spikeletData{1,i},obj.spikeletData{2,i},'.g');
                %                     end
                
                if ~isempty(obj.intStimAmp)
                    titleString = sprintf('Sweep #%d, LED stimulation', sweepN);
                    title(titleString, 'FontSize', 8);
                end
                xlabel('Time (s)','FontSize', 8);
                ylabel('Spike rate (Hz)','FontSize', 8);
                
                if plotSweep ~= 0
                    break;
                end
                
                maxGeneSpikeRate = max([obj.baselineGeneSpikeRate(sweepN),...
                    obj.LEDOFFGeneSpikeRate(sweepN),...
                    obj.LEDONGeneSpikeRate(sweepN),...
                    obj.AStimGeneSpikeRate(sweepN)]);
                
                xlim([obj.recordTime(sweepN,1), obj.recordTime(sweepN,end)]);
                ylim([0, maxGeneSpikeRate*3]);

            
                
                %                 end
                %             else
                %                 figure
                %                 hold on
                %                 if ~isempty(obj.baseSpike{1,plotSweep})
                %
                %                     plot(obj.recordTime(plotSweep,1) - 1/obj.SR + obj.baseSpike{1,plotSweep},obj.baselineIndSpikeRate{plotSweep},'.b');
                %
                %                     baselineTimeStart = obj.recordTime(plotSweep,1) + (obj.baselineDataPoint{plotSweep}(1) - 1)/obj.SR;
                %                     baselineTimeEnd = obj.recordTime(plotSweep,1) + (obj.baselineDataPoint{plotSweep}(end) - 1)/obj.SR;
                %
                %                     line([baselineTimeStart, baselineTimeEnd],...
                %                         [obj.baselineGeneSpikeRate(plotSweep), obj.baselineGeneSpikeRate(plotSweep)], 'color', 'b')
                %
                %                 end
                %
                %                 if ~isempty(obj.LEDONSpike{1,plotSweep})
                %                     for LEDN = 1:size(obj.LEDONSpike,3)
                %                         plot(obj.recordTime(plotSweep,1) - 1/obj.SR + obj.LEDONSpike{1,plotSweep, LEDN},obj.LEDONIndSpikeRate{plotSweep, LEDN},'.r');
                %
                %                         LEDONTimeStart = obj.recordTime(plotSweep,1) + (obj.LEDONDataPoint{plotSweep, LEDN}(1) - 1)/obj.SR;
                %                         LEDONTimeEnd = obj.recordTime(plotSweep,1) + (obj.LEDONDataPoint{plotSweep, LEDN}(end) - 1)/obj.SR;
                %
                %                         line([LEDONTimeStart, LEDONTimeEnd],...
                %                             [obj.LEDONGeneSpikeRate(plotSweep, LEDN), obj.LEDONGeneSpikeRate(plotSweep, LEDN)], 'color', 'r')
                %                     end
                %                 end
                %
                %                 if ~isempty(obj.LEDOFFSpike{1,plotSweep})
                %                     for LEDN = 1:size(obj.LEDONSpike,3)
                %                         plot(obj.recordTime(plotSweep,1) - 1/obj.SR + obj.LEDOFFSpike{1,plotSweep, LEDN},obj.LEDOFFIndSpikeRate{plotSweep, LEDN},'.k');
                %
                %                         LEDOFFTimeStart = obj.recordTime(plotSweep,1) + (obj.LEDOFFDataPoint{plotSweep, LEDN}(1) - 1)/obj.SR;
                %                         LEDOFFTimeEnd = obj.recordTime(plotSweep,1) + (obj.LEDOFFDataPoint{plotSweep, LEDN}(end) - 1)/obj.SR;
                %
                %                         line([LEDOFFTimeStart, LEDOFFTimeEnd],...
                %                             [obj.LEDOFFGeneSpikeRate(plotSweep, LEDN), obj.LEDOFFGeneSpikeRate(plotSweep, LEDN)], 'color', 'k')
                %                     end
                %                 end
                %
                %                 if ~isempty(obj.AStimSpike{1,plotSweep})
                %                     plot(obj.recordTime(plotSweep,1) - 1/obj.SR + obj.AStimSpike{1,plotSweep},obj.AStimIndSpikeRate{plotSweep},'.k');
                %
                %                     AStimTimeStart = obj.recordTime(plotSweep,1) + (obj.AStimDataPoint{plotSweep}(1) - 1)/obj.SR;
                %                     AStimTimeEnd = obj.recordTime(plotSweep,1) + (obj.AStimDataPoint{plotSweep}(end) - 1)/obj.SR;
                %
                %                     line([AStimTimeStart, AStimTimeEnd],...
                %                         [obj.AStimGeneSpikeRate(plotSweep), obj.AStimGeneSpikeRate(plotSweep)], 'color', 'k')
                %                 end
                %
                %
                %
                %                 %                     if ~isempty(obj.spikeletData{1,i})
                %                 %                         hold on
                %                 %                         plot(obj.spikeletData{1,i},obj.spikeletData{2,i},'.g');
                %                 %                     end
                %
                %                 if ~isempty(obj.intStimAmp)
                %                     titleString = sprintf('Sweep #%d, current = %3.1fpA', plotSweep, obj.intStimAmp(plotSweep));
                %                     title(titleString, 'FontSize', 8);
                %                 end
                %                 xlabel('Time (s)','FontSize', 8);
                %                 ylabel('Spike rate (Hz)','FontSize', 8);
                %
            end
        end
            

        
        
        function plotSPSSpikeRate(obj, varargin)
            
        end
        
        
        function plotHFSSpikeRate(obj, varargin)
            plotSweep = 0;
            HFSAmp = [];
            
            subplotCol = 11;
            subplotRow = 5;
            
            
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    if strcmp(varargin{1}{i}, 'sweepN')
                        plotSweep = varargin{1}{i + 1};
                    end
                    
                    if strcmp(varargin{1}{i}, 'HFSAmp')
                        HFSAmp = varargin{1}{i + 1};
                    end
                end
            end
            
            obj.getSpikeRate;
%             figure
            
            if obj.ndim == 3
%                 figure
                for sweepN = 1:size(obj.HFSSpike, 2)
                    if plotSweep ~= 0
                        sweepN = plotSweep;
                    end
                    
                    if obj.extStimStartPt(sweepN) == 0
                        continue
                    end
                    
                    if (obj.ndim == 3)&(plotSweep == 0)
                        sweepMod = mod(sweepN, subplotCol*subplotRow);
                        if sweepMod == 0
                            figure
                            sweepMod = 1;
                        end
                        subplot(subplotRow,subplotCol,sweepMod);

                    end
                    %                     plot(obj.recordTime(sweepN,:),obj.patch(sweepN,:),'.k');
                    hold on
                    if ~isempty(obj.baseSpike{1,sweepN})
                        
                        plot(obj.recordTime(sweepN,1) - 1/obj.SR + obj.baseSpike{1,sweepN}',obj.baselineIndSpikeRate{sweepN},'.b');
                        
                        baselineTimeStart = obj.recordTime(sweepN,1) + (obj.baselineDataPoint{sweepN}(1) - 1)/obj.SR;
                        baselineTimeEnd = obj.recordTime(sweepN,1) + (obj.baselineDataPoint{sweepN}(end) - 1)/obj.SR;
                        
                        line([baselineTimeStart, baselineTimeEnd],...
                            [obj.baselineGeneSpikeRate(sweepN), obj.baselineGeneSpikeRate(sweepN)], 'color', 'b')
                        
                    end
                    
                    if ~isempty(obj.HFSSpike{1,sweepN})
                        
                        plot(obj.recordTime(sweepN,1) - 1/obj.SR + obj.HFSSpike{1,sweepN}',obj.HFSIndSpikeRate{sweepN},'.r');
                        
                        HFSTimeStart = obj.recordTime(sweepN,1) + (obj.HFSDataPt{sweepN}(1) - 1)/obj.SR;
                        HFSTimeEnd = obj.recordTime(sweepN,1) + (obj.HFSDataPt{sweepN}(end) - 1)/obj.SR;
                        
                        line([HFSTimeStart, HFSTimeEnd],...
                            [obj.HFSGeneSpikeRate(sweepN), obj.HFSGeneSpikeRate(sweepN)], 'color', 'r');
                        
%                         plot(obj.recordTime(sweepN,1) - 1/obj.SR + obj.HFSGapSpike{1,sweepN}',obj.HFSGapIndSpikeRate{sweepN},'.k');
%                         
%                         HFSGapTimeStart = obj.recordTime(sweepN,1) + (obj.HFSGapDataPoint{sweepN}(1) - 1)/obj.SR;
%                         HFSGapTimeEnd = obj.recordTime(sweepN,1) + (obj.HFSGapDataPoint{sweepN}(end) - 1)/obj.SR;
%                         
%                         line([HFSGapTimeStart, HFSGapTimeEnd],...
%                             [obj.HFSGapGeneSpikeRate(sweepN), obj.HFSGapGeneSpikeRate(sweepN)], 'color', 'k')
                        
                    end
                    
                    
                    if ~isempty(obj.AStimSpike{1,sweepN})
                        plot(obj.recordTime(sweepN,1) - 1/obj.SR + obj.AStimSpike{1,sweepN}',obj.AStimIndSpikeRate{sweepN},'.m');
                        
                        AStimTimeStart = obj.recordTime(sweepN,1) + (obj.AStimDataPoint{sweepN}(1) - 1)/obj.SR;
                        AStimTimeEnd = obj.recordTime(sweepN,1) + (obj.AStimDataPoint{sweepN}(end) - 1)/obj.SR;
                        
                        line([AStimTimeStart, AStimTimeEnd],...
                            [obj.AStimGeneSpikeRate(sweepN), obj.AStimGeneSpikeRate(sweepN)], 'color', 'm')
                    end
                    
                    
                    
%                     if ~isempty(obj.intStimAmp)
                    if ~isempty(HFSAmp)
                        titleString = sprintf('Sweep #%d, HFS Amp: %d\\muA', sweepN, HFSAmp(sweepN));
                        title(titleString, 'FontSize', 8);
                    end
                    xlabel('Time (s)','FontSize', 8);
                    ylabel('Spike rate (Hz)','FontSize', 8);
                    
                    if plotSweep ~= 0
                        break;
                    end
                    
                    
                    maxGeneSpikeRate = max([obj.baselineGeneSpikeRate(sweepN),...
                        obj.HFSGeneSpikeRate(sweepN),...
                        obj.AStimGeneSpikeRate(sweepN)]);
                    
                    xlim([obj.recordTime(sweepN,1), obj.recordTime(sweepN,end)]);
                    if maxGeneSpikeRate ~= 0
                        ylim([0, maxGeneSpikeRate*3]);
                    end
                end
            else
%                 figure
%                 hold on
                
                if ~isempty(obj.baseSpike{1})
                    
                    plot(obj.dataTime(1) - 1/obj.SR + obj.baseSpike{1}',obj.baselineIndSpikeRate{1},'.b');
                    
                    baselineTimeStart = obj.dataTime(1) + (obj.baselineDataPoint{1}(1) - 1)/obj.SR;
                    baselineTimeEnd = obj.dataTime(1) + (obj.baselineDataPoint{1}(end) - 1)/obj.SR;
                    
                    line([baselineTimeStart, baselineTimeEnd],...
                        [obj.baselineGeneSpikeRate(1), obj.baselineGeneSpikeRate(1)], 'color', 'b')
                    
                end
                
                if ~isempty(obj.AStimSpike{1})
                    plot(obj.dataTime(1) - 1/obj.SR + obj.AStimSpike{1}',obj.AStimIndSpikeRate{1},'.m');
                    
                    AStimTimeStart = obj.dataTime(1) + (obj.AStimDataPoint{1}(1) - 1)/obj.SR;
                    AStimTimeEnd = obj.dataTime(1) + (obj.AStimDataPoint{1}(end) - 1)/obj.SR;
                    
                    line([AStimTimeStart, AStimTimeEnd],...
                        [obj.AStimGeneSpikeRate(1), obj.AStimGeneSpikeRate(1)], 'color', 'm')
                end
                
                for HFSN = 1:size(obj.HFSSpike, 2)

                    
                    if ~isempty(obj.HFSSpike{1,HFSN})
                        
                        plot(obj.dataTime(1) - 1/obj.SR + obj.HFSSpike{1,HFSN}',obj.HFSIndSpikeRate{HFSN},'.r');
                        
                        HFSTimeStart = obj.dataTime(1) + (obj.HFSDataPt{HFSN}(1) - 1)/obj.SR;
                        HFSTimeEnd = obj.dataTime(1) + (obj.HFSDataPt{HFSN}(end) - 1)/obj.SR;
                        
                        line([HFSTimeStart, HFSTimeEnd],...
                            [obj.HFSGeneSpikeRate(HFSN), obj.HFSGeneSpikeRate(HFSN)], 'color', 'r');
                        
                        if HFSN ~= size(obj.HFSSpike, 2)
                            plot(obj.dataTime(1) - 1/obj.SR + obj.HFSGapSpike{1,HFSN}',obj.HFSGapIndSpikeRate{HFSN},'.k');
                            
                            HFSGapTimeStart = obj.dataTime(1) + (obj.HFSGapDataPoint{HFSN}(1) - 1)/obj.SR;
                            HFSGapTimeEnd = obj.dataTime(1) + (obj.HFSGapDataPoint{HFSN}(end) - 1)/obj.SR;
                            
                            line([HFSGapTimeStart, HFSGapTimeEnd],...
                                [obj.HFSGapGeneSpikeRate(HFSN), obj.HFSGapGeneSpikeRate(HFSN)], 'color', 'k')
                        end
                        
                    end
                    
                    if ~isempty(HFSAmp)
                        titleString = sprintf('Sweep #%d, HFS Amp: %d\\muA', HFSAmp);
                        title(titleString, 'FontSize', 8);
                    end
                    xlabel('Time (s)','FontSize', 8);
                    ylabel('Spike rate (Hz)','FontSize', 8);
                    
                    if plotSweep ~= 0
                        break;
                    end
                    
                    maxGeneSpikeRate = max([obj.baselineGeneSpikeRate(1),...
                        obj.HFSGeneSpikeRate(HFSN),...
                        obj.AStimGeneSpikeRate(1)]);
                    
                    xlim([obj.dataTime(1,1), obj.dataTime(1,end)]);
                    ylim([0, maxGeneSpikeRate*3]);
                    
                end
                
            end
        end
            
        
        
        function [spikeletAmp, spikeletDataPoint] = sortSpike(obj,tempPatch,spikeAmp, spikeDataPoint,spikeCutOff)
            % sorting spikes into spikelet or spikes, also eliminate false
            % positive spikelets
            
            spikeletAmp = spikeAmp(spikeAmp < spikeCutOff);
            spikeletDataPoint = spikeDataPoint(spikeAmp < spikeCutOff);
            
            testSpikeletInterval1 = spikeletDataPoint + 100;
            testSpikeletInterval2 = spikeletDataPoint - 100;
            testSpikeletInterval2(testSpikeletInterval2 <= 0) = 1;
            testSpikeletInterval1(testSpikeletInterval1 >= length(tempPatch)) = length(tempPatch);
            testRealSpikelet = ...
                ((spikeletAmp - tempPatch(testSpikeletInterval1)) > 0.2)&((spikeletAmp - tempPatch(testSpikeletInterval2)) > 0.2);
            % if the difference between a spikelet and its
            % previous 10 data and next point is not bigger than 0.1mV,
            % then it is a fake spikelet (false positive from
            % the findpeaks function)
            spikeletAmp = spikeletAmp(testRealSpikelet);
            spikeletDataPoint = spikeletDataPoint(testRealSpikelet);
            
            
            
            
        end
        
        function saveImage(obj, saveName)
            h = gcf;
            set(h,'position',get(0,'screensize'))
            [pathstr,name,ext] = fileparts(obj.path);
            saveFullName = fullfile(pathstr,saveName);
            saveas(h,saveFullName,'tif');
        end
        
        function plotVCInputR(obj)
            
            plot(obj.dataTime(:),obj.patch(:),'b');
            title('VC input resistance');
            xlabel('Time (s)');
            ylabel('I (pA)');
            
            legend('VC recording', 'Location','NorthWest');
            
        end
        
        
        function plotSpike(obj, varargin)
            
            
            if ~isempty(strfind(obj.protocol, 'IC AP family'))
                obj.plotIntSpike(varargin);
                
            elseif ~isempty(strfind(obj.protocol, 'LED'))
                obj.plotLightSpike(varargin);
                
            elseif ~isempty(strfind(obj.protocol, 'SPS'))
                obj.plotSPSSpike(varargin);
                
            elseif ~isempty(strfind(obj.protocol, 'HFS'))
                obj.plotHFSSpike(varargin);
                
            end
            
            
        end
        
        function plotIntSpike(obj, varargin)
            
            plotSweep = 0;
            
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    if strcmp(varargin{1}{i}, 'sweepN')
                        plotSweep = varargin{1}{i + 1};
                    end
                    
                end
            end
            
%             if isempty(obj.spikeData)
%                 obj.getSpikes;
%             end
            figure;
            
            if plotSweep == 0
                for i = 1:obj.nSweeps
                    subplot(3,ceil(obj.nSweeps/3),i);
                    plot(obj.recordTime(i,:),obj.patch(i,:),'.k');
                    
                    if ~isempty(obj.baseSpike{1,i})
                        hold on
                        plot(obj.recordTime(i,1) - 1/obj.SR + obj.baseSpike{1,i},obj.baseSpike{2,i},'.c');
                        %                     else
                        %                         displayString = sprintf('No spikes detected for Sweep #%d', i);
                        %                         display(displayString);
                    end
                    
                    if ~isempty(obj.stimSpike{1,i})
                        hold on
                        plot(obj.recordTime(i,1) - 1/obj.SR + obj.stimSpike{1,i},obj.stimSpike{2,i},'.r');
                    end
                    
                    if ~isempty(obj.AStimSpike{1,i})
                        hold on
                        plot(obj.recordTime(i,1) - 1/obj.SR + obj.AStimSpike{1,i},obj.AStimSpike{2,i},'.m');
                    end
                    
                    
                    
                    %                     if ~isempty(obj.spikeletData{1,i})
                    %                         hold on
                    %                         plot(obj.spikeletData{1,i},obj.spikeletData{2,i},'.g');
                    %                     end
                    
                    if ~isempty(obj.intStimAmp)
                        titleString = sprintf('Sweep #%d, I = %3.0fpA', i, obj.intStimAmp(i));
                        title(titleString, 'FontSize', 10,'FontWeight', 'bold');
                    end
                    xlabel('Time (s)','FontSize', 8,'FontWeight', 'bold');
                    ylabel('Vm (mV)','FontSize', 8,'FontWeight', 'bold');
                    
                end
            else
                plot(obj.recordTime(plotSweep,:),obj.patch(plotSweep,:),'.b');
                
                if ~isempty(obj.spikeData{1,plotSweep})
                    hold on
                    plot(obj.spikeData{1,plotSweep},obj.spikeData{2,plotSweep},'.r');
                else
                    %                     displayString = sprintf('No spikes detected for Sweep #%d', plotType);
                    %                     display(displayString);
                end
                
                if ~isempty(obj.spikeletData{1,plotSweep})
                    hold on
                    plot(obj.spikeletData{1,plotSweep},obj.spikeletData{2,plotSweep},'.g');
                end
                
                if ~isempty(obj.intStimAmp)
                    titleString = sprintf('Sweep #%d, I = %3.0fpA', plotSweep, obj.intStimAmp(plotSweep));
                    title(titleString,'FontSize', 25);
                end
                xlabel('Time (s)','FontSize', 25);
                ylabel('Vm (mV)','FontSize', 25);
                set(gca, 'FontSize', 25);
                
                
            end
            
        end
        
        
        function plotLightSpike(obj,varargin)
            
            if isempty(obj.LEDONSpike)
                obj.getLightSpike(10,0.01);
                
            end
            
            
            plotSweep = 0;
            
            for i = 1:length(varargin)
                if strcmp(varargin{1}{i}, 'sweepN')|strcmp(varargin{1}{i}, 'stimN')
                    plotSweep = varargin{1}{i + 1};
                end
                
            end
            
            
            
            
            if obj.ndim == 3
                
                figure;
                if plotSweep == 0
                    for i = 1:obj.nSweeps
                        subplot(3,ceil(obj.nSweeps/3),i);
                        plot(obj.recordTime(i,:),obj.patch(i,:),'.b');
                        
                        hold on
                        obj.plotLED(obj.patch(i,:), i);
                        
                        plot(obj.recordTime(i,1) - 1/obj.SR + obj.baseSpike{1,i}, obj.baseSpike{2,i}, '.c');
                        plot(obj.recordTime(i,1) - 1/obj.SR + obj.AStimSpike{1,i}, obj.AStimSpike{2,i}, '.m');
                        
                        for LEDNum = 1:size(obj.LEDONDataPoint,2)
                            plot(obj.recordTime(i,1) - 1/obj.SR + obj.LEDONSpike{1, i, LEDNum},obj.LEDONSpike{2, i, LEDNum},'.r');
                            plot(obj.recordTime(i,1) - 1/obj.SR + obj.LEDOFFSpike{1, i, LEDNum},obj.LEDOFFSpike{2, i, LEDNum},'.k');
                            
                        end
                        
                        
                        %                 if ~isempty(obj.spikeletData)
                        %                     if ~isempty(obj.spikeletData{1})
                        %                         hold on
                        %                         plot(obj.spikeletData{1},obj.spikeletData{2},'.g');
                        %                     end
                        %                 end
                        title('Episodic LED (membrane voltage) recording');
                        xlabel('Time (s)');
                        ylabel('Vm (mV)');
                        
                        %                         legend('IC recording', 'LED', 'Location','NorthWest');
                    end
                else
                    figure
                    plot(obj.recordTime(plotSweep,:),obj.patch(plotSweep,:),'.b');
                    hold on
                    obj.plotLED(obj.patch(plotSweep,:), plotSweep);
                    
                    plot(obj.recordTime(plotSweep,1) - 1/obj.SR + obj.baseSpike{1,plotSweep}, obj.baseSpike{2,plotSweep}, '.c');
                    plot(obj.recordTime(plotSweep,1) - 1/obj.SR + obj.AStimSpike{1,plotSweep}, obj.AStimSpike{2,plotSweep}, '.m');
                    
                    for LEDNum = 1:size(obj.LEDONDataPoint,2)
                        plot(obj.recordTime(plotSweep,1) - 1/obj.SR + obj.LEDONSpike{1, plotSweep, LEDNum},obj.LEDONSpike{2, plotSweep, LEDNum},'.r');
                        plot(obj.recordTime(plotSweep,1) - 1/obj.SR + obj.LEDOFFSpike{1, plotSweep, LEDNum},obj.LEDOFFSpike{2, plotSweep, LEDNum},'.k');
                        
                    end
                    
                    
                    %                 if ~isempty(obj.spikeletData)
                    %                     if ~isempty(obj.spikeletData{1})
                    %                         hold on
                    %                         plot(obj.spikeletData{1},obj.spikeletData{2},'.g');
                    %                     end
                    %                 end
                    title('Episodic LED (membrane voltage) recording');
                    xlabel('Time (s)');
                    ylabel('Vm (mV)');
                    
                    
                end
                
                
            elseif obj.ndim == 2
                plot(obj.dataTime(:),obj.patch(:),'.b');
                
                hold on
                obj.plotLED(obj.patch);
                
                
                
                plot(obj.LEDONSpike{1},obj.LEDONSpike{2},'.r');
                plot(obj.LEDOFFSpike{1},obj.LEDOFFSpike{2},'.r');
                
                
                
                
                %                 if ~isempty(obj.spikeletData)
                %                     if ~isempty(obj.spikeletData{1})
                %                         hold on
                %                         plot(obj.spikeletData{1},obj.spikeletData{2},'.g');
                %                     end
                %                 end
                title('Gap free IC (membrane voltage) recording');
                xlabel('Time (s)');
                ylabel('Vm (mV)');
                
                legend('IC recording', 'LED', 'Location','NorthWest');
                
            end
        end
        
        function plotSPSSpike(obj, varargin)
            %             if isempty(obj.spikeData)
            %                 obj.getSpike;
            %             end
            
            close all;
            plotDataPt = 1:obj.nSamples;
            
            plotAdjTime = 0;
            % when this variable is one, plot with the time as the
            % reference point (0)
            
            stimN = 0;
            
            ylimVal = [-120, 50];
            
            stimAmp = 0:5:120;
            stimRep = 20;
            
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    
                    if strcmp(varargin{i}, 'plotDataPt')
                        plotDataPt = varargin{i + 1};
                    end
                    
                    
                    if strcmp(varargin{i}, 'stimN')
                        stimN = varargin{i + 1};
                    end
                    
                    if strcmp(varargin{i}, 'ylimVal')
                        ylimVal = varargin{i + 1};
                    end
                    
                    if strcmp(varargin{i}, 'stimAmp')
                        stimAmp = varargin{i + 1};
                    end
                    
                    if strcmp(varargin{i}, 'stimRep')
                        stimRep = varargin{i + 1};
                    end
                    
                end
            end
            
            if stimN ~= 0
                plotSweep = stimN;
                
            else
                plotSweep = 1:obj.nSweeps;
                
            end
            
            plotTime = obj.recordTime(1,plotDataPt)*1000;
            % Time converted to ms
                
            stimAmpI = 0;
            
            for i = 1:length(plotSweep)
                
                stimMod = mod(i, stimRep);
                
                if stimMod == 1
                    stimAmpI = stimAmpI + 1;
                end
                
                if obj.extStimStartPt(plotSweep(i)) == 0
                    continue
                end
                
                plotMod = mod(i, 25);
                if (plotMod == 1)&(length(plotSweep) > 1)
                    figure
                    set(gcf,'units','normalized','outerposition',[0 0 1 1])
                end
                
                if plotMod == 0
                    plotMod = 25;
                end
                if length(plotSweep) > 1
                    subplot(5,5,plotMod);
                end
                
                
                adjPlotTime = plotTime - obj.extStimStartPt(plotSweep(i))/obj.SR*1000;
                
                %                                 plotDataPt = 4900:5500;
                %                                 plot(obj.recordTime(1,:),obj.patch(i,:),'.k');
                
                if plotAdjTime == 1
                    plot(adjPlotTime, obj.patch(plotSweep(i),plotDataPt),'.k');
                else
                    plot(plotTime, obj.patch(plotSweep(i),plotDataPt),'.k');
                end
                
                if ~isempty(obj.baseSpike{1,plotSweep(i)})
                    hold on
                    baseTime = obj.baseSpike{1,plotSweep(i)}*1000;
                    adjBaseTime = baseTime - obj.extStimStartPt(plotSweep(i))/obj.SR*1000;
                    if plotAdjTime == 1
                        plot(adjBaseTime, obj.baseSpike{2,plotSweep(i)},'.c');
                    else
                        plot(baseTime, obj.baseSpike{2,plotSweep(i)},'.c');
                    end
                    %                     else
                    %                         displayString = sprintf('No spikes detected for Sweep #%d', i);
                    %                         display(displayString);
                end
                
                if ~isempty(obj.AStimSpike{1,plotSweep(i)})
                    hold on
                    AStimTime = obj.AStimSpike{1,plotSweep(i)}*1000;
                    adjAStimTime = AStimTime - obj.extStimStartPt(plotSweep(i))/obj.SR*1000;
                    %                     plot(obj.AStimSpike{1,plotSweep(i)},obj.AStimSpike{2,plotSweep(i)},'.m');
                    if plotAdjTime == 1
                        plot(adjAStimTime,obj.AStimSpike{2,plotSweep(i)},'.b','MarkerSize', 30);
                    else
                        plot(AStimTime,obj.AStimSpike{2,plotSweep(i)},'.b','MarkerSize', 30);
                    end
                end
                
                if ~isempty(obj.SPSSpike)
                    if i <= size(obj.SPSSpike,2)
                        if ~isempty(obj.SPSSpike{1,plotSweep(i)})
                            hold on
                            
                            stimTime = obj.SPSSpike{1,plotSweep(i)}*1000;
                            adjStimTime = stimTime - obj.extStimStartPt(plotSweep(i))/obj.SR*1000;
                            if plotAdjTime == 1
                                plot(adjStimTime,obj.SPSSpike{2,plotSweep(i)},'.r','MarkerSize', 20);
                            else
                                plot(stimTime,obj.SPSSpike{2,plotSweep(i)},'.r','MarkerSize', 20);
                            end
                        end
                    end
                end
                
                if ~isempty(obj.stimPosPeak{1,plotSweep(i)})
                    hold on
%                     baseTime = obj.baseSpike{1,plotSweep(i)}*1000;
                    cathStimTime = obj.stimPosPeak{1,plotSweep(i)}*1000;
%                     timeInterval = 1000/obj.SR;
%                     endTime = (obj.extStimEndPt(plotSweep(i)) - obj.extStimStartPt(plotSweep(i)))/obj.SR*1000;
                    adjCathStimTime = cathStimTime  - obj.extStimStartPt(plotSweep(i))/obj.SR*1000;
                    if plotAdjTime == 1
                        plot(adjCathStimTime, obj.stimPosPeak{2,plotSweep(i)},'.g');
                    else
                        plot(cathStimTime, obj.stimPosPeak{2,plotSweep(i)},'.g');
                    end
                end
                
                if ~isempty(obj.stimNegPeak{1,plotSweep(i)})
                    hold on
%                     baseTime = obj.baseSpike{1,plotSweep(i)}*1000;
                    anodStimTime = obj.stimNegPeak{1,plotSweep(i)}*1000;
%                     timeInterval = 1000/obj.SR;
%                     endTime = (obj.extStimEndPt(plotSweep(i)) - obj.extStimStartPt(plotSweep(i)))/obj.SR*1000;
                    adjAnodStimTime = anodStimTime  - obj.extStimStartPt(plotSweep(i))/obj.SR*1000;
                    if plotAdjTime == 1
                        plot(adjAnodStimTime, obj.stimNegPeak{2,plotSweep(i)}, '.','color', [1, 0.5, 0]);
                    else
                        plot(anodStimTime, obj.stimNegPeak{2,plotSweep(i)}, '.','color', [1, 0.5, 0]);
                    end
                end
                
                if plotAdjTime == 1
                    xlim([adjPlotTime(1), adjPlotTime(end)]);
                else
                    xlim([plotTime(1), plotTime(end)]);
                end
               
                ylim([ylimVal(1), ylimVal(2)]);
                
                %                 if ~isempty(obj.SPSSpikelet{1,i})
                %                     hold on
                %                     plot(obj.SPSSpikelet{1,i},obj.SPSSpikelet{2,i},'.g');
                %                 end
                
                %                 if ~isempty(obj.intStimAmp)
                %                     titleString = sprintf('Sweep #%d, current = %dpA', i, obj.intStimAmp(i));
                %                     title(titleString, 'FontSize', 8);
                %                 end
                xlabel('Time (ms)','FontSize', 8);
                ylabel('Vm (mV)','FontSize', 8);
                
               
                
                titleString = sprintf('StimAmp: %d (uA)', stimAmp(stimAmpI)); 
                title(titleString);
                
                
                
            end
        end
            
        function plotHFSSpike(obj, varargin)
%             if isempty(obj.stimSpike)
%                 obj.getHFSSpike;
%             end
            
            plotSweep = 0;
            plotSig = obj.HFSfilteredSig;
            HFSAmp = [];
            
            plotDataPt = 1:obj.nSamples;
            
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    if strcmp(varargin{i}, 'sweepN')
                        plotSweep = varargin{1}{i + 1};
                    end
                    
                    if strcmp(varargin{i}, 'plotType')
                        if strcmp(varargin{i + 1}, 'raw');
                            plotSig = obj.patch;
                        else
                            plotSig = obj.HFSfilteredSig;
                        end
                    end
                    
                    if strcmp(varargin{i}, 'HFSAmp')
                        HFSAmp = varargin{i + 1};
                    end
                    
                    if strcmp(varargin{i}, 'plotDataPt')
                        plotDataPt = varargin{i + 1};
                    end
                    
                end
            end
            figure;

            
            if obj.ndim == 2
                
                
                plot(obj.dataTime,plotSig,'b');
                %                   plot(obj.dataTime,obj.HFSfilteredSig,'.b');
                % plot filtered signal or not
                
                
                hold on
%                 obj.plotExtStim(obj.patch);
                
                
                
                plot(obj.baseSpike{1},obj.baseSpike{2},'.m');
                
                
                plot(obj.AStimSpike{1},obj.AStimSpike{2},'.k');
                
                for HFSStimN = 1:size(obj.HFSSpike,2)
                    plot(obj.HFSSpike{1, HFSStimN} - 1/obj.SR,obj.HFSSpike{2, HFSStimN},'.r');
                    if HFSStimN <= size(obj.HFSGapSpike,2)
                        plot(obj.HFSGapSpike{1, HFSStimN} - 1/obj.SR,obj.HFSGapSpike{2, HFSStimN},'.k');
                    end
                    
                end
                
                %                 if ~isempty(obj.spikeletData)
                %                     if ~isempty(obj.spikeletData{1})
                %                         hold on
                %                         plot(obj.spikeletData{1},obj.spikeletData{2},'.g');
                %                     end
                %                 end
                title('Gap free IC (membrane voltage) recording');
                xlabel('Time (s)');
                ylabel('Vm (mV)');
                
%                 legend('IC recording', 'LED', 'Location','NorthWest');
                
                
            else
                for stimN = 1:obj.nSweeps
                    subplot(5,ceil(obj.nSweeps/5),stimN);
                    hold on
                    plotSigSeg = plotSig(stimN,plotDataPt);
                    plot(obj.recordTime(1,plotDataPt),plotSigSeg,'b');
                    %                   plot(obj.dataTime,obj.HFSfilteredSig,'.b');
                    % plot filtered signal or not

                    %                 obj.plotExtStim(obj.patch);
                    
                    
                    
                    plot(obj.baseSpike{1,stimN},obj.baseSpike{2,stimN},'.m',...
                        'MarkerSize', 20);
                    
                    
                    plot(obj.AStimSpike{1,stimN},obj.AStimSpike{2,stimN},'.k',...
                        'MarkerSize', 20);
                    
                    
                    plot(obj.HFSSpike{1, stimN} - 1/obj.SR,obj.HFSSpike{2, stimN},'.r');

                    
                    %                 if ~isempty(obj.spikeletData)
                    %                     if ~isempty(obj.spikeletData{1})
                    %                         hold on
                    %                         plot(obj.spikeletData{1},obj.spikeletData{2},'.g');
                    %                     end
                    %                 end
                    if ~isempty(HFSAmp)
                        titleString = sprintf('HFS Amplitde: %d\\muA', HFSAmp(stimN));
                    else
                        titleString = 'HFS stimulation';
                    end
                    title(titleString);
                    xlabel('Time (s)');
                    xlim([obj.recordTime(1,plotDataPt(1)), obj.recordTime(1,plotDataPt(end))]);
                    ylabel('Vm (mV)');
                    
                end
                
            end
        end
        
        
        
        
        
        
        
        
        
        
        
        function plotLED(obj,data , i)
            if (obj.ndim == 2)
                plot(obj.dataTime, obj.LED*2*max(diff(data)) + 0.95*median(data),'color',[1,0.5,0]);
            else
                plot(obj.recordTime(i,:), obj.LED(i,:)*2*max(diff(data)) + 0.95*median(data),'color',[1,0.5,0]);
            end
            
        end
        
        
        function plotExtStim(obj,data,i)
            if (obj.ndim == 2)
                plot(obj.dataTime, obj.extStim*2*max(diff(data)) + 0.95*median(data),'color',[1,0.5,0]);
            else
                plot(obj.recordTime(i,:), obj.extStim(i,:)*2*max(diff(data)) + 0.95*median(data),'color',[1,0.5,0]);
            end
            
        end
        
        function calcDF(obj)
            % nPrevBaseFrame = number of frames before stimulation used for
            % baseline fluorescence
            imSignal = obj.imSig;
            data = imSignal.smData(7,:);
            for i = 1:obj.nSweeps
                obj.F0(i) = mean(data(obj.baselineFrame(i,:)));
                obj.DF(i,:) = data(obj.sweepFrame(i,:)) - obj.F0(i);
                obj.DFF(i,:) = obj.DF(i,:)/obj.F0(i)*100;
                % Fluorescent change represented into percentage
                
                
            end
        end
        
        function assignFPS(obj, adFrame)
            % adFrames is the additional frames added before and after the
            % recorded sweeps frames
            
            imSignal = obj.imSig;
            
            if ~isempty(imSignal.framesPerSweep)
                obj.framesPerSweep = imSignal.framesPerSweep;
            else
                obj.framesPerSweep = obj.recordTime(1,end)*imSignal.SR;
            end
            
            intSweepFrame = ceil(obj.framesPerSweep);
            % carry the number of frames in each sweep into the closest integer
            
            intSTSTFrame = ceil(obj.STST*imSignal.SR);
            
            if nargin < 2
                obj.sweepFrame = zeros(obj.nSweeps, (intSweepFrame + 1));
            else
                obj.sweepFrame = zeros(obj.nSweeps, (intSweepFrame + 1 + 2*adFrame));
            end
            
            for i = 1:obj.nSweeps
                
                
                % the frame number in each sweep (changed over
                % different i)
                
                
                
                
                if i == 1
                    obj.sweepFrame(1,:) = 1:(intSweepFrame + 2*adFrame + 1);
                else
                    if nargin < 2
                        obj.sweepFrame(i,:) = intSTSTFrame*(i - 1):intSTSTFrame*(i - 1) + intSweepFrame;
                    else
                        obj.sweepFrame(i,:) = (intSTSTFrame*(i - 1) - adFrame):(intSTSTFrame*(i - 1) + intSweepFrame + adFrame);
                    end
                end
                obj.sweepFrameTime(i,:) = obj.sweepFrame(i,:)/imSignal.SR;
                % convert the frame number in each sweep into the
                % corresponding time
                
            end
        end
        
        function phasePlot(obj,nSweep)
            % plotting the membrane voltage phase plot
            plot(obj.patch(nSweep,1:end - 1)), diff(obj.patch(nSweep,1:end - 1));
            
            if ~isempty(obj.stimInitAmp)
                titleString = sprintf('Sweep #%d, current = %dpA', nSweep, obj.intStimAmp(nSweep));
                title(titleString);
            end
            xlabel('Vm (mV)');
            ylabel('diff(Vm) (mV)');
        end
        
        function plotOverInt(obj,sub)
            % subplot nSweeps of results into n number of subplots
            
            imSignal = obj.imSig;
            
            figure(2)
            
            if (nargin > 1)&(sub == 1)
                for i = 1:obj.nSweeps
                    data = obj.DFF(i,:);
                    %                     data = obj.sweepFrame(i,:);
                    
                    %                     frameTime = obj.sweepFrameTime(i,:);
                    frameTime = obj.sweepFrameTime(i,:);
                    
                    subplot(3,ceil(obj.nSweeps/3),i);
                    %                     imObj.framesPerSweep = obj.recordTime(1,end)*imObj.SR;
                    
                    % frames per sweep
                    
                    % obj.recordTime(1,end) represents the sweeping time
                    %                     sweepEntries = (floor((i - 1)*imObj.framesPerSweep + 1):ceil(i*imObj.framesPerSweep));
                    
                    %                     plot(SweepFrameTimes,imSignal.convInt(7,sweepFrame),'b');
                    
                    
                    plot(frameTime,data,'-g','LineWidth', 3);
                    hold on
                    
                    stimTimePlot = (i - 1)*obj.STST + obj.stimTime;
                    stimTimeAmp = ...
                        min(data) - 0.1*(max(data) - min(data))*ones(1,2);
                    
                    xlabel('Time (s)');
                    %                     ylabel('Fluorescence (AU)');
                    ylabel('DF/F0 (%)');
                    
                    currentStimAmp = obj.stimInitAmp + obj.stimDAmp*(i - 1);
                    if strcmp(obj.stimType, 'IC')
                        titleString = sprintf('Sweep #%d, Current = %dpA', i, currentStimAmp);
                    else
                        titleString = sprintf('Sweep #%d, Voltage = %dmV', i, currentStimAmp);
                    end
                    title(titleString);
                    
                    plot(stimTimePlot, stimTimeAmp,'-k','LineWidth', 3);
                    
                    
                    
                    %                     plot(obj.spikeData{1,i},data((i*ceil(obj.spikeData{1,1}*imSignal.SR)),'+r');
                    spikeDataPointInSweep  = ceil(obj.spikeData{1,i}*imSignal.SR) - obj.sweepFrame(i,1) + 1;
                    
                    for j = 1:length(obj.spikeData{1,i})
                        if obj.spikeData{1,i}(j) < stimTimePlot(2)
                            % if the spikes occur within the stimulation
                            % time, it is the usual depolarisation spikes.
                            % Otherwise, if the spikes occurs after
                            % stimulation, the spike is
                            % after-hyperpolarisation spike, which is a
                            % feature of the OFF RGCs.
                            plot(obj.spikeData{1,i}(j),data(spikeDataPointInSweep(j)),'+r');
                        else
                            plot(obj.spikeData{1,i}(j),data(spikeDataPointInSweep(j)),'+b');
                        end
                    end
                    
                    
                    
                    hold off
                end
            end
        end
        
        
    end
    
    
    
    
end








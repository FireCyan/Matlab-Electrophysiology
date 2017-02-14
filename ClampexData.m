% Copyright 2017 Chih-Yu (John) Yang
% This class is for loading my Digidata1440A recordings from retinal
% ganglion cell patching. This class may be helpful for other
% electrophysiologists who want to manipulate their data in Matlab rather
% than Axon softwares. I have also uploaded a script showing how to construct
% a Clampex object and do anaylsis. 

% Because I do artificial stimulation on retinal ganglion cells, I have
% channels that record external stimulating electrode trigger and LED
% triggers, which may not be relevant for conventional
% electrophysiologists.

% Need to download abdload.mat, which is on Matlab File Exchange. I have
% also included in my repository, but just saying this for copyright
% purpose, abfload.mat is written by Harald Hentschke, not me.

classdef ClampexData < hgsetget
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    % handling of Clampex data
    
    properties
        
        %% Basic properties
        path = '';
        wholeData = [];
        cDate = '';
        % nSamples x nChannels x nSweeps
        
        gapFreeData = [];
        SR = ''; % sampling rate
        rsSR = ''; % resampled sampling rate
        headerStruct = [];
        nChannels = '';
        nSamples = '';
        nSweeps = ''; % number of sweeps in the recordings
        ndim = ''; % ndimensions. ndim = 2 means gap-free recording. ndim = 3 means sweep recording
        nLine = '';
        
        sweepFrames = [];
        % number of frames in a sweep. This is for CaImage frame number
        % when using episodic recording mode
        

        LEDSeqSETime = [];
        % The time when LED ON and OFF sequence starts and end. 1st rows
        % shows the sequence starts time and 2nd row shows the sequence end
        % times.
        
        %% Patched neuron
        patch = [];
        patchNo = '';
        patchGapFree = [];
        rsPatch = '';
        
        %% Frame trigger
        frame = [];
        frameNo = '';
        frameGapFree = [];;
        rsFrame = '';
        STSFrameATime = ''; % one frame acquisition time (average), should be fairly consistent though
        frameTime = [];
        frameProcessTime = [];
        % The time it takes to for processing and sending the data to the computer
        
        % the frame acquiring time, can be seen from the print screen shot,
        % but it is best to deduce the accurate frame time by the high
        % sampling rate recording. It will be 2 x nSamples array where the
        % 1st row contains the time the particular frame starts imaging and
        % 2nd row contains the time when that frame imaging ends.
        
        %% Line trigger
        line = [];
        lineNo = '';
        lineGapFree = [];
        rsLine = '';
        lineTime = [];
        % the starting time of each line scan. It is organised into to
        % nframes x (nlines per frame)
        lineGap = '';
        % the gap between each line acquiring time
        %         lineAStartTime = [];
        lineATime = '';
        %         linePerFrameTime = '';
        linesPerFrame = '';
        lineFrameTime = '';
        % Using the line triggers to deduce the frame time
        
        %% Extracellular stimulation recording
        extStim = [];
        extStimNo = '';
        
        %% Intracellular stimulation recording
        intStim = [];
        intStimNo = '';
        intStimGapFree = [];
        
        
        %% LED stimulation recording
        LED = [];
        LEDNo = '';
        LEDGapFree = [];;
        rsLED = '';
%         LEDOnSETime = [];
%         % LED ON start and end time
%         LEDOffSETime = [];
%         % LED OFF start and end time
%         LEDOnTime = [];
%         
%         LEDOffTime = [];
        
        
        stimFlag = 0;
        % when stimFlag = 0, it means the imaging is a spontaneous
        % recording. If stimFlag = 1, it means the imaging is accompanied
        % with stimulation.
        
        
        %% Time points
        dataTime = [];
        % total dataTime including the time when digitiser is not recording
        % between episodic gaps
        rsDataTime = [];
        % resampled dataTime
        recordTime = [];
        % the time when the digitsier is recording. It is organised into
        % nSweeps x nSamples
        nRS = '';
        STST = ''; % sweep start-to-start time
        
        gapFreeTime = [];
        % The gap free times of data point recorded in sweep mode 
        % and converted to gap free mode.
        
        
        %         syncImIndex = [];
        
        testTime = 2;
        nFrameTest = 5;
        nLineTest = 100;
        nLEDTest = 1;
        
        STSTrigger = [];
        
        STSTriggerNo = '';
        
                
        STSStartPt = [];
        STSEndPt = [];
        STSTotalPt;
        
        extStimStartPt = [];
        extStimEndPt = [];
        
    end
    
    
    
    methods (Static)
        
        function filePath = findFile(folderPath,num)
            % This static function is mainly for automatically anaylse the
            % pClamp data. It won't output the object instance in the workspace
            % but will save the analysed results and objects into the folder.
            % If the user want to look at the instance and its properties,
            % he/she should use the initialisation function
            
            % STST is start-to-start time. This value needs to be input in
            % by the user if pClamp recording is in episodic mode.
            % Otherwise accurate timing of each timing cannot be calculated
            if num > 99
                findFile = sprintf('\\*%d.abf', num);
            elseif num > 9
                findFile = sprintf('\\*0%d.abf', num);
            else
                findFile = sprintf('\\*00%d.abf', num);
            end
            
               
            fileStruct = dir([folderPath,findFile]);
            filePath = fullfile(folderPath,fileStruct.name);
        end
        
        function obj = loadobj(s)
            if isstruct(s)
                obj = ClampexData;
                obj = reload(obj,s);
            end
        end
    end
    
    
    

    
    methods
        
        
        function obj = ClampexData(path,assignNo,abfEnd, varargin)
            
            deleteSweep = [];
            cDate = '';
            
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    if strcmp(varargin{i}, 'deleteSweep')
                        deleteSweep = varargin{i + 1};
                    end
                    if strcmp(varargin{i}, 'cDate')
                        cDate = varargin{i + 1};
                    end
                    
                end
            end
            
            
            obj.path = path;
            [obj.wholeData,ST,obj.headerStruct] = abfload(obj.path,'stop',abfEnd);
            obj.SR = 1/(ST*10^-6); % convert to s
            obj.rsSR = obj.SR;
            
            obj.cDate = cDate;
            
            obj.ndim = ndims(obj.wholeData);
            
            if ~isempty(deleteSweep)
                obj.assignChannel(assignNo, 'deleteSweep', deleteSweep);
            else
                obj.assignChannel(assignNo)
            end
        end
        
        function s = saveobj(obj)
            s = obj;
        end
        
        function obj = reload(obj, s)
        end
        
        function assignChannel(obj,assignNo,varargin)
            
            deleteSweep = [];
            
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    if strcmp(varargin{i}, 'deleteSweep')
                        deleteSweep = varargin{i + 1};
                    end
                end
            end
            
            %% handling sweep recordings
            
            if obj.ndim == 3
                
                if ~isempty(assignNo)
                    for i = 1:length(assignNo)
%                         recChNames = {};
%                         recChNames{i} = obj.headerStruct.recChNames{assignNo(i)};
                        for j = 1:size(obj.wholeData,3)
                            tempData(:,i,j) = obj.wholeData(:,assignNo(i),j)';
                        end
                        %                         tempRecChNames{assignNo(i)} = recChNames{i};
                    end
                    obj.wholeData = tempData;
                    
                    %                     recChNames = tempRecChNames;
                else
                    recChNames = obj.headerStruct.recChNames;
                end
                               
                obj.wholeData(:,:,deleteSweep) = [];
                
                if length(obj.headerStruct.sweepStartInPts) > 1
                    obj.STST = mean(diff(obj.headerStruct.sweepStartInPts))/obj.SR;
                else 
                    obj.STST = 0;
                end
                [obj.nSamples, obj.nChannels, obj.nSweeps] = size(obj.wholeData);
                
                obj.dataTime = 1/obj.SR:1/obj.SR:(obj.STST*obj.nSweeps);
                obj.rsDataTime = 1/obj.SR:1/obj.SR:(obj.STST*obj.nSweeps);
                
                for i = 1:obj.nSweeps
                    obj.recordTime(i,:) = ...
                        ((i - 1)*obj.STST + 1/obj.SR):1/obj.SR:((i - 1)*obj.STST + obj.nSamples/obj.SR);
                end
                
                
                
                for i = 1:size(obj.wholeData,2)
                    if strcmp(recChNames{i},'V0')|strcmp(recChNames{i},'IN 0')|strcmp(recChNames{i},'I0')
                        % Patch recording
                        obj.patch = zeros(obj.nSweeps,obj.nSamples);
                        obj.patch = reshape(obj.wholeData(:,i,:),obj.nSamples,obj.nSweeps)';
                        obj.patchNo = i;
                        obj.patchGapFree = zeros(2,round(obj.nSweeps*obj.nSamples));
                        
                        if obj.STST ~= 0
                            for j = 1:obj.nSweeps
                                obj.patchGapFree(1,round((j - 1)*obj.nSamples + 1):round(j*obj.nSamples)) =...
                                    ((j - 1)*obj.STST + 1/obj.SR):1/obj.SR:((j - 1)*obj.STST + obj.nSamples/obj.SR);
                                
                                obj.patchGapFree(2,round((j - 1)*obj.nSamples + 1):round(j*obj.nSamples)) =...
                                    obj.patch(j,:);
                                % data amplitude
                            end
                        end
                        
                    elseif strcmp(recChNames{i},'Ext_stim')|strcmp(recChNames{i},'IN 4')
                        % External stimulating electrode recording
                        obj.extStim = zeros(obj.nSweeps,obj.nSamples);
                        obj.extStim = reshape(obj.wholeData(:,i,:),obj.nSamples,obj.nSweeps)';
                        obj.extStimNo = i;
                        obj.stimFlag = 1;
                        
                        for stimN = 1:size(obj.extStim,1)
%                             [~, extStimStartLocs] =...
%                                 findpeaks(obj.extStim(stimN,:), 'MINPEAKHEIGHT', 2,'MINPEAKDISTANCE',round(0.0005*obj.SR));
                            
                            extStimStartLocs = find(diff(obj.extStim(stimN,:)) > 2);
                            
%                             [~, extStimEndLocs] =...
%                                 findpeaks(-diff(obj.extStim(stimN,:)), 'MINPEAKHEIGHT', 2,'MINPEAKDISTANCE',round(0.0005*obj.SR));
                            
                            extStimEndLocs = find(diff(-obj.extStim(stimN,:)) > 2);
                            
                            if ~isempty(extStimStartLocs)
                                if extStimStartLocs(1) < 100
                                    % if extStimLocs(1) is less than 100 data
                                    % points, then the trigger is for initiating
                                    % digitiser
                                    obj.extStimStartPt(stimN) = extStimStartLocs(2);
                                    
                                else
                                    obj.extStimStartPt(stimN) = extStimStartLocs(1);
                                end
                                
                                obj.extStimEndPt(stimN) = extStimEndLocs(end);
                            else
                                obj.extStimStartPt(stimN) = 0;
                                obj.extStimEndPt(stimN) = 0;
                            end
                                
                            
                            
                            
                        end
                        

                        % 1st one is for starting Digitiser
                        
                    elseif strcmp(recChNames{i},'Line_trig')|strcmp(recChNames{i},'IN 6')|strcmp(recChNames{i},'LVAL')
                        % Trigger out for line scanning from two-photon laser microscope
                        obj.line = zeros(obj.nSweeps,obj.nSamples);
                        obj.line = reshape(obj.wholeData(:,i,:),obj.nSamples,obj.nSweeps)';
                        obj.lineNo = i;
                    elseif strcmp(recChNames{i},'Frame_trig')|strcmp(recChNames{i},'IN 7')|strcmp(recChNames{i},'FVAL')
                        % Trigger out for frame from two-photon laser microscope
                        obj.frame = zeros(obj.nSweeps,obj.nSamples);
                        obj.frame = reshape(obj.wholeData(:,i,:),obj.nSamples,obj.nSweeps)';
                        obj.frameNo = i;
                        
                        % Thorlab camera triggers (trigger up = start of
                        % acqusition, trigger down = end of acquisition)
                        for frameI = 1:size(obj.frame,1)
                            frameAStartTime = obj.dataTime(diff(obj.frame(frameI,:)) > 1.5) + 1/obj.SR;
                            % frame acquisition starting time
                            if isempty(frameAStartTime)
                                continue
                            end
                            tooClose = diff(frameAStartTime) < 1.5/obj.SR;
                            frameAStartTime(tooClose) = [];
                            frameAEndTime = obj.dataTime(diff(obj.frame(frameI,:)) < -1);
                            tooClose = diff(frameAEndTime) < 1.5/obj.SR;
                            % if the two frame found is too close, ie, less
                            % than 1 sampling point, then chances is the
                            % same frame time is being counted twice
                            frameAEndTime(tooClose) = [];
                            %                             frameAEndTime(end) = [];
                            %                             obj.meanFrameATime = mean(frameAEndTime - frameAStartTime);
                            
%                             frameAStartTime(end) = [];
%                             frameAEndTime(end) = [];
                            
                            shiftFrameATime = frameAStartTime;
                            shiftFrameATime(1) = [];
                            %                             shiftFrameATime(end) = [];
                            
                            obj.STSFrameATime{frameI} = median(shiftFrameATime - frameAStartTime(1:(end - 1)));
                            obj.frameTime{frameI} = frameAEndTime;
                            
                            nAStartTime = length(frameAStartTime);
                            nAEndtime = length(frameAEndTime);
                            obj.frameProcessTime{frameI} = mean(frameAEndTime - frameAStartTime);
                            obj.sweepFrames(frameI) = length(frameAEndTime);
                        end
                        
                        
                        
                    elseif strcmp(recChNames{i},'LED')|strcmp(recChNames{i},'IN 5')
                        % LED (on = high, off = low) recording
                        obj.LED = zeros(obj.nSweeps,obj.nSamples);
                        obj.LED = reshape(obj.wholeData(:,i,:),obj.nSamples,obj.nSweeps)';
                        obj.LEDNo = i;
                        obj.stimFlag = 1;
                    elseif strcmp(recChNames{i},'I1')|strcmp(recChNames{i},'IN 1')
                        obj.intStim = zeros(obj.nSweeps,obj.nSamples);
                        obj.intStim = reshape(obj.wholeData(:,i,:),obj.nSamples,obj.nSweeps)';
                        obj.intStimNo = i;
                        obj.stimFlag = 1;
                        
                    end
                end
                
                
                obj.convertToGapFree;
                
                
                if ~isempty(obj.LEDGapFree)
                    LEDOnStartTime = obj.dataTime(diff(obj.gapFreeTime(:)) > 2.5) + 1/obj.SR;
                    LEDOnEndTime = obj.dataTime(diff(obj.gapFreeTime(:)) < -2.5);
                    
                    LEDOffStartTime = LEDOnEndTime + 1/obj.SR;
                    LEDOffEndTime = LEDOnStartTime - 1/obj.SR;
                end

                
                
                
                %% handling normal gap free recording
            elseif obj.ndim == 2
                [obj.nSamples, obj.nChannels] = size(obj.wholeData);
                obj.nSweeps = 1;
                obj.dataTime = 1/obj.SR:1/obj.SR:(obj.nSamples*1/obj.SR);
                obj.rsDataTime = 1/obj.SR:1/obj.SR:(obj.nSamples*1/obj.SR);
                
                recChNames = obj.headerStruct.recChNames;
                
                obj.recordTime = ...
                        (1/obj.SR):(1/obj.SR):(obj.nSamples/obj.SR);
                
                if ~isempty(assignNo)
                    for i = 1:length(assignNo)
%                         recChNames{i} = obj.headerStruct.recChNames{assignNo(i)};
                        tempData(:,i) = obj.wholeData(:,assignNo(i))';
%                         tempRecChNames{assignNo(i)} = recChNames{i};
                    end
                    obj.wholeData = tempData;
%                     recChNames = tempRecChNames;
%                 else
%                     recChNames = obj.headerStruct.recChNames;
                end
                
                for i = 1:length(recChNames)
                    if strcmp(recChNames{i},'V0')|strcmp(recChNames{i},'IN 0')|strcmp(recChNames{i},'I0')
                        obj.patch = zeros(1,obj.nSamples);
                        obj.patch = obj.wholeData(:,i)';
                        obj.patchNo = i;
                        
                    elseif strcmp(recChNames{i},'I1')|strcmp(recChNames{i},'IN 1')
                        obj.intStim = zeros(obj.nSweeps,obj.nSamples);
                        obj.intStim = reshape(obj.wholeData(:,i,:),obj.nSamples,obj.nSweeps)';
                        obj.intStimNo = i;
                        obj.stimFlag = 1;
                        
                        
                    elseif strcmp(recChNames{i},'Ext_stim')|strcmp(recChNames{i},'IN 4')
                        % There are two modes since 20150707 experiment,
                        % one is recording trigger (as from 20150707), the
                        % other is recording the voltage waveform (before
                        % 20150707.
                        obj.extStim = zeros(1,obj.nSamples);
                        obj.extStim = obj.wholeData(:,i)';
                        obj.extStimNo = i;
                        obj.stimFlag = 1;
                        [extStimPeak, extStimStartLocs] =...
                            findpeaks(obj.extStim, 'MINPEAKHEIGHT', 2,'MINPEAKDISTANCE',round(0.0005*obj.SR));
                        
                        if extStimStartLocs(1) < 100
                            % if extStimLocs(1) is less than 100 data
                            % points, then the trigger is for initiating
                            % digitiser
                            obj.extStimStartPt = extStimStartLocs(2:end);
                        else
                            obj.extStimStartPt = extStimStartLocs;
                        end
                        % 1st one is for starting Digitiser
                        
                        
                        
                    elseif strcmp(recChNames{i},'LED')|strcmp(recChNames{i},'IN 5')
                        obj.LED = zeros(1,obj.nSamples);
                        obj.LED = obj.wholeData(:,i)';
                        obj.LEDNo = i;
                        obj.stimFlag = 1;
                    elseif strcmp(recChNames{i},'ExtSync2')
                        obj.STSTrigger = zeros(1,obj.nSamples);
                        obj.STSTrigger = obj.wholeData(:,i)';
                        obj.STSTriggerNo = i;
                        [STSPeak, STSLocs] =...
                            findpeaks(obj.STSTrigger, 'MINPEAKHEIGHT', 2,'MINPEAKDISTANCE',round(0.5*obj.SR));
                        if ~isempty(STSPeak)
                            if STSLocs(1) < 100
                                % If less than 100 data points, then this
                                % trigger is for initialising something, not
                                % for time stamping stimulation or start of
                                % stimulation.
                                STSLocs(1) = [];
                            end
                            
                            obj.STSStartPt = STSLocs;
                            meanSTSDiff = round(mean(diff(STSLocs)));
                            obj.STSEndPt = [STSLocs(2:end), STSLocs(end) + meanSTSDiff];
                            STSTemp = obj.STSEndPt - obj.STSStartPt;
                            obj.STSTotalPt = min(STSTemp);
                        end
                        

                        
                       
                    elseif strcmp(recChNames{i},'Line_trig')|strcmp(recChNames{i},'IN 6')|strcmp(recChNames{i},'LVAL')
                        obj.line = zeros(1,obj.nSamples);
                        obj.line = obj.wholeData(:,i)';
                        obj.lineNo = i;
                        
                        %                             lineAStartTime = 1/obj.SR;
                        
                        lineAStartTime = obj.dataTime(diff(obj.line) > 2) + 1/obj.SR;
                        % frame acquisition starting time
                        tooClose = diff(lineAStartTime) < 1.5/obj.SR;
                        lineAStartTime(tooClose) = [];
                        lineAEndTime = obj.dataTime(diff(obj.line) < -2);
                        tooClose = diff(lineAEndTime) < 1.5/obj.SR;
                        % if the two frame found is too close, ie, less
                        % than 1 sampling point, then chances is the
                        % same frame time is being counted twice
                        
                        lineAEndTime(tooClose) = [];
                        %                             frameAEndTime(end) = [];
                        lineLengthAvg = mean(diff(lineAStartTime));
                        
                        tooFar = find(diff(lineAStartTime) > 1.5*lineLengthAvg);
                        shiftLineAST = circshift(lineAStartTime, [0 -1]);
                        shiftLineAST(end) = [];
                        
                        lineAStartTime(end) = [];
                        if ~isempty(tooFar)
                            obj.linesPerFrame = tooFar(1);
                            obj.lineATime = median(shiftLineAST(1:(tooFar(1) - 1)) - lineAStartTime(1:(tooFar(1) - 1)));
                            %                     obj.linePerFrameTime = tooFar(1)*obj.lineATime;
                            obj.lineGap = lineAStartTime(tooFar(1) + 1) - lineAStartTime(tooFar(1)) - obj.lineATime;
                            %                             obj.meanLineATime = mean(lineAEndTime - lineAStartTime);
                            obj.lineTime = lineAStartTime;
                        end
                        
                        
                    elseif strcmp(recChNames{i},'Frame_trig')|strcmp(recChNames{i},'IN 7')|strcmp(recChNames{i},'FVAL')
                        
                        obj.frame = zeros(1,obj.nSamples);
                        obj.frame = obj.wholeData(:,i)';
                        obj.frameNo = i;
                        
                        frameAStartTime = [];
                        % The first frame acquisition starts from the
                        % very beginning
                        if (strcmp(obj.cDate, '20150707'))
                            % For scientific camera frame, where trigger up
                            % means start of acquisition and trigger down
                            % means end of acquisition
                            frameAStartTime = [frameAStartTime, obj.dataTime(diff(obj.frame) > 1.5) + 1/obj.SR];
                            % frame acquisition starting time
                            tooClose = diff(frameAStartTime) < 1.5/obj.SR;
                            frameAStartTime(tooClose) = [];
                            frameAEndTime = obj.dataTime(diff(obj.frame) < -1);
                            tooClose = diff(frameAEndTime) < 1.5/obj.SR;
                            % if the two frame found is too close, ie, less
                            % than 1 sampling point, then chances is the
                            % same frame time is being counted twice
                            frameAEndTime(tooClose) = [];
                            %                             frameAEndTime(end) = [];
                            %                             obj.meanFrameATime = mean(frameAEndTime - frameAStartTime);
                            
                            frameAStartTime(end) = [];
                            frameAEndTime(end) = [];
                            
                            shiftFrameATime = frameAStartTime;
                            shiftFrameATime(1) = [];
%                             shiftFrameATime(end) = [];
                            
                            obj.STSFrameATime = median(shiftFrameATime - frameAStartTime(1:(end - 1)));
                            obj.frameTime = frameAEndTime;
                            
                            nAStartTime = length(frameAStartTime);
                            nAEndtime = length(frameAEndTime);
                            obj.frameProcessTime = mean(frameAEndTime - frameAStartTime); 
                            
                        else
                            % For Lavision 2P frame, where trigger down
                            % means start of acquisition and trigger up
                            % means end of acquisition
                            frameAStartTime = [frameAStartTime, obj.dataTime(diff(obj.frame) < -1) + 1/obj.SR];
                            % frame acquisition starting time
                            tooClose = diff(frameAStartTime) < 1.5/obj.SR;
                            frameAStartTime(tooClose) = [];
                            frameAEndTime = obj.dataTime(diff(obj.frame) > 1.5);
                            tooClose = diff(frameAEndTime) < 1.5/obj.SR;
                            % if the two frame found is too close, ie, less
                            % than 1 sampling point, then chances is the
                            % same frame time is being counted twice
                            frameAEndTime(tooClose) = [];
                            
                            if ~isempty(obj.lineTime)
                                % 2P starts scanning without the first
                                % frame trigger down, so assign first line
                                % time as frame time
                                frameAStartTime = [obj.lineTime(1), frameAStartTime];
                            end
                            
%                             if 
                            %                             frameAEndTime(end) = [];
                            %                             obj.meanFrameATime = mean(frameAEndTime - frameAStartTime);
                            shiftFrameATime = circshift(frameAStartTime, [0, -1]);
                            shiftFrameATime(end) = [];
                            
                            obj.STSFrameATime = median(shiftFrameATime - frameAStartTime(1:(end - 1)));
                            obj.frameTime = frameAEndTime;
                            
                            nAStartTime = length(frameAStartTime);
                            nAEndtime = length(frameAEndTime);
                            obj.frameProcessTime = mean(frameAStartTime - frameAEndTime(1:nAStartTime)); 
                            
                            
                            
                        end
                        
                    end
                end
            end
        end
        
        
        function displayData(obj) % to deduce what type of data is recorded in each channel by inspection
            
            for i = 1:obj.nChannels
                subplot(obj.nChannels, 1,i);
                if obj.ndim == 3
                    if i == obj.LEDNo
                        plot(obj.recordTime(1,:),obj.LED(1,:));
                        ylabel('LED ON/OFF (0/5V)');
                        xlabel('Time (s)');
                    elseif i == obj.patchNo
                        plot(obj.recordTime(1,:),obj.patch(1,:));
                        ylabel('Vm (mV)');
                        xlabel('Time (s)');
                    elseif i == obj.frameNo
                        plot(obj.recordTime(1,:),obj.frame(1,:));
                        ylabel('Frame A/R (0/5V)');
                        xlabel('Time (s)');
                    elseif i ==  obj.lineNo
                        plot(obj.recordTime(1,:),obj.line(1,:));
                        ylabel('Line A/R (5/0V)');
                        xlabel('Time (s)');
                    else
                        xlabel('Corrupted signal');
                    end
                elseif obj.ndim == 2
                    if i == obj.LEDNo
                        plot(obj.dataTime, obj.LED);
                        ylabel('LED ON/OFF (0/5V)');
                        xlabel('Time (s)');
                    elseif i == obj.patchNo
                        plot(obj.dataTime, obj.patch);
                        ylabel('Vm (mV)');
                        xlabel('Time (s)');
                    elseif i == obj.frameNo
                        plot(obj.dataTime, obj.frame);
                        ylabel('Frame A/R (0/5V)');
                        xlabel('Time (s)');
                    elseif i == obj.lineNo
                        plot(obj.dataTime, obj.line);
                        ylabel('Line A/R (5/0V)');
                        xlabel('Time (s)');
                    else
                        xlabel('Corrupted signal');
                    end
                end
            end
        end
        
        
        function convertToGapFree(obj)
            if ~isempty(obj.LED)
                obj.LEDGapFree = zeros(1,round(obj.nSweeps*obj.nSamples));
            end
            
            if ~isempty(obj.frame)
                obj.frameGapFree = zeros(1,round(obj.nSweeps*obj.nSamples));
                %                 obj.frameGapFree(1:round(1/obj.fps*obj.SR):end) = 5;
            end
            
            if ~isempty(obj.line)
                obj.lineGapFree = zeros(1,round(obj.nSweeps*obj.nSamples));
            end
            
            if ~isempty(obj.intStim)
            obj.intStimGapFree = zeros(1,round(obj.nSweeps*obj.nSamples));
            end
            




            

            for i = 1:obj.nSweeps
                if isempty(obj.gapFreeTime) 
                    obj.gapFreeTime = zeros(1,round(obj.nSweeps*obj.nSamples));                    
                end
                
                obj.gapFreeTime(round((i - 1)*obj.nSamples + 1):round(i*obj.nSamples)) =...
                    ((i - 1)*obj.STST + 1/obj.SR):1/obj.SR:((i - 1)*obj.STST + obj.nSamples/obj.SR);
                if ~isempty(obj.LED)
                    obj.LEDGapFree(round((i - 1)*obj.nSamples + 1):round(i*obj.nSamples)) =...
                        obj.LED(i,:);
                    
                    if isempty(obj.gapFreeTime)
                        
                        
                    end
                    
%                     LEDOnStartTime = obj.dataTime(diff(obj.LED(:)) > 2.5) + 1/obj.SR;
%                     LEDOnEndTime = obj.dataTime(diff(obj.LED(:)) < -2.5);
%                     
%                     LEDOffStartTime = LEDOnEndTime + 1/obj.SR;
%                     LEDOffEndTime = LEDOnStartTime - 1/obj.SR;
%                     
%                     
%                     for j = 1:length(LEDOnStartTime)
%                         obj.LEDOnSETime(1,j) = LEDOnStartTime(j);
%                         obj.LEDOnSETime(2,j) = LEDOnEndTime(j);
%                         obj.LEDOnTime = [obj.LEDOnTime, LEDOnStartTime(j):1/obj.SR:LEDOnEndTime(j)];
%                         if j ~= length(LEDOnStartTime)&((LEDOffEndTime(j + 1) - LEDOffStartTime(j)) < 1)
%                             obj.LEDOffSETime(1,j) = LEDOffStartTime(j);
%                             obj.LEDOffSETime(2,j) = LEDOffEndTime(j);
%                             obj.LEDOffTime = [obj.LEDOffTime, LEDOffStartTime(j):1/obj.SR:LEDOffEndTime(j + 1)];
%                         end
%                     end
%                     
%                     obj.LEDOnTime(2,:) = 5;
%                     obj.LEDOffTime(2,:) = 0;
                    
                end
                
                if ~isempty(obj.intStim)
                    obj.intStimGapFree(round((i - 1)*obj.nSamples + 1):round(i*obj.nSamples)) =...
                        obj.intStim(i,:);
                end
                
                if ~isempty(obj.frame)
                    obj.frameGapFree(round((i - 1)*obj.nSamples + 1):round(i*obj.nSamples)) =...
                        obj.frame(i,:);
                end
                
                if ~isempty(obj.line)
                    %                    obj.lineGapFree(round(((i - 1)*sweepT*obj.SR) + 1):round(((i - 1)*sweepT*obj.SR) + obj.nSamples)) =...
                    %                         obj.line(:,i);
                    obj.lineGapFree(round((i - 1)*obj.nSamples + 1):round(i*obj.nSamples)) =...
                        obj.line(i,:);
                end
                
            end
        end
        
        
        function rsElec(obj,newSR)
            if ~isempty(obj.LED)
                obj.rsLED = resample(obj.LED, newSR, round(obj.SR));
            end
            
            if ~isempty(obj.patch)
                obj.rsPatch = resample(obj.patch, newSR, round(obj.SR));
            end
            if ~isempty(obj.frame)
                obj.rsFrame = resample(obj.frame, newSR, round(obj.SR));
            end
            if ~isempty(obj.line)
                obj.rsLine = resample(obj.line, newSR, round(obj.SR));
            end
            
            obj.nRS = ceil(size(obj.wholeData,1)*newSR/round(obj.SR));
            
            obj.rsSR = newSR;
            obj.rsDataTime = 1/obj.rsSR:1/obj.rsSR:(obj.nRS*1/obj.rsSR);
        end
        
        
    end
    
    
end


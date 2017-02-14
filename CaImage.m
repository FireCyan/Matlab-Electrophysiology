% Copyright 2017 Chih-Yu (John) Yang
% CaImage is for analysing images from Ca2+ imaging. I would recommend draw
% region-of-interest (ROI) by ImageJ, then load the ROIs here to do the
% analysis.


classdef CaImage < hgsetget
    % change from hgsetget to handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    % 20140910 Make all functions more modular, ie, use the same variable to
    % store the average values obtained by different methods
    
    properties
        %% Image properties
        path = ''; % can be either a path or a image matrix
        folderPath = '';
        frames = 0;
        referenceIm = [];
        frame = 1;
        width = 0;
        height = 0;
        SR = ''; % original sampling rate
        rsSR = ''; % resampled sampling rate
        rsSampleRate = '';
        color = 'Gray';
        cmax = 1000;
        cmin = 0;
        frameShift = [0 0]; % first entry is the frame, second is the shifted row, third is shifted column
        % Reference frame for translational correction
        %         image = [];
        imHandle = '';
        patch = '';
        % patch data for easier manipulation
        
        initialData = [];
        % Store referenced images into matrix for faster access and
        % manipulation
        
%         RBBgSubData = [];

        cDate = '';
        % experiment date
        pClampNo = [];
        imageNo = [];
        
        sweepFrames = [];
        % number of frames per stimulation sweep. This should be
        % automatically determined by the number of frame triggers
        % recorded. However, in cases such as 20160412 one, the last two
        % frame triggers were not recorded, so have to manually put in this
        % value.
        
        %% ROIs
        
        ROI = [];
        ROICentre = [];
        sigROI = {};
        % 1s row is the BW mask (0 for all entries and 1 for circled
        % region). 2nd row is the circled points's position (rounded to
        % integer). 3rd row is the centre position (rounded to integer). If
        % the ROI method is convolution, then 2nd and 3rd values will be
        % the same.
        
        electrodeROI = [];
        
        tempSigROI = {};
        % To store ROIs that have been temporarily deleted for analysis
        bgROI = {};
        % background ROI. Note that it is defined as the 1st selected sigROI.
        
%         ROIN = [];
        % current ROI number for getting the frame time (using line time to
        % get a more accurate timing of the scanning
        
        refROI = {};
        % Reference ROI for mapping onto different recordings that imaged
        % the same imaging field
        
        refVertMov = 0;
        refHoriMov = 0;
        % The translation movement to map the reference ROI onto this
        % CaObj. ith column corresponds to ith reference ROI.
        
        refROINum = 0;
        % reference ROI number used
        
        refROIDate = [];
        % Date that reference ROI number is made
        
        refROITime = [];
        % time that reference ROI number is made
        
        varMat = [];
        % matrix containing the variance of each pixel
        
        %% Extracellular stimulation

        % number of stimulation pulses/bursts
        stimNo = [];
        stimOnFrameCell = {};
        % similar to LEDOnFrameCell, where the 1st row contains the frame
        % number and 2nd row contains the frame time during stimulation.
        % Additionally, the data points and data times are contained in row
        % 3 and row 4 respectively.
        stimOffFrameCell = {};
        % similar to LEDOffFrameCell, where the 1st row contains the frame
        % number and 2nd row contains the frame time after stimulation.
        % Additionally, the data points and data times are contained in row
        % 3 and row 4 respectively.
        
        stimCoupleWindow = 0.2;
        % 200ms as the determination of whether a Ca response is couple to
        % electrical response or not. Additional frames are to be added
        % after stimFrame for finding out the actual maximum amplitude of
        % calcium signal.
        stimFrame = {};
        % 20140906 definition
        % 1st row cells, 1st row = baseline frames, 2nd row = baseline frame time
        % The baseline value will be 100ms before the onset of the stimulation
        % (approximately 10 frames for 100Hz, 2 frames for 20Hz). If the
        % imageing frequency is too slow (<20Hz), then 500ms before the
        % onset should be used.
        % 20141001 additional
        % 4th row cells
        % 1st row = baseline data points(for Clampex data)
        % 2nd row = baseline data time (for Clampex data)
        
        % 2nd row cells, 1st row = onset frames, 2nd row = onset frame time
        % The onset frames will simply be the time when the stimulation is
        % on.
        % 20141001 additional
        % 5th row cells
        % 1st row = onset data points(for Clampex data)
        % 2nd row = onset data time (for Clampex data)
        
        % 3rd row cells, 1st row = offset frames, 2nd row = offset frame time
        % Offset response should be at least 5
        % second, if in cases where the offset is disrupted by another
        % stimulation, then just use the time between the stimulations as
        % the offset time.
        % 20141001 additional
        % 6th row
        % 1st row = offset data points(for Clampex data)
        % 2nd row = offset data time (for Clampex data)
        
        % 7th row = baseline + onset + offset, frames and frame time
        % 8th row = baseline + onset + offset, data pts and data time
        
        stimFilePath
        
        segStimFrame = {};
        
        segStimTime = {};
        
        stimSegFrameBs = {};
        stimSegFrameTimeBs = {};
        stimSegFrameON = {};
        stimSegFrameTimeON = {};
        stimSegFrameOFF = {};
        stimSegFrameTimeOFF = {};
        stimSegFrameWhole = {};
        stimSegFrameTimeWhole = {};
        
        DFFSegSigBs = {};
        DFFSegSigON = {};
        DFFSegSigOFF = {};
        
        riseFrame = {};
        % Frames at which the signal is rising. It is organised into nROI x
        % nStim cell format
        decayFrame = {};
        % Frames at which the signal is decaying. It is organised into nROI x
        % nStim cell format
        stimVar = 0;
        % The parameter that is being varied.
        % 1 = amplitude, 3 = duration and 13 = reps.
        
        extStimSegAmp = [];
        % Extracellular stimulation amplitude for each segmented signal
        
        extStimSegDur = [];
        % Extracellular stimulation duration for each segmented signal
        
        extStimSegFreq = [];
        % Extracellular stimulation frequency for each segmented signal
        
        extStimSegRep = [];
        % Extracellular stimulation repeat for each segmented signal
        
        extStimSegBurstInt = [];
        % Extracellular stimulation burst interval for each segmented
        % signal
        
        extStimVar = '';
        % The variable in extracellular stimulation
        
        
        %% Stimulation properties
        % For properties that stay the same throughout all stimulations,
        % only 1 value is given. If the properties change in different
        % stimnulation number, then the variable becomes a vector that
        % holds values for each stimulation pulse/burst
        stimPro = [];
        % the stimulation protocol number applied
        duration = [];
        reps = [];
        amp = [];
        freq = [];
        
        
        
        %% LED stimulation
        nLED = [];
        % number of LED stimulations (every ON and OFF counted as one).
        bgLEDRemoveInt = [];
        % LED artefact removed convoluted background intensities
        LEDOnFrameCell = {};
        % LED ON cell form
        LEDOnMFrame = [];
        % LED ON matrix form
        % I divided ON frame into middle part (M) and edge part (Edge) because the
        % edges sometimes have incomplete illumination and causes incorrect
        % LED ON intensity calculation
        LEDOnEdgeFrame = [];
        % The edges of the LED On time. The values of edges are calculated
        % as the average of its neighbours (+2 and -2 points).
        LEDOffFrameCell = {};
        % LED OFF cell form
        LEDOffFrame = [];
        % LED OFF matrix form
        LEDStimTime = [];
        
                %% LED stimulation
       
        LEDONFrame = {};
        % LED ON cell form
        
        LEDOFFFrame = {};
        % LED OFF cell form
        
        lightONTime = 0;
        lightOFFTime = 0;
        
        ONInc = [];
        % Light ON increase in response indicator
        ONDec = [];
        % Light ON decrease in response indicator
        OFFInc = [];
        % Light OFF increase in response indicator
        OFFDec = [];
        % Light OFF decrease in response indicator
        
        lightPVal = '';
        
        RGCType = {};
        % row array where column number is the ROI number 
        
        % Light ON or OFF response indicator, inc means increase, dec means
        % decrease, in ON or OFF segment
        
%         LEDOnMFrame = [];
        % LED ON matrix form
        % I divided ON frame into middle part (M) and edge part (Edge) because the
        % edges sometimes have incomplete illumination and causes incorrect
        % LED ON intensity calculation
%         LEDOnEdgeFrame = [];
        % The edges of the LED On time. The values of edges are calculated
        % as the average of its neighbours (+2 and -2 points).

%         LEDOffFrame = [];
        % LED OFF matrix form
        
        
        
        %% Extracted signals
        rawSig = [];
        % the values of each signal over time, in averaged form (if using
        % convolution method, it's the value of the local maximal points.
        % If using manual ROI selection method, rawSig is the average of
        % the pixels within the manually drawn boundary.
        rawSegSig = {};
        % nROI x nStim cell
        % DF/F0 signal

        
        SGSig =[];
        % Savitzky-Golay filtered signal
        SGSegSig = {};
        % Savitzky-Gola filtered signal organised into nRO x nStim cell
        % format.
        
        normSig = [];
        % normSig = normalise rawSig = (rawSig - baseline)/baseline
        
        LEDReSig = [];
        % LED artefact removed signal
        sigBl = [];
        % baseline value (found by median) of each cell
        bgVal = [];
        % the background values for subtraction for each frame,1 x nFrames
        bgBl = [];
        % baseline value of the background
        
        %% Segmented signals
        % Signals that are segmented into corresponding stimulation
        % signals. They are organised into cell format and the cell matrix
        % is organised into nROI x nStim, each entry's row and column
        % represents the ROI number and stimulation number of that
        % segmented signal
        blSig = [];
        % baseline signal
        stimSig = [];
        % signal during stimulation time
        offSig = [];
        % signal after stimulation
        riseSig = [];
        % signal that is rising from baseline to peak
        decaySig = [];
        % signal from peak to baseline.
        wholeSig = [];
        % overall signal in one stimulation
        stimAndOffSig = [];
        % stimulation and offset signals combined
        
        checkSig = [];
        
        blStd = []
        % baseline standard deviation
        
        
        
        %% Calcium signal properties
        %         calSig = [];
        %         % calSignal is in structure form, where each stimulation represents
        %         % one instance.
        %         % Calcium signal is divided into 3 sections: 1. Baseline values
        %         % (baseline frames) 2. Rise time values (n frames before the peak)
        %         % 3. Decay times values (m frames after the peak)
        %         % Calcium signals, either due to LED stimulation or electrical
        %         % stimulation. The method is to find the stimulation time first
        %         % (for LED, it will be divided into LED ON or OFF frames. For
        %         % electrical stim, it will be stim onset and offset frames).
        %         stimSegSig = {};
        %         % Stimulation coupled signals in segmented form. Each column
        %         % contains the
        
        % Each variable is a nROI x nStim matrix
        riseTime = [];
        % full rise time from 0% to 100%
        halfRiseFrame
        % half rise frame from 0% to 50%
        halfRiseTime = [];
        % half rise time from 0% to 50%
        rise10to90Time = [];
        % rise time from 10% to 90%
        
        decayTime = [];
        aboveThreshTime = [];
        % The calcium signal threshold is defined to be 3 times the
        % baseline noise level. This is to find when the cell starts
        % responding and potentially find how many stimulation repetition
        % is needed for calcium signal to occur.
        maxAmp = [];
        maxAmpLocs = [];
        offsetTime = [];
        blTime = [];
        blNoiseLevel = [];
        
        blMean = [];
        peakCaAmp = [];
        peakFall = [];
        HDT = [];
        bsReturn = [];
        
        % noise level during baseline, calculated as standard deviation
        stimCoupledRes = [];
        % A logical array indicating whether there is an electrical
        % stimlation coupled response or not. A response is determined when
        % signal maximum is 3 times the value of baseline noise (std value)

        %% Save and load
        gROI = [];
        % responsive ROIs to the stimulation
        lsROI = [];
        % Cells that are responsive initially to laser scanning
        
        
        %% Image processing
        mask = []; % static segmentation method adapted from Tomek
        translationCorrectFlag = 0;
        initialTwoDData = [];
        backgroundSubtracted = [];
        temporalFiltered = [];
        %% Convolution method
        convIm = [];
        % Convoluted image from convolution function
        convR = '';
        % convolution item size
        %% Intracellular stimulation (Not many experiments yet 20141002)
        APstimTime = [];
        framesPerSweep;
        % number of frames recorded for each pClamp sweep
        % This is also for Thorlab camera where the number of frames for each
        % sweep can be set. This has to be manually assigned.
        
        intStimAmp = [];
        % For recording intracellular stimulation amplitudes
        
        intStimDur = [];
        % For recording intracellular stimulation durations
        
        
        %% Fit data parameters
        riseFitExpCoef = {};
        % Exponential 2 fitting curve parameters for rising fluorescence
        % The variable is in the form of nROI x nStim
        % where a structured variable is in each entry and the structured
        % variable contains a b c and d, the coefficients of the
        % exponential 2 fit
        decayFitExpCoef = {};
        % Exponential 2 fitting curve parameters for decaying fluorescence
        % The variable is in the form of nROI x nStim
        % where a structured variable is in each entry and the structured
        % variable contains a b c and d, the coefficients of the
        % exponential 2 fit
        
        riseFitLinCoef = {};
        
        RDFitExpCoef = {};
        % fitting both rise and decay time simultaneously with 2
        % exponentials
        RDFitExpGOF = {};
        
        fitModelCoef = {};
        % The model coefficients fitted by Calvin Eiber method for fitting
        % Calcium signals. Organised into cell format with nROI x nStim entries
        fitModelR2 = [];
        % R2 values for Calvin Eiber Calcium signal fitting. Organised with
        % nROI x nStim entries
        CaModel = [];
        % Calvin Eiber Ca signal model
        
        riseFitExpGOF = {}
        % Exponential 2 fitting curve Goodness of fit for rising fluorescence
        % The variable is in the form of nROI x nStim
        % where a structured variable is in each entry and the structured
        % variable contains SSE, R2, dfe, adjusted R2 and RMSE
        decayFitExpGOF = {}
        % Exponential 2 fitting curve Goodness of fit for decaying fluorescence
        % The variable is in the form of nROI x nStim
        % where a structured variable is in each entry and the structured
        % variable contains
        
        decayTau = [];
        % Fitted decay singals with time constants (tau) recorded; the time
        % constant values are stroed in the usua nROI x nStim way.
        
        % The fitting is either exp or exp2, and it should be figured out
        % by the number of coefficients
        
        segFitInfo = {};
        fitR2 = [];
        
        %% Unsorted
        plotData = [];
        % A temporary variable that contains the data for last plot. This
        % is used for rescaling stimulus to fit well with the calcium
        % signals
        nRS = '' % number of resampled average intensities data point
        
        cellPixelPercent = 0;
        % percentage of cell occupying the imaging field
        
        
        segSigEntry = {};
        
        %         stimTime = {};
        %         % Extracellular stimulation voltage recorded by Clampex. Both long
        %         % pulse stimulation or burst of several short pulses stimulation
        %         % are separated by certain set amount of interval duration.
        %         % each pulse/burst is stored in a column
        %         % 1st row contains the data point of stimulation start time
        
        
        alpha = [];
        threshF = [];
        minGrowth = [];
        
        coupledResponse = [];
        
        segSigPeakVal = [];
        segSigPeakFrame = [];
        % Frame at which peak occurs counting from stimluation trigger
        segSigHDT = [];
        
        nCoupledSig = [];
        
        
        avgSig;
        avgSigSTD; % standard deviation
        avgSigSE; % standard error
        % These are for episodic recordings of averaged trial recordings
        
        avgSegSig;
        avgSegSTD;
        avgCoupledSig;
        
        avgSigPeakVal;
        avgSigHDT;
        avgSigPeakFrame;
        
        cROIN = 1;
        cStimN = 1;;
        % current pointing ROI number and stimulation number
        
        
        fontSize = 14;
        legendSize = 14;
        tickLabelSize = 14;
        fontWeight = 'bold';
        
        %   Copyright 2014 John Yang
        
        %% Raw signal from the original picture will go through several signal processing
        %  1. Rolling ball background subtraction to correct uneven illumination
        %  2. Background subtraction
        RBBgSubData = [];
        % Background subtraction by rolling ball background subtraction.
        % The algorithm is: First, do top-hat filtering with a ball
        % structure (by strel), second, subtract the original by the ball
        % filtered image
        RBBgSubSig = [];
        bgSubSig = [];
        % Simple background subtraction
        
        stimType;
        
        
%         framesPerSweep;
        
        % For patched RGC calcium signal analysis
        patchedCaBlMean = [];
        patchedCaSigDFF = [];
        patchedCaArea = [];
%         patchedCaAreaPerAmp = [];
        patchedFrameTime = [];

            
        electrodeXShift = 0;
        electrodeYShift = 0;
        % These are for electrode position adjustment for electrode outside
        % of the screen, but still able to locate the position of the electrode
    end
    
    properties (Dependent)
        image;
        % nROI x nStim cell
        fps;
        segImageTime;        
        % The timing for each frame
        
        %  3. Bleaching correction
        %  4. DF/F0 signal
        %  5. DF/F0 signal in segments
        %  This processed signals will be in order and subsequent function
        %  will call the previous process signals only to make this
        %  processing a production line like.

        bgSubBlCorrSig = [];
        
        rawBlCorrSig = [];
        DFFSig = [];
        DFFSegSig = {};
        
        
        

        
        
        frameTime;
        
        % Image time steps        
        cellFrameTime;        
        % The corrected image time for each cell, since every line is
        % scanned at different time
        
        nStim;
        
        CEDistance;
        % cell to electrode distance;
        
        coupledSigProb;

        nROI;
                
        segFrame;
        stimStartFrame; 
        blSegFrame;
        
        segFrameTime;
    end
    

    
    properties (Access = private)
        map = 'Gray';
    end
    
    %     methods (Abstract)
    %     end
    
    
    methods (Static)
        function filePath = findFile(folderPath,num,varargin)
            % This static function is mainly for automatically anaylse the
            % pClamp data. It won't output the object instance in the workspace
            % but will save the analysed results and objects into the folder.
            % If the user want to look at the instance and its properties,
            % he/she should use the initialisation function
            
            % STST is start-to-start time. This value needs to be input in
            % by the user if pClamp recording is in episodic mode.
            % Otherwise accurate timing of each timing cannot be calculated
            
            dyeNum = [];
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'dyeNum')
                    dyeNum = varargin{i + 1};
                end
            end
            
            if ~isempty(dyeNum)
                if num < 10
                    findFile = sprintf('\\0%d*', num);
                elseif num >= 10
                    findFile = sprintf('\\%d*', num);
                end
                
                tempFileStruct = dir([folderPath,findFile]);
                for i = 1:length(tempFileStruct)
                    if tempFileStruct(i).isdir == 1
                        tempFilePath = fullfile(folderPath,tempFileStruct(i).name);
                        break
                    end
                end
                
                
                
                findDye = sprintf('\\*PMT [%d]*', dyeNum);
                fileStruct = dir([tempFilePath,findDye]);
                filePath = fullfile(tempFilePath,fileStruct.name);
            else
                
                findFile = sprintf('\\*%d*.tif', num);
                fileStruct = dir([folderPath,findFile]);
                for i = 1:length(fileStruct)
                    filePath{i} = fullfile(folderPath,fileStruct(i).name);
                end
            end
            
        end
    end
    
    methods
        function obj = CaImage(path, pClampNo, imageNo)
            obj.path = path;
            [obj.folderPath,~,~] = fileparts(path);
            
            if nargin > 1
                obj.pClampNo = pClampNo;
            elseif nargin > 2
                obj.pClampNo = pClampNo;
                obj.imageNo = imageNo;
            end
            
            
            if ischar(path)
                
                temp = imfinfo(obj.path);
                obj.frames = length(temp);
                obj.width = temp(1).Width;
                obj.height = temp(1).Height;
                obj.frameShift = zeros(obj.frames,2);
            else
                obj.frames = size(obj.path,3);
                obj.width = size(obj.path,2);
                obj.height = size(obj.path,1);
            end
            
        end
        
        function val = get.image(obj)
            if ischar(obj.path)
                val = imread(obj.path, obj.frame);
            else
                val = obj.path(:,:,obj.frame);
            end
            
        end
        
        %         function val = get.imageSR(obj)
        %             if ischar(obj.pathSR)
        %                 val = imread(obj.pathSR, obj.frame);
        %             else
        %                 val = obj.pathSR(:,:,obj.frame);
        %             end
        %
        %         end
        
        function draw(obj, HHandle)
            if (nargin > 1)
                obj.imHandle = imagesc(obj.image, 'Parent', HHandle);
            elseif ~isempty(obj.ROI)
                set(obj.imHandle,'CData',obj.image);
            else
                obj.imHandle = imagesc(obj.image);
            end
            colormap(obj.map);
            caxis([obj.cmin, obj.cmax])
            colorbar;
            colorMatrix = ['y', 'm', 'c', 'r', 'b', 'w'];
            if ~isempty(obj.ROI)
                hold on
                for j = 1:length(obj.ROI)
                    b = obj.ROI{j};
                    t = obj.ROICentre{j};
                    colorMod = mod(j,length(colorMatrix));
                    if colorMod == 0
                        colorMod = length(colorMatrix);
                    end
                    plot(b(:,2),b(:,1),colorMatrix(colorMod),'LineWidth',1);
                    text(t(1), t(2), num2str(t(3)), 'Color',colorMatrix(colorMod),...
                        'horizontalalignment','c', 'verticalalignment','m')
                    
                end
            end
            hold off
        end
        
        function fps = get.fps(obj)
            if isa(obj.patch,'patchData')
                
                if length(obj.patch.STSFrameATime) == 1
                    fps = 1/obj.patch.STSFrameATime;
                else
                    fps = 1/obj.patch.STSFrameATime{1};
                end
                
            else
                fps = obj.patch;
            end
            
        end
        
        function frameTime = get.frameTime(obj)
            if isa(obj.patch,'patchData')
%                 if (~isempty(obj.patch.frameTime))&(obj.patch.ndim ~= 3)&(length(obj.patch.frameTime) == obj.frames)
                if (~isempty(obj.patch.frameTime))
                    frameTime = obj.patch.frameTime;
%                 else
%                     frameTime = 1/obj.fps:1/obj.fps:(obj.frames*1/obj.fps);
                end
            else
                
                frameTime = 1/obj.fps:1/obj.fps:(obj.frames*1/obj.fps);
            end
        end
%         function setSR(obj, patch)
%             if isa(patch,'patchData')
%                 obj.fps = 1/patch.STSFrameATime;
%                 if (~isempty(patch.frameTime))&(patch.ndim ~= 3)&(length(patch.frameTime) == obj.frames)
%                     obj.frameTime = patch.frameTime;
%                 else
%                     obj.frameTime = 1/obj.fps:1/obj.fps:(obj.frames*1/obj.fps);
%                 end
%             else
%                 obj.fps = patch;
%                 % patch is actually the frame rate in this case
%                 obj.frameTime = 1/obj.fps:1/obj.fps:(obj.frames*1/obj.fps);
%             end
%         end
%         
        
        
        
        
        function set.color(obj, val)
            m = 256;
            obj.map = zeros(m, 3);
            switch val
                case 'Gray'
                    obj.map = 'Gray';
                    obj.color = 'Gray';
                case 'Red'
                    obj.map(:,1) = linspace(0,1,m);
                    obj.color = 'Red';
                case 'Green'
                    obj.map(:,2) = linspace(0,1,m);
                    obj.color = 'Green';
                case 'Blue'
                    obj.map(:,3) = linspace(0,1,m);
                    obj.color = 'Blue';
                otherwise
                    obj.map = 'Gray';
                    obj.color = 'Gray';
            end
            
        end
        
        function findcValue(obj)
            originalFrame = obj.frame;
            for i = 1: obj.frames
                obj.frame = i;
                currentCmax = max(max(obj.image));
                currentCmin = min(min(obj.image));
                if obj.cmax < currentCmax
                    obj.cmax = currentCmax;
                end
                if obj.cmin > currentCmin
                    obj.cmin = currentCmin;
                end
            end
            obj.frame = originalFrame;
            
        end
        
        function findcValueSR(obj)
            originalFrame = obj.frame;
            for i = 1: obj.frames
                obj.frame = i;
                currentCmax = max(max(obj.imageSR));
                currentCmin = min(min(obj.imageSR));
                if obj.cmax < currentCmax
                    obj.cmax = currentCmax;
                end
                if obj.cmin > currentCmin
                    obj.cmin = currentCmin;
                end
            end
            obj.frame = originalFrame;
            
        end
        
        
        function plotCollapsed(obj)
            tempField = zeros(obj.height,obj.width);
            
            for i = 1:obj.frames
                obj.frame = i;
                tempField = tempField + double(obj.image());
            end
            
            imagesc(tempField);
            
            set(gca,'YDir','normal')
            colormap(obj.map);
            caxis([5000 8000])
            colorbar;
            
        end
        %         function frameS(obj,frame,rowShift,colShift)
        %             obj.frameShift = [frame,rowShift,colShift];
        %         end
        function field = output3D(obj)
            field = zeros( obj.height,obj.width,obj.frames);
            
            for i = 1:obj.frames
                obj.frame = i;
                field(:,:, i) = double(obj.image());
            end
        end
        
        function transCorr(obj, fileType, saveIm)
            
            % fileType: 1 for path, 2 for matlab data and 3 for handle
            % save: 1 if the user wants to save the image in tif form
            
            [obj.path, obj.frameShift] = CellPreprocessing2(obj.path, fileType,saveIm);
            
            obj.translationCorrectFlag = 1;
            
        end
        
        function imageShift(obj, shift) % shift is in 2(x and y image)x number of frames
            
            tempShiftIm = zeros(obj.width, obj.height, obj.frames);
            for i = size(shift,1)
                obj.frame = i;
                tempOrgIm = obj.image;
                if (shift(i,1) <= 0)&&(shift(i,2) <= 0)
                    tempShiftIm(1:(end + shift(i,1)),1:(end + shift(i,2)),i) = tempOrgIm(1:(end + shift(i,1)),1:(end + shift(i,2)));
                    
                elseif (shift(i,1) >= 0)&&(shift(i,2) <= 0)
                    tempShiftIm((1 + shift(i,1)):end,1:(end + shift(i,2)),i) = tempOrgIm((1 + shift(i,1)):end,1:(end + shift(i,2)));
                    
                elseif (shift(i,1) <= 0)&&(shift(i,2) >= 0)
                    tempShiftIm(1:(end + shift(i,1)),(1 + shift(i,2)):end,i) = tempOrgIm(1:(end + shift(i,1)),(1 + shift(i,2)):end);
                    
                elseif (shift(i,1) >= 0)&&(shift(i,2) >= 0)
                    tempShiftIm((1 + shift(i,1)):end,(1 + shift(i,2)):end,i) = tempOrgIm((1 + shift(i,1)):end,(1 + shift(i,2)):end);
                    
                end
                
            end
            
            obj.path = tempShiftIm;
            
            obj.translationCorrectFlag = 1;
            
            obj.frameShift = shift;
            
        end
        
        function thresholding(obj,frame,threshValue,otherMask)
            if nargin < 4
                obj.frame = frame;
                obj.mask = obj.image > threshValue;
                usedMask = obj.mask;
            else
                usedMask = otherMask;
            end
            
            obj.rawSig = zeros(1, obj.frames);
            for i = 1:obj.frames
                obj.frame = i;
                obj.rawSig(i) = sum(sum(obj.image(usedMask)))/sum(sum(usedMask));
            end
        end
        
        function temporalFilter(obj,filterType,avgFrame)
            % filterType 1. average 2. min 3. median
            obj.temporalFiltered = zeros(obj.width,obj.height, obj.frames);
            InitialMatrix(obj)
            
            for i = 1:obj.frames
                switch filterType
                    case 'average'
                        if (i <= avgFrame)
                            obj.temporalFiltered(:,:,i) = sum(obj.InitialData(:,:,i:avgFrame+i-1),3)/avgFrame;
                        elseif ((i + avgFrame) >= obj.frames)
                            obj.temporalFiltered(:,:,i) = sum(obj.InitialData(:,:,i-avgFrame:i),3)/avgFrame;
                        else
                            obj.temporalFiltered(:,:,i) = sum(obj.InitialData(:,:,(i - avgFrame):(i + avgFrame)),3)/(2*avgFrame+1);
                        end
                    case 'min'
                        if (i <= avgFrame)
                            obj.temporalFiltered(:,:,i) = min(obj.InitialData(:,:,i:avgFrame+i-1), [], 3);
                        elseif ((i + avgFrame) >= obj.frames)
                            obj.temporalFiltered(:,:,i) = min(obj.InitialData(:,:,i-avgFrame:i), [],3);
                        else
                            obj.temporalFiltered(:,:,i) = min(obj.InitialData(:,:,(i - avgFrame):(i + avgFrame)), [],3);
                        end
                        
                    case 'median'
                        if (i <= avgFrame)
                            obj.temporalFiltered(:,:,i) = median(obj.InitialData(:,:,i:avgFrame+i-1),3);
                        elseif ((i + avgFrame) >= obj.frames)
                            obj.temporalFiltered(:,:,i) = median(obj.InitialData(:,:,i-avgFrame:i),3);
                        else
                            obj.temporalFiltered(:,:,i) = median(obj.InitialData(:,:,(i - avgFrame):(i + avgFrame)),3);
                        end
                end
            end
        end
        
        
        function readStim(obj,year)
            if strcmp(obj.patch.stimType, 'intStim')
                obj.intStimAmp = obj.patch.intStimAmp;
                % For recording intracellular stimulation amplitudes
                
                obj.intStimDur = obj.patch.intStimDur;
                % For recording intracellular stimulation durations
            elseif strcmp(obj.patch.stimType, 'extStim')
                
                if nargin > 1
                    if year == 2014
                        
                        %                 findFile = sprintf('\\%d*', num);
                        findFile = '\*.txt';
                        % txt file contains the ext stimulation information
                        
                        tempFileStruct = dir([obj.folderPath,findFile]);
                        
                        %                 for i = 1:length(tempFileStruct)
                        %                     if tempFileStruct(i).isdir == 1
                        %                         tempFilePath = fullfile(folderPath,tempFileStruct(i).name);
                        %                         break
                        %                     end
                        %                 end
                        
                        extStimFilePath = fullfile(obj.folderPath,tempFileStruct.name);
                        
                        obj.readExtStim(extStimFilePath);
                        
                    else
                        obj.readExtStim(obj.stimFilePath);
                        obj.nStim = length(obj.extStimSegAmp);
                    end
                end
                
            end
        end
        
        function readExtStim(obj,fn)
            fid = fopen(fn);
            opts = struct;
            
            while ~feof(fid)
                
                lineIn = strtrim(fgetl(fid));
                
                if isempty(lineIn) || ~max(lineIn == ':') || lineIn(1) == '%'
                    continue
                end
                
                [n,v] = strtok(lineIn,':');
                [v,u] = strtok(v(2:end),'#');
                
                n = lower(strtrim(n));
                if strfind(n,'forread')
                    varVal = str2num(strtrim(v));
                    varVal(varVal == 0) = [];
                    if strfind(n,'amp')
                        obj.extStimVar = 'amp';
                    elseif strfind(n,'dur')
                        obj.extStimVar = 'dur';
                    elseif strfind(n,'rep')
                        obj.extStimVar = 'rep';
                    end
                else
                    v = eval(strtrim(v));
                end
                
                
                opts.(n) = v;
                if ~isempty(u),
                    u = strtrim(u(2:end));
                    opts.([n '_unit']) = u;
                end
            end
            
            fclose(fid);
            clear lineIn n u v fid
            %
            % Fix units to default units
            %
            
            if strcmp(opts.amplitude_unit, 'ma')
                opts.amplitude = opts.amplitude * 1000;
            end
            opts.amplitude_unit = 'ua';
            
            if strcmp(opts.duration_unit, 'ms')
                opts.duration = opts.duration * 10^3;
            elseif strcmp(opts.duration_unit, 's')
                opts.duration = opts.duration * 10^6;
            end
            opts.duration_unit = 'us';
            opts.duration = round(opts.duration);
            
            if strcmp(opts.burst_interval_unit, 'ms')
                opts.burst_interval = opts.burst_interval * 10^3;
            elseif strcmp(opts.burst_interval_unit, 's')
                opts.burst_interval = opts.burst_interval * 10^6;
            end
            opts.burst_interval_unit = 'us';
            opts.burst_interval = round(opts.burst_interval);
            
            if isempty(varVal)
                varVal = ones(1,opts.burst_repeat);
            end
            
            obj.extStimSegAmp = opts.amplitude(1)*ones(1,length(varVal));
            obj.extStimSegDur = opts.duration(1)*ones(1,length(varVal));
            obj.extStimSegFreq = opts.frequency(1)*ones(1,length(varVal));
            obj.extStimSegRep = opts.repeat(1)*ones(1,length(varVal));
            obj.extStimSegBurstInt = opts.burst_interval(1)*ones(1,length(varVal));
            
            if strcmp(obj.extStimVar, 'amp')
                obj.extStimSegAmp = varVal;
            elseif strcmp(obj.extStimVar, 'dur')
                obj.extStimSegDur = varVal;
            elseif strcmp(obj.extStimVar, 'rep')
                obj.extStimSegRep = varVal;
            end
            
            %             if size(opts.amplitude,2) > 1
            %
            %             elseif size(opts.duration,2) > 1
            %
            %             elseif size(opts.repeat,2) > 1
            %
            %             else
            %
            %             end
            
        end
        
        function calSigSort(obj, ROIN)
            % removal artefacts (for example, LED artefact) and then get
            % the calcium signals from stimFrame (sorted into baseline,
            % onset and offset frames.
            
            % Calcium signal is divided into 3 sections: 1. Baseline values
            % (baseline frames) 2. Rise time values (n frames before the peak)
            % 3. Decay times values (m frames after the peak)
            % Calcium signals, either due to LED stimulation or electrical
            % stimulation. The method is to find the stimulation time first
            % (for LED, it will be divided into LED ON or OFF frames. For
            % electrical stim, it will be stim onset and offset frames).
            
            obj.calSig = [];
            
            obj.stimFrameSort(ROIN);
            obj.LEDRemoval(ROIN);
            stimData = obj.LEDReSig;
            
            stimData = smooth(stimData, 5);
            
            for stimNo = 1:obj.nLED
                
                obj.calSig(stimNo).bsSig = stimData(obj.stimFrame{1,stimNo}(1,:));
                % baseline frame signal
                obj.calSig(stimNo).stimSig = stimData(obj.stimFrame{2,stimNo}(1,:));
                % stim onset frame signal
                obj.calSig(stimNo).offsetSig = stimData(obj.stimFrame{3,stimNo}(1,:));
                % stim offset frame signal
                
                obj.calSig(stimNo).overallSig = [obj.calSig(stimNo).bsSig',...
                    obj.calSig(stimNo).stimSig', obj.calSig(stimNo).offsetSig'];
                
                obj.calSig(stimNo).bsAvg = mean(obj.calSig(stimNo).bsSig);
                % baseline average value before the onset of the
                % stimulation
                
                if length(obj.calSig(stimNo).bsSig) >= 3
                    obj.calSig(stimNo).bsNoiseSD = std(obj.calSig(stimNo).bsSig);
                elseif length(obj.calSig(stimNo).bsSig) >= 2
                    obj.calSig(stimNo).bsNoiseSD = abs(diff(obj.calSig(stimNo).bsSig));
                end
                % standard deviation as described in the paper: noise std
                % at baseline time
                
                %                 [obj.calSig(stimNo).peakAmp, obj.calSig(stimNo).peakLoc] = ...
                %                     max(obj.calSig(stimNo).overallSig);
                
                
                %                 peakThresh = obj.calSig(stimNo).bsAvg + 3*obj.calSig(stimNo).bsNoiseSD;
                %                 [obj.calSig(stimNo).peakAmp, obj.calSig(stimNo).peakLoc] = ...
                %                     findpeaks(obj.calSig(stimNo).overallSig, 'MINPEAKDISTANCE', 10, 'MINPEAKHEIGHT', peakThresh);
                %                 % maybe should look for local peaks
                %
                %                 if isempty(obj.calSig(stimNo).peakAmp)
                
                
                [obj.calSig(stimNo).peakAmp, obj.calSig(stimNo).peakLoc] = ...
                    max(obj.calSig(stimNo).overallSig);
                
                
                % Finding the peak in each stim calcium signal
                obj.calSig(stimNo).overallPeakLoc = obj.stimFrame{1,stimNo}(1,1) + obj.calSig(stimNo).peakLoc - 1;
                % the peak frame number in overall frames
                
                
                if length(obj.calSig(stimNo).bsSig) >= 2
                    obj.calSig(stimNo).checkSig = (obj.calSig(stimNo).peakAmp(1) - obj.calSig(stimNo).bsAvg) > 3*obj.calSig(stimNo).bsNoiseSD;
                else
                    obj.calSig(stimNo).checkSig = -1;
                end
                % CheckSig is to determine whether a valid signal occurs or
                % not. A signal is defined by Xiaowe Chen 2011 by the
                % following definition: A fluorescent change was accepted
                % as a calcium signal when its amplitude (peak value) was
                % three times larger than the standard deviation of the
                % noise that was determined for a period of 100ms just
                % before auditory stimulation.
                % 0 = no signal, 1 = signal is present and -1 =
                % undetermined (too less baseline data to get noise level).
                
                
                
                obj.calSig(stimNo).DFFAllSig = ...
                    (obj.calSig(stimNo).overallSig - obj.calSig(stimNo).bsAvg)/obj.calSig(stimNo).bsAvg;
                
                
                
                obj.calSig(stimNo).sigOnLatency = obj.frameTime(obj.calSig(stimNo).overallPeakLoc(1)) - obj.stimFrame{2,stimNo}(2,1);
                % the latency between the onset of the stimulation
                % and the peak response
                obj.calSig(stimNo).sigOffLatency = obj.frameTime(obj.calSig(stimNo).overallPeakLoc(1)) - obj.stimFrame{2,stimNo}(2,end);
                % the latency between the offset of the stimulation
                % and the peak response
                
                obj.calSig(stimNo).sigArea = trapz(obj.calSig(stimNo).overallSig*obj.fps);
                % sigArea = fluorescent signal area. It is computed with
                % trapzoidal area method, with sampling rate as the time
                % step.
                
            end
            
        end
        
        
        
%         function saveRollingBallBgSub(obj, savePath, R)
%             
%             if nargin < 3
%                 SE = strel('disk', 20, 0);
%             else
%                 SE = strel('disk', R, 0);
%             end
%             
%             obj.RBBgSubData = zeros(obj.height,obj.width,'uint16');
%             TiffObj = Tiff(obj.path,'r'); % open file for reading
%             
%             
%             TiffObj.setDirectory(1);
%             tempIm = TiffObj.read();
%             
%             imFilt = imopen(tempIm,SE);
%           
%             tempRBBgSubData = tempIm - imFilt;
%             
%             imwrite(tempRBBgSubData, savePath);
%             
%             clear TiffObj
%             
%         end
        
        
        
        function rollingBallBgSub(obj, varargin)
            % background subtraction by rolling ball method
            % R = radius size
            
            R = 20;
            savePath = '';
            
            for i = 1:length(varargin)
                
                if strcmp(varargin{i}, 'R')
                    R = varargin{i + 1};
                elseif strcmp(varargin{i}, 'savePath')
                    savePath = varargin{i + 1};
                
                elseif strcmp(varargin{i}, 'framesToCompute')
                    framesToCompute = varargin{i + 1};
                end

            end
          
%             if nargin < 2
%                 SE = strel('disk', 20, 0);
%             else
%                 SE = strel('disk', R, 0);
%             end
            
            SE = strel('disk', R, 0);
            
            frameN = length(framesToCompute);

            obj.RBBgSubData = zeros(obj.height,obj.width,frameN,'uint16'); % 3D matrix containing 16bits values
            TiffObj = Tiff(obj.path,'r'); % open file for reading

            TiffObj.setDirectory(1);
            tempIm = TiffObj.read();
            imFilt = imopen(tempIm,SE);
            
            for frameI = 1:length(framesToCompute)
                frameI
                TiffObj.setDirectory(framesToCompute(frameI));
                tempIm = TiffObj.read();
                obj.initialData(:,:,frameI) = tempIm;
                
                
                
                obj.RBBgSubData(:,:,frameI) = tempIm - imFilt;
                
                if ~isempty(savePath)
                    imwrite(obj.RBBgSubData(:,:,frameI), savePath, 'tiff', 'WriteMode', 'append');
                end
                
%                 obj.RBBgSubData(:,:,frameN) = imtophat(tempIm,SE);
            end   

            
            clear TiffObj
        
        end
        
        function bgROISub(obj)
            obj.bgSubSig = zeros(obj.nROI,obj.frames);
                      
            for ROIN = 1:obj.nROI    
                sig = obj.RBBgSubSig(ROIN,:);
                obj.bgSubSig(ROIN,:) = sig - obj.bgVal;
            end
        end

        
        
        
        function bgSubSig = bgSub(obj,bgsType,numFrames)
            % background subtraction with 1.avg 2.median 3.avg over
            % numFrames 4 weighted avg over numFrames
            %             obj.backgroundSubtracted = zeros(obj.width,obj.height, obj.frames);
            %             initialMatrix(obj);
            bgSubSig = zeros(size(sig,1), size(sig,2));
            
            switch bgsType
                case 'bgROI'
                    % Subtract the average value calculated from the bgROI
                    for i = 1:obj.frames
                        obj.frame = i;
                        avgBg = mean(mean(obj.image));
                        bgSubSig(:,i) = sig(:,i) - avgBg;
                    end
                
                case 'average'
                    %                     avgPix = mean(obj.initialData,3);
                    for i = 1:obj.frames
                        obj.frame = i;
                        avgBg = mean(mean(obj.image(obj.bgROI{1,1})));
                        bgSubSig(:,i) = sig(:,i) - avgBg;
                    end
                    
                case 'median'
                    medianPix = median(obj.initialData,3);
                    for i = 1:obj.frames
                        obj.backgroundSubtracted(:,:,i) = obj.initialData(:,:,i) - medianPix;
                    end
                    
                case 'avgOverFrames'
                    for i = 1:obj.frames
                        if (i <= numFrames)
                            obj.backgroundSubtracted(:,:,i) = sum(obj.initialData(:,:,i:numFrames+i-1),3)/numFrames;
                        elseif ((i + numFrames) >= obj.frames)
                            obj.backgroundSubtracted(:,:,i) = sum(obj.initialData(:,:,i-numFrames:i),3)/numFrames;
                        else
                            obj.backgroundSubtracted(:,:,i) = sum(obj.initialData(:,:,(i - numFrames):(i + numFrames)),3)/(2*numFrames+1);
                        end
                    end
                    
                    %                 case 'weightedAvgFrames'
                    %                     w = ones(numFrames-1,1)*15/(100*(numFrames-1));
                    %                     weights = [w;0.2;0.3;0.2;w];
                    %
                    %                     for i = 1:obj.frames
                    %                         if (i <= numFrames)
                    %                             weightedInitialData = weights * InitialData(:,:,
                    %                             obj.backgroundSubtracted(:,:,i) = sum(obj.InitialData(:,:,i:numFrames+i-1),3)/numFrames;
                    %                         elseif ((i + numFrames) >= obj.frames)
                    %                             obj.backgroundSubtracted(:,:,i) = sum(obj.InitialData(:,:,i-numFrames:i),3)/numFrames;
                    %                         else
                    %                             obj.backgroundSubtracted(:,:,i) = sum(obj.InitialData(:,:,(i - numFrames):(i + numFrames)),3)/(2*numFrames+1);
                    %                         end
                    %                     end
                    %
            end
            
        end
        
        function thresholdBgs(obj,thresholdType)
            initialMatrix(obj);
            
            switch thresholdType
                case 'avg'
                    meanThreshold = mean(mean(mean(obj.initialData,3),2),1) + mean(mean(std(obj.initialData,3),2),1);
                    
                    for i = 1:obj.frames
                        for j = 1:obj.width
                            for k = 1:obj.height
                                if (obj.initialData(k,j,i) < meanThreshold)
                                    obj.initialData(k,j,i) = 0;
                                else
                                end
                            end
                        end
                    end
                    
                case 'mode'
                    modeThreshold = mode(mean(mean(obj.initialData,3),2),1) + mean(mean(std(obj.initialData,3),2),1);
                    
                    for i = 1:obj.frames
                        for j = 1:obj.width
                            for k = 1:obj.height
                                if (obj.initialData(k,j,i) < modeThreshold)
                                    obj.initialData(k,j,i) = 0;
                                else
                                end
                            end
                        end
                    end
                    
            end
        end
        
        function normRawSig(obj)
            obj.sigBl = zeros(size(obj.rawSig,1),1);
            if obj.cellPixelPercent > 0.3
                obj.bgBl = median(obj.bgVal);
            else
                obj.bgBl = zeros(size(obj.rawSig,1),1);
            end
            
            for ROIN = 1:size(obj.rawSig,1)
                obj.sigBl(ROIN) = median(obj.rawSig(ROIN,:));
                nSig = obj.rawSig(ROIN,:) - obj.sigBl(ROIN);
                
                if obj.cellPixelPercent > 0.3
                    nBg = obj.bgVal - obj.bgBl;
                else
                    obj.bgBl(ROIN) = median(obj.bgVal(ROIN,obj.sigROI{3,ROIN}(1),:));
                    nBg = obj.bgVal(ROIN,obj.sigROI{3,ROIN}(1) - obj.bgBl(ROIN));
                end
                
                obj.normSig(ROIN,:) = (nSig - nBg)./obj.sigBl(ROIN);
                %                 maxInt = max(nSig);
                %                 findMax
                
            end
        end
        
        function CEDistance = imageJElectrode(obj, imageJROIPath)
            % Get the ROI for stimulating electrode
            [electrodeROI] = ReadImageJROI(imageJROIPath);
            BW = zeros(obj.height, obj.width);
            %             obj.electrodeROI{2,1} = electrodeROI{1, 1}.mnCoordinates;
            circleBound = electrodeROI.vnRectBounds;
            %                 BW = createMask(sROI{1, i}.mnCoordinates);
            % vnRectBounds(1) = height low value
            % vnRectBounds(2) = width low value
            % vnRectBounds(3) = height high value
            % vnRectBounds(4) = width high value
            circleD = abs((circleBound(3) - circleBound(1))/2);
            
            hCentre = round((circleBound(3) + circleBound(1))/2);
            wCentre = round((circleBound(2) + circleBound(4))/2);
            
            BW = mat2gray(obj.circ(BW, hCentre, wCentre, circleD));
            
            
            BW = imdilate(BW,strel('square',3)); %# dilation
            BW = imfill(BW,'holes');             %# fill inside silhouette
            BW = imerode(BW,strel('square',3));  %# erode
            
            %                 BW = imfill(BWTemp,'hole');
            obj.electrodeROI{1,1} = BW;
            pos = regionprops(BW,'ConvexHull');
            stat = regionprops(BW,'centroid');
            
            obj.electrodeROI{2,1} = round(pos.ConvexHull);
            
            obj.electrodeROI{3,1} = [round(stat.Centroid(2)),round(stat.Centroid(1))];
            %             obj.getIntMROI;
            
            CEDistance = obj.CEDistance;
            
        end
        
        function CEDistance = get.CEDistance(obj)
           % calculated electrode to cell distance
           if isempty(obj.electrodeROI)
               %                display('Electrode ROI has not been assinged yet');
               CEDistance = [];
           else
               
               CEDistance = zeros(1, obj.nROI);
               elecCentre = obj.electrodeROI{3,1};
               % 1st element is the y (height), second is the x (width)
               
               for i = 1:obj.nROI
                   
                   if ~isempty(obj.sigROI)&~isempty(obj.sigROI{3,i})

                   cellCentre = obj.sigROI{3,i};
                   CEDistance(1,i) = sqrt((elecCentre(1) + obj.electrodeYShift - cellCentre(1))^2 + ...
                       (elecCentre(2) + obj.electrodeXShift - cellCentre(2))^2);
                   end
               end
               
           end
        end
        
        
        function imageJROI(obj, imageJROIPath, varargin)
            % Update the ROI from the ROI sets drawn in ImageJ
            % must specify type as either sig (signal) or bg (background)
            % to update sigROI or bgROI respectively
            
            [sROI] = ReadImageJROI(imageJROIPath);
            
%             type = 'sig';
            bgNum = 0;
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'type')
                    type = varargin{i + 1};
                end
                
                if strcmp(varargin{i}, 'bgNum')
                    % the ROI number for background ROI
                    bgNum = varargin{i + 1};
                end
            end
            
            
            
            if strcmp(type, 'sig')
        %% Signal ROIs
                obj.sigROI = cell(3,length(sROI) - 1);
                
                sROICount = 1;
                
                for i = 1:length(sROI)
                    if i == bgNum;
                        continue;
                    end
                    BWTemp = zeros(obj.height, obj.width);
                    BWTemp = mat2gray(BWTemp);
                    obj.sigROI{2,sROICount} = sROI{1, i}.mnCoordinates;
                    %                 BW = createMask(sROI{1, i}.mnCoordinates);
                    
                    for j = 1:size(sROI{1, i}.mnCoordinates,1)
                        if sROI{1, i}.mnCoordinates(j,1) == 0
                            sROI{1, i}.mnCoordinates(j,1) = 1;
                        end
                        
                        if sROI{1, i}.mnCoordinates(j,2) == 0
                            sROI{1, i}.mnCoordinates(j,2) = 1;
                        end
                        BWTemp(sROI{1, i}.mnCoordinates(j,2),sROI{1, i}.mnCoordinates(j,1)) = 1;
                        
                        
                    end
                    
                    BW = imdilate(BWTemp,strel('square',10)); %# dilation
                    BW = imfill(BW,'holes');             %# fill inside silhouette
                    BW = imerode(BW,strel('square',10));  %# erode
                    
                    %                 BW = imfill(BWTemp,'hole');
                    obj.sigROI{1,sROICount} = BW;
                    stat = regionprops(BW,'centroid');
                    
                    obj.sigROI{3,sROICount} = [round(stat.Centroid(2)),round(stat.Centroid(1))];
                    
                    sROICount = sROICount + 1;
                end
            elseif strcmp(type, 'bg')
        %% Background ROI (default is the first entry of the ROI set
                if length(sROI) == 1
                    bgROI = sROI;
                else
                    bgROI = sROI{1, bgNum};
                end
                
                obj.bgROI = cell(3,length(bgROI));
                
                
                
                BWTemp = zeros(obj.height, obj.width);
                BWTemp = mat2gray(BWTemp);
                obj.bgROI{2,1} = bgROI.mnCoordinates;
                %                 BW = createMask(sROI{1, i}.mnCoordinates);
                
                for j = 1:size(bgROI.mnCoordinates,1)
                    if bgROI.mnCoordinates(j,1) == 0
                        bgROI.mnCoordinates(j,1) = 1;
                    end
                    
                    if bgROI.mnCoordinates(j,2) == 0
                        bgROI.mnCoordinates(j,2) = 1;
                    end
                    BWTemp(bgROI.mnCoordinates(j,2),bgROI.mnCoordinates(j,1)) = 1;
                    
                    
                end
                
                BW = imdilate(BWTemp,strel('square',10)); %# dilation
                BW = imfill(BW,'holes');             %# fill inside silhouette
                BW = imerode(BW,strel('square',10));  %# erode
                
                %                 BW = imfill(BWTemp,'hole');
                obj.bgROI{1,1} = BW;
                stat = regionprops(BW,'centroid');
                
                obj.bgROI{3,1} = [round(stat.Centroid(2)),round(stat.Centroid(1))];
                
            end
            
            
            
            %             obj.getIntMROI;
        end
        
        
            
        
        function mROI(obj)
            [x,y,b] = ginput(1);
            if b == 'c'
                obj.sigROI = {};
                obj.bgROI = {};
                obj.nROI = 0;
            end
            %             obj.plotROI('var')
            obj.plotROI;
            while b ~= 'e'
                
                if b == 'c'
                    obj.sigROI = {};
                    obj.bgROI = {};
                end
                
                
                h = imfreehand;
                %                 pos = h.getPosition;
                %                 [pos(:,1), pos(:,2)] = checkBound(pos(:,1),pos(:,2),obj.height, obj.width);
                BW = createMask(h);
                stat = regionprops(BW,'centroid');
                pos = regionprops(BW,'ConvexHull');
                
                if isempty(obj.bgROI)
                    obj.bgROI{1,end + 1} = BW;
                    obj.bgROI{2,end} = round(pos.ConvexHull);
                    obj.bgROI{3,end} = [round(stat.Centroid(2)),round(stat.Centroid(1))];
                else
                    obj.sigROI{1,end + 1} = BW;
                    obj.sigROI{2,end} = round(pos.ConvexHull);
                    obj.sigROI{3,end} = [round(stat.Centroid(2)),round(stat.Centroid(1))];
                end
                
%                 obj.plotROI('max');
%                 obj.plotROI('var');
                obj.plotROI;
                delete(h);
                [x,y,b] = ginput(1);
                obj.nROI = size(obj.sigROI,2);
                
            end
            
            emptyEntry = cellfun(@isempty,obj.sigROI);
            sumEmp = sum(emptyEntry,1);
            emptyCol = find(sumEmp > 0);
            emptyEntry(:,emptyCol) = 1;
            obj.sigROI(emptyEntry) = [];
            obj.sigROI = reshape(obj.sigROI,3,[]);
            % This is to get rid of any empty entries when sometimes the
            % imfreehand was not executed properly by the user
            
            obj.nROI = size(obj.sigROI,2);
            obj.getIntMROI;
            
        end
        
        
        function sortROI(obj)
            
            sigROITemp = [obj.sigROI, obj.tempSigROI];
            sigROICt = zeros(obj.nROI, 2);
            
            
            for i = 1:length(sigROITemp)
                sigROICt(i,:) = sigROITemp{3,i};
            end
            
            [temp, sortI] = sortrows(sigROICt, [2,1]);
            
            obj.organiseROI(sigROITemp, sortI);
        end
        
        function remapROI(obj, remapROIMat)
            %             nROI = obj.nROI;
%             assignedSigROI = [];
            
            sigROITemp = obj.sigROI;
            rawSigTemp = obj.rawSig;
            
            RBBgSubSigTemp = obj.RBBgSubSig;
            %             if ~exist('sortI')
            %                 sortI = 1:size(sigROITemp,2);
            %             end
            obj.sigROI = cell(3,size(remapROIMat,2));
            obj.rawSig = zeros(size(remapROIMat,2),size(rawSigTemp,2));
            
            obj.RBBgSubSig = zeros(size(remapROIMat,2),size(RBBgSubSigTemp,2));
            
            
            obj.sigROI = cell(3, max((remapROIMat(2,:))));
            %% 20150923 version
            for i = 1:size(remapROIMat,2)
                if ~isnan(remapROIMat(2,i))
                    if i ==65
                        i
                    end
                    obj.sigROI{1,i} = sigROITemp{1,remapROIMat(2,i)};
                    obj.sigROI{2,i} = sigROITemp{2,remapROIMat(2,i)};
                    obj.sigROI{3,i} = sigROITemp{3,remapROIMat(2,i)};
                    
%                     assignedSigROI = [assignSigROI, remapROIMat(2,i)];
                    if ~isempty(obj.rawSig)
                        obj.rawSig(i,:) = rawSigTemp(remapROIMat(2,i),:);
                    end
                    
                    if ~isempty(obj.RBBgSubSig)
                        obj.RBBgSubSig(i,:) = RBBgSubSigTemp(remapROIMat(2,i),:);
                    end
                end
            end
            
            allROI = 1:size(sigROITemp,2);
            missingROIs = setdiff(allROI,remapROIMat(2,:));
            
            for i = 1:length(missingROIs)
                if isempty(obj.sigROI(missingROIs(i)))
                    obj.sigROI{1,missingROIs(i)} = sigROITemp{1,missingROIs(i)};
                    obj.sigROI{2,missingROIs(i)} = sigROITemp{2,missingROIs(i)};
                    obj.sigROI{3,missingROIs(i)} = sigROITemp{3,missingROIs(i)};
                    
                    if ~isempty(obj.rawSig)
%                         obj.rawSig(missingROIs(i),:) = rawSigTemp(remapROIMat(2,missingROIs(i)),:);
                         obj.rawSig(missingROIs(i),:) = rawSigTemp(missingROIs(i),:);
                    end
                    
                    if ~isempty(obj.RBBgSubSig)
                        obj.RBBgSubSig(missingROIs(i),:) = RBBgSubSigTemp(missingROIs(i),:);
                    end
                else
                    obj.sigROI{1,end + 1} = sigROITemp{1,missingROIs(i)};
                    obj.sigROI{2,end} = sigROITemp{2,missingROIs(i)};
                    obj.sigROI{3,end} = sigROITemp{3,missingROIs(i)};
                    
                    if ~isempty(obj.rawSig)
                        obj.rawSig(end + 1,:) = rawSigTemp(missingROIs(i),:);
                    end
                    
                    if ~isempty(obj.RBBgSubSig)
                        obj.RBBgSubSig(end + 1,:) = RBBgSubSigTemp(missingROIs(i),:);
                    end
                end
            end
            
%             for i = 1:length(sigROITemp)
%                 if ~isnan(remapROIMat(2,i))
%                     obj.sigROI{1,remapROIMat(2,i)} = sigROITemp{1,i};
%                     obj.sigROI{2,remapROIMat(2,i)} = sigROITemp{2,i};
%                     obj.sigROI{3,remapROIMat(2,i)} = sigROITemp{3,i};
%                     
%                     
%                     obj.rawSig(remapROIMat(2,i),:) = rawSigTemp(i,:);
%                     
%                     obj.RBBgSubSig(remapROIMat(2,i),:) = RBBgSubSigTemp(i,:);
%                 end
%             end
%             obj.organiseROI(sigROITemp,remapROIMat(2,:));
        end
        
        function organiseROI(obj, sigROITemp, sortI)
            
            rawSigTemp = obj.rawSig;
            
            RBBgSubSigTemp = obj.RBBgSubSig;
%             if ~exist('sortI')
%                 sortI = 1:size(sigROITemp,2);
%             end
            
            for i = 1:length(sigROITemp)
                
                %                 if nargin > 1
                %                     testSkip = sum(skipROI == i);
                %                 else
                %                     testSkip = 0;
                %                 end
                
                %                 if testSkip == 0
                obj.sigROI{1,i} = sigROITemp{1,sortI(i)};
                obj.sigROI{2,i} = sigROITemp{2,sortI(i)};
                obj.sigROI{3,i} = sigROITemp{3,sortI(i)};
                
                if ~isempty(obj.rawSig)                
                    obj.rawSig(i,:) = rawSigTemp(sortI(i),:);
                end
                
                if ~isempty(obj.RBBgSubSig)
                    obj.RBBgSubSig(i,:) = RBBgSubSigTemp(sortI(i),:);
                end
                %                 else
                %                     obj.sigROI{1,i} = [];
                %                     obj.sigROI{2,i} = [];
                %                     obj.sigROI{3,i} = [];
                %
                %                     obj.tempSigROI{1,end + 1} = sigROITemp{1,sortI(i)};
                %                     obj.tempSigROI{2,i} = sigROITemp{2,sortI(i)};
                %                     obj.tempSigROI{3,i} = sigROITemp{3,sortI(i)};

            end
        end

        
        function reMROI(obj,ROIN)
            obj.plotROI;
            h = imfreehand;
            
            BW = createMask(h);
                stat = regionprops(BW,'centroid');
                pos = regionprops(BW,'ConvexHull');
            
            obj.sigROI{1,ROIN} = BW;
            obj.sigROI{2,ROIN} = round(pos.ConvexHull);
            obj.sigROI{3,ROIN} = [round(stat.Centroid(2)),round(stat.Centroid(1))];
            
            obj.plotROI;
        end
        
        function storeGROI(obj, ROIN)
            for i = 1:length(ROIN)
                obj.gROI(end + 1) = ROIN(i);
            end
        end
        
        function storeLSROI(obj,ROIN)
            for i = 1:length(ROIN)
                obj.lsROI(end + 1) = ROIN(i);
            end
        end
        
        function saveMROI(obj,saveName)
            [pathstr,~,~] = fileparts(obj.path);
            sigROI = obj.sigROI;
            bgROI = obj.bgROI;
            gROI = obj.gROI;
            lsROI = obj.lsROI;
            
            % mROIByConv parameters
            alpha = obj.alpha;
            minGrowth = obj.minGrowth;
            threshF = obj.threshF;
            
            
            currTime = clock;
            currTime = strcat(num2str(currTime(4)),num2str(currTime(5)));
            refNum = obj.refROINum;
            currDate = date;
            
            if nargin > 1
                saveString = sprintf('%s', saveName);
            else
                saveString = sprintf('%s %s  refROINum %d',...
                    currDate, currTime, refNum);
            end
            
            if (obj.refROINum ~= 0)&isnumeric(obj.refROINum)
                vertMove = obj.refVertMov;
                horiMove = obj.refHoriMov;
                refROINum = obj.refROINum;
                refROIDate = obj.refROIdate;
                refROITime = obj.refROITime;
                save([pathstr,saveString], 'sigROI', 'bgROI','gROI', 'lsROI', 'alpha', 'minGrowth', 'threshF',...
                    'vertMove','horiMove','refROINum', 'refROIDate','refROITime');
            else
                save([pathstr,saveString], 'sigROI', 'bgROI','gROI', 'lsROI','alpha', 'minGrowth', 'threshF');
            end
            
            
        end
        
        function moveROI(obj, horiMove, vertMove)
            
            
            for i = 1:size(obj.sigROI,2)
                obj.sigROI{2,i}(:,1) = obj.sigROI{2,i}(:,1) + vertMove;

                obj.sigROI{2,i}(:,2) = obj.sigROI{2,i}(:,2) + horiMove;
                
                BWTemp = zeros(obj.height, obj.width);
                BWTemp = mat2gray(BWTemp);

                for j = 1:size(obj.sigROI{2,i},1)
                    if obj.sigROI{2,i}(j,1) < 1
                        obj.sigROI{2,i}(j,1) = 1;
                    elseif obj.sigROI{2,i}(j,1) > obj.width
                        obj.sigROI{2,i}(j,1) = obj.width;
                    end
                    
                    if obj.sigROI{2,i}(j,2) < 1
                        obj.sigROI{2,i}(j,2) = 1;
                    elseif obj.sigROI{2,i}(j,2) > obj.height
                        obj.sigROI{2,i}(j,2) = obj.height;
                    end
                    BWTemp(obj.sigROI{2,i}(j,2),obj.sigROI{2,i}(j,1)) = 1;
                    
                end
                BW = imdilate(BWTemp,strel('square',3)); %# dilation
                BW = imfill(BW,'holes');             %# fill inside silhouette
                BW = imerode(BW,strel('square',3));  %# erode
    
              
                obj.sigROI{1,i} = BW;
                stat = regionprops(BW,'centroid');
            
                obj.sigROI{3,i} = [round(stat.Centroid(2)),round(stat.Centroid(1))];

                
                
            end
            
%             obj.getIntMROI;
        end
        
        function saveRefROI(obj,num)
            % Need to number the ref ROI because there may be more than 1
            % imaging field.
            [pathstr,~,~] = fileparts(obj.path);
            parts = strsplit(pathstr, '\');
            parentDir = [];
            for i = 1:(length(parts) - 1)
                if i == 1
                    parentDir = [parentDir, parts{i}];
                else
                    parentDir = [parentDir,'\' parts{i}];
                end
            end
            refROI = obj.sigROI;
            pClampNo = obj.pClampNo;
            imageNo = obj.imageNo;
            
            if nargin < 2
                findROI = sprintf('\\*.mat');
                ROIStruct = dir([parentDir,findROI]);
                num = length(ROIStruct) - 1;
                % If num is not updated (incremented) by the user, then the
                % default will assume it is updating the existing reference
                % ROI and will save a copy of the updated RefROI with
                % different date and time but the same num.
            end
            
            currTime = clock;
            currTime = strcat(num2str(currTime(4)),num2str(currTime(5)));
            
            saveString = sprintf('/%s %s referenceROINum %d', date, currTime, num);
            save([parentDir,saveString], 'refROI','pClampNo','imageNo');
        end
        
        function loadRefROI(obj, num)
            % Need to number the ref ROI because there may be more than 1
            % imaging field.
            [pathstr,~,~] = fileparts(obj.path);
            parts = strsplit(pathstr, '\');
            parentDir = [];
            for i = 1:(length(parts) - 1)
                if i == 1
                    parentDir = [parentDir, parts{i}];
                else
                    parentDir = [parentDir,'\' parts{i}];
                end
            end
            
            if nargin > 1
                findROI = sprintf('\\*%d.mat', num);
            else
                findROI = sprintf('\\*.mat');
            end
            ROIStruct = dir([parentDir,findROI]);
            
            S = [ROIStruct(:).datenum].';
            [~,SI] = sort(S);
            ROIStruct = ROIStruct(SI);
            ROIPath = fullfile(parentDir,ROIStruct(end).name);
            
            refROIName = ROIStruct(end).name;
            refROINameSplit = strsplit(refROIName, ' ');
            tempMat = load(ROIPath);
            
            obj.refROI = tempMat.refROI;
            obj.refROIDate = refROINameSplit{1};
            obj.refROITime = refROINameSplit{2};
            obj.refROINum = refROINameSplit{end}(1);
            
        end
        
        
        function loadMROI(obj,ROIName, noInitDat)
            [pathstr,name,ext] = fileparts(obj.path);
            
            if nargin < 2
                findROI = sprintf('\\*.mat');
                ROIStruct = dir([pathstr,findROI]);
                S = [ROIStruct(:).datenum].';
                [~,SI] = sort(S);
                ROIStruct = ROIStruct(SI);
                ROIPath = fullfile(pathstr,ROIStruct(end).name);
            else
                ROIPath = fullfile(pathstr,ROIName);
                
                
            end
            

            
            tempMat = load(ROIPath);
            
            obj.sigROI = tempMat.sigROI;
            obj.bgROI = tempMat.bgROI;
            
            if isfield(tempMat,'gROI') == 1
                obj.gROI = tempMat.gROI;
            end
            
            if isfield(tempMat, 'lsROI') == 1
                obj.lsROI = tempMat.lsROI;
            end
            
            if isfield(tempMat, 'horiMove') == 1
                obj.refHoriMove = tempMat.horiMove;
                obj.refVertMove = tempMat.vertMove;
            end
            
            if isfield(tempMat, 'refROINum') == 1
                obj.refROINum = tempMat.refROINum;
            end
            
            if isfield(tempMat, 'refROIDate') == 1
                obj.refROINum = tempMat.refROINum;
            end
            
            
            obj.nROI = size(obj.sigROI,2);
            if (nargin < 3)|(noInitDat == 0)
                obj.getIntMROI;
                % do not load initialData matrix if the user just wants to
                % have a look at the object or for other purposes that do
                % not require initialData matrix.
            end
        end
        
        function mRefROI(obj, num, refSigROI)
            % This function is used when the user has a reference ROIs (in
            % sigROI format) and would like to apply these on the Ca image with
            % ROI movement. This is called when several recordings of the
            % same imaging field (similar cells) are to be analysed and
            % give the benefit of the same ROI number for the same cell for
            % comparison and analysis.
            
            % Method: write a reference ROI template and then record the
            % movement required to map the ROI template for each recording.
            % Later on when additional ROIs are added to a recording, they
            % can be conveniently added to other recordings as well by
            % looking at the movement values.
            CaObj.sigROI = {};
            
            vertMov = 0;
            horiMov = 0;
            
            BWTemp = zeros(obj.height,obj.width);
            
            if nargin == 3
                
            elseif nargin < 2
                obj.loadRefROI;
                refSigROI = obj.refROI;
            elseif nargin < 3
                obj.loadRefROI(num)
                refSigROI = obj.refROI;
                obj.refROINum = num;
            end
            
            for i = 1:length(refSigROI)
                BWTemp = BWTemp + refSigROI{1,i};
            end
            
            [hI, wI] = find(BWTemp ~= 0);
            
            
            [x,y,b] = ginput(1);
            
            if isa(obj.image,'uint16')
                logVal = log10(double(obj.image));
                maxIntensity = max(max(logVal));
                minVal = maxIntensity - 1.6;
                photo = (logVal - minVal )/(maxIntensity - minVal);
            end
            
            imagesc(photo);
            %             obj.plotROI
            while b ~= 'e'
                if b == 'd'
                    break
                end
                [x,y,b] = ginput(1);
                
                if b == 30
                    % 30 => upArrow
                    hI = hI - 1;
                    vertMov = vertMov - 1;
                elseif b == 31
                    % 31 => downArrow
                    hI = hI + 1;
                    vertMov = vertMov + 1;
                elseif b == 28;
                    % 28 => leftArrow
                    wI = wI - 1;
                    horiMov = horiMov - 1;
                elseif b == 29;
                    % 29 => rightArrow
                    wI = wI + 1;
                    horiMov = horiMov + 1;
                end
                
                exceedBound = find((hI < 1)|(wI < 1)|(hI > obj.height)|(wI > obj.width));
                hITemp = hI;
                wITemp = wI;
                hITemp(exceedBound) = [];
                wITemp(exceedBound) = [];
                
                if isa(obj.image,'uint16')
                    logVal = log10(double(obj.image));
                    maxIntensity = max(max(logVal));
                    minVal = maxIntensity - 1.6;
                    photo = (logVal - minVal )/(maxIntensity - minVal);
                end
                
                pbordered = cat(3,photo,photo,photo);
                
                BW = zeros(obj.height,obj.width);
                for i = 1:length(hITemp)
                    BW(hITemp(i),wITemp(i)) = 1;
                    pbordered(hITemp(i),wITemp(i),1) = 1;
                    pbordered(hITemp(i),wITemp(i),2) = 0;
                    pbordered(hITemp(i),wITemp(i),3) = 0;
                end
                
                
                imshow(pbordered);
            end
            
            if b == 'd'
                % cancel
            else
                
                obj.refVertMov = vertMov;
                obj.refHoriMov = horiMov;
                
                for i = 1:length(refSigROI)
                    
                    BWTemp = refSigROI{1,i};
                    [hI, wI] = find(BWTemp ~= 0);
                    hI = hI + vertMov;
                    wI = wI + horiMov;
                    exceedBound = find((hI < 1)|(wI < 1)|(hI > obj.height)|(wI > obj.width));
                    hITemp = hI;
                    wITemp = wI;
                    hITemp(exceedBound) = [];
                    wITemp(exceedBound) = [];
                    
                    BW = zeros(obj.height, obj.width);
                    for j = 1:1:length(hITemp)
                        BW(hITemp(j),wITemp(j)) = 1;
                    end
                    
                    stat = regionprops(BW,'centroid');
                    pos = regionprops(BW,'ConvexHull');
                    
                    if ~isempty(hITemp)
                        obj.sigROI{1,i} = BW;
                        obj.sigROI{2,i} = round(pos.ConvexHull);
                        obj.sigROI{3,i} = [round(stat.Centroid(2)),round(stat.Centroid(1))];
                    else
                        obj.sigROI{1,i} = BW;
                        obj.sigROI{2,i} = [];
                        obj.sigROI{3,i} = [];
                    end
                    
                end
                obj.nROI = size(obj.sigROI,2);
            end
        end
        
        
        function dates = querryMROI(obj)
            [pathstr,name,ext] = fileparts(obj.path);
            
            dates = {};
            
            findROI = sprintf('\\*.mat');
            ROIStruct = dir([pathstr,findROI]);
            
            % sort the ROIStruct
            S = [ROIStruct(:).datenum].';
            [~,SI] = sort(S);
            ROIStruct = ROIStruct(SI);
            
            for i = 1:length(ROIStruct)
                dates{end + 1} = ROIStruct(i).date;
            end
            
        end
        
        function mROIByVarM(obj, alpha, threshF, minGrowth)
            % Use variance matrix to find the ROIs (It is not very good in
            % terms of getting good shape of cells becaus sometimes it is
            % only part of the cells responding.
            
            % alpha is a coefficient for threshold intensity for finding
            % background value
            %             if nargin < 4
            %                 alpha = 1;
            %             end
            warning('off','all');
            %             distance = CalDistance(convR);
            % Generating distance matrix for the making of coupling strength matrix
            
            if isempty(obj.initialData)
                obj.initialMatrix;
            end
            %             obj.convR = convR;
            %             circItem = circConvItem(convR);
            %             rectItem = rectConvItem(1,convR);
            
            
            varMat = zeros(obj.width, obj.height);
            
            
            for row = 1:obj.height
                for col = 1:obj.width
                    varMat(col,row) = var(double(obj.initialData(col, row,:)));
                end
            end
            
            obj.sigROI = {};
            obj.bgROI = {};
            
            % 1s row is the BW mask (0 for all entries and 1 for circled
            % region). 2nd row is the circled points' position (rounded to
            % integer). 3rd row is the centre position (rounded to integer). If
            % the ROI method is convolution, then 2nd and 3rd values will be
            % the same.
            
            meanVarMat = mean(mean(varMat));
            thresh = alpha * meanVarMat;
            % randomly select a region tha has lo variance as the
            % background ROI
            bgVarMat = varMat < thresh;
            
            while true
                bgBW = zeros(obj.height, obj.width);
                rRow = round((obj.width - (0.2*obj.width))*rand(1,1)) + round(0.2*obj.width);
                rCol = round((obj.height - (0.2*obj.height))*rand(1,1)) + round(0.2*obj.height);
                % The scalar (obj.width - round(0.1*obj.width)) and
                % + round(0.1*obj.width)is to make sure that the
                % background ROI is not going to be too close to the field
                % edges.
                
                rRowAll = (rRow - round(0.02*obj.width)):(rRow + round(0.02*obj.width));
                rColAll = (rCol - round(0.02*obj.height)):(rCol + round(0.02*obj.height));
                bgBW(rRowAll,rColAll) = 1;
                if sum(sum(bgBW(~bgVarMat))) == 0
                    break
                end
            end
            
            
            bgBW = logical(bgBW);
            
            pos = regionprops(bgBW,'ConvexHull');
            
            obj.bgROI{1,end + 1} = bgBW;
            obj.bgROI{2,end} = round(pos.ConvexHull);
            obj.bgROI{3,end} = [rRow,rCol];
            
            
            
            while true
                
                [maxTest, maxI] = max(varMat(:));
                thresh
                maxTest
                if maxTest < thresh
                    break
                end
                [roiCCol, roiCRow] = ind2sub(size(varMat),maxI);
                opts.threshold = maxTest*threshF;
                %                 opts.threshold = threshF;
                opts.minGrowth = minGrowth;
                
                minThresh = 2*thresh;
                
                maskTemp = obj.fillROI(varMat, roiCCol, roiCRow ,opts, minThresh);
                pos = regionprops(maskTemp,'ConvexHull');
                obj.sigROI{1, end + 1} = maskTemp;
                nROI = size(obj.sigROI,2);
                obj.sigROI{2, nROI} = round(pos.ConvexHull);
                obj.sigROI{3, nROI}= [roiCCol, roiCRow];
                varMat(maskTemp) = 0;
                
            end
            
            obj.nROI = nROI;
            obj.alpha = alpha;
            obj.threshF = threshF;
            obj.minGrowth = minGrowth;
            
            %             obj.getIntMROI;
            
            %             padCircLayer = [size(circItem,1),size(circItem,2)];
            %             padImCirc = padarray(obj.initialData(:,:,refFrame), [padCircLayer(1), padCircLayer(2)],'circular');
            %             convImCircTemp = conv2(padImCirc,circItem,'same');
            %
            %             padRectLayer = [size(rectItem,1),size(rectItem,2)];
            %             padImRect = padarray(obj.initialData(:,:,refFrame), [padRectLayer(1), padRectLayer(2)],'circular');
            %             convImRectTemp = conv2(padImRect,rectItem,'same');
            %
            %             convImCirc = convImCircTemp(padCircLayer(1) + 1:end - padCircLayer(1), padCircLayer(2) + 1:end - padCircLayer(2));
            %             convImRect = convImRectTemp(padRectLayer(1) + 1:end - padRectLayer(1), padRectLayer(2) + 1:end - padRectLayer(2));
            %
            %
            %
            %             if nargin < 3
            %                 thresh = mean(mean(convImCirc));
            %             else
            %                 thresh = alpha*mean(mean(convImCirc));
            %             end
            
            
        end
        
        function mROIByInt(obj, alpha, threshF, minGrowth)
            % Use image intensity to draw ROIs.
            
            % alpha is a coefficient for threshold intensity for finding
            % background value
            %             if nargin < 4
            %                 alpha = 1;
            %             end
            warning('off','all');
            %             distance = CalDistance(convR);
            % Generating distance matrix for the making of coupling strength matrix
            
            if isempty(obj.initialData)
                obj.initialMatrix;
            end
            %             obj.convR = convR;
            %             circItem = circConvItem(convR);
            %             rectItem = rectConvItem(1,convR);
            
            
            
            % 1s row is the BW mask (0 for all entries and 1 for circled
            % region). 2nd row is the circled points' position (rounded to
            % integer). 3rd row is the centre position (rounded to integer). If
            % the ROI method is convolution, then 2nd and 3rd values will be
            % the same.
            
            imMat = zeros(obj.width, obj.height);
            %             imMatSqrt = zeros(obj.width, obj.height);
            
            for row = 1:obj.height
                for col = 1:obj.width
                    imMat(col,row) = max(double(obj.initialData(col, row,:)));
                    imMatLog(col,row) = max(log10(double(obj.initialData(col, row,:))));
                end
            end
            
            
            meanInt = mean(mean(imMat));
            thresh = alpha * meanInt;
            % randomly select a region tha has lo variance as the
            % background ROI
            bgImMat = meanInt < thresh;
            
            while true
                bgBW = zeros(obj.height, obj.width);
                rRow = round((obj.width - (0.2*obj.width))*rand(1,1)) + round(0.2*obj.width);
                rCol = round((obj.height - (0.2*obj.height))*rand(1,1)) + round(0.2*obj.height);
                % The scalar (obj.width - round(0.1*obj.width)) and
                % + round(0.1*obj.width)is to make sure that the
                % background ROI is not going to be too close to the field
                % edges.
                
                rRowAll = (rRow - round(0.02*obj.width)):(rRow + round(0.02*obj.width));
                rColAll = (rCol - round(0.02*obj.height)):(rCol + round(0.02*obj.height));
                bgBW(rRowAll,rColAll) = 1;
                if sum(sum(bgBW(~bgImMat))) == 0
                    break
                end
            end
            
            
            bgBW = logical(bgBW);
            
            pos = regionprops(bgBW,'ConvexHull');
            
            obj.bgROI{1,end + 1} = bgBW;
            obj.bgROI{2,end} = round(pos.ConvexHull);
            obj.bgROI{3,end} = [rRow,rCol];
            
            obj.alpha = alpha;
            obj.threshF = threshF;
            obj.minGrowth = minGrowth;
            
            while true
                
                [maxTest, maxI] = max(imMat(:));
                thresh
                maxTest
                if maxTest < thresh
                    break
                end
                [roiCCol, roiCRow] = ind2sub(size(imMat),maxI);
                opts.threshold = maxTest*threshF;
                %                 opts.threshold = threshF;
                opts.minGrowth = minGrowth;%Hello, john. Hello CH or KJ?; Walk home safe
                
                minThresh = meanInt;
                % the minimum threshold at which the fill algorithm should
                % stop
                maskTemp = obj.fillROI(imMat, roiCCol, roiCRow ,opts, minThresh);
                pos = regionprops(maskTemp,'ConvexHull');
                obj.sigROI{1, end + 1} = maskTemp;
                nROI = size(obj.sigROI,2);
                obj.sigROI{2, nROI} = round(pos.ConvexHull);
                obj.sigROI{3, nROI}= [roiCCol, roiCRow];
                imMat(maskTemp) = 0;
                obj.nROI = nROI;
            end
            
            obj.nROI = nROI;
            
        end
        
        function mROIByLog(obj, alpha, threshF, minGrowth)
            % Use image intensity to draw ROIs.
            
            % alpha is a coefficient for threshold intensity for finding
            % background value
            %             if nargin < 4
            %                 alpha = 1;
            %             end
            warning('off','all');
            %             distance = CalDistance(convR);
            % Generating distance matrix for the making of coupling strength matrix
            
            if isempty(obj.initialData)
                obj.initialMatrix;
            end
            %             obj.convR = convR;
            %             circItem = circConvItem(convR);
            %             rectItem = rectConvItem(1,convR);
            
            
            
            % 1s row is the BW mask (0 for all entries and 1 for circled
            % region). 2nd row is the circled points' position (rounded to
            % integer). 3rd row is the centre position (rounded to integer). If
            % the ROI method is convolution, then 2nd and 3rd values will be
            % the same.
            
            imMatLog = zeros(obj.width, obj.height);
            %             imMatSqrt = zeros(obj.width, obj.height);
            
            for row = 1:obj.height
                for col = 1:obj.width
                    imMatLog(col,row) = max(log10(double(obj.initialData(col, row,:))));
                end
            end
            
            obj.sigROI = {};
            obj.bgROI = {};
            
            
            meanInt = mean(mean(imMatLog));
            thresh = alpha * meanInt;
            % randomly select a region tha has lo variance as the
            % background ROI
            bgImMat = meanInt < thresh;
            
            while true
                bgBW = zeros(obj.height, obj.width);
                rRow = round((obj.width - (0.2*obj.width))*rand(1,1)) + round(0.2*obj.width);
                rCol = round((obj.height - (0.2*obj.height))*rand(1,1)) + round(0.2*obj.height);
                % The scalar (obj.width - round(0.1*obj.width)) and
                % + round(0.1*obj.width)is to make sure that the
                % background ROI is not going to be too close to the field
                % edges.
                
                rRowAll = (rRow - round(0.02*obj.width)):(rRow + round(0.02*obj.width));
                rColAll = (rCol - round(0.02*obj.height)):(rCol + round(0.02*obj.height));
                bgBW(rRowAll,rColAll) = 1;
                if sum(sum(bgBW(~bgImMat))) == 0
                    break
                end
            end
            
            
            bgBW = logical(bgBW);
            
            pos = regionprops(bgBW,'ConvexHull');
            
            obj.bgROI{1,end + 1} = bgBW;
            obj.bgROI{2,end} = round(pos.ConvexHull);
            obj.bgROI{3,end} = [rRow,rCol];
            
            obj.alpha = alpha;
            obj.threshF = threshF;
            obj.minGrowth = minGrowth;
            
            while true
                
                [maxTest, maxI] = max(imMatLog(:));
                thresh
                maxTest
                if maxTest < thresh
                    break
                end
                [roiCCol, roiCRow] = ind2sub(size(imMatLog),maxI);
                opts.threshold = maxTest*threshF;
                %                 opts.threshold = threshF;
                opts.minGrowth = minGrowth;%Hello, john. Hello CH or KJ?; Walk home safe
                
                minThresh = meanInt;
                % the minimum threshold at which the fill algorithm should
                % stop
                maskTemp = obj.fillROI(imMatLog, roiCCol, roiCRow ,opts, minThresh);
                pos = regionprops(maskTemp,'ConvexHull');
                obj.sigROI{1, end + 1} = maskTemp;
                nROI = size(obj.sigROI,2);
                obj.sigROI{2, nROI} = round(pos.ConvexHull);
                obj.sigROI{3, nROI}= [roiCCol, roiCRow];
                imMatLog(maskTemp) = 0;
                obj.nROI = nROI;
            end
            
            %             obj.nROI = nROI;
            
        end
        
        function getStimCoupled(obj)
            coupledMatrix = zeros(obj.nROI, length(obj.patch.STSStartPt));
            
            for i = 1:obj.nROI
                rawSigTemp = obj.rawBlCorrSig(i,:);
                for j = 1:length(obj.patch.STSStartPt)
                    STSStartFrame = round(obj.patch.STSStartPt(j)/obj.patch.SR*obj.fps);
                    STSEndFrame = round(obj.patch.STSEndPt(j)/obj.patch.SR*obj.fps);
                    
                    stimFrameTemp =...
                        find((obj.patch.extStimStartPt > obj.patch.STSStartPt(j))&(obj.patch.extStimStartPt < obj.patch.STSEndPt(j)));
                    stimStartFrame = floor(obj.patch.extStimStartPt(stimFrameTemp(1))/obj.patch.SR*obj.fps);
                    % The first trigger of extStimStartPt is to start
                    % camera imaging, not for indcation of stim time
                    
                    rawSegSig = rawSigTemp(STSStartFrame:STSEndFrame);
                    
                    bsF = median(rawSegSig);
                    DF = rawSegSig - bsF;
                    DFFSigTemp = DF/bsF*100;

                    % Turn into DFF signal at each segment for coupled
                    % repsonse check
                    
                    
                    DFFBaseSig = DFFSigTemp(1:(stimStartFrame - STSStartFrame + 1));
                    blSTD = std(DFFBaseSig);
                    
                    DFFStimSig = DFFSigTemp((stimStartFrame - STSStartFrame + 2):(STSEndFrame - STSStartFrame + 1));
                    sigSTD = std(DFFStimSig);
                    
                    smDFFStimSig = smooth(DFFStimSig, 5);
                    if ~isempty(smDFFStimSig > 3*blSTD)
                        [extStimFramePeak, extStimFrameLocs] =...
                            findpeaks(DFFStimSig, 'MINPEAKHEIGHT', 3*blSTD,'MINPEAKDISTANCE',round(5));
                    end
                    stimLatency = extStimFrameLocs/obj.fps; 
                    % get the latency of the calcium signals
                    
%                     if stimLatency < 0.2 % If stimulation latency is less than 100ms, the calcium signal is considered as a coupled response
                    if ~isempty(stimLatency)
                        coupledMatrix(i,j) = stimLatency(1);
                        if i == 55
                        i
                        end
                        j
                        obj.fitSig(i,j, obj.frameTime(STSStartFrame:STSEndFrame))
                    end
%                     end
                    
                    
                    
                end
            end
            
            obj.coupledResponse = coupledMatrix;
        end
        
        
        function dROI(obj, ROIN)
            % delete ROI
            ROIN = sort(ROIN,'descend');
            for i = 1:length(ROIN)
                obj.tempSigROI(:,end + 1) = obj.sigROI(:,ROIN(i));
                obj.sigROI(:,ROIN(i)) = [];
                obj.nROI = obj.nROI - 1;
            end
        end
        
        function aROI(obj, ROIInfo,ROIN)
            obj.sigROI(:,(ROIN + 1:end + 1)) = obj.sigROI(:,(ROIN:end));
            obj.sigROI(:,ROIN) = ROIInfo;
            obj.nROI = length(obj.sigROI);
        end
        
        function mask = fillROI(obj, data, maxR, maxC ,opts, minThresh)
            % The funciton uses the input maximum row and maximum column
            % as the centre point of the cell and then expand from this point
            % to fill the ROI. Developed by Calvin Eiber
            
            % parse inputs
            
            % keep track of the difference between original point and newly
            % found points. If the differences became too great, then it
            % means it is out of the ROI range.
            x = round(maxR); y = round(maxC);
            
            % default configuration
            if ~exist('opts','var'), opts = struct(); end
            %             if isfield(opts, 'threshold'), T = opts.threshold; else T = 100; end
            
            if isfield(opts,'threshold')
                T = opts.threshold;
            else
                T = 100;
            end
            if isfield(opts, 'minGrowth')
                if isempty(opts.minGrowth)
                    opts.minGrowth = 0.1;
                end
            end
            
            mask = false(size(data));
            mask(x,y) = true;
            
            TRecord = []; % threshold record
            TRecord(1) = T;
            
            while true
                
                old_mask = mask; % expand mask by 1 voxel in each direction
                mask(1:end-2,:) = mask(1:end-2,:) | old_mask(2:end-1,:);
                mask(3:end, :) = mask(3:end, :) | old_mask(2:end-1,:);
                mask(:,1:end-2) = mask(:,1:end-2) | old_mask(:,2:end-1);
                mask(:, 3:end) = mask(:, 3:end) | old_mask(:,2:end-1);
                %                 mask(:,:,1:end-2) = mask(:,:,1:end-2) | old_mask(:,:,2:end-1);
                %                 mask(:, :, 3:end) = mask(:, :, 3:end) | old_mask(:,:,2:end-1);
                
                new_mask = mask & ~old_mask;
                %                 threshVal = mean(data(old_mask))*T;
                keepers = data(new_mask) > T; % keep only sufficiently bright voxels
                newMaskVal = data(new_mask);
                mask(new_mask) = keepers;
                
                % If we are updating threshold on each iteration,  do it.
                %                 if isfield(opts, 'thresholdFcn'),
                % update threshold
                %                     if length(keepers) < 30
                %                     T = obj.threshF*mean(newMaskVal(keepers));
                T = obj.threshF*mean(data(mask));
                TRecord(end + 1) = T;
                diffTR = diff(TRecord);
                meanDiffTR = mean(diffTR);
                %                     firstDiffTR = TRecord(1) - TRecord(2);
                %                     lastTDiff = TRecord(end - 1) - TRecord(end);
                
                if T > 1.5*minThresh
                    diffTRFac = 0.7;
                else
                    diffTRFac = 0.95;
                end
                
                if length(TRecord) > 2
                    if TRecord(end) < diffTRFac*meanDiffTR
                        break
                    end
                end
                %                 end
                % If we kept 10% [opts.MinGrowth] or fewer of the points on the
                % perimeter, terminate ROI growth
                sum(keepers)
                if T < minThresh
                    break
                end
                %                 opts.minGrowth * sum(new_mask(:))
                if (sum(keepers) <= opts.minGrowth * sum(new_mask(:)))
                    break
                end
                if (sum(keepers) == 0)
                    break
                end
                if (opts.minGrowth * sum(new_mask(:)) == 0)
                    break
                end
            end
        end
        
        
        
        function getIntMROI(obj,varargin)
            %             if isempty(obj.initialData)
            %                 obj.initialMatrix;
            %             end
            
            startFrame = 1;
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'startFrame')
                    startFrame = varargin{i + 1};
                end
            end
            
            analyseBlockN = 500;
            % Need to break down the number of frames to be analysed
            % otherwise memory will be full
            
            obj.rawSig = zeros(size(obj.sigROI,2), obj.frames);
            obj.RBBgSubSig = zeros(size(obj.sigROI,2), obj.frames);
            
%             endFrame = startFrame;
            
            while startFrame < obj.frames
                obj.RBBgSubData = [];
                obj.initialData = [];
                
                if (startFrame + analyseBlockN - 1) < obj.frames
                    framesToCompute = startFrame:(startFrame + analyseBlockN - 1);
                else
                    framesToCompute = startFrame:obj.frames;
                end
                
                frameN = length(framesToCompute);
                
            
%             if isempty(obj.RBBgSubData)
                obj.rollingBallBgSub('framesToCompute',framesToCompute);
%             end
            
            
            
            rawTempImage = reshape(obj.initialData,[obj.height, obj.width*frameN]);
            rawTempImage = double(rawTempImage);
            
            
            RBBgSubTempImage = reshape(obj.RBBgSubData,[obj.height, obj.width*frameN]);
            RBBgSubTempImage = double(RBBgSubTempImage);
            
            

            
            for ROIN = 1:size(obj.sigROI,2)
                ROIN;
                sigROIMat = repmat(obj.sigROI{1,ROIN}, 1,frameN);
                
%                 rawSigTemp = tempImage.*sigROIMat;

                pixelNum = sum(sum(logical(obj.sigROI{1,ROIN})));
                
                %% RBBgSub signal
                
                RBBgSubSigTemp = RBBgSubTempImage(logical(sigROIMat));
                
                RBBgSubSegSigTemp = reshape(RBBgSubSigTemp, [pixelNum, frameN]);
                
                obj.RBBgSubSig(ROIN,framesToCompute) = sum(RBBgSubSegSigTemp,1)/pixelNum;
                
                
                %% Raw signal
                
                rawSigTemp = rawTempImage(logical(sigROIMat));
                
                rawSegSigTemp = reshape(rawSigTemp, [pixelNum, frameN]);
                
                obj.rawSig(ROIN,framesToCompute) = sum(rawSegSigTemp,1)/pixelNum;

                
                clear sigROIMat pixelNum rawSigTemp rawSigSegTemp
                
            end
            
            
            for frameI = 1:length(framesToCompute)
                tempBgROIVal = obj.RBBgSubData(:,:,frameI);
                
                tempBgROIVal = tempBgROIVal((logical(obj.bgROI{1})));
                
                obj.bgVal(framesToCompute) = mean(tempBgROIVal);
                
            end
            
%             obj.rawSig(ROIN,frameN) = sum(sum(tempImage(logical(obj.sigROI{1,ROIN}))))/sum(sum(logical(obj.sigROI{1,ROIN})));
            
            
            %             for ROIN = 1:size(obj.sigROI,2)
            %                 for frameN = 1:obj.frames
            %
            % %                     tempImage = obj.initialData(:,:,frameN);
            %    obj.rawSig(ROIN,frameN) = sum(sum(tempImage(logical(obj.sigROI{1,ROIN}))))/sum(sum(logical(obj.sigROI{1,ROIN})));
            %                 end
            %             end
            %             obj.getDFFSig;
            
            startFrame = startFrame + analyseBlockN;
            end
        end
        
        
        function getIntByConv(obj, convR,alpha)
            % alpha is a coefficient for threshold intensity for finding
            % background value
            warning('off','all');
            %             distance = CalDistance(convR);
            % Generating distance matrix for the making of coupling strength matrix
            
            obj.convR = convR;
            circItem = circConvItem(convR);
            rectItem = rectConvItem(1,convR);
            
            obj.sigROI = {};
            obj.rawSig = [];
            
            obj.bgROI = {};
            obj.bgVal = [];
            
            
            for i = 1:obj.frames
                i
                %                 obj.frame = i;
                padCircLayer = [size(circItem,1),size(circItem,2)];
                padImCirc = padarray(obj.initialData(:,:,i), [padCircLayer(1), padCircLayer(2)],'circular');
                convImCircTemp = conv2(padImCirc,circItem,'same');
                
                padRectLayer = [size(rectItem,1),size(rectItem,2)];
                padImRect = padarray(obj.initialData(:,:,i), [padRectLayer(1), padRectLayer(2)],'circular');
                convImRectTemp = conv2(padImRect,rectItem,'same');
                
                convImCirc = convImCircTemp(padCircLayer(1) + 1:end - padCircLayer(1), padCircLayer(2) + 1:end - padCircLayer(2));
                convImRect = convImRectTemp(padRectLayer(1) + 1:end - padRectLayer(1), padRectLayer(2) + 1:end - padRectLayer(2));
                
                if i == 1
                    % This is for finding the local maximum, which
                    % represents the average intensities of a particular
                    % cell after convolution. The local maximum found are
                    % further tested with their intensities to make sure
                    % they are clear cells for analysis
                    
                    obj.frame = i;
                    
                    if nargin < 3
                        thresh = mean(mean(convImCirc));
                    else
                        thresh = alpha*mean(mean(convImCirc));
                    end
                    
                    BW = imregionalmax(convImCirc);
                    %                     obj.localMaxF = BW;
                    
                    BW(convImCirc < thresh) = 0;
                    obj.sigROI{1,1} = BW;
                    [localMX, localMY] = find(BW == 1);
                    
                    nCell = length(localMX);
                    
                    for ROIN = 1:nCell
                        obj.sigROI{3,ROIN}(1) = localMX(ROIN);
                        obj.sigROI{3,ROIN}(2) = localMY(ROIN);
                    end
                    
                    
                    cellAreaCheck = (convImCirc > thresh);
                    cellAreaPixelNumber = sum(sum(cellAreaCheck));
                    
                    obj.cellPixelPercent = cellAreaPixelNumber/(obj.height*obj.width);
                    
                    %                     lineMin = min(convImRect,[],2);
                    %                     lineMax = max(convImRect,[],2);
                    lineMax = max(convImCirc,[],2);
                    if (obj.cellPixelPercent > 0.3)&(obj.width < 30)
                        
                        diffLineMax = diff(lineMax);
                        edge = find(diffLineMax > 0.8*(max(diffLineMax)));
                        lowIntLineMax = find(lineMax < 0.2*(max(lineMax)));
                        
                        lowIntAndEdge = ismember(lowIntLineMax, edge);
                        if any(lowIntAndEdge);
                            findTemp = find(lowIntAndEdge == 1);
                            bgPtTemp = lowIntLineMax(findTemp(end));
                            obj.bgROI{3,1}(1) = bgPtTemp(end);
                        else
                            lowIntLineMaxTemp = lowIntLineMax;
                            lowIntLineMaxTemp(lowIntLineMaxTemp > round(obj.height/2)) = obj.height - lowIntLineMaxTemp(lowIntLineMaxTemp > round(obj.height/2));
                            
                            
                            % find the line that is farthest away from the
                            % edge
                            [~,maxI] = max(lowIntLineMaxTemp);
                            obj.bgROI{3,1} = lowIntLineMax(maxI);
                        end
                    end
                end
                
                %                 [lineMin,lineMinI] = min(convImRect,[],2);
                %                 lineMax = max(convImRect,[],2);
                
                if (obj.cellPixelPercent > 0.3)&(obj.width < 30)
                    %                     lineMin = min(convImRect,[],2);
                    [lineMin, minPos] = min(convImCirc,[],2);
                    obj.bgVal(i) = lineMin(obj.bgROI{3,1}(1));
                    obj.bgROI{3,1}(end + 1) = minPos(obj.bgROI{3,1}(1));
                else
                    %                     [lineMin,lineMinI] = min(convImRect,[],2);
                    [lineMin,lineMinI] = min(convImCirc,[],2);
                    
                    
                    
                    for ROIN = 1:nCell
                        obj.bgROI{3,ROIN} = [obj.sigROI{3,ROIN}(1),lineMinI(obj.sigROI{3,ROIN}(1))];
                        obj.bgVal(ROIN,i) = lineMin(obj.sigROI{3,ROIN}(1));
                    end
                    
                end
                %
                %                 obj.frame = i;
                %                 obj.bgVal(i) =  sum(sum(obj.image(obj.bgBW)))/sum(sum(obj.bgBW));
                
                for ROIN = 1:size(obj.sigROI, 2)
                    obj.convIm(:,:,i) = convImCirc;
                    obj.rawSig(ROIN,i) = convImCirc(obj.sigROI{3,ROIN}(1),obj.sigROI{3,ROIN}(2));
                    % convInt values, not corrected by background
                    % subtraction
                    %                     if cellPixelPercent > 0.3
                    %                         obj.convIntSubBg(j,i) = obj.convInt(j,i) - obj.bgVal(i);
                    %                     else
                    %                         obj.convIntSubBg(j,i) = obj.convInt(j,i) - obj.bgLineVal(obj.convPos(1,j),i);
                    %                     end
                    
                    obj.sigROI{1,1}(obj.sigROI{3,ROIN}(1),obj.sigROI{3,ROIN}(2)) = ROIN;
                end
            end
        end
        
        
        function plotROI(obj, plotMethod, ROIN)
            
            % plot manual ROIs on the image
            if nargin > 2
                figure
                ROIForPlot = ROIN;
            else
                gcf
                clf
                ROIForPlot = 1:obj.nROI;
            end
            
            if isa(obj.image,'uint16')
                if (nargin < 2)|strcmp(plotMethod, 'org')
                    logVal = log10(double(obj.image));
                    
                    maxIntensity = max(max(logVal));
                    minVal = maxIntensity - 1.6;
                    photo = (logVal - minVal )/(maxIntensity - minVal);
                elseif nargin > 1
                    if strcmp(plotMethod, 'max')
                        zProject = zeros(obj.height,obj.width);
                        for i = 1:obj.frames
                            obj.frame = 1;
                            zProject = max(zProject,double(obj.image));
                        end
                        
                        logVal = log10(zProject);
%                         photo = logVal;
                    elseif strcmp(plotMethod, 'var')
                        
                        if isempty(obj.varMat)
                            varMat = zeros(obj.height, obj.width);
                            if isempty(obj.initialData)
                                obj.initialMatrix;
                            end
                            
                            for row = 1:obj.height
                                for col = 1:obj.width
                                    varMat(row,col) = var(double(obj.initialData(row, col,:)));
                                end
                            end
                            obj.varMat = varMat;
                        end
                        
                        logVal = log10(double(obj.varMat));
                        maxIntensity = max(max(logVal));
                        minVal = maxIntensity - 4.5;
                        photo = (logVal - minVal )/(maxIntensity - minVal);
                        
                        
                    elseif strcmp(plotMethod, 'diff')
                        
                        
                        diffMat = zeros(obj.height, obj.width);
                        if isempty(obj.initialData)
                            obj.initialMatrix;
                        end
                        
                        for row = 1:obj.height
                            for col = 1:obj.width
                                tempDiff = diff(double(obj.initialData(row, col,:)));
                                diffMat(row,col) = sum(tempDiff.^2);
                                
                            end
                        end

                        
                        logVal = log10(double(diffMat));
                        maxIntensity = max(max(logVal));
                        minVal = maxIntensity - 4.5;
                        photo = (logVal - minVal )/(maxIntensity - minVal);
                    end
                    
                end
                
                pbordered = cat(3,photo,photo,photo);
                [x,y] = size(photo);
                
                % background point
                if ~isempty(obj.bgROI)
                    if size(obj.bgROI, 1) > 3
                        pbordered(obj.bgROI{3,1}(1),obj.bgROI{3,1}(2),1) = 0; % draw in red
                        pbordered(obj.bgROI{3,1}(1),obj.bgROI{3,1}(2),2) = 0;
                        pbordered(obj.bgROI{3,1}(1),obj.bgROI{3,1}(2),3) = 1;
                    end
                    
                    if size(obj.sigROI, 1) > 0
                        if ~isempty(obj.sigROI{2,1})
                            index = (obj.bgROI{2,1}(:,1) - 1)*obj.height + obj.bgROI{2,1}(:,2);
                            pbordered(index) = rand(1);
                            pbordered(index + obj.width*obj.height) = rand(1);
                            if index < obj.width*obj.height
                                pbordered(index + 2*obj.width*obj.height) = rand(1);
                            end
                        end
                    end
                end
                
                for ROIN = 1:size(ROIForPlot, 2)
%                     ROIN
                    if ~isempty(obj.sigROI{3,ROIForPlot(ROIN)})
                        pbordered(obj.sigROI{3,ROIForPlot(ROIN)}(1),obj.sigROI{3,ROIForPlot(ROIN)}(2),1) = 1; % draw in red
                        pbordered(obj.sigROI{3,ROIForPlot(ROIN)}(1),obj.sigROI{3,ROIForPlot(ROIN)}(2),2) = 0;
                        pbordered(obj.sigROI{3,ROIForPlot(ROIN)}(1),obj.sigROI{3,ROIForPlot(ROIN)}(2),3) = 0;
                    end
                    if ~isempty(obj.sigROI{2,ROIForPlot(ROIN)})
                        index = (obj.sigROI{2,ROIForPlot(ROIN)}(:,1) - 1)*obj.height + obj.sigROI{2,ROIForPlot(ROIN)}(:,2);
                        pbordered(index) = rand(1);
                        pbordered(index + obj.width*obj.height) = rand(1);
                        if index < obj.width*obj.height
                            pbordered(index + 2*obj.width*obj.height) = rand(1);
                        end
                    end
                    
                end
                
                imshow(pbordered);
                if nargin < 2
                    for ROIN = 1:size(ROIForPlot,2)
                        if ~isempty(obj.sigROI{3,ROIForPlot(ROIN)})
                            text(obj.sigROI{3,ROIForPlot(ROIN)}(2),obj.sigROI{3,ROIForPlot(ROIN)}(1),num2str(ROIForPlot(ROIN)),'color',[1,0.6,0.1],'FontSize', 20, 'FontWeight', 'bold');
                            %                       text(obj.sigROI{3,ROIN}(2),obj.sigROI{3,ROIN}(1),num2str(ROIN),'color','c','FontSize', 20);
                        end
                    end
                end
                
                %             axis xy
                
                set(gca, 'visible', 'on', 'XGrid', 'on', 'YGrid', 'on');
                set(gcf,'position',get(0,'screensize'))
                
            end
            
            axis equal 
            axis tight
        end
        
        
        function getIntensities(obj,mask)
            
            if nargin > 1
                obj.mask = mask;
            end
            
            nROI = max(max(abs(obj.mask)));
            %             obj.sumInt = zeros(nROI,obj.frames);
            obj.rawSig = zeros(nROI, obj.frames);
            %This one-way pass is slightly faster than the n-way find above...
            elements = zeros(1,nROI); %How many pixels are there in the image for various neurons.
            for frameN = 1:obj.frames
                obj.frame = frameN;
                frameN
                for ROIN = 1:nROI
                    tempROIMask = (obj.mask == -ROIN);
                    if frameN == 1
                        elements(ROIN) = sum(sum(tempROIMask));
                    end
                    obj.rawSig(ROIN,frameN) = sum(sum((obj.image(tempROIMask))));
                end
                
            end
            
            
            for ROIN = 1:nROI
                if (elements(ROIN) > 0)
                    obj.rawSig(ROIN,:) = obj.rawSig(ROIN,:)./elements(ROIN);
                end
            end
        end
       
        
        function getStimFrameSort(obj, ROIN, offSetTime, blTime)
            % offSetTime and blTime are user defined time after and before
            % stimulation for offset signal and baseline signal
            % calculations respectively.
            % If they are not given, default values are given
            
            obj.stimFrame = {};
            
            bsFrameForOverlap = 5;
            % If baseline frame overlaps with previous stim offset time,
            % then just get 5 frames before stimulation as baseline frame.
            
            if (nargin < 3)|isempty(offSetTime)|(offSetTime == 0)
                offSetTime = 5;
            end
            % 5 second as the offset time after the end of
            % stimulation. The user can change this value here.
            
            if (nargin < 4)|isempty(blTime)|(blTime == 0)
                blTime = 1;
            end
            % 0.2 second as the baseline time before the start of
            % stimulation.
            
            
            if isempty(obj.stimOnFrameCell)
                getStimFrame(obj,ROIN)
            end
            
            if ~isempty(obj.patch.lineATime)
                obj.cellImageTime = obj.frameTime(1:obj.frames) - (obj.patch.lineATime)*(obj.height - obj.sigROI{3,ROIN}(1));
            else
                obj.cellImageTime = obj.frameTime(1:obj.frames);
            end
            
            if obj.patch.ndim == 3
                dataTimeTemp = obj.patch.gapFreeTime;
            else
                dataTimeTemp = obj.patch.dataTime;
            end
            
            for stimN = 1:obj.nStim
                
                obj.stimFrame{2,stimN}(1,:) = obj.stimOnFrameCell{1,stimN};
                obj.stimSegFrameON{stimN} = obj.stimOnFrameCell{1,stimN};
                obj.stimFrame{2,stimN}(2,:) = obj.stimOnFrameCell{2,stimN};
                obj.stimSegFrameTimeON{stimN} = obj.stimOnFrameCell{2,stimN};
                
                obj.stimFrame{5,stimN}(1,:) = obj.stimOnFrameCell{3,stimN};
                obj.stimFrame{5,stimN}(2,:) = obj.stimOnFrameCell{4,stimN};
                % Stimulation (extracellular stimulation) onset frames
                
                
                currentStimStart = obj.stimOnFrameCell{2,stimN}(1);
                blTimeCheck = currentStimStart - blTime;
                % baseline time check to see if it overlaps with
                % previous stim's offset time
                
                if stimN == 1
                    if blTimeCheck < obj.frameTime(1)
                        
                        blTemp = obj.cellImageTime < currentStimStart;
                        obj.stimFrame{1, stimN}(1,:) = find(blTemp);
                        obj.stimFrame{1, stimN}(2,:) = obj.cellImageTime(blTemp);
                        obj.stimFrame{4, stimN}(1,:) = 1:(obj.stimOnFrameCell{3,stimN}(1) - 1);
                        obj.stimFrame{4, stimN}(2,:) = dataTimeTemp(obj.stimFrame{4,stimN}(1,:));
                    else
                        blPts = ceil(blTime*obj.patch.SR);
                        blTemp = (obj.cellImageTime > blTimeCheck)&(obj.cellImageTime < currentStimStart);
                        obj.stimFrame{1, stimN}(1,:) = find(blTemp);
                        obj.stimFrame{1, stimN}(2,:) = obj.cellImageTime(blTemp);
                        obj.stimFrame{4, stimN}(1,:) = (obj.stimOnFrameCell{3,stimN}(1) - 1 - blPts):(obj.stimOnFrameCell{3,stimN}(1) - 1);
                        obj.stimFrame{4, stimN}(2,:) = dataTimeTemp(obj.stimFrame{4, stimN}(1,:));
                    end
                else
                    if blTimeCheck > offEndTime
                        %                             % if baseline time does not overlap with offset
                        %                             % time
                        blPts = ceil(blTime*obj.patch.SR);
                        % Convert baseline time into baseline data
                        % points by timing the baseline time with
                        % Clampex sampling rate
                        blTemp = (obj.cellImageTime > blTimeCheck)&(obj.cellImageTime < currentStimStart);
                        obj.stimFrame{1, stimN}(1,:) = find(blTemp);
                        obj.stimFrame{1, stimN}(2,:) = obj.cellImageTime(blTemp);
                        %                         stimN
                        obj.stimFrame{4, stimN}(1,:) = (obj.stimOnFrameCell{3,stimN}(1) - 1 - blPts):(obj.stimOnFrameCell{3,stimN}(1) - 1);
                        obj.stimFrame{4, stimN}(2,:) = dataTimeTemp(obj.stimFrame{4, stimN}(1,:));
                    else
                        blTimeTemp = find(obj.cellImageTime < currentStimStart);
                        blTimeTemp = blTimeTemp(end);
                        blPts = ceil(bsFrameForOverlap/obj.fps*obj.patch.SR);
                        % Converting baseline frames into baseline data
                        % points by dividing number of frames by image
                        % sampling rate first (to get the baseline time
                        % when overlapping) then times this number with
                        % Clampex sampling rate
                        obj.stimFrame{1, stimN}(1,:) = (blTimeTemp - bsFrameForOverlap):1:blTimeTemp;
                        obj.stimFrame{1, stimN}(2,:) = obj.cellImageTime(obj.stimFrame{1, stimN}(1,:));
                        obj.stimFrame{4, stimN}(1,:) = (obj.stimOnFrameCell{3,stimN}(1) - 1 - blPts):(obj.stimOnFrameCell{3,stimN}(1) - 1);
                        obj.stimFrame{4, stimN}(2,:) = dataTimeTemp(obj.stimFrame{4, stimN}(1,:));
                    end
                end
                
                obj.stimSegFrameBs{stimN} = obj.stimFrame{1, stimN}(1,:);
                if isempty(obj.stimSegFrameBs{stimN})
                    obj.stimSegFrameBs{stimN} = 1;
                end
                
                obj.stimSegFrameTimeBs{stimN} = obj.stimFrame{1, stimN}(2,:);
                
                
                if stimN ~= obj.nStim
                    % Test if another stimulation occurs before the 5s
                    % decay time. If so, set the offset time between the
                    % end of current stimlation onset and the start of next
                    % stimulation onset
                    currentStimEnd = obj.stimOnFrameCell{2,stimN}(end);
                    nextStimStart = obj.stimOnFrameCell{2,stimN + 1}(1);
                else
                    currentStimEnd = obj.stimOnFrameCell{2,stimN}(end);
                    nextStimStart = obj.frameTime(end);
                end
                
                
                offTimeCheck = nextStimStart - currentStimEnd;
                offEndTime = currentStimEnd + offSetTime;
                % the time at which the offset time ends
                
                if offTimeCheck < offSetTime
                    offSetTemp = (obj.cellImageTime > currentStimEnd)&(obj.cellImageTime < nextStimStart);
                    offSetPts = offSetTime*obj.patch.SR;
                    obj.stimFrame{3,stimN}(1,:) = find(offSetTemp);
                    obj.stimFrame{3,stimN}(2,:) = obj.cellImageTime(offSetTemp);
                    
                    if stimN ~= obj.nStim
                        obj.stimFrame{6,stimN}(1,:) = (obj.stimOnFrameCell{3,stimN}(end) + 1):(obj.stimOnFrameCell{3,stimN + 1}(1) - 1);
                    else
                        obj.stimFrame{6,stimN}(1,:) = (obj.stimOnFrameCell{3,stimN}(end) + 1):(obj.patch.nSamples*obj.patch.nSweeps);
                    end
                    
                    obj.stimFrame{6,stimN}(2,:) = dataTimeTemp(obj.stimFrame{6, stimN}(1,:));
                else
                    offSetPts = offSetTime*obj.patch.SR;
                    offSetFrame = find((obj.cellImageTime > currentStimEnd)&(obj.cellImageTime < offEndTime));
                    obj.stimFrame{3,stimN}(1,:) = offSetFrame;
                    obj.stimFrame{3,stimN}(2,:) = obj.cellImageTime(offSetFrame);
                    
                    endDataPoint = (obj.stimOnFrameCell{3,stimN}(end) + offSetPts);
                    
                    if endDataPoint < (obj.patch.nSamples*obj.patch.nSweeps)
                        obj.stimFrame{6,stimN}(1,:) = (obj.stimOnFrameCell{3,stimN}(end) + 1):(obj.stimOnFrameCell{3,stimN}(end) + offSetPts);
                    else
                        obj.stimFrame{6,stimN}(1,:) = (obj.stimOnFrameCell{3,stimN}(end) + 1):(obj.patch.nSamples*obj.patch.nSweeps);
                    end
                    obj.stimFrame{6,stimN}(2,:) = dataTimeTemp(obj.stimFrame{6, stimN}(1,:));
                end
                
                obj.stimSegFrameOFF{stimN} = obj.stimFrame{3,stimN}(1,:);
                obj.stimSegFrameTimeOFF{stimN} = obj.stimFrame{3,stimN}(2,:);
                
                obj.stimFrame{7,stimN}(1,:) = [obj.stimFrame{1, stimN}(1,:), obj.stimFrame{2, stimN}(1,:), obj.stimFrame{3, stimN}(1,:)];
                obj.stimFrame{7,stimN}(2,:) = [obj.stimFrame{1, stimN}(2,:), obj.stimFrame{2, stimN}(2,:), obj.stimFrame{3, stimN}(2,:)];
                obj.stimFrame{8,stimN}(1,:) = [obj.stimFrame{4, stimN}(1,:), obj.stimFrame{5, stimN}(1,:), obj.stimFrame{6, stimN}(1,:)];
                obj.stimFrame{8,stimN}(2,:) = [obj.stimFrame{4, stimN}(2,:), obj.stimFrame{5, stimN}(2,:), obj.stimFrame{6, stimN}(2,:)];
                
                obj.stimSegFrameWhole{stimN} = obj.stimFrame{7,stimN}(1,:);
                obj.stimSegFrameTimeWhole{stimN} = obj.stimFrame{7,stimN}(2,:);
            end
        end
        
        
        function rawBlCorrSig = get.rawBlCorrSig(obj) % photobleach correction
%             frame = 1:size(obj.rawSig,2);

            if isempty(obj.bgSubSig)
               obj.bgROISub; 
            end
            frame = 1:size(obj.bgSubSig,2);

            %             for ROIIndex = 1:size(obj.rawSig,1)
            initialVal = obj.bgSubSig(obj.cROIN, 1);
            finalVal = obj.bgSubSig(obj.cROIN, end);
            diffFluo = abs(finalVal - initialVal);
            %                     checkPhotoBleach
            fluo = obj.bgSubSig(obj.cROIN, :);
            if diffFluo > 10
                
                rawBlCorrSig = detrend(fluo) + fluo(1);
                
%                 smFluo = smooth(fluo,20);
%                 p = polyfit(frame(1:end/2),smFluo(1:end/2)',1);
                
                %                     yFit =  p(1)*frame + p(2);
                
                %                     photobleachCorrSig = fluo - yFit;
%                 rawBlCorrSig = fluo - p(1)*frame;
            else
                rawBlCorrSig = fluo;
            end
            
            %             end
        end
        
                function sortRGCLightType(obj,varargin)
            % This function is to determine whether a RGC is an ON, OFF or
            % ON-OFF RGC or did not respond to light stimulation. This
            % function gets the light ON and OFF segment and then
            % determines by comparing the baseline and peaks to determine
            % whether they respond to light ON or OFF.
            
            obj.ONInc = zeros(1,obj.nROI);
            obj.ONDec = zeros(1,obj.nROI);
            obj.OFFInc = zeros(1,obj.nROI);
            obj.OFFDec = zeros(1,obj.nROI);
            
            obj.RGCType = cell(1, obj.nROI);
            
            allROI = 1:obj.nROI;

            pTestVal = 0.025;
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'ROIN')
                    % Check plotType
                    allROI = varargin{i + 1};
                end
                
                if strcmp(varargin{i}, 'pVal')
                    % Check plotType
                    pTestVal = varargin{i + 1};
                end
                
                
            end
                
            obj.lightPVal = pTestVal;
            
            if (strcmp(obj.stimType, 'light'))|(strcmp(obj.patch.stimType, 'light'))
                
                for i = 1:length(allROI)
                    
                    obj.cROIN = allROI(i);
                    
                    avgSegSig = obj.getAvgSegSig;
                    segSigSTD = obj.getAvgSegSTD;
                    avgONSegSig = avgSegSig{1};
                    avgOFFSegSig = avgSegSig{2};
                    
                    obj.lightONTime = length(avgONSegSig)/obj.fps;
                    obj.lightOFFTime = length(avgOFFSegSig)/obj.fps;
                    
                    avgONSegSTD = segSigSTD{1};
                    avgOFFSegSTD = segSigSTD{2};
                    
                    coupledTime = 1;
                    % The calcium signal should rise to significant value
                    % coupledTime after stimulation
                    
                    % Note that the end of OFF segment can be considered as the
                    % baseline of ON segment and the end of ON segment can be
                    % considered as the baseline of OFF segment.
                    %% Getting all the essential values for ON and OFF determination
                    
                    coupledFrame = round(coupledTime*obj.fps);
                    
                    LEDONStartVal = avgONSegSig(1);
                    LEDONStartSTD = avgONSegSTD(1);
                    
                    LEDONEndVal = avgONSegSig(end);
                    LEDONEndSTD = avgONSegSig(end);
                    
                    LEDOFFStartVal = avgOFFSegSig(1);
                    LEDOFFStartSTD = avgOFFSegSTD(1);
                    
                    LEDOFFEndVal = avgOFFSegSig(end);
                    LEDOFFEndSTD = avgOFFSegSTD(end);
                    if coupledFrame < length(avgONSegSig)
                        [LEDONPeak, LEDONPeakLoc] = max(avgONSegSig(1:coupledFrame));
                    else
                        [LEDONPeak, LEDONPeakLoc] = max(avgONSegSig);
                    end
                    LEDONPeakSTD = avgONSegSTD(LEDONPeakLoc);
                    
                    
                    if coupledFrame < length(avgONSegSig)
                        [LEDOFFPeak, LEDOFFPeakLoc] = max(avgOFFSegSig(1:coupledFrame));
                    else
                        [LEDOFFPeak, LEDOFFPeakLoc] = max(avgOFFSegSig);
                    end
                    LEDOFFPeakSTD = avgOFFSegSTD(LEDOFFPeakLoc);
                    
%                     LEDONBlVal = mean(avgOFFSegSig((end - 5):end));
%                     % ON's baseline is considered to be the end of OFF segment
%                     LEDOFFBlVal = mean(avgONSegSig((end - 5):end));
%                     % OFF's baseline is considered to be the end of ON segment
%                     
%                     LEDONBlSTD = std(avgOFFSegSig((end - 5):end));
%                     LEDOFFBlSTD = std(avgONSegSig((end - 5):end));
                    

                    if i == 25
                        i
                    end
                    
                    %% Determining ON response
                    %                     ONInitBl = LEDONStartVal - LEDONBlVal;
                    %                     ONLastInit = LEDONEndVal - LEDONStartVal;
                    
                    
                    %                     if (ONLastInit < -3)&(LEDONPeakLoc == 1)

                    
                    % On segment initial value minus baseline
                    %                     ONPeakbl = LEDONPeak - LEDONBlVal;
                    %                     ONPeakInit = LEDONPeak - LEDONStartVal;
                    
                    tValONInc = obj.tTest(LEDONPeak, LEDONStartVal, LEDONPeakSTD, LEDONStartSTD);
                    pValONInc = 1 - tcdf(tValONInc, obj.patch.nStim - 1);
                    
                    tValONDec = obj.tTest(LEDONStartVal, LEDONEndVal, LEDONStartSTD, LEDONEndSTD);
                    pValONDec = 1 - tcdf(tValONDec, obj.patch.nStim - 1);
                    
                    
%                     frameONDiff = diff(avgONSegSig(1:LEDONPeakLoc));
%                     monoIncreaseONCheck = all(frameONDiff > 0);
                    
                    %                     if (ONPeakbl > 3*LEDONBlSTD)&(ONPeakInit > 5)
                    if (pValONInc < pTestVal)&(tValONInc > 0)
                        % If the peak is a true calcium signal, then we
                        % would expect all frames between initial frame and
                        % peak frame should be monotonically increasing. So
                        % a true calcium has to also satisfy this
                        % condition.
                        
                        % Check if two frames after the peak are also
                        % significantly different so to eliminate the
                        % possibility of having a high noise peak. Only
                        % apply to ON because there is bleaching problem
                        
%                         if (LEDONPeakLoc + 2) < length(avgONSegSig)
%                             t2ValONInc = obj.tTest(avgONSegSig(LEDONPeakLoc + 2), LEDONStartVal, avgONSegSTD(LEDONPeakLoc + 2), LEDONStartSTD);
%                             p2ValONInc = 1 - tcdf(t2ValONInc, obj.patch.nStim - 1);
%                             
%                             t1ValONInc = obj.tTest(avgONSegSig(LEDONPeakLoc + 1), LEDONStartVal, avgONSegSTD(LEDONPeakLoc + 1), LEDONStartSTD);
%                             p1ValONInc = 1 - tcdf(t1ValONInc, obj.patch.nStim - 1);
                        if (LEDONPeakLoc + 1) < length(avgONSegSig)
                            
                            t1ValONInc = obj.tTest(avgONSegSig(LEDONPeakLoc + 1), LEDONStartVal, avgONSegSTD(LEDONPeakLoc + 1), LEDONStartSTD);
                            p1ValONInc = 1 - tcdf(t1ValONInc, obj.patch.nStim - 1);
                            
%                             p2ValONInc = 0;
                        elseif (LEDONPeakLoc > 2)
                            p1ValONInc = 0;
                            
                            
                        else                            
                            p1ValONInc = 1;
                            % If the frame is too short, then render
                            % detection of significant fluorescent change
                            % invalid.
                            
                        end
                        
                        if (p1ValONInc < pTestVal)
                            obj.ONInc(allROI(i)) = 1;
                        end
                        
                    end
                        
                    if (pValONDec < pTestVal)&(tValONDec > 0)&(obj.ONInc(allROI(i)) ~= 1)
                        % If there is no peak signal (no light ON
                        % response), then check if the signal is decreasing
                        % significantly, which can be used to identify OFF
                        % response if OFF stim is too short.
                        if length(avgONSegSig) > 1
                            
                            t1ValONDec = obj.tTest(avgONSegSig(2), LEDONEndVal, avgONSegSTD(2), LEDONEndSTD);
                            p1ValONDec = 1 - tcdf(t1ValONDec, obj.patch.nStim - 1);
                            
                            % Compare the second value to the last one to
                            % make sure that significant difference in
                            % first and last frame values is not due to
                            % noise.
                            
                            %                             p2ValONInc = 0;
                                                       
                        else
                            p1ValONDec = 1;
                            % If the frame is too short, then render
                            % detection of significant fluorescent change
                            % invalid.
                            
                        end
                        
                        if (p1ValONDec < pTestVal)
                            obj.ONDec(allROI(i)) = 1;
                        end
                        % If ON initial value is lower than baseline, then that
                        % means ON causes decrease in cell response
                        
                    end
                    
                    %% Determining OFF response
                    %                     OFFLastInit = LEDOFFEndVal - LEDOFFStartVal;
                    %                     OFFInitBl = LEDOFFStartVal - LEDOFFBlVal;
                    
                    tValOFFInc = obj.tTest(LEDOFFPeak, LEDOFFStartVal, LEDOFFPeakSTD, LEDOFFStartSTD);
                    pValOFFInc = 1 - tcdf(tValOFFInc, obj.patch.nStim - 1);
                    
                    tValOFFDec = obj.tTest(LEDOFFStartVal, LEDOFFEndVal, LEDOFFStartSTD, LEDOFFEndSTD);
                    pValOFFDec = 1 - tcdf(tValOFFDec, obj.patch.nStim - 1);
                    
                    %                     frameOFFDiff = diff(avgOFFSegSig(1:LEDOFFPeakLoc));
                    %                     monoIncreaseOFFCheck = all(frameOFFDiff > 0);
                    
                    %                     if (OFFLastInit < -3)&(LEDOFFPeakLoc == 1)
                    if (pValOFFInc < pTestVal)&(tValOFFInc > 0)
                        
                         
%                         if (LEDOFFPeakLoc + 2) < length(avgOFFSegSig)
%                             t2ValOFFInc = obj.tTest(avgOFFSegSig(LEDOFFPeakLoc + 2), LEDOFFStartVal, avgOFFSegSTD(LEDOFFPeakLoc + 2), LEDOFFStartSTD);
%                             p2ValOFFInc = 1 - tcdf(t2ValOFFInc, obj.patch.nStim - 1);
%                             
%                             t1ValOFFInc = obj.tTest(avgOFFSegSig(LEDOFFPeakLoc + 1), LEDOFFStartVal, avgOFFSegSTD(LEDOFFPeakLoc + 1), LEDOFFStartSTD);
%                             p1ValOFFInc = 1 - tcdf(t1ValOFFInc, obj.patch.nStim - 1);
                        if (LEDOFFPeakLoc + 1) < length(avgOFFSegSig)
                            
                            t1ValOFFInc = obj.tTest(avgOFFSegSig(LEDOFFPeakLoc + 1), LEDOFFStartVal, avgOFFSegSTD(LEDOFFPeakLoc + 1), LEDOFFStartSTD);
                            p1ValOFFInc = 1 - tcdf(t1ValOFFInc, obj.patch.nStim - 1);
                            
%                             p2ValOFFInc = 0;
                        elseif (LEDOFFPeakLoc > 2)
                            t1ValOFFInc = obj.tTest(avgOFFSegSig(LEDOFFPeakLoc - 1), LEDOFFStartVal, avgOFFSegSTD(LEDOFFPeakLoc - 1), LEDOFFStartSTD);
                            p1ValOFFInc = 1 - tcdf(t1ValOFFInc, obj.patch.nStim - 1);
                            
                        else
                            p1ValOFFInc = 1;
                            % If the frame is too short, then render
                            % detection of significant fluorescent change
                            % invalid.
                            
                        end
                        
                        
                        if (p1ValOFFInc < pTestVal)
                           obj.OFFInc(allROI(i)) = 1;
                        end
                        
                        % If ON initial value is lower than baseline, then that
                        % means ON causes decrease in cell response
                    end
                    
                    if (pValOFFDec < pTestVal)&(tValOFFDec > 0)&(obj.OFFInc(allROI(i)) ~= 1)
                        % If there is significant decrease in DFF, then
                        % mark the Dec indicator as 1
                        
                        if length(avgOFFSegSig) > 1
                            
                            t1ValOFFDec = obj.tTest(avgOFFSegSig(2), LEDOFFEndVal, avgOFFSegSTD(2), LEDOFFEndSTD);
                            p1ValOFFDec = 1 - tcdf(t1ValOFFDec, obj.patch.nStim - 1);
                            
                            % Compare the second value to the last one to
                            % make sure that significant difference in
                            % first and last frame values is not due to
                            % noise.
                            
%                             p2ValOFFInc = 0;
                        else
                            p1ValOFFDec = 1;
                            % If the frame is too short, then render
                            % detection of significant fluorescent change
                            % invalid.
                            
                        end
                        
                        if (p1ValOFFDec < pTestVal)
                            obj.OFFDec(allROI(i)) = 1;
                        end
                        
                    end
                    
                    % On segment initial value minus baseline
                    %                     OFFPeakBl = LEDOFFPeak - LEDOFFBlVal;
                    % light bleaching will cause OFFPeakBl almost always be
                    % negative
                    %                     OFFPeakInit = LEDOFFPeak - LEDOFFStartVal;
                    
                    
                    
                    %                     if (OFFPeakInit > 5)
                    
                      % If the difference between ON segment peak and initial
                        % value is 3 times larger than the difference between
                        % the initial value and baseline, then it is reasonable
                        % to assume that there is light ON increase.

                    
                    
                    if obj.ONInc(allROI(i))&obj.OFFInc(allROI(i))
                        obj.RGCType{i} = 'ON-OFF';
                    elseif obj.ONInc(allROI(i))|((obj.OFFDec(allROI(i)))&(length(avgONSegSig) < 3))
                        % Light ON response if the cell has increase in
                        % fluorescence in ON frames, or, if the ON frames
                        % are shorter than 3 frames, there is decrease in
                        % fluorescence in OFF frames.
                        obj.RGCType{i} = 'ON';
                    elseif obj.OFFInc(allROI(i))|((obj.ONDec(allROI(i)))&(length(avgOFFSegSig) < 3))
                        % Light OFF response if the cell has increase in
                        % fluorescence in OFF frames, or, if the OFF frames
                        % is shorter than 3 frames, there is decrease in
                        % fluorescence in ON frames.
                        obj.RGCType{i} = 'OFF';
                    else
                        obj.RGCType{i} = 'NaN';
                    end
                    
                end
                
            end
        end
        
        
        function RBBgSubBlCorrSig = get.bgSubBlCorrSig(obj) % photobleach correction
%             frame = 1:size(obj.rawSig,2);

%             if isempty(obj.bgSubSig)
%                obj.bgROISub; 
%             end
            
            frame = 1:size(obj.RBBgSubSig,2);

            %             for ROIIndex = 1:size(obj.rawSig,1)
            initialVal = obj.RBBgSubSig(obj.cROIN, 1);
            finalVal = obj.RBBgSubSig(obj.cROIN, end);
            diffFluo = abs(finalVal - initialVal);
            %                     checkPhotoBleach
            fluo = obj.RBBgSubSig(obj.cROIN, :);
            if diffFluo > 10
                
                RBBgSubBlCorrSig = detrend(fluo) + fluo(1);
                
%                 smFluo = smooth(fluo,20);
%                 p = polyfit(frame(1:end/2),smFluo(1:end/2)',1);
                
                %                     yFit =  p(1)*frame + p(2);
                
                %                     photobleachCorrSig = fluo - yFit;
%                 rawBlCorrSig = fluo - p(1)*frame;
            else
                RBBgSubBlCorrSig = fluo;
            end
            
            %             end
        end
        
        function getLEDFrame(obj)
            % This is a simple function that converts patchData object's
            % LEDONTime and LEDOFFTime (Digitiser form) into LEDONFrame and
            % LEDOFFFrame (image form)
            for i = 1:length(obj.patch.LEDONTime)
                %% Light ON frames
                firstFrameI = find(obj.frameTime > obj.patch.LEDONTime{i}(1));
                lightONFirstFrame = firstFrameI(1);
                
                lastFrameI = find(obj.frameTime < obj.patch.LEDONTime{i}(end));
                lightONLastFrame = lastFrameI(end);
                
                obj.LEDONFrame{i} = lightONFirstFrame:lightONLastFrame;
                
                %% Light OFF frames
                firstFrameI = find(obj.frameTime > obj.patch.LEDOFFTime{i}(1));
                lightOFFFirstFrame = firstFrameI(1);
                
                lastFrameI = find(obj.frameTime < obj.patch.LEDOFFTime{i}(end));
                lightOFFLastFrame = lastFrameI(end);
                
                obj.LEDOFFFrame{i} = lightOFFFirstFrame:lightOFFLastFrame;
                
            end
            
        end

        
        function LEDRemoval(obj,ROIN)
            % Removing LED artefact for background of local maxima
            
            if size(obj.bgVal,1) > 2
                bgLEDOnIntAvg = mean(obj.bgVal(ROIN,obj.LEDOnMFrame));
            else
                bgLEDOnIntAvg = mean(obj.bgVal(obj.LEDOnMFrame));
                %                 realLEDOnFrame = obj.LEDOnFrame;
            end
            
            blFrame = 1:obj.frames;
            blFrame(ismember(blFrame,obj.LEDOnMFrame)) = [];
            blFrame(ismember(blFrame,obj.LEDOnEdgeFrame)) = [];
            
            if size(obj.bgVal,1) > 2
                bgLEDOffIntAvg = mean(obj.bgVal(ROIN,blFrame));
            else
                bgLEDOffIntAvg = mean(obj.bgVal(1,blFrame));
            end
            
            bgLEDCorrVal = bgLEDOnIntAvg - bgLEDOffIntAvg;
            
            obj.bgLEDRemoveInt = obj.bgVal(1,:);
            %                 obj.bgLEDRemoveInt = obj.bgLineVal(cellPos,:);
            if ~isempty(obj.LEDOnMFrame)
                obj.bgLEDRemoveInt(obj.LEDOnMFrame) = obj.bgLEDRemoveInt(obj.LEDOnMFrame) -...
                    padarray(bgLEDCorrVal,[size(obj.LEDOnMFrame,1) - 1,size(obj.LEDOnMFrame,2) - 1],'replicate','post');
            end
            
            %                 bgLEDCorrVal = bgLEDOnIntAvg - bgLEDOffIntAvg;
            for LEDNo = 1:length(obj.LEDOnFrameCell)
                % bring down those intensities partially affected by
                % LED by taking average values of its baseline before
                % LED
                if ~isempty(obj.LEDOnFrameCell{1,LEDNo})
                    if (obj.LEDOnFrameCell{1,LEDNo}(1) > 2)
                        obj.bgLEDRemoveInt(obj.LEDOnFrameCell{1,LEDNo}(1)) = obj.bgLEDRemoveInt(obj.LEDOnFrameCell{1,LEDNo}(1) - 2);
                    end
                    if (obj.LEDOnFrameCell{1,LEDNo}(end) < obj.frames - 2)
                        obj.bgLEDRemoveInt(obj.LEDOnFrameCell{1,LEDNo}(end)) = obj.bgLEDRemoveInt(obj.LEDOnFrameCell{1,LEDNo}(end) + 2);
                    end
                end
            end
            
            %             cellPos = obj.sigROI{3,ROIN}(1);
            
            
            % Removing LED artefact for local maxima
            
            LEDOnIntAvg = mean(obj.rawSig(ROIN,obj.LEDOnMFrame),2);
            LEDOffIntAvg = mean(obj.rawSig(ROIN,blFrame),2);
            obj.LEDReSig = obj.rawSig(ROIN,:);
            
            LEDCorrVal = LEDOnIntAvg - LEDOffIntAvg;
            
            if ~isempty(obj.LEDOnMFrame)
                obj.LEDReSig(obj.LEDOnMFrame) = obj.LEDReSig(obj.LEDOnMFrame) -...
                    padarray(LEDCorrVal,[0,size(obj.LEDOnMFrame,2) - 1],'replicate','post');
            end
            
            % for those uncorrected LED frame or false corrected frame,
            % simply use diff to find them and change their value to the
            % average of their neighbours
            
            for i = 1:length(obj.LEDOnFrameCell)
                % bring down those intensities partially affected by
                % LED by taking average values of its baseline before
                % LED
                if ~isempty(obj.LEDOnFrameCell{1,i})
                    if (obj.LEDOnFrameCell{1,i}(1) > 2)
                        obj.LEDReSig(obj.LEDOnFrameCell{1,i}(1)) = obj.LEDReSig(obj.LEDOnFrameCell{1,i}(1) - 2);
                    end
                    if (obj.LEDOnFrameCell{1,i}(end) < obj.frames - 2)
                        obj.LEDReSig(obj.LEDOnFrameCell{1,i}(end)) = obj.LEDReSig(obj.LEDOnFrameCell{1,i}(end) + 2);
                    end
                end
            end
            
        end
        
        
        
        
        function anaCaSig(obj,smMethod, smFac)
            % Calculating the calcium signal properties, ie, baseline noise
            % level, rise time, decay time and maximum amplitude
            
            for ROINo = 1:obj.nROI
                stimFrameSort(obj, ROINo)
                for stimNo = 1:obj.nStim
                    obj.getSegSig;
                    % Extracting the frames for a particular stimulation
                    wholeSigFrame = obj.stimFrame{7,stimNo}(1,:);
                    if isempty(obj.stimFrame{1,stimNo})
                        continue
                    end
                    blSigFrame = obj.stimFrame{1,stimNo}(1,:) - obj.stimFrame{1,stimNo}(1,1) + 1;
                    stimSigFrame = obj.stimFrame{2,stimNo}(1,:) - obj.stimFrame{1,stimNo}(1,1) + 1;
                    offSigFrame = obj.stimFrame{3,stimNo}(1,:) - obj.stimFrame{1,stimNo}(1,1) + 1;
                    stimCoupledWindowFrame = stimSigFrame(1):(stimSigFrame(end) + ceil(obj.stimCoupleWindow*obj.fps));
                    
                    DFF = obj.DFFSegSig{ROINo,stimNo};
                    obj.wholeSig{ROINo, stimNo} =  DFF;
                    
                    if nargin < 2
                        data =  DFF;
                    elseif smMethod == 'm'
                        data = smooth(DFF,smFac);
                    elseif smMethod == 'e'
                        data = smoothts(DFF,'e', smFac);
                    elseif smMethod == 's'
                        data = smooth(DFF,'sgolay',smFac);
                    end
                    
                    obj.blSig{ROINo, stimNo} = data(blSigFrame);
                    obj.stimSig{ROINo, stimNo} =  data(stimSigFrame);
                    obj.offSig{ROINo, stimNo} =  data(offSigFrame);
                    stimCoupledSig = data(stimCoupledWindowFrame);
                    % stimulation frames plus 200ms of offset frame for
                    % finding signal maximum amplitude and also the frame at
                    % which signal is 3 times above baseline noise level
                    
                    
                    obj.blNoiseLevel(ROINo, stimNo) = std(obj.blSig{ROINo, stimNo});
                    obj.blMean(ROINo, stimNo) = mean(obj.blSig{ROINo, stimNo});
                    
                    
                    %% Calcium signal: maximum amplitude
                    [peakValue,obj.maxAmpLocs(ROINo, stimNo)] = ...
                        max(stimCoupledSig);
                    framesForMax =...
                        (obj.maxAmpLocs(ROINo, stimNo) - 1):(obj.maxAmpLocs(ROINo, stimNo) + 6);
                    % The maximum is calculated as the mean value in a time
                    % window around the peak response amplitde (1 frame
                    % before and 6 frames after the peak).
                    framesForMax(framesForMax < 1) = [];
                    framesForMax(framesForMax > length(stimCoupledSig)) = [];
                    obj.maxAmp(ROINo, stimNo) = mean(stimCoupledSig(framesForMax));
                    
                    obj.maxAmpLocs(ROINo, stimNo) = obj.maxAmpLocs(ROINo, stimNo) + stimSigFrame(1) - 1;
                    % correcting the max amplitude frame number. Before I
                    % was getting the max amplitude from stimCoupledSig, so
                    % the frame was counted from 1. However, in overall
                    % segmented signal, the max amplitude location is
                    % referred to the frame number counting from the
                    % baseline frame, so need to correct it here.
                    
                    
                    
                    
                    %                    obj.decayFrame{ROINo, stimNo} = obj.maxAmpLocs(ROINo, stimNo);
                    
                    obj.stimCoupledRes(ROINo, stimNo) = obj.maxAmp(ROINo, stimNo) > 3*obj.blNoiseLevel(ROINo, stimNo);
                    aboveThreshTemp = find(stimCoupledSig > 3*obj.blNoiseLevel(ROINo, stimNo));
                    if ~isempty(aboveThreshTemp)
                        obj.aboveThreshTime(ROINo, stimNo) = aboveThreshTemp(1);
                    end
                    
                    %% Calcium signal: rise time
                    obj.riseTime(ROINo, stimNo) = obj.maxAmpLocs(ROINo, stimNo)/obj.fps;
                    
                    halfMax = 0.5*obj.maxAmp(ROINo, stimNo);
                    halfMaxRiseTemp = find(stimCoupledSig > halfMax);
                    %                    obj.halfRiseFrame(ROINo, stimNo) =  halfMaxRiseTemp(1);
                    %                    obj.halfRiseTime(ROINo, stimNo) = obj.halfRiseFrame(ROINo, stimNo)/obj.fps;
                    % the rise time is defined to be
                    obj.riseFrame{ROINo, stimNo} = stimSigFrame(1):obj.maxAmpLocs(ROINo, stimNo);
                    % The rising frame is from the start of stimulation to
                    % the peak of calcium signal
                    obj.riseSig{ROINo, stimNo} = data(obj.riseFrame{ROINo, stimNo})';
                    
                    
                    %% Calcium signal: decay time
                    %                    decayTemp = find(obj.offSig{ROINo, stimNo} > halfMax);
                    % decayTemp is the frame at which the peak fluorescence
                    % has dropped to half of the peak value.
                    %                    obj.decayTime(ROINo, stimNo) = decayTemp(end)/obj.fps;
                    
                    peakToReturnBlDiff = peakValue - DFF(end);
                    
                    %                    backToBlFrame =...
                    %                        find((obj.offSig{ROINo, stimNo} - DFF(end)) < 0.1*peakToReturnBlDiff);
                    obj.decayFrame{ROINo, stimNo} = (obj.maxAmpLocs(ROINo, stimNo) + 1):offSigFrame(end);
                    obj.decaySig{ROINo, stimNo} = data(obj.decayFrame{ROINo, stimNo})';
                    % The decay frame is from the peak fluorescence frame to
                    % the frame at which the fluoresence values dropped
                    % to 10% of the peak - baseline (returning baseline)
                    % value.
                    
                    
                    
                    %                    obj.maxAmp(ROINo, stimNo) = max(obj.stimSig{ROINo, stimNo}) - obj.blMean(ROINo, stimNo);
                    
                    %                        obj.maxAmp(ROINo, stimNo) = 0;
                    
                    % Note the baseline mean is calculated for each stim
                    % segment because previous baseline is calculated as the
                    % mean of the whole signal, which is not accurate to
                    % describe each stim segments' baseline.
                    
                    %                    obj.stimCoupledRes(ROINo, stimNo) = obj.maxAmp(ROINo, stimNo) > 3*obj.blNoiseLevel(ROINo, stimNo);
                    
                    %                    riseTimeLocsTemp = ...
                    %                        find((obj.stimAndOffSig{ROINo, stimNo}- obj.blMean(ROINo, stimNo)) > 0.9*obj.maxAmp(ROINo, stimNo));
                    %                    riseTimeLocs(ROINo, stimNo) = riseTimeLocsTemp(1);
                    
                    
                end
            end
        end
        
        
        function labelCoupling(obj, ROIN, coupleLength)
            
            obj.cellImageTime = obj.frameTime(1:obj.frames) - (obj.patch.lineATime)*(obj.height - obj.sigROI{3,ROIN}(1));
            
            coupledPts = [];
            coupledPtTime = [];
            coupledDiffTime = [];
            
            for i = 1:length(obj.patch.stimStartTime)
                stimStart = obj.patch.stimStartTime(i);
                stimEnd = obj.patch.stimEndTime(i);
                tempStartPt = find(obj.cellImageTime > stimStart);
                
                tempEndPt = find(obj.cellImageTime < stimEnd);
                
                if (tempStartPt(1) > coupleLength)&(tempStartPt(1) < obj.frames - (coupleLength - 2))
                    tempStartCPts = (tempStartPt(1) - coupleLength):(tempStartPt(1) + (coupleLength - 1));
                elseif (tempStartPt(1) > coupleLength)
                    tempStartCPts = (tempStartPt(1) - coupleLength):obj.frames;
                elseif (tempStartPt(1) < obj.frames - (coupleLength - 2))
                    tempStartCPts = 1:(tempStartPt(1) + (coupleLength - 1));
                else 
                    tempStartCPts = 1:obj.frames;
                end
                
                
                if (tempEndPt(end) > (coupleLength - 1))&(tempEndPt(end) < obj.frames - (coupleLength - 1))
                    tempEndCPts = (tempEndPt(end) - (coupleLength - 1)):(tempEndPt(end) + (coupleLength));
                elseif (tempEndPt(end) > (coupleLength - 1))
                    tempEndCPts = (tempEndPt(end) - (coupleLength - 1)):obj.frames;
                elseif (tempEndPt(end) < obj.frames - (coupleLength - 1))
                    tempEndCPts = 1:(tempEndPt(end) + (coupleLength));
                else 
                    tempEndCPts = 1:obj.frames;
                end
                
                
                        
                coupledPts = [coupledPts tempStartCPts tempEndCPts];
                
                
                coupledDiffTime = [coupledDiffTime, (obj.cellImageTime(tempStartCPts) - stimStart), obj.cellImageTime(tempEndCPts) - stimEnd];
                
            end
            
            coupledPtTime = obj.cellImageTime(coupledPts);
            
            coupledDiffTime = coupledDiffTime*1000;
            % convert to millisecond.
            
            gca
            
            labelpoints (coupledPtTime, obj.DFFSig(ROIN,coupledPts), coupledDiffTime);
            
            
        end
        
%         function plotMain(obj,plotType, ROIN, stimN, smMethod,smFac)
        function plotMain(obj, varargin)
            % plot all signals in \DeltaF/F0 form
            % if ROIN = 0, plot background ROI signal
            
            % If 'plotType', then next arg should be:
            % 'DFF' (default), 'raw','DFFSeg, 'rawSeg', 'fitSeg'

            
            % If 'smooth', then next 2 args should be:
            % 1. 'm' (default), 'e', 's'
            % 2. number of smoothing points
            

            
%             gcf;
%             clf;
            
            
            
            
            fitData = [];
            
            ROIN = 1:obj.nROI;
            plotType = 'DFF';
            smMethod = 'm';
            smFac = 5;
            plotSTD = 0;
            stimN = 0;
            
            if nargin > 1 
                for i = 1:length(varargin)
                    
                    if strcmp(varargin{i}, 'plotType')
                        % Check plotType
                        checkArgLength = length(varargin{i});
                        plotType = varargin{i + 1};
                        
                        
%                         if checkArgLength == 1
%                             plotType = varargin{i + 1};
%                             ROIN = 1:obj.nROI;
%                         else
%                             plotType = varargin{i + 1}(1);
%                             if strcmp(plotType, 'DFF')|strcmp(plotType, 'raw');
%                                 ROIN = varargin{i + 1}(2:end);
%                             else
%                                 ROIN = varargin{i + 1}(2:(end - 1));
%                                 stimN = varargin{i + 1}(end);
%                             end
%                         end
                        
                    elseif strcmp(varargin{i}, 'ROIN')
                        ROIN = varargin{i + 1};
%                         obj.cROIN = ROIN;
                        
                    elseif strcmp(varargin{i}, 'stimN')
                        stimN = varargin{i + 1};
                        obj.cStimN = stimN;
                        
                    elseif strcmp(varargin{i}, 'smooth')
                        if length(varargin{i + 1}) == 1
                            % only smooth method, no number of smoothing
                            % points
                            smMethod = varargin{i + 1};
                            smFac = 5;
                        else
                            smMethod = varargin{i + 1}{1};
                            smFac = varargin{i + 1}{2};
                        end
%                         
%                         
%                     elseif strcmp(varargin{i}, 'plotSTD')
%                         plotSTD = 1;
                        
                    end
                    
                end
            end
            
            plotLoc = 1:obj.nROI;
            
            for i = 1:length(ROIN)
                if ROIN(i) == 0
                    rawBgSig = obj.bgVal(1,:);
                    
                    bgBsF = median(rawBgSig);
                    bgDF = rawBgSig - bgBsF;
                    bgDFF = bgDF/bgBsF*100;
                end
                
                
                
                if length(ROIN) > 1
                    %                     subplot(3,ceil(length(ROIN)/3),i);
                    if plotLoc(i) > 16
                        plotLoc(i:end) = plotLoc(i:end) - 16;
                        figure;
                    end
                    subplot(4,4,plotLoc(i));
                end
                
                obj.cROIN = ROIN(i);
                
                frameTime = obj.cellFrameTime;
                if (strcmp(plotType, 'DFF'))|(strcmp(plotType, 'fit'))
                    if ROIN(i) ~= 0
                        obj.cROIN = ROIN(i);
                        plotData = obj.DFFSig;
                    else
                        plotData = bgDFF;
                    end
                    
                elseif strcmp(plotType, 'rawSubBlSub')
                    %% Based on raw signal with bleaching correction
                    if ROIN(i) ~= 0
                        obj.cROIN = ROIN(i);
                        plotData = obj.rawBlCorrSig;
                    else
                        plotData = rawBgSig;
                    end
                    
                elseif strcmp(plotType, 'raw')
                    %% Based on raw signal
                    if ROIN(i) ~= 0
                        plotData = obj.rawSig(ROIN(i),:);
                    else
                        plotData = rawBgSig;
                    end
                    
                elseif strcmp(plotType, 'bgSubBlSub')
                    
                    %% Based on background subtraction and bleaching corerction
                    if ROIN(i) ~= 0
                        obj.cROIN = ROIN(i);
                        plotData = obj.bgSubBlCorrSig;
                    else
                        plotData = rawBgSig;
                    end
                    
                elseif strcmp(plotType, 'bgSub')
                    %% Based on background subtraction signal
                    if ROIN(i) ~= 0
                        if isempty(obj.bgSubSig)
                           obj.bgROISub; 
                        end
                        plotData = obj.bgSubSig(ROIN(i),:);
                    else
                        plotData = rawBgSig;
                    end

                end

                %                 elseif (strcmp(plotType, 'rawSeg'))|(strcmp(plotType, 'DFFSeg')|strcmp(plotType, 'fitSeg'))
                if stimN ~= 0
                    frameTime = obj.segFrameTime;
                    segFrame = obj.segFrame;
                    
                    blSegFrame = obj.blSegFrame;
                    
                    plotData = plotData(segFrame);
                    
                    if strcmp(plotType, 'fit')
                        obj.cROIN = ROIN(i);
                        obj.cStimN = stimN;
                        fitData = obj.fitSig;
                    end
                    
%                     if strcmp(plotType, 'rawSeg')
%                         plotData = obj.rawSig(ROIN(i),segFrame);
%                     elseif strcmp(plotType, 'DFFSeg')
%                         plotData = obj.getDFFSegSig(ROIN(i), stimN);
%                     end
%                     if strcmp(plotType, 'fitSeg')
%                         fitData = obj.fitSig(ROIN(i), stimN, frameTime);
%                         plotData = obj.getDFFSegSig(ROIN(i), stimN);
%                     end
                    
                    
                    
                end
                
                
%                 if nargin > 5
                if obj.patch.ndim == 2
                    smData = obj.getSmData(plotData, smMethod, smFac);
                else
                    for sweepI = 1:length(plotData)
                        smData{sweepI} = ...
                            obj.getSmData(plotData{sweepI}, smMethod, smFac);
                    end
                end
                   
%                 elseif nargin > 4
%                     [smData, smFac] = obj.getSmData(plotData, smMethod);
%                 else
%                     [smData, smFac] = obj.getSmData(plotData);
%                 end

                if smMethod == 'm'
                    smSigPrint = sprintf('Smoothed trace (%d moving smooth)', smFac);
                elseif smMethod == 'e'
                    smSigPrint = sprintf('Smoothed trace (exponential smooth with alpha = %2.1f)', smFac);
                elseif smMethod == 's'
                    smSigPrint = sprintf('Smoothed trace (Savitzky-Golay filter with degree = %d)', smFac);
                end
                
                legendText = {};
                
                if obj.patch.ndim == 2
                    recordN = 1;
                else
                    recordN = length(plotData);
                end
                
                % Plot Patched ROI membrane voltage
                
                for sweepI = 1:recordN
                    
                    obj.cStimN = sweepI;
                    if obj.patch.ndim == 3
                        modVal = mod(sweepI,25);
                        if modVal == 1
                            figure;
                        elseif modVal == 0
                            modVal = 25;
                        end
                       subplot(5,5,modVal)
                    end
                
                    hold on;
                
                %% Plot intracellular stimulation time
                if ~isempty(obj.patch.intStim)
                    intraStimH = obj.plotIntraStim(smData);
                    legendText{2,end + 1} = 'Int stim time';
                    legendText{1,end} = intraStimH;
                end
                
                %% plot patched cell membrane potential
                if ~isempty(obj.patch.patch)
                    intraPatch = obj.plotPatch(smData);
                    legendText{2,end + 1} = 'patched membrane voltage';
                    legendText{1,end} = intraPatch;
                end
                
                %% Plot LED stimulation
                if ~isempty(obj.patch.LED)
                    LEDH = obj.plotLED(smData);
                    legendText{2, end + 1} = 'LED';
                    legendText{1,end} = LEDH;
                end
                
                %% Plot extracellular stimulation information
                if ~isempty(obj.patch.extStim)
                    if ~isempty(obj.stimFilePath)
                        obj.readExtStim(obj.stimFilePath);
                        obj.nStim = length(obj.stimOnFrameCell);
                    end
%                     obj.getStimFrameSort(ROIN(i))
                    if obj.patch.ndim == 2
                        extStimH = obj.plotExtStim(smData);
                        legendText{2,end + 1} = 'Ext stim (rescaled)';
                        legendText{1,end} = extStimH;
                    else
                        extStimH = obj.plotExtStim(smData{sweepI});
                        legendText{2,end + 1} = 'Ext stim (rescaled)';
                        legendText{1,end} = extStimH;
                    end
                end
                
                
                if ~isempty(obj.patch.extStim)
                    if isempty(obj.extStimVar)
                        xlString = sprintf('Time (s); Frames per second: %4.2f', obj.fps);
                    else
                        if strcmp(obj.extStimVar, 'amp')
                            xlString = ...
                                sprintf('Time (s); Frames per second: %4.2f; Stim Duration: %dus; Stim Freq: %3.0fHz; Stim reps: %d',...
                                obj.fps, obj.extStimSegDur(1), obj.extStimSegFreq(1), obj.extStimSegRep(1));
                        elseif strcmp(obj.extStimVar, 'dur')
                            xlString = ...
                                sprintf('Time (s); Frames per second: %4.2f; Stim Amplitude: %duA; Stim Freq: %3.0fHz; Stim reps: %d',...
                                obj.fps, obj.extStimSegAmp(1), obj.extStimSegFreq(1), obj.extStimSegRep(1));
                        elseif strcmp(obj.extStimVar, 'rep')
                            xlString = ...
                                sprintf('Time (s); Frames per second: %4.2f; Stim Amplitude: %duA; Stim Duration: %dus; Stim Freq: %3.0fHz; ',...
                                obj.fps, obj.extStimSegAmp(1), obj.extStimSegDur(1), obj.extStimSegFreq(1));
                        end
                    end
                elseif ~isempty(obj.patch.LED)
                    xlString = sprintf('Time (s); Frames per second: %4.2f', obj.fps);
                end
                
                if obj.patch.ndim == 2
                DFFH = plot(frameTime,plotData,'.g');
                
                smDFFH = plot(frameTime,smData,'color', [1,0,0]);
                if ~isempty(fitData)
                    fitDFFH = plot(frameTime((blSegFrame(end) + 1):end),fitData,'color', [0,0,0]);
                end
                else
                    DFFH = plot(frameTime{sweepI},plotData{sweepI},'.g');
                
                    smDFFH = plot(frameTime{sweepI},smData{sweepI},'color', [1,0,0]);
                end
                
                obj.plotData = plotData;
                 
%                 if plotSTD == 1
%                    obj.plotSTD; 
%                 end
                
                
                if ~isempty(DFFH)
                    legendTemp = cell(2,size(legendText,2) + 1);
                    legendTemp(:,2:end) = legendText;
                    if strcmp(plotType, 'DFF')
                        legendTemp{2,1} = '\DeltaF/F0 trace';
                    elseif strcmp(plotType, 'raw')
                        legendTemp{2,1} = 'Raw trace';
                    elseif strcmp(plotType, 'rawBlSub')
                        legendTemp{2,1} = 'Raw bleach subtracted trace';
                    end
                    legendTemp{1,1} = DFFH;
                    legendText = legendTemp;
                end
                
                if ~isempty(smDFFH)
                    legendTemp = cell(2,size(legendText,2) + 1);
                    legendTemp(:,2:end) = legendText;
                    legendTemp{2,1} = smSigPrint;
                    legendTemp{1,1} = smDFFH;
                    
                    legendText = legendTemp;
                    
                end
                
                if ~isempty(obj.patch.LED)
                    titleString = sprintf('ROI number = %d with LED stimulation. Plot type: %s', ROIN(i), plotType);
                elseif ~isempty(obj.patch.intStim)
                    titleString = sprintf('ROI number = %d with intracellular stimulation. Plot type: %s', ROIN(i), plotType);
                elseif ~isempty(obj.patch.extStim)
                    titleString = sprintf('ROI number = %d with extracellular stimulation. Plot type: %s', ROIN(i), plotType);
                end
                title(titleString,'FontSize', obj.fontSize,'FontWeight', obj.fontWeight);
                
                hold off
                xlabel(xlString,'FontSize', obj.fontSize,'FontWeight', obj.fontWeight);
                
                if strcmp(plotType, 'DFF')
                    ylabel('\DeltaF/F0 (%)','FontSize', obj.fontSize,'FontWeight', obj.fontWeight);
                else 
                    ylabel('Fluorescence (AU)','FontSize', obj.fontSize,'FontWeight', obj.fontWeight);
                end
                
                
                if (length(ROIN) == 1)&(obj.patch.ndim == 2)
%                     hForLegend = cell2mat(legendText(1,:));
                    
                    hLegend = legend([legendText{1,:}], legendText{2,:});
%                     hLegend = legend(hForLegend,legendText(2,:),'Location','NorthWest');
                    set(hLegend, 'FontSize', obj.legendSize,'FontWeight', obj.fontWeight);
                end
                set(gca,'FontSize', obj.tickLabelSize,'FontWeight', obj.fontWeight);
                end
            end
        end
        
        
        function anaSweepCalSig(obj,varargin)
            
            
            blPreTime = 0.5;
            AStimPostTime = 0.3;
            
%             framesPerSweep = 95;
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'blPreTime');
                    blPreTime = varargin{i + 1};
                end
                
                if strcmp(varargin{i}, 'AStimPostTime');
                    AStimPostTime = varargin{i + 1};
                end
                
%                 if strcmp(varargin{i}, 'framesPerSweep');
%                     framesPerSweep = varargin{i + 1};
%                 end
            end
            
            blPrePt = blPreTime*obj.patch.SR;
            AStimPostPt = AStimPostTime*obj.patch.SR;
            
            rawSubBgSig = obj.rawSig - obj.bgVal;
            
            obj.patchedCaSigDFF = zeros(obj.patch.nSweeps, obj.framesPerSweep);
            
            obj.patchedFrameTime = [];
            for stimI = 1:obj.patch.nSweeps
                
%                 rawSubBgSig = obj.rawSig - obj.bgVal;
                
                % Getting the frame times
                frameTrigger = obj.patch.frame(stimI,:);
                frameDiff = diff(frameTrigger);
                
                framePt = find(frameDiff < -1) + 1;
                obj.patchedFrameTime(stimI,:) = framePt/obj.patch.SR;
                 
                
                
                % Getting baseline calcium signal. This is by average
                % frames 0.5s before the stimulation
                
                stimStartPt = obj.patch.extStimStartPt(stimI);
                
                if stimStartPt == 0
                    continue
                end
                blPt = stimStartPt - blPrePt;
                blFrame = (stimI - 1)*obj.framesPerSweep + find((framePt > blPt)&(framePt < stimStartPt));
                
                FZero = mean(rawSubBgSig(blFrame));
                obj.patchedCaBlMean(stimI) = FZero;
                
                % Getting calcium signal. This is done by calculating the
                % area (with baseline as 0) from the time of stimulation
                % till 1 s after the offset of the stimulation
                
                sigFrame = ((stimI - 1)*obj.framesPerSweep + 1):stimI*obj.framesPerSweep;
                
                sigDFF = (rawSubBgSig(sigFrame) - FZero)/FZero;
                obj.patchedCaSigDFF(stimI,:) = sigDFF;
                
                stimEndPt = obj.patch.extStimEndPt(stimI);
                
                AStimPt = stimEndPt + AStimPostPt;
                
                CaSigFrame = find((framePt > stimStartPt)&(framePt < AStimPt));
                

                
                obj.patchedCaArea(stimI) = sum(sigDFF(CaSigFrame));
                
            end
            
            
        end
        
        function plotStim(obj,varargin)
            %% this is to simply plot the electrophysiological and imaging data recorded simultaneously
            
            
            ROII = 1;
            
            stimN = 1;
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'ROII');
                    ROII = varargin{i + 1};
                end
                
                if strcmp(varargin{i}, 'stimN');
                    stimN = varargin{i + 1};
                end
            end
            
            startFrameI = (obj.framesPerSweep*(stimN - 1) + 2);
            % Get rid of the first frame because it is usually dimmer than
            % others due to the initial openning of the imaging LED
            endFrameI = (obj.framesPerSweep*stimN);
            
            frameForPlot = startFrameI:endFrameI;
            
%             plot(obj.frameTime(frameForPlot), obj.rawSig(ROII, frameForPlot));
            
            hold on

            
%             plot(obj.patch.recordTime(stimN,:), obj.patch.patch(stimN,:));

            frameIs = 2:95;
            frameTrigger = obj.patch.frame(stimN,:);
            frameDiff = diff(frameTrigger);
            
            framePt = find(frameDiff < -1) + 1;
            frameTime = framePt/obj.patch.SR;
            
            [haxes,hline1,hline2] = plotyy(frameTime(frameIs), obj.rawSig(ROII, frameForPlot),...
                obj.patch.recordTime(1,:), obj.patch.patch(stimN,:));
            
%             set(hline1, 'LineStyle', '.');
            
            
        end
        
        function plotSweepStimAvg(obj, varargin)
            % This is for recording that used episodic sweep mode, each
            % trials were repeated 3 times for HFS
            ROIN = 1:obj.nROI; 
            % If ROI number is not specified, plot all ROI signals
            
            
            % This parameters are for HFS recordings
            stimRep = 3;
                        
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'ROIN');
                    ROIN = varargin{i + 1};
                end
                
                if strcmp(varargin{i}, 'stimRep');
                    stimRep = varargin{i + 1};
                end
            end

            
            for ROII = 1:length(ROIN)
                figure;
                obj.cROIN = ROIN(ROII);
                
                plotI = 1;
                
                if isempty(obj.avgSig)
                    obj.getAvgSig;
                end
                
                if isempty(obj.avgSig{ROIN(ROII)})
                    obj.getAvgSig('ROIN',ROIN(ROII));
                end
                
                if isempty(obj.avgSigPeakVal)
                    
                    obj.calcAvgSigProperties;
                end
                
                for stimI = 1:obj.patch.nSweeps
                    subplot(5,5,plotI);
                    obj.cStimN = stimI;
                    
                    if plotI < size(obj.avgSig, 2)
                        if isempty(obj.avgSig{ROIN(ROII),plotI})
                            continue
                        end
                    else
                        continue
                    end
                                        
                    modVal = mod((stimI - 1), stimRep);
                    
                    if modVal == 0
                        if stimI == 1
                            plot(obj.avgSig{ROIN(ROII),plotI});
                        else
                            shadedErrorBar((1:length(obj.avgSig{ROIN(ROII),plotI}))/obj.fps,obj.avgSig{ROIN(ROII),plotI},...
                                obj.avgSigSE{ROIN(ROII),plotI},'-k', 0);
                            hold on
                            
                            peakCaSigTime = (obj.stimStartFrame + obj.avgSigPeakFrame(ROIN(ROII),plotI) - 1)/...
                                obj.fps;
                            
                            plot(peakCaSigTime, obj.avgSigPeakVal(ROIN(ROII),plotI),'.r');
                            
                        end
                        plotI = plotI + 1;
                    end
                                       
                end
            end
            
        end
        
        function getAvgSig(obj,varargin)
            % This is for recording that used episodic sweep mode, each
            % trials were repeated 3 times for HFS
            ROIN = 1:obj.nROI;
            % If ROI number is not specified, plot all ROI signals
            
            % This parameters are for HFS recordings
            
            stimRep = 3;
            
            obj.avgSig = cell(length(ROIN), (obj.patch.nSweeps - 1)/stimRep);
            
            % If the stimulation number is not specified, plot all
            % stimulation averaged Ca signals
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'ROIN');
                    ROIN = varargin{i + 1};
                end
                
                if strcmp(varargin{i}, 'stimRep');
                    stimRep = varargin{i + 1};
                end
            end
            
            
            
            for ROII = 1:length(ROIN)
                obj.cROIN = ROIN(ROII);
                ampI = 1;
                
                for stimI = 1:obj.patch.nSweeps
                    obj.cStimN = stimI;
                    
                    
                    sweepSig = obj.DFFSegSig;
                    
                    if isempty(sweepSig)
                        continue;
                    end
                    sweepSTD = std(sweepSig);
                    
                    modVal = mod((stimI - 1), stimRep);
                    
                    if stimI ~= 1
                        % Averaging the sweeps by trial repetitions. the
                        % first entry is the spontaneous recording so is not
                        % included for averaging.
                        if modVal == 1
                            avgSig = zeros(1, length(sweepSig));
                            avgSigSTD = zeros(1, length(sweepSig));
                            % Start averaging with the first trial (sweep 2,
                            % 5, 8,....)
                        end  
                        
                        if length(sweepSig) < length(avgSig)
                            sweepSig((end + 1):length(avgSig)) = 0;
                        end
                        
                        avgSig = avgSig + sweepSig;
                        avgSigSTD = avgSigSTD + sweepSTD;
                    end
                    
                                        
                    if modVal == 0
                        if stimI == 1
%                             plot(sweepSig);
                            obj.avgSig{ROIN(ROII),ampI} = sweepSig;
                        else
                            avgSigSE = avgSigSTD/sqrt(stimRep);
%                             shadedErrorBar((1:length(avgSig))/obj.fps,avgSig, avgSigSE,'-k', 0);
                            obj.avgSig{ROIN(ROII),ampI} = avgSig;
                            obj.avgSigSTD{ROIN(ROII),ampI} = avgSigSTD;
                            obj.avgSigSE{ROIN(ROII),ampI} = avgSigSE;
                        end
                        ampI = ampI + 1;
                    end
                    
                end
            end
            
        end
        
        function plotStimAvg(obj, varargin)
            % This is for recording that used gap-free recording mode where
            % the stimulation segments were decided by triggers
            
            ROIN = 1:obj.nROI;
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'smooth');
                    smoothMethod = varargin{i + 1}{1};
                    smFac = varargin{i + 1}{2};
                end
                if strcmp(varargin{i}, 'ROIN');
                    ROIN = varargin{i + 1};
                end
            end
            
%             obj.cROIN = ROIN;
            if isempty(obj.avgSigPeakVal)
               obj.calcAvgSigProperties; 
            end
            
            for ROII = 1:length(ROIN)
                
                modI = mod(ROII, 25);
                
                if modI == 1
                    figure;
                end
                
                if modI == 0
                    modI = 25;
                end
                    
                
                if length(ROIN) > 1
                    subplot(5,5, modI)
                end
                
                obj.cROIN = ROIN(ROII);
                
                avgSegSig = obj.getAvgSegSig;
                avgSegSTD = obj.getAvgSegSTD;
                
                %             avgSegSTD = [avgSegSig + avgSegSTD; (avgSegSig - avgSegSTD)];
                
                shadedErrorBar((1:length(avgSegSig))/obj.fps,avgSegSig, avgSegSTD,'-k', 0);
                hold on
                
                adjpeakFrame = obj.avgSigPeakFrame(ROIN(ROII))+ obj.stimStartFrame - 1;
                plot(adjpeakFrame/obj.fps , obj.avgSigPeakVal(ROIN(ROII)), '.r', 'MarkerSize', 20);
                
                extStimTime = (obj.patch.extStimStartPt(5) - obj.patch.STSStartPt(5))/obj.patch.SR;
                
                line([extStimTime, extStimTime], [-5, 10], 'Color','r');
                
                if length(ROIN) == 1
                    xlString = sprintf('Time (s); Frames per second: %4.2f', obj.fps);
                    titleString = sprintf('ROI number = %d, Averaged all calcium signal \\DeltaF/F0 with extracellular stimulation', ROIN);
                    
                    title(titleString,'FontSize', obj.fontSize,'FontWeight', obj.fontWeight);
                    
                    xlabel(xlString,'FontSize', obj.fontSize,'FontWeight', obj.fontWeight);
                    ylabel('\DeltaF/F0 (%)','FontSize', obj.fontSize,'FontWeight', obj.fontWeight);
                    set(gca,'FontSize', obj.tickLabelSize,'FontWeight', obj.fontWeight);
                    
                end
            end
            
        end
        
        function plotStimCoupledAvg(obj, ROIN, varargin)
            
            obj.cROIN = ROIN;
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'smooth');
                    smoothMethod = varargin{i + 1}{1};
                    smFac = varargin{i + 1}{2};
                end
            end
            
            [plotSig, ~, stimCoupledI] = obj.getAvgCoupledSig;
            

            
            maxPeak = max(plotSig);
            stringStimCoupledI = num2str(stimCoupledI);
            avgCoupledSigSTD = obj.getAvgCoupledSigSTD;
            
            if sum(plotSig) == 0
                plotSig = obj.getAvgSegSig;
                avgCoupledSigSTD = obj.getAvgSegSTD;
            end
            
            % avgSegSTD = [avgSegSig + avgSegSTD; (avgSegSig - avgSegSTD)];
            
            shadedErrorBar((1:length(plotSig))/obj.fps,plotSig, avgCoupledSigSTD,'-k', 0);
            
            obj.cStimN = 1;
            STF = obj.stimStartFrame/obj.fps;
            
            h = get(gca, 'Children');
            
            line([STF, STF], [-1, maxPeak],'LineWidth',2,'Color',[0 0 1]);
            
            SDLegendString = sprintf('%s 1 SD', setstr(177));
            legend(h([3,4]),{'Averaged calcium signal \DeltaF/F0', SDLegendString});
            
            xlString = sprintf('Time (s); Stimulation coupled segments: [%s]; Frames per second: %4.2f', stringStimCoupledI, obj.fps);
            titleString = sprintf('ROI number = %d. Averaged coupled calcium signal \\DeltaF/F0 with extracellular stimulation', ROIN);
            
            title(titleString,'FontSize', obj.fontSize,'FontWeight', obj.fontWeight);
            
            xlabel(xlString,'FontSize', obj.fontSize,'FontWeight', obj.fontWeight);
            ylabel('\DeltaF/F0 (%)','FontSize', obj.fontSize,'FontWeight', obj.fontWeight);
            set(gca,'FontSize', obj.tickLabelSize,'FontWeight', obj.fontWeight);
            
        end
        
   
        
        
        function plotCalvinFit(obj, ROINo, stimNo)
            
            for ROIN = 1:length(ROINo)
                for stimN = 1:length(stimNo)
                    obj.plotSegSig(ROIN, stimN);
                    hold on;
                    
                    caModel = @(x,p) (p(5) + p(1) .* (x > p(2)) .* (1-exp( -(x - p(2))/p(3) )) .* ...
                        exp( -(x - p(2))/p(4) ));
                    
                    dataTime = [obj.stimSegFrameTimeON{stimN}, obj.stimSegFrameTimeOFF{stimN}];
                    
                    
                    
                    T = dataTime - dataTime(1);
                    %             if ~isnan(obj.segFitInfo{segSigN,1})
                    %                 plot(3880:(3879+length(caModel(T,obj.segFitInfo{segSigN,1}))), caModel(T,obj.segFitInfo{segSigN,1}),'k-')
                    %             end
                    %             if size(obj.segFitInfo, 1) < segSigN
                    %                 fitCalvinMethod(obj,segSigN);
                    %             else
                    if isempty(obj.segFitInfo{ROIN,stimN})
                        obj.fitCalvinMethod;
                    end
                    %             end
                    if ~isnan(obj.segFitInfo{ROIN,stimN})
                        plot(dataTime, caModel(T,obj.segFitInfo{ROIN,stimN}),'r-')
                    end
                end
            end
        end
        
%         function plotROIInfo(obj, plotSigType, plotType, ROIN, stimN)
        function plotROIImagesc(obj, varargin)
            % imagesc figure that plots the calcium signal peak, HDT of all ROIs
            % with ROI locations.
            ROIN = 1:obj.nROI;
            stimN = 0;
            
            
            plotSigType = 'avgSig';
            plotProperty = 'peak';

            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'sigType')
                    plotSigType = varargin{i + 1};
                    
                elseif strcmp(varargin{i}, 'property')
                    plotProperty = varargin{i + 1};
                    
                elseif strcmp(varargin{i}, 'ROIN')
                    ROIN = varargin{i + 1};
                    
                elseif strcmp(varargin{i}, 'stimN')
                    stimN = varargin{i + 1};
                    
                end
            end
                
            
            if (length(ROIN) > obj.nROI)|(ROIN == 0)
                ROIN((obj.nROI + 1):end) = [];
                % If ROIN is more than the number of ROIs, then get whole
                % ROI properties and clear the excess.
            end
            
            imEdge = zeros(obj.height, obj.width);
            
            %             obj.getSigProperties(plotSigType,ROIN,stimN);
            if strcmp(plotSigType, 'DFF')|strcmp(plotSigType, 'smDFF')|strcmp(plotSigType, 'fitDFF')
                % segSignal can be defined as DFF, smDFF or fitDFF
                if isempty(obj.segSigPeakVal)
                    obj.calcSegSigProperties(plotSigType);
                end
                
                peakVal = obj.segSigPeakVal(:,stimN);
                HDT = obj.segSigHDT(:,stimN);
                latency = obj.segSigPeakFrame/obj.fps;
            elseif strcmp(plotSigType, 'avgSig')
                if isempty(obj.avgSigPeakVal)
                    obj.calcAvgSigProperties;
                    
                end
                peakVal = obj.avgSigPeakVal;
                HDT = obj.avgSigHDT;
                latency = obj.avgSigPeakFrame/obj.fps;
            end
            
            BWTotal = zeros(obj.height, obj.width);
            
            if strcmp(plotProperty, 'peak')
                plotSig = peakVal;
                %                 BWTotal = -ones(obj.height, obj.width);
                %                 cmin = -1;
            elseif strcmp(plotProperty, 'HDT')
                plotSig = HDT;
                %                 BWTotal = -0.1*ones(obj.height, obj.width);
                %                 cmin = -0.1;
            elseif strcmp(plotProperty, 'latency')
                plotSig = latency;
                %                 BWTotal = -0.01*ones(obj.height, obj.width);
                %                 cmin = -0.01;
            end
            
            if ~isempty(obj.electrodeROI)
                % plot the stimulating electrode if it is included in the
                % object
                electrodeI = (obj.electrodeROI{1,1} > 0);
                if size(electrodeI) == size(BWTotal)
                    BWTotal(electrodeI) = -1;
                end
            end
            
            for i = 1:length(ROIN)
                if ~isempty(obj.sigROI{1,ROIN(i)})
                    checkOverlap = (BWTotal > 0)&(obj.sigROI{1,ROIN(i)} > 0);
                    BWTotal(checkOverlap) = 0;
                    % Overlapped ROI will have the earlier one dropped down to
                    % 0 so not to sum the peak value.
                    imfull = (obj.sigROI{1,ROIN(i)} > 0);
                    imEdge = imEdge + edge(imfull);
                    
                    BWTotal(imfull) = 0;
                    % make sure that ROI interior value is property value
                    % instead of whatever values that has been assigned to
                    % sigROI BW.
                    
                    BWTotal = BWTotal + obj.sigROI{1,ROIN(i)}*plotSig(ROIN(i));
                end
            end
            

            
            imEdge(imEdge > 0) = 1;
            imEdge = logical(imEdge);
            cmax = max(max(BWTotal));
            
            BWTotal(imEdge) = cmax/2;
            % Draw the edges of the ROI for easier visual inspection
            
            BWTotal(BWTotal == -1) = cmax;
            
            if cmax == 0
                cmax = 0.1
            end
            
            clims = [0, cmax];
            
          
%             surf(BWTotal, clims);
%             contour(BWTotal);
            imagesc(BWTotal, clims);
            colorbar;
%             BWTotal = uint8(obj.addPlotROINum(BWTotal));
%             imshow(BWTotal);
            obj.plotROINum;
            
            axis equal
            axis tight
        end
        
        
   
        
        function plotROINum(obj, imOrg)
            gca;
            boxColor = 'red';
            for i = 1:obj.nROI
                if ~isempty(obj.sigROI{3,i})
                    textPos = obj.sigROI{3,i};
                    ROINum = sprintf('%d',i);
                    %                 im = insertText(imOrg,textPos,ROINum,'BoxColor',boxColor);
                    text(textPos(2), textPos(1), ROINum, 'Color', [1, 0.6, 0],...
                        'FontWeight', 'bold', 'FontSize', 14);
                end
            end
        end
        
        function LEDH = plotLED(obj,data)
            rescaleLED = obj.patch.LED;
            maxLED = max(rescaleLED);
            medianLED = median(rescaleLED);
            rescaleLED = rescaleLED - maxLED;
            rescaleLED = rescaleLED/maxLED;
            rescaleLED = rescaleLED*0.5*max((data) - median(data)) + median(data);
            minRescaleLED = min(rescaleLED);
            maxRescaleLED = max(rescaleLED);
            
            %             LEDOnTime = obj.patch.LEDOnTime;
            %             LEDOnDataTime = ceil(LEDOnTime*obj.patch.SR);
            %             plot(obj.patch.dataTime, obj.patch.LED*0.1*max(diff(data)) + min(data),'.','color',[1,0.5,0]);
            gcf
            hold on
            for i = 1:length(obj.patch.stimStartTime)
                LEDH = plot([obj.patch.stimStartTime(i), obj.patch.stimEndTime(i)],...
                    [minRescaleLED, minRescaleLED], '-b', 'LineWidth', 4);
                % orange colour 'color',[1,0.5,0]
            end
        end
        
        function extStimH = plotExtStim(obj,data,dataPts)
            
            fontSize = 15;
            reductionF = 1;
            % Reduce the number of datapoints drawn to increase speed
            % recrorded external stimulation waveform
            rescaleExtStim = obj.patch.extStim(obj.cStimN,:);
            maxExtStim = max(rescaleExtStim);
%             medianExtStim = median(rescaleExtStim);
%             rescaleExtStim = rescaleExtStim - medianExtStim;
            rescaleExtStim = rescaleExtStim/maxExtStim;
            rescaleExtStim = rescaleExtStim*0.5*max((data) - median(data)) + median(data);
            
            if nargin < 3
                if obj.patch.ndim == 2
                    if obj.patch.dataTime(end) > (obj.frameTime(end) + 20)
                        extStimH = plot(obj.patch.dataTime(1:reductionF:ceil(obj.frameTime(end)*obj.patch.SR)),...
                            rescaleExtStim(1:reductionF:ceil(obj.frameTime(end)*obj.patch.SR)), 'color',[0.5,0.5,1]);
                    else
                        extStimH = plot(obj.patch.dataTime(1:reductionF:end), rescaleExtStim(1:reductionF:end), 'color',[0.5,0.5,1]);
                    end
                else
                        extStimH = plot(obj.patch.dataTime(1:reductionF:ceil(obj.frameTime{obj.cStimN}(end)*obj.patch.SR)),...
                            rescaleExtStim(1:reductionF:ceil(obj.frameTime{obj.cStimN}(end)*obj.patch.SR)), 'color',[0.5,0.5,1]);
                    
                end
            else
                if obj.patch.dataTime(end) > (obj.frameTime(end) + 20)
                    extStimH = plot(obj.patch.dataTime(dataPts(1):reductionF:ceil(obj.frameTime(end)*...
                        obj.patch.SR)), rescaleExtStim(dataPts(1):reductionF:ceil(obj.frameTime(end)*...
                        obj.patch.SR)), 'color',[0.5,0.5,1]);
                else
                    extStimH = plot(obj.patch.dataTime(dataPts),...
                        rescaleExtStim(dataPts), 'color',[0.5,0.5,1]);
                end
            end
            
            for stimN = 1:obj.nStim
                
                if strcmp(obj.extStimVar, 'amp')
                    stimString = sprintf('%duA', obj.extStimSegAmp(stimN));
                elseif strcmp(obj.extStimVar, 'dur')
                    stimString = sprintf('%dus', obj.extStimSegDur(stimN));
                elseif strcmp(obj.extStimVar, 'rep')
                    stimString = sprintf('%dpulses', obj.extStimSegRep(stimN));
                end
                
                yAxisTemp = ylim;
                
                
%                 text(obj.stimSegFrameTimeON{stimN}(1), yAxisTemp(1) + 0.1*diff(yAxisTemp) , stimString, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'FontSize', fontSize);
            end
        end
        
        
        
        function intraStimH = plotIntraStim(obj,data)
            gca
            fontSize = 15;
            patch = obj.patch;
            for i = 1:size(patch.recordTime, 1)
                intraStimH = line([patch.stimTime(i,1),patch.stimTime(i,2)], [min(data),min(data)],'Color','k','LineWidth',5);
                intStimAmpString = sprintf('%dpA', patch.stimInitAmp + i*patch.stimDAmp);
                text(patch.stimTime(i,1), min(data)-0.5,intStimAmpString, 'Color','k','FontSize',fontSize);
                %     plot(patch.stimTime(i,:),-70,'-','LineWidth',20);
            end
        end
        
        function intraReH = plotPatch(obj,data)
            % intraRe = intracellular voltage recording
            if obj.patch.ndim == 2
                %                 plot(obj.patch.dataTime, (obj.patch.patch - median(obj.patch.patch))*0.03*max(diff(data)) + 0.95*median(data),'color','b');
                if nargin > 1
                    intraReH = plot(obj.patch.dataTime, (obj.patch.patch - median(obj.patch.patch))*0.05*max(diff(data)) + 0.5*min(data),'color','b');
                else
                    intraReH = plot(obj.patch.dataTime, obj.patch.patch,'color','b');
                end
            elseif obj.patch.ndim == 3
                if nargin > 1
                    %                 plot(obj.patch.patchGapFree(1,:), (obj.patch.patchGapFree(2,:) - median(obj.patch.patchGapFree(2,:)))*0.03*max(diff(data)) + 0.8*median(data),'color','b');
                    intraReH = plot(obj.patch.patchGapFree(1,:), (obj.patch.patchGapFree(2,:) - median(obj.patch.patchGapFree(2,:)))*0.3*max(diff(data)) + 0.5*min(data),'color','b');
                else
                    intraReH = plot(obj.patch.patchGapFree(1,:), obj.patch.patchGapFree(2,:),'color','b');
                end
            end
        end
        
        %         function rsIm(obj, clampObj)
        %             if ~isempty(obj.averageInt)
        %                 obj.rsAvgInt = resample(obj.averageInt, clampObj.rsSR, round(obj.fps));
        %
        %             end
        %
        %             obj.nRS = ceil(size(obj.averageInt,1)*clampObj.rsSR/round(obj.fps));
        %             obj.rsSR = clampObj.rsSR;
        %             obj.frameTime = 1/obj.rsSR:1/obj.rsSR:(obj.nRS*1/obj.rsSR);
        % %             obj.syncImIndex = zeros(1, length(imFrameTime));
        % %             for i = 1:length(imFrameTime)
        % %                 obj.syncImIndex(i) = imFrameTime/
        % %
        % %
        % %             end
        %         end
        
        function saveImage(obj, saveName)
            h = gcf;
            set(h,'position',get(0,'screensize'))
            [pathstr,name,ext] = fileparts(obj.path);
            saveNameSpec = sprintf('_Clampex_%d_2P_%d', obj.pClampNo, obj.imageNo);
            saveName = strcat(saveName,saveNameSpec);
            saveFullName = fullfile(pathstr,saveName);
            saveas(h,saveFullName,'tif');
            saveas(h,saveFullName,'fig');
        end
        
        function changeXLabel(obj,labelName)
            h = gcf;
            ax = findall(h, 'type', 'axes');
            
            xlStringTemp = sprintf('Time (s); Frames per second: %4.2f; Stim Protocol: ', obj.fps);
            xlStringAll = strcat(xlStringTemp, labelName);
            
            for nAxe = 1:length(ax)
                xlabel(ax(nAxe),xlStringAll);
            end
            
            
            
            
        end
        
        
        
        function initialMatrix(obj)
            %             for i = 1:obj.frames
            %                 obj.frame = i;
            %                 obj.initialData(:,:,i) = obj.image;
            %             end
            obj.initialData = zeros(obj.height,obj.width,obj.frames,'uint16'); % 3D matrix containing 16bits values
            TiffObj = Tiff(obj.path,'r'); % open file for reading
            for i = 1:obj.frames
                TiffObj.setDirectory(i);
                obj.initialData(:,:,i) = TiffObj.read();
            end
        end
        
        function initialTwoDMatrix(obj)
            for i = 1:obj.frames
                obj.frame = i;
                obj.initialTwoDData(:,i) = reshape(obj.image,obj.height*obj.width,1);
            end
        end
        
        function cellFrameTime = get.cellFrameTime(obj)
            if obj.patch.ndim == 2
                if ~isempty(obj.patch.lineATime)
                    cellFrameTime = obj.frameTime(1:obj.frames) - (obj.patch.lineATime)*(obj.height - obj.sigROI{3,obj.cROIN}(1));
                else
                    cellFrameTime = obj.frameTime(1:obj.frames);
                end
            else
                cellFrameTime = obj.frameTime;
            end
            %             temp = obj.DFFSig;
        end
        
        
        function DFFSig = get.DFFSig(obj)
%             DFFSig = zeros(obj.nROI,size(obj.rawSig,2));
%             rawBlCorrSig = obj.rawBlCorrSig;

            RBBgSubBlCorrSig = obj.bgSubBlCorrSig;
            
%             bgSubSig = obj.bgSubSig;
%             rawSignal = rawBlCorrSig;
            
            blFrames = [];
            
            if obj.patch.ndim == 2
                for i = 1:obj.nStim
                    obj.cStimN  = i;
                    blFrames = [blFrames, obj.blSegFrame];
                end
                
                bsF = median(RBBgSubBlCorrSig(blFrames));
                DF = RBBgSubBlCorrSig - bsF;
                DFFSig = DF/abs(bsF)*100;
            else
                for i = 1:obj.patch.nSweeps
                    obj.cStimN  = i;
                    DFFSig{i} = obj.DFFSegSig;
                end
                
            end
            

            %             end
        end
        
        function [smData smFac] = getSmData(obj, inData, smMethod, smFac)
            if nargin < 4
                smFac = 5;
            end
            
            if nargin < 3
                smData = smooth(inData,smFac);
            elseif smMethod == 'm'
                smData = smooth(inData,smFac);
            elseif smMethod == 'e'
                smData = smoothts(inData,'e',smFac);
            elseif smMethod == 's'
                smData = smooth(inData,'sgolay',smFac);
            end
        end
        
        function DFFSegSig = get.DFFSegSig(obj)
            if obj.patch.ndim == 2
                segFrame = obj.segFrame;
                %             rawSegSig = obj.rawBlCorrSig(segFrame);
              
                bgSubBlCorrSig = obj.bgSubBlCorrSig(segFrame);
                
                blFrame = obj.blSegFrame;
                
                bsF = median(bgSubBlCorrSig(blFrame));
                DF = bgSubBlCorrSig - bsF;
                DFFSegSig = DF/bsF*100;
            else
                
                if isempty(obj.sweepFrames)
                    obj.sweepFrames = obj.patch.sweepFrames;
                    % If some imaging frames were not recorded by
                    % digitiser, need to manually input the image frame
                    % number (as opposed to trigger recorded frames, which
                    % is decided by ClampexData)
                end
                initialFrameNum = obj.sweepFrames*(obj.cStimN - 1) + 1;
                
                if (strcmp(obj.cDate, '20160412'))&(obj.cStimN > 40)
                    initialFrameNum = initialFrameNum - 3;
                    % drop frame correction
                end
                
                if obj.segFrame == 0
                    obj.cStimN = obj.cStimN + 1;
                    % Usually I have the first recording without
                    % stimulation, but still want to get the calcium
                    % signals, so just simply get the stimulation frames
                    % from the next sweep
                end
                
                if obj.cStimN <= length(obj.patch.sweepFrames)
                segFrame = ...
                    initialFrameNum:(initialFrameNum + obj.patch.sweepFrames(obj.cStimN) - 1);
                else
                    DFFSegSig = [];
                    return;
                end
                
                if (strcmp(obj.cDate, '20160412'))&(obj.cStimN == 40)
                    segFrame((end - 2):end) = [];
                    % drop frame correction
                end
                % rawSegSig = obj.rawBlCorrSig(segFrame);
                
                if obj.cStimN == 72
                obj.cStimN
                end
                
                if segFrame(end) > length(obj.bgSubBlCorrSig)
                    DFFSegSig = [];
                    return;
                end
                bgSubBlCorrSig = obj.bgSubBlCorrSig(segFrame);
                
                blFrame = obj.blSegFrame;
                
                bsF = median(bgSubBlCorrSig(blFrame));
                DF = bgSubBlCorrSig - bsF;
                DFFSegSig = DF/bsF*100;
                
            end
        end
        
        function segFrameTime = get.segFrameTime(obj)
            segFrame = obj.segFrame;
            segFrameTime = obj.frameTime(segFrame);
        end
        
        function segFrame = get.segFrame(obj)
            if obj.patch.ndim == 2
                dataPt = [obj.patch.STSStartPt(obj.cStimN) obj.patch.STSEndPt(obj.cStimN)];
            else
                dataPt = [obj.patch.extStimStartPt(obj.cStimN) obj.patch.extStimEndPt(obj.cStimN)];
            end
            startFrame = floor((dataPt(1)/obj.patch.SR - obj.frameTime(1))*obj.fps);
            endFrame = ceil((dataPt(end)/obj.patch.SR - obj.frameTime(1))*obj.fps); 
            % Need to check when is the first frame time
            
            if startFrame < 0
                startFrame = 1;
            end
            segFrame = startFrame:endFrame;
        end

       
        function stimStartFrame = get.stimStartFrame(obj)
            
            if obj.patch.ndim == 2
                STSStartFrame = round(obj.patch.STSStartPt(obj.cStimN)/obj.patch.SR*obj.fps);
                
                stimFrameTemp =...
                    find((obj.patch.extStimStartPt > obj.patch.STSStartPt(obj.cStimN))&(obj.patch.extStimStartPt < obj.patch.STSEndPt(obj.cStimN)));
                stimStartFrameTemp = floor(obj.patch.extStimStartPt(stimFrameTemp(1))/obj.patch.SR*obj.fps);
                % The first trigger of extStimStartPt is to start
                % camera imaging, not for indcation of stim time
                
                stimStartFrame = stimStartFrameTemp - STSStartFrame;
            else
               
%                 stimStartFrame = floor(obj.patch.extStimStartPt(obj.cStimN)/obj.patch.SR*obj.fps);        

                if length(obj.frameTime) < obj.cStimN
                    stimStartFrame = [];
                    return
                   
                end
                if obj.patch.extStimStartPt(obj.cStimN) > obj.frameTime{obj.cStimN}(end)*obj.patch.SR
                    stimStartFrame = [];
                    return
                   
                end
                
                testDiff = obj.patch.extStimStartPt(obj.cStimN) - obj.frameTime{obj.cStimN}*obj.patch.SR;
                
                testDiff(testDiff > 0 ) = 1000000;
                
                [~, stimStartFrame] = min(abs(testDiff));
                
                
                % The first trigger of extStimStartPt is to start
                % camera imaging, not for indcation of stim time
                
            end
            
        end

        
        function blSegFrame = get.blSegFrame(obj)
            stimStartFrame = obj.stimStartFrame;
            blSegFrame = 1:stimStartFrame;
            
        end
       
        
        function fitData = fitSig(obj, smFac)
            % Calvin method, with other things trimmed
           
            % modelling only the stimulation onset and offset time, not
            % including the baseline time.
                        
%             obj.cROIN = ROIN;
%             obj.cStimN = stimN;
            DFFSegSigFull = obj.DFFSegSig;
            % with baseline
            stimStartFrame = obj.stimStartFrame;
            
%             blMean = obj.DFFsegBlMean(ROIN, stimN);
            
            DFFSegSig = DFFSegSigFull(stimStartFrame:end);
            
            if nargin < 5
                smoothData = smooth(DFFSegSig, 5)';
            else
                smoothData = smooth(DFFSegSig, smFac)';
            end
            
            vars = {'A', 't_0', '\tau_O_N', '\tau_O_F_F', 'c'};
            
            caModel = @(x,p) (p(5) + p(1) .* (x > p(2)) .* (1-exp( -(x - p(2))/p(3) )) .* ...
                exp( -(x - p(2))/p(4) ));
            
%             if blMean > -5
                P0 = [25 0 0.4 0.8 -1];
                % This is the original values used until 21/11/2014
%             else
%                 P0 = [blMean 0 0.4 0.8 -1];
                % 100uA fitting for 20141010 Clampex 20 couldn't be
                % fitted well with the original P(0) = 25, but fine
                % with blMean
%             end
            
            %                 P0 = [blMean 0 0.4 0.8 -1];
            O = optimset('MaxFunEvals',1E5,'MaxIter',1E5);
            
            %             P = zeros(length(W)-1, length(P0));
            %             V = zeros(length(W)-1, 1);
            segFrameTime = obj.segFrameTime;
            T = segFrameTime(stimStartFrame:end);
            T = T - T(1);
%             T = T(stimStartFrame:end);
            % getting rid of the time offset
            
            
            fitFcn = @(p) sum((smoothData - caModel(T, p)).^2);
            
            [obj.segFitInfo{obj.cROIN,obj.cStimN},~,~] = fminsearch(fitFcn,P0,O);
            
            fitData = caModel(T, obj.segFitInfo{obj.cROIN,obj.cStimN});
            
            obj.fitR2(obj.cROIN,obj.cStimN) = 1 - (var(smoothData -...
                caModel(T,obj.segFitInfo{obj.cROIN,obj.cStimN}))/ ...
                var(smoothData));
            % R square value (1 - SSres/SStot where SSres is sum(yi -
            % fi) squared and SStot is sum(yi - ymean) squared
            
            
            % RMS of the fit function minus the experimental values.
            
        end
        
        function calcAvgSigProperties(obj, varargin)
            
            
            ROIN = 1:obj.nROI;
            % If ROI number is not specified, plot all ROI signals
            
            
            % This parameters are for HFS recordings
            stimRep = 3;
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'ROIN');
                    ROIN = varargin{i + 1};
                end
                
                if strcmp(varargin{i}, 'stimRep');
                    stimRep = varargin{i + 1};
                end
            end
            
            nStimAmp = (obj.patch.nSweeps - 1)/stimRep;
            
            if obj.patch.ndim == 2
                % Gap-free recording mode
                obj.avgSigPeakVal = zeros(obj.nROI, 1);
                obj.avgSigHDT  = zeros(obj.nROI, 1);
                obj.avgSigPeakFrame = zeros(obj.nROI, 1);
                
                obj.cStimN = 1;
                
                for i = 1:obj.nROI
                    obj.cROIN = i;
%                     if i == 76
%                         i
%                     end
%                     [avgCoupledSig] = obj.getAvgCoupledSig;
                    [avgSegSig] = obj.getAvgSegSig;
                    stimStartFrame = obj.stimStartFrame;
                    
%                     avgStimSig = avgCoupledSig(stimStartFrame:end);
                    avgStimSig = avgSegSig(stimStartFrame:end);
                    
                    blFrame = obj.blSegFrame;
%                     blMean = mean(avgCoupledSig(blFrame));
                    blMean = mean(avgSegSig(blFrame));
%                     blSTD = std(avgCoupledSig(blFrame));
                    blSTD = std(avgSegSig(blFrame));
                    
                    [peakVal, peakFrame] = findpeaks(avgStimSig, 'MINPEAKHEIGHT', 3*blSTD,'MINPEAKDISTANCE',round(0.5*obj.fps));
                    if ~isempty(peakFrame)
                        latency = peakFrame(1)/obj.fps;
                    end
                    
                    if peakVal < 5
                        latencyCheckVal = 0.5;
                    else
                        latencyCheckVal = 1.5;
                        % For strong peak, the latency may be prolonged (or
                        % because of summation). So make the latency check
                        % value longer
                    end
                    
                    if ~isempty(peakVal)&(latency < latencyCheckVal)
                        % The signal is only considered as stimulation
                        % coupled if it has more than 3 times the baseline
                        % STD value and less than 500ms (0.5) latency
                        
                        
                        peakCaAmp = peakVal(1) - blMean;
                        % blMean is the mean value of the baseline frames obtained
                        % from the data (not from the fit data)
                        
                        bsReturn = mean(avgStimSig((end - 10):end));
                        
                        % Get the returned baseline value by averaging the values
                        % of the last 10 frames
                        
                        peakFall = peakVal(1) - bsReturn;
                        halfPeakFall = 0.5*peakFall + bsReturn;
                        
                        %             stimStartFrame = obj.getStimStartFrame(stimN);
                        
                        belowHalfPeakFall = find(avgStimSig(peakFrame(1):end) < halfPeakFall);
                        if ~isempty(belowHalfPeakFall)
                            halfDecayTime = belowHalfPeakFall(1)/obj.fps;
                            % belowHalfPeakFall(1) is the number of frames required
                            % to fall below half the (peak - returned baseline). So
                            % the decay time is the number of frames required
                            % divided by the sampling time.
                        else
                            halfDecayTime = NaN;
                        end
                        
                        obj.avgSigPeakVal(i) = peakCaAmp;
                        obj.avgSigHDT(i) = halfDecayTime;
                        obj.avgSigPeakFrame(i) = peakFrame(1);
                    else
                        obj.avgSigPeakVal(i) = 0;
                        obj.avgSigHDT(i) = 0;
                        obj.avgSigPeakFrame(i) = 0;
                    end
                end
                
            else
                obj.avgSigPeakVal = zeros(obj.nROI, nStimAmp);
                obj.avgSigHDT  = zeros(obj.nROI, nStimAmp);
                obj.avgSigPeakFrame = zeros(obj.nROI, nStimAmp);
                
                if isempty(obj.avgSig)
                    obj.getAvgSig;
                elseif isempty(obj.avgSig{1,1})
                    obj.getAvgSig;
                end
                
                for ROII = 1:length(ROIN)
                    
                    obj.cROIN = ROIN(ROII);
                    
                    for stimAmpI = 1:size(obj.avgSig,2)
                        
                        obj.cStimN = stimAmpI;
                        
                        if isempty(obj.avgSig{ROIN(ROII),stimAmpI})
                            continue
                        end
                        
%                         if stimAmpI == 22
%                             stimAmpI
%                         end
                        
                        avgSig = obj.avgSig{ROIN(ROII), stimAmpI};
                        stimStartFrame = obj.stimStartFrame;
                        if isempty(stimStartFrame)
                            continue;
                        end
                        
                        avgStimSig = avgSig(stimStartFrame:end);
                        
                        %                         blFrame = avgSig(1:stimStartFrame);
                        blMean = mean(avgSig(1:stimStartFrame));
                        blSTD = std(avgSig(1:stimStartFrame));
                        
                        caSigThresh = 3*blSTD/sqrt(stimRep);
                        % To be considered as Ca signal, the signal has to 
                        % be greater than 3 times the standard error 
                        % (standard error = standard deviation/sqrt(repetitions))
                        
                        [peakVal, peakFrame] = findpeaks(avgStimSig, 'MINPEAKHEIGHT', caSigThresh,'MINPEAKDISTANCE',round(0.5*obj.fps));
                        
                        [~, maxI] = max(peakVal);
                        
                        if ~isempty(peakFrame)
                            latency = peakFrame(maxI)/obj.fps;
                        end
                        
                        if peakVal < 5
                            latencyCheckVal = 0.5;
                        else
                            latencyCheckVal = 1.5;
                            % For strong peak, the latency may be prolonged (or
                            % because of summation). So make the latency check
                            % value longer
                        end
                        
                        if ~isempty(peakVal)&(latency < latencyCheckVal)
                            % The signal is only considered as stimulation
                            % coupled if it has more than 3 times the baseline
                            % STD value and less than 500ms (0.5) latency
                            
                            
                            peakCaAmp = peakVal(maxI) - blMean;
                            % blMean is the mean value of the baseline frames obtained
                            % from the data (not from the fit data)
                            
                            bsReturn = mean(avgStimSig((end - 10):end));
                            
                            % Get the returned baseline value by averaging the values
                            % of the last 10 frames
                            
                            peakFall = peakVal(maxI) - bsReturn;
                            halfPeakFall = 0.5*peakFall + bsReturn;
                            
                            %             stimStartFrame = obj.getStimStartFrame(stimN);
                            
                            belowHalfPeakFall = find(avgStimSig(peakFrame(maxI):end) < halfPeakFall);
                            if ~isempty(belowHalfPeakFall)
                                halfDecayTime = belowHalfPeakFall(maxI)/obj.fps;
                                % belowHalfPeakFall(1) is the number of frames required
                                % to fall below half the (peak - returned baseline). So
                                % the decay time is the number of frames required
                                % divided by the sampling time.
                            else
                                halfDecayTime = NaN;
                            end
                            
                            obj.avgSigPeakVal(ROIN(ROII), stimAmpI) = peakCaAmp;
                            obj.avgSigHDT(ROIN(ROII), stimAmpI) = halfDecayTime;
                            obj.avgSigPeakFrame(ROIN(ROII), stimAmpI) = peakFrame(maxI);
                        else
                            obj.avgSigPeakVal(ROIN(ROII), stimAmpI) = 0;
                            obj.avgSigHDT(ROIN(ROII), stimAmpI) = 0;
                            obj.avgSigPeakFrame(ROIN(ROII), stimAmpI) = 0;
                        end
                    end
                end
            end
        end

        
        function calcSegSigProperties(obj,dataType,smFac,smMethod)
            % Calculate the properties based on each stimulation segment
%             function calcSigProperties(obj,dataType,ROIN,stimN, smFac, smMethod)
            obj.segSigPeakVal = zeros(obj.nROI, obj.nStim);
            obj.segSigHDT  = zeros(obj.nROI, obj.nStim);
            obj.segSigPeakFrame = zeros(obj.nROI, obj.nStim);
            obj.nCoupledSig = zeros(obj.nROI, 1);
            
%             if length(ROIN) > obj.nROI
%                 ROIN((obj.nROI + 1):end) = [];
%                 % If ROIN is more than the number of ROIs, then get whole
%                 % ROI properties and clear the excess.
%             end
%                         
%             if length(stimN) > obj.nStim
%                 ROIN((obj.nStim + 1):end) = [];
%                 % If stimN is more than the number of stims, then get whole
%                 % stim properties and clear the excess.
%             end
            
            
            for i = 1:obj.nROI
                for j = 1:obj.nStim
                    
                    obj.cROIN = i;
                    obj.cStimN = j;
                    
                    stimStartFrame = obj.stimStartFrame;
                    
                    
                    DFFSegSig = obj.DFFSegSig;
                    DFFSegSigAfterStim = DFFSegSig(stimStartFrame:end);
                    
                    blFrame = obj.blSegFrame;
                    blMean = mean(DFFSegSig(blFrame));
                    blSTd = std(DFFSegSig(blFrame));
                    
                    if strcmp(dataType, 'DFF')
                        anaData = DFFSegSigAfterStim;
                    elseif strcmp(dataType, 'smDFF')
                        if nargin > 3
                            anaData = obj.getSmData(DFFSegSigAfterStim, smFac, smMethod);
                        elseif nargin > 2
                            anaData = obj.getSmData(DFFSegSigAfterStim, smFac);
                        else
                            anaData = obj.getSmData(DFFSegSigAfterStim);
                        end
                    elseif strcmp(dataType, 'fitDFF')
                        if nargin > 4
                            anaData = obj.fitSig(smFac);
                        else
                            anaData = obj.fitSig;
                        end
                    end
                    
                    [peakVal, peakFrame] = max(anaData);
                    latency = peakFrame/obj.fps;
                    
                    if (peakVal > 3*blSTd)&(latency < 0.5)
                        % The signal is only considered as stimulation
                        % coupled if it has more than 3 times the baseline
                        % STD value and less than 500ms (0.5) latency

                        
                        peakCaAmp = peakVal - blMean;
                        % blMean is the mean value of the baseline frames obtained
                        % from the data (not from the fit data)
                        
                        bsReturn = mean(anaData((end - 10):end));
                        
                        % Get the returned baseline value by averaging the values
                        % of the last 10 frames
                        
                        peakFall = peakVal - bsReturn;
                        halfPeakFall = 0.5*peakFall + bsReturn;
                        
                        %             stimStartFrame = obj.getStimStartFrame(stimN);
                        
                        belowHalfPeakFall = find(anaData(peakFrame:end) < halfPeakFall);
                        if ~isempty(belowHalfPeakFall)
                            halfDecayTime = belowHalfPeakFall(1)/obj.fps;
                            % belowHalfPeakFall(1) is the number of frames required
                            % to fall below half the (peak - returned baseline). So
                            % the decay time is the number of frames required
                            % divided by the sampling time.
                        else
                            halfDecayTime = NaN;
                        end
                        
                        obj.segSigPeakVal(i, j) = peakCaAmp;
                        obj.segSigHDT(i, j) = halfDecayTime;
                        obj.segSigPeakFrame(i, j) = peakFrame;
                    else
                        obj.segSigPeakVal(i, j) = 0;
                        obj.segSigHDT(i, j) = 0;
                        obj.segSigPeakFrame(i, j) = 0;
                        
%                         obj.segSigPeakVal(ROIN(i), stimN(j)) = peakCaAmp;
%                         obj.segSigHDT(ROIN(i), stimN(j)) = halfDecayTime;
%                         obj.segSigPeakFrame(ROIN(i), stimN(j)) = peakFrame;
%                     else
%                         obj.segSigPeakVal(ROIN(i), stimN(j)) = 0;
%                         obj.segSigHDT(ROIN(i), stimN(j)) = 0;
%                         obj.segSigPeakFrame(ROIN(i), stimN(j)) = 0;
                    
                    end
                end
                obj.nCoupledSig(i) = sum(obj.segSigPeakVal(i, :) > 0);
            end
            
        end
        
        function [avgCoupledSig, wholeCoupledSig, stimCoupledI] = getAvgCoupledSig(obj)
            % avgCoupledSig = mean over stimulation coupled signals
            % wholeCoupledSig = sum of stimulation coupled signals
            % stimCoupledI = stimulation coupled segments
            minFrame = round(obj.patch.STSTotalPt/obj.patch.SR*obj.fps);
            wholeCoupledSig = zeros(obj.nStim, minFrame);
            
            stimCoupledI = [];
            
            for i = 1:obj.nStim
                
                obj.cStimN = i;                
                DFFSegSig = obj.DFFSegSig;
                
                stimStartFrame = obj.stimStartFrame;
                anaData = DFFSegSig(stimStartFrame:end);
                
                blFrame = obj.blSegFrame;
                blMean = mean(DFFSegSig(blFrame));
                blSTD = std(DFFSegSig(blFrame));
                
                [peakVal, peakFrame] = findpeaks(anaData, 'MINPEAKHEIGHT', 3*blSTD,'MINPEAKDISTANCE',round(0.5*obj.fps));
                
                if ~isempty(peakFrame)
                    latency = peakFrame(1)/obj.fps;
                else 
                    latency = 100;
                end
                
                if isempty(peakVal)
                    peakVal = 0;
                end

                
                if (peakVal(1) > 3*blSTD)&(latency < 0.6)
                    stimCoupledI = [stimCoupledI, i];
                    wholeCoupledSig(i,:) = DFFSegSig(1:minFrame);
                end                
            end
            avgCoupledSig = mean(wholeCoupledSig,1);
        end
        
        function avgCoupledSigSTD = getAvgCoupledSigSTD(obj)
            [~, wholeCoupledSig] = obj.getAvgCoupledSig;
            avgCoupledSigSTD = std(wholeCoupledSig,1,1);            
        end
        
        function avgSegSig = getAvgSegSig(obj)
            % Average every point of each stimulation section to give an
            % averaged segment signal
%             avgTime = 3;
            minFrame = round(obj.patch.STSTotalPt/obj.patch.SR*obj.fps);
            
            avgSegSig = zeros(obj.nStim, minFrame);
            
                for i = 1:obj.nStim
                    
                    obj.cStimN = i;
                    
%                     stimStartFrame = obj.stimStartFrame;
                    
                    DFFSegSig = obj.DFFSegSig;
                    
                    if minFrame <= length(DFFSegSig)
                        % 20150707 recording 9 to 19 has the first stim
                        % segment started at the wrong time, so will lack 1
                        % frame
                        avgSegSig(i,:) = DFFSegSig(1:minFrame);
                    else
                        DFFSegSigTemp = DFFSegSig;
                        DFFSegSig = zeros(1, minFrame);
                        nDiffFrame = minFrame - length(DFFSegSigTemp);
                        startFrame = nDiffFrame + 1;
                        avgSegSig(i, startFrame:end) = DFFSegSigTemp;
                    end
                    
                end
  
             avgSegSig = mean(avgSegSig,1);
        end
        
        
        function avgSegSTD = getAvgSegSTD(obj)
%             avgTime = 3;
%             avgFrame = round(avgTime*obj.fps);
            minFrame = round(obj.patch.STSTotalPt/obj.patch.SR*obj.fps);
            
%             avgSegSTD = zeros(obj.nROI, avgFrame);
            
            tempSegSig = zeros(obj.nStim, minFrame);
            

            
            for i = 1:obj.nStim
                obj.cStimN = i;
                
%                 stimStartFrame = obj.stimStartFrame;
                

                DFFSegSig = obj.DFFSegSig;
                
                if minFrame <= length(DFFSegSig)
                    % 20150707 recording 9 to 19 has the first stim
                    % segment started at the wrong time, so will lack 1
                    % frame
                    tempSegSig(i,:) = DFFSegSig(1:minFrame);
                else
                    DFFSegSigTemp = DFFSegSig;
                    DFFSegSig = zeros(1, minFrame);
                    nDiffFrame = minFrame - length(DFFSegSigTemp);
                    startFrame = nDiffFrame + 1;
                    tempSegSig(i, startFrame:end) = DFFSegSigTemp;
                end
            end
            
            
            
            avgSegSTD = std(tempSegSig,1,1);            
        end
     

        
        function fitExpRiseSig(obj,ROINo, stimNo, smMethod, smFac)
            %             obj.anaCaSig(smMethod, smFac);
            data = obj.riseSig{ROINo, stimNo};
            data = obj.maxAmp(ROINo, stimNo) - data;
            dataTime = obj.frameTime(obj.riseFrame{ROINo, stimNo});
            dataTime = dataTime - dataTime(1);
            %             [fitresult, gof] = obj.fitExp2(dataTime,data);
            [fitresult, gof] = obj.fitExp(dataTime,data');
            obj.riseFitExpCoef{ROINo, stimNo} = fitresult;
            obj.riseFitExpGOF{ROINo, stimNo} = gof;
        end
        
        function fitExp2RiseSig(obj,ROINo, stimNo, smMethod, smFac)
            %             obj.anaCaSig(smMethod, smFac);
            data = obj.riseSig{ROINo, stimNo};
            data = obj.maxAmp(ROINo, stimNo) - data;
            dataTime = obj.frameTime(obj.riseFrame{ROINo, stimNo});
            dataTime = dataTime - dataTime(1);
            %             [fitresult, gof] = obj.fitExp2(dataTime,data);
            [fitresult, gof] = obj.fitExp2(dataTime,data');
            obj.riseFitExpCoef{ROINo, stimNo} = fitresult;
            obj.riseFitExpGOF{ROINo, stimNo} = gof;
        end
        
        function fitExp2RDSig(obj,ROINo, stimNo, smMethod, smFac)
            %             obj.anaCaSig(smMethod, smFac);
            data = [obj.riseSig{ROINo, stimNo}', obj.decaySig{ROINo, stimNo}'];
            %             data = data/obj.max(data);
            dataTime = [obj.frameTime(obj.riseFrame{ROINo, stimNo}), obj.frameTime(obj.decayFrame{ROINo, stimNo})];
            dataTime = dataTime - dataTime(1);
            %             [fitresult, gof] = obj.fitExp2(dataTime,data);
            [fitresult, gof] = obj.fitExp2(dataTime,data);
            obj.RDFitExpCoef{ROINo, stimNo} = fitresult;
            obj.RDFitExpGOF{ROINo, stimNo} = gof;
        end
        
        function fitLinRiseSig(obj,ROINo, stimNo, smMethod, smFac)
            %             obj.anaCaSig(smMethod, smFac);
            data = obj.riseSig{ROINo, stimNo};
            dataTime = obj.frameTime(obj.riseFrame{ROINo, stimNo});
            dataTime = dataTime - dataTime(1);
            %             [fitresult, gof] = obj.fitExp2(dataTime,data);
            obj.riseFitLinCoef{ROINo, stimNo} = polyfit(dataTime,data',1);
            
        end
        
        
        function fitExpDecaySig(obj,ROINo, stimNo, smMethod, smFac)
            %             obj.anaCaSig(smMethod, smFac);
            data = obj.decaySig{ROINo, stimNo};
            %             data = data/obj.max(data);
            dataTime = obj.frameTime(obj.decayFrame{ROINo, stimNo});
            dataTime = dataTime - dataTime(1);
            %             [fitresult, gof] = obj.fitExp2(dataTime,data);
            [fitresult, gof] = obj.fitExp(dataTime,data');
            obj.decayFitExpCoef{ROINo, stimNo} = fitresult;
            obj.decayFitExpGOF{ROINo, stimNo} = gof;
        end
        
    
        
        
        function fitDFFSeg(obj,smMethod, smFac)
            obj.anaCaSig(smMethod, smFac);
            for ROINo = 1:obj.nROI
                for stimNo = 1:obj.nStim
                    if obj.stimCoupledRes(ROINo, stimNo) == 1
                        obj.fitCalvinMethod(ROINo, stimNo);
                    end
                end
            end
        end
        
        function [fitresult, gof] = fitExp2(obj,dataTime,data)
            [xData, yData] = prepareCurveData(dataTime, data);
            
            % Set up fittype and options.
            ft = fittype('exp2');
            opts = fitoptions('Method', 'NonlinearLeastSquares');
            opts.Display = 'Off';
            opts.StartPoint = [dataTime(1) dataTime(2) data(1) data(2)];
            
            % Fit model to data.
            [fitresult, gof] = fit(xData, yData, ft, opts);
            
            % Plot fit with data.
            %             figure('Name','untitled fit 1' );
            %             h = plot(fitresult, xData, yData);
            %             legend(h, 'data vs. dataTime', 'untitled fit 1', 'Location', 'NorthEast' );
            %             % Label axes
            %             xlabel('dataTime' );
            %             ylabel('data');
            %             grid on
        end
        
        function [fitresult, gof] = fitExp(obj,dataTime,data)
            [xData, yData] = prepareCurveData(dataTime, data);
            
            % Set up fittype and options.
            ft = fittype('exp1');
            opts = fitoptions('Method', 'NonlinearLeastSquares');
            opts.Display = 'Off';
            opts.StartPoint = [dataTime(1) dataTime(2)];
            %             opts.StartPoint = [1.29001320957929 2.01231412434391 -2603802574908.06 -25.8726883356248];
            
            % Fit model to data.
            [fitresult, gof] = fit(xData, yData, ft, opts);
            
            % Plot fit with data.
            %             figure('Name','untitled fit 1' );
            %             h = plot(fitresult, xData, yData);
            %             legend(h, 'data vs. dataTime', 'untitled fit 1', 'Location', 'NorthEast' );
            %             % Label axes
            %             xlabel('dataTime' );
            %             ylabel('data');
            %             grid on
        end
        
        
        function nStim = get.nStim(obj)
            nStim = length(obj.patch.STSStartPt);
            
            
        end
        
        function BW = circ(obj, BW, hCentre,wCentre,r)
            
            th = 0:pi/50:2*pi;
            
            hUnit = round(r * sin(th) + hCentre);
            
            wUnit = round(r * cos(th) + wCentre);
            
%             BW = zeros(obj.height, obj.width);
            
            for i = 1:length(th)
               BW(hUnit(i),wUnit(i)) = 1; 
            end
            
        end
        
        
%         function coupledSigProb = get.coupledSigProb(obj)
%             if isempty(obj.nCoupledSig)
%                 %                display('Number of coupled signals have not been calculated yet');
%                 coupledSigProb = [];
%             else
%                 coupledSigProb = obj.nCoupledSig/obj.nStim;
%             end
%         end
        
        function set.cROIN(obj,ROIN)
            obj.cROIN = ROIN;
            
        end
        
        function set.cStimN(obj,stimN)
            obj.cStimN = stimN;
        end
        
        function stimN = get.cStimN(obj)
            stimN = obj.cStimN;            
        end
        
        function ROIN = get.cROIN(obj)
            ROIN = obj.cROIN;
            
        end
        
        function nROI = get.nROI(obj)
           nROI = size(obj.sigROI, 2); 
            
        end
        
    end
    
end

function [TraceGradData,DockingData] =...
    Det_and_Apply_Threshold_Filters(FigureHandles,UniversalData,Options,TraceRunMedian,CurrTraceCropped)

% Here we apply the three tests (gradient test, difference test, spike test)
% to our current trace and identify frame numbers which test positive in 
% the various tests. This information is then relayed back (in TraceGradData)
% and analyzed later on to classify each trace as a fusion event, etc.

% Define and calculate some variables and values before we start

    ProbableDockingEvent = 'n';
        % The trigger for a probable docking event is set to y only if the type
        % of fusion data is correct and if the event meets certain criteria
    MaxIntensity = max(TraceRunMedian.Trace);
    
    %Calc the gradient of the running median
        GradTraceRunMed = gradient(TraceRunMedian.Trace);
         
    %Set all grad trace before Binding to zero (don't see events
    %before Binding).
    if strcmp(Options.NoFusionAllowedBeforepHDrop,'y')
        BindFrameNumIdx = find(TraceRunMedian.FrameNumbers==UniversalData.StandardBindFrameNum);
        GradTraceRunMed(1:BindFrameNumIdx) = 0;
    end

    % Set grad trace near the end to zero (otherwise occasionally get 
    % funny edge effects for a small percentage of events where the 
    % last frame or two is erroneously identified as a fusion event)
    if strcmp(Options.ClipFramesAtEnd,'y')
        GradTraceRunMed(end-Options.NumFramesToClip:end) = 0;
    end
    
% -------------------- Gradient Test -------------------- 

% Define gradient filters (i.e. thresholds), etc. which will be used to
% identify potential fusion events by the gradient tests (positive and
% negative)

    if strcmp(Options.TypeofFusionData, 'TetheredVesicle')
        STDFilterFactor = Options.GradientThreshold;
        STDFilterFactorDifferencePos = Options.DifferenceTraceThresholdPos;
        STDFilterFactorDifferenceNeg = Options.DifferenceTraceThresholdNeg;
        PositiveGradFilter = STDFilterFactor*std(GradTraceRunMed);
        RangeToFilterPositive = Options.NumFramesBetweenGradientEvents;
        RangeToFilterNegative = Options.NumFramesBetweenGradientEvents;
        VirusStopbyIndex = find(TraceRunMedian.FrameNumbers==...
            max([UniversalData.StandardBindFrameNum + 1,Options.FrameToStartAnalysis + Options.RunMedHalfLength,...
            min(TraceRunMedian.FrameNumbers)]));
                %Added 250731 - added min(TraceRunMedian.FrameNumbers) to
                %account for when we want the frame to start the analysis
                %to be greater than 1. Currently it is hardcoded at 1. So
                %we choose instead the first frame number in the running median of
                %the trace.
            % This determines where we will start analyzing for fusion events
            if isempty(VirusStopbyIndex)
                ErrorinCalculatingVirusStopbyIndex;
            end
    elseif strcmp(Options.TypeofFusionData, 'SLBSelfQuench')
        % This is a special case in which we need to decide whether there
        % is a docking event at the beginning. If so, then we only look for
        % fusion events in the rest of the trace. To detect whether there
        % is a docking event, we will look for a very large intensity jump
        % near the beginning of the trace which remains high for the
        % remainder
        STDFilterFactor = Options.GradientThreshold;
        STDFilterFactorDifference = Options.DifferenceTraceThreshold;
        VirusStopbyIndex = find(TraceRunMedian.FrameNumbers==UniversalData.FrameAllVirusStoppedBy);
        TruncatedGradTraceDock = GradTraceRunMed(BindFrameNumIdx:VirusStopbyIndex);
        FrameNumbersDockTrace =TraceRunMedian.FrameNumbers(BindFrameNumIdx:VirusStopbyIndex);
        MaxGradientBeforeStop = max(TruncatedGradTraceDock);
        EarlyMedianIntensity = median(TraceRunMedian.Trace(1:UniversalData.StandardBindFrameNum));
        RangeToFilterDock = 5;
        DockGradientFilter = MaxGradientBeforeStop - 0.5*MaxGradientBeforeStop;
        DockFilteredRunMed =  TruncatedGradTraceDock > DockGradientFilter;
        
        if MaxGradientBeforeStop > 0.15*MaxIntensity && EarlyMedianIntensity < 1500
            % Probably There Is a Docking Event
            ProbableDockingEvent = 'y';    
            PositiveGradFilter = STDFilterFactor*std(GradTraceRunMed(VirusStopbyIndex:end));
        else
            PositiveGradFilter = STDFilterFactor*std(GradTraceRunMed(BindFrameNumIdx:end));
        end
        RangeToFilterPositive = Options.NumFramesBetweenGradientEvents;
        RangeToFilterNegative = Options.NumFramesBetweenGradientEvents;
        
    end

    % Now we define the negative gradient filter (which is the same regardless 
    % of data type, but really only relevant for fusion to supported bilayers).
        
        NegativeGradFilter = -STDFilterFactor*std(GradTraceRunMed(VirusStopbyIndex:end));
    
    % Apply gradient filters
        NegFilteredGradTrace = GradTraceRunMed < NegativeGradFilter;
        PosFilteredGradTrace = GradTraceRunMed > PositiveGradFilter;
            
    % Filter out adjacent test results that are higher than the cutoff (a
    % single event can sometimes be above the cutoff for several frames on
    % either side, and we don't want to identify each of those frames as separate events)
        FilterOption='First Only';
        [NegFilteredGradTrace] = Filter_Adjacent_Gradient_Results(RangeToFilterNegative,NegFilteredGradTrace,FilterOption,GradTraceRunMed);
        [PosFilteredGradTrace] = Filter_Adjacent_Gradient_Results(RangeToFilterPositive,PosFilteredGradTrace,FilterOption,GradTraceRunMed);
        if strcmp(ProbableDockingEvent, 'y')
            FilterOption='Max Then First';
            [DockFilteredRunMed] = Filter_Adjacent_Gradient_Results(RangeToFilterDock,DockFilteredRunMed,FilterOption,TruncatedGradTraceDock);
        end
        
% -------------------- Difference Test --------------------
    % Calculate difference trace (difference between current point and the
    % nth previous point). This will be used to identify either a slow
    % increase in fluorescence or a slow decrease in fluorescence which may
    % indicate an unusual fusion event.
        NumberPointsBack = Options.NumberFramesBackToSubtract;
        DifferenceTrace = zeros(length(TraceRunMedian.Trace),1);
        DifferenceTrace(VirusStopbyIndex+NumberPointsBack+1:end) =  TraceRunMedian.Trace(VirusStopbyIndex+NumberPointsBack+1:end) -...
            TraceRunMedian.Trace(VirusStopbyIndex+1:end-NumberPointsBack);

        DifferenceFilterPos = STDFilterFactorDifferencePos* std(GradTraceRunMed(VirusStopbyIndex:end));
        DifferenceFilterNeg = -STDFilterFactorDifferenceNeg* std(GradTraceRunMed(VirusStopbyIndex:end));
        
        FilteredDiffTracePos = DifferenceTrace > DifferenceFilterPos;
        FilteredDiffTraceNeg = DifferenceTrace < DifferenceFilterNeg;
        RangeToFilterDifference = NumberPointsBack;

    % Analyze the filtered difference traces to pull out metrics such 
    % as number of clusters of frames above the cutoff value, size 
    % of frame clusters, etc.
        [DiffNegClusterData] = Analyze_Trace_Clusters(FilteredDiffTraceNeg,5);
        [DiffPosClusterData] = Analyze_Trace_Clusters(FilteredDiffTracePos,5);
    
% -------------------- Spike Test --------------------
    % Calculate difference between raw data values and a running median
    % that is of a wider range, this is to identify transient spikes in the
    % data which might indicate a self quenched fusion event that returns to a similar
    % intensity value, or something bright that flies by quickly out of focus
        RunMedHalfLength = Options.RunMedHalfLengthSpike;

%         StartIdx = RunMedHalfLength + 1;
%         EndIdx = length(TraceRunMedian.Trace)-RunMedHalfLength;

        StartIdx = Options.FrameToStartAnalysis;
        EndIdx = length(TraceRunMedian.Trace);

        TraceRunMedianWiderRange = zeros(length(StartIdx:EndIdx),1);
        %We also set up a vector with the actual frame numbers
        %corresponding to each index position.
        SpikeFrameNumbers = TraceRunMedian.FrameNumbers(StartIdx:EndIdx);

%         for n = StartIdx:EndIdx
%             TraceRunMedianWiderRange(n-RunMedHalfLength) = median(TraceRunMedian.Trace(n-RunMedHalfLength:n+RunMedHalfLength));
%         end

        
        TraceLength = length(TraceRunMedian.FrameNumbers);
        OldTrace = TraceRunMedian.Trace;
        TraceRunMedianWiderRange = zeros(size(TraceRunMedian.FrameNumbers));
    
        for n = 1:TraceLength
            if n - RunMedHalfLength < 1
                TraceRunMedianWiderRange(n) = median(OldTrace(1:n+RunMedHalfLength));
            elseif n + RunMedHalfLength > TraceLength
                TraceRunMedianWiderRange(n) = median(OldTrace(n-RunMedHalfLength:end));
            else
                TraceRunMedianWiderRange(n) = median(OldTrace(n-RunMedHalfLength:n+RunMedHalfLength));
            end
            
        end

        SpikeTrace = TraceRunMedian.Trace - TraceRunMedianWiderRange; 
        STDFilterFactorSpike = Options.SpikeThreshold;
        SpikeFilter = STDFilterFactorSpike* std(SpikeTrace(VirusStopbyIndex:end));
        FilteredSpikeTrace = SpikeTrace > SpikeFilter;
        RangeToFilterSpike = Options.NumFramesBetweenSpikeEvents;
        
    % Filter out adjacent test results that are higher than the cutoff (a
    % single event can sometimes be above the cutoff for several frames on
    % either side, and we don't want to identify each of those frames as separate events)
        FilterOption='First Only';
        [FilteredSpikeTrace] = Filter_Adjacent_Gradient_Results(...
            RangeToFilterSpike,FilteredSpikeTrace,FilterOption,[]);
    
% -------------------- Plot --------------------
    %Plot gradient & filters
    if strcmp(Options.DisplayFigures,'y')
        set(0,'CurrentFigure',FigureHandles.GradientWindow)
        cla
        plot(TraceRunMedian.FrameNumbers,GradTraceRunMed,'r-');
        hold on
%         plot(TraceRunMedian.FrameNumbers,DifferenceTrace ,'g--');
%         plot(SpikeFrameNumbers,SpikeTrace,'m--');
        plot(TraceRunMedian.FrameNumbers,ones(1,length(TraceRunMedian.FrameNumbers))*NegativeGradFilter,'g--');
        plot(TraceRunMedian.FrameNumbers,ones(1,length(TraceRunMedian.FrameNumbers))*PositiveGradFilter,'g--');
        if strcmp(ProbableDockingEvent, 'y')
            plot(TraceRunMedian.FrameNumbers,ones(1,length(TraceRunMedian.FrameNumbers))*DockGradientFilter,'b--');
        end
        TitleInfo = strcat('Gradient test',...
            '; Event =',num2str(UniversalData.TraceNumber),'/',num2str(UniversalData.NumTraces),...
            '; STD=',num2str(std(GradTraceRunMed(VirusStopbyIndex:end))));
        title(TitleInfo)
        xlabel('Frame Number');
        hold off
        Draw_Events_On_Plot(NegFilteredGradTrace,TraceRunMedian.FrameNumbers,FigureHandles.GradientWindow,'r-')
        Draw_Events_On_Plot(PosFilteredGradTrace,TraceRunMedian.FrameNumbers,FigureHandles.GradientWindow,'k-')
        
        set(0,'CurrentFigure',FigureHandles.DiagnosticWindow)
        cla
        hold on
        plot(TraceRunMedian.FrameNumbers,DifferenceTrace ,'g-');
        plot(TraceRunMedian.FrameNumbers,ones(1,length(TraceRunMedian.FrameNumbers))*DifferenceFilterNeg,'g--');
        plot(TraceRunMedian.FrameNumbers,ones(1,length(TraceRunMedian.FrameNumbers))*DifferenceFilterPos,'g--');
            
        plot(SpikeFrameNumbers,SpikeTrace,'m-');
        plot(SpikeFrameNumbers,ones(1,length(SpikeFrameNumbers))*SpikeFilter,'m--');
        title('Difference (green) and Spike (magenta) tests')
        xlabel('Frame Number');
        ylabel('Intensity Difference (AU)');
        hold off
        
        Draw_Events_On_Plot(FilteredDiffTraceNeg,TraceRunMedian.FrameNumbers,FigureHandles.DiagnosticWindow,'k:')
        Draw_Events_On_Plot(FilteredDiffTracePos,TraceRunMedian.FrameNumbers,FigureHandles.DiagnosticWindow,'b:')
        Draw_Events_On_Plot(FilteredSpikeTrace,SpikeFrameNumbers,FigureHandles.DiagnosticWindow,'r:')
 
    end
    

% -------------------- Data Management --------------------   
    if strcmp(Options.TypeofFusionData,'SLBSelfQuench')
        [DockingData]=Compile_Docking_Data(ProbableDockingEvent,DockFilteredRunMed,VirusStopbyIndex,...
            UniversalData.FrameAllVirusStoppedBy,FrameNumbersDockTrace,RangeToFilterDock);
    else 
        DockingData = [];
    end
    
    % Compile all of the trace and gradient data into one structure
        TraceGradData.CroppedRawTrace = CurrTraceCropped;
        TraceGradData.TraceRunMedian = TraceRunMedian;
        TraceGradData.GradTraceRunMed = GradTraceRunMed;
        TraceGradData.PositiveGradFilter = PositiveGradFilter;
        TraceGradData.NegativeGradFilter = NegativeGradFilter;
        TraceGradData.NegFilteredGradTrace = NegFilteredGradTrace;
        TraceGradData.PosFilteredGradTrace= PosFilteredGradTrace;
        TraceGradData.RangeToFilterPositive = RangeToFilterPositive;
        TraceGradData.RangeToFilterNegative = RangeToFilterNegative;

        TraceGradData.SpikeTrace = SpikeTrace;
        TraceGradData.SpikeFilter = SpikeFilter;
        TraceGradData.RangeToFilterSpike = RangeToFilterSpike;
        TraceGradData.FilteredSpikeTrace= FilteredSpikeTrace;
        TraceGradData.SpikeFrameNumbers = SpikeFrameNumbers;

        TraceGradData.DifferenceTrace = DifferenceTrace;
        TraceGradData.DifferenceFilterNeg = DifferenceFilterNeg;
        TraceGradData.DifferenceFilterPos = DifferenceFilterPos;
        TraceGradData.DiffPosClusterData = DiffPosClusterData;
        TraceGradData.DiffNegClusterData = DiffNegClusterData;
        TraceGradData.DiffTraceFrameNumbers = TraceRunMedian.FrameNumbers;
        TraceGradData.FilteredDiffTracePos = FilteredDiffTracePos;
        TraceGradData.FilteredDiffTraceNeg = FilteredDiffTraceNeg;
        TraceGradData.RangeToFilterDifference= RangeToFilterDifference;
end

function [DockingData]=Compile_Docking_Data(ProbableDockingEvent,DockFilteredRunMed,VirusStopbyIndex,VirusStopbyFrameNumber,...
    FrameNumbersDockTrace,RangeToFilterDock)
DockingData.VirusStopbyIndex= VirusStopbyIndex;
DockingData.VirusStopbyFrameNumber= VirusStopbyFrameNumber;
DockingData.ProbableDockingEvent = ProbableDockingEvent; 
DockingData.DockFilteredRunMed =DockFilteredRunMed;
DockingData.FrameNumbersDockTrace=FrameNumbersDockTrace;
DockingData.RangeToFilterDock= RangeToFilterDock;
end

function Draw_Events_On_Plot(FilteredTrace,FrameNumbers,FigureHandle,LineType)
    set(0,'CurrentFigure',FigureHandle)
    hold on
    NumberofEvents = sum(FilteredTrace);
    EventFrameNumbers = FrameNumbers(FilteredTrace);
    LineToPlot = ylim;
    for j= 1:NumberofEvents
        XToPlot = [EventFrameNumbers(j), EventFrameNumbers(j)];
        plot(XToPlot,LineToPlot,LineType)
    end
    hold off
end
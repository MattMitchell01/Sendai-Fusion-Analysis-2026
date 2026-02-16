function [TraceRunMedian,FigureHandles] = ...
    Run_Med_And_Plot(CurrTraceCropped,FigureHandles,UniversalData,Options,CurrentVirusData)
    
%Calc the running median, plot running median and where binding
%occurred
    RunMedHalfLength = Options.RunMedHalfLength;

        StartIdx = Options.FrameToStartAnalysis;
        EndIdx = length(CurrTraceCropped.Trace);
        %TraceRunMedian.Trace = zeros(length(StartIdx:EndIdx),1);
        
        %We also set up a vector with the actual frame numbers and time
        %vector corresponding to each index position.
        TraceRunMedian.FrameNumbers = CurrTraceCropped.FrameNumbers(StartIdx:EndIdx);
        TraceRunMedian.TimeVector = CurrTraceCropped.TimeVector(StartIdx:EndIdx);

%     for n = StartIdx:EndIdx
%         TraceRunMedian.Trace(n-RunMedHalfLength) = median(CurrTraceCropped.Trace(n-RunMedHalfLength:n+RunMedHalfLength));
%     end


    TraceLength = length(CurrTraceCropped.FrameNumbers);
    OldTrace = CurrTraceCropped.Trace;
    TraceRunMedian.Trace = zeros(size(CurrTraceCropped.FrameNumbers));

    for n = 1:TraceLength
        if n - RunMedHalfLength < 1
            TraceRunMedian.Trace(n) = median(OldTrace(1:n+RunMedHalfLength));
        elseif n + RunMedHalfLength > TraceLength
            TraceRunMedian.Trace(n) = median(OldTrace(n-RunMedHalfLength:end));
        else
            TraceRunMedian.Trace(n) = median(OldTrace(n-RunMedHalfLength:n+RunMedHalfLength));
        end
        
    end


%Plot running median
    if strcmp(Options.DisplayFigures,'y')
        set(0,'CurrentFigure',FigureHandles.TraceWindow)
        cla
        hold on
        % Plot green if a good virus, and red if a bad virus (only if we are analyzing bad viruses)
        if strcmp(CurrentVirusData.IsVirusGood,'y')
            % plot(CurrTraceCropped.FrameNumbers,CurrTraceCropped.Trace,'b-')
            
            plot(TraceRunMedian.FrameNumbers,TraceRunMedian.Trace,'g-');
        else
            plot(TraceRunMedian.FrameNumbers,TraceRunMedian.Trace,'r-');
        end
        %Plot a line where the bind time occurred
            set(0,'CurrentFigure',FigureHandles.TraceWindow)
            hold on
                LineToPlot = ylim;
                XToPlot = [UniversalData.StandardBindFrameNum, UniversalData.StandardBindFrameNum];
            plot(XToPlot,LineToPlot,'m--')
        xlabel('Frame Number');
        ylabel('Intensity (AU)');
        hold off
    end
end
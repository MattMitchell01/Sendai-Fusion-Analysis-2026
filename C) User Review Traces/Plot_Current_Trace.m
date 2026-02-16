function [FigureHandles] = Plot_Current_Trace(FigureHandles,CurrentVirusData,UniversalData,...
    CurrTrace_Corrected,PlotCounter,CurrentTraceNumber,Options)

FusionData = CurrentVirusData.FusionData;
DockingData = CurrentVirusData.DockingData;
TraceGradData = CurrentVirusData.TraceGradData;
ChangedByUser = CurrentVirusData.ChangedByUser;

% Plot the current trace to the current subplot axis
    set(0,'CurrentFigure',FigureHandles.MasterWindow)
    set(FigureHandles.MasterWindow,'CurrentAxes',FigureHandles.SubHandles(PlotCounter));
    cla
    
    hold on

if strcmp(FusionData.Designation,'2 Fuse')
    
    plot(CurrTrace_Corrected,'r-')  
    LineToPlot = ylim;
    set(gca,'XTickLabel',''); 
    set(gca,'YTickLabel','');
    
    XToPlot = [FusionData.FuseFrameNumbers(1), FusionData.FuseFrameNumbers(1)];
    plot(XToPlot,LineToPlot,'g--')
    XToPlot = [FusionData.FuseFrameNumbers(2), FusionData.FuseFrameNumbers(2)];
    plot(XToPlot,LineToPlot,'m--')

elseif strcmp(FusionData.Designation,'1 Fuse')
    
    % HERE IS WHERE YOU CHANGE THE FILTER CONDITION. YOU WOULD ONLY CHANGE
    % THE PART AFTER THE && SIGN
    if strcmp(Options.ApplyFilter,'y') && FusionData.FuseFrameNumbers(1) > 960  
        plot(CurrTrace_Corrected,'y-')
            % We plot as a yellow trace if the filter condition is true
    else
        plot(CurrTrace_Corrected,'b-')
    end

    LineToPlot = ylim;
    set(gca,'XTickLabel',''); 
    set(gca,'YTickLabel','');
    
    if strcmp(ChangedByUser,'Incorrect Designation-Changed') || strcmp(ChangedByUser,'Correct Designation, Incorrect Wait Time')
        % No wait time to plot if this designation
    else
        XToPlot = [FusionData.FuseFrameNumbers(1), FusionData.FuseFrameNumbers(1)];
        plot(XToPlot,LineToPlot,'g--')
    end
    
elseif strcmp(FusionData.Designation,'No Fusion')
    
    plot(CurrTrace_Corrected,'k-')
    LineToPlot = ylim;
    set(gca,'XTickLabel',''); 
    set(gca,'YTickLabel','');
    
elseif strcmp(FusionData.Designation,'Slow')
    
    plot(CurrTrace_Corrected,'r-')
    LineToPlot = ylim;
    set(gca,'XTickLabel',''); 
    set(gca,'YTickLabel','');
    
    TraceRunMedian = TraceGradData.TraceRunMedian;
    DiffPosClusterData = TraceGradData.DiffPosClusterData;
    SlowFusePosFrameNumbers = TraceRunMedian.FrameNumbers(DiffPosClusterData.ClusterStartIndices);
    for d = 1:length(SlowFusePosFrameNumbers)
        XToPlot = [SlowFusePosFrameNumbers(d), SlowFusePosFrameNumbers(d)];
        plot(XToPlot,LineToPlot,'c--')
    end
    
elseif strcmp(FusionData.Designation,'Unbound')
    
    plot(CurrTrace_Corrected,'c-')
    LineToPlot = ylim;
    set(gca,'XTickLabel',''); 
    set(gca,'YTickLabel','');
    
else
    
    plot(CurrTrace_Corrected,'k.')
    LineToPlot = ylim;
    set(gca,'XTickLabel',''); 
    set(gca,'YTickLabel','');
end

    Title = strcat(num2str(PlotCounter),'-',FusionData.Designation,'-',num2str(CurrentTraceNumber),...
        '-ID:',num2str(CurrentVirusData.VirusIDNumber),'-X',num2str(round(CurrentVirusData.Coordinates(1))),...
        ',Y',num2str(round(CurrentVirusData.Coordinates(2))));
    title(Title);
    
    if strcmp(Options.ShowBindFrame,'y')
        XToPlot = [UniversalData.StandardBindFrameNum, UniversalData.StandardBindFrameNum];
        plot(XToPlot,LineToPlot,'k--')
    end
    
    if strcmp(Options.ExpandXAxis,'y')
        xlim([-30 length(CurrTrace_Corrected)+30])
    end
    

if strcmp(DockingData.IsMobile ,'y')
    XToPlot = [DockingData.StopFrameNum, DockingData.StopFrameNum];
    plot(XToPlot,LineToPlot,'b--')
end
     
hold off
drawnow
end
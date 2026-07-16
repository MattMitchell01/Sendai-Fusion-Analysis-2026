function [AllResults] = Calculate_Bootstrap_Confidence_Intervals(CumX,CumYNorm,FigureHandles,CurrentColor,...
    AllResults,WaitingTimeList,Options)

    NumberBootstraps = Options.NumberBootstraps;
    
    SourceDistribution = WaitingTimeList;
    NumberDataPoints = length(SourceDistribution);
    IndexMatrix = randi(NumberDataPoints,[NumberDataPoints,NumberBootstraps]);
    BootstrapMatrix = SourceDistribution(IndexMatrix);
    BootstrapMatrix  = sort(BootstrapMatrix);
    
    NumberCDFPoints = length(CumX);
    BootstrapCDFMatrix = zeros(NumberCDFPoints,NumberBootstraps);
    
    for j = 1:NumberBootstraps
        CurrentBootstrapVector = BootstrapMatrix(:,j);
        for i = 1:NumberCDFPoints
            CDFValue = find(CurrentBootstrapVector == CumX(i),1,'last');
            if ~isempty(CDFValue)
                BootstrapCDFMatrix(i,j) = CDFValue;
            else
                if i == 1
                    BootstrapCDFMatrix(i,j) = 0;
                else
                    BootstrapCDFMatrix(i,j) = BootstrapCDFMatrix(i-1,j);
                end
            end
        end
    end
    
    ConfidenceInterval = Options.ConfidenceInterval;
    PercentileRange = [100 - ConfidenceInterval,ConfidenceInterval];
    CDFPercentiles = prctile(BootstrapCDFMatrix,PercentileRange,2);
    CDFPercentilesNorm = CDFPercentiles/NumberDataPoints;
    
    set(0,'CurrentFigure',FigureHandles.BootstrapWindow)
        hold on
        plot(CumX,CumYNorm,CurrentColor.DataPoints);
        plot(CumX,CDFPercentilesNorm,CurrentColor.FitLine,'LineWidth',2)
%         plot(WaitingTimeList,(1:length(WaitingTimeList))/length(WaitingTimeList),'mx');
        drawnow
end
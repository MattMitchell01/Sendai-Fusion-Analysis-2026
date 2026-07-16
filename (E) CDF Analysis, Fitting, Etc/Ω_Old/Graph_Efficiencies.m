function [] = Graph_Efficiencies(FigureHandles,ResultsReport,Options)

    NumberFiles = length(ResultsReport);
    
    for j=  1:NumberFiles
        XLabels{j,1} = ResultsReport(j).Name(1:min(length(ResultsReport(j).Name),10));
        BarData(j,1) = ResultsReport(j).PercentFuse1*100;%
    
        if strcmp(Options.ShowPercentAnyFuseAnd1FuseTogether,'y')
            BarData(j,2) = ResultsReport(j).PercentAnyFusion*100;
        end
    
        NumberVirus(j) = ResultsReport(j).NumVirusFuse1_Total;
        
        
        if strcmp(Options.BootstrapEfficiencies,'y')
             [BootstrapError(j,1:2),BootstrapData(j).FusePercents] =...
                 Calculate_Bootstrap_Error(NumberVirus(j),ResultsReport(j).PercentFuse1,Options);
            
             BootstrapData(j).PercentFuse = ResultsReport(j).PercentFuse1;
             BootstrapData(j).TotalNumberAnalyzed = round(NumberVirus(j)/ResultsReport(j).PercentFuse1);
             BootstrapData(j).NumberVirusFused = NumberVirus(j);
        
        end
        
    end
    
    if strcmp(Options.BootstrapEfficiencies,'y')
        Calculate_P_Value(BootstrapData,ResultsReport,Options);
        Error = BootstrapError*100;
    else
        Error = zeros(size(BarData));
    end
    
    
    set(0,'CurrentFigure',FigureHandles.EfficiencyWindow)
    
    BarHandle = barwitherr(Error, BarData);
    % BarHandle(1).BarWidth = 0.5;
    set(gca,'XTickLabel',XLabels)
    ylabel('% Fusion (Lipid Mixing)')
    %set(BarHandle(1),'FaceColor','y');

    BarHandle.FaceColor = 'flat';
    
    for i = 1:size(BarHandle.CData,1)
        InputInfo.FileNumber = i;
        [CurrentColor] = Choose_Color(InputInfo);

        BarHandle.CData(i,:) = CurrentColor.FaceValue;
    
    end


    CurrentAxes = gca;
    CurrentAxes.XTickLabelRotation=45;
end

function [BootstrapError,FusePercents] = Calculate_Bootstrap_Error(NumberFuseEvents,PercentFuse1, Options)

NumberBootstraps = Options.EfficiencyNumberBootstraps;
ConfidenceInterval = Options.EfficiencyConfidenceInterval;
TotalNumberVirus = round(NumberFuseEvents/PercentFuse1);

SourceDistribution = zeros(TotalNumberVirus,1);
SourceDistribution(1:NumberFuseEvents) = 1;

IndexMatrix = randi(TotalNumberVirus,[TotalNumberVirus,NumberBootstraps]);
BootstrapMatrix = SourceDistribution(IndexMatrix);

FuseMatrix = BootstrapMatrix == 1;
% NoFuseMatrix = BootstrapMatrix == 0;

FuseCounts = sum(FuseMatrix,1);
FusePercents = FuseCounts/TotalNumberVirus;
% NoFuseCounts = sum(NoFuseMatrix,1);
% 
% figure
% histogram(FusePercents)

PercentileRange = [100 - ConfidenceInterval,ConfidenceInterval];
BootstrapErrorCounts = prctile(FuseCounts,PercentileRange);
BootstrapError = prctile(FusePercents,PercentileRange);

BootstrapError = [PercentFuse1 - BootstrapError(1), BootstrapError(2) - PercentFuse1]

end

function [] = Calculate_P_Value(BootstrapData,ResultsReport,Options)
   
NumberFiles = length(ResultsReport);
PValueMatrix = cell(NumberFiles +1,NumberFiles +1);
PValueMatrix{1,1} = 'p, Efficiencies';

NumberBootstraps = Options.EfficiencyNumberBootstraps;

for j=  1:NumberFiles
    Name = ResultsReport(j).Name(1:min(length(ResultsReport(j).Name),8));
    PValueMatrix{j+1,1} = Name;
    PValueMatrix{1,j+1} = Name;
    
    CurrentTotalNumberVirus = BootstrapData(j).TotalNumberAnalyzed;
    CurrentNumberVirusFused = BootstrapData(j).NumberVirusFused;
    CurrentPercentFuse = BootstrapData(j).PercentFuse;
    
%     CurrentBootstrapData = BootstrapData(j).FusePercents;
%     CurrentMedian = median(CurrentBootstrapData);
%     CurrentMax = max(CurrentBootstrapData);
    
    for k=  1:NumberFiles
    
        TotalNumberVirusToCompare = BootstrapData(k).TotalNumberAnalyzed;
        NumberVirusFusedToCompare = BootstrapData(k).NumberVirusFused;
        PercentFuseToCompare = BootstrapData(k).PercentFuse;
    
        % Generate the distribution for the null hypothesis
        TotalNumberVirusNull = TotalNumberVirusToCompare + CurrentTotalNumberVirus;
        NumberVirusFusedNull = NumberVirusFusedToCompare + CurrentNumberVirusFused;
        
        IndexMatrix1 = randi(TotalNumberVirusNull,[CurrentTotalNumberVirus,NumberBootstraps]);
        IndexMatrix2 = randi(TotalNumberVirusNull,[TotalNumberVirusToCompare,NumberBootstraps]);
        
        BootstrapMatrix1 = IndexMatrix1 <= NumberVirusFusedNull;
        BootstrapMatrix2 = IndexMatrix2 <= NumberVirusFusedNull;

        FusePercents1 = sum(BootstrapMatrix1,1)/CurrentTotalNumberVirus *100;
        FusePercents2 = sum(BootstrapMatrix2,1)/TotalNumberVirusToCompare *100;
        
        NullDistribution = FusePercents1 - FusePercents2;
        
%         figure
%         histogram(NullDistribution)
%         histogram(FusePercents1)
        
        % Determine P value using null hypothesis distribution that the two fusion efficiencies are the same
        TestStatistic = abs(CurrentPercentFuse - PercentFuseToCompare) *100;
        
            % Multiply by 2 to make it a two-tailed P value
            PValue = 2*(sum(NullDistribution > TestStatistic))/length(NullDistribution);
% Old method to bootstrap P value (incorrect I believe)
%         BootstrapDataToCompare = BootstrapData(k).FusePercents;
%         MedianToCompare = median(BootstrapDataToCompare);
%         MaxToCompare= max(BootstrapDataToCompare);
%         MinToCompare= min(BootstrapDataToCompare);
%         
%         if CurrentMedian > MedianToCompare
%             ComparisonData = CurrentBootstrapData - BootstrapDataToCompare;
%             FilteredData = ComparisonData <=0;
%             FilteredData = CurrentBootstrapData<=MaxToCompare; 
%             FilteredData = CurrentBootstrapData<=MedianToCompare; 
%         else
%             ComparisonData = CurrentBootstrapData - BootstrapDataToCompare;
%             FilteredData = ComparisonData >=0;
%             FilteredData = CurrentBootstrapData>=MinToCompare;
%             FilteredData = CurrentBootstrapData>=MedianToCompare;
%         end
%         PValue = sum((FilteredData))/NumberBootstraps;
        
        PValueMatrix{k+1,j+1} = PValue;
    end
end

PValueMatrix

end
function [AllDataToAnalyze,CDFReports] = Randomness_Parameter_Analysis(AllDataToAnalyze,CDFReports, Options)

    NumberDataFiles = length(AllDataToAnalyze);
    
    for FileNumber = 1:NumberDataFiles
        SortedWaitTimeList = AllDataToAnalyze(FileNumber).CDFData.SortedWaitTimeList;
        
        % Shift data to account for dead time
        if strcmp(Options.AccountForDeadTime,'y') 
            SortedWaitTimeList = SortedWaitTimeList - AllDataToAnalyze(FileNumber).CDFData.DeadTime;
        end
        
        MeanWaitingTimeSquared = (mean(SortedWaitTimeList))^2;
        MeanOfSquared = mean(SortedWaitTimeList.^2);
        
        RandomnessParameter = (MeanOfSquared - MeanWaitingTimeSquared)/MeanWaitingTimeSquared;
        Nmin = 1/RandomnessParameter;
        
        [NMinPercentiles,RandParamPercentiles, MeanNMin, MeanRandParam] = Bootstrap_Rand_Param(SortedWaitTimeList,Options);
        
        AllDataToAnalyze(FileNumber).RandomnessParameter = RandomnessParameter;
        AllDataToAnalyze(FileNumber).Nmin = Nmin;

        % Add to CDFReport
        CDFReports(FileNumber).RandParam_Mean = MeanRandParam;
        CDFReports(FileNumber).RandParam_ConfInts = RandParamPercentiles;
        CDFReports(FileNumber).Nmin_Mean = MeanNMin;
        CDFReports(FileNumber).Nmin_ConfInts = NMinPercentiles;
        CDFReports(FileNumber).RPNmin_ConfIntOptionUsed = Options.ConfInt_RandParam;
        CDFReports(FileNumber).RPNMin_DeadTimeAccountedFor = Options.AccountForDeadTime;
    end

end


function [NMinPercentiles,RandParamPercentiles, MeanNMin, MeanRandParam] = Bootstrap_Rand_Param(SortedWaitTimeList,Options)

    NumberBootstraps = Options.NumberBootstraps_RandParam;
    
    SourceDistribution = SortedWaitTimeList;
    NumberDataPoints = length(SourceDistribution);
    IndexMatrix = randi(NumberDataPoints,[NumberDataPoints,NumberBootstraps]);
    BootstrapMatrix = SourceDistribution(IndexMatrix);
    %BootstrapMatrix  = sort(BootstrapMatrix);
    
    RandParamList = zeros(1,NumberBootstraps);
    NminList = zeros(1,NumberBootstraps);

    for j = 1:NumberBootstraps
        CurrentBootstrapVector = BootstrapMatrix(:,j);
        CurrMeanWaitingTimeSquared = (mean(CurrentBootstrapVector))^2;
        CurrMeanOfSquared = mean(CurrentBootstrapVector.^2);
        
        CurrRandomnessParameter = (CurrMeanOfSquared - CurrMeanWaitingTimeSquared)/CurrMeanWaitingTimeSquared;
        CurrNmin = 1/CurrRandomnessParameter;

        RandParamList(j) = CurrRandomnessParameter;
        NminList(j) = CurrNmin;
    end
    
    % Calc confidence intervals
        ConfidenceInterval = Options.ConfInt_RandParam;
        PercentileRange = [100 - ConfidenceInterval,ConfidenceInterval];
        NMinPercentiles = prctile(NminList,PercentileRange,2);
        MeanRandParam = mean(RandParamList);
        MeanNMin = mean(NminList);
        RandParamPercentiles = prctile(RandParamList,PercentileRange,2);
    
    %set(0,'CurrentFigure',FigureHandles.BootstrapWindow)
        %hold on
        % Diagnostic plots
            % figure(47)
            % hist(RandParamList)
            % title(strcat("Rand Param, mean = ", num2str(MeanRandParam)));
            % 
            % figure(48)
            % hist(NminList)
            % title(strcat("Nmin, mean = ", num2str(MeanNMin)));
            % drawnow
end
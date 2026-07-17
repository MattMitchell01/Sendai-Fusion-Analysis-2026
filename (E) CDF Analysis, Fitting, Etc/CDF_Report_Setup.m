function [CDFReports] = CDF_Report_Setup(AllDataToAnalyze)
    % Here we set up the basic results report. Other items will be added in
    % the different modules depending on analysis options you chose.

    NumberDataFiles = length(AllDataToAnalyze);

    for FileNumber = 1:NumberDataFiles
            CurrCDFData = AllDataToAnalyze(FileNumber).CDFData;

            CDFReports(FileNumber).DataName = AllDataToAnalyze(FileNumber).CDFData.DataLabelForPlot;
            CDFReports(FileNumber).FusionTriggerType = AllDataToAnalyze(FileNumber).CDFData.FusionTriggerType;
            if AllDataToAnalyze(FileNumber).CDFData.DeadTime ~=0
                CDFReports(FileNumber).DeadTime = AllDataToAnalyze(FileNumber).CDFData.DeadTime;
            else
                CDFReports(FileNumber).DeadTime = "Not used - set at zero";
            end
            CDFReports(FileNumber).MeanpHtoFuse= mean(AllDataToAnalyze(FileNumber).CDFData.SortedWaitTimeList);
            CDFReports(FileNumber).Efficiency = num2str(CurrCDFData.UsefulInfo.NumberFusedDataPoints) + "/" + ...
                        num2str(CurrCDFData.UsefulInfo.NumberTotalAnalyzed) + "; " + ...
                        num2str(CurrCDFData.Efficiency*100,'%.1f') + "%";
            CDFReports(FileNumber).IncludedInEfficiencyCalculation = AllDataToAnalyze(FileNumber).CDFData.Options.ToIncludeInEfficiencyCalculation;
            CDFReports(FileNumber).Unbinding= num2str(CurrCDFData.UsefulInfo.NumberUnbound) + "/" + ...
                        num2str(CurrCDFData.UsefulInfo.NumberVirusTotalInFractUnboundCalc) + "; " + num2str(CurrCDFData.UsefulInfo.FractionUnbound*100,'%.1f') + "%";
            
            if isfield(CurrCDFData.UsefulInfo, 'NumberOther')
                CDFReports(FileNumber).NumberOther = CurrCDFData.UsefulInfo.NumberOther;
            elseif isfield(CurrCDFData.UsefulInfo,'NumberSlow_Other')
                display("WARNING: Old File Format - make sure you want to do this")
                CDFReports(FileNumber).NumberSlow_Other = CurrCDFData.UsefulInfo.NumberSlow_Other;
            else
                display("Incorrect File Format??")
                StopProgram
            end
            
            

    end

end
function [FigureHandles] = Plot_CDFs(AllDataToAnalyze,Options,FigureHandles)
    NumberDataFiles = length(AllDataToAnalyze);

    for FileNumber = 1:NumberDataFiles

        CurrCDFData = AllDataToAnalyze(FileNumber).CDFData;
        
        % Plot current CDF data in CDF window
            [CurrentColor] = Choose_Color('Data',FileNumber);
            
            set(0,'CurrentFigure',FigureHandles.CDFWindow)
            hold on
                CumXToPlot = CurrCDFData.CumX;
    
                if strcmp(Options.NormalizeCDFsInComparisonPlot,'y')
                    CumYToPlot = CurrCDFData.CumYNorm;
                else
                    CumYToPlot = CurrCDFData.CumYNorm.*CurrCDFData.Efficiency;
                end
    
                plot(CumXToPlot,CumYToPlot, CurrentColor.DataPoints);

        % Set up legend for CDF window
            LegendInfoCDF{1,FileNumber} = CurrCDFData.DataLabelForPlot + "; " +...
                num2str(CurrCDFData.UsefulInfo.NumberFusedDataPoints) + "/" + ...
                num2str(CurrCDFData.UsefulInfo.NumberTotalAnalyzed) + "; " + ...
                num2str(CurrCDFData.Efficiency*100,'%.1f') + "%" ;
    end

    % Show legend, title, xylabels in CDF window
        set(0,'CurrentFigure',FigureHandles.CDFWindow)
                
            if strcmp(CurrCDFData.FusionTriggerType,'Binding')
                xlabel('Waiting Time (min)');
            elseif strcmp(CurrCDFData.FusionTriggerType,'pH')
                xlabel('Waiting Time (sec)');
            else
                xlabel('Waiting Time (data type undefined)');
            end
    
            if strcmp(Options.NormalizeCDFsInComparisonPlot,'y')
                CDFPlotTitle = "Data Window: CDFs Normalized";
                ylabel('Normalized Prob of Lipid Mixing');
                ylim([0 1]);
            else
                CDFPlotTitle = "Data Window: CDFs NOT Normalized";
                ylabel('Prop of Lipid Mixing Events');
            end
    
            legend(LegendInfoCDF,'Location','southeast','Interpreter', 'none');
            title(CDFPlotTitle, 'Interpreter', 'none')


end
function [AllDataToAnalyze,CDFReports] = Fit_CDFs_And_Plot_Fits(AllDataToAnalyze,CDFReports, FigureHandles, Options)

    % Determine number of fits and number of files
    if strcmp(Options.RunFit_SingleFitMultipleCDFs,'y') && strcmp(Options.RunFit_MultipleFitsSingleCDF,'n')
        NumFitsToPerform = 1;
        NumberDataFilesToFit = length(AllDataToAnalyze);
    elseif strcmp(Options.RunFit_SingleFitMultipleCDFs,'n') && strcmp(Options.RunFit_MultipleFitsSingleCDF,'y')
        NumFitsToPerform = length(Options.FitToPerform);
        NumberDataFilesToFit = 1;
    else
        disp("!!!!ERROR: Set up your fit options correctly please.")
        StopProgramNow
    end

    NumThingsPlotted_FitWind = 0;
    NumThingsPlotted_ResidWind = 0;


    for FileNumber = 1:NumberDataFilesToFit
        
        CurrCDFData = AllDataToAnalyze(FileNumber).CDFData;
        
        % Plot current CDF data in CDF window
            [CurrentColor] = Choose_Color('Data',FileNumber);
            
            set(0,'CurrentFigure',FigureHandles.FitWindow)
            hold on
                CumXToPlot = CurrCDFData.CumX;
                CumYToPlot = CurrCDFData.CumYNorm;
                plot(CumXToPlot,CumYToPlot, CurrentColor.DataPoints);
                NumThingsPlotted_FitWind = NumThingsPlotted_FitWind + 1;

            % Record legend info for CDF window
                LegendInfo_Fit{1,NumThingsPlotted_FitWind} = CurrCDFData.DataLabelForPlot;

        % Run and Plot Fit Solutions
            for FitNumber = 1:NumFitsToPerform
                
                FitType = Options.FitToPerform(FitNumber).FitType;

                % Send info to sorter, which then runs the fit script
                    [AllDataToAnalyze,CDFReports,FitLine] = Fit_Type_Sorter(AllDataToAnalyze,CDFReports, Options, FitType, FileNumber);

                % Plot Current Fit Solution
                    if strcmp(Options.RunFit_SingleFitMultipleCDFs,'y')
                        [CurrentColor] = Choose_Color('Single Fit',FileNumber);
                    elseif strcmp(Options.RunFit_MultipleFitsSingleCDF,'y')
                        [CurrentColor] = Choose_Color('MultipleFits',FitNumber);
                    end
                
                    set(0,'CurrentFigure',FigureHandles.FitWindow)
                    hold on
                        plot(FitLine.FitXValues,FitLine.FitYValues, CurrentColor.FitLine, 'LineWidth',1.5);
                        NumThingsPlotted_FitWind = NumThingsPlotted_FitWind + 1;

                    % Record legend info for fit window
                        LegendInfo_Fit{1,NumThingsPlotted_FitWind} = FitType;

                % Plot Residuals
                    set(0,'CurrentFigure',FigureHandles.ResidualsWindow)
                    hold on
                        Residuals = CumYToPlot - FitLine.FitYValues;
                        plot(FitLine.FitXValues,Residuals, CurrentColor.ResidualPoints);
                        NumThingsPlotted_ResidWind = NumThingsPlotted_ResidWind + 1;

                    % Record legend info for fit window
                        if strcmp(Options.RunFit_SingleFitMultipleCDFs,'y')
                            LegendInfo_Resid{1,NumThingsPlotted_ResidWind} = CurrCDFData.DataLabelForPlot;
                        elseif strcmp(Options.RunFit_MultipleFitsSingleCDF,'y')
                            LegendInfo_Resid{1,NumThingsPlotted_ResidWind} = FitType;
                        end
            end

        
    end

    % Show legend, title, xylabels in fit window
        set(0,'CurrentFigure',FigureHandles.FitWindow)  
            if strcmp(CurrCDFData.FusionTriggerType,'Binding')
                LabelX = 'Waiting Time (min)';
            elseif strcmp(CurrCDFData.FusionTriggerType,'pH')
                LabelX = 'Waiting Time (sec)';
            else
                LabelX = 'Waiting Time (data type undefined)';
            end

            xlabel(LabelX);
            ylabel('Normalized Prob of Lipid Mixing');
            ylim([0 1]);
            title("Fit Window")
          
            legend(LegendInfo_Fit,'Location','southeast','Interpreter', 'none');
            

    % Show legend, title, xylabels in residual window
        set(0,'CurrentFigure',FigureHandles.ResidualsWindow)  
            xlabel(LabelX);
            ylabel('Residuals');
            legend(LegendInfo_Resid,'Location','southeast','Interpreter', 'none');
            title("Residuals Plot")

            drawnow
end
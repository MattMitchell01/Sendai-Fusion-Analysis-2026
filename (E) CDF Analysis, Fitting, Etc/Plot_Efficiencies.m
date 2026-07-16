function [FigureHandles] = Plot_Efficiencies(AllDataToAnalyze,Options,FigureHandles)
    
    %if strcmp(Options.SimpleComparisonOnly,'y')             
         NumberDataFiles = length(AllDataToAnalyze);
    
         for FileNumber = 1:NumberDataFiles
            CurrCDFData = AllDataToAnalyze(FileNumber).CDFData;
            
            % Plot current efficiency data in efficiency window
                [CurrentColor] = Choose_Color('Data',FileNumber);
                
                set(0,'CurrentFigure',FigureHandles.EfficiencyWindow)
                hold on
                    x_bar = string(CurrCDFData.DataLabelForPlot);
                    y_bar = CurrCDFData.Efficiency;
        
                    bar(x_bar,y_bar, 'FaceColor',CurrentColor.FaceValue);
         end

         ylabel("Efficiency (Fraction)")
         title("Efficiency Plot")
     
    %elseif strcmp(Options.SimpleComparisonOnly,'n')
        % To add: bootstrap errors, etc

    %end
end


       
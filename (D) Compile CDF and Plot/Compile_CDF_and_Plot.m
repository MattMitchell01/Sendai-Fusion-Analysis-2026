function [ErrorFlag,DataToSave,CDFData] = Compile_CDF_and_Plot(DataFileName,...
    DataFilePath,FigureHandles,Options,DataToSave)
       
    % Extract useful info: filename, etc.
        TextFilename = DataFileName;
            IdxOfDot = find(TextFilename=='.');
            TextFilenameWODot = TextFilename(1:IdxOfDot-1);
        UsefulInfo.Filename = TextFilename;
        UsefulInfo.Name = TextFilenameWODot + "-CDF"; %This will be the combo filename if multiple files are combined
        UsefulInfo.Filepath = DataFilePath;

        if strcmp(Options.ExcludeOutlierPoints,'y')
            UsefulInfo.Name = UsefulInfo.Name + "-EDITED";
        end

        if Options.TimeCutoffLow ~= 0
            if Options.DeadTime ~=0
                UsefulInfo.Name = UsefulInfo.Name + "-DEAD" + num2str(Options.DeadTime);
            else
                UsefulInfo.Name = UsefulInfo.Name + "-BegTRIMMED";
            end
        end

        if ~isnan(Options.TimeCutoffHigh) 
            UsefulInfo.Name = UsefulInfo.Name + "-EndTRIMMED";
        end


    % Extract data, compile non-normalized CDF, calculate efficiency
    [ErrorFlag,SortedWaitTimeList,CumX,CumY,Efficiency,UsefulInfo] = Extract_Data(DataToSave,...
        UsefulInfo,Options);

    %Compile normalized CDF data
        CumYDecay = max(CumY)-CumY;
        CumYDecayNorm = CumYDecay/max(CumY);
        CumYNorm = CumY/max(CumY);

    %Record the CDF data
        CDFData.CumX = CumX;
        CDFData.CumY = CumY;
        CDFData.CumYNorm = CumYNorm;
        CDFData.OtherCumYVersions.CumYDecay = CumYDecay;
        CDFData.OtherCumYVersions.CumYDecayNorm = CumYDecayNorm;
        CDFData.SortedWaitTimeList = SortedWaitTimeList;
        CDFData.Name = UsefulInfo.Name;
        CDFData.Efficiency = Efficiency;
        CDFData.FusionTriggerType = Options.FusionTrigger;
        CDFData.DeadTime = Options.DeadTime;
        CDFData.UsefulInfo = UsefulInfo;
        CDFData.Options = Options;
   
 % Plot the data
        
    % CDF
        set(0,'CurrentFigure',FigureHandles.CDFWindow)
        plot(CumX,CumYNorm, 'bo');
        LegendInfo = num2str(CDFData.UsefulInfo.NumberFusedDataPoints) + "/" + ...
            num2str(CDFData.UsefulInfo.NumberTotalAnalyzed) + "; " + num2str(CDFData.Efficiency*100,'%.1f') + "%" ;
    
        legend(LegendInfo,'Location','southeast','Interpreter', 'none');
            title(CDFData.Name,'Interpreter', 'none')
            xlabel('Waiting Time (min)');
            ylabel('Normalized Prob of Lipid Mixing');
            %ylim([0 1]);
            xlim([0 max(CumX)]);

    % Bar graph of efficiency and unbinding
        set(0,'CurrentFigure',FigureHandles.EffUnbWindow)
        x_bar = ["Efficiency" "Unbinding"];
        y_bar = [CDFData.Efficiency UsefulInfo.FractionUnbound];
        bar(x_bar, y_bar);
        ylabel("Fraction")
        title(CDFData.Name,'Interpreter', 'none')

       

% Display stats
    disp("   ---------Report---------")
    disp("   Data = " + CDFData.Name)
    disp("      Efficiency = " + num2str(CDFData.UsefulInfo.NumberFusedDataPoints) + "/" + ...
            num2str(CDFData.UsefulInfo.NumberTotalAnalyzed) + "; " + num2str(CDFData.UsefulInfo.PercentFuse1*100,'%.1f') + "%" )
    disp("      Unbinding = " + num2str(CDFData.UsefulInfo.NumberUnbound) + "/" + ...
            num2str(CDFData.UsefulInfo.NumberVirusTotalInFractUnboundCalc) + "; " + num2str(UsefulInfo.FractionUnbound*100,'%.1f') + "%")
    disp("           Reminder: Unbinding events are typically excluded from efficiency calc")
    disp("   -----------------------")

end
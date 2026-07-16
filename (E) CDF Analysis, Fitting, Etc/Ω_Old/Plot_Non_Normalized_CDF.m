function Plot_Non_Normalized_CDF(FigureHandles,AllResults,XLimit,NumberDataFiles,LegendInfoFuse1)

    set(0,'CurrentFigure', FigureHandles.NotNormalizedFusionWindow)
    hold on

    for w = 1:NumberDataFiles
        InputInfo.FileNumber = w;
        [CurrentColor] = Choose_Color(InputInfo);

        NotNormCumY = AllResults(w).CDFData.CumYNorm .* AllResults(w).ResultsReport.PercentFuse1;
        plot(AllResults(w).CDFData.CumX, NotNormCumY, CurrentColor.DataPoints)

    end

    legend(LegendInfoFuse1,'Location','southeast');
    xlabel('Waiting Time (min)');
    ylabel('Fraction of Virions Fused');
    xlim([0 XLimit]);

end
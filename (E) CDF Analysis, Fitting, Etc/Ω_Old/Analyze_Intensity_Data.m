function [FigureHandles] = Analyze_Intensity_Data(FigureHandles,Compiled_Fuse1_Data,UsefulInfo,XLimitForPlot)

% Change color as we go along
[CurrentColor]=Choose_Color(UsefulInfo);

FigureHandles.IntensityWindow = figure(4);
FigureHandles.IntensityCDFWind = figure(5);
set(FigureHandles.IntensityWindow, 'Position', [562   1    560   420]);
set(FigureHandles.IntensityCDFWind, 'Position', [-5   1    560   420]);

for b = 1:length(Compiled_Fuse1_Data)
    pHtoFList(b) =  Compiled_Fuse1_Data(b).pHtoFuse1Time;
    IntensityList(b) = Compiled_Fuse1_Data(b).IntensityJumpFuse1;
end

    IdxToUse = pHtoFList > UsefulInfo.TimeCutoff;
    pHtoFList = pHtoFList(IdxToUse);
    IntensityList = IntensityList(IdxToUse);
[SortedFusionList,IndexSortedList] = sort(pHtoFList);
SortedIntensityList = IntensityList(IndexSortedList);
Y_CDF = 1:length(SortedFusionList);

for k = 1:length(SortedIntensityList)
    if k == 1
        Y_Intensity(k) = SortedIntensityList(k);
    else
        Y_Intensity(k) = SortedIntensityList(k) + Y_Intensity(k-1);
    end
end
Y_IntensityNorm = Y_Intensity/max(Y_Intensity);
Y_CDFNorm = Y_CDF/max(Y_CDF);

set(0,'CurrentFigure',FigureHandles.IntensityWindow)
hold on
plot(pHtoFList,IntensityList, CurrentColor.DataPoints);
xlabel('Waiting Time (s)');
ylabel('Intensity Jump');
xlim([0 XLimitForPlot]);

set(0,'CurrentFigure',FigureHandles.IntensityCDFWind)
hold on
plot(SortedFusionList,Y_IntensityNorm, 'ks');
plot(SortedFusionList,Y_CDFNorm, CurrentColor.DataPoints);
xlabel('Waiting Time (s)');
ylabel('ProportionFused');
xlim([0 XLimitForPlot]);

end
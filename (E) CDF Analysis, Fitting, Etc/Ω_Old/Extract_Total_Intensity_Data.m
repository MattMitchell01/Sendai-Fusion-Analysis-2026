function [TimePoints,IntensityValuesRunMedian,MeanFusionTimeRough,NumberDataPoints] =...
    Extract_Total_Intensity_Data(InputData,FileNumber,FigureHandles,CurrentColor)

if FileNumber == 1.753
    TimeInterval = 1.152;
else
    TimeInterval = .288;
end



if isfield(InputData,'OtherDataToSave')
    PHdropFrameNum = InputData.OtherDataToSave.PHdropFrameNum;
    AverageVideoIntensity= InputData.OtherDataToSave.AverageVideoIntensity;
    TotalVideoIntensity= InputData.OtherDataToSave.TotalVideoIntensity;
elseif isfield(InputData,'Other_Data_To_Save')
    PHdropFrameNum = InputData.Other_Data_To_Save.PHdropFrameNum;
    AverageVideoIntensity= InputData.Other_Data_To_Save.AverageVideoIntensity;
    TotalVideoIntensity= InputData.Other_Data_To_Save.TotalVideoIntensity;
end

%     PHdropFrameNum = PHdropFrameNum + 4;
% PHdropFrameNum = 1;
% TimeInterval = 5;

IntensityValues = TotalVideoIntensity(PHdropFrameNum:end);
    IntensityValues = IntensityValues';

TimePoints = (1:length(IntensityValues))*TimeInterval;
NumberDataPoints = length(IntensityValues);

    RunMedHalfLength = 1;
        StartIdx = RunMedHalfLength + 1;
        EndIdx = length(IntensityValues)-RunMedHalfLength;
    IntensityValuesRunMedian = zeros(1,length(IntensityValues));
        IntensityValuesRunMedian(1:StartIdx) = mean(IntensityValues(1:StartIdx));
        IntensityValuesRunMedian(EndIdx:end) = mean(IntensityValues(EndIdx:end));
        
    for n = StartIdx:EndIdx
        IntensityValuesRunMedian(n) = mean(IntensityValues(n-RunMedHalfLength:n+RunMedHalfLength));
    end
    
% Set Initial Baseline to Zero
    IntensityValuesRunMedian = IntensityValuesRunMedian - median(IntensityValuesRunMedian(1:10));
    
    DistanceToMidpoint = abs(IntensityValuesRunMedian - max(IntensityValuesRunMedian)/2);
    Index50Per = DistanceToMidpoint == min(DistanceToMidpoint);
    MeanFusionTimeRough = TimePoints(Index50Per);
    MeanFusionTimeRough = MeanFusionTimeRough(1);

    
% Plot raw intensity values baselined to zero
    AverageVideoIntensity = AverageVideoIntensity - AverageVideoIntensity(1);
    RawTimePoints = (1:length(AverageVideoIntensity))*TimeInterval;
    set(0,'CurrentFigure',FigureHandles.BootstrapWindow)
    hold on
    plot(RawTimePoints,AverageVideoIntensity,CurrentColor.FitLine)
    xlabel('Time (s)');
    ylabel('Average Intensity');
    
end
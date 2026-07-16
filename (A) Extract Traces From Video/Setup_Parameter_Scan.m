function [FileOptions,NumberOfParameters] = Setup_Parameter_Scan(FileOptions,i,CurrentOptions)

    CurrentOptions.QuickScanThreshold = 'y';
    InputRange = CurrentOptions.ScanRange;
    ParameterRangeDiff = InputRange(2) - InputRange(1);
    NumberIntervals = 4;
    for b=  1:NumberIntervals+1
        ParameterRange(b) = InputRange(1) + (b-1)*ParameterRangeDiff/NumberIntervals;
    end
    
    
    NumberOfParameters = length(ParameterRange);
    CurrentLabel = CurrentOptions.ExtraLabel;
    for j= 1:NumberOfParameters
        CurrentOptions.Threshold = ParameterRange(j);
        CurrentOptions.ExtraLabel =  strcat(CurrentLabel,...
            ';TH=', num2str(ParameterRange(j)));
        FileOptions(i).Parameter(j).Options = CurrentOptions;
    end
    
    if strcmp(CurrentOptions.QuickScanThreshold,'y')
        % Note: only 5 thresholds can be displayed at once
        
        ZoomImage = 'y';
        Quick_Scan_Threshold(FileOptions,NumberOfParameters,ZoomImage)
    end
end
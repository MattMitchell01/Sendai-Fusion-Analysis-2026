function [ErrorFlag,SortedWaitTimeList,CumX,CumY,Efficiency,UsefulInfo] = Extract_Data(DataToSave,...
        UsefulInfo,Options)
%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%

    TypeOfInputData = 'Normal CDF-Improved Analysis';
    ErrorFlag = false;

    if strcmp(TypeOfInputData,'Normal CDF-Improved Analysis')
        % Keeping this if statement just to remind Bob of the way that legacy
        % data was handled in the old program (before 250915).
    
        AnalyzedTraceData = DataToSave.CombinedAnalyzedTraceData;
        NumberofTraces = length(AnalyzedTraceData);
        NumberFuse1 = 0;
        NumberFuse1ToPlot = 0;
        NumberFuse2 = 0;
        NumberNoFuse = 0;
        NumberUnbound = 0;
        NumberOther = 0;
        WaitTimeList = [];
        
        for k = 1:NumberofTraces
            CurrentFusionData = AnalyzedTraceData(k).FusionData;
            if strcmp(AnalyzedTraceData(k).ChangedByUser,'Incorrect Designation-Not Changed')
                % The designation is wrong, but has not been corrected, so we will skip it.
                disp('   ***Some designations have been previously flagged as incorrect, but have not yet been changed, and so have been skipped.')
                
            elseif strcmp(CurrentFusionData.Designation,'Other')
                % This trace was flagged as hard to classify, so we will ignore it
                NumberOther = NumberOther + 1;
                
            else
                if strcmp(CurrentFusionData.Designation,'No Fusion')
                    NumberNoFuse = NumberNoFuse + 1;
                elseif strcmp(CurrentFusionData.Designation,'1 Fuse')
                    if strcmp(AnalyzedTraceData(k).ChangedByUser,'Not analyzed') ||...
                            strcmp(AnalyzedTraceData(k).ChangedByUser,'Reviewed, Fuse frame chosen by user')
                        NumberFuse1 = NumberFuse1 + 1;
                        NumberFuse1ToPlot = NumberFuse1ToPlot + 1;

                        if strcmp(Options.FusionTrigger,'Binding')
                            WaitTimeList(NumberFuse1ToPlot) = CurrentFusionData.BindtoFusionTime(1);
                        elseif strcmp(Options.FusionTrigger,'pH')
                            WaitTimeList(NumberFuse1ToPlot) = CurrentFusionData.pHtoFusionTime(1);
                        else
                            error('Incorrect data type or fusion trigger chosen -- pHtoFusionTime no longer exists in the new schema.');
                        end

                    else
                        % Note: in this case we can't Necessarily Trust The Wait Time
                        NumberFuse1 = NumberFuse1 + 1;
                    end
    
                elseif strcmp(CurrentFusionData.Designation,'2 Fuse')
                    NumberFuse2 = NumberFuse2 + 1;
                elseif strcmp(CurrentFusionData.Designation,'Unbound')
                    NumberUnbound = NumberUnbound + 1;
                else
                    NumberOther = NumberOther + 1;
                end
            end
        end
        
        IdxToUse = WaitTimeList > Options.TimeCutoffLow; 

        if ~isnan(Options.TimeCutoffHigh)
            Index2 = WaitTimeList<Options.TimeCutoffHigh;
            IdxToUse =IdxToUse  +Index2;
            IdxToUse = IdxToUse ==2;
        end
        WaitTimeList = WaitTimeList(IdxToUse);
        WaitTimeList = WaitTimeList';
        SortedWaitTimeList = sort(WaitTimeList);
    
        %Change to cumulative distribution function (not multiple
        %data points at repeated time points)
        if ~isempty(SortedWaitTimeList)
            [CumX, CumY] = Generate_Prop_Cum(SortedWaitTimeList);
        else
            CumX = 0; 
            CumY = 0;
        end
        
        if strcmp(Options.ToIncludeInEfficiencyCalculation,'Fuse1 and No Fusion')
            NumberTotalAnalyzed = NumberFuse1 + NumberNoFuse;
            
        elseif strcmp(Options.ToIncludeInEfficiencyCalculation,'Fuse1, Fuse2, and No Fusion')
            NumberTotalAnalyzed = NumberFuse1 + NumberFuse2 + NumberNoFuse;
            
        elseif strcmp(Options.ToIncludeInEfficiencyCalculation,'Fuse1, Fuse2, Unbound, and No Fusion')
            NumberTotalAnalyzed = NumberFuse1 + NumberFuse2 + NumberNoFuse + NumberUnbound;

        elseif strcmp(Options.ToIncludeInEfficiencyCalculation,'All')
            NumberTotalAnalyzed = NumberFuse1 + NumberFuse2 + NumberNoFuse + NumberUnbound + NumberOther;
            
        else
            disp('   Error: No ToIncludeInEfficiencyCalculation defined')
            ErrorFlag = true;
            return
        end
            
        UsefulInfo.NumberTotalAnalyzed = NumberTotalAnalyzed;
        UsefulInfo.NumberFusedDataPoints = length(SortedWaitTimeList);
        UsefulInfo.MeanFusion1Time = mean(SortedWaitTimeList);
    
        if ~isnan(Options.TimeCutoffHigh)
            UsefulInfo.PercentFuse1 = length(SortedWaitTimeList)/NumberTotalAnalyzed;
            UsefulInfo.TotalNumberFuse1 = length(SortedWaitTimeList);
        else
            UsefulInfo.PercentFuse1 = NumberFuse1/NumberTotalAnalyzed;
            UsefulInfo.TotalNumberFuse1 = NumberFuse1;
        end
        
        UsefulInfo.NumberUnbound = NumberUnbound;
        UsefulInfo.NumberVirusTotalInFractUnboundCalc = NumberFuse1 + NumberNoFuse + NumberUnbound;
        UsefulInfo.FractionUnbound = NumberUnbound./UsefulInfo.NumberVirusTotalInFractUnboundCalc;
        UsefulInfo.NumberOther = NumberOther;
        
        if strcmp(Options.ToIncludeInEfficiencyCalculation,'Fuse1 and No Fusion')
            UsefulInfo.PercentAnyFusion = UsefulInfo.PercentFuse1;
        else
            UsefulInfo.PercentAnyFusion = (NumberFuse1 + NumberFuse2)/NumberTotalAnalyzed;
        end

        % Define efficiency
        Efficiency = UsefulInfo.PercentFuse1;
        
    else
        % Shouldn't be possible to get here - error would have occurred
        % earlier.
    end

end

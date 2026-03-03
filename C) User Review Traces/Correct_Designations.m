function [RerunThisRound, CorrectedAnalysisData, ErrorCounter] = Correct_Designations(IncorrectPlotIndices,...
    PreviousAnalysisData,CurrentTraceRange,CorrectedAnalysisData, ErrorCounter,Options,UniversalData,FigureHandles)

    RerunThisRound = 'n';
    NumberIncorrect = length(IncorrectPlotIndices);

    for j = 1:NumberIncorrect
        CurrentIndex = floor(IncorrectPlotIndices(j));
        MaxRoundPlots = length(CurrentTraceRange);

        if isnan(CurrentIndex) || CurrentIndex < 1 || CurrentIndex > MaxRoundPlots
            fprintf('Invalid plot number %g for this round. Valid plot numbers are 1 to %d.\n', ...
                IncorrectPlotIndices(j), MaxRoundPlots);
            RerunThisRound = 'y';
            break
        end

        TraceNumberIndex = CurrentTraceRange(CurrentIndex);

        if  strcmp(Options.QuickModeNoCorrection,'y')
            CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Incorrect Designation-Not Changed';
            
        elseif strcmp(Options.QuickModeNoCorrection,'n')
            
            CurrentDesignationCode = rem(IncorrectPlotIndices(j), 1);
            if CurrentDesignationCode == 0.0
                if strcmp(PreviousAnalysisData(TraceNumberIndex).Designation, 'No Fusion')
                    % There has been a mistake, re-run the last round
                    RerunThisRound = 'y';
                    break
                else

                    CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Incorrect Designation-Changed';
                    CorrectedAnalysisData(TraceNumberIndex).Designation = 'No Fusion';
                    CorrectedAnalysisData(TraceNumberIndex).FusionData.Designation = 'No Fusion';
                end
            elseif CurrentDesignationCode > 0.09 && CurrentDesignationCode < 0.11
                if strcmp(PreviousAnalysisData(TraceNumberIndex).Designation, '1 Fuse')
                    % There has been a mistake, re-run the last round
                    RerunThisRound = 'y';
                    break
                else
                    CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Incorrect Designation-Changed';
                    CorrectedAnalysisData(TraceNumberIndex).Designation = '1 Fuse';
                    CorrectedAnalysisData(TraceNumberIndex).FusionData.Designation = '1 Fuse';

                    if strcmp(Options.FixWaitTime,"y")
                        [CorrectedAnalysisData] = FixWaitTime(CorrectedAnalysisData,UniversalData,FigureHandles,TraceNumberIndex,Options);
                    end

                end
            elseif CurrentDesignationCode > 0.11 && CurrentDesignationCode < 0.13
                if ~strcmp(PreviousAnalysisData(TraceNumberIndex).Designation, '1 Fuse')
                    % There has been a mistake, re-run the last round
                    RerunThisRound = 'y';
                    break
                else
                    CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Correct Designation, Incorrect Wait Time';
                    CorrectedAnalysisData(TraceNumberIndex).Designation = '1 Fuse';
                    CorrectedAnalysisData(TraceNumberIndex).FusionData.Designation = '1 Fuse';

                    if strcmp(Options.FixWaitTime,"y")
                        [CorrectedAnalysisData] = FixWaitTime(CorrectedAnalysisData,UniversalData,FigureHandles,TraceNumberIndex,Options);
                    end

                end
            elseif CurrentDesignationCode > 0.19 && CurrentDesignationCode < 0.21
                if strcmp(PreviousAnalysisData(TraceNumberIndex).Designation, '2 Fuse') ...
                        || strcmp(PreviousAnalysisData(TraceNumberIndex).Designation, 'Slow')
                    % There has been a mistake, re-run the last round
                    RerunThisRound = 'y';
                    break
                else
                    CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Incorrect Designation-Changed';
                    CorrectedAnalysisData(TraceNumberIndex).Designation = 'Slow';
                    CorrectedAnalysisData(TraceNumberIndex).FusionData.Designation = 'Slow';
                end
            elseif CurrentDesignationCode > 0.29 && CurrentDesignationCode < 0.31
                    if strcmp(PreviousAnalysisData(TraceNumberIndex).Designation, 'Unbound')
                        % There has been a mistake, re-run the last round
                        RerunThisRound = 'y';
                        break
                    else
                        CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Incorrect Designation-Changed';
                        CorrectedAnalysisData(TraceNumberIndex).Designation = 'Unbound';
                        CorrectedAnalysisData(TraceNumberIndex).FusionData.Designation = 'Unbound';
                    end
            elseif CurrentDesignationCode > 0.89 && CurrentDesignationCode < 0.91
                    CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Incorrect Designation-Changed';
                    CorrectedAnalysisData(TraceNumberIndex).Designation = 'Strange-Ignore';
                    CorrectedAnalysisData(TraceNumberIndex).FusionData.Designation = 'Strange-Ignore';
            else
                % There has been an error, re-run the last round to avoid crash
                RerunThisRound = 'y';
                break
            end
        end
    end
    
    if  strcmp( RerunThisRound,'n')
        ErrorCounter = ErrorCounter + NumberIncorrect;
    end
end

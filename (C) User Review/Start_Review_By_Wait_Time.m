function Start_Review_By_Wait_Time(SelectedFilePath)
%
% ------------------------------------------------------------------------
% Written by Matthew D. Mitchell, Rawle Lab, Williams College, 2026.
% New script developed as part of the 2026 overhaul and expansion of the
% Sendai fusion analysis pipeline originally created by Bob Rawle
% (Kasson Lab, University of Virginia, 2016; Rawle et al., Disentangling
% Viral Membrane Fusion from Receptor Binding Using Synthetic DNA-Lipid
% Conjugates, Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% pipeline subsequently updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
% Start_Review_By_Wait_Time  Targeted re-review of '1 Fuse' traces inside
% a reviewer-chosen wait-time window (minutes) -- for chasing down a jump
% in Part D's CDF back to the specific traces causing it.
%
% SelectedFilePath is an existing Part C User Review File (<name>-Revd.mat)
% -- picked by An_Easy_Start_ReviewByWaitTime.m and handed in directly.
% Unlike Start_User_Review.m, this never bootstraps a fresh file from raw
% Part B output -- it only ever operates on an already-reviewed file.
%
% For each window: Find_Traces_In_Wait_Window.m selects the matching
% traces (mirroring Extract_Data.m's own CDF-inclusion filter exactly), 4
% at a time in a 2x2 grid, using the SAME PlotNumber.Code correction
% workflow as Part C (Correct_Designations.m / Get_Designation_Code_Table.m),
% and saves back into the SAME file -- ReviewQueue/CurrentQueuePosition/
% ReviewMeta are carried through untouched, since this tool never advances
% Part C's own Low/Medium/High resume position.

    disp('====================================')
    disp('   Review By Wait Time -- targeted re-review for a CDF jump')

    if ~isfile(SelectedFilePath)
        error('Start_Review_By_Wait_Time:FileNotFound', 'File not found: %s', SelectedFilePath);
    end

    Loaded = load(SelectedFilePath, 'DataToSave', 'ReviewQueue', 'CurrentQueuePosition', 'ReviewMeta');
    DataToSave           = Loaded.DataToSave;
    ReviewQueue          = Loaded.ReviewQueue;
    CurrentQueuePosition = Loaded.CurrentQueuePosition;
    ReviewMeta           = Loaded.ReviewMeta;
    UserReviewFilePath   = SelectedFilePath;

    Options = Setup_Options_User_Review();

    % Same soft-guard pattern Start_User_Review.m uses -- a dataset that
    % predates ClipWidth just quietly loses the running-median dead-zone
    % crop / H2 overlay inputs Plot_Current_Trace.m optionally uses, not
    % an error.
    VideoTimeVector = DataToSave.OtherDataToSave.VideoTimeVector;
    ClipWidth = [];
    if isfield(DataToSave.OtherDataToSave, 'ClipWidth')
        ClipWidth = DataToSave.OtherDataToSave.ClipWidth;
    end
    ShowFocusDotsGlobal = strcmp(Options.ShowFocusDotsAlways, 'y');

    % Four traces at a time, 2x2 -- this tool's own grid, independent of
    % Options.Low/Medium/High (those are baked into ReviewQueue segments
    % at Part C preprocessing time and have nothing to do with this
    % ad hoc, wait-time-filtered pass).
    BatchRows = 2;
    BatchCols = 2;

    while true
        fprintf('\n   --- New wait-time window ---\n');
        RawWindow = strtrim(input(['   Enter wait time window in minutes as "Low High" ' ...
            '(e.g. 12 13), or blank/q to stop: '], 's'));
        if isempty(RawWindow) || strcmpi(RawWindow, 'q')
            break
        end

        WindowVals = str2double(strsplit(RawWindow));
        if numel(WindowVals) ~= 2 || any(isnan(WindowVals)) || WindowVals(1) > WindowVals(2)
            fprintf('   Could not parse -- enter two numbers, low then high, e.g. "12 13".\n');
            continue
        end
        LowMinutes = WindowVals(1);
        HighMinutes = WindowVals(2);

        FilteredIndices = Find_Traces_In_Wait_Window(DataToSave.CombinedAnalyzedTraceData, LowMinutes, HighMinutes);
        fprintf('   Found %d trustworthy "1 Fuse" trace(s) with wait time in [%.4g, %.4g] min.\n', ...
            numel(FilteredIndices), LowMinutes, HighMinutes);
        if isempty(FilteredIndices)
            continue
        end

        [DataToSave, Quit] = Review_Filtered_Batch(FilteredIndices, DataToSave, VideoTimeVector, ClipWidth, ...
            Options, BatchRows, BatchCols, ShowFocusDotsGlobal, UserReviewFilePath, ReviewQueue, ...
            CurrentQueuePosition, ReviewMeta);
        if Quit
            disp(' ')
            disp('   Review-by-wait-time session paused (q). Progress saved.')
            disp('====================================')
            return
        end
    end

    disp(' ')
    disp('   Review-by-wait-time session ended.')
    disp('====================================')
end

function [DataToSave, Quit] = Review_Filtered_Batch(FilteredIndices, DataToSave, VideoTimeVector, ClipWidth, ...
    Options, BatchRows, BatchCols, ShowFocusDotsGlobal, UserReviewFilePath, ReviewQueue, CurrentQueuePosition, ReviewMeta)
% Review_Filtered_Batch  Round-by-round review loop over one wait-time
% window's FilteredIndices -- same shape as Process_Review_Segment /
% Plot_Review_Batch in Start_User_Review.m, simplified: no tiers, no
% H1/H4 special-casing, no 'j'/'s' commands (a narrow wait-time window's
% filtered list doesn't need cross-segment jumping or manual focus-jump
% registration -- those are Low/Medium/High-queue concepts, not needed
% here since this tool has no queue of its own to navigate).

    BatchSize = BatchRows * BatchCols;
    NumFiltered = numel(FilteredIndices);
    Quit = false;

    FigureHandles = Create_Master_Window(BatchRows, BatchCols);
    if strcmp(Options.FixWaitTime, 'y')
        FigureHandles.FixWaitPlot = figure(2);
    end

    BatchStart = 1;
    while BatchStart <= NumFiltered
        BatchEnd = min(BatchStart + BatchSize - 1, NumFiltered);
        CurrentTraceRange = FilteredIndices(BatchStart:BatchEnd);

        fprintf('\n--- Wait-time review -- trace(s) %d-%d of %d matching this window ---\n', ...
            BatchStart, BatchEnd, NumFiltered);

        % Unconditionally clear every subplot slot before drawing, same
        % reasoning as Plot_Review_Batch in Start_User_Review.m -- a
        % smaller final batch would otherwise leave a larger earlier
        % batch's stale plots sitting in the unused slots.
        set(0, 'CurrentFigure', FigureHandles.MasterWindow)
        for k = 1:numel(FigureHandles.SubHandles)
            set(FigureHandles.MasterWindow, 'CurrentAxes', FigureHandles.SubHandles(k));
            cla
            title('')
        end

        for k = 1:numel(CurrentTraceRange)
            GlobalIdx = CurrentTraceRange(k);
            Overlay = [];
            if ShowFocusDotsGlobal
                Overlay = struct('ShowFocusDots', true);
            end
            Plot_Current_Trace(FigureHandles, DataToSave.CombinedAnalyzedTraceData(GlobalIdx), ...
                VideoTimeVector, ClipWidth, k, GlobalIdx, Options, Overlay, []);
        end

        PreviousAnalysisData = DataToSave.CombinedAnalyzedTraceData;
        ErrorCounter = 0;
        RerunThisRound = 'y';
        GoBack = false;

        while strcmp(RerunThisRound, 'y')
            fprintf('\nEnter PlotNumber.Code corrections (comma-separated):\n');
            fprintf('  (blank) = all correct, move to next group\n');
            fprintf('  b       = back to the previous group\n');
            fprintf('  q       = quit and save\n');
            RawInput = strtrim(input('> ', 's'));

            if strcmpi(RawInput, 'q')
                Save_User_Review_File(UserReviewFilePath, DataToSave, ReviewQueue, CurrentQueuePosition, ReviewMeta);
                Quit = true;
                return
            end

            if strcmpi(RawInput, 'b')
                if BatchStart == 1
                    fprintf('Already at the first group of matching traces -- nothing to go back to.\n');
                    continue
                end
                BatchStart = max(1, BatchStart - BatchSize);
                GoBack = true;
                RerunThisRound = 'n';
                continue
            end

            if isempty(RawInput)
                IncorrectPlotIndices = [];
            else
                Tokens = strtrim(strsplit(RawInput, ','));
                Tokens = Tokens(~cellfun(@isempty, Tokens));

                if isempty(Tokens)
                    IncorrectPlotIndices = [];
                else
                    if any(~contains(Tokens, '.'))
                        fprintf(['Each entry must be PlotNumber.Code (e.g. 2.100) -- ' ...
                            '"%s" has no decimal code attached. Try again.\n'], ...
                            strjoin(Tokens(~contains(Tokens, '.')), ', '));
                        RerunThisRound = 'y';
                        continue
                    end

                    IncorrectPlotIndices = str2double(Tokens);
                    if any(isnan(IncorrectPlotIndices))
                        fprintf('Could not parse one or more entries -- try again.\n');
                        RerunThisRound = 'y';
                        continue
                    end
                end
            end

            [RerunThisRound, CorrectedAnalysisData, ErrorCounter] = Correct_Designations(IncorrectPlotIndices, ...
                PreviousAnalysisData, CurrentTraceRange, DataToSave.CombinedAnalyzedTraceData, ErrorCounter, ...
                Options, VideoTimeVector, FigureHandles, ClipWidth, []);

            DataToSave.CombinedAnalyzedTraceData = CorrectedAnalysisData;
        end

        if GoBack
            continue
        end

        BatchStart = BatchEnd + 1;
        Save_User_Review_File(UserReviewFilePath, DataToSave, ReviewQueue, CurrentQueuePosition, ReviewMeta);
    end

    fprintf('\n   All %d matching trace(s) in this window reviewed.\n', NumFiltered);
end

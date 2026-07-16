function Start_User_Review(SelectedFilePath, MasterAnalysisFolder)
%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
% Start_User_Review  Bootstrap entry point for Part C -- User Review Traces.
%
% SelectedFilePath is either a raw Part B output .mat file (first run) or
% an existing User Review File (resume) -- picked by
% An_Easy_Start_UserReview.m and handed in directly.
%
% MasterAnalysisFolder is An_Easy_Start_UserReview.m's DataLocation_UserReview
% -- the shared folder (across every experiment/dataset) where
% Key_Traces_Collection.mat and Performance_Log.mat live and accumulate.
% See Extract_Future_Key_Traces.m / Update_Performance_Log.m.
%
% Bootstraps/resumes the User Review File, prints a segment summary,
% then runs the round-by-round review loop for the Low and Medium tiers
% plus all six High subgroups (H1, H2a, H2b, H3, H4, H5).
% See the design notes on the Part C roadmap.

disp('====================================')
disp('Initiating User Review Trace Program.......')

close all

Options = Setup_Options_User_Review();
Options
fprintf('   Low:         Rows = %d, Cols = %d\n', Options.Low.Rows,    Options.Low.Cols);
fprintf('   Medium:      Rows = %d, Cols = %d\n', Options.Medium.Rows, Options.Medium.Cols);
fprintf('   High (H2a-H5): Rows = %d, Cols = %d\n', Options.High.Rows, Options.High.Cols);
if isfield(Options, 'HighOverrides') && isfield(Options.HighOverrides, 'H1')
    fprintf('   High/H1:     Rows = %d, Cols = %d, TracesPerBatch = %d\n', Options.HighOverrides.H1.Rows, ...
        Options.HighOverrides.H1.Cols, Options.HighOverrides.H1.TracesPerBatch);
end
AreOptionsOk = strtrim(input('   Options look good? (Enter or y to proceed, n to cancel): ','s'));
if strcmpi(AreOptionsOk,'n')
    disp('   Options no good? Then change em and run the program again!')
    disp('   Program terminated.')
    disp('====================================')
    return
end

if strcmp(Options.ReviewOldDataMode, 'y')
    Options.ReportPerformanceLog = 'n';
    Options.ExtractKeyTraces = 'n';
    disp('   Review Old Data Mode: auto-disabling ReportPerformanceLog/ExtractKeyTraces')
    disp('   so this legacy-dataset session does not write into the shared master tracking files.')
end

UserReviewFilePath = Get_User_Review_File_Path(SelectedFilePath, Options);
IsResuming = isfile(UserReviewFilePath);

if IsResuming

    disp('   Found existing review file -- resuming.')
    disp(['   ' UserReviewFilePath])

    Loaded = load(UserReviewFilePath, 'DataToSave', 'ReviewQueue', 'CurrentQueuePosition', 'ReviewMeta');
    DataToSave           = Loaded.DataToSave;
    ReviewQueue          = Loaded.ReviewQueue;
    CurrentQueuePosition = Loaded.CurrentQueuePosition;
    ReviewMeta           = Loaded.ReviewMeta;

    LiveOptionsUsed = struct('ReviewFolderName', Options.ReviewFolderName, ...
        'Low', Options.Low, 'Medium', Options.Medium, 'High', Options.High, ...
        'HighOverrides', Get_High_Overrides(Options));
    if ~isequal(LiveOptionsUsed, ReviewMeta.PreprocessingOptionsUsed)
        disp('   NOTE: current Options.ReviewFolderName/Low/Medium/High/HighOverrides differ from')
        disp('   what this review file was preprocessed with.')
    end

    % Backfill for a User Review File saved before FinalHarvestComplete
    % existed -- same getfield_or_empty-style tolerance Refresh_Segment_Grids
    % already uses for legacy segment structs, rather than requiring a fresh
    % re-preprocess just to gain this one flag.
    if ~isfield(ReviewMeta, 'FinalHarvestComplete')
        ReviewMeta.FinalHarvestComplete = false;
    end
    if ~isfield(ReviewMeta, 'ManualFocusSubtractIndices')
        ReviewMeta.ManualFocusSubtractIndices = [];
    end

    % Grid/layout metadata (Rows/Cols/TracesPerBatch/NumBatches) is purely a
    % display/batching concern -- it has zero bearing on which traces belong
    % to a segment (StartIndex/EndIndex/NumTraces/Tier/SubgroupName) or where
    % CurrentQueuePosition points, so it's always safe to re-derive from the
    % CURRENT Options on every resume, even if trace bucketing/order was
    % frozen at first-run preprocessing time. This is what actually lets an
    % in-progress review file pick up a newly added/changed grid (e.g. H1's
    % override) without losing any review progress -- previously this was
    % only ever a passive warning, which left an already-baked segment
    % permanently stuck on its old grid.
    [ReviewQueue, GridsChanged] = Refresh_Segment_Grids(ReviewQueue, Options);
    if GridsChanged
        disp('   Grid/layout settings (Rows/Cols/batch size) refreshed from current Options for')
        disp('   one or more segments -- trace bucketing/order and review progress are untouched.')
    end

else

    disp('   No existing review file found -- running preprocessing...')

    Raw = load(SelectedFilePath, 'DataToSave');
    CombinedAnalyzedTraceData = Raw.DataToSave.CombinedAnalyzedTraceData;

    [SortedTraceData, ReviewQueue] = Build_Review_Queue(CombinedAnalyzedTraceData, Options);
    CurrentQueuePosition = 1;

    DataToSave = struct();
    DataToSave.CombinedAnalyzedTraceData = SortedTraceData;
    DataToSave.OtherDataToSave           = Raw.DataToSave.OtherDataToSave;
    DataToSave.ReviewOptions             = Options;

    ReviewMeta = struct();
    ReviewMeta.SourceFilePath           = SelectedFilePath;
    ReviewMeta.CreatedTimestamp         = char(datetime('now'));
    ReviewMeta.SchemaVersion            = 1;
    ReviewMeta.PreprocessingOptionsUsed = struct('ReviewFolderName', Options.ReviewFolderName, ...
        'Low', Options.Low, 'Medium', Options.Medium, 'High', Options.High, ...
        'HighOverrides', Get_High_Overrides(Options));
    ReviewMeta.FinalHarvestComplete = false;
    ReviewMeta.ManualFocusSubtractIndices = [];

    Save_User_Review_File(UserReviewFilePath, DataToSave, ReviewQueue, CurrentQueuePosition, ReviewMeta);

    disp('   Preprocessing complete. Review file created:')
    disp(['   ' UserReviewFilePath])

end

% ---- Segment summary ----
NumTraces = numel(DataToSave.CombinedAnalyzedTraceData);
fprintf('\n   %d traces total. Queue position: %d of %d.\n\n', NumTraces, CurrentQueuePosition, NumTraces);
fprintf('   %-14s %-14s %10s\n', 'Tier', 'Subgroup', 'NumTraces');
fprintf('   %-14s %-14s %10s\n', repmat('-',1,14), repmat('-',1,14), repmat('-',1,10));
for k = 1:numel(ReviewQueue.Segments)
    seg = ReviewQueue.Segments(k);
    fprintf('   %-14s %-14s %10d\n', seg.Tier, seg.SubgroupName, seg.NumTraces);
end

% ---- Where to start this session ----
% Offered on BOTH a fresh run and a resume -- a reviewer may want to jump
% straight to a specific tier (L/M/H) or trace number even on the very
% first run, not just when picking back up later.
CurrentQueuePosition = Prompt_Start_Position(ReviewQueue, CurrentQueuePosition, NumTraces, IsResuming);

% ---- Round-by-round review loop -- Low, Medium, and all of High. ----
% Process_Review_Segment is driven entirely by the Segment row it's given,
% so wiring up Medium here was just looping it over more segments, not
% rewriting it -- see the design notes' roadmap, Step 4. Every High subgroup now has
% its own diagnostic overlay (Plot_Trace_Pair_H1.m for H1;
% Plot_Current_Trace.m's DiagnosticOverlay.ShowFocusDots for H2a/H2b/H3/H5;
% DiagnosticOverlay.TitleSuffix via Get_H4_Failed_Test_Label.m for H4).
% ismember (not strcmp) preserves ReviewQueue.Segments' existing physical
% order -- Low, then Medium's 1 Fuse -> 2 Fuse -> Unbound subgroups.
%
% OtherDataToSave.UniversalData no longer exists (PART_B_HANDOFF_NOTES.md
% section 1a) -- VideoTimeVector is its replacement for bind-time
% computation (Fix_Fusion_Wait_Time.m/Fix_Unbind_Wait_Time.m). ClipWidth
% (section 5) is a new passthrough field, soft-guarded since it's only
% needed by optional display enhancements (the H2a/H2b focus-jump
% correction and the running-median dead-zone crop, both in
% Plot_Current_Trace.m) rather than a core workflow -- missing ClipWidth
% just means those features quietly stay off, not an error.
VideoTimeVector = DataToSave.OtherDataToSave.VideoTimeVector;
ClipWidth = [];
if isfield(DataToSave.OtherDataToSave, 'ClipWidth')
    ClipWidth = DataToSave.OtherDataToSave.ClipWidth;
end

% Manual focus-jump correction ('s' round command, Process_Review_Segment)
% needs a dataset-wide global<->clipped conversion for registration-time
% validation: GlobalFindingImageFrame is numerically identical for every
% trace in this file (one shared finding-image scan locates every particle
% at once, not one per virion), so Global = clippedIdx + GlobalFindingImageFrame
% is a single, dataset-wide conversion, not a per-trace one.
% GlobalFocusFrameNumbers is Part A's own record of every raw/global frame
% a focus/refocus event occurred (Find_And_Analyze_Particles.m ->
% OtherDataToSave.GlobalFocusFrameNumbers = Options.FocusFrameNumbers).
% Both soft-guarded the same way ClipWidth is above -- a dataset that
% predates either field just quietly disables the 's' command, not an error.
GlobalFindingImageFrame = [];
if isfield(DataToSave.OtherDataToSave, 'GlobalFindingImageFrame')
    GlobalFindingImageFrame = DataToSave.OtherDataToSave.GlobalFindingImageFrame;
end
GlobalFocusFrameNumbers = [];
if isfield(DataToSave.OtherDataToSave, 'GlobalFocusFrameNumbers')
    GlobalFocusFrameNumbers = DataToSave.OtherDataToSave.GlobalFocusFrameNumbers;
end
if isempty(GlobalFindingImageFrame) || isempty(GlobalFocusFrameNumbers)
    disp('   NOTE: GlobalFindingImageFrame/GlobalFocusFrameNumbers not found in this dataset')
    disp('   (older Part A/B output) -- the manual focus-jump correction command (s) will be')
    disp('   unavailable this session.')
end

% Session-wide, growing, persisted list -- seeded from ReviewMeta (already
% backfilled to [] above if this is an older file), carried through every
% Process_Review_Segment call below, and written back into ReviewMeta
% before every subsequent Save_User_Review_File call.
ManualFocusSubtractIndices = ReviewMeta.ManualFocusSubtractIndices;

ActiveSegments = ReviewQueue.Segments(ismember({ReviewQueue.Segments.Tier}, {'Low','Medium'}) | ...
    ismember({ReviewQueue.Segments.SubgroupName}, {'H1','H2a','H2b','H3','H4','H5'}));

% while, not for -- 'b' at the very first batch of a segment can cross
% back into the previous segment's last batch (see Process_Review_Segment),
% which needs s to move BACKWARD, not just forward.
s = 1;
while s >= 1 && s <= numel(ActiveSegments)
    Segment = ActiveSegments(s);

    if CurrentQueuePosition > Segment.EndIndex
        s = s + 1;   % already fully reviewed in a prior session
        continue
    end
    if CurrentQueuePosition < Segment.StartIndex
        s = s - 1;   % a cross-segment 'b' moved us into an earlier segment
        continue
    end

    [DataToSave, CurrentQueuePosition, Quit, ManualFocusSubtractIndices] = Process_Review_Segment(ActiveSegments, s, DataToSave, ...
        VideoTimeVector, ClipWidth, Options, UserReviewFilePath, ReviewQueue, ReviewMeta, CurrentQueuePosition, ...
        GlobalFindingImageFrame, GlobalFocusFrameNumbers, ManualFocusSubtractIndices);

    if Quit
        disp(' ')
        disp('   Review session paused (q). Progress saved -- re-run to resume where you left off.')
        disp('====================================')
        return
    end
end

disp(' ')
disp('   Review complete -- Low, Medium, and all of High (H1-H5).')

% ---- Final key-trace harvest + performance-log update ----
% Runs only on this normal-completion path, never on the Quit early-return
% above -- an interrupted session hasn't finished being reviewed, so neither
% a 20% sample nor a full performance tally would be meaningful yet.
% ReviewMeta.FinalHarvestComplete guards against re-harvesting a second 20%
% sample (or double-counting this session's TP/FN/FP into the cumulative
% log) if Start_User_Review is ever re-run against an already-fully-reviewed
% User Review File.
if ~ReviewMeta.FinalHarvestComplete
    % One last confirmation before the harvest (and the cross-session
    % Performance_Log.mat/Key_Traces_Collection.mat writes it triggers) --
    % gives the reviewer a chance to back out and keep reviewing/correcting
    % instead. Blank Enter or 'y' proceeds, only 'n' holds off, matching this
    % file's existing Options-confirmation convention above. Declining just
    % ends the session without harvesting -- CurrentQueuePosition is already
    % past every segment, so re-running lands right back at this same
    % prompt, askable again as many times as needed.
    ReadyToFinish = strtrim(input(['   Have you reviewed all traces and are ready to finish? ' ...
        '(Enter or y to proceed, n to hold off): '], 's'));
    if strcmpi(ReadyToFinish, 'n')
        disp('   Okay -- holding off. Re-run to pick up here again whenever you are ready to finish.')
        disp('====================================')
        return
    end

    disp(' ')
    disp('   Calculating performance and preparing to show results -- please wait, do not quit...')
    if ~strcmp(Options.ReportPerformanceLog, 'y') && ~strcmp(Options.ExtractKeyTraces, 'y')
        disp('   (Note: ReportPerformanceLog and ExtractKeyTraces are both set to ''n''.')
        disp('   Skipping performance-log reporting and key-trace extraction this session.)')
    elseif ~strcmp(Options.ReportPerformanceLog, 'y')
        disp('   (Note: Options.ReportPerformanceLog is set to ''n''.')
        disp('   Skipping performance-log reporting this session.)')
    elseif ~strcmp(Options.ExtractKeyTraces, 'y')
        disp('   (Note: Options.ExtractKeyTraces is set to ''n''.')
        disp('   Skipping key-trace extraction this session.)')
    end

    EligibleMask = arrayfun(@(t) ~strcmp(t.AlgoDesignation, 'Other') && ~strcmp(t.Designation, 'Other'), ...
        DataToSave.CombinedAnalyzedTraceData);
    EligibleTraces = DataToSave.CombinedAnalyzedTraceData(EligibleMask);

    if strcmp(Options.ReportPerformanceLog, 'y')
        Update_Performance_Log(EligibleTraces, MasterAnalysisFolder, ReviewMeta.SourceFilePath);
    end
    if strcmp(Options.ExtractKeyTraces, 'y')
        Extract_Future_Key_Traces(EligibleTraces, DataToSave.OtherDataToSave, MasterAnalysisFolder);
    end

    ReviewMeta.FinalHarvestComplete = true;
    ReviewMeta.ManualFocusSubtractIndices = ManualFocusSubtractIndices;
    Save_User_Review_File(UserReviewFilePath, DataToSave, ReviewQueue, CurrentQueuePosition, ReviewMeta);
else
    disp(' ')
    disp('   Final key/performance harvest already completed for this dataset -- skipping to avoid double-counting.')
end

% ---- Final review summary: designation tally + finding-image virus map ----
% Runs only on this normal-completion path (every segment reached its end),
% never on the Quit early-return above -- an interrupted session still has
% unreviewed traces, so a full tally wouldn't be meaningful yet. Mirrors
% Part B's own end-of-analysis step (Analyze_Trace_Data.m ->
% Display_All_Analyzed_Viruses.m), but against DataToSave.CombinedAnalyzedTraceData
% -- the reviewer-corrected data, not the raw algorithm output.
disp(' ')
disp('   --- Final designation summary ---')
Print_Designation_Summary(DataToSave.CombinedAnalyzedTraceData)

if isfield(DataToSave.OtherDataToSave, 'FindingImage')
    minImageShow = 400;   % Part A defaults (Setup_Options_Extract_Traces.m),
    maxImageShow = 800;   % used as fallback if OtherDataToSave.Options is absent
    if isfield(DataToSave.OtherDataToSave, 'Options')
        if isfield(DataToSave.OtherDataToSave.Options, 'MinImageShow')
            minImageShow = DataToSave.OtherDataToSave.Options.MinImageShow;
        end
        if isfield(DataToSave.OtherDataToSave.Options, 'MaxImageShow')
            maxImageShow = DataToSave.OtherDataToSave.Options.MaxImageShow;
        end
    end
    Display_Reviewed_Viruses(DataToSave.OtherDataToSave.FindingImage, DataToSave.CombinedAnalyzedTraceData, ...
        minImageShow, maxImageShow, DataToSave.OtherDataToSave.GlobalFindingImageFrame);
else
    fprintf('   WARNING: OtherDataToSave.FindingImage not found -- skipping final virus overlay figure.\n');
end

disp('====================================')

end

function [DataToSave, CurrentQueuePosition, Quit, ManualFocusSubtractIndices] = Process_Review_Segment(ActiveSegments, SegmentIdx, DataToSave, ...
    VideoTimeVector, ClipWidth, Options, UserReviewFilePath, ReviewQueue, ReviewMeta, CurrentQueuePosition, ...
    GlobalFindingImageFrame, GlobalFocusFrameNumbers, ManualFocusSubtractIndices)
% Process_Review_Segment  Round-by-round review loop for ONE ReviewQueue
% segment (ActiveSegments(SegmentIdx)). Tier-agnostic -- driven entirely by
% the Segment row's own Tier/SubgroupName/StartIndex/EndIndex/Rows/Cols;
% carries no notion of "Low" itself, so Low/Medium/High all share this
% same loop body. Takes the full ActiveSegments array (not just its own
% row) so 'b' at the very first batch of a segment can look BACKWARD past
% it -- into the previous non-empty segment's last batch -- which a
% single Segment row has no way to know about on its own.
%
% GlobalFindingImageFrame/GlobalFocusFrameNumbers (both [] if this dataset
% predates them) and ManualFocusSubtractIndices are for the 's' round
% command (manual focus-jump correction, see Plot_Current_Trace.m /
% Apply_Manual_Focus_Corrections.m) -- ManualFocusSubtractIndices is now a
% genuine OUTPUT of this function too (grows via 's'), unlike every other
% input here, which only ever flows in.

    Segment = ActiveSegments(SegmentIdx);
    Quit = false;
    BatchSize = Get_Segment_Batch_Size(Segment);
    NumTraces = numel(DataToSave.CombinedAnalyzedTraceData);   % for the 'j' jump prompt below
    ManualFocusAvailable = ~isempty(GlobalFindingImageFrame) && ~isempty(GlobalFocusFrameNumbers);

    FigureHandles = Create_Master_Window(Segment.Rows, Segment.Cols);
    if strcmp(Options.FixWaitTime,'y')
        FigureHandles.FixWaitPlot = figure(2);
    end

    % One caption for the whole grid, pinned to the BOTTOM of the figure --
    % text doesn't change round to round within a segment, so set it once
    % here rather than in the per-round loop. An annotation textbox (not
    % sgtitle -- which sits at the top and collided with the top row of
    % subplot titles) is a figure child, not an axes, so it's unaffected by
    % the per-round unconditional cla() over FigureHandles.SubHandles, and
    % Create_Master_Window's per-segment clf() clears any leftover one from
    % a previous segment automatically.
    annotation(FigureHandles.MasterWindow, 'textbox', [0 0 1 0.035], ...
        'String', Get_Segment_Title_Text(Segment), 'EdgeColor', 'none', 'FitBoxToText', 'off', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 11);

    % A User Review File preprocessed before the H1 two-pane layout existed
    % has its H1 segment (if any) baked with the old shared grid (no
    % TracesPerBatch field, Rows != 2) -- Plot_Trace_Pair_H1 requires the
    % new shape, so fall back to the plain single-pane view for a segment
    % that doesn't have it rather than erroring out and losing review
    % progress. Re-run against a fresh User Review File to pick up the new
    % layout. Only Rows==2 is checked (not Cols, which is a tunable "traces
    % per round" knob for this design, not a fixed shape requirement).
    IsH1PairLayout = strcmp(Segment.SubgroupName, 'H1') && isfield(Segment, 'TracesPerBatch') ...
        && Segment.Rows == 2;
    if strcmp(Segment.SubgroupName, 'H1') && ~IsH1PairLayout
        fprintf(['NOTE: this H1 segment was preprocessed before the two-pane H1 layout existed ' ...
            '(grid %dx%d) -- showing the plain single-pane view for these traces. Re-run against ' ...
            'a fresh User Review File to get the full zoom-pane layout.\n'], Segment.Rows, Segment.Cols);
    end

    % H2a/H2b/H3/H5: same single-pane Plot_Current_Trace layout as
    % Low/Medium, just with every focus event drawn on top (red/black/blue
    % dot triple, same marker Fix_Fusion_Wait_Time.m/Fix_Unbind_Wait_Time.m
    % already use in the blown-up picker) -- see Plot_Focus_Dot_Markers.m.
    % No pane-layout change needed here, unlike H1, so this is just an
    % extra DiagnosticOverlay argument on the ordinary call. H4 gets the
    % same focus-dot overlay too, but needs its own branch below since it
    % also carries a per-trace TitleSuffix -- deliberately excluded here.
    %
    % Options.ShowFocusDotsAlways ('y' by default) extends the dot overlay
    % to every other segment too (Low, Medium, H1's two panes) -- with it
    % on, NeedsFocusDotOverlay is true unconditionally and the plain 'else'
    % branch below never actually runs; it's kept only as the fallback for
    % when a reviewer turns the global option off.
    ShowFocusDotsGlobal = strcmp(Options.ShowFocusDotsAlways, 'y');
    NeedsFocusDotOverlay = ShowFocusDotsGlobal || ismember(Segment.SubgroupName, {'H2a', 'H2b', 'H3', 'H5'});

    % H2a/H2b: subtract the focus-induced jump near the candidate/actual
    % fusion frame so the reviewer sees the trace as if that jump hadn't
    % happened (Get_H2_Focus_Jump_Trigger_Frame.m, applied inside
    % Plot_Current_Trace.m). Requires ClipWidth to convert Part B's
    % Review_PriorityData coordinates -- on a dataset that predates that
    % field, this quietly declines rather than erroring, same as every
    % other optional overlay's graceful-degrade convention in this file.
    ApplyH2Correction = ismember(Segment.SubgroupName, {'H2a', 'H2b'}) ...
        && strcmp(Options.ApplyH2FocusJumpCorrection, 'y') && ~isempty(ClipWidth);

    while CurrentQueuePosition <= Segment.EndIndex
        BatchStart = CurrentQueuePosition;
        BatchEnd   = min(BatchStart + BatchSize - 1, Segment.EndIndex);
        CurrentTraceRange = BatchStart:BatchEnd;

        % Brief label, NOT Get_Segment_Title_Text (that's the long
        % descriptive caption for the figure annotation) -- Low is just
        % "Low", Medium adds its subgroup ("Medium: 1 Fuse"), High's
        % SubgroupName is already a short code ("H1", "H2a", ...) on its
        % own.
        if strcmp(Segment.Tier, 'Low')
            ShortLabel = 'Low';
        elseif strcmp(Segment.Tier, 'Medium')
            ShortLabel = sprintf('Medium: %s', Segment.SubgroupName);
        else
            ShortLabel = Segment.SubgroupName;
        end
        fprintf('\n--- %s -- traces %d-%d (%d-%d / %d) ---\n', ...
            ShortLabel, BatchStart, BatchEnd, Segment.StartIndex, Segment.EndIndex, NumTraces);

        % Clearing + drawing the whole batch is factored into Plot_Review_Batch
        % (below) so it can be called a SECOND time, mid-round, right after
        % the 's' command registers a new manual focus-jump index -- see the
        % 's' branch further down. Every unconditional-clear/gcf-reset
        % reasoning documented there still applies on every call, not just
        % this first one.
        Plot_Review_Batch(FigureHandles, Segment, DataToSave, CurrentTraceRange, VideoTimeVector, ClipWidth, ...
            Options, IsH1PairLayout, ShowFocusDotsGlobal, NeedsFocusDotOverlay, ApplyH2Correction, ManualFocusSubtractIndices);

        PreviousAnalysisData = DataToSave.CombinedAnalyzedTraceData;
        ErrorCounter = 0;
        RerunThisRound = 'y';
        GoBack = false;

        while strcmp(RerunThisRound,'y')
            fprintf('\nEnter PlotNumber.Code corrections (comma-separated):\n');
            fprintf('  (blank) = all correct, move to next group\n');
            fprintf('  b       = back to the previous group\n');
            fprintf('  j       = jump to a trace number or section\n');
            if ManualFocusAvailable
                fprintf('  s       = subtract a focus jump from all traces\n');
            end
            fprintf('  q       = quit and save\n');
            RawInput = strtrim(input('> ','s'));

            if strcmpi(RawInput,'q')
                ReviewMeta.ManualFocusSubtractIndices = ManualFocusSubtractIndices;
                Save_User_Review_File(UserReviewFilePath, DataToSave, ReviewQueue, CurrentQueuePosition, ReviewMeta);
                Quit = true;
                return
            end

            if strcmpi(RawInput,'b')
                PrevBatchStart = max(Segment.StartIndex, BatchStart - BatchSize);
                if PrevBatchStart == BatchStart
                    % Already at this segment's first batch -- try crossing
                    % back into the previous non-empty segment's last batch
                    % instead of refusing outright.
                    PrevSeg = Find_Previous_Nonempty_Segment(ActiveSegments, SegmentIdx);
                    if isempty(PrevSeg)
                        fprintf('Already at the very first group of traces -- nothing to go back to.\n');
                        continue
                    end
                    PrevSegBatchSize = Get_Segment_Batch_Size(PrevSeg);
                    CurrentQueuePosition = PrevSeg.StartIndex + (PrevSeg.NumBatches - 1) * PrevSegBatchSize;
                    fprintf('Crossing back into %s / %s (last group of that segment).\n', PrevSeg.Tier, PrevSeg.SubgroupName);
                    % Not saved, same as within-segment 'b' -- a crash right
                    % after this resumes from the last real save, not here.
                    return   % outer loop steps SegmentIdx back down to find it
                end
                CurrentQueuePosition = PrevBatchStart;
                GoBack = true;
                RerunThisRound = 'n';
                continue
            end

            if strcmpi(RawInput,'j')
                NewPosition = Prompt_Jump_Target(ReviewQueue, NumTraces, '   Jump to:');
                if isempty(NewPosition)
                    % Cancelled (blank entry inside the jump prompt) --
                    % redo this same round prompt, nothing has changed.
                    continue
                end
                CurrentQueuePosition = NewPosition;
                % Not saved -- same reasoning as 'b' above: a jump is a
                % navigational move, not a completed round. Always return,
                % even for a target inside THIS segment, so the outer
                % while(s) loop in Start_User_Review.m re-locates the
                % right segment and re-enters fresh (new figure, re-printed
                % caption) -- the same path 'b' already uses when it
                % crosses a segment boundary, reused here unconditionally
                % rather than adding a separate same-segment fast path.
                return
            end

            if strcmpi(RawInput,'s')
                if ~ManualFocusAvailable
                    fprintf(['Manual focus-jump correction is unavailable for this dataset -- ' ...
                        'GlobalFindingImageFrame/GlobalFocusFrameNumbers not found (see note at session start).\n']);
                    continue
                end

                [ManualFocusSubtractIndices, AnyRegistered] = Prompt_Manual_Focus_Correction( ...
                    ManualFocusSubtractIndices, GlobalFindingImageFrame, GlobalFocusFrameNumbers);

                if AnyRegistered
                    % Meaningful reviewer work, not a pure navigational move
                    % like 'b'/'j' -- saved immediately (also carried
                    % forward into every later normal save via ReviewMeta).
                    ReviewMeta.ManualFocusSubtractIndices = ManualFocusSubtractIndices;
                    Save_User_Review_File(UserReviewFilePath, DataToSave, ReviewQueue, CurrentQueuePosition, ReviewMeta);

                    % Re-draw the CURRENT batch immediately so the reviewer
                    % sees the correction applied right away, then fall back
                    % into this same round's prompt -- RerunThisRound is
                    % untouched (still 'y'), so the while test above loops
                    % straight back to it.
                    Plot_Review_Batch(FigureHandles, Segment, DataToSave, CurrentTraceRange, VideoTimeVector, ClipWidth, ...
                        Options, IsH1PairLayout, ShowFocusDotsGlobal, NeedsFocusDotOverlay, ApplyH2Correction, ...
                        ManualFocusSubtractIndices);
                end
                continue
            end

            if isempty(RawInput)
                IncorrectPlotIndices = [];
            else
                Tokens = strtrim(strsplit(RawInput, ','));

                % Drop empty tokens from a trailing/leading/doubled comma
                % (e.g. "5.100, 9.100," -- typing a hanging comma after the
                % last entry, so the next code can just be typed after it
                % while looking at the plots, is a normal habit this should
                % tolerate rather than reject) -- strsplit(',') on a
                % trailing comma produces a trailing '' token that would
                % otherwise fail the '.' check below.
                Tokens = Tokens(~cellfun(@isempty, Tokens));

                if isempty(Tokens)
                    % Input was comma(s) only (e.g. "," or ",,,") -- same as
                    % a blank entry.
                    IncorrectPlotIndices = [];
                else
                    % Require an actual PlotNumber.Code format -- checked on
                    % the RAW STRING, not the parsed value, since
                    % str2double('1') and str2double('1.0') are numerically
                    % identical. Without this, a bare plot number typo (no
                    % code at all) would silently parse as ".0" (rem(1,1)==0
                    % -> MilliCode 0 -> 'No Fusion') and get applied with
                    % zero confirmation.
                    if any(~contains(Tokens, '.'))
                        fprintf(['Each entry must be PlotNumber.Code (e.g. 5.100) -- ' ...
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
                Options, VideoTimeVector, FigureHandles, ClipWidth, ManualFocusSubtractIndices);

            DataToSave.CombinedAnalyzedTraceData = CorrectedAnalysisData;
        end

        if GoBack
            % Jump back one group -- don't advance CurrentQueuePosition.
            % The outer while re-plots starting at the earlier position
            % using the current (already-corrected) DataToSave, so any
            % fixes made earlier show up live: new designation, new
            % color, new fuse-frame line, etc.
            continue
        end

        CurrentQueuePosition = BatchEnd + 1;
        ReviewMeta.ManualFocusSubtractIndices = ManualFocusSubtractIndices;
        Save_User_Review_File(UserReviewFilePath, DataToSave, ReviewQueue, CurrentQueuePosition, ReviewMeta);
    end
end

function Plot_Review_Batch(FigureHandles, Segment, DataToSave, CurrentTraceRange, VideoTimeVector, ClipWidth, ...
    Options, IsH1PairLayout, ShowFocusDotsGlobal, NeedsFocusDotOverlay, ApplyH2Correction, ManualFocusSubtractIndices)
% Plot_Review_Batch  Draws CurrentTraceRange into FigureHandles' subplot
% grid -- the per-round drawing pass Process_Review_Segment's outer while
% loop always ran exactly once at the top of each batch, now also callable
% a SECOND time, mid-round, right after the 's' command registers a new
% ManualFocusSubtractIndices entry, so the reviewer sees the correction
% applied to the CURRENT batch immediately. Every parameter here is
% something the caller already has in scope -- MATLAB local functions
% don't close over caller variables, so nothing here can be implicitly
% picked up; it all has to come in as an argument.

    % Unconditionally clear EVERY subplot slot before drawing anything,
    % every single call -- not just the slots this round happens to need.
    % Create_Master_Window's clf only runs once per SEGMENT, not once per
    % BATCH (or per re-draw), so without this, a partial/smaller batch
    % would leave whatever a larger earlier batch in this same segment drew
    % still sitting in the unused slots. Clearing everything up front,
    % unconditionally, removes any dependency on correctly computing "which
    % slots are left over" -- there is nothing left to get wrong.
    % set(0,'CurrentFigure',...) is required here, not just
    % set(MasterWindow,'CurrentAxes',...) -- bare cla/title('') below
    % resolve via gca, which depends on the SESSION-WIDE current figure
    % (gcf), not on which figure's CurrentAxes property was last set.
    % Without this line, if Fix_Fusion_Wait_Time.m/Fix_Unbind_Wait_Time.m
    % opened FixWaitPlot (figure 2) during this round's corrections, gcf is
    % left pointing at figure 2 and never reset -- so this entire loop
    % silently clears figure 2's axes instead, leaving any subplot slot the
    % next draw doesn't happen to redraw stuck showing stale content from a
    % previous batch.
    set(0,'CurrentFigure',FigureHandles.MasterWindow)
    for k = 1:numel(FigureHandles.SubHandles)
        set(FigureHandles.MasterWindow, 'CurrentAxes', FigureHandles.SubHandles(k));
        cla
        title('')
    end

    for k = 1:numel(CurrentTraceRange)
        GlobalIdx = CurrentTraceRange(k);
        if IsH1PairLayout
            H1Overlay = [];
            if ShowFocusDotsGlobal
                H1Overlay = struct('ShowFocusDots', true);
            end
            Plot_Trace_Pair_H1(FigureHandles, DataToSave.CombinedAnalyzedTraceData(GlobalIdx), ...
                VideoTimeVector, ClipWidth, k, GlobalIdx, Options, Segment, H1Overlay, ManualFocusSubtractIndices);
        elseif strcmp(Segment.SubgroupName, 'H4')
            Overlay = struct('ShowFocusDots', true, ...
                'TitleSuffix', Get_H4_Failed_Test_Label(DataToSave.CombinedAnalyzedTraceData(GlobalIdx)));
            Plot_Current_Trace(FigureHandles, DataToSave.CombinedAnalyzedTraceData(GlobalIdx), ...
                VideoTimeVector, ClipWidth, k, GlobalIdx, Options, Overlay, ManualFocusSubtractIndices);
        elseif ApplyH2Correction
            Overlay = struct('ShowFocusDots', true, 'ApplyFocusJumpCorrection', true);
            Plot_Current_Trace(FigureHandles, DataToSave.CombinedAnalyzedTraceData(GlobalIdx), ...
                VideoTimeVector, ClipWidth, k, GlobalIdx, Options, Overlay, ManualFocusSubtractIndices);
        elseif NeedsFocusDotOverlay
            Plot_Current_Trace(FigureHandles, DataToSave.CombinedAnalyzedTraceData(GlobalIdx), ...
                VideoTimeVector, ClipWidth, k, GlobalIdx, Options, struct('ShowFocusDots', true), ManualFocusSubtractIndices);
        else
            Plot_Current_Trace(FigureHandles, DataToSave.CombinedAnalyzedTraceData(GlobalIdx), ...
                VideoTimeVector, ClipWidth, k, GlobalIdx, Options, [], ManualFocusSubtractIndices);
        end
    end
end

function CurrentQueuePosition = Prompt_Start_Position(ReviewQueue, CurrentQueuePosition, NumTraces, IsResuming)
% Prompt_Start_Position  Ask the reviewer where to start this session --
% offered both the first time a dataset is ever opened (right after
% preprocessing, CurrentQueuePosition == 1) and on every later resume, so
% jumping straight to a section or trace number doesn't require having
% reviewed anything yet.
%
% Blank keeps CurrentQueuePosition as-is (trace 1 on a fresh run, wherever
% left off when resuming). Anything else is handled by the shared
% Prompt_Jump_Target/Resolve_Jump_Target helpers below -- also used by
% Process_Review_Segment's in-round 'j' command, so both places accept
% the exact same jump codes and can never drift out of sync.

    Segment = Find_Segment_For_Position(ReviewQueue, CurrentQueuePosition);
    if isempty(Segment)
        % Past every segment -- the whole dataset has already been reviewed.
        % Still offer the jump prompt rather than silently returning: a
        % reviewer needs a way to go back and revisit/correct specific
        % traces even after finishing (the final-harvest guard already
        % prevents a completed dataset from being double-harvested, so
        % going back in to make more corrections is always safe). Blank
        % Enter keeps CurrentQueuePosition as-is and proceeds straight
        % through to the completion/harvest step exactly as before.
        fprintf('\n   This entire dataset (%d traces) has already been reviewed.\n', NumTraces);
        NewPosition = Prompt_Jump_Target(ReviewQueue, NumTraces, ...
            '   Press Enter to finish, or jump to a trace/section to review again:');
        if ~isempty(NewPosition)
            CurrentQueuePosition = NewPosition;
        end
        return
    end

    if IsResuming
        PositionInSegment = CurrentQueuePosition - Segment.StartIndex + 1;
        fprintf('\n   You left off at trace %d of %d in %s / %s.\n', ...
            PositionInSegment, Segment.NumTraces, Segment.Tier, Segment.SubgroupName);
        fprintf('   (overall queue position %d of %d)\n', CurrentQueuePosition, NumTraces);
        HeaderText = '   Press Enter to continue from there, or jump to:';
    else
        fprintf('\n   Ready to begin at trace 1 of %d (%s / %s).\n', ...
            NumTraces, Segment.Tier, Segment.SubgroupName);
        HeaderText = '   Press Enter to start from the beginning, or jump to:';
    end

    NewPosition = Prompt_Jump_Target(ReviewQueue, NumTraces, HeaderText);
    if ~isempty(NewPosition)
        CurrentQueuePosition = NewPosition;
    end
end

function NewPosition = Prompt_Jump_Target(ReviewQueue, NumTraces, HeaderText)
% Prompt_Jump_Target  Interactive jump-target prompt, shared by
% Prompt_Start_Position (session start) and Process_Review_Segment's
% in-round 'j' command. HeaderText is the caller-specific framing line
% printed once before the shared Print_Jump_Options list, so each call
% site can phrase "why you're seeing this" differently while showing the
% exact same set of valid codes.
%
% Loops on an invalid entry (with an explanatory message) until either a
% valid target is entered (returns that position) or the reviewer cancels
% with a blank entry (returns []).

    NewPosition = [];
    while true
        fprintf('\n%s\n', HeaderText);
        Print_Jump_Options(NumTraces);
        RawInput = strtrim(input('> ','s'));

        if isempty(RawInput)
            return   % cancelled
        end

        [Candidate, Success, Message] = Resolve_Jump_Target(ReviewQueue, NumTraces, RawInput);
        if ~Success
            fprintf('   %s\n', Message);
            fprintf('   Try again (blank to cancel).\n');
            continue
        end

        JumpSegment = Find_Segment_For_Position(ReviewQueue, Candidate);
        fprintf('   Jumping to trace %d (%s / %s).\n', Candidate, JumpSegment.Tier, JumpSegment.SubgroupName);
        NewPosition = Candidate;
        return
    end
end

function Print_Jump_Options(NumTraces)
% Print_Jump_Options  The shared, multi-line list of valid jump-target
% codes -- one line per option rather than a single long sentence, so a
% reviewer can scan it at a glance. Shared by every caller of
% Prompt_Jump_Target so the offered codes can never drift out of sync
% between the session-start prompt and the in-round 'j' prompt.

    fprintf('     1-%d              a trace number\n', NumTraces);
    fprintf('     L / M / H        start of Low / Medium / High\n');
    fprintf('     H1 / H2A / H2B   start of that High subgroup\n');
    fprintf('     H3 / H4 / H5\n');
end

function [NewPosition, Success, Message] = Resolve_Jump_Target(ReviewQueue, NumTraces, RawInput)
% Resolve_Jump_Target  Parse a jump-target string into a queue position.
% Pure lookup, no I/O -- shared by every Prompt_Jump_Target call site.
% Case-insensitive. Accepts:
%   L / M / H                      -> start of that tier's first segment
%   H1 / H2A / H2B / H3 / H4 / H5   -> start of that specific High subgroup
%   a whole number 1-NumTraces      -> that exact queue position
% Success=false (with an explanatory Message) for anything else, including
% a section with zero traces in this dataset or an out-of-range number.

    NewPosition = [];
    Success = false;
    Message = '';

    Upper = upper(RawInput);

    TierLetters = struct('L','Low', 'M','Medium', 'H','High');
    if isfield(TierLetters, Upper)
        TargetTier = TierLetters.(Upper);
        TierSegments = ReviewQueue.Segments(strcmp({ReviewQueue.Segments.Tier}, TargetTier));
        if isempty(TierSegments)
            Message = sprintf('No %s-tier traces in this dataset.', TargetTier);
            return
        end
        NewPosition = TierSegments(1).StartIndex;
        Success = true;
        return
    end

    % SubgroupCodes/SubgroupNames are parallel lists -- ReviewQueue.Segments'
    % own SubgroupName casing is 'H2a'/'H2b', not the all-caps code a
    % reviewer types.
    SubgroupCodes = {'H1','H2A','H2B','H3','H4','H5'};
    SubgroupNames = {'H1','H2a','H2b','H3','H4','H5'};
    SubgroupIdx = find(strcmp(Upper, SubgroupCodes), 1);
    if ~isempty(SubgroupIdx)
        TargetSubgroup = SubgroupNames{SubgroupIdx};
        SubgroupSegments = ReviewQueue.Segments(strcmp({ReviewQueue.Segments.SubgroupName}, TargetSubgroup));
        if isempty(SubgroupSegments)
            Message = sprintf('No %s traces in this dataset.', TargetSubgroup);
            return
        end
        NewPosition = SubgroupSegments(1).StartIndex;
        Success = true;
        return
    end

    Requested = str2double(RawInput);
    if isnan(Requested) || Requested ~= floor(Requested)
        Message = 'Not a whole number, tier letter, or High subgroup code.';
        return
    end
    if Requested < 1 || Requested > NumTraces
        Message = sprintf('Out of range -- valid trace numbers are 1 to %d.', NumTraces);
        return
    end

    NewPosition = Requested;
    Success = true;
end

function [ManualFocusSubtractIndices, AnyRegistered] = Prompt_Manual_Focus_Correction(ManualFocusSubtractIndices, ...
    GlobalFindingImageFrame, GlobalFocusFrameNumbers)
% Prompt_Manual_Focus_Correction  Interactive loop for the round-prompt
% 's' command -- registers one or more reviewer-identified focus-jump
% frames (entered in clipped coordinates, same as the master grid's
% x-axis) into the session-wide, growing ManualFocusSubtractIndices list.
% Validated ONCE per entry, dataset-wide (Resolve_Manual_Focus_Index) --
% not per-trace, since FrameNumFound is numerically identical for every
% trace in one dataset (confirmed equal to
% OtherDataToSave.GlobalFindingImageFrame), so
% Global = EnteredIdx + GlobalFindingImageFrame is a single, dataset-wide
% conversion. A trace only actually gets the subtraction applied where the
% index is in-bounds for that trace's own clipped length -- see
% Apply_Manual_Focus_Corrections.m -- so no further per-trace validity
% check is needed here.
%
% Loops back for another entry after each successful (or duplicate)
% registration; exits only on 'N' (blank is not a shortcut -- it falls
% through to the "not a whole number" reprompt like any other bad entry).
%
% AnyRegistered is true iff at least one NEW index was actually added this
% call -- lets the caller (Process_Review_Segment) skip the immediate
% save/replot when the reviewer types 's' and then immediately backs out.

    AnyRegistered = false;
    while true
        fprintf('\nEnter a focus-frame index as you see it on the plots (n to stop):\n');
        RawInput = strtrim(input('> ','s'));
        if strcmpi(RawInput,'n')
            return
        end

        EnteredIdx = str2double(RawInput);
        if isnan(EnteredIdx) || EnteredIdx ~= floor(EnteredIdx)
            fprintf('Not a whole number -- try again.\n');
            continue
        end

        [Valid, Message] = Resolve_Manual_Focus_Index(EnteredIdx, GlobalFindingImageFrame, GlobalFocusFrameNumbers);
        if ~Valid
            fprintf('%s\n', Message);
            continue   % re-prompt for the SAME entry
        end

        if ismember(EnteredIdx, ManualFocusSubtractIndices)
            fprintf('Frame %d is already registered -- skipping duplicate.\n', EnteredIdx);
        else
            ManualFocusSubtractIndices = sort([ManualFocusSubtractIndices, EnteredIdx]);
            AnyRegistered = true;
            fprintf('Frame %d Registered. The jump due to focus at frame %d will now be subtracted out on every trace.\n', ...
                EnteredIdx, EnteredIdx);
        end
    end
end

function [Valid, Message] = Resolve_Manual_Focus_Index(EnteredIdx, GlobalFindingImageFrame, GlobalFocusFrameNumbers)
% Resolve_Manual_Focus_Index  Pure validation for one manual focus-jump
% index, entered in the SAME clipped coordinates the master grid's x-axis
% uses. Validated ONCE, dataset-wide (not per-trace) -- see
% Prompt_Manual_Focus_Correction's header for why.

    Global = EnteredIdx + GlobalFindingImageFrame;
    if ismember(Global, GlobalFocusFrameNumbers)
        Valid = true;
        Message = '';
    else
        Valid = false;
        Message = sprintf('Frame %d was not a focus frame.', EnteredIdx);
    end
end

function Segment = Find_Segment_For_Position(ReviewQueue, Position)
% Find_Segment_For_Position  Which ReviewQueue.Segments row Position
% currently falls in -- [] if Position is past every segment.

    Segment = [];
    for k = 1:numel(ReviewQueue.Segments)
        s = ReviewQueue.Segments(k);
        if Position >= s.StartIndex && Position <= s.EndIndex
            Segment = s;
            return
        end
    end
end

function TitleText = Get_Segment_Title_Text(Segment)
% Get_Segment_Title_Text  One-line caption text describing what's being
% reviewed in this segment -- shown in the bottom annotation strip below
% the whole subplot grid so a reviewer always knows which tier/subgroup/
% rule they're looking at.

    if strcmp(Segment.Tier, 'Low')
        TitleText = 'Low Priority';
        return
    end
    if strcmp(Segment.Tier, 'Medium')
        TitleText = sprintf('Medium Priority: %s', Segment.SubgroupName);
        return
    end

    switch Segment.SubgroupName
        case 'H1'
            TitleText = 'H1: Traces where the fusion/unbound event fell within 20 frames of the trace start or end';
        case 'H2a'
            TitleText = 'H2a: Traces where a focus event was the only reason a fusion call was vetoed (near-miss No Fusion)';
        case 'H2b'
            TitleText = 'H2b: Traces where a focus event occurred within 2 frames of the fusion frame';
        case 'H3'
            TitleText = 'H3: Traces where the winning fusion cluster was unusually small';
        case 'H4'
            TitleText = 'H4: Traces that narrowly failed exactly one scoring test';
        case 'H5'
            TitleText = 'H5: Traces designated Other (ambiguous)';
        otherwise
            TitleText = sprintf('%s / %s', Segment.Tier, Segment.SubgroupName);
    end
end

function [ReviewQueue, Changed] = Refresh_Segment_Grids(ReviewQueue, Options)
% Refresh_Segment_Grids  Re-derive every segment's Rows/Cols/TracesPerBatch/
% NumBatches from the CURRENT Options via Resolve_Grid_Config.m -- the same
% lookup Build_Review_Queue.m used at first-run preprocessing time, so a
% resumed queue's grid can never permanently drift from what live Options
% now say. Only grid/layout fields are touched; StartIndex/EndIndex/
% NumTraces/Tier/SubgroupName (trace bucketing/order) and the caller's
% CurrentQueuePosition are never read or written here.

    Changed = false;
    for k = 1:numel(ReviewQueue.Segments)
        Seg = ReviewQueue.Segments(k);
        GridConfig = Resolve_Grid_Config(Seg.Tier, Seg.SubgroupName, Options);

        if isfield(GridConfig, 'TracesPerBatch')
            TracesPerBatch = GridConfig.TracesPerBatch;
        else
            TracesPerBatch = GridConfig.Rows * GridConfig.Cols;
        end
        NumBatches = ceil(Seg.NumTraces / max(TracesPerBatch, 1));

        if Seg.Rows ~= GridConfig.Rows || Seg.Cols ~= GridConfig.Cols || ...
                ~isequal(getfield_or_empty(Seg, 'TracesPerBatch'), TracesPerBatch) || Seg.NumBatches ~= NumBatches
            Changed = true;
        end

        ReviewQueue.Segments(k).Rows           = GridConfig.Rows;
        ReviewQueue.Segments(k).Cols           = GridConfig.Cols;
        ReviewQueue.Segments(k).TracesPerBatch = TracesPerBatch;
        ReviewQueue.Segments(k).NumBatches     = NumBatches;
    end
end

function Value = getfield_or_empty(S, FieldName)
% getfield_or_empty  S.(FieldName) if it exists, else [] -- lets
% Refresh_Segment_Grids compare against a legacy segment struct that
% predates the TracesPerBatch field entirely without an isfield branch at
% every comparison site.

    if isfield(S, FieldName)
        Value = S.(FieldName);
    else
        Value = [];
    end
end

function HighOverrides = Get_High_Overrides(Options)
% Get_High_Overrides  Options.HighOverrides if present, else struct() --
% keeps the PreprocessingOptionsUsed isequal comparison well-defined
% regardless of whether the live/loaded Options ever set HighOverrides.

    if isfield(Options, 'HighOverrides')
        HighOverrides = Options.HighOverrides;
    else
        HighOverrides = struct();
    end
end

function BatchSize = Get_Segment_Batch_Size(Segment)
% Get_Segment_Batch_Size  Traces consumed per round for this segment.
% TracesPerBatch defaults to Rows*Cols for every segment except H1 (2
% physical axes per trace -- one in each of its fixed 2 rows) -- see
% Build_Review_Queue.m.

    if isfield(Segment, 'TracesPerBatch')
        BatchSize = Segment.TracesPerBatch;
    else
        BatchSize = Segment.Rows * Segment.Cols;
    end
end

function PrevSeg = Find_Previous_Nonempty_Segment(ActiveSegments, SegmentIdx)
% Find_Previous_Nonempty_Segment  The nearest segment before ActiveSegments
% (SegmentIdx) that actually has traces in it -- [] if SegmentIdx is 1 or
% every earlier segment is empty. Skips empty segments (e.g. a dataset
% with zero '2 Fuse' traces) rather than landing on a segment with no
% valid last batch to go back to.

    PrevSeg = [];
    for k = SegmentIdx-1 : -1 : 1
        if ActiveSegments(k).NumTraces > 0
            PrevSeg = ActiveSegments(k);
            return
        end
    end
end

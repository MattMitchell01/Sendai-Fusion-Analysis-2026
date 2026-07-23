function [FigureHandles] = Plot_Current_Trace(FigureHandles, CurrentVirusData, VideoTimeVector, ClipWidth, ...
    PlotCounter, CurrentTraceNumber, Options, DiagnosticOverlay, ManualFocusSubtractIndices, ClearedIndexRanges)
%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
% Plot_Current_Trace  Draw one trace into the master grid window.
%
% Tier-agnostic: which grid cell it lands in comes from
% PlotCounter/FigureHandles.SubHandles (set by Create_Master_Window from
% that segment's Rows/Cols) -- nothing here branches on Tier, so Low,
% Medium, and High all share this function.
%
% VideoTimeVector is accepted for call-site-shape consistency with other
% callers of this parameter (see PART_B_HANDOFF_NOTES.md section 1a) but
% is not read here -- the only thing this function used to need from
% OtherDataToSave.UniversalData (StandardBindFrameNum) always equals
% FrameNumFound, so it's hardcoded below instead.
%
% ClipWidth (Part B's OtherDataToSave.ClipWidth, or [] if a dataset
% predates that field) is used two ways: (1) when Options.UseRunMed=='y',
% the visible x-range excludes the first/last ClipWidth frames -- Part B
% can never report a fusion/unbound event there (see "dead zone" crop
% below); (2) passed through to Get_H2_Focus_Jump_Trigger_Frame.m for the
% H2a/H2b focus-jump-removal overlay.
%
% DiagnosticOverlay (optional, default []) is a forward-compat extension
% point for High-tier overlays, carrying whichever optional fields a given
% subgroup needs:
%   .ShowFocusDots -- draws every focus event as a red/black dot pair
%     via Plot_Focus_Dot_Markers.m (the same marker
%     Fix_Fusion_Wait_Time.m/Fix_Unbind_Wait_Time.m already draw in the
%     blown-up picker figure). Used by H2a, H2b, H3, and H5.
%   .TitleSuffix -- an extra line appended to the subplot title. Used by
%     H4 to name which scoring test the trace narrowly failed (see
%     Get_H4_Failed_Test_Label.m).
%   .ApplyFocusJumpCorrection -- H2a/H2b only: subtracts the focus-induced
%     step (value_at_focus_frame - value_at_frame_before) from every
%     index from the triggering focus frame onward, so the reviewer sees
%     the trace as if that one focus jump hadn't happened. Display-only;
%     never touches saved data. See Get_H2_Focus_Jump_Trigger_Frame.m.
% H1 doesn't use this argument at all (its overlay is a pane-layout
% difference, see Plot_Trace_Pair_H1.m).
%
% ManualFocusSubtractIndices (optional, default []) -- session-wide,
% growing list of reviewer-registered clipped-coordinate focus-jump
% frames (Start_User_Review.m's 's' round command), applied via
% Apply_Manual_Focus_Corrections.m to EVERY trace where in-bounds, in
% every tier/segment -- unlike DiagnosticOverlay.ApplyFocusJumpCorrection
% above (H2a/H2b only, auto-detected), this is universal and manually
% curated, and can compose with that H2 correction on the same trace.
%
% ClearedIndexRanges (optional, default []) -- session-wide, growing Nx2
% list of reviewer-registered [Start,End] clipped-coordinate index ranges
% (Start_User_Review.m's 'c' round command), NaN'd out via
% Apply_Cleared_Index_Ranges.m so no dots/line are drawn there, on EVERY
% trace where in-bounds. Purely a display exclusion -- never touches saved
% data, unlike the focus-jump corrections above. Applied AFTER Run_Med
% smoothing (see below), not before, so it can't be interpolated across.
%
% Plots the CLIPPED trace (Clip_Trace_For_Review) since that's the
% coordinate system FuseFrameNumbers/etc. are meaningful in. Part B
% stores FusionData.FuseFrameNumbers in GLOBAL coordinates (see
% Analyze_Trace_Data.m in Part B), so any frame drawn against the clipped
% trace here must be converted via clippedIdx = globalIdx - FrameNumFound
% -- clip(1) corresponds to global frame FrameNumFound+1.

if nargin < 8
    DiagnosticOverlay = [];
end
if nargin < 9
    ManualFocusSubtractIndices = [];
end
if nargin < 10
    ClearedIndexRanges = [];
end

FusionData    = CurrentVirusData.FusionData;
UnboundData   = CurrentVirusData.UnboundData;
ChangedByUser = CurrentVirusData.ChangedByUser;
FrameNumFound = CurrentVirusData.FrameNumFound;

[ClippedTrace, ClippedFocusFrames] = Clip_Trace_For_Review(CurrentVirusData);

% Focus-jump correction (H2a/H2b only) runs BEFORE Run_Med smoothing --
% Part B's own Seek_Fusion.m computes its Test 5 focus-jump ratio against
% the RAW trace specifically because "focus events are 1-2 frame spikes
% that the median smoother suppresses"; computing the jump size after
% smoothing here would underestimate it and leave a smeared residual.
if ~isempty(DiagnosticOverlay) && isfield(DiagnosticOverlay,'ApplyFocusJumpCorrection') && DiagnosticOverlay.ApplyFocusJumpCorrection
    TriggerFrame = Get_H2_Focus_Jump_Trigger_Frame(CurrentVirusData, ClippedFocusFrames, ClipWidth);
    if ~isempty(TriggerFrame) && TriggerFrame >= 2 && TriggerFrame <= numel(ClippedTrace)
        FocusJump = ClippedTrace(TriggerFrame) - ClippedTrace(TriggerFrame - 1);
        ClippedTrace(TriggerFrame:end) = ClippedTrace(TriggerFrame:end) - FocusJump;
        DiagnosticOverlay.TitleSuffix = 'Focus-jump corrected';
    end
end

% Manual, reviewer-curated focus-jump correction -- universal (every
% tier/segment), independent of and additive to the H2a/H2b-only automatic
% correction above (a trace can get both). Runs BEFORE Run_Med for the
% same reason: a 1-2 frame focus spike is what the median smoother would
% suppress/underestimate.
%
% Deliberately does NOT add a title suffix for this -- per the user, an
% "(Manual focus correction applied)" line was cluttering/distorting the
% subplot title. The correction is still applied to the plotted trace
% itself; only the title annotation was removed. (H2's 'Focus-jump
% corrected' and H4's failed-test label are unaffected -- those stay.)
if ~isempty(ManualFocusSubtractIndices)
    ClippedTrace = Apply_Manual_Focus_Corrections(ClippedTrace, ManualFocusSubtractIndices);
end

if strcmp(Options.UseRunMed, 'y')
    ClippedTrace = Run_Med(ClippedTrace, Options);
end

% Reviewer-curated "clear dots" exclusion ('c' round command) -- runs
% AFTER Run_Med, not before, unlike the focus-jump corrections above:
% Run_Med.m's median(...,'omitnan') would otherwise interpolate right
% across a NaN'd gap using real neighboring values and silently un-clear
% it. Applying this last, right before the trace is drawn, guarantees the
% excluded indices reach plot() as NaN -- MATLAB draws a broken line there
% natively, never an interpolated jump across the gap.
%
% Deliberately does NOT add a title suffix for this either, same reasoning
% as the manual focus correction above -- "(Dots cleared in registered
% range(s))" was cluttering/distorting the subplot title. The exclusion
% itself is still applied to the plotted trace.
if ~isempty(ClearedIndexRanges)
    ClippedTrace = Apply_Cleared_Index_Ranges(ClippedTrace, ClearedIndexRanges);
end

% Focus-event frames are NOT gapped out here -- every index is plotted,
% focused or not. (An earlier version NaN'd them out; reverted per user.)

% Visible/plotted range: when running-median smoothing is on, the first
% and last ClipWidth frames are a dead zone -- Part B can never report a
% fusion/unbound event there (see PART_B_HANDOFF_NOTES.md section 5), so
% there's no review value in showing them. This is computed BEFORE
% plotting and actually trims what gets drawn (not just xlim) so those
% frames are never part of the line itself -- otherwise Options.ExpandXAxis's
% +/-30 padding (much larger than ClipWidth) re-reveals them regardless of
% xlim, which is exactly the bug this replaced. Run_Med above already
% smoothed using the FULL, untrimmed trace as context, so the boundary
% values retained here were computed with genuine neighbor data reaching
% into the dead zone, not an artificially shrunk window (mirrors Part B's
% own Calculate_Sliding_Window.m, which smooths the full array first and
% only then drops ClipWidth off each end). Explicit x-coordinates
% (VisibleStart:VisibleEnd, not the default 1:length) keep this aligned
% with every marker/overlay below, which are all computed in true clipped-
% index space.
VisibleStart = 1;
VisibleEnd = numel(ClippedTrace);
if strcmp(Options.UseRunMed, 'y') && ~isempty(ClipWidth) && ClipWidth > 0 ...
        && (numel(ClippedTrace) - 2*ClipWidth) >= 1
    VisibleStart = ClipWidth + 1;
    VisibleEnd = numel(ClippedTrace) - ClipWidth;
end
VisibleX = VisibleStart:VisibleEnd;
VisibleY = ClippedTrace(VisibleStart:VisibleEnd);

% Plot the current trace to the current subplot axis
set(0,'CurrentFigure',FigureHandles.MasterWindow)
NumSubHandles = length(FigureHandles.SubHandles);
PlotCounterToUse = min(PlotCounter, NumSubHandles);
set(FigureHandles.MasterWindow,'CurrentAxes',FigureHandles.SubHandles(PlotCounterToUse));
cla
hold on

if strcmp(FusionData.Designation,'2 Fuse')

    plot(VisibleX,VisibleY,'r-')
    LineToPlot = ylim;
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');

    if strcmp(ChangedByUser,'Incorrect Designation-Changed') || strcmp(ChangedByUser,'Correct Designation, Incorrect Wait Time')
        % Just re-designated, fuse frame(s) not picked yet -- nothing to plot.
        % (Same guard as the 1 Fuse branch below -- previously missing here,
        % so a stale FuseFrameNumbers left over from a prior designation
        % could get drawn against the freshly-corrected 2 Fuse trace.)
    elseif ~isempty(FusionData.FuseFrameNumbers)
        ClippedFuse1 = FusionData.FuseFrameNumbers(1) - FrameNumFound;
        plot([ClippedFuse1, ClippedFuse1], LineToPlot, 'g--')
        if numel(FusionData.FuseFrameNumbers) > 1
            % Dark green, not magenta/red -- the trace itself is already
            % plotted red for 2 Fuse, so the second fuse-frame marker uses
            % a distinct shade of green instead of competing with it.
            ClippedFuse2 = FusionData.FuseFrameNumbers(2) - FrameNumFound;
            plot([ClippedFuse2, ClippedFuse2], LineToPlot, '--', 'Color', [0 0.5 0])
        end
    end

elseif strcmp(FusionData.Designation,'1 Fuse')

    plot(VisibleX,VisibleY,'b-')
    LineToPlot = ylim;
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');

    if strcmp(ChangedByUser,'Incorrect Designation-Changed') || strcmp(ChangedByUser,'Correct Designation, Incorrect Wait Time')
        % Just re-designated, fuse frame not picked yet -- nothing to plot.
    elseif ~isempty(FusionData.FuseFrameNumbers)
        ClippedFuse1 = FusionData.FuseFrameNumbers(1) - FrameNumFound;
        plot([ClippedFuse1, ClippedFuse1], LineToPlot, 'g--')
    end

elseif strcmp(FusionData.Designation,'No Fusion')

    plot(VisibleX,VisibleY,'k-')
    LineToPlot = ylim;
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');

elseif strcmp(FusionData.Designation,'Unbound')

    plot(VisibleX,VisibleY,'c-')
    LineToPlot = ylim;
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');

    if strcmp(ChangedByUser,'Incorrect Designation-Changed') || strcmp(ChangedByUser,'Correct Designation, Incorrect Wait Time')
        % Just re-designated, unbind frame not picked yet -- nothing to plot.
        % (Previously missing the repick-mode ('.33') check here, unlike
        % the 1 Fuse branch's equivalent guard.)
    elseif isfield(UnboundData,'UnboundFrameNumbers') && ~isempty(UnboundData.UnboundFrameNumbers)
        ClippedUnbind = UnboundData.UnboundFrameNumbers(1) - FrameNumFound;
        plot([ClippedUnbind, ClippedUnbind], LineToPlot, 'g--')
    end

elseif strcmp(FusionData.Designation,'Other')

    plot(VisibleX,VisibleY,'m-')
    LineToPlot = ylim;
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');

else

    plot(VisibleX,VisibleY,'k.')
    LineToPlot = ylim;
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
end

if ~isempty(DiagnosticOverlay) && isfield(DiagnosticOverlay, 'ShowFocusDots') && DiagnosticOverlay.ShowFocusDots
    Plot_Focus_Dot_Markers(ClippedTrace, ClippedFocusFrames);
end

% Kept deliberately terse: plot number, designation, virus ID. No global
% trace number (redundant with the plot number within a round) and no
% X/Y coordinates -- clutters the subplot grid.
Title = sprintf('%d - %s (ID %d)', PlotCounter, FusionData.Designation, ...
    CurrentVirusData.VirusIDNumber);
if ~isempty(DiagnosticOverlay) && isfield(DiagnosticOverlay, 'TitleSuffix')
    % H4 uses this to name which scoring test the trace narrowly failed
    % (see Get_H4_Failed_Test_Label.m) -- a second title line, not baked
    % into the format string above, so this stays a generic extension
    % point rather than something only H4 can use.
    Title = sprintf('%s\n%s', Title, DiagnosticOverlay.TitleSuffix);
end
title(Title);

if strcmp(Options.ShowBindFrame,'y')
    % StandardBindFrameNum == FrameNumFound always (Part B convention --
    % see PART_B_HANDOFF_NOTES.md section 1a, UniversalData no longer
    % exists), so the bind frame is always clipped-coordinate 0.
    ClippedBindFrame = 0;
    plot([ClippedBindFrame, ClippedBindFrame], LineToPlot, 'k--')
end

% VisibleStart/VisibleEnd (computed above, before plotting) already
% trimmed the drawn line itself, so this padding is now purely cosmetic
% breathing room around real plotted data -- it can no longer re-reveal
% the dead zone the way it used to when only xlim (not the line) was
% cropped.
if strcmp(Options.ExpandXAxis,'y')
    xlim([VisibleStart - 30, VisibleEnd + 30])
else
    xlim([VisibleStart, VisibleEnd])
end

hold off
drawnow
end

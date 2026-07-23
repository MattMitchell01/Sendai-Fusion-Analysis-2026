function [ClippedTrace, ClippedFrameNums] = Plot_Trace_With_Focus_Markers(TraceStruct, TitleStr, ClipWidth, Options, ManualFocusSubtractIndices, ClearedIndexRanges)
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
% Plot_Trace_With_Focus_Markers  Draws the clipped trace + bind-frame line
% + focus-frame dot-pair overlay into whatever figure/axes is currently
% current (caller sets that via set(0,'CurrentFigure',...) first, same
% convention as the callers below).
%
% Shared by Fix_Fusion_Wait_Time.m and Fix_Unbind_Wait_Time.m so this exact
% look ("the raw trace + focus frames, like when picking a fusion/unbind
% frame") only needs to be maintained in one place.
%
% Returns the clipped trace and its frame-number vector so callers that
% do their own ginput/picking don't have to re-derive them. ClippedFrameNums
% is already restricted to the same dead-zone-excluded range this figure
% plots (see below), so a caller doing nearest-frame ginput matching against
% it can never resolve a click to an unscoreable frame.
%
% ManualFocusSubtractIndices -- session-wide, reviewer-registered list from
% Start_User_Review.m's 's' round command (see Plot_Current_Trace.m for the
% full explanation) -- applied here too so this figure never shows a raw
% trace the master grid doesn't.
%
% ClearedIndexRanges -- session-wide, reviewer-registered list of
% [Start,End] ranges from Start_User_Review.m's 'c' round command (see
% Plot_Current_Trace.m for the full explanation) -- applied here too so
% this figure never shows a gap the master grid doesn't. Applied AFTER
% Run_Med (see below), same reasoning as Plot_Current_Trace.m: Run_Med's
% median(...,'omitnan') would otherwise interpolate right across a NaN'd
% gap using real neighboring values and silently un-clear it.

    [ClippedTrace, ClippedFocusFrames] = Clip_Trace_For_Review(TraceStruct);

    % Same automatic H2a/H2b focus-jump-removal correction the master grid
    % applies (Plot_Current_Trace.m) -- gated on the trace's OWN frozen
    % Review_PriorityData.Rule (via Get_H2_Focus_Jump_Trigger_Frame.m,
    % which already returns [] for any non-H2a/H2b trace), not on which
    % Segment/Tier happened to call this function -- this picker has no
    % Segment in scope at all, so keying off the trace's own data is both
    % simpler and guarantees this can never disagree with the master grid
    % about whether a given trace is jump-corrected. Runs BEFORE the manual
    % correction and Run_Med below, same order Plot_Current_Trace.m uses
    % and for the same reason: a 1-2 frame focus spike is exactly what the
    % median smoother would suppress/underestimate if applied first.
    FocusJumpCorrected = false;
    if strcmp(Options.ApplyH2FocusJumpCorrection, 'y') && ~isempty(ClipWidth)
        TriggerFrame = Get_H2_Focus_Jump_Trigger_Frame(TraceStruct, ClippedFocusFrames, ClipWidth);
        if ~isempty(TriggerFrame) && TriggerFrame >= 2 && TriggerFrame <= numel(ClippedTrace)
            FocusJump = ClippedTrace(TriggerFrame) - ClippedTrace(TriggerFrame - 1);
            ClippedTrace(TriggerFrame:end) = ClippedTrace(TriggerFrame:end) - FocusJump;
            FocusJumpCorrected = true;
        end
    end

    % Same manual, reviewer-curated correction the master grid applies
    % (Plot_Current_Trace.m) -- for consistency, this picker figure should
    % never show a raw trace the grid itself doesn't. Runs BEFORE Run_Med,
    % same reasoning as Plot_Current_Trace.m: a 1-2 frame focus spike is
    % what the median smoother would suppress/underestimate if applied first.
    % No title annotation for this -- per the user, "(Manual focus
    % correction applied)" was cluttering/distorting the title. The
    % correction itself is still applied to the plotted trace.
    if ~isempty(ManualFocusSubtractIndices)
        ClippedTrace = Apply_Manual_Focus_Corrections(ClippedTrace, ManualFocusSubtractIndices);
    end

    % Same running-median smoothing as the master grid (Plot_Current_Trace.m),
    % gated on the same Options.UseRunMed/Options.RunMedHalfLength -- one
    % setting controls smoothing everywhere in Part C, so a trace never
    % looks different here than it does in the grid it was flagged from.
    % Smoothing only changes VALUES, never the trace's length/index
    % alignment, so the click-matching against ClippedFrameNums below (and
    % the frame numbers stored back into FusionData/UnboundData) are
    % unaffected either way.
    if strcmp(Options.UseRunMed, 'y')
        ClippedTrace = Run_Med(ClippedTrace, Options);
    end

    % Reviewer-curated "clear dots" exclusion ('c' round command) -- same
    % NaN-gap treatment as Plot_Current_Trace.m, applied here too so this
    % picker figure is never gap-free where the master grid isn't. No title
    % annotation for this either, same reasoning as the manual correction
    % above -- the exclusion itself is still applied to the plotted trace.
    if ~isempty(ClearedIndexRanges)
        ClippedTrace = Apply_Cleared_Index_Ranges(ClippedTrace, ClearedIndexRanges);
    end

    % Same dead-zone exclusion as the master grid (Plot_Current_Trace.m):
    % Part B can never report a fusion/unbound event in the first/last
    % ClipWidth frames (PART_B_HANDOFF_NOTES.md section 5), so a reviewer
    % should never be able to pick one here either -- otherwise this figure
    % could offer a click target the master grid never showed at all.
    % Gated on Options.UseRunMed the same way the master grid's crop is, so
    % the two stay in lockstep regardless of that setting.
    VisibleStart = 1;
    VisibleEnd = numel(ClippedTrace);
    if strcmp(Options.UseRunMed, 'y') && ~isempty(ClipWidth) && ClipWidth > 0 ...
            && (numel(ClippedTrace) - 2*ClipWidth) >= 1
        VisibleStart = ClipWidth + 1;
        VisibleEnd = numel(ClippedTrace) - ClipWidth;
    end

    % StandardBindFrameNum == FrameNumFound always (Part B convention --
    % see PART_B_HANDOFF_NOTES.md section 1a, UniversalData no longer
    % exists), so the bind frame is always clipped-coordinate 0.
    ClippedBindFrame = 0;
    ClippedFrameNums = VisibleStart:VisibleEnd;

    cla
    hold on

    % Same designation-color palette as Plot_Current_Trace.m -- by the
    % time this figure opens, TraceStruct.FusionData.Designation already
    % holds the JUST-CORRECTED value (Correct_Designations.m sets it
    % before calling Fix_Fusion_Wait_Time/Fix_Unbind_Wait_Time), so a trace flagged
    % from No Fusion to 1 Fuse shows up blue here immediately, not black.
    switch TraceStruct.FusionData.Designation
        case '1 Fuse'
            TraceColor = 'b-';
        case '2 Fuse'
            TraceColor = 'r-';
        case 'Unbound'
            TraceColor = 'c-';
        case 'Other'
            TraceColor = 'm-';
        otherwise
            TraceColor = 'k-';   % 'No Fusion' and any unrecognized designation
    end

    plot(ClippedFrameNums, ClippedTrace(VisibleStart:VisibleEnd), TraceColor)
    if FocusJumpCorrected
        title(sprintf('%s\nFocus-jump corrected', TitleStr))
    else
        title(TitleStr)
    end

    LineToPlot = ylim;
    plot([ClippedBindFrame, ClippedBindFrame], LineToPlot, 'k--')

    % Focus-frame overlay: red dot before / black dot on each focus event,
    % on the actual (un-gapped) trace values, so the reviewer can see the
    % jump when picking a frame. Shared with Plot_Current_Trace.m's
    % High-tier H2a/H2b overlay -- see Plot_Focus_Dot_Markers.m.
    Plot_Focus_Dot_Markers(ClippedTrace, ClippedFocusFrames);
end

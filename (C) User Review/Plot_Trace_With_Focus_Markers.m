function [ClippedTrace, ClippedFrameNums] = Plot_Trace_With_Focus_Markers(TraceStruct, TitleStr, ClipWidth, Options, ManualFocusSubtractIndices)
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
% + focus-frame dot-triple overlay into whatever figure/axes is currently
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

    % Deliberately does NOT apply Run_Med here, even when Options.UseRunMed
    % is 'y' -- this figure is for picking an exact frame (fuse/unbind), so
    % it always shows the RAW trace, unlike the master grid
    % (Plot_Current_Trace.m), which does smooth when that option is on.
    [ClippedTrace, ClippedFocusFrames] = Clip_Trace_For_Review(TraceStruct);

    % Same manual, reviewer-curated correction the master grid applies
    % (Plot_Current_Trace.m) -- for consistency, this picker figure should
    % never show a raw trace the grid itself doesn't. Order relative to
    % the dead-zone crop below doesn't matter -- this only changes values,
    % never the trace's length.
    ManualCorrectionApplied = false;
    if ~isempty(ManualFocusSubtractIndices)
        [ClippedTrace, ManualCorrectionApplied] = Apply_Manual_Focus_Corrections(ClippedTrace, ManualFocusSubtractIndices);
    end

    % Same dead-zone exclusion as the master grid (Plot_Current_Trace.m):
    % Part B can never report a fusion/unbound event in the first/last
    % ClipWidth frames (PART_B_HANDOFF_NOTES.md section 5), so a reviewer
    % should never be able to pick one here either -- otherwise this figure
    % could offer a click target the master grid never showed at all.
    % Gated on Options.UseRunMed the same way the master grid's crop is, so
    % the two stay in lockstep regardless of that setting; only the
    % smoothing itself (never applied here) differs between the two, by
    % design, for precise click-picking.
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
    if ManualCorrectionApplied
        title(sprintf('%s\n(Manual focus correction applied)', TitleStr))
    else
        title(TitleStr)
    end

    LineToPlot = ylim;
    plot([ClippedBindFrame, ClippedBindFrame], LineToPlot, 'k--')

    % Focus-frame overlay: red dot before / black dot on / blue dot after
    % each focus event, on the actual (un-gapped) trace values, so the
    % reviewer can see the jump when picking a frame. Shared with
    % Plot_Current_Trace.m's High-tier H2a/H2b overlay -- see
    % Plot_Focus_Dot_Markers.m.
    Plot_Focus_Dot_Markers(ClippedTrace, ClippedFocusFrames);
end

function [ClippedTrace, Applied] = Apply_Cleared_Index_Ranges(ClippedTrace, ClearedIndexRanges)
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
% Apply_Cleared_Index_Ranges  NaNs out every reviewer-registered "clear
% dots" index range (Start_User_Review.m's 'c' round command) in
% ClippedTrace, clamped to what's actually in-bounds for THIS trace's
% (shorter or longer) clipped length -- shared by Plot_Current_Trace.m
% (master grid) and Plot_Trace_With_Focus_Markers.m (picker figure) so the
% exact same gap is visible in both places.
%
% ClearedIndexRanges is an Nx2 array of [Start, End] inclusive-range rows
% (registered dataset-wide, session-wide -- every trace shares the same
% clipped-coordinate meaning for a given index, same as
% ManualFocusSubtractIndices). A row that falls entirely outside this
% trace's length is silently a no-op for that trace.
%
% Purely a display transform -- NaN in a MATLAB plot() call draws a broken
% line (gap), never an interpolated jump across the gap, so no extra
% plotting logic is needed beyond setting these values to NaN before the
% trace is drawn. Must be applied AFTER any running-median smoothing
% (Run_Med.m uses median(...,'omitnan'), which would otherwise interpolate
% right across a NaN gap using real neighboring values and silently
% un-clear it) -- see Plot_Current_Trace.m for where this is called.
%
% Applied is true iff at least one row in ClearedIndexRanges actually
% overlapped this trace's in-bounds length -- lets callers add a
% title-suffix line only on traces the clearing actually touched.

    Applied = false;
    for k = 1:size(ClearedIndexRanges, 1)
        Start = max(1, ClearedIndexRanges(k,1));
        End   = min(numel(ClippedTrace), ClearedIndexRanges(k,2));
        if Start <= End
            ClippedTrace(Start:End) = NaN;
            Applied = true;
        end
    end
end

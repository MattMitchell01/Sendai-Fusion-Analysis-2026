function [ClippedTrace, Applied] = Apply_Manual_Focus_Corrections(ClippedTrace, ActiveIndices)
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
% Apply_Manual_Focus_Corrections  Subtracts a reviewer-curated focus-jump
% artifact from ClippedTrace at every index in ActiveIndices that's
% actually in-bounds for THIS trace's (shorter or longer) clipped length --
% shared by Plot_Current_Trace.m (master grid) and
% Plot_Trace_With_Focus_Markers.m (picker figure) so the exact same
% correction is visible in both places.
%
% ActiveIndices was already validated dataset-wide, once, at registration
% time (Start_User_Review.m's 's' round command, Resolve_Manual_Focus_Index)
% -- every trace in one dataset shares the same GlobalFindingImageFrame, so
% a clipped-coordinate index means the same real moment for every trace.
% This function only needs a per-trace BOUNDS check, not a re-validation
% against GlobalFocusFrameNumbers.
%
% Identical math to the existing H2a/H2b auto-correction
% (Plot_Current_Trace.m / Get_H2_Focus_Jump_Trigger_Frame.m): the value at
% the jump minus the value just before it, subtracted from that index
% onward. When more than one index is in-bounds for a given trace, each is
% applied in ascending order, sequentially -- each subsequent jump is
% computed against the ALREADY-partially-corrected trace.
%
% Applied is true iff at least one index in ActiveIndices was actually
% in-bounds (>=2 and <= numel(ClippedTrace)) for this trace -- lets
% callers add a title-suffix line only on traces the correction actually
% touched.

    Applied = false;
    SortedIndices = sort(ActiveIndices);
    for k = 1:numel(SortedIndices)
        Idx = SortedIndices(k);
        if Idx >= 2 && Idx <= numel(ClippedTrace)
            Jump = ClippedTrace(Idx) - ClippedTrace(Idx - 1);
            ClippedTrace(Idx:end) = ClippedTrace(Idx:end) - Jump;
            Applied = true;
        end
    end
end

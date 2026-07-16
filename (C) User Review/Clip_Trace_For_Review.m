function [ClippedTrace, ClippedFocusFrames] = Clip_Trace_For_Review(TraceStruct)
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
% Clip_Trace_For_Review  Reproduce Part B's global -> clipped trace/frame
% conversion for display purposes.
%
% WARNING: THIS MUST MATCH Analyze_Trace_Data.m IN PART B EXACTLY. If that
% logic changes, change it here too.
%
% Revised per PART_B_HANDOFF_NOTES.md (Part A/B breaking change): Part A's
% StartAnalysisFrameNumber is now FrameNumToFindParticles + 1 (previously
% equal to it), so Trace_BackSub/TimeVector no longer contain the finding
% image at all -- Trace_BackSub's own native index 1 is ALREADY global
% frame FrameNumFound+1. No slicing is needed (or correct) to exclude the
% finding image; it was never included in the first place under this
% convention. (This supersedes an earlier revision of this file that did
% Trace_BackSub(2:end) under the OLD convention, where index 1 truly was
% the finding image -- that convention no longer matches what Part A/B
% actually produce.)
%
% Because ClippedTrace(1) == global frame FrameNumFound+1, global
% (raw-frame-number-equivalent) <-> clipped conversion is a plain
% clippedIdx = globalIdx - FrameNumFound / globalIdx = clippedIdx +
% FrameNumFound -- no extra +/-1. This formula itself is unchanged from
% before; only the clip boundary moved.
% ClippedFocusFrames is FocusFrameNumbers_Shifted converted into that same
% clipped coordinate system, filtered to the valid range.

    FrameNumFound = TraceStruct.FrameNumFound;
    ClippedTrace  = TraceStruct.Trace_BackSub;

    FocusGlobal = [];
    if isfield(TraceStruct, 'focusframenumbers') && ~isempty(TraceStruct.focusframenumbers)
        ff = TraceStruct.focusframenumbers;
        FocusGlobal = ff(~isnan(ff(:)))';
    elseif isfield(TraceStruct, 'FocusFrameNumbers_Shifted') && ~isempty(TraceStruct.FocusFrameNumbers_Shifted)
        ff = TraceStruct.FocusFrameNumbers_Shifted;
        % FocusFrameNumbers_Shifted now comes out of Part A already in
        % clipped coordinates (since StartAnalysisFrameNumber-1 ==
        % FrameNumFound under the new convention) -- reconstruct global
        % with +FrameNumFound, NOT +1 (PART_B_HANDOFF_NOTES.md section 2).
        FocusGlobal = ff(~isnan(ff(:)))' + FrameNumFound;
    end

    ClippedFocusFrames = FocusGlobal - FrameNumFound;
    ClippedFocusFrames = ClippedFocusFrames(ClippedFocusFrames >= 1 & ClippedFocusFrames <= numel(ClippedTrace));
end

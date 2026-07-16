function [TriggerFocusFrame_Clipped] = Get_H2_Focus_Jump_Trigger_Frame(CurrentVirusData, ClippedFocusFrames, ClipWidth)
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
% Get_H2_Focus_Jump_Trigger_Frame  Identify the single focus-event frame
% (Part C clipped coordinates) responsible for an H2a/H2b High-priority
% flag, for Plot_Current_Trace.m's focus-jump-removal overlay. Returns []
% if none can be reliably identified.
%
% Review_PriorityData.Score15Frames / .TriggerFusionFrame /
% .TriggerFocusFrame (set by Part B's Assign_Review_Priority.m) live in
% Seek_Fusion's internally-trimmed ScoreArray-index coordinate space, NOT
% Part C's own clippedIdx = globalIdx - FrameNumFound convention used by
% Clip_Trace_For_Review.m. Verified conversion (PART_B_HANDOFF_NOTES.md
% section 6, cross-checked against Analyze_Trace_Data.m's own
% fusionFramesGlobal = fusionFramesClipped + fs + ClipWidth line):
%   PartC_clippedIdx = value + ClipWidth
%
% H2b: Review_PriorityData.TriggerFocusFrame is Part B's own direct
% record of which focus frame tripped the Test 5 veto near the actual
% fusion frame. Converted and cross-checked against ClippedFocusFrames
% (must land within +/-1 of a real entry) before use -- a defensive
% canary against any coordinate-space mismatch, rather than trusting the
% converted value blindly.
%
% H2a: Part B does NOT record which focus frame tripped the veto for each
% Score15Frames candidate (frames where Tests 1-4 passed but Test 5
% failed) -- only the candidate frame(s) themselves. Score15Frames is
% clustered by contiguous run (gap of exactly 1 = same cluster, any
% larger gap = a separate, unrelated near-miss region). If more than one
% cluster exists, this function declines to guess which one matters and
% returns [] -- silently picking the wrong one would apply a real-looking
% but incorrect correction. With exactly one cluster, the nearest actual
% ClippedFocusFrames entry to any frame in that cluster is used --
% justified because H2a's own membership criterion guarantees a
% genuinely nearby focus event exists (Test 5's search window is small
% and bounded), without needing to duplicate Seek_Fusion's exact
% ratio-maximization window logic into Part C.

    TriggerFocusFrame_Clipped = [];

    if isempty(ClippedFocusFrames)
        return
    end
    if ~isfield(CurrentVirusData, 'Review_PriorityData')
        return
    end
    PriorityData = CurrentVirusData.Review_PriorityData;
    if ~isfield(PriorityData, 'Rule')
        return
    end

    if strcmp(PriorityData.Rule, 'H2b')
        if ~isfield(PriorityData, 'TriggerFocusFrame') || isempty(PriorityData.TriggerFocusFrame)
            return
        end
        Candidate = PriorityData.TriggerFocusFrame(1) + ClipWidth;
        Distances = abs(ClippedFocusFrames(:) - Candidate);
        [MinDist, MatchIdx] = min(Distances);
        if MinDist <= 1
            % Snap to the actual matched ClippedFocusFrames entry, not the
            % converted Candidate value -- they should coincide exactly
            % (the +ClipWidth conversion is exact, not approximate), but if
            % the ±1 tolerance ever catches a real discrepancy, the entry
            % Part C itself already computed as a genuine focus frame is
            % the more trustworthy of the two.
            TriggerFocusFrame_Clipped = ClippedFocusFrames(MatchIdx);
        end

    elseif strcmp(PriorityData.Rule, 'H2a')
        if ~isfield(PriorityData, 'Score15Frames') || isempty(PriorityData.Score15Frames)
            return
        end
        SortedCandidates = sort(PriorityData.Score15Frames(:)');
        Gaps = diff(SortedCandidates);
        NumClusters = 1 + sum(Gaps > 1);
        if NumClusters > 1
            return   % multiple disjoint near-miss regions -- decline rather than guess
        end

        CandidatesClipped = SortedCandidates + ClipWidth;
        Distances = abs(ClippedFocusFrames(:) - CandidatesClipped(:)');
        [~, LinearIdx] = min(Distances(:));
        [FocusRow, ~] = ind2sub(size(Distances), LinearIdx);
        TriggerFocusFrame_Clipped = ClippedFocusFrames(FocusRow);
    end
end

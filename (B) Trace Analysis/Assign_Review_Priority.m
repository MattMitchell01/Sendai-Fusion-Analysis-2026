function [Review_Priority, Review_PriorityData] = Assign_Review_Priority( ...
    designation, focusClipped, sfData, ubData, Assign_Review_Priority_Options)
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
% Assign_Review_Priority  Assigns a review priority (Low / Medium / High) to a trace.
%
% Three-pass logic — priority can only be raised, never lowered:
%   Pass 1 — Assign everything Low (default).
%   Pass 2 — Raise qualifying traces to Medium.
%   Pass 3 — Raise qualifying traces to High, in ordered rules.
%             First matching rule wins (early exit). This ensures that
%             more specific rules take precedence over general ones.
%
% Inputs:
%   designation   - string; final classification from Classify_Trace
%   focusClipped  - focus frame indices in the SAME coordinate space as
%                   sfData/ubData's ScoreArray/SmoothedTrace/FuseFramesClipped
%                   (i.e. Seek_Fusion/Seek_Unbounds's internally-trimmed
%                   space, index 1 == clip(ClipWidth+1) -- NOT the raw
%                   focusClipped computed against the full `clip` in
%                   Analyze_Trace_Data.m; that must be shifted by -ClipWidth
%                   first). [] if none.
%
%   sfData   - struct built in Analyze_Trace_Data step 4c; fields:
%                .ScoreArray, .DiffTrace, .FracTrace, .T5Ratio,
%                .SmoothedTrace, .FusionClusterCount, .ClusterSizes,
%                .FuseFramesClipped, .FuseFramesGlobal, .Options
%
%   ubData   - struct built in Analyze_Trace_Data step 4c; fields:
%                .ScoreArray, .DiffTrace, .FracTrace,
%                .SmoothedTrace, .UnboundClusterCount,
%                .UnboundFramesClipped, .UnboundFramesGlobal, .Options
%
%   Assign_Review_Priority_Options  - struct from Parameters_Assign_Review_Priority
%
% Outputs:
%   Review_Priority      - 'Low', 'Medium', or 'High'
%   Review_PriorityData  - struct recording which rule triggered the assigned Review_Priority

% =========================================================================
% PASS 1 — Assign Low to every trace (default)
% =========================================================================

Review_Priority = 'Low';

Review_PriorityData         = struct();
Review_PriorityData.Pass    = 1;
Review_PriorityData.Reason  = 'Default — no elevating criteria met';

% =========================================================================
% PASS 2 — Raise to Medium
%
% Traces with a positive event designation are elevated to Medium.
% The reason recorded specifies which positive designation was assigned:
%   '1 Fuse'  — algorithm detected a single fusion event
%   '2 Fuse'  — algorithm detected two fusion events
%   'Unbound' — algorithm detected an unbinding event
% =========================================================================

if ismember(designation, {'1 Fuse', '2 Fuse', 'Unbound'})
    Review_Priority              = 'Medium';
    Review_PriorityData.Pass     = 2;
    Review_PriorityData.Reason   = ['Positive event designation: ' designation];
end

% =========================================================================
% PASS 3 — Raise to High (ordered rules — first match wins)
%
% Rules are evaluated in order. The first rule that fires assigns High and
% returns immediately — more specific rules must come before general ones.
% =========================================================================

boundaryZone = Assign_Review_Priority_Options.H1.BoundaryZone;

% -------------------------------------------------------------------------
% H1 — Event frame within boundary zone of trace edge
%
% A fusion or unbound frame that falls in the first or last 20 frames of
% the clipped trace is likely unreliable: the algorithm may have latched
% onto an artefact at the start before the virus fully settles, or onto
% signal decay at the very end of the recording. Any such trace gets the
% highest review priority regardless of how cleanly the other tests passed.
%
% Applies to: '1 Fuse', '2 Fuse', 'Unbound'
% For '2 Fuse': fires if ANY fusion frame is in the boundary zone.
% -------------------------------------------------------------------------

if ismember(designation, {'1 Fuse', '2 Fuse'}) && ~isempty(sfData.FuseFramesClipped)

    clipLength   = numel(sfData.SmoothedTrace);
    triggerFrame = [];
    edgeLocation = '';

    for f = sfData.FuseFramesClipped
        if f <= boundaryZone
            triggerFrame = f;
            edgeLocation = 'start';
            break;
        elseif f >= (clipLength - boundaryZone + 1)
            triggerFrame = f;
            edgeLocation = 'end';
            break;
        end
    end

    if ~isempty(triggerFrame)
        Review_Priority                    = 'High';
        Review_PriorityData.Pass           = 3;
        Review_PriorityData.Rule           = 'H1';
        Review_PriorityData.Reason         = 'H1: Fusion frame within boundary zone of trace edge';
        Review_PriorityData.TriggerFrame   = triggerFrame;
        Review_PriorityData.EdgeLocation   = edgeLocation;
        return;
    end

end

if strcmp(designation, 'Unbound') && ~isempty(ubData.UnboundFramesClipped)

    clipLength   = numel(ubData.SmoothedTrace);
    triggerFrame = [];
    edgeLocation = '';

    for f = ubData.UnboundFramesClipped
        if f <= boundaryZone
            triggerFrame = f;
            edgeLocation = 'start';
            break;
        elseif f >= (clipLength - boundaryZone + 1)
            triggerFrame = f;
            edgeLocation = 'end';
            break;
        end
    end

    if ~isempty(triggerFrame)
        Review_Priority                    = 'High';
        Review_PriorityData.Pass           = 3;
        Review_PriorityData.Rule           = 'H1';
        Review_PriorityData.Reason         = 'H1: Unbound frame within boundary zone of trace edge';
        Review_PriorityData.TriggerFrame   = triggerFrame;
        Review_PriorityData.EdgeLocation   = edgeLocation;
        return;
    end

end

% -------------------------------------------------------------------------
% H2a — 'No Fusion' trace where T5 focus veto was the only failing test
%
% At least one frame in sfData.ScoreArray scored 15 (T1+T2+T3+T4 all
% passed) but not 31 (T5 fired because a focus event in the look-forward
% window had a large enough jump ratio). The trace would have been
% classified as fusion if not for the focus veto. Worth High review because
% the veto may be a false positive, or genuine fusion coincided with a focus
% event.
%
% Applies to: 'No Fusion' only. ('1 Fuse'/'2 Fuse' traces already have
% score-31 clusters — score-15 frames alongside them are not the deciding
% factor for those designations.)
% -------------------------------------------------------------------------

if strcmp(designation, 'No Fusion')

    score15Frames = find(sfData.ScoreArray == 15);

    if ~isempty(score15Frames)
        Review_Priority                   = 'High';
        Review_PriorityData.Pass          = 3;
        Review_PriorityData.Rule          = 'H2a';
        Review_PriorityData.Reason        = 'H2a: No Fusion but T5 focus veto was the only failing test';
        Review_PriorityData.Score15Frames = score15Frames;
        return;
    end

end

% -------------------------------------------------------------------------
% H2b — '1 Fuse' or '2 Fuse' with a focus event near the fusion frame
%
% A focus event in focusClipped falls within FocusProximityWindow frames
% of one of the recorded fusion frames (sfData.FuseFramesClipped — the
% specific frames returned by FindFusionFrame, not all score-31 frames).
% T5 did not veto the fusion candidate (ratio was below threshold), but the
% proximity of a focus artifact still makes the fusion frame timing suspect.
%
% For '2 Fuse': fires if ANY fusion frame is close to a focus event.
% -------------------------------------------------------------------------

if ismember(designation, {'1 Fuse', '2 Fuse'}) && ...
        ~isempty(focusClipped) && ~isempty(sfData.FuseFramesClipped)

    focusWindow = Assign_Review_Priority_Options.H2.FocusProximityWindow;
    triggerFusionFrame = [];
    triggerFocusFrame  = [];

    for fusionFrame = sfData.FuseFramesClipped
        nearFocus = focusClipped( ...
            abs(focusClipped - fusionFrame) <= focusWindow);
        if ~isempty(nearFocus)
            triggerFusionFrame = fusionFrame;
            triggerFocusFrame  = nearFocus(1);
            break;
        end
    end

    if ~isempty(triggerFusionFrame)
        Review_Priority                        = 'High';
        Review_PriorityData.Pass               = 3;
        Review_PriorityData.Rule               = 'H2b';
        Review_PriorityData.Reason             = 'H2b: Fusion frame within FocusProximityWindow of a focus event';
        Review_PriorityData.TriggerFusionFrame = triggerFusionFrame;
        Review_PriorityData.TriggerFocusFrame  = triggerFocusFrame;
        return;
    end

end

% -------------------------------------------------------------------------
% H3 — '1 Fuse' or '2 Fuse' with a small score-31 cluster
%
% A short cluster means only a handful of frames passed all five SeekFusion
% tests. This is more consistent with a coincidental noise spike than a
% genuine dequench rise. Empirically ~40% of false positives fall into this
% bucket while only ~10% of true fusions do (Observed June 2026 by Matt).
%
% ClusterSizes comes from Find_Fusion_Frame, one value per cluster.
% For '2 Fuse': fires if ANY cluster is small.
% -------------------------------------------------------------------------

if ismember(designation, {'1 Fuse', '2 Fuse'}) && ~isempty(sfData.ClusterSizes)

    maxSmall  = Assign_Review_Priority_Options.H3.MaxSmallClusterSize;
    smallMask = sfData.ClusterSizes <= maxSmall;

    if any(smallMask)
        Review_Priority                     = 'High';
        Review_PriorityData.Pass            = 3;
        Review_PriorityData.Rule            = 'H3';
        Review_PriorityData.Reason          = 'H3: Fusion cluster has small cluster size';
        Review_PriorityData.ClusterSizes    = sfData.ClusterSizes;
        Review_PriorityData.SmallClusterIdx = find(smallMask);
        return;
    end

end

% -------------------------------------------------------------------------
% H4 — 'No Fusion' trace with a one-test-failure frame (either algorithm)
%
% At least one frame in sfData.ScoreArray or ubData.ScoreArray reached all
% but one test of a full candidate score. Worth High review: a borderline
% parameter threshold may be causing a genuine event to be missed.
%
% SeekFusion one-test-failure scores (T1 required; exactly one of T2–T4 failed;
% T5 one-failure score = 15 already caught by H2a above):
%   Score 29 = T1+T3+T4+T5: T2 (fractional rise) failed
%   Score 27 = T1+T2+T4+T5: T3 (sustained elevation) failed
%   Score 23 = T1+T2+T3+T5: T4 (pre-rise level) failed
%
% SeekUnbounds one-test-failure scores (T1 required; exactly one of T2–T4 failed):
%   Score 13 = T1+T3+T4: T2 (fractional drop) failed
%   Score 11 = T1+T2+T4: T3 (sustained floor) failed
%   Score 7  = T1+T2+T3: T4 (pre-drop elevation) failed
%
% Both algorithms are checked. Either or both can trigger this rule.
% Applies to: 'No Fusion' only.
% -------------------------------------------------------------------------

if strcmp(designation, 'No Fusion')

    sfH4.Triggered     = false;
    sfH4.Score29Frames = find(sfData.ScoreArray == 29); % T2 (fractional rise) failed
    sfH4.Score27Frames = find(sfData.ScoreArray == 27); % T3 (sustained elevation) failed
    sfH4.Score23Frames = find(sfData.ScoreArray == 23); % T4 (pre-rise level) failed
    if ~isempty(sfH4.Score29Frames) || ~isempty(sfH4.Score27Frames) || ~isempty(sfH4.Score23Frames)
        sfH4.Triggered = true;
    end

    ubH4.Triggered     = false;
    ubH4.Score13Frames = find(ubData.ScoreArray == 13); % T2 (fractional drop) failed
    ubH4.Score7Frames  = find(ubData.ScoreArray == 7);  % T4 (pre-drop elevation) failed

    % T3 (sustained floor) near-misses are filtered by violation count.
    % Only frames where t3ViolationCount <= UB_T3MaxViolations are true near-misses.
    ub_t3MaxViol       = Assign_Review_Priority_Options.H4.UB_T3MaxViolations;
    score11Candidates  = find(ubData.ScoreArray == 11);
    if ~isempty(score11Candidates) && isfield(ubData, 'T3ViolationCount')
        violCounts         = ubData.T3ViolationCount(score11Candidates);
        score11Candidates  = score11Candidates(violCounts <= ub_t3MaxViol);
    end
    ubH4.Score11Frames = score11Candidates;

    if ~isempty(ubH4.Score13Frames) || ~isempty(ubH4.Score11Frames) || ~isempty(ubH4.Score7Frames)
        ubH4.Triggered = true;
    end

    if sfH4.Triggered || ubH4.Triggered
        triggered = {};
        if sfH4.Triggered,  triggered{end+1} = 'SeekFusion';   end
        if ubH4.Triggered,  triggered{end+1} = 'SeekUnbounds'; end

        Review_Priority                      = 'High';
        Review_PriorityData.Pass             = 3;
        Review_PriorityData.Rule             = 'H4';
        Review_PriorityData.Reason           = ['H4: No Fusion — one-test failure in: ' strjoin(triggered, ', ')];
        Review_PriorityData.SeekFusionH4     = sfH4;
        Review_PriorityData.SeekUnboundsH4   = ubH4;
        return;
    end

end

% -------------------------------------------------------------------------
% H5 — Designation is 'Other'
%
% Any trace that was designated as 'Other' is automatically flagged for High review. 
%
% Applies to: 'Other' only.
% -------------------------------------------------------------------------

if strcmp(designation, 'Other')
    Review_Priority              = 'High';
    Review_PriorityData.Pass     = 3;
    Review_PriorityData.Rule     = 'H5';
    Review_PriorityData.Reason   = 'H5: Designation is Other';
    return;
end

% [Additional High rules (H6, ...) to be implemented here]

end

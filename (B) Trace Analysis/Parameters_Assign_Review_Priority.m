function Options = Parameters_Assign_Review_Priority()
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
% Parameters_Assign_Review_Priority  Parameters for Assign_Review_Priority.
%
% To override a value:
%   Options = Parameters_Assign_Review_Priority();
%   Options.H1.BoundaryZone = 30;
%
% Loaded automatically as Options.AssignReviewPriority by Setup_Options.m,
% the overall pipeline options entry point used by Analyze_Trace_Data.m.


% =========================================================================
% H1 — BOUNDARY ZONE (Pass 3, highest priority rule)
%
% Fusion and unbound events that occur very near the start or end of the
% clipped trace are flagged as High priority. Events this close to the
% edge are unreliable: there may not be enough baseline before the event
% to confirm it was a genuine rise/drop, or the recording may have ended
% before the signal fully resolved.
%
% BoundaryZone: number of frames from either edge of the clipped trace
% within which an event frame triggers High priority.
% =========================================================================

Options.H1.BoundaryZone = 20;


% =========================================================================
% H2b — FOCUS PROXIMITY WINDOW (Pass 3, High priority)
%
% A '1 Fuse' or '2 Fuse' trace is flagged when a focus event falls very
% close to the recorded fusion frame. T5 may not have vetoed the fusion
% candidate (the ratio was below threshold), but a nearby focus event still
% raises doubt about whether the fusion frame localization is reliable.
%
% FocusProximityWindow: one-sided half-window in frames. A focus event
% triggers the flag if it falls within this many frames of a fusion frame
% on either side. Total window width = 2 × FocusProximityWindow + 1.
% e.g. FocusProximityWindow = 2  →  window of 5 frames centred on the
% fusion frame: [fusionFrame−2, fusionFrame−1, fusionFrame, fusionFrame+1, fusionFrame+2]
% =========================================================================

Options.H2.FocusProximityWindow = 2;


% =========================================================================
% H3 — SMALL FUSION CLUSTER SIZE (Pass 3, High priority)
%
% Fusion traces where any score-31 cluster spans only a few frames are
% flagged as High priority. A short cluster is a weak signal: only a
% handful of frames passed all five SeekFusion tests, which is more
% consistent with a coincidental noise spike than a real dequench rise.
% Empirically ~40% of false positives fall into this category.
%
% MaxSmallClusterSize: clusters with this many frames or fewer trigger the
% flag. e.g. 3 = a cluster of 3 consecutive score-31 frames.
% =========================================================================

Options.H3.MaxSmallClusterSize = 3;


% =========================================================================
% H4 — ONE-TEST-FAILURE NEAR-MISS (Pass 3, High priority)
%
% Controls how selective H4 is for SeekUnbounds T3 failures (score 11:
% T1+T2+T4 passed but T3 — sustained post-drop floor — failed).
%
% T3 scans FloorLength (500) frames after the drop and checks whether they
% all stay below the postDropCeiling. A trace that barely failed (e.g. 3
% frames above the ceiling out of 500) is a genuine near-miss. A trace
% that clearly failed (e.g. 300 frames above) is not a near-miss and should
% not be flagged as High priority.
%
% UB_T3MaxViolations: score-11 frames are only counted as a near-miss when
% t3ViolationCount at that frame is <= this value. Frames exceeding this
% threshold are excluded from ubH4 triggering.
%
% Tune this value based on observed H4 counts:
%   Lower = stricter (fewer High traces, only very close near-misses flagged)
%   Higher = more permissive (more High traces)
% =========================================================================

Options.H4.UB_T3MaxViolations = 25;


end

function Options = Parameters_Seek_Fusion()
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
% Parameters_Seek_Fusion
% Returns the default parameters structure for Seek_Fusion.
%
% To override a value, call this function and then modify the field:
%   Options = Parameters_Seek_Fusion();
%   Options.Test1.Threshold = 300;
%   [algoClass, scoreArray] = Seek_Fusion(trace, Options);
%
% Loaded automatically as Options.SeekFusion by Setup_Options.m, the overall
% pipeline options entry point used by Analyze_Trace_Data.m.


% =========================================================================
% DISPLAY
% =========================================================================

Options.Display.DoPlot = false;


% =========================================================================
% PRE-PROCESSING
% =========================================================================

Options.Preprocess.Method      = 'median';
Options.Preprocess.WindowWidth = 25;
Options.Preprocess.RaiseToZero = 1;



% =========================================================================
% TEST 1 — LookForward Rise
% Bit score: +1
%
% At frame i (pre-rise), look LookForward frames ahead:
%   diff = smoothed(i+LookForward) - smoothed(i)
% If this value exceeds Threshold (positive), frame i scores +1.
% Gates Tests 2 and 3 via score==3 requirement.
% =========================================================================

Options.Test1.LookForward = 15;   % frames to look forward
Options.Test1.Threshold = 200;  % raw intensity rise threshold (intensity units)


% =========================================================================
% TEST 2 — Fractional LookForward Rise
% Bit score: +2
%
% Runs on all frames where smoothed(i+LookForward) > smoothed(i).
% frac = (smoothed(i+LookForward) - smoothed(i)) / max(abs(smoothed(i)), 1)
% Denominator is smoothed(i) — the pre-rise level before the jump.
% Passes if frac > threshold, where threshold is looked up from
% Thresholds based on which 1000-unit bin smoothed(i) falls into:
%
%   Bin  1: ref in [    0,  1000)  → Thresholds(1)
%   Bin  2: ref in [ 1000,  2000)  → Thresholds(2)
%   ...
%   Bin 15: ref in [14000, 15000+) → Thresholds(15)
% =========================================================================

%                      0    1k   2k   3k   4k   5k   6k   7k   8k   9k  10k  11k  12k  13k  14k+
Options.Test2.Thresholds = [0.15,0.20,0.15,0.12,0.15,0.10,0.15,0.08,0.15,0.08,0.15,0.15,0.05,0.15,0.08];


% =========================================================================
% TEST 3 — Post-Rise Sustained Elevation
% Bit score: +4
%
% Only applied to frames that score exactly 3 (T1+T2 both fired).
% riseFloor = mean(smoothed(i), smoothed(i+LookForward))
% Condition: min(smoothed(i+LookForward : i+LookForward+CeilingLength)) >= riseFloor
% Passes when the trace stays at or above the midpoint between pre- and
% post-rise levels for CeilingLength frames after the rise. Confirms sustained.
% =========================================================================

Options.Test3.CeilingLength = 90;   % frames after the rise to check for sustained elevation


% =========================================================================
% TEST 4 — Pre-Rise Sustained Low
% Bit score: +8
%
% Only applied when T1+T2+T3 have all fired (score == 7).
% Uses the same midpoint formula as T3, scaled by PreRiseFraction:
%
%   riseFloor = mean(smoothed(i), smoothed(i+LookForward))   ← same as T3
%   ceiling   = PreRiseFraction * riseFloor
%
% PreRiseFraction = 1.0  →  ceiling = midpoint of the step (mirrors T3 exactly)
% PreRiseFraction < 1.0  →  ceiling lower (more strict — rejects more FPs)
% PreRiseFraction > 1.0  →  ceiling higher (more lenient — fewer TPs lost)
%
% Passes when max(smoothed(i-PreRiseLength : i-1)) < ceiling.
% Auto-passes at the trace boundary (fewer than 2 pre-rise frames available).
% =========================================================================

Options.Test4.PreRiseLength   = 150;   % frames before candidate to check
Options.Test4.PreRiseFraction = 1.000; % ceiling = fraction × midpoint (1.0 = same as T3 floor)


% =========================================================================
% CLASSIFICATION
% =========================================================================
%
% MinEventSeparation — minimum gap (in frames) between two clusters of
% high-scoring frames required to count them as separate fusion events
% (i.e. to call 2 Fuse instead of 1 Fuse).
%
% WHY THIS SHOULD EQUAL LookForward:
%   A single fusion event peaking at frame p produces a cluster of
%   high-scoring pre-rise frames spanning roughly p-LookForward to p-1
%   (each frame fires because it looks forward into the same rise).
%   Two events peaking at p1 and p2 produce clusters whose gap is
%   approximately (p2 - LookForward) - p1 = p2 - p1 - LookForward.
%   Setting MinEventSeparation = LookForward means a full look-forward
%   window of silence must exist between clusters before we believe they
%   are independent — protecting against noise that splits one event's
%   cluster into two pieces with a small gap. There is no physical
%   motivation for setting this to any value other than LookForward.
%   It is exposed as a parameter only for experimentation; do not change
%   it without a clear empirical reason.
% =========================================================================

Options.Classification.MinEventSeparation = Options.Test1.LookForward;


% =========================================================================
% TEST 5 — Focus Frame Veto
% Bit score: +16
%
% Applied only to frames that scored exactly 15 (passed all of T1–T4).
% For each candidate frame i, the window [i, i+LookForward] is examined.
% Any focus frame whose clipped index falls in that window is checked:
%
%   focusJump = trace(fIdx) - trace(fIdx-1)        % raw trace, not smoothed
%   windowRange = max(smoothed(i:i+LookForward)) - min(smoothed(i:i+LookForward))
%
%   If focusJump > FocusJumpThreshold * windowRange → frame is vetoed.
%
% FocusJumpThreshold is a fraction (0–1). A value of 0.5 means: if the
% intensity jump at the focus frame accounts for more than 50% of the
% total intensity swing in the detection window, the candidate is rejected
% as a focus artefact rather than a true fusion event.
%
% If no focus indices are passed to SeekFusion, all candidates auto-pass.
%
% CALLER RESPONSIBILITY — index alignment:
%   focusframenumbers in the trace struct are GLOBAL frame numbers.
%   Seek_Fusion receives a clipped trace starting at global frame
%   FrameNumFound+1 (the finding image itself, FrameNumFound, is excluded).
%   Convert before calling:  clippedIdx = globalFocusFrame - FrameNumFound
%   Failure to do this will silently apply the veto at the wrong frames.
% =========================================================================

Options.Test5.FocusJumpThreshold = 0.425;

% FocusLookExtra — extends the focus-frame search window beyond i+LookForward.
%
% Normally T5 searches for focus frames only inside [i, i+LookForward].
% This parameter adds extra frames to the RIGHT of that window:
%   focus search window = [i, i + LookForward + FocusLookExtra]
%
% Motivation: the median smoother uses a half-width of 12 frames. A focus
% spike at frame f = i+LookForward+1 (just outside the detection window)
% can still elevate smoothed(i+LookForward) through the forward half of the
% median window, causing T1/T2 to fire while T5 never sees the spike.
% A small extension (1–5 frames) catches these without a large regression.
%
% NOTE: windowRange for the ratio calculation is still computed over the
% standard [i, i+LookForward] window. Only the focus-frame search is extended.
%
% Set to 0 to reproduce the original T5 behaviour (no extension).
Options.Test5.FocusLookExtra = 1;


end

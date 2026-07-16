function Options = Parameters_Seek_Unbounds()
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
% Parameters_Seek_Unbounds
% Returns the default parameters structure for Seek_Unbounds.
%
% To override a value, call this function and then modify the field:
%   Options = Parameters_Seek_Unbounds();
%   Options.Test1.Threshold = -0.5;
%   [algoClass, scoreArray] = Seek_Unbounds(trace, Options);
%
% Loaded automatically as Options.SeekUnbounds by Setup_Options.m, the overall
% pipeline options entry point used by Analyze_Trace_Data.m.


% =========================================================================
% DISPLAY
% =========================================================================

% Plot a diagnostic figure after classification
Options.Display.DoPlot = false;


% =========================================================================
% PRE-PROCESSING
% Applied to the raw trace before any tests are run.
% =========================================================================

% Smoothing method passed to Calculate_Sliding_Window: 'median' or 'mean'
Options.Preprocess.Method = 'median';

% Full width of the sliding window (frames)
% The window extends floor(Width/2) frames on each side, including center
Options.Preprocess.WindowWidth = 3; %20

% Shift the smoothed trace up so its minimum value equals zero
% 1 = apply shift, 0 = leave as-is
Options.Preprocess.RaiseToZero = 1;


% =========================================================================
% TEST 1 — Lookback Drop (diffPast)
% Bit score: +1
%
% At each frame i, computes:
%   diffPast(i) = smoothed(i) - smoothed(i - LookBack)
% If this value falls below Threshold, the frame scores +1.
% =========================================================================

% Number of frames to look back in the diffPast computation
Options.Test1.LookBack = 8;  % R12 best

% Raw intensity drop threshold. Frames where diffPast < Threshold score +1.
% Negative values detect intensity drops (unbound behavior).
Options.Test1.Threshold = -175;  % R12 best


% =========================================================================
% TEST 2 — Fractional lookback drop
% Bit score: +2
%
% Only applied to frames that passed Test 1 (score == 1) and where
% smoothed(i) < smoothed(i - LookBack).
% Passes if EITHER condition holds:
%   fracDrop = (smoothed(i) - smoothed(i-LookBack)) / smoothed(i-LookBack) < Threshold
%   OR raw step = smoothed(i-LookBack) - smoothed(i) > StepThreshold
% =========================================================================

% Fractional drop threshold (negative). Frames where fracDrop < Threshold score +2.
Options.Test2.Threshold = -0.60;  % R12 best

% Raw intensity step threshold. Frames where the absolute drop exceeds this
% value also score +2, regardless of the fractional threshold.
Options.Test2.StepThreshold = 5000;  % R12 best (effectively disables absolute-step fallback)


% =========================================================================
% TEST 3 — Post-drop floor check
% Bit score: +4
%
% Only applied to frames that already score 3 (passed Tests 1 and 2).
% postDropCeiling = mean(smoothedTrace(i), smoothedTrace(i - LookBack)).
% Condition: max(smoothedTrace(i+1 : i+FloorLength)) <= postDropCeiling.
% Passes when the trace stays at or below that ceiling for the next
% FloorLength frames, confirming the drop is sustained.
% =========================================================================

% Number of frames forward from index i to check for non-recovery
Options.Test3.FloorLength = 500;  % R12 best



% =========================================================================
% TEST 4 — Pre-drop sustained elevation check
% Bit score: +8
%
% Only applied to frames that passed T1+T2 (bitand(score,3)==3).
% Confirms the trace was stably elevated BEFORE the drop, rejecting:
%   - Landing events (brief elevation then drop)
%   - End-of-trace drops (trace drops at the very end, T3 trivially passes)
%
% For each candidate frame i:
%   refFrame     = i - LookBack        (the pre-drop reference frame)
%   winStart     = max(refFrame - PreElevationLength, 1)
%   floorCeiling = mean(smoothed(i), smoothed(refFrame))   (same midpoint as T3)
%
%   Passes if: min(smoothed(winStart : refFrame)) >= floorCeiling
%
% Classification: score == 15 (T1+T2+T3+T5 = 1+2+4+8, all four tests) = 'Unbound'
%                 score == 11 (T1+T2+T5, T3 failed) = near-miss / landing candidate
% =========================================================================

% Number of frames to look back from the pre-drop reference frame.
% Window shrinks automatically at the start of the trace.
Options.Test5.PreElevationLength = 20;


% =========================================================================
% CLASSIFICATION
% =========================================================================
%
% MinEventSeparation — minimum gap (in frames) between two clusters of
% high-scoring frames required to count them as separate unbinding events.
% Set equal to LookBack for the same reason as SeekFusion: a single event
% produces a cluster spanning ~LookBack frames, so two clusters must be
% separated by at least LookBack frames to be considered independent.
% =========================================================================

Options.Classification.MinEventSeparation = Options.Test1.LookBack;


end

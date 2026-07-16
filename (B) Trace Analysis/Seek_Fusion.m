function [scoreArray, diffTrace, fracTrace, t5Ratio, smoothedTrace] = Seek_Fusion(trace, Options, focusIndices, clipWidth)
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
% Seek_Fusion  Scores each frame of a trace using a binary accumulator.
%             Uses forward-looking tests: frame i is the pre-rise low point,
%             and we look LookForward frames ahead to confirm a rise.
%             Each frame accumulates points when it passes detection tests.
%             Test bit scores are powers of 2 (1, 2, 4, 8, 16) for extensibility.
%             Classification (1 Fuse / 2 Fuse / No Fusion) is NOT done here —
%             that is the responsibility of the caller.
%
%             No NaN processing is performed internally. Pass a raw trace —
%             do not NaN focus frames before calling. Windowing via Calculate_Sliding_Window
%             is symmetric at boundaries (same as Seek_Unbounds). Runs off the
%             median by default (Options.Preprocess.Method = 'median').
%
%             The first/last clipWidth frames of `trace` are used only as
%             smoothing context (see Calculate_Sliding_Window) and are
%             dropped before any test runs — they never appear in any
%             output. All outputs are therefore length numel(trace) - 2*clipWidth,
%             and index 1 of every output corresponds to trace(clipWidth+1).
%
% Inputs:
%   trace        - Raw trace vector (clipped to frameStart; no NaN substitution)
%   Options      - Options struct from Parameters_Seek_Fusion.m
%   focusIndices - (optional) Focus frame indices in the SAME (pre-trim)
%                  coordinates as `trace`. Convert from global before calling:
%                      clippedIdx = globalFocusFrame - FrameNumFound
%                  (Seek_Fusion internally shifts these by -clipWidth to
%                  match the trimmed trace.) If empty or omitted, Test 5
%                  passes all candidates automatically.
%   clipWidth    - (optional) Frames to drop from each end after smoothing.
%                  Must match the value used elsewhere in the same pipeline
%                  run (Seek_Unbounds, and the +clipWidth global-frame
%                  conversion in Analyze_Trace_Data.m) — it is a single
%                  shared dead-zone parameter, not per-algorithm. Defaults
%                  to 5 if omitted, for standalone/dev usage only.
%
% Outputs:
%   scoreArray   - Integer vector (length numel(trace) - 2*clipWidth); each
%                  frame's accumulated score. Score 31 (all 5 tests pass) =
%                  fusion candidate frame.
%   diffTrace    - T1 key value: smoothed(i+LookForward) - smoothed(i) at each frame.
%                  The raw forward rise used by Test 1.
%   fracTrace    - T2 key value: fractional rise at each frame (NaN where not a rise).
%                  The normalized rise used by Test 2.
%   t5Ratio      - T5 key value: max(focusJump / windowRange) at each T1-passing frame.
%                  NaN where no focus frame is in the detection window, or T1 did not fire.
%   smoothedTrace - Preprocessed trace used by all tests.
%
% Usage:
%   Options = Parameters_Seek_Fusion();
%   [scoreArray, diffTrace, fracTrace, t5Ratio, smoothedTrace] = Seek_Fusion(trace, Options);
%   [scoreArray, diffTrace, fracTrace, t5Ratio, smoothedTrace] = Seek_Fusion(trace, Options, focusIndices);
%   [scoreArray, diffTrace, fracTrace, t5Ratio, smoothedTrace] = Seek_Fusion(trace, Options, focusIndices, clipWidth);

if nargin < 3 || isempty(focusIndices)
    focusIndices = [];
end
if nargin < 4 || isempty(clipWidth)
    clipWidth = 5;
end

% ---- Pre-processing: smooth the trace (full trace incl. clip-zone context) --
smoothedTrace = Calculate_Sliding_Window(trace, ...
    Options.Preprocess.WindowWidth, ...
    Options.Preprocess.Method, ...
    Options.Preprocess.RaiseToZero, ...
    clipWidth);

% Drop the same first/last clipWidth frames from the raw trace and focus
% indices so everything stays index-aligned with smoothedTrace.
trace = trace(clipWidth+1 : end-clipWidth);
if ~isempty(focusIndices)
    focusIndices = focusIndices - clipWidth;
    focusIndices = focusIndices(focusIndices >= 1 & focusIndices <= numel(trace));
end

n          = numel(smoothedTrace);
scoreArray = zeros(size(trace));
lf         = Options.Test1.LookForward;

% ---- Gate architecture --------------------------------------------------
% T1 is the sole gate for all downstream tests. Once T1 passes at frame i,
% T2, T3, T4, and T5 all evaluate at that frame independently of one another.
% No downstream test requires any other downstream test to have passed first.
% A frame reaches score 31 (fusion candidate) only when T1 passes and T2,
% T3, T4, and T5 also all pass at that frame.

% ---- Test 1 — LookForward Rise, bit score = 1 ------------------------
% At frame i (pre-rise), look lf frames ahead to detect a rise.
% diff = smoothed(i+lf) - smoothed(i)
% Passes if diff > Threshold. Score stored at frame i (the rise start).
% NaN at the tail (last lf frames) — boundary only, not focus-related.
diffTrace = NaN(size(smoothedTrace));
for i = 1 : n - lf
    diffTrace(i) = smoothedTrace(i + lf) - smoothedTrace(i);
    if diffTrace(i) > Options.Test1.Threshold
        scoreArray(i) = scoreArray(i) + 1;
    end
end

% ---- Test 2 — Fractional LookForward Rise, bit score = 2 -------------
% Runs on every frame that passed Test 1 (bitand(score,1)==1).
% frac = (smoothed(i+lf) - smoothed(i)) / max(abs(smoothed(i)), 1)
% Denominator is smoothed(i) — the pre-rise level (before the jump).
% Threshold is looked up from Options.Test2.Thresholds based on the
% 1000-unit bin that smoothed(i) falls into (bin 1=[0,1000), bin 15=[14000,∞)).
fracTrace = NaN(size(smoothedTrace));
for i = 1 : n - lf
    if bitand(scoreArray(i), 1) == 1
        ref = smoothedTrace(i);
        if smoothedTrace(i + lf) > ref
            fracTrace(i) = (smoothedTrace(i + lf) - ref) / max(abs(ref), 1);
            bin = min(floor(max(ref, 0) / 1000) + 1, 15);
            if fracTrace(i) > Options.Test2.Thresholds(bin)
                scoreArray(i) = scoreArray(i) + 2;
            end
        end
    end
end

% ---- Test 3 — Post-Rise Sustained Elevation, bit score = 4 -----------
% Runs on every frame that passed Test 1 (bitand(score,1)==1).
% riseFloor = mean(smoothed(i), smoothed(i+lf))
% Passes when the trace stays at or above riseFloor for CeilingLength frames
% starting from i+lf. Confirms the rise is sustained — not a transient spike.
for i = 1 : n - lf - 1
    if bitand(scoreArray(i), 1) == 1
        riseFloor = mean([smoothedTrace(i), smoothedTrace(i + lf)]);
        windowEnd = min(i + lf + Options.Test3.CeilingLength, n);
        if i + lf <= n && min(smoothedTrace(i+lf:windowEnd)) >= riseFloor
            scoreArray(i) = scoreArray(i) + 4;
        end
    end
end

% ---- Test 4 — Pre-Rise Sustained Low, bit score = 8 ------------------
% Runs on every frame that passed Test 1 (bitand(score,1)==1).
% Checks that the PreRiseLength frames immediately before frame i are all
% below a ceiling derived from the same midpoint as T3's riseFloor:
%
%   riseFloor = mean([smoothed(i), smoothed(i+lf)])   ← same formula as T3
%   ceiling   = PreRiseFraction * riseFloor
%
% PreRiseFraction = 1.0  →  ceiling = midpoint of the step (mirrors T3)
% PreRiseFraction < 1.0  →  ceiling lower (stricter — rejects more FPs)
% PreRiseFraction > 1.0  →  ceiling higher (more lenient — fewer TPs lost)
%
% Passes when max(smoothed(preStart:i-1)) < ceiling.
% Confirms the trace was genuinely low before the rise (not already elevated).
% Auto-passes when fewer than 2 pre-rise frames are available (trace boundary).
for i = 1 : n - lf - 1
    if bitand(scoreArray(i), 1) == 1
        riseFloor = mean([smoothedTrace(i), smoothedTrace(i + lf)]);
        ceiling   = Options.Test4.PreRiseFraction * riseFloor;
        preStart  = max(1, i - Options.Test4.PreRiseLength);
        if preStart <= i - 1
            if max(smoothedTrace(preStart : i-1)) < ceiling
                scoreArray(i) = scoreArray(i) + 8;
            end
        else
            scoreArray(i) = scoreArray(i) + 8;
        end
    end
end

% ---- Test 5 — Focus Frame Veto, bit score = 16 -----------------------
% Runs on every frame that passed Test 1 (bitand(score,1)==1).
% For each such frame i:
%   1. Compute windowRange = max(smoothed(i:i+lf)) - min(smoothed(i:i+lf))
%      (always the STANDARD detection window — independent of FocusLookExtra)
%   2. Find any focus indices in [i, i+lf+FocusLookExtra] (extended window).
%      FocusLookExtra=0 reproduces original [i, i+lf] behaviour.
%   3. For each focus index fIdx, compute the jump into the focus frame
%      using the RAW trace (not smoothed). Focus events are 1-2 frame
%      spikes that the median smoother suppresses:
%        focusJump = trace(fIdx) - trace(fIdx-1)
%   4. If focusJump > FocusJumpThreshold * windowRange → T5 fails; its bit
%      is not added. The frame's score stays at its current value (T1 ±
%      whatever combination of T2/T3/T4 passed at that frame).
%   5. If no focus in window, or all jumps are below threshold → T5 passes,
%      score += 16. The frame reaches 31 only if T2, T3, and T4 also passed.
%
% If focusIndices was not provided (empty), all candidates auto-pass.
%
% NOTE on index alignment:
%   focusframenumbers are GLOBAL frame numbers.
%   The caller must convert: clippedIdx = globalFocusFrame - FrameNumFound
%   before passing focusIndices (same pre-trim coordinates as `trace` —
%   Seek_Fusion applies the -clipWidth shift internally). Wrong indices
%   silently produce bad results.

candidates = find(bitand(scoreArray, 1) == 1);
t5Ratio    = NaN(size(trace));

if isempty(focusIndices)
    scoreArray(candidates) = scoreArray(candidates) + 16;
else
    for k = 1:numel(candidates)
        i           = candidates(k);
        winEnd      = min(i + lf, n);
        focusWinEnd = min(i + lf + Options.Test5.FocusLookExtra, n);
        winRange    = max(smoothedTrace(i:winEnd)) - min(smoothedTrace(i:winEnd));

        focusInWin = focusIndices(focusIndices >= i & focusIndices <= focusWinEnd);

        if isempty(focusInWin)
            scoreArray(i) = scoreArray(i) + 16;
        else
            maxRatio = 0;
            for fi = 1:numel(focusInWin)
                fIdx = focusInWin(fi);
                if fIdx > 1
                    focusJump = trace(fIdx) - trace(fIdx - 1);
                    if focusJump > 0 && winRange > 0
                        maxRatio = max(maxRatio, focusJump / winRange);
                    end
                end
            end
            t5Ratio(i) = maxRatio;

            if maxRatio <= Options.Test5.FocusJumpThreshold
                scoreArray(i) = scoreArray(i) + 16;
            end
        end
    end
end

% ---- Optional plot ------------------------------------------------------
if Options.Display.DoPlot
    meanTrace = Calculate_Sliding_Window(trace, Options.Preprocess.WindowWidth, 'mean', Options.Preprocess.RaiseToZero, 0);
    xFrames   = 1:numel(trace);

    figure(2);

    yyaxis left;
    hold on;
    plot(xFrames, trace,         '-', 'Color', 'k',         'DisplayName', 'Raw Trace');
    plot(xFrames, smoothedTrace, '-', 'Color', [0.2 0.6 1], 'LineWidth', 2, ...
        'DisplayName', sprintf('%s (w=%d)', Options.Preprocess.Method, Options.Preprocess.WindowWidth));
    plot(xFrames, meanTrace,     '-', 'Color', 'g',         'LineWidth', 2, 'DisplayName', 'Mean');
    ylabel('Intensity (Background Subtracted)');

    yyaxis right;
    plot(xFrames, diffTrace,  '-', 'Color', [1 0.5 0],   'LineWidth', 2, ...
        'DisplayName', sprintf('T1 Diff (fwd=%d)', lf));
    plot(xFrames, fracTrace,  '-', 'Color', [0 0.7 0.9], 'LineWidth', 2, ...
        'DisplayName', sprintf('T2 Frac (fwd=%d)', lf));
    plot(xFrames, scoreArray, '-', 'Color', [0.5 0 0.8], 'LineWidth', 2, 'DisplayName', 'Score');
    ylabel('Diff / Frac / Score');
    hold off;

    legend('show', 'Location', 'best');
    xlabel('Frame');
    title(sprintf('Seek_Fusion  |  Score-31 frames: %d  |  T1 Thresh: %d', sum(scoreArray >= 31), Options.Test1.Threshold));
end

end

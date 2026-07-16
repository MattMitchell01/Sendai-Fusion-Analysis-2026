function [scoreArray, diffTrace, fracTrace, smoothedTrace, t3ViolationCount] = Seek_Unbounds(trace, Options, clipWidth)
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
% Seek_Unbounds  Scores each frame of a trace using a binary accumulator.
%               Each frame accumulates points when it passes detection tests.
%               Test bit scores are powers of 2 (1, 2, 4, 8) for extensibility.
%               Classification (Unbound / No Fusion) is NOT done here —
%               that is the responsibility of the caller.
%
%               The first/last clipWidth frames of `trace` are used only as
%               smoothing context (see Calculate_Sliding_Window) and are
%               dropped before any test runs — they never appear in any
%               output. All outputs are therefore length numel(trace) - 2*clipWidth,
%               and index 1 of every output corresponds to trace(clipWidth+1).
%
% Inputs:
%   trace     - Raw trace vector (clipped to frameStart+1; no focus-frame zeroing needed)
%   Options   - Options struct from Parameters_Seek_Unbounds.m
%   clipWidth - (optional) Frames to drop from each end after smoothing.
%               Must match the value used elsewhere in the same pipeline
%               run (Seek_Fusion, and the +clipWidth global-frame conversion
%               in Analyze_Trace_Data.m) — it is a single shared dead-zone
%               parameter, not per-algorithm. Defaults to 5 if omitted, for
%               standalone/dev usage only.
%
% Outputs:
%   scoreArray      - Integer vector (length numel(trace) - 2*clipWidth); each
%                     frame's accumulated score. Score 15 (all four tests
%                     pass) = unbound candidate frame.
%   diffTrace       - T1 key value: smoothed(i) - smoothed(i-LookBack) at each frame.
%                     The raw lookback drop used by Test 1.
%   fracTrace       - T2 key value: fractional drop at each frame (NaN where not a drop).
%                     The normalized drop used by Test 2.
%   smoothedTrace   - Preprocessed trace used by all tests.
%   t3ViolationCount - Integer vector (same size as trace).
%                     At frames where T1 passed and T3 was evaluated: the number of
%                     frames in the post-drop window [i+1, i+FloorLength] that exceeded
%                     postDropCeiling. 0 means T3 passed; >0 means it failed by that
%                     many frames. NaN at all frames where T1 did not pass (T3 not run).
%                     Used by Assign_Review_Priority H4 to distinguish genuine near-misses
%                     (very few violations) from clear T3 failures (many violations).
%
% Usage:
%   Options = Parameters_Seek_Unbounds();
%   [scoreArray, diffTrace, fracTrace, smoothedTrace, t3ViolationCount] = Seek_Unbounds(trace, Options);
%   [scoreArray, diffTrace, fracTrace, smoothedTrace, t3ViolationCount] = Seek_Unbounds(trace, Options, clipWidth);

if nargin < 3 || isempty(clipWidth)
    clipWidth = 5;
end

% ---- Pre-processing: smooth the trace (full trace incl. clip-zone context) --
smoothedTrace = Calculate_Sliding_Window(trace, ...
    Options.Preprocess.WindowWidth, ...
    Options.Preprocess.Method, ...
    Options.Preprocess.RaiseToZero, ...
    clipWidth);

% Drop the same first/last clipWidth frames from the raw trace so it stays
% index-aligned with smoothedTrace.
trace = trace(clipWidth+1 : end-clipWidth);

% ---- Score accumulator --------------------------------------------------
scoreArray       = zeros(size(trace));
t3ViolationCount = NaN(size(trace));

% ---- Gate architecture --------------------------------------------------
% T1 is the sole gate for all downstream tests. Once T1 passes at frame i,
% T2, T3, and T4 all evaluate at that frame independently of one another.
% No downstream test requires any other downstream test to have passed first.
% A frame reaches score 15 (unbound candidate) only when T1 passes and T2,
% T3, and T4 also all pass at that frame.

% ---- Test 1 — Lookback drop (diffPast), bit score = 1 ---------------
% Raw intensity drop over LookBack frames. Gates T2, T3, and T4 — all three
% require T1 to have passed (bitand(score,1)==1) and are otherwise independent.
diffTrace = diffPast(smoothedTrace, Options.Test1.LookBack);

for i = 1 : numel(diffTrace)
    if ~isnan(diffTrace(i)) && diffTrace(i) < Options.Test1.Threshold
        scoreArray(i) = scoreArray(i) + 1;
    end
end

% ---- Test 2 — Fractional lookback drop, bit score = 2 ---------------
% Runs on every frame that passed Test 1 (bitand(score,1)==1). Passes if EITHER:
%   fracTrace(i) = (smoothed(i) - smoothed(i-LookBack)) / smoothed(i-LookBack) < Threshold
%   OR the raw step (smoothed(i-LookBack) - smoothed(i)) > StepThreshold
fracTrace = NaN(size(smoothedTrace));
for i = Options.Test1.LookBack + 1 : numel(smoothedTrace)
    if bitand(scoreArray(i), 1) == 1
        ref = smoothedTrace(i - Options.Test1.LookBack);
        if smoothedTrace(i) < ref
            step          = ref - smoothedTrace(i);
            fracTrace(i)  = (smoothedTrace(i) - ref) / ref;
            if (ref ~= 0 && fracTrace(i) < Options.Test2.Threshold) || step > Options.Test2.StepThreshold
                scoreArray(i) = scoreArray(i) + 2;
            end
        end
    end
end

% ---- Test 3 — Sustained post-drop floor, bit score = 4 ---------------
% Runs on every frame that passed Test 1 (bitand(score,1)==1).
% postDropCeiling = mean(smoothedTrace(i), smoothedTrace(i - LookBack)) —
%   the midpoint between the post-drop level and the pre-drop level.
% Passes when the smoothed trace stays at or below that ceiling for the
% next FloorLength frames, confirming the drop is sustained.
% t3ViolationCount(i) records how many frames in the window exceeded the
% ceiling — 0 when T3 passes, >0 when it fails.
n3 = numel(smoothedTrace);

for i = Options.Test1.LookBack + 1 : n3 - 1
    if bitand(scoreArray(i), 1) == 1
        postDropCeiling        = mean([smoothedTrace(i), smoothedTrace(i - Options.Test1.LookBack)]);
        windowEnd              = min(i + Options.Test3.FloorLength, n3);
        window                 = smoothedTrace(i+1 : windowEnd);
        nViolations            = sum(window > postDropCeiling);
        t3ViolationCount(i)    = nViolations;
        if nViolations == 0
            scoreArray(i) = scoreArray(i) + 4;
        end
    end
end

% ---- Test 4 — Pre-drop sustained elevation check, bit score = 8 ------
% Runs on every frame that passed Test 1 (bitand(score,1)==1).
% refFrame = i - LookBack (the pre-drop reference frame).
% winStart = max(refFrame - PreElevationLength, 1) — shrinks at trace start.
% elevationFloor = mean(smoothed(i), smoothed(refFrame)) — same midpoint as T3,
%   but used here as a FLOOR: the pre-drop window must stay above it.
% Passes if min(smoothed(winStart:refFrame)) >= elevationFloor.
for i = Options.Test1.LookBack + 1 : numel(smoothedTrace) - 1
    if bitand(scoreArray(i), 1) == 1
        refFrame       = i - Options.Test1.LookBack;
        winStart       = max(refFrame - Options.Test5.PreElevationLength, 1);
        elevationFloor = mean([smoothedTrace(i), smoothedTrace(refFrame)]);
        if min(smoothedTrace(winStart : refFrame)) >= elevationFloor
            scoreArray(i) = scoreArray(i) + 8;
        end
    end
end

% ---- Optional plot ------------------------------------------------------
if Options.Display.DoPlot
    meanTrace = Calculate_Sliding_Window(trace, Options.Preprocess.WindowWidth, 'mean', Options.Preprocess.RaiseToZero, 0);
    xFrames   = 1:numel(trace);

    figure(1);

    % Left y-axis: raw trace, smoothed trace, mean trace
    yyaxis left;
    hold on;
    plot(xFrames, trace,         '-', 'Color', 'k',         'DisplayName', 'Raw Trace');
    plot(xFrames, smoothedTrace, '-', 'Color', [1 0.4 0.7], 'LineWidth', 2, ...
        'DisplayName', sprintf('%s (w=%d)', Options.Preprocess.Method, Options.Preprocess.WindowWidth));
    plot(xFrames, meanTrace,     '-', 'Color', 'g',         'LineWidth', 2, 'DisplayName', 'Mean');
    ylabel('Intensity (Background Subtracted)');

    % Right y-axis: diffTrace, fracTrace (Test 2), and scoreArray
    yyaxis right;
    plot(xFrames, diffTrace,  '-', 'Color', [1 0.5 0],   'LineWidth', 2, ...
        'DisplayName', sprintf('T1 Diff (back=%d)', Options.Test1.LookBack));
    plot(xFrames, fracTrace,  '-', 'Color', [0 0.7 0.9], 'LineWidth', 2, ...
        'DisplayName', sprintf('T2 Frac (back=%d)', Options.Test1.LookBack));
    plot(xFrames, scoreArray, '-', 'Color', [0.5 0 0.8], 'LineWidth', 2, 'DisplayName', 'Score');
    ylabel('Diff / Frac / Score');
    hold off;

    legend('show', 'Location', 'best');
    xlabel('Frame');
    title(sprintf('Seek_Unbounds  |  Score==15 frames: %d  |  Threshold: %.2f', sum(scoreArray == 15), Options.Test1.Threshold));
end

end

% -------------------------------------------------------------------------
function result = diffPast(trace, back)
% diffPast  Raw intensity change looking back `back` frames.
%   result(i) = trace(i) - trace(i-back)
%   Indices where i-back < 1 are set to NaN.

n      = numel(trace);
result = NaN(size(trace));

for i = back + 1 : n
    result(i) = trace(i) - trace(i - back);
end

end

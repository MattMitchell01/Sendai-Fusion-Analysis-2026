function [unboundFrames, clusterSizes, clusterFrames] = Find_Unbound_Frame(trace, scoreArray, Options, preSmoothed)
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
% Find_Unbound_Frame  Localizes the exact unbinding frame for each event
%                     detected by Seek_Unbounds.
%
% Called after Seek_Unbounds on any trace classified as 'Unbound'.
% For each cluster of high-scoring frames, uses a step-function fit on the
% smoothed trace to find the frame index where the intensity drop actually
% occurs — i.e. the unbinding frame.
%
% Algorithm:
%   1. Find all frames with score >= 15 (score 15 = all four tests pass).
%   2. Cluster them using a 5-frame separation threshold (matching Bob's
%      NumFramesBetweenDifferentClusters in his original algorithm).
%   3. For each cluster, anchor at clusterStart and search the window
%      [clusterStart - LookBack, clusterStart] for the split point k
%      that maximizes (mean of smoothed pre-k) - (mean of smoothed post-k).
%      This is the step-function fit in the drop direction: the k that best
%      divides the window into a high pre-region and a low post-region is
%      the unbinding frame.
%
% Inputs:
%   trace       - Raw clipped trace (same vector passed to Seek_Unbounds; clip(1)
%                 == global frame FrameNumFound+1). Unused when preSmoothed
%                 is supplied (the normal case from Analyze_Trace_Data.m).
%   scoreArray  - Score array output from Seek_Unbounds (same length as trace)
%   Options     - Seek_Unbounds options struct from Parameters_Seek_Unbounds.m
%                 (uses Preprocess.WindowWidth, Preprocess.Method,
%                  Preprocess.RaiseToZero, Test1.LookBack)
%   preSmoothed - (optional) Pre-computed smoothed trace. If supplied,
%                 smoothing is skipped (avoids redundant Calculate_Sliding_Window calls).
%
% Outputs:
%   unboundFrames - Row vector of unbinding frame indices in CLIPPED coordinates.
%                   Empty if scoreArray has no score==15 frames.
%   clusterSizes  - Row vector, one element per cluster: number of score-15
%                   frames in that cluster.
%   clusterFrames - NaN-padded numeric matrix (nClusters × maxClusterSize).
%                   Row c lists clipped-coordinate indices of every score-15
%                   frame in cluster c; trailing NaNs fill shorter rows.
%
% Coordinate note:
%   Output indices share whatever coordinate space `scoreArray`/`preSmoothed`
%   were given in. When called from Analyze_Trace_Data.m, that space is
%   Seek_Unbounds's internally-trimmed trace (index 1 = global frame
%   FrameNumFound + ClipWidth + 1, since Seek_Unbounds drops the first/last
%   ClipWidth frames -- see Calculate_Sliding_Window). To convert to global
%   frame numbers, use the same offset Analyze_Trace_Data.m uses:
%       globalFrame = unboundFrames + FrameNumFound + ClipWidth

% Smooth the trace — use pre-smoothed if provided
if nargin < 4 || isempty(preSmoothed)
    smoothed = Calculate_Sliding_Window(trace, ...
        Options.Preprocess.WindowWidth, ...
        Options.Preprocess.Method, ...
        Options.Preprocess.RaiseToZero);
else
    smoothed = preSmoothed;
end

lb = Options.Test1.LookBack;

% ---- Step 1: find all unbound candidate frames --------------------------
highIdx = find(scoreArray >= 15);
if isempty(highIdx)
    unboundFrames = [];
    clusterSizes  = [];
    clusterFrames = {};
    return;
end

% ---- Step 2: cluster candidate frames -----------------------------------
% Use the same MinEventSeparation as Analyze_Trace_Data so that the number
% of clusters found here matches the ubClusterCount used for classification.
clusterSep       = Options.Classification.MinEventSeparation;
gaps             = diff(highIdx);
splitPoints      = find(gaps > clusterSep);
clusterStartsIdx = [1, splitPoints + 1];
clusterEndsIdx   = [splitPoints, numel(highIdx)];
clusterStarts    = highIdx(clusterStartsIdx);
clusterSizes     = clusterEndsIdx - clusterStartsIdx + 1;

nClusters     = numel(clusterStarts);
unboundFrames = zeros(1, nClusters);
clusterFrames = cell(1, nClusters);

% ---- Step 3: step-function fit within each cluster's window -------------
% Frame i with score>=15 means the drop occurred between i-LookBack and i.
% So the search window is [clusterStart - LookBack, clusterStart].
% We look for split point k that maximizes preMean - postMean (drop direction).
for c = 1 : nClusters
    cs       = clusterStarts(c);
    winStart = max(cs - lb, 1);
    winEnd   = cs;

    bestK    = winStart + 1;   % fallback if window is 1 frame
    bestStep = -Inf;

    for k = winStart + 1 : winEnd
        preMean  = mean(smoothed(winStart : k-1));
        postMean = mean(smoothed(k        : winEnd));
        stepH    = preMean - postMean;   % positive when there is a drop
        if stepH > bestStep
            bestStep = stepH;
            bestK    = k;
        end
    end

    unboundFrames(c)  = bestK;
    clusterFrames{c}  = highIdx(clusterStartsIdx(c) : clusterEndsIdx(c));
end

end

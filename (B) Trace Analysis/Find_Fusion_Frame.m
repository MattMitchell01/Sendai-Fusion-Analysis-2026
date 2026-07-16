function [fusionFrames, clusterSizes, clusterFrames] = Find_Fusion_Frame(trace, scoreArray, Options, preSmoothed)
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
% Find_Fusion_Frame  Localizes the exact fusion frame for each event detected
%                    by Seek_Fusion.
%
% Called after Seek_Fusion on any trace classified as '1 Fuse' or '2 Fuse'.
% For each cluster of high-scoring frames, uses a step-function fit on the
% smoothed trace to find the frame index where the intensity rise actually
% occurs — i.e. the fusion frame.
%
% Algorithm:
%   1. Find all frames with score >= 31.
%   2. Cluster them using a 5-frame separation threshold (Bob's value,
%      matching NumFramesBetweenDifferentClusters in his original algorithm).
%   3. For each cluster, anchor at clusterStart and search the window
%      [clusterStart, clusterStart + LookForward] for the split point k
%      that maximizes (mean of smoothed post-k) - (mean of smoothed pre-k).
%      This is the step-function fit: the k that best divides the window
%      into a low pre-region and a high post-region is the fusion frame.
%
% Inputs:
%   trace      - Raw clipped trace (same vector passed to SeekFusion; clip(1)
%                == global frame FrameNumFound+1). Unused when preSmoothed
%                is supplied (the normal case from Analyze_Trace_Data.m).
%   scoreArray - Score array output from SeekFusion (same length as trace)
%   Options    - Seek_Fusion options struct from Parameters_Seek_Fusion.m
%                (uses Preprocess.WindowWidth, Preprocess.Method,
%                 Preprocess.RaiseToZero, Test1.LookForward)
%
% Outputs:
%   fusionFrames - Row vector of fusion frame indices in CLIPPED coordinates.
%                  1 element for a 1 Fuse trace, 2 for a 2 Fuse trace.
%                  Empty if scoreArray has no score-31 frames.
%   clusterSizes  - Row vector, one element per cluster: the number of
%                   score-31 frames in that cluster. Used downstream by
%                   Assign_Review_Priority (H3) to flag weak-signal events.
%   clusterFrames - Cell array {1 × nClusters}. Each cell contains a row
%                   vector of clipped-coordinate frame indices for every
%                   score-31 frame in that cluster.
%                   e.g. {[271 272 273], [652 654 655]} for a 2-Fuse trace.
%
% Coordinate note:
%   Output indices share whatever coordinate space `scoreArray`/`preSmoothed`
%   were given in. When called from Analyze_Trace_Data.m, that space is
%   Seek_Fusion's internally-trimmed trace (index 1 = global frame
%   FrameNumFound + ClipWidth + 1, since Seek_Fusion drops the first/last
%   ClipWidth frames -- see Calculate_Sliding_Window). To convert to global
%   frame numbers, use the same offset Analyze_Trace_Data.m uses:
%       globalFrame = fusionFrames + FrameNumFound + ClipWidth

% Smooth the trace — use pre-smoothed if provided (avoids redundant computation
% when SeekFusion has already smoothed the same trace with the same parameters)
if nargin < 4 || isempty(preSmoothed)
    smoothed = Calculate_Sliding_Window(trace, ...
        Options.Preprocess.WindowWidth, ...
        Options.Preprocess.Method, ...
        Options.Preprocess.RaiseToZero);
else
    smoothed = preSmoothed;
end

n  = numel(smoothed);
lf = Options.Test1.LookForward;

% ---- Step 1: find all fusion candidate frames ---------------------------
highIdx = find(scoreArray >= 31);
if isempty(highIdx)
    fusionFrames  = [];
    clusterSizes  = [];
    clusterFrames = {};
    return;
end

% ---- Step 2: cluster candidate frames -----------------------------------
% Use the same MinEventSeparation as Analyze_Trace_Data so that the number
% of clusters found here matches the sfClusterCount used for classification.
clusterSep       = Options.Classification.MinEventSeparation;
gaps             = diff(highIdx);
splitPoints      = find(gaps > clusterSep);
clusterStartsIdx = [1, splitPoints + 1];          % indices into highIdx
clusterEndsIdx   = [splitPoints, numel(highIdx)]; % indices into highIdx
clusterStarts    = highIdx(clusterStartsIdx);
clusterSizes     = clusterEndsIdx - clusterStartsIdx + 1;

nClusters     = numel(clusterStarts);
fusionFrames  = zeros(1, nClusters);
clusterFrames = cell(1, nClusters);

% ---- Step 3: step-function fit within each cluster's window -------------
for c = 1 : nClusters
    cs     = clusterStarts(c);
    winEnd = min(cs + lf, n);

    bestK    = min(cs + 1, winEnd);   % fallback if window is 1 frame
    bestStep = -Inf;

    for k = cs + 1 : winEnd
        preMean  = mean(smoothed(cs : k-1));
        postMean = mean(smoothed(k  : winEnd));
        stepH    = postMean - preMean;
        if stepH > bestStep
            bestStep = stepH;
            bestK    = k;
        end
    end

    fusionFrames(c)  = bestK;
    clusterFrames{c} = highIdx(clusterStartsIdx(c) : clusterEndsIdx(c));
end

end

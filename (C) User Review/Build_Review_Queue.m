function [SortedTraceData, ReviewQueue] = Build_Review_Queue(CombinedAnalyzedTraceData, Options)
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
% Build_Review_Queue  Physically reorder traces into review groups and
% record per-group batch/grid metadata.
%
% Buckets every trace by Review_Priority into Low -> Medium{1 Fuse,2 Fuse,
% Unbound} -> High{H1,H2a,H2b,H3,H4,H5}, in that fixed order, then returns
% the trace struct array PHYSICALLY REORDERED to match -- same variable
% shape/fields as raw Part B output, only the trace order differs. Any
% trace that doesn't match an expected Designation/Rule bucket (schema
% drift, a future rule like H6, etc.) is routed to a trailing
% "Uncategorized" segment with a warning, never silently dropped.

    N = numel(CombinedAnalyzedTraceData);
    Priorities   = {CombinedAnalyzedTraceData.Review_Priority};
    Designations = {CombinedAnalyzedTraceData.Designation};

    LowIdx    = find(strcmp(Priorities, 'Low'));
    MediumIdx = find(strcmp(Priorities, 'Medium'));
    HighIdx   = find(strcmp(Priorities, 'High'));

    % ---- Review Old Data Mode: optionally filter to disagreement traces
    % only, BEFORE the Low/Medium/High/H1-H5 bucketing below -- so a
    % disagreement trace still lands in exactly the segment it would in a
    % live run (e.g. an edge-of-trace disagreement still gets H1). Nothing
    % is dropped from the dataset: traces excluded by this filter are
    % appended to their own trailing segment further down, never lost.
    if strcmp(Options.ReviewOldDataMode, 'y') && strcmp(Options.ReviewOldDataOnlyDifferences, 'y')
        IsFlagged = arrayfun(@(t) isfield(t, 'IsLegacyDisagreement') && t.IsLegacyDisagreement, ...
            CombinedAnalyzedTraceData);
        FlagIdx = find(IsFlagged);
    else
        FlagIdx = 1:N;
    end
    LowIdx    = intersect(LowIdx, FlagIdx);
    MediumIdx = intersect(MediumIdx, FlagIdx);
    HighIdx   = intersect(HighIdx, FlagIdx);

    % ---- Medium: subdivide by Designation, fixed order ----
    MediumOrder   = {'1 Fuse', '2 Fuse', 'Unbound'};
    MediumBuckets = cell(1, numel(MediumOrder));
    for k = 1:numel(MediumOrder)
        MediumBuckets{k} = MediumIdx(strcmp(Designations(MediumIdx), MediumOrder{k}));
    end
    UnmatchedMedium = setdiff(MediumIdx, [MediumBuckets{:}]);

    % ---- High: subdivide by Review_PriorityData.Rule, fixed order ----
    % Guarded isfield access -- Low/Medium traces (and, in principle, any
    % future rule-less High trace) have no .Rule field at all.
    HighOrder = {'H1', 'H2a', 'H2b', 'H3', 'H4', 'H5'};
    Rules = repmat({''}, 1, N);
    for i = HighIdx
        if isfield(CombinedAnalyzedTraceData(i).Review_PriorityData, 'Rule')
            Rules{i} = CombinedAnalyzedTraceData(i).Review_PriorityData.Rule;
        end
    end
    HighBuckets = cell(1, numel(HighOrder));
    for k = 1:numel(HighOrder)
        HighBuckets{k} = HighIdx(strcmp(Rules(HighIdx), HighOrder{k}));
    end
    UnmatchedHigh = setdiff(HighIdx, [HighBuckets{:}]);

    Unmatched = [UnmatchedMedium, UnmatchedHigh];

    % ---- Assemble master order + per-segment grid metadata ----
    Segments = struct('Tier', {}, 'SubgroupName', {}, 'StartIndex', {}, 'EndIndex', {}, ...
                       'NumTraces', {}, 'Rows', {}, 'Cols', {}, 'TracesPerBatch', {}, 'NumBatches', {});
    MasterOrder = [];

    % GridConfig for every segment is resolved via Resolve_Grid_Config.m --
    % the same lookup Start_User_Review.m's Refresh_Segment_Grids reuses on
    % every resume, so a preprocessed-then-resumed queue can never drift
    % out of sync with this initial assembly.
    [MasterOrder, Segments] = AppendSegment(MasterOrder, Segments, 'Low', 'All', LowIdx, ...
        Resolve_Grid_Config('Low', 'All', Options));
    for k = 1:numel(MediumOrder)
        [MasterOrder, Segments] = AppendSegment(MasterOrder, Segments, 'Medium', MediumOrder{k}, MediumBuckets{k}, ...
            Resolve_Grid_Config('Medium', MediumOrder{k}, Options));
    end
    for k = 1:numel(HighOrder)
        SubgroupName = HighOrder{k};
        [MasterOrder, Segments] = AppendSegment(MasterOrder, Segments, 'High', SubgroupName, HighBuckets{k}, ...
            Resolve_Grid_Config('High', SubgroupName, Options));
    end
    if ~isempty(Unmatched)
        [MasterOrder, Segments] = AppendSegment(MasterOrder, Segments, 'Uncategorized', 'Unmatched', Unmatched, ...
            Resolve_Grid_Config('Uncategorized', 'Unmatched', Options));
        warning('Build_Review_Queue:UnmatchedTraces', ...
            '%d trace(s) did not match any expected Designation/Rule bucket -- routed to a trailing Uncategorized segment.', ...
            numel(Unmatched));
    end

    % ---- Review Old Data Mode: trailing segment for traces excluded by the
    % FlagIdx filter above (i.e. everything NOT flagged .IsLegacyDisagreement
    % when Options.ReviewOldDataOnlyDifferences='y'). Tier 'NotInReviewScope'
    % doesn't match Start_User_Review.m's ActiveSegments filter
    % (Low/Medium/H1-H5 only), so this segment is walked past automatically,
    % the same mechanism that already skips 'Uncategorized' today -- no
    % changes needed there. Resolve_Grid_Config's 'otherwise' branch already
    % falls back to Options.Low's grid for any unrecognized Tier.
    ExcludedIdx = setdiff(1:N, FlagIdx);
    if ~isempty(ExcludedIdx)
        [MasterOrder, Segments] = AppendSegment(MasterOrder, Segments, 'NotInReviewScope', ...
            'ReviewOldDataFiltered', ExcludedIdx, Resolve_Grid_Config('NotInReviewScope', 'ReviewOldDataFiltered', Options));
    end

    assert(numel(MasterOrder) == N && numel(unique(MasterOrder)) == N, ...
        'Build_Review_Queue:BadPermutation', 'ReviewQueue must account for every trace exactly once.');

    SortedTraceData      = CombinedAnalyzedTraceData(MasterOrder);
    ReviewQueue.Segments = Segments;
end

function [MasterOrder, Segments] = AppendSegment(MasterOrder, Segments, Tier, SubgroupName, Idx, GridConfig)
    Idx = Idx(:)';   % preserve original file order within the group
    n = numel(Idx);

    % TracesPerBatch defaults to Rows*Cols (1 trace = 1 physical axis, true
    % for every segment except H1, which uses 2 physical axes per trace --
    % one in each of its fixed 2 rows -- see Options.HighOverrides.H1 /
    % Plot_Trace_Pair_H1.m).
    if isfield(GridConfig, 'TracesPerBatch')
        TracesPerBatch = GridConfig.TracesPerBatch;
    else
        TracesPerBatch = GridConfig.Rows * GridConfig.Cols;
    end

    Seg.Tier           = Tier;
    Seg.SubgroupName   = SubgroupName;
    Seg.StartIndex     = numel(MasterOrder) + 1;
    Seg.EndIndex       = numel(MasterOrder) + n;
    Seg.NumTraces      = n;
    Seg.Rows           = GridConfig.Rows;
    Seg.Cols           = GridConfig.Cols;
    Seg.TracesPerBatch = TracesPerBatch;
    Seg.NumBatches     = ceil(n / max(TracesPerBatch, 1));

    Segments(end+1) = Seg;
    MasterOrder = [MasterOrder, Idx];
end

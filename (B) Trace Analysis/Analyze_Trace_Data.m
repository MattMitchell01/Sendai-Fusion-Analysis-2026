function saveFilePath = Analyze_Trace_Data(dataFilePath)
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
% Analyze_Trace_Data  New Part B — classifies lipid-mixing traces using
%                     Seek_Unbounds and Seek_Fusion, then localizes fusion
%                     timing with Find_Fusion_Frame.
%
% Replacement for Bob's Start_Trace_Analysis / Analyze_Current_Data_Set.
% Input: raw Program A output file (InputTraceData format).
% Output: C/D-compatible DataToSave.CombinedAnalyzedTraceData.
%
% Inputs:
%   dataFilePath  - Full path to the .mat file from Program A. The output
%                   save path is derived automatically from this (see SAVE
%                   section below) -- this function takes no saveFilePath.
%
% Pipeline per trace:
%   1.  clip = Trace_BackSub as-is (Part A already excludes the finding image;
%       clip(1) == global frame FrameNumFound+1), convert focus indices to
%       clipped coords
%   2.  Seek_Unbounds (smooths internally; drops first/last ClipWidth frames)
%   3.  Seek_Fusion   (smooths internally; drops first/last ClipWidth frames)
%   4.  Find_Fusion_Frame (4b: Find_Unbound_Frame)
%   5.  Compute timing fields (bind-to-fusion, bind-to-unbound)
%   6.  Detect_Landing [NOT IMPLEMENTED -- local var only, not saved]
%   7.  Classify_Trace → final designation
%   8.  Assign_Review_Priority → Low/Medium/High + H1-H5 rule data
%   9.  Build output struct and save

addpath(genpath(fileparts(mfilename('fullpath'))));

fprintf('====================================\n');
fprintf('Analyze_Trace_Data\n');
fprintf('  Input:  %s\n', dataFilePath);
fprintf('====================================\n');

% =========================================================================
% LOAD OPTIONS
% =========================================================================
% Setup_Options.m is the single overall/pipeline options entry point —
% assembles the per-algorithm parameter structs (Parameters_Seek_Fusion.m,
% Parameters_Seek_Unbounds.m, Parameters_Assign_Review_Priority.m) alongside
% pipeline-wide settings (Display.DoPlot, ClipWidth).
Options = Setup_Options();

sfOpts = Options.SeekFusion;
sfOpts.Display.DoPlot = Options.Display.DoPlot;

ubOpts = Options.SeekUnbounds;
ubOpts.Display.DoPlot = Options.Display.DoPlot;

Assign_Review_Priority_Options = Options.AssignReviewPriority;

% ClipWidth — frames dropped from each end of every trace before scoring
% (their values are still used as real smoothing context; see
% Calculate_Sliding_Window). Single shared value for Seek_Fusion AND
% Seek_Unbounds, sourced once from Setup_Options.m — NOT duplicated into
% either algorithm's own parameters — so both detectors and the
% global-frame conversion below always agree on the same dead-zone boundary.
ClipWidth = Options.ClipWidth;
% hlOpts    = Setup_Options_Detect_Landing();   % [DISABLED — testing pipeline]

% =========================================================================
% LOAD INPUT DATA (Part A format only)
% =========================================================================
fprintf('Loading data...\n');
raw = load(dataFilePath);

% Accept both field names used by different versions of Program A
if isfield(raw, 'VirusDataToSave')
    Traces = raw.VirusDataToSave;
elseif isfield(raw, 'InputTraceData')
    Traces = raw.InputTraceData;
else
    error('Analyze_Trace_Data: could not find trace data. Expected VirusDataToSave or InputTraceData. Got fields: %s', ...
        strjoin(fieldnames(raw), ', '));
end

% StandardBindTime: accept both field names used by different Program A versions.
% Guard against placeholder value 0 (means not set).
fileStandardBindTime = [];
if isfield(raw, 'OtherDataToSave') && isfield(raw.OtherDataToSave, 'StandardBindTime') ...
        && raw.OtherDataToSave.StandardBindTime ~= 0
    fileStandardBindTime = raw.OtherDataToSave.StandardBindTime;
elseif isfield(raw, 'OtherImportedData') && isfield(raw.OtherImportedData, 'StandardBindTime') ...
        && raw.OtherImportedData.StandardBindTime ~= 0
    fileStandardBindTime = raw.OtherImportedData.StandardBindTime;
end

% FrameNumFound (fs) -- Part A's Find_And_Analyze_Particles.m does NOT save
% this per-trace (VirusDataToSave(n) has no such field); it's a single
% file-level value, since all traces in one Program A file share the same
% pH-drop/finding-image frame. Saved unconditionally as
% OtherDataToSave.GlobalFindingImageFrame -- read once here rather than per
% trace inside the loop below.
if ~isfield(raw, 'OtherDataToSave') || ~isfield(raw.OtherDataToSave, 'GlobalFindingImageFrame')
    error('Analyze_Trace_Data: could not find OtherDataToSave.GlobalFindingImageFrame (the file-level finding-image frame number) in this Program A output file.');
end
fs = raw.OtherDataToSave.GlobalFindingImageFrame;

% VideoTimeVector (raw/global-frame-indexed, milliseconds) is required below
% for standardBindTime/BindtoFusionTime/BindtoUnboundTime: tr.TimeVector is
% natively indexed starting at GlobalStartAnalysisFrame (fs+1), so it cannot
% be indexed with global frame numbers like fs/fusionFramesGlobal directly.
if ~isfield(raw.OtherDataToSave, 'VideoTimeVector')
    error('Analyze_Trace_Data: could not find OtherDataToSave.VideoTimeVector (the file-level, raw-frame-indexed time vector) in this Program A output file.');
end

nTotal = numel(Traces);
fprintf('  %d traces loaded.\n', nTotal);

% Keep only viruses marked as good. Any value other than 'y' is skipped.
if isfield(Traces, 'IsVirusGood')
    goodMask = strcmp({Traces.IsVirusGood}, 'y');
    Traces   = Traces(goodMask);
    fprintf('  %d good virions (IsVirusGood==''y''), %d skipped.\n\n', ...
        sum(goodMask), nTotal - sum(goodMask));
else
    fprintf('  WARNING: IsVirusGood field not found — analyzing all traces.\n\n');
end

nTr = numel(Traces);

% =========================================================================
% PER-TRACE PIPELINE
% =========================================================================
AnalyzedTraceData = [];
n1F = 0; n2F = 0; nNF = 0; nUB = 0; nOT = 0;
nLow = 0; nMed = 0; nHigh = 0;
nH1 = 0; nH2a = 0; nH2b = 0; nH3 = 0; nH4 = 0; nH5 = 0;
nH4_SF29 = 0; nH4_SF27 = 0; nH4_SF23 = 0;   % SF: T2/T3/T4 failed
nH4_UB13 = 0; nH4_UB11 = 0; nH4_UB7  = 0;   % UB: T2/T3/T4 failed
nH4_DoubleHit = 0;                            % traces with both an SF and UB near-miss

for i = 1:nTr

    tr        = Traces(i);
    fullTrace = tr.Trace_BackSub;
    % fs (FrameNumFound) is loop-invariant -- read once above from
    % OtherDataToSave.GlobalFindingImageFrame, not per-trace. Timing lookups
    % (standardBindTime/BindtoFusionTime/BindtoUnboundTime, step 5 below) use
    % the file-level OtherDataToSave.VideoTimeVector, not tr.TimeVector --
    % see step 5's comment for why.

    % ------------------------------------------------------------------
    % 1. CLIP and focus indices
    % ------------------------------------------------------------------
    % Part A now excludes BOTH the raw time-zero frame(s) AND the finding
    % image from Trace_BackSub before it ever reaches Part B.
    % StartAnalysisFrameNumber is now FrameNumToFindParticles + 1 (no longer
    % equal to it, as it used to be by convention), and Trace_BackSub is
    % built starting exactly at StartAnalysisFrameNumber (Part A's
    % Create_Video_Matrix_Auto_Time_Vector.m / Find_And_Analyze_Particles.m).
    % So Trace_BackSub's own native index 1 is ALREADY global frame
    % FrameNumFound+1 -- the finding image itself is never present in
    % Trace_BackSub at all, and no slicing is needed (or correct) to
    % exclude it.
    %
    % clip(1) == the frame immediately after the finding image, exactly as
    % it always has -- only HOW that alignment is achieved changed (Part A
    % now does it upstream, instead of Part B slicing off native index 1).
    % The global (raw-frame-number-equivalent) <-> clipped conversion is
    % therefore unchanged: clippedIdx = globalIdx - FrameNumFound /
    % globalIdx = clippedIdx + FrameNumFound -- no +/-1.
    %

    clip = fullTrace;

    % Note: `clip` here is the FULL trace as Part A now delivers it -- do NOT
    % additionally trim ClipWidth off it here. ClipWidth is a SEPARATE,
    % additional dead zone (first/last 5 frames never scored) applied
    % INSIDE Seek_Fusion/Seek_Unbounds themselves, after their own internal
    % smoothing has used those frames as real context (see
    % Calculate_Sliding_Window). Passing the full, untrimmed `clip` in is
    % what lets that smoothing context exist at all; the +ClipWidth offset
    % is applied only where fusionFramesClipped/unboundFramesClipped are
    % converted back to global frames, below.

    focusGlobal = [];
    if isfield(tr, 'focusframenumbers') && ~isempty(tr.focusframenumbers)
        ff = tr.focusframenumbers;
        focusGlobal = ff(~isnan(ff(:)))';
    elseif isfield(tr, 'FocusFrameNumbers_Shifted') && ~isempty(tr.FocusFrameNumbers_Shifted)
        ff = tr.FocusFrameNumbers_Shifted;
        % FocusFrameNumbers_Shifted = Options.FocusFrameNumbers - (StartAnalysisFrameNumber-1).
        % StartAnalysisFrameNumber-1 now equals FrameNumFound (see above), so
        % this field comes out of Part A ALREADY in clipped coordinates
        % (== focusGlobal - fs) -- it is NOT the old "0-indexed global" value
        % (that reading only held under the OLD Part A convention, where
        % StartAnalysisFrameNumber was fixed independent of FrameNumFound).
        % Add fs back here so the single focusClipped = focusGlobal - fs
        % conversion below recovers the already-clipped value exactly.
        focusGlobal = ff(~isnan(ff(:)))' + fs;
    end

    focusClipped = [];
    if ~isempty(focusGlobal)
        focusClipped = focusGlobal - fs;
        focusClipped = focusClipped(focusClipped >= 1 & focusClipped <= numel(clip));
    end

    % ------------------------------------------------------------------
    % 2. SEEK UNBOUNDS (runs on all traces)
    % ------------------------------------------------------------------
    [ubScore, ubDiff, ubFrac, ubSmoothed, ubT3ViolCount] = Seek_Unbounds(clip, ubOpts, ClipWidth);

    % Count unbound clusters: score-15 frames grouped by MinEventSeparation.
    % score 15 = T1+T2+T3+T5 (all four tests). ubClusterCount == 0 → no unbinding,
    % == 1 → Unbound, >= 2 → Other
    highUbIdx = find(ubScore >= 15);
    if isempty(highUbIdx)
        ubClusterCount = 0;
    else
        ubGaps         = diff(highUbIdx);
        ubClusterCount = sum(ubGaps > ubOpts.Classification.MinEventSeparation) + 1;
    end

    % ------------------------------------------------------------------
    % 3. SEEK FUSION (runs on all traces)
    % ------------------------------------------------------------------
    [sfScore, sfDiff, sfFrac, sfT5, sfSmoothed] = ...
        Seek_Fusion(clip, sfOpts, focusClipped, ClipWidth);

    % Count fusion clusters: score-31 frames grouped by MinEventSeparation.
    % sfClusterCount == 0 → No Fusion, == 1 → 1 Fuse, >= 2 → 2 Fuse
    highSfIdx = find(sfScore >= 31);
    if isempty(highSfIdx)
        sfClusterCount = 0;
    else
        sfGaps         = diff(highSfIdx);
        sfClusterCount = sum(sfGaps > sfOpts.Classification.MinEventSeparation) + 1;
    end

    % ------------------------------------------------------------------
    % 4. FIND FUSION FRAME (fusion traces only)
    % ------------------------------------------------------------------
    fusionFramesClipped  = [];
    fusionFramesGlobal   = [];
    fusionClusterSizes   = [];
    fusionClusterFrames  = {};

    if sfClusterCount > 0
        [fusionFramesClipped, fusionClusterSizes, fusionClusterFrames] = Find_Fusion_Frame(clip, sfScore, sfOpts, sfSmoothed);
        % sfScore/sfSmoothed (and therefore fusionFramesClipped) are indexed
        % relative to Seek_Fusion's internally-trimmed trace, whose index 1
        % is clip(ClipWidth+1) == global frame fs+1+ClipWidth.
        fusionFramesGlobal = fusionFramesClipped + fs + ClipWidth;
    end

    % ------------------------------------------------------------------
    % 4b. FIND UNBOUND FRAME (unbound traces only)
    % ------------------------------------------------------------------
    unboundFramesClipped  = [];
    unboundFramesGlobal   = [];
    unboundClusterSizes   = [];
    unboundClusterFrames  = zeros(0, 0);

    if ubClusterCount >= 1
        [unboundFramesClipped, unboundClusterSizes, unboundClusterFrames] = ...
            Find_Unbound_Frame(clip, ubScore, ubOpts, ubSmoothed);
        % Same ClipWidth offset as the fusion side above — ubScore/ubSmoothed
        % are indexed relative to Seek_Unbounds's internally-trimmed trace.
        unboundFramesGlobal = unboundFramesClipped + fs + ClipWidth;
    end

    % ------------------------------------------------------------------
    % 4c. PACKAGE ALGORITHM RESULTS INTO STRUCTS
    % Built once here and reused in step 8 (Assign_Review_Priority) and
    % step 9 (output struct) — no redundant intermediate structs.
    % ------------------------------------------------------------------
    sfData.ScoreArray          = sfScore;
    sfData.DiffTrace           = sfDiff;
    sfData.FracTrace           = sfFrac;
    sfData.T5Ratio             = sfT5;
    sfData.SmoothedTrace       = sfSmoothed;
    sfData.FusionClusterCount  = sfClusterCount;
    sfData.ClusterSizes        = fusionClusterSizes;
    sfData.ClusterFrames       = fusionClusterFrames;
    sfData.FuseFramesClipped   = fusionFramesClipped;
    sfData.FuseFramesGlobal    = fusionFramesGlobal;
    sfData.Options             = sfOpts;

    ubData.ScoreArray            = ubScore;
    ubData.DiffTrace             = ubDiff;
    ubData.FracTrace             = ubFrac;
    ubData.SmoothedTrace         = ubSmoothed;
    ubData.T3ViolationCount      = ubT3ViolCount;
    ubData.UnboundClusterCount   = ubClusterCount;
    ubData.ClusterSizes          = unboundClusterSizes;
    ubData.ClusterFrames         = unboundClusterFrames;
    ubData.UnboundFramesClipped  = unboundFramesClipped;
    ubData.UnboundFramesGlobal   = unboundFramesGlobal;
    ubData.Options               = ubOpts;

    % ------------------------------------------------------------------
    % 5. TIMING FIELDS
    % ------------------------------------------------------------------
    standardBindFrameNum = fs;

    % tr.TimeVector (timeVec) is natively indexed starting at
    % GlobalStartAnalysisFrame == fs+1 -- the finding image itself (global
    % frame fs) and everything before it are excluded upstream by Part A (see
    % the frame-indexing design notes). fs, fusionFramesGlobal, and
    % unboundFramesGlobal are all GLOBAL frame numbers, so indexing timeVec
    % with them directly reads the wrong frame's time entirely (or throws an
    % out-of-bounds error once a global frame number exceeds numel(timeVec),
    % which happens for any event detected near the end of the clip, since
    % fusionFramesGlobal/unboundFramesGlobal carry a +fs+ClipWidth offset).
    % Use the file-level VideoTimeVector instead, which IS indexed by true
    % raw/global frame number (milliseconds; convert to minutes to match
    % timeVec's units).
    if ~isempty(fileStandardBindTime)
        standardBindTime = fileStandardBindTime;
    else
        standardBindTime = raw.OtherDataToSave.VideoTimeVector(fs) / (1000 * 60);
    end

    if ~isempty(fusionFramesGlobal)
        fuseFrameNumbers      = fusionFramesGlobal;
        bindToFusionNumFrames = fuseFrameNumbers - standardBindFrameNum;
        bindToFusionTime      = raw.OtherDataToSave.VideoTimeVector(fuseFrameNumbers) / (1000 * 60) - standardBindTime;
    else
        fuseFrameNumbers      = [];
        bindToFusionNumFrames = [];
        bindToFusionTime      = [];
    end

    % Mirrors the fusion-timing block above exactly, using Find_Unbound_Frame's
    % own result (unboundFramesGlobal, step 4b) -- this IS Part B's definitive
    % unbind-frame call, computed with the same rigor as Find_Fusion_Frame's,
    % not just a diagnostic candidate.
    if ~isempty(unboundFramesGlobal)
        unboundFrameNumbers    = unboundFramesGlobal;
        bindToUnboundNumFrames = unboundFrameNumbers - standardBindFrameNum;
        bindToUnboundTime      = raw.OtherDataToSave.VideoTimeVector(unboundFrameNumbers) / (1000 * 60) - standardBindTime;
    else
        unboundFrameNumbers    = [];
        bindToUnboundNumFrames = [];
        bindToUnboundTime      = [];
    end

    % Frozen record of Part B's own bind-timing call, alongside the frame
    % indices sfData/ubData already carry (FuseFramesGlobal/FuseFramesClipped,
    % UnboundFramesGlobal/UnboundFramesClipped) -- Part C never writes to
    % SeekFusionData/SeekUnboundsData, so these values are frozen simply by
    % never being touched, the same invariant as AlgoDesignation, just
    % enforced structurally (by which struct gets mutated) instead of by a
    % naming convention within a struct that also gets actively corrected.
    sfData.BindtoFusionNumFrames  = bindToFusionNumFrames;
    sfData.BindtoFusionTime       = bindToFusionTime;
    ubData.BindtoUnboundNumFrames = bindToUnboundNumFrames;
    ubData.BindtoUnboundTime      = bindToUnboundTime;

    % ------------------------------------------------------------------
    % 6. DETECT LANDING  [NOT IMPLEMENTED]
    % ------------------------------------------------------------------
    % Landing detection is not implemented — insufficient ground-truth data.
    % landingData is kept as an empty struct only because Classify_Trace
    % inspects it; it is not saved in the output struct.
    landingData = struct();

    % ------------------------------------------------------------------
    % 7. CLASSIFY — assign final designation
    % ------------------------------------------------------------------
    designation = Classify_Trace(ubClusterCount, sfClusterCount, landingData);

    % ------------------------------------------------------------------
    % 8. ASSIGN REVIEW PRIORITY
    % ------------------------------------------------------------------
    % H2b compares focus indices directly against sfData.FuseFramesClipped,
    % which lives in Seek_Fusion's internally-trimmed coordinate space (index
    % 1 == clip(ClipWidth+1)) -- NOT the same space as focusClipped (indexed
    % against the full, untrimmed clip). Re-derive focus indices in that same
    % trimmed space here so the comparison isn't silently off by ClipWidth frames.
    focusClippedTrimmed = focusClipped - ClipWidth;
    focusClippedTrimmed = focusClippedTrimmed(focusClippedTrimmed >= 1 & focusClippedTrimmed <= numel(sfSmoothed));

    [reviewPriority, reviewPriorityData] = Assign_Review_Priority( ...
        designation, focusClippedTrimmed, sfData, ubData, Assign_Review_Priority_Options);

    % Human-readable one-line label, e.g. 'High H1' -- so a reviewer doesn't
    % have to dig into Review_PriorityData.Rule just to see which rule fired.
    % Review_Priority/Review_PriorityData themselves are untouched (several
    % places, including the summary counters below, rely on Review_Priority
    % being exactly 'Low'/'Medium'/'High').
    if strcmp(reviewPriority, 'High') && isfield(reviewPriorityData, 'Rule')
        reviewPriorityLabel = [reviewPriority ' ' reviewPriorityData.Rule];
    else
        reviewPriorityLabel = reviewPriority;
    end

    % ------------------------------------------------------------------
    % 9. BUILD OUTPUT STRUCT
    % ------------------------------------------------------------------
    out = tr;   % pass all input fields through
    if isfield(out, 'IgnoreFrameNumbers_Shifted')
        out = rmfield(out, 'IgnoreFrameNumbers_Shifted');
    end

    % FrameNumFound is NOT one of tr's real passthrough fields (Part A never
    % saves it per-trace -- see the fs lookup above) but IS a documented
    % top-level output field, since every downstream consumer (Parts C/D/E,
    % this file's own timing math) expects it per-trace. Add it explicitly.
    out.FrameNumFound = fs;

    % Top-level fields
    out.AlgoDesignation     = designation;   % frozen record of Part B's original call -- Part C corrections
                                              % update Designation/FusionData.Designation but must NEVER touch
                                              % this, so a future improved Part B can be graded against what
                                              % the reviewer actually found.
    out.Designation         = designation;   % mirrors FusionData.Designation for C/D/E compatibility
    out.ChangedByUser       = 'Not analyzed';
    out.Review_Priority      = reviewPriority;
    out.Review_PriorityLabel = reviewPriorityLabel;   % e.g. 'Low', 'Medium', 'High H1'
    out.Review_PriorityData  = reviewPriorityData;

    % FusionData — D-required fields, kept deliberately minimal: Designation,
    % FuseFrameNumbers, BindtoFusionTime only. No frozen "Algo" copies here
    % either -- the algorithm's own original call (frame indices AND
    % derived bind-timing) lives exclusively in SeekFusionData
    % (FuseFramesGlobal/FuseFramesClipped/BindtoFusionNumFrames/
    % BindtoFusionTime), which Part C never writes to, so it's frozen
    % structurally rather than by an Algo-prefixed duplicate sitting in the
    % same struct FixWaitTime.m actively corrects.
    %
    % BindtoFusionNumFrames and StandardBindFrameNum are deliberately NOT
    % included here (removed -- see the design notes on output
    % fields): StandardBindFrameNum is a pure duplicate of the
    % already-present top-level FrameNumFound (always equal to it), and
    % BindtoFusionNumFrames is a frame-count difference with no known
    % consumer -- only the TIME value (BindtoFusionTime) matters, e.g. for
    % Part D's CDF plots. Both are still computed and kept on
    % SeekFusionData below (frozen record), just not duplicated here.
    out.FusionData.Designation           = designation;   % primary field read by C/D/E
    out.FusionData.FuseFrameNumbers      = fuseFrameNumbers;
    out.FusionData.BindtoFusionTime      = bindToFusionTime;

    % UnboundData — mirrors FusionData's shape exactly, but for the unbind
    % event, so Unbound-specific data doesn't get mixed into FusionData (a
    % trace can only ever be one designation at a time; keeping the two
    % separate makes that obvious rather than implied). Find_Unbound_Frame
    % (step 4b) is Part B's own definitive unbind-frame call, computed with
    % the same rigor as Find_Fusion_Frame's -- non-empty whenever
    % Seek_Unbounds detected an unbind event at all (unboundClusterCount >=
    % 1, including on traces ultimately designated 'Other' for having 2+
    % such events), empty only if none was detected. Same as FusionData,
    % no Algo-prefixed copies here -- the frozen original lives in
    % SeekUnboundsData instead. Also same as FusionData, BindtoUnboundNumFrames
    % and StandardBindFrameNum are deliberately not duplicated here.
    out.UnboundData.UnboundFrameNumbers    = unboundFrameNumbers;
    out.UnboundData.BindtoUnboundTime      = bindToUnboundTime;

    % Full algorithm diagnostics — top-level fields
    out.SeekFusionData   = sfData;
    out.SeekUnboundsData = ubData;

    % Count by final designation
    if strcmp(designation, 'Unbound'),      nUB = nUB + 1;
    elseif strcmp(designation, '1 Fuse'),   n1F = n1F + 1;
    elseif strcmp(designation, '2 Fuse'),   n2F = n2F + 1;
    elseif strcmp(designation, 'Other'),    nOT = nOT + 1;
    else,                                   nNF = nNF + 1;
    end

    % Count by review priority
    if strcmp(reviewPriority, 'High'),        nHigh = nHigh + 1;
    elseif strcmp(reviewPriority, 'Medium'),  nMed  = nMed  + 1;
    else,                                     nLow  = nLow  + 1;
    end

    % Count High-priority rule triggers
    if strcmp(reviewPriority, 'High') && isfield(reviewPriorityData, 'Rule')
        switch reviewPriorityData.Rule
            case 'H1',  nH1  = nH1  + 1;
            case 'H2a', nH2a = nH2a + 1;
            case 'H2b', nH2b = nH2b + 1;
            case 'H3',  nH3  = nH3  + 1;
            case 'H4'
                nH4 = nH4 + 1;
                sfScores = sfData.ScoreArray;
                ubScores = ubData.ScoreArray;
                hasSF = any(sfScores == 29) || any(sfScores == 27) || any(sfScores == 23);
                hasUB = any(ubScores == 13) || any(ubScores == 11) || any(ubScores == 7);
                if any(sfScores == 29), nH4_SF29 = nH4_SF29 + 1; end
                if any(sfScores == 27), nH4_SF27 = nH4_SF27 + 1; end
                if any(sfScores == 23), nH4_SF23 = nH4_SF23 + 1; end
                if any(ubScores == 13), nH4_UB13 = nH4_UB13 + 1; end
                if any(ubScores == 11), nH4_UB11 = nH4_UB11 + 1; end
                if any(ubScores ==  7), nH4_UB7  = nH4_UB7  + 1; end
                if hasSF && hasUB,      nH4_DoubleHit = nH4_DoubleHit + 1; end
            case 'H5',  nH5  = nH5  + 1;
        end
    end

    if isempty(AnalyzedTraceData)
        AnalyzedTraceData = out;
    else
        AnalyzedTraceData(end+1) = out; %#ok<AGROW>
    end

    if mod(i, 100) == 0
        fprintf('  Processed %d / %d traces...\n', i, nTr);
    end
end

% =========================================================================
% SAVE
% =========================================================================

DataToSave.CombinedAnalyzedTraceData = AnalyzedTraceData;

% OtherDataToSave: ClipWidth + EffectiveDeadTime first, then this run's
% algorithm parameters, then an exact copy-through of everything Part A saved
% under its own OtherDataToSave (GlobalTimeZeroFrame, VideoTimeVector, etc.).
DataToSave.OtherDataToSave.ClipWidth = ClipWidth;

% EffectiveDeadTime: real-world time (MINUTES, matching every other timing
% field this pipeline saves) of the first global video frame that actually
% gets scored by Seek_Fusion/Seek_Unbounds once ClipWidth's dead zone is
% excluded -- GlobalStartAnalysisFrame is clip(1)'s global frame, so
% clip(ClipWidth+1) (the first surviving frame) is GlobalStartAnalysisFrame +
% ClipWidth. NaN for old-format Part A files that don't carry these fields.
if isfield(raw, 'OtherDataToSave') ...
        && isfield(raw.OtherDataToSave, 'GlobalStartAnalysisFrame') ...
        && isfield(raw.OtherDataToSave, 'VideoTimeVector')
    effectiveDeadTimeGlobalFrame = raw.OtherDataToSave.GlobalStartAnalysisFrame + ClipWidth;
    DataToSave.OtherDataToSave.EffectiveDeadTime = ...
        raw.OtherDataToSave.VideoTimeVector(effectiveDeadTimeGlobalFrame) / (1000 * 60);   % ms -> minutes
else
    DataToSave.OtherDataToSave.EffectiveDeadTime = NaN;   % old-format Part A file
end

DataToSave.OtherDataToSave.SeekFusionParameters     = sfOpts;
DataToSave.OtherDataToSave.SeekUnboundsParameters    = ubOpts;
DataToSave.OtherDataToSave.ReviewPriorityParameters  = Assign_Review_Priority_Options;

if isfield(raw, 'OtherDataToSave')
    passthroughFields = fieldnames(raw.OtherDataToSave);
    for k = 1:numel(passthroughFields)
        DataToSave.OtherDataToSave.(passthroughFields{k}) = raw.OtherDataToSave.(passthroughFields{k});
    end
end

% Derive label: strip _ExtractedTraces suffix, append -AnalyzedTraces (matches 26.0 Bob_Style_Save + ExtraLabel)
[inputDir, inputName, ~] = fileparts(dataFilePath);
idx = strfind(inputName, '_ExtractedTraces');
if ~isempty(idx)
    label = [inputName(1:idx(1)-1) '-AnalyzedTraces'];
else
    label = [inputName '-AnalyzedTraces'];
end

% Save into Analysis/ one level up from the input file's folder (matches 26.0 BobStyleSave='y')
% e.g. ExptFolder/Traces/file.mat  →  ExptFolder/Analysis/label.mat
exptFolder = fileparts(inputDir);
analysisSaveDir = fullfile(exptFolder, 'Analysis');
if exist(analysisSaveDir, 'dir') == 0
    mkdir(analysisSaveDir);
end
saveFilePath = fullfile(analysisSaveDir, [label '.mat']);
save(saveFilePath, 'DataToSave');

% =========================================================================
% DISPLAY — finding image with a box drawn around every analyzed trace,
% colored by final Designation. See Display_All_Analyzed_Viruses.m
% (adapted from Bob's original of the same name).
% =========================================================================
if isfield(raw, 'OtherDataToSave') && isfield(raw.OtherDataToSave, 'FindingImage')

    minImageShow = 400;   % Part A defaults (Setup_Options_Extract_Traces.m),
    maxImageShow = 800;   % used as fallback if OtherDataToSave.Options is absent
    if isfield(raw.OtherDataToSave, 'Options')
        if isfield(raw.OtherDataToSave.Options, 'MinImageShow')
            minImageShow = raw.OtherDataToSave.Options.MinImageShow;
        end
        if isfield(raw.OtherDataToSave.Options, 'MaxImageShow')
            maxImageShow = raw.OtherDataToSave.Options.MaxImageShow;
        end
    end

    % fs (== OtherDataToSave.GlobalFindingImageFrame, read once above) is
    % the same file-level finding-image frame number Options.FrameNumToFindParticles
    % holds -- use it directly rather than re-deriving it from firstTrace
    % (which, like every trace, has no FrameNumFound field of its own).
    Display_All_Analyzed_Viruses(raw.OtherDataToSave.FindingImage, AnalyzedTraceData, ...
        minImageShow, maxImageShow, fs);
else
    fprintf('  WARNING: OtherDataToSave.FindingImage not found -- skipping virus overlay figure.\n');
end

% =========================================================================
% SUMMARY
% =========================================================================
nTotal = n1F + n2F + nNF + nUB + nOT;
fprintf('\n====================================\n');
fprintf('Analyze_Trace_Data complete.\n');
fprintf('  Total:      %d\n', nTotal);
fprintf('  1 Fuse:     %d  (%.1f%%)\n', n1F, 100*n1F/max(nTotal,1));
fprintf('  2 Fuse:     %d  (%.1f%%)\n', n2F, 100*n2F/max(nTotal,1));
fprintf('  No Fusion:  %d  (%.1f%%)\n', nNF, 100*nNF/max(nTotal,1));
fprintf('  Unbound:    %d  (%.1f%%)\n', nUB, 100*nUB/max(nTotal,1));
fprintf('  Other:      %d  (%.1f%%)\n', nOT, 100*nOT/max(nTotal,1));
fprintf('\n');
fprintf('  Review Priority:\n');
fprintf('    High:     %d  (%.1f%%)\n', nHigh, 100*nHigh/max(nTotal,1));
fprintf('      H1 (event near edge):          %d\n', nH1);
fprintf('      H2a (focus-vetoed No Fusion):  %d\n', nH2a);
fprintf('      H2b (focus near fusion frame): %d\n', nH2b);
fprintf('      H3 (small cluster):            %d\n', nH3);
fprintf('      H4 (one-test-failure miss):    %d\n', nH4);
fprintf('        Fusion near-misses (SF):\n');
fprintf('          T2 failed - fractional rise  (score 29): %d\n', nH4_SF29);
fprintf('          T3 failed - sustained elev.  (score 27): %d\n', nH4_SF27);
fprintf('          T4 failed - pre-rise check   (score 23): %d\n', nH4_SF23);
fprintf('        Unbound near-misses (UB):\n');
fprintf('          T2 failed - fractional drop  (score 13): %d\n', nH4_UB13);
fprintf('          T3 failed - sustained floor  (score 11): %d\n', nH4_UB11);
fprintf('          T4 failed - pre-drop elev.   (score  7): %d\n', nH4_UB7);
if nH4_DoubleHit > 0
    fprintf('        SF+UB double-hit traces: %d\n', nH4_DoubleHit);
end
fprintf('      H5 (Other designation):        %d\n', nH5);
fprintf('    Medium:   %d  (%.1f%%)\n', nMed,  100*nMed /max(nTotal,1));
fprintf('    Low:      %d  (%.1f%%)\n', nLow,  100*nLow /max(nTotal,1));
fprintf('  Saved to: %s\n', saveFilePath);
fprintf('====================================\n');

end

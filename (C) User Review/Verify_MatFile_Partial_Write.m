function Verify_MatFile_Partial_Write()
%
% ------------------------------------------------------------------------
% Written by Matthew D. Mitchell, Rawle Lab, Williams College, 2026.
% ------------------------------------------------------------------------
%
% Verify_MatFile_Partial_Write  Standalone diagnostic -- confirms that
% matfile() partial writes into a nested struct-array field
% (DataToSave.CombinedAnalyzedTraceData, where DataToSave itself stays a
% single top-level MAT variable) actually work as expected on THIS
% MATLAB installation, BEFORE Save_User_Review_File_Fast.m is ever wired
% into the live Part C review loop.
%
% Run this yourself (type Verify_MatFile_Partial_Write at the command
% line, or press Run in the editor) and report back the PASS/FAIL result
% and the two printed timings. Nothing here touches any real User Review
% File -- it builds and deletes its own synthetic, throwaway test file in
% the system temp folder. Safe to run repeatedly.
%
% What it checks, mirroring exactly what Save_User_Review_File_Fast.m
% will do in production:
%   1. Build a synthetic DataToSave shaped like a real User Review File --
%      a modest CombinedAnalyzedTraceData struct array, plus an
%      OtherDataToSave carrying a large dummy "image" and a large dummy
%      time vector (stand-ins for the real FindingImage/VideoTimeVector).
%   2. Save it once as -v7.3, then time a FULL resave (today's approach)
%      as a baseline.
%   3. Modify ONE trace element in memory, then patch ONLY that element
%      on disk via a writable matfile handle, timed separately.
%   4. Reload fresh from disk and verify with isequaln: the patched
%      element matches exactly what was intended, EVERY other trace
%      element is byte-identical to before the patch, and
%      OtherDataToSave/ReviewOptions are completely untouched.
%
% A FAIL here means the partial-write approach is not safe to use as
% planned on this MATLAB installation -- stop and reconsider before
% touching Save_User_Review_File.m/Start_User_Review.m.

    TestFilePath = fullfile(tempdir, 'Verify_MatFile_Partial_Write_TestFile.mat');
    if isfile(TestFilePath)
        delete(TestFilePath);
    end

    fprintf('=== Verify_MatFile_Partial_Write ===\n\n');

    % ---- Build a synthetic DataToSave, shaped like the real thing ----
    NumTraces = 50;
    CombinedAnalyzedTraceData = repmat(struct( ...
        'Designation', 'No Fusion', ...
        'ChangedByUser', 'Not analyzed', ...
        'FusionData', struct('Designation', 'No Fusion', 'FuseFrameNumbers', [], 'BindtoFusionTime', []), ...
        'UnboundData', struct('UnboundFrameNumbers', [], 'BindtoUnboundTime', [])), 1, NumTraces);

    OtherDataToSave = struct();
    OtherDataToSave.FindingImage    = uint16(randi(65535, 1024, 1024));   % stand-in for the real finding image
    OtherDataToSave.VideoTimeVector = (1:200000) * 33.4;                  % stand-in for a long video's time vector

    ReviewOptions = struct('UseRunMed', 'y', 'RunMedHalfLength', 1);

    DataToSave = struct('CombinedAnalyzedTraceData', CombinedAnalyzedTraceData, ...
        'OtherDataToSave', OtherDataToSave, 'ReviewOptions', ReviewOptions);
    ReviewQueue = struct('Segments', struct('Tier', {'Low'}, 'StartIndex', {1}, 'EndIndex', {NumTraces}));
    CurrentQueuePosition = 1;
    ReviewMeta = struct('SchemaVersion', 1, 'ManualFocusSubtractIndices', [], 'ClearedIndexRanges', zeros(0,2));

    fprintf('Saving initial synthetic file (-v7.3)...\n');
    save(TestFilePath, 'DataToSave', 'ReviewQueue', 'CurrentQueuePosition', 'ReviewMeta', '-v7.3');

    OriginalDataToSave = DataToSave;   % untouched in-memory copy, for comparison later

    % ---- Baseline: time a FULL resave (today's approach) ----
    tic;
    save(TestFilePath, 'DataToSave', 'ReviewQueue', 'CurrentQueuePosition', 'ReviewMeta', '-v7.3');
    FullSaveTime = toc;
    fprintf('Full resave time:    %.3f sec\n', FullSaveTime);

    % ---- Modify ONE trace element in memory (mirrors what a single ----
    % ---- round's PlotNumber.Code correction would do)               ----
    TestIdx = 27;
    ModifiedTrace = DataToSave.CombinedAnalyzedTraceData(TestIdx);
    ModifiedTrace.Designation = '1 Fuse';
    ModifiedTrace.ChangedByUser = 'Incorrect Designation-Changed';
    ModifiedTrace.FusionData.Designation = '1 Fuse';
    ModifiedTrace.FusionData.FuseFrameNumbers = 12345;
    ModifiedTrace.FusionData.BindtoFusionTime = 4.2;

    % ---- Targeted patch: only this ONE element, via matfile ----------
    % ---- (exactly the pattern Save_User_Review_File_Fast.m uses)    ----
    tic;
    m = matfile(TestFilePath, 'Writable', true);
    m.DataToSave.CombinedAnalyzedTraceData(1, TestIdx) = ModifiedTrace;
    m.ReviewQueue = ReviewQueue;
    m.CurrentQueuePosition = CurrentQueuePosition;
    m.ReviewMeta = ReviewMeta;
    PartialSaveTime = toc;
    fprintf('Partial patch time:  %.3f sec\n\n', PartialSaveTime);

    % ---- Reload fresh from disk and verify ----
    Reloaded = load(TestFilePath, 'DataToSave');
    ReloadedTraces = Reloaded.DataToSave.CombinedAnalyzedTraceData;

    Pass = true;

    if isequaln(ReloadedTraces(TestIdx), ModifiedTrace)
        fprintf('[PASS] Patched element (%d) matches the intended modified value.\n', TestIdx);
    else
        fprintf('[FAIL] Patched element (%d) does NOT match the intended modified value.\n', TestIdx);
        Pass = false;
    end

    UntouchedIdx = setdiff(1:NumTraces, TestIdx);
    if isequaln(ReloadedTraces(UntouchedIdx), OriginalDataToSave.CombinedAnalyzedTraceData(UntouchedIdx))
        fprintf('[PASS] Every OTHER trace element (%d of them) is untouched.\n', numel(UntouchedIdx));
    else
        fprintf('[FAIL] One or more OTHER trace elements changed unexpectedly.\n');
        Pass = false;
    end

    if isequaln(Reloaded.DataToSave.OtherDataToSave, OriginalDataToSave.OtherDataToSave)
        fprintf('[PASS] OtherDataToSave (image + time vector) is completely untouched.\n');
    else
        fprintf('[FAIL] OtherDataToSave changed -- this would mean the partial write touched data it should not have.\n');
        Pass = false;
    end

    if isequaln(Reloaded.DataToSave.ReviewOptions, OriginalDataToSave.ReviewOptions)
        fprintf('[PASS] ReviewOptions is completely untouched.\n');
    else
        fprintf('[FAIL] ReviewOptions changed unexpectedly.\n');
        Pass = false;
    end

    fprintf('\n');
    if Pass
        fprintf('=== Result: PASS ===\n');
    else
        fprintf('=== Result: FAIL ===\n');
    end
    fprintf('Speedup this run: %.1fx (full resave / partial patch)\n', FullSaveTime / max(PartialSaveTime, eps));

    delete(TestFilePath);
end

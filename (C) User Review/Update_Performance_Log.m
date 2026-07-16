function Update_Performance_Log(EligibleTraces, MasterAnalysisFolder, SourceFilePath)
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
% Update_Performance_Log  Accumulates TP/Miss/FP counts (AlgoDesignation vs.
% final reviewer Designation) from this session into a running cross-session
% log saved as Performance_Log.mat in MasterAnalysisFolder, then prints a
% console report of the updated cumulative percentages.
%
% EligibleTraces has already had every trace where AlgoDesignation=='Other'
% or Designation=='Other' excluded by the caller (Start_User_Review.m) --
% this function does no further filtering.
%
% Categories: Fusion (either side is '1 Fuse' or '2 Fuse'), 1 Fuse, 2 Fuse,
% Unbound -- confirmed with the user; 'No Fusion' is intentionally not
% scored here. Formulas match Part B's own current scorer (Score_Part_B.m,
% ../Part (B)/_Dev_B) Lipid Mixing Trace Analysis/Score_Part_B.m):
%   TP    = final==cat & algo==cat
%   Miss  = final==cat & algo~=cat   (missed detection -- "FN" elsewhere)
%   FP    = final~=cat & algo==cat   (false alarm)
%   PercentDetected   = 100 * TP / (TP + Miss)
%   FalsePositiveRate = 100 * FP / TP
% Cumulative totals are always recomputed from raw counts (never by
% averaging percentages across sessions), so the running percentages stay
% mathematically correct regardless of how session sizes vary.
%
% PerformanceLog.OverallResults is a MATLAB table (one row per category,
% named columns TP/Miss/FP/PercentDetected/FalsePositiveRate) rather than a
% struct array -- a struct array's fields each collapse into one flattened
% cross-category list in the MATLAB Variable Editor/Workspace preview (e.g.
% TP shows as "[135,128,7,24]" with no visible link back to which category
% each number belongs to). A table renders as a proper one-row-per-category
% spreadsheet grid instead. Kept deliberately minimal (just these 6 columns,
% per the user) -- earlier fraction-string/denominator display columns
% (DetectedFraction, TrueEventCount, AlgoCallCount, FalseAlarmFraction) were
% tried and then dropped again in favor of this simpler set.
%
% PerformanceLog.ContributingFileLog is a struct array, one entry per
% completed review session, recording which source file and timestamp
% contributed which session-only TP/Miss/FP counts to the OverallResults
% cumulative totals above -- an audit trail behind the running totals, not
% itself the results.
%
% (Field names: OverallResults was called Categories, and ContributingFileLog
% was called History, in earlier versions -- renamed for clarity since
% neither old name described what the field actually holds. Older saved
% logs are migrated in place on load, see below.)
%
% BREAKING SCHEMA CHANGE from earlier versions (struct-array Categories;
% then FN instead of Miss; then extra fraction-string columns) -- no
% backward-compat shim for the table SHAPE, consistent with this project's
% established policy (see the design notes' Algo*-field/FusionData-trim migration
% notes). Retired individual COLUMNS are still cleaned up automatically on
% load (see below) so a reviewer doesn't have to manually delete/regenerate
% the file just to drop one field.
%
% PerformanceLog.Priority_Group_Performance is a second, independent
% breakdown -- same TP/FP arithmetic, but scoped to each review-priority
% group (Low, Medium, H1, H2a, H2b, H3, H4) instead of the whole dataset,
% so a reviewer can see which specific tier/rule the algorithm is actually
% struggling on. H5 is intentionally omitted (unspecified stat shape, and
% its traces are almost always AlgoDesignation=='Other', which
% EligibleTraces has already excluded before this function runs).

LogFilePath = fullfile(MasterAnalysisFolder, 'Performance_Log.mat');

CategoryNames = {'Fusion'; '1 Fuse'; '2 Fuse'; 'Unbound'};
NumCategories = numel(CategoryNames);

% Priority_Group_Performance groups: 'Missed' groups (Low/H2a/H4) report how
% many traces got CORRECTED INTO each category; 'FP' groups (Medium/H1/H2b/
% H3) report stats for traces that STARTED in each category (algo's own
% bucket assignment).
GroupNames      = {'Low', 'Medium', 'H1', 'H2a', 'H2b', 'H3', 'H4'};
GroupShapes     = {'Missed', 'FP', 'FP', 'Missed', 'FP', 'FP', 'Missed'};
GroupCategories = {{'1 Fuse','2 Fuse','Unbound'}, {'1 Fuse','2 Fuse','Unbound'}, {'1 Fuse','2 Fuse','Unbound'}, ...
                    {'1 Fuse','2 Fuse'}, {'1 Fuse','2 Fuse'}, {'1 Fuse','2 Fuse','Unbound'}, {'1 Fuse','2 Fuse','Unbound'}};
NumGroups = numel(GroupNames);

if isfile(LogFilePath)
    Loaded = load(LogFilePath, 'PerformanceLog');
    PerformanceLog = Loaded.PerformanceLog;

    % Older saves used 'Categories'/'History' as the field names -- rename
    % in place so an existing cumulative log carries forward under the
    % current names instead of requiring a manual delete/regenerate.
    if isfield(PerformanceLog, 'Categories') && ~isfield(PerformanceLog, 'OverallResults')
        PerformanceLog.OverallResults = PerformanceLog.Categories;
        PerformanceLog = rmfield(PerformanceLog, 'Categories');
    end
    if isfield(PerformanceLog, 'History') && ~isfield(PerformanceLog, 'ContributingFileLog')
        PerformanceLog.ContributingFileLog = PerformanceLog.History;
        PerformanceLog = rmfield(PerformanceLog, 'History');
    end

    % An older save used 'FN' as the column name -- rename in place so the
    % cumulative counts carry forward under the current name.
    if ismember('FN', PerformanceLog.OverallResults.Properties.VariableNames)
        PerformanceLog.OverallResults.Properties.VariableNames{ ...
            strcmp(PerformanceLog.OverallResults.Properties.VariableNames, 'FN')} = 'Miss';
    end

    % Retired display columns from earlier versions -- strip whichever of
    % these are present rather than requiring a manual delete/regenerate.
    RetiredColumns = {'FPperTP', 'PercentFalsePositive', 'DetectedFraction', ...
        'TrueEventCount', 'AlgoCallCount', 'FalseAlarmFraction'};
    for k = 1:numel(RetiredColumns)
        if ismember(RetiredColumns{k}, PerformanceLog.OverallResults.Properties.VariableNames)
            PerformanceLog.OverallResults.(RetiredColumns{k}) = [];
        end
    end
else
    PerformanceLog = struct();
    PerformanceLog.OverallResults = table(CategoryNames, zeros(NumCategories,1), zeros(NumCategories,1), zeros(NumCategories,1), ...
        'VariableNames', {'Category','TP','Miss','FP'});
    PerformanceLog.ContributingFileLog = struct('SourceFilePath', {}, 'Timestamp', {}, ...
        'TP', {}, 'Miss', {}, 'FP', {});
end

% Additive migration -- an older Performance_Log.mat won't have this field
% yet; initialize at zero rather than requiring a fresh file.
if ~isfield(PerformanceLog, 'Priority_Group_Performance')
    PerformanceLog.Priority_Group_Performance = struct();
    for g = 1:NumGroups
        PerformanceLog.Priority_Group_Performance.(GroupNames{g}) = ...
            Init_Group_Performance(GroupShapes{g}, GroupCategories{g});
    end
end

AlgoDesig  = {EligibleTraces.AlgoDesignation};
FinalDesig = {EligibleTraces.Designation};

SessionTP   = zeros(1, NumCategories);
SessionMiss = zeros(1, NumCategories);
SessionFP   = zeros(1, NumCategories);

for c = 1:NumCategories
    CategoryName = CategoryNames{c};
    if strcmp(CategoryName, 'Fusion')
        IsFinal = ismember(FinalDesig, {'1 Fuse', '2 Fuse'});
        IsAlgo  = ismember(AlgoDesig,  {'1 Fuse', '2 Fuse'});
    else
        IsFinal = strcmp(FinalDesig, CategoryName);
        IsAlgo  = strcmp(AlgoDesig,  CategoryName);
    end

    SessionTP(c)   = sum(IsFinal & IsAlgo);
    SessionMiss(c) = sum(IsFinal & ~IsAlgo);
    SessionFP(c)   = sum(~IsFinal & IsAlgo);

    PerformanceLog.OverallResults.TP(c)   = PerformanceLog.OverallResults.TP(c)   + SessionTP(c);
    PerformanceLog.OverallResults.Miss(c) = PerformanceLog.OverallResults.Miss(c) + SessionMiss(c);
    PerformanceLog.OverallResults.FP(c)   = PerformanceLog.OverallResults.FP(c)   + SessionFP(c);
end

% Recomputed fresh from the just-updated cumulative TP/Miss/FP every time --
% pure display columns, not additional stored state.
TP   = PerformanceLog.OverallResults.TP;
Miss = PerformanceLog.OverallResults.Miss;
FP   = PerformanceLog.OverallResults.FP;

PercentDetected   = nan(NumCategories, 1);
FalsePositiveRate = nan(NumCategories, 1);
for c = 1:NumCategories
    if (TP(c) + Miss(c)) > 0
        PercentDetected(c) = 100 * TP(c) / (TP(c) + Miss(c));
    end
    if TP(c) > 0
        FalsePositiveRate(c) = 100 * FP(c) / TP(c);
    end
end
PerformanceLog.OverallResults.PercentDetected   = PercentDetected;
PerformanceLog.OverallResults.FalsePositiveRate = FalsePositiveRate;

% ---- Priority_Group_Performance ----
% Group membership is recovered from Review_Priority/Review_PriorityData.Rule
% -- both frozen at Part B time and never touched by any Part C correction
% (see the design notes on save strategy) -- not from ReviewQueue.Segments, so this
% works regardless of file/trace order.
NumEligible = numel(EligibleTraces);
GroupKeyOfTrace = repmat({''}, 1, NumEligible);
for i = 1:NumEligible
    switch EligibleTraces(i).Review_Priority
        case 'Low'
            GroupKeyOfTrace{i} = 'Low';
        case 'Medium'
            GroupKeyOfTrace{i} = 'Medium';
        case 'High'
            if isfield(EligibleTraces(i).Review_PriorityData, 'Rule')
                GroupKeyOfTrace{i} = EligibleTraces(i).Review_PriorityData.Rule;
            end
    end
end

% A repick (.11/.22/.33) can only apply when the code's designation already
% matches the trace's live designation (Apply_Designation_Code.m's
% 'repick'-vs-'new' guard) -- so Designation==AlgoDesignation plus one of
% these ChangedByUser stamps unambiguously means "designation confirmed
% correct, but the frame/time was re-picked."
RepickStamps = {'Reviewed, Fuse frame chosen by user', ...
    'Reviewed, Unbind frame chosen by user', 'Correct Designation, Incorrect Wait Time'};
IsRepicked            = ismember({EligibleTraces.ChangedByUser}, RepickStamps);
IsCorrectDesignation  = strcmp(FinalDesig, AlgoDesig);

for g = 1:NumGroups
    GroupName  = GroupNames{g};
    Categories = GroupCategories{g};
    InGroup    = strcmp(GroupKeyOfTrace, GroupName);

    GroupPerf = PerformanceLog.Priority_Group_Performance.(GroupName);
    GroupPerf.NumTracesTotal = GroupPerf.NumTracesTotal + sum(InGroup);
    GroupPerf.TP             = GroupPerf.TP + sum(InGroup & IsCorrectDesignation);
    if GroupPerf.NumTracesTotal > 0
        GroupPerf.OverallAccuracy = 100 * GroupPerf.TP / GroupPerf.NumTracesTotal;
    end

    for c = 1:numel(Categories)
        Category   = Categories{c};
        IsAlgoCat  = strcmp(AlgoDesig, Category);
        IsFinalCat = strcmp(FinalDesig, Category);
        if strcmp(GroupShapes{g}, 'Missed')
            GroupPerf.MissedByCategory.Missed(c) = GroupPerf.MissedByCategory.Missed(c) + ...
                sum(InGroup & IsFinalCat & ~IsAlgoCat);
        else
            GroupPerf.ByCategory.NumTracesInBucket(c) = GroupPerf.ByCategory.NumTracesInBucket(c) + ...
                sum(InGroup & IsAlgoCat);
            GroupPerf.ByCategory.FP(c) = GroupPerf.ByCategory.FP(c) + ...
                sum(InGroup & IsAlgoCat & ~IsFinalCat);
            GroupPerf.ByCategory.DesignationCorrectButTimeDelta(c) = ...
                GroupPerf.ByCategory.DesignationCorrectButTimeDelta(c) + ...
                sum(InGroup & IsAlgoCat & IsFinalCat & IsRepicked);
        end
    end

    PerformanceLog.Priority_Group_Performance.(GroupName) = GroupPerf;
end

FileLogEntry = struct('SourceFilePath', SourceFilePath, 'Timestamp', char(datetime('now')), ...
    'TP', SessionTP, 'Miss', SessionMiss, 'FP', SessionFP);
PerformanceLog.ContributingFileLog(end+1) = FileLogEntry;

save(LogFilePath, 'PerformanceLog');

fprintf('\n   --- Cumulative performance log (%s) ---\n', LogFilePath);
disp(PerformanceLog.OverallResults(:, {'Category','TP','Miss','FP','PercentDetected','FalsePositiveRate'}))
disp('   (Priority_Group_Performance breakdown saved to PerformanceLog -- see the .mat file.)')

end

function GroupPerf = Init_Group_Performance(Shape, Categories)
% Zeroed starting struct for one Priority_Group_Performance group.
Categories = Categories(:);
NumCats = numel(Categories);
GroupPerf.NumTracesTotal = 0;
GroupPerf.TP = 0;
GroupPerf.OverallAccuracy = NaN;
if strcmp(Shape, 'Missed')
    GroupPerf.MissedByCategory = table(Categories, zeros(NumCats,1), ...
        'VariableNames', {'Category','Missed'});
else
    GroupPerf.ByCategory = table(Categories, zeros(NumCats,1), zeros(NumCats,1), zeros(NumCats,1), ...
        'VariableNames', {'Category','NumTracesInBucket','FP','DesignationCorrectButTimeDelta'});
end
end

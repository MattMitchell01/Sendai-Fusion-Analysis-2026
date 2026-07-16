function Score_Part_B()
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
% Score_Part_B  Compares a Part B algorithm output file against a
%               human-reviewed key file and reports per-class metrics.
%
% Opens two file-picker dialogs:
%   1. ALGORITHM OUTPUT file  —  Part B .mat output (new or Bob's)
%   2. KEY file               —  same dataset, designations set by human reviewer
%
% Both files must contain DataToSave.CombinedAnalyzedTraceData with a
% .Designation field. Traces are matched by position (same order assumed).
%
% Scored classes: 1 Fuse, 2 Fuse, Unbound, No Fusion
% Skipped from all metrics if KEY has: '' (unreviewed), 'Slow', 'Other'
%
% Metrics per class:
%   KeyCount  — traces in the key with that class
%   AlgoCount — traces the algorithm assigned that class
%   TP        — algo and key agree
%   Recall    — TP / KeyCount  ("how many real X did we catch")
%   FP        — algo said X but key disagrees
%   FP/TP     — false-alarm rate relative to true positives

fprintf('\n==============================================================\n');
fprintf('  SCORE PART B\n');
fprintf('==============================================================\n\n');

% =========================================================================
% FILE SELECTION
% =========================================================================
fprintf('Step 1: Select the ALGORITHM OUTPUT file (Part B .mat)...\n');
[algoFile, algoDir] = uigetfile('*.mat', 'Select Algorithm Output File');
if isequal(algoFile, 0)
    fprintf('Cancelled.\n');  return;
end

fprintf('Step 2: Select the KEY file (human-reviewed .mat)...\n');
[keyFile, keyDir] = uigetfile('*.mat', 'Select Key File');
if isequal(keyFile, 0)
    fprintf('Cancelled.\n');  return;
end

algoPath = fullfile(algoDir, algoFile);
keyPath  = fullfile(keyDir,  keyFile);

fprintf('\n  Algo : %s\n', algoPath);
fprintf('  Key  : %s\n\n', keyPath);

% =========================================================================
% LOAD
% =========================================================================
fprintf('Loading...\n');
algoData = loadTraces(algoPath);
keyData  = loadTraces(keyPath);
if isempty(algoData) || isempty(keyData),  return;  end

nAlgo = numel(algoData);
nKey  = numel(keyData);
fprintf('  Algorithm output : %d traces\n', nAlgo);
fprintf('  Key              : %d traces\n\n', nKey);

if nAlgo ~= nKey
    fprintf('  WARNING: counts differ. Scoring first %d (matched by position).\n\n', ...
        min(nAlgo, nKey));
end
nMatch = min(nAlgo, nKey);

% =========================================================================
% EXTRACT AND FILTER
% =========================================================================
algoDesig = extractDesig(algoData, nMatch);
keyDesig  = extractDesig(keyData,  nMatch);

skipMask = strcmp(keyDesig, '') | strcmp(keyDesig, 'Slow') | strcmp(keyDesig, 'Other');
nSkip    = sum(skipMask);
nScored  = nMatch - nSkip;

fprintf('  Matched  : %d\n', nMatch);
fprintf('  Skipped  : %d  (key is empty / Slow / Other)\n', nSkip);
fprintf('  Scored   : %d\n\n', nScored);

algoD = algoDesig(~skipMask);
keyD  = keyDesig(~skipMask);

% =========================================================================
% PER-CLASS METRICS
% =========================================================================
classes = {'1 Fuse', '2 Fuse', 'Unbound', 'No Fusion'};
nC      = numel(classes);
R       = struct('Class', classes, 'KeyN', {0}, 'AlgoN', {0}, ...
                 'TP', {0}, 'FP', {0}, 'FN', {0}, ...
                 'Recall', {NaN}, 'FPperTP', {NaN});

for c = 1:nC
    cls      = classes{c};
    isKey    = strcmp(keyD,  cls);
    isAlgo   = strcmp(algoD, cls);

    R(c).KeyN   = sum(isKey);
    R(c).AlgoN  = sum(isAlgo);
    R(c).TP     = sum(isKey  &  isAlgo);
    R(c).FP     = sum(~isKey &  isAlgo);
    R(c).FN     = sum(isKey  & ~isAlgo);

    if R(c).KeyN  > 0,  R(c).Recall  = R(c).TP / R(c).KeyN;   end
    if R(c).TP    > 0,  R(c).FPperTP = R(c).FP / R(c).TP;     end
end

nCorrect    = sum(strcmp(algoD, keyD));
overallRecall = nCorrect / max(nScored, 1);

% =========================================================================
% REPORT
% =========================================================================
div = '  ────────────────────────────────────────────────────';

fprintf('\n============================================================\n');
fprintf('  PART B ASSESSMENT  —  %d traces scored  (%d skipped)\n', nScored, nSkip);
fprintf('  Algo : %s\n', algoPath);
fprintf('  Key  : %s\n', keyPath);
fprintf('============================================================\n');

% --- Detection rate -------------------------------------------------------
fprintf('\n  DETECTION RATE  (of the real events, how many did we catch?)\n');
fprintf('%s\n', div);
for c = 1:nC
    r = R(c);
    if r.KeyN == 0
        fprintf('  %-10s  no examples in key\n', r.Class);
    else
        fprintf('  %-10s  %d / %d caught  (%s)    missed %d\n', ...
            r.Class, r.TP, r.KeyN, pctStr(r.Recall), r.FN);
    end
end
fprintf('%s\n', div);

% --- False alarms ---------------------------------------------------------
fprintf('\n  FALSE ALARMS  (algo said X, but key says it was something else)\n');
fprintf('%s\n', div);
for c = 1:nC
    r = R(c);
    if r.FP == 0
        fprintf('  %-10s  no false alarms\n', r.Class);
    elseif r.TP > 0
        fprintf('  %-10s  %d false alarms  (%.1f wrong for every 1 correct)\n', ...
            r.Class, r.FP, r.FP / r.TP);
    else
        fprintf('  %-10s  %d false alarms  (0 correct — all wrong)\n', r.Class, r.FP);
    end
end
fprintf('%s\n', div);

% --- Biggest errors -------------------------------------------------------
% Collect all error pairs (key→algo) and sort by count descending
errorPairs  = {};
errorCounts = [];
for c = 1:nC
    for d = 1:nC
        if c == d,  continue;  end
        n = sum(strcmp(keyD, classes{c}) & strcmp(algoD, classes{d}));
        if n > 0
            errorPairs{end+1}  = sprintf('%-10s called "%s"', classes{c}, classes{d}); %#ok<AGROW>
            errorCounts(end+1) = n; %#ok<AGROW>
        end
    end
    % also count key=c called Other/unrecognised
    nOther = sum(strcmp(keyD, classes{c}) & ~ismember(algoD, classes));
    if nOther > 0
        errorPairs{end+1}  = sprintf('%-10s called "Other"', classes{c}); %#ok<AGROW>
        errorCounts(end+1) = nOther; %#ok<AGROW>
    end
end

if ~isempty(errorCounts)
    [~, idx] = sort(errorCounts, 'descend');
    fprintf('\n  WHERE THE ERRORS ARE  (sorted by frequency)\n');
    fprintf('%s\n', div);
    for k = 1:numel(idx)
        fprintf('  %4d  %s\n', errorCounts(idx(k)), errorPairs{idx(k)});
    end
    fprintf('%s\n', div);
end

% --- Overall --------------------------------------------------------------
fprintf('\n  Overall: %d / %d correct  (%s)\n\n', ...
    nCorrect, nScored, pctStr(overallRecall));

end


% =========================================================================
% LOCAL HELPERS
% =========================================================================

function data = loadTraces(filePath)
data = [];
try
    raw = load(filePath);
catch e
    fprintf('  ERROR loading file: %s\n', e.message);
    return;
end

% Standard Part B output
if isfield(raw, 'DataToSave') && isfield(raw.DataToSave, 'CombinedAnalyzedTraceData')
    data = raw.DataToSave.CombinedAnalyzedTraceData;
    return;
end

% Flat struct at top level
if isfield(raw, 'CombinedAnalyzedTraceData')
    data = raw.CombinedAnalyzedTraceData;
    return;
end

% Search one level deep for any struct array with a Designation field
fnames = fieldnames(raw);
for k = 1:numel(fnames)
    top = raw.(fnames{k});
    if isstruct(top) && isfield(top, 'Designation')
        data = top;
        fprintf('  Note: loaded from top-level field "%s"\n', fnames{k});
        return;
    end
    if isstruct(top) && ~isempty(fieldnames(top))
        subnames = fieldnames(top);
        for j = 1:numel(subnames)
            sub = top.(subnames{j});
            if isstruct(sub) && isfield(sub, 'Designation')
                data = sub;
                fprintf('  Note: loaded from "%s.%s"\n', fnames{k}, subnames{j});
                return;
            end
        end
    end
end

fprintf('  ERROR: could not find trace data with a Designation field in:\n');
fprintf('         %s\n', filePath);
fprintf('  Top-level fields: %s\n', strjoin(fnames, ', '));
end


function desig = extractDesig(data, n)
% Pull .Designation from each trace; fall back to FusionData.Designation.
desig = cell(n, 1);
for i = 1:n
    d = '';
    if isfield(data(i), 'Designation') && ~isempty(data(i).Designation)
        d = data(i).Designation;
    elseif isfield(data(i), 'FusionData') && isfield(data(i).FusionData, 'Designation') ...
            && ~isempty(data(i).FusionData.Designation)
        d = data(i).FusionData.Designation;
    end
    if isnumeric(d),  d = '';  end   % guard against NaN/0 placeholders
    desig{i} = d;
end
end


function s = pctStr(v)
if isnan(v),  s = '     n/a';
else,         s = sprintf('%6.1f%%', 100 * v);
end
end


function s = ratioStr(v)
if isnan(v),  s = '    n/a';
else,         s = sprintf('%7.3f', v);
end
end

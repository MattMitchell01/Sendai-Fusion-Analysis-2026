function Flag_Legacy_Disagreements(bOutputPath, oldReviewedPath)
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
% Flag_Legacy_Disagreements  Compares the New B Output's raw algorithmic
% call against the Old User Review File's final human-reviewed designation,
% and flags every disagreement directly on the New B Output file (in place).
%
% Matching key: VirusIDNumber (verified unique within a single video file,
% and reliably shared between a raw Program A file and its corresponding
% reviewed file -- see Part (R) plan notes).
%
% Comparison field: .AlgoDesignation (New B Output's frozen, pre-review
% call) vs. the Old User Review File's .Designation (final human-reviewed
% designation) -- NOT top-level .Designation, so the comparison stays
% pinned to the frozen field even if this file is re-saved later by
% something else.
%
% The Old User Review File can contain 'Slow', which New B Output never
% outputs -- this needs no special-casing: any 'Slow'-labeled trace is
% automatically a disagreement under plain string comparison, correctly
% routing it to human review under the new pipeline's real designation set.
%
% Adds, for EVERY matched trace (agreements included, not just
% disagreements):
%   .IsLegacyDisagreement          logical
%   .LegacyComparison.OldDesignation
%   .LegacyComparison.NewAlgoDesignation
%   .LegacyComparison.Match        logical (NaN if unmatched)
%
% Only an agree/disagree percentage is reported -- deliberately no
% category-by-category breakdown (e.g. a confusion matrix). Seeing which
% specific categories the algorithm confuses before reviewing the
% disagreements themselves would bias a reviewer's expectations going in.
%
% Inputs:
%   bOutputPath     - path to a New B Output .mat (DataToSave.CombinedAnalyzedTraceData)
%   oldReviewedPath - path to the corresponding Old User Review File .mat
%                     (DataToSave.CombinedAnalyzedTraceData, old-format Designation)
%
% Saves DataToSave back to bOutputPath, overwriting it in place.

fprintf('Flag_Legacy_Disagreements:\n');
fprintf('  New B Output: %s\n', bOutputPath);
fprintf('  Old User Review File: %s\n', oldReviewedPath);

BLoaded   = load(bOutputPath, 'DataToSave');
DataToSave = BLoaded.DataToSave;
BTraces   = DataToSave.CombinedAnalyzedTraceData;

OldLoaded = load(oldReviewedPath, 'DataToSave');
OldTraces = OldLoaded.DataToSave.CombinedAnalyzedTraceData;

% ---- Build VirusIDNumber -> Designation lookup from the Old User Review File ----
OldIDs   = [OldTraces.VirusIDNumber];
OldDesigs = {OldTraces.Designation};
if numel(unique(OldIDs)) ~= numel(OldIDs)
    error('Flag_Legacy_Disagreements: VirusIDNumber is not unique in %s -- cannot build a reliable lookup.', ...
        oldReviewedPath);
end
OldMap = containers.Map(num2cell(OldIDs), OldDesigs);

BIDs = [BTraces.VirusIDNumber];
UnmatchedInB   = [];
UsedOldIDs     = false(1, numel(OldIDs));

nMatch = 0;
nMismatch = 0;
nUnmatched = 0;

for i = 1:numel(BTraces)
    vid = BIDs(i);
    algoDesig = BTraces(i).AlgoDesignation;

    if isKey(OldMap, vid)
        oldDesig = OldMap(vid);
        isMatch  = strcmp(oldDesig, algoDesig);

        BTraces(i).LegacyComparison.OldDesignation    = oldDesig;
        BTraces(i).LegacyComparison.NewAlgoDesignation = algoDesig;
        BTraces(i).LegacyComparison.Match              = isMatch;
        BTraces(i).IsLegacyDisagreement                = ~isMatch;

        if isMatch
            nMatch = nMatch + 1;
        else
            nMismatch = nMismatch + 1;
        end

        oldPos = find(OldIDs == vid, 1);
        UsedOldIDs(oldPos) = true;
    else
        BTraces(i).LegacyComparison.OldDesignation    = '';
        BTraces(i).LegacyComparison.NewAlgoDesignation = algoDesig;
        BTraces(i).LegacyComparison.Match              = NaN;
        BTraces(i).IsLegacyDisagreement                = false;   % nothing to compare against

        UnmatchedInB(end+1) = vid; %#ok<AGROW>
        nUnmatched = nUnmatched + 1;
    end
end

UnmatchedInOld = OldIDs(~UsedOldIDs);

DataToSave.CombinedAnalyzedTraceData = BTraces;
save(bOutputPath, 'DataToSave');

% ---- Report ----
fprintf('\n  --- Comparison summary ---\n');
fprintf('  Matched traces:    %d\n', nMatch + nMismatch);
fprintf('    Agreements:      %d (%.1f%%)\n', nMatch, 100*nMatch/max(nMatch+nMismatch,1));
fprintf('    Disagreements:   %d (%.1f%%)\n', nMismatch, 100*nMismatch/max(nMatch+nMismatch,1));
fprintf('  Unmatched (in New B Output, no Old User Review File match): %d\n', nUnmatched);
if ~isempty(UnmatchedInB)
    fprintf('    VirusIDNumbers: %s\n', mat2str(UnmatchedInB));
end
if ~isempty(UnmatchedInOld)
    fprintf('  WARNING: %d VirusIDNumber(s) in Old User Review File have no New B Output counterpart ', ...
        numel(UnmatchedInOld));
    fprintf('(possible video-pair mismatch): %s\n', mat2str(UnmatchedInOld));
end

fprintf('\n  Flagged file saved: %s\n', bOutputPath);

end

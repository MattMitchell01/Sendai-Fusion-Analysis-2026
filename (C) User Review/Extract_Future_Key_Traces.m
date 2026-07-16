function Extract_Future_Key_Traces(EligibleTraces, OtherDataToSave, MasterAnalysisFolder)
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
% Extract_Future_Key_Traces  Randomly samples 20% of this session's eligible
% traces and appends them to a shared, cross-experiment collection
% (Key_Traces_Collection.mat, in MasterAnalysisFolder) -- future ground-truth
% key material.
%
% EligibleTraces has already had every trace where AlgoDesignation=='Other'
% or Designation=='Other' excluded by the caller (Start_User_Review.m) --
% this function samples from that pool as-is, no further filtering.
%
% Each sampled trace gets a new OtherOptions field (no space -- MATLAB
% struct field names can't contain spaces) attached LAST, holding a copy of
% this session's DataToSave.OtherDataToSave (GlobalFindingImageFrame,
% Options, VideoTimeVector, etc. -- everything EXCEPT FindingImage, see
% below). Traces pulled out of their original per-file struct array and
% pooled with traces from other experiments/assays would otherwise lose
% that file-level context -- attaching it per-trace keeps every pooled
% trace self-describing on its own.
%
% FindingImage is deliberately dropped before attaching: it's a full
% camera-resolution image (e.g. 1024x1024 uint16, ~2MB), not useful for
% whatever this collection is actually used for, and since it's identical
% for every trace sampled from the same session, attaching it per-trace
% means every single sampled trace pays that ~2MB cost again -- this was
% confirmed to be the dominant reason Key_Traces_Collection.mat had grown
% far larger than any single review file (a review file stores hundreds of
% full traces but only ONE copy of the finding image for the whole file;
% this collection was storing one copy PER SAMPLED TRACE).

KeyFilePath = fullfile(MasterAnalysisFolder, 'Key_Traces_Collection.mat');

SampleSize = round(0.20 * numel(EligibleTraces));
if SampleSize == 0
    fprintf('\n   No traces sampled for the key-trace collection (0%% of %d eligible traces).\n', ...
        numel(EligibleTraces));
    return
end

SampledIdx = randperm(numel(EligibleTraces), SampleSize);
NewKeyTraces = EligibleTraces(SampledIdx);

OtherOptionsToSave = OtherDataToSave;
if isfield(OtherOptionsToSave, 'FindingImage')
    OtherOptionsToSave = rmfield(OtherOptionsToSave, 'FindingImage');
end
[NewKeyTraces.OtherOptions] = deal(OtherOptionsToSave);

if isfile(KeyFilePath)
    Loaded = load(KeyFilePath, 'KeyTraceData');
    KeyTraceData = [Loaded.KeyTraceData, NewKeyTraces];
else
    KeyTraceData = NewKeyTraces;
end

save(KeyFilePath, 'KeyTraceData');

fprintf('\n   Added %d traces to the key-trace collection (%s).\n', SampleSize, KeyFilePath);
fprintf('   Collection now holds %d traces total.\n', numel(KeyTraceData));

end

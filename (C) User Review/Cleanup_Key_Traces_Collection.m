function Cleanup_Key_Traces_Collection(KeyFilePath)
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
% Cleanup_Key_Traces_Collection  One-time retroactive fix for
% Key_Traces_Collection.mat files saved before Extract_Future_Key_Traces.m
% stopped attaching FindingImage to every sampled trace's OtherOptions.
%
% Every trace in the collection carries a full copy of its source
% session's OtherOptions (formerly including FindingImage, a full
% camera-resolution image -- e.g. 1024x1024 uint16, ~2MB). Since every
% trace sampled from the same session carried an IDENTICAL image copy,
% this was the dominant reason the collection file had grown far larger
% than any single review file. Extract_Future_Key_Traces.m now excludes
% FindingImage when attaching OtherOptions to newly-sampled traces; this
% script strips it from traces already saved in an existing collection
% file, so past sessions' contributions shrink too, not just future ones.
%
% Backs up the original file (append '-PreCleanupBackup' to the filename)
% before overwriting, since this touches real accumulated data -- if
% anything looks wrong after running this, the untouched original is
% still right there to restore from.
%
% KeyFilePath (optional): full path to Key_Traces_Collection.mat. Defaults
% to a uigetfile prompt if not given.

    if nargin < 1 || isempty(KeyFilePath)
        [FileName, FileDir] = uigetfile('*.mat', 'Select Key_Traces_Collection.mat');
        if isequal(FileName, 0)
            disp('Cancelled -- no file selected.');
            return
        end
        KeyFilePath = fullfile(FileDir, FileName);
    end

    if ~isfile(KeyFilePath)
        error('Cleanup_Key_Traces_Collection: file not found: %s', KeyFilePath);
    end

    fprintf('Loading %s ...\n', KeyFilePath);
    Loaded = load(KeyFilePath, 'KeyTraceData');
    KeyTraceData = Loaded.KeyTraceData;

    OriginalInfo = dir(KeyFilePath);
    fprintf('Loaded %d traces. Current file size: %.1f MB.\n', ...
        numel(KeyTraceData), OriginalInfo.bytes / 1e6);

    NumWithImage = 0;
    for k = 1:numel(KeyTraceData)
        if isfield(KeyTraceData(k), 'OtherOptions') && isfield(KeyTraceData(k).OtherOptions, 'FindingImage')
            KeyTraceData(k).OtherOptions = rmfield(KeyTraceData(k).OtherOptions, 'FindingImage');
            NumWithImage = NumWithImage + 1;
        end
    end

    if NumWithImage == 0
        fprintf('No traces carried a FindingImage field -- nothing to clean up.\n');
        return
    end

    [Dir, Name, Ext] = fileparts(KeyFilePath);
    BackupPath = fullfile(Dir, [Name '-PreCleanupBackup' Ext]);
    if isfile(BackupPath)
        error(['Cleanup_Key_Traces_Collection: backup path already exists (%s) -- ' ...
            'refusing to overwrite a prior backup. Move/rename it first if you want to re-run this.'], BackupPath);
    end
    fprintf('Backing up original to %s ...\n', BackupPath);
    copyfile(KeyFilePath, BackupPath);

    fprintf('Stripped FindingImage from %d of %d traces. Saving...\n', NumWithImage, numel(KeyTraceData));
    save(KeyFilePath, 'KeyTraceData');

    NewInfo = dir(KeyFilePath);
    fprintf('Done. File size: %.1f MB -> %.1f MB.\n', OriginalInfo.bytes / 1e6, NewInfo.bytes / 1e6);
    fprintf('Original backed up at: %s\n', BackupPath);
end

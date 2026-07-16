function newFilePath = Convert_Legacy_ProgramA_Output(oldTracesPath, outputFolder)
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
% Convert_Legacy_ProgramA_Output  Converts a legacy (pre-migration) Program A
% output file into the format current Part B (Analyze_Trace_Data.m) expects.
%
% Legacy files use the OLD frame-indexing convention: StartAnalysisFrameNumber
% == FrameNumToFindParticles, so Trace_BackSub's native index 1 IS the
% finding-image frame itself, and FocusFrameNumbers_Shifted is offset by
% (FrameNumToFindParticles-1). Current Part B assumes the NEW convention:
% StartAnalysisFrameNumber == FrameNumToFindParticles+1, so clip(1) ==
% FrameNumFound+1 and FocusFrameNumbers_Shifted is offset by FrameNumFound.
% See the Part B design notes on frame indexing for the full history.
%
% This function does NOT modify Analyze_Trace_Data.m or any Part B algorithm
% file -- it only reshapes the INPUT so the unmodified pipeline produces
% correct output (including a non-NaN OtherDataToSave.EffectiveDeadTime,
% required for Part D compatibility).
%
% Input:
%   oldTracesPath - full path to a legacy Traces/*.mat file (VirusDataToSave +
%                   OtherDataToSave, old convention).
%   outputFolder  - folder to save the converted file into. Required, not
%                   defaulted to oldTracesPath's own folder -- Run_Part_R_Driver.m
%                   deliberately keeps the legacy input (Old Traces/) and the
%                   converted output (Traces/) in separate folders, so this
%                   function never writes back into the folder it read from.
%
% Output:
%   newFilePath - path to the newly saved, current-convention-compatible file,
%                 written into outputFolder. The original file is never modified.

fprintf('Convert_Legacy_ProgramA_Output: %s\n', oldTracesPath);

Raw = load(oldTracesPath);

if isfield(Raw, 'VirusDataToSave')
    Traces = Raw.VirusDataToSave;
else
    error('Convert_Legacy_ProgramA_Output: could not find VirusDataToSave in %s.', oldTracesPath);
end

if ~isfield(Raw, 'OtherDataToSave') || ~isfield(Raw.OtherDataToSave, 'FrameNumToFindParticles')
    error('Convert_Legacy_ProgramA_Output: could not find OtherDataToSave.FrameNumToFindParticles in %s.', oldTracesPath);
end
if ~isfield(Raw.OtherDataToSave, 'VideoTimeVector')
    error('Convert_Legacy_ProgramA_Output: could not find OtherDataToSave.VideoTimeVector in %s.', oldTracesPath);
end
if ~isfield(Raw.OtherDataToSave, 'FocusFrameNumbers')
    error('Convert_Legacy_ProgramA_Output: could not find OtherDataToSave.FocusFrameNumbers in %s.', oldTracesPath);
end

fs = Raw.OtherDataToSave.FrameNumToFindParticles;
videoTimeVector = Raw.OtherDataToSave.VideoTimeVector;
focusGlobal = Raw.OtherDataToSave.FocusFrameNumbers;

% -------------------------------------------------------------------------
% VALIDATION -- confirm the OLD-convention assumption holds for THIS file
% before blindly slicing. Fail loudly rather than silently mis-convert.
% -------------------------------------------------------------------------
nTr = numel(Traces);

if isfield(Traces, 'FrameNumFound')
    frameNumFoundValues = unique([Traces.FrameNumFound]);
    if numel(frameNumFoundValues) ~= 1
        error(['Convert_Legacy_ProgramA_Output: per-trace FrameNumFound is not uniform across ' ...
            'all traces in %s (found %d distinct values). Expected a single shared value.'], ...
            oldTracesPath, numel(frameNumFoundValues));
    end
    if frameNumFoundValues(1) ~= fs
        error(['Convert_Legacy_ProgramA_Output: per-trace FrameNumFound (%d) does not match ' ...
            'file-level OtherDataToSave.FrameNumToFindParticles (%d) in %s.'], ...
            frameNumFoundValues(1), fs, oldTracesPath);
    end
else
    error('Convert_Legacy_ProgramA_Output: VirusDataToSave has no FrameNumFound field in %s.', oldTracesPath);
end

for i = 1:nTr
    tr = Traces(i);
    if numel(tr.Trace_BackSub) ~= numel(tr.TimeVector)
        error(['Convert_Legacy_ProgramA_Output: trace %d has mismatched Trace_BackSub (%d) and ' ...
            'TimeVector (%d) lengths in %s.'], i, numel(tr.Trace_BackSub), numel(tr.TimeVector), oldTracesPath);
    end
    expectedLen = numel(videoTimeVector) - fs + 1;
    if numel(tr.Trace_BackSub) ~= expectedLen
        error(['Convert_Legacy_ProgramA_Output: trace %d Trace_BackSub length (%d) does not match ' ...
            'the expected OLD-convention length (numel(VideoTimeVector)-FrameNumFound+1 = %d) in %s. ' ...
            'This file may not use the assumed legacy convention -- do not convert blindly.'], ...
            i, numel(tr.Trace_BackSub), expectedLen, oldTracesPath);
    end
end

fprintf('  Validation passed: %d traces, uniform FrameNumFound = %d, OLD-convention length holds.\n', nTr, fs);

% -------------------------------------------------------------------------
% PER-TRACE CONVERSION
% -------------------------------------------------------------------------
% New-pipeline formula: FocusFrameNumbers_Shifted = FocusFrameNumbers - fs
% (StartAnalysisFrameNumber-1 == fs under the NEW convention). Same value
% reused for every trace, matching how Part A computes it once per file.
newFocusFrameNumbersShifted = focusGlobal - fs;

for i = 1:nTr
    tr = Traces(i);

    % Drop native index 1 (the finding-image frame) so index 1 of the
    % converted trace becomes global frame fs+1, matching current Part B's
    % clip(1) == FrameNumFound+1 assumption.
    tr.Trace_BackSub = tr.Trace_BackSub(2:end);
    tr.TimeVector    = tr.TimeVector(2:end);
    if isfield(tr, 'Trace') && ~isempty(tr.Trace)
        tr.Trace = tr.Trace(2:end);
    end

    tr.FocusFrameNumbers_Shifted = newFocusFrameNumbersShifted;

    % FrameNumFound, Coordinates, Eccentricity, Area, FullFilePath,
    % StreamFilename, BoxAroundVirus, VirusIDNumber, BindingTime,
    % IgnoreFrameNumbers_Shifted, IsVirusGood, ReasonVirusFailed all pass
    % through unchanged -- Analyze_Trace_Data.m overwrites FrameNumFound
    % from the file-level value itself and never reads the per-trace copy.

    Traces(i) = tr;
end

% -------------------------------------------------------------------------
% FILE-LEVEL OtherDataToSave -- synthesize the fields current Part B needs
% -------------------------------------------------------------------------
NewOther = struct();
NewOther.GlobalFindingImageFrame  = fs;
NewOther.GlobalStartAnalysisFrame = fs + 1;   % makes EffectiveDeadTime compute for real, not NaN

% Passthrough, unchanged convention/units.
passthroughFields = {'VideoTimeVector', 'FindingImage', 'Options', 'FocusFrameNumbers', ...
    'ThresholdsUsed', 'RoughBackground', 'TotalVideoIntensity', 'AverageVideoIntensity', ...
    'StandardBindTime'};
for k = 1:numel(passthroughFields)
    fieldName = passthroughFields{k};
    if isfield(Raw.OtherDataToSave, fieldName)
        NewOther.(fieldName) = Raw.OtherDataToSave.(fieldName);
    end
end

% GlobalTimeZeroFrame and XYOffsets are deliberately NOT fabricated here --
% neither Part B (Analyze_Trace_Data.m) nor Part D reads either field
% (verified by direct code inspection), and legacy files carry no equivalent
% source data to derive them from. Inventing a placeholder value would be
% worse than omitting the field entirely.

VirusDataToSave  = Traces;
OtherDataToSave  = NewOther; %#ok<NASGU>

% -------------------------------------------------------------------------
% SAVE -- into outputFolder, clearly-marked new filename. Original file
% (still sitting wherever oldTracesPath pointed) is never touched.
% -------------------------------------------------------------------------
[~, oldName, oldExt] = fileparts(oldTracesPath);
idx = strfind(oldName, '_ExtractedTraces');
if isempty(idx)
    error(['Convert_Legacy_ProgramA_Output: expected "_ExtractedTraces" in the input filename ' ...
        '(%s) so the converted name and Part B''s own label-derivation logic stay in sync.'], oldName);
end
newName = [oldName(1:idx(1)-1) '_LegacyConverted' oldName(idx(1):end)];
if exist(outputFolder, 'dir') ~= 7
    mkdir(outputFolder);
end
newFilePath = fullfile(outputFolder, [newName oldExt]);

save(newFilePath, 'VirusDataToSave', 'OtherDataToSave');

fprintf('  Saved converted file: %s\n', newFilePath);

end

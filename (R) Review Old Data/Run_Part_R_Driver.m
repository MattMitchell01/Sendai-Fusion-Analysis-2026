function Run_Part_R_Driver()
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
% Run_Part_R_Driver  Single-video driver for Part R -- Review Old Data.
%
% Run ONCE per experiment folder (the folder labeled with the six-digit
% date + description, e.g. 250613_SomeExperimentName -- the one directly
% containing Traces/ and Analysis/). On that one run, it:
%
%   0. Restructures the experiment folder so old and new data are clearly
%      separated (see "Folder restructuring" below) -- this only runs once;
%      it errors out rather than silently repeating if the folder has
%      already been restructured.
%   1. Convert_Legacy_ProgramA_Output  -- reshape the legacy Program A
%      output you pick into current Part B's expected format, saved into
%      the new Traces/ folder
%   2. Analyze_Trace_Data              -- run the UNMODIFIED current Part B
%   3. Flag_Legacy_Disagreements       -- compare against the Old User
%      Review File's designations, flag every disagreement in place
%
% and prints an agreement-rate summary for that video. You pick both files
% yourself, one at a time; nothing is searched for automatically.
%
% Folder restructuring (step 0, happens immediately after you pick the
% experiment folder, before either file picker):
%   <ExptFolder>/Traces/    -- renamed to <ExptFolder>/Old Traces/ (raw
%                              legacy Program A output, left untouched)
%   <ExptFolder>/Traces/    -- NOT created here, deliberately -- creating it
%                              empty right away just leaves a confusing
%                              empty folder sitting next to Old Traces/
%                              while you're still picking which raw file to
%                              convert. Convert_Legacy_ProgramA_Output
%                              creates it lazily, the moment it actually has
%                              a converted file to save into it.
%   <ExptFolder>/Analysis/  -- renamed to <ExptFolder>/Old Analysis/ (so
%                              Old Analysis/AnalysisReviewed/ holds the
%                              legacy human-reviewed ground truth)
%   <ExptFolder>/Analysis/  -- also NOT created here, same reasoning;
%                              Analyze_Trace_Data creates it automatically
%                              the first time it runs, since the old one is
%                              now out of the way, and it fills with ONLY
%                              New B Output
%
% After restructuring, you pick:
%   1. The raw legacy Traces file for the video you want to reprocess --
%      now lives in <ExptFolder>/Old Traces/.
%   2. The matching Analysis/AnalysisReviewed file for the SAME video
%      number -- now lives in <ExptFolder>/Old Analysis/AnalysisReviewed/.
%      This is the Old User Review File to compare against. Click Cancel on
%      this second pick if no reviewed file exists for this video --
%      Convert+Analyze still run, just without the disagreement-flagging
%      step.
%
% (An earlier version picked the raw Traces file directly, with no
% experiment-folder restructuring step, and no root-folder batch mode
% before that. Restructuring was added so old legacy data and newly
% generated Part B/C data can never be confused for each other by folder
% name alone.)

ThisFolder = fileparts(mfilename('fullpath'));
addpath(genpath(ThisFolder));

PartBFolder = fullfile(fileparts(ThisFolder), 'Part (B)', '_Dev_B) Lipid Mixing Trace Analysis');
if exist(PartBFolder, 'dir') == 7
    addpath(genpath(PartBFolder));
else
    error('Run_Part_R_Driver: could not find Part B folder at expected path: %s', PartBFolder);
end

fprintf('====================================\n');
fprintf('Run_Part_R_Driver\n');
fprintf('====================================\n');
disp('  Run this ONCE per experiment folder. You will be prompted three times,')
disp('  in order -- a message will tell you exactly what to pick each time:')
disp('    1) The EXPERIMENT folder')
disp('    2) The RAW TRACES file to convert (from Old Traces/, after restructuring)')
disp('    3) The OLD ANALYSIS (reviewed) file to compare against (Cancel if none exists)')
disp(' ')

disp('  >> STEP 1 of 3: Now select the EXPERIMENT folder -- the one labeled with')
disp('     the six-digit date + description (e.g. 250613_SomeExperimentName)')
disp('     that directly contains Traces/ and Analysis/.')
ExptFolder = uigetdir(pwd, 'Select the EXPERIMENT folder (contains Traces/ and Analysis/ directly)');
if isequal(ExptFolder, 0)
    disp('  Cancelled -- no experiment folder selected.');
    return
end
disp(' ')

% ---- Step 0: one-time folder restructuring ----
TracesFolder      = fullfile(ExptFolder, 'Traces');
OldTracesFolder   = fullfile(ExptFolder, 'Old Traces');
AnalysisFolder    = fullfile(ExptFolder, 'Analysis');
OldAnalysisFolder = fullfile(ExptFolder, 'Old Analysis');

if exist(OldTracesFolder, 'dir') == 7 || exist(OldAnalysisFolder, 'dir') == 7
    error(['Run_Part_R_Driver: %s already contains ''Old Traces'' or ''Old Analysis'' -- ' ...
        'this experiment folder has already been restructured by Part R. Run_Part_R_Driver ' ...
        'only restructures a given experiment folder once.'], ExptFolder);
end
if exist(TracesFolder, 'dir') ~= 7
    error('Run_Part_R_Driver: no Traces/ folder found in %s.', ExptFolder);
end

fprintf('  Restructuring %s:\n', ExptFolder);
movefile(TracesFolder, OldTracesFolder);
fprintf('    Traces/ -> Old Traces/  (a new Traces/ will be created once the converted file is saved)\n');

if exist(AnalysisFolder, 'dir') == 7
    movefile(AnalysisFolder, OldAnalysisFolder);
    fprintf('    Analysis/ -> Old Analysis/  (Analyze_Trace_Data will create a fresh Analysis/)\n');
else
    fprintf('    No existing Analysis/ folder found -- nothing to rename.\n');
end
disp(' ')

% ---- Step 1: pick the raw Traces file (now under Old Traces/) ----
disp('  >> STEP 2 of 3: Now select the RAW TRACES file for this video that you')
disp('     want to convert -- it now lives in the Old Traces/ folder just')
disp('     created above. (There is no Traces/ folder yet -- it is not')
disp('     created until the converted file is saved into it.)')
[TracesFile, TracesFilePicked] = uigetfile(fullfile(OldTracesFolder, '*.mat'), ...
    'Step 1: Select the raw legacy Traces file for this video (in Old Traces/)');
if isequal(TracesFile, 0)
    disp('  Cancelled -- no Traces file selected. (The folder restructuring above already happened;');
    disp('  re-run and pick a file to continue, no need to redo the restructuring.)');
    return
end
TracesPath = fullfile(TracesFilePicked, TracesFile);
VideoLabel = Extract_Leading_Number(TracesFile);
fprintf('  Traces file: %s (video %s)\n\n', TracesPath, VideoLabel);

% ---- Step 2: pick the matching Old User Review File (now under Old Analysis/AnalysisReviewed/) ----
OldReviewedStartFolder = fullfile(OldAnalysisFolder, 'AnalysisReviewed');
if exist(OldReviewedStartFolder, 'dir') ~= 7
    OldReviewedStartFolder = OldAnalysisFolder;
end
disp(' ')
disp('  >> STEP 3 of 3: Now select the OLD ANALYSIS (reviewed) file for the SAME')
disp('     video to compare against -- it lives in Old Analysis/AnalysisReviewed/.')
disp('     Click Cancel if no reviewed file exists for this video.')
[ReviewedFile, ReviewedFilePicked] = uigetfile(fullfile(OldReviewedStartFolder, '*.mat'), ...
    'Step 2: Select the matching Old User Review File (Cancel if none exists)');
if isequal(ReviewedFile, 0)
    disp('  No Old User Review File selected -- Flag_Legacy_Disagreements will be skipped.');
    disp(' ');
    ReviewedPath = '';
else
    ReviewedPath = fullfile(ReviewedFilePicked, ReviewedFile);
    ReviewedLabel = Extract_Leading_Number(ReviewedFile);
    if ~strcmp(ReviewedLabel, VideoLabel)
        warning('Run_Part_R_Driver:VideoNumberMismatch', ...
            ['Selected reviewed file''s leading video number (%s) does not match the Traces ' ...
             'file''s (%s) -- double check you picked the matching file.'], ReviewedLabel, VideoLabel);
    end
    fprintf('  Old User Review File: %s\n\n', ReviewedPath);
end

try
    ConvertedPath = Convert_Legacy_ProgramA_Output(TracesPath, TracesFolder);
    BOutputPath   = Analyze_Trace_Data(ConvertedPath);
catch ME
    error('Run_Part_R_Driver: error processing video %s: %s', VideoLabel, ME.message);
end

% Flag_Legacy_Disagreements already prints its own comparison summary
% (Matched/Agreements/Disagreements/Unmatched) and the saved file path --
% not repeated here, so results are reported exactly once.
if isempty(ReviewedPath)
    fprintf(['  No reviewed file was selected -- skipping Flag_Legacy_Disagreements. This New B ' ...
        'Output has no .IsLegacyDisagreement field; review it in Part C with ' ...
        'Options.ReviewOldDataOnlyDifferences=''n'', not ''y''.\n\n']);
    fprintf('  New B Output saved: %s\n', BOutputPath);
else
    Flag_Legacy_Disagreements(BOutputPath, ReviewedPath);
end

fprintf('====================================\n');
fprintf('Run_Part_R_Driver complete -- video %s\n', VideoLabel);
fprintf('====================================\n');
fprintf('  Ready for Part C review from: %s\n', AnalysisFolder);

end

% =============================================================================
function NumStr = Extract_Leading_Number(FileName)
% Extract_Leading_Number  Leading digit run of a filename, used only as a
% sanity-check that the two manually-picked files (Traces + Reviewed) are
% for the same video number -- not used to search for or match files.

    Tokens = regexp(FileName, '^(\d+)', 'tokens', 'once');
    if isempty(Tokens)
        NumStr = '';
    else
        NumStr = Tokens{1};
    end
end

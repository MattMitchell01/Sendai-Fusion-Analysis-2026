%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
% An_Easy_Start_UserReview
%
% Entry point for the new Part C -- User Review Traces.
% BEFORE RUNNING: check Setup_Options_User_Review.m.
%
% Select a Part B output file (or a previously-reviewed file, to resume)
% and the resolved path is handed straight to Start_User_Review.

disp('====================================');
disp('   Part C  —  User Review Traces');

% ── Data location: comment one line, uncomment the other ─.──────────────────
% Matt recommends that you manually set the DataLocation to the Master Analysis
% folder on whatever computer you are using that way grabbing the individual
% data files is easier.
%
% IMPORTANT: this is no longer just a file-picker convenience default -- it
% is also the shared MASTER ANALYSIS FOLDER where Key_Traces_Collection.mat
% and Performance_Log.mat live and accumulate across EVERY dataset/experiment
% you ever review, not just this one. Make sure this points to one stable,
% correct, shared location (not a per-experiment subfolder) -- pointing it
% somewhere different between sessions will fragment your key-trace
% collection and performance log across multiple files instead of building
% up one running record. See Extract_Future_Key_Traces.m / Update_Performance_Log.m.

%DataLocation_UserReview = uigetdir([], 'Choose parent data folder');
DataLocation_UserReview = '/Users/littlem/Research/analysis';
% ───────────────────────────────────────────────────────────────────────────

disp('   Choose your Part B output file (or a previously-reviewed file, to resume)...');
[inputFile, inputDir] = uigetfile('*.mat', 'Select Part B Output File', DataLocation_UserReview);
if isequal(inputFile, 0); disp('   Cancelled.'); return; end
SelectedFilePath = fullfile(inputDir, inputFile);

disp('   Starting review...');
Start_User_Review(SelectedFilePath, DataLocation_UserReview);

%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
% An_Easy_Start_TraceAnalysis
%
% Entry point for the new Part B pipeline.
% BEFORE RUNNING: check Parameters_Seek_Fusion.m and Parameters_Seek_Unbounds.m.
%
% Output is saved automatically next to your data:
%   ExptFolder/Traces/103_ExtractedTraces.mat  →  ExptFolder/Analysis/103-AnalyzedTraces.mat

disp('====================================');
disp('   Part B  —  Lipid Mixing Trace Analysis');

% ── Data location: comment one line, uncomment the other ───────────────────
%DataLocation_TraceAnalysis = uigetdir([], 'Choose parent data folder');
DataLocation_TraceAnalysis = '/Users/littlem/Research/analysis';
% ───────────────────────────────────────────────────────────────────────────

disp('   Choose your Program A output file (.mat)...');
[inputFile, inputDir] = uigetfile('*.mat', 'Select Program A Output File', DataLocation_TraceAnalysis);
if isequal(inputFile, 0); disp('   Cancelled.'); return; end
dataFilePath = fullfile(inputDir, inputFile);

disp('   Running analysis...');
Analyze_Trace_Data(dataFilePath);

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
% An_Easy_Start_ReviewByWaitTime
%
% Entry point for targeted re-review by wait-time window -- for chasing a
% jump/artifact spotted in Part D's CDF back to the specific traces
% causing it, without re-running the full Low/Medium/High Part C queue.
%
% Select an EXISTING Part C User Review File (<name>-Revd.mat, produced by
% Start_User_Review.m) -- this tool never bootstraps a fresh one from raw
% Part B output. Corrections you make here are saved back into that same
% file, so a subsequent Part D run picks them up immediately.

disp('====================================');
disp('   Review By Wait Time  —  targeted re-review for a CDF jump');

disp('   Choose the User Review File to re-review (e.g. *-Revd.mat)...');
[inputFile, inputDir] = uigetfile('*.mat', 'Select User Review File (-Revd.mat)');
if isequal(inputFile, 0); disp('   Cancelled.'); return; end
SelectedFilePath = fullfile(inputDir, inputFile);

disp('   Starting wait-time review...');
Start_Review_By_Wait_Time(SelectedFilePath);

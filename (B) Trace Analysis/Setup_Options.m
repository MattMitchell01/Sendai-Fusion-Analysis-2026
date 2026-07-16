function Options = Setup_Options()
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
% Setup_Options  Overall, pipeline-level options for Analyze_Trace_Data.
%                Assembles the per-algorithm parameter structs (from
%                Parameters_Seek_Fusion.m, Parameters_Seek_Unbounds.m,
%                Parameters_Assign_Review_Priority.m) alongside settings
%                that apply across the whole pipeline rather than to any
%                single algorithm.
%


% =========================================================================
% DISPLAY
%
% Single pipeline-wide toggle. Analyze_Trace_Data.m applies this to every
% sub-algorithm's own Display.DoPlot after loading it (replacing whatever
% each algorithm's own Parameters_*.m file defaults to), so there is one
% place to turn diagnostic plotting on/off for an entire run.
% =========================================================================

Options.Display.DoPlot = false;


% =========================================================================
% CLIP WIDTH
%
% Frames dropped from each end of every trace by Seek_Fusion/Seek_Unbounds
% after smoothing (their real values are still used as smoothing context --
% see Calculate_Sliding_Window). A single shared value for BOTH detectors --
% lives here, not duplicated in either algorithm's own parameters, so both
% detectors and the +ClipWidth global-frame conversion in
% Analyze_Trace_Data.m always agree on the same dead-zone boundary.
% =========================================================================

Options.ClipWidth = 5;


% =========================================================================
% PER-ALGORITHM PARAMETERS
% =========================================================================

Options.SeekFusion           = Parameters_Seek_Fusion();
Options.SeekUnbounds         = Parameters_Seek_Unbounds();
Options.AssignReviewPriority = Parameters_Assign_Review_Priority();

end

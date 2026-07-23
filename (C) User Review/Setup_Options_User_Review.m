function [Options] = Setup_Options_User_Review()
%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
    
    % ------Easy reference guide for prompt codes-----
        % To enter at prompt: PlotNumber.DesignationCode
        % See "Reference Guide For Prompt Codes.rtf" for the full table.
        % DesignationCode as follows:
        % .0   = No Fusion
        % .1   = 1 Fuse
        % .11  = 1 Fuse designation is already correct, but fuse frame/wait time is wrong
        % .2   = 2 Fuse (click twice: once per fusion time)
        % .22  = 2 Fuse designation is already correct, but fuse frames/wait times are wrong
        % .3   = Unbound event (see sharp decrease; click the unbind frame)
        % .33  = Unbound designation is already correct, but unbind frame/wait time is wrong
        % .9   = Other / hard to classify

    % ------Options to double check-------
        Options.Label = '-Revd';
        Options.ApplyFilter = 'n';
            % Currently unused -- the reference tool's dataset-specific
            % yellow-trace filter branch (1 Fuse plotted yellow if
            % FuseFrameNumbers(1) > 960) was dropped when Plot_Current_Trace.m
            % was adapted to the new Part B schema. Left here in case a
            % similar per-dataset filter is needed again later.

    % -------Options you are less likely to change regularly-------
        Options.ReviewFolderName = 'User_Revd_Analysis';
            % User Review File folder name. Treat as effectively constant --
            % changing it between sessions on a dataset that already has a
            % User Review File will cause Start_User_Review to no longer
            % recognize that file (see Get_User_Review_File_Path.m).

        Options.Low    = struct('Rows', 4, 'Cols', 5);
        Options.Medium = struct('Rows', 3, 'Cols', 4);
        Options.High   = struct('Rows', 3, 'Cols', 4);
            % Shared grid for every High subgroup EXCEPT H1 (H2a, H2b, H3,
            % H4, H5 all use this one uniform setting -- not a per-subgroup
            % override each, on purpose, since none of them need a layout
            % different from the others). Smaller than Low/Medium's grids
            % since H2a/H2b's traces also carry a full focus-dot overlay
            % (see Plot_Current_Trace.m's DiagnosticOverlay.ShowFocusDots),
            % a busier plot than the plain single-line view -- a starting
            % point, tune Rows/Cols here if it needs adjusting.
            % Baked into ReviewQueue.Segments at first-run preprocessing time --
            % editing these later only affects newly-preprocessed datasets.

        Options.HighOverrides.H1 = struct('Rows', 2, 'Cols', 4, 'TracesPerBatch', 4);
            % H1 is the one High subgroup that needs its own custom
            % layout instead of the shared Options.High: Rows is fixed at
            % 2 (row 1 = full-trace overview, row 2 = the same trace's
            % zoomed first/last-100-frame view, directly below it in the
            % same column) -- Cols is the tunable "traces per round" knob,
            % starting at 5. TracesPerBatch (5) tells the review loop each
            % round consumes Cols traces, not Rows*Cols (2 physical axes
            % per trace here, not 1) -- see Build_Review_Queue.m and
            % Plot_Trace_Pair_H1.m.

        Options.QuickModeNoCorrection = 'n';
    
        Options.FixWaitTime = 'y'; %allow user to manually fix wait time when they choose 1 Fuse, 2 Fuse or Unbound.

        Options.UseRunMed = 'y'; % show running median instead of raw trace
        Options.RunMedHalfLength = 1; % num of data points on either side to include in running median (window of 3: 1 each side)

        Options.ShowBindFrame = 'n'; %Will draw a black dashed line
        Options.ExpandXAxis = 'y'; %will expand the x-axis so the beginning/end of the trace aren't obscured by y-axis

        Options.ShowFocusDotsAlways = 'y';
            % Global override: draw the red/black focus-event dot
            % overlay (Plot_Focus_Dot_Markers.m) on every trace in every
            % segment -- Low, Medium, and all six High subgroups -- not
            % just H2a/H2b/H3/H4/H5, which already always show it
            % regardless of this flag. Set to 'n' to fall back to the old
            % per-subgroup behavior (only H2a/H2b/H3/H4/H5 get dots).

        Options.ApplyH2FocusJumpCorrection = 'y';
            % H2a/H2b only: subtract the focus-induced jump near the
            % candidate/actual fusion frame from every index onward, so
            % the reviewer sees the trace as if that jump hadn't happened
            % (Get_H2_Focus_Jump_Trigger_Frame.m, applied inside
            % Plot_Current_Trace.m). Display-only -- never touches saved
            % data. Requires Part B's OtherDataToSave.ClipWidth passthrough
            % field; quietly does nothing on a dataset that predates it.
            % Set to 'n' to fall back to the plain focus-dot overlay only.

        Options.GrabExampleTrace = 'n'; %Use to grab example trace.
            % The rest of user review will be ignored. It will grab the trace of the starting trace number.

        Options.ReportPerformanceLog = 'y';
            % At the end of a completed review session, update the cumulative
            % TP/FN/FP performance log (Performance_Log.mat, in the master
            % analysis folder -- see An_Easy_Start_UserReview.m's
            % DataLocation_UserReview) comparing this session's AlgoDesignation
            % vs. final Designation. See Update_Performance_Log.m.

        Options.ExtractKeyTraces = 'y';
            % At the end of a completed review session, harvest a random 20%
            % of this session's traces (excluding anything designated 'Other'
            % by either the algorithm or the reviewer) into the shared
            % cross-experiment Key_Traces_Collection.mat, in the same master
            % analysis folder. See Extract_Future_Key_Traces.m.

        Options.ReviewOldDataMode = 'y';
            % 'y' when reviewing a legacy dataset reprocessed through Part R
            % (Convert_Legacy_ProgramA_Output.m / Flag_Legacy_Disagreements.m).
            % Auto-disables ReportPerformanceLog/ExtractKeyTraces above (see
            % Start_User_Review.m) so a legacy-dataset review session never
            % writes into the shared cross-dataset master tracking files used
            % for live-pipeline QA.

        Options.ReviewOldDataOnlyDifferences = 'y';
            % Only meaningful when ReviewOldDataMode='y'. 'y' filters the
            % review queue down to traces with .IsLegacyDisagreement==true
            % (set by Flag_Legacy_Disagreements.m) before Low/Medium/High/
            % H1-H5 bucketing -- excluded traces are still preserved in the
            % saved file, just routed to a trailing segment the review loop
            % never walks. See Build_Review_Queue.m.

end

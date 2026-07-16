function [FigureHandles] = Plot_Trace_Pair_H1(FigureHandles, CurrentVirusData, VideoTimeVector, ClipWidth, ...
    PlotCounter, CurrentTraceNumber, Options, Segment, DiagnosticOverlay, ManualFocusSubtractIndices)
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
% Plot_Trace_Pair_H1  H1's two-row trace display: row 1 is the full-trace
% overview for every trace in the batch, row 2 is the same traces' zoomed
% first/last-100-frame view, directly below its own column (whichever edge
% triggered the H1 flag). Reuses Plot_Current_Trace unchanged for all
% designation-coloring/marker logic in both panes -- only the axes each
% call targets, and the zoom pane's xlim/title, differ.
%
% Rows is fixed at 2 by this design (row 1 = full, row 2 = zoom); Cols is
% the tunable "traces per round" knob (Segment.TracesPerBatch == Cols).
% PlotCounter is this trace's 1-based ordinal position within the current
% batch (1..Segment.Cols) -- same meaning as everywhere else, and what the
% reviewer types back at the PlotNumber.Code prompt. tight_subplot fills
% row-major, so row 1 occupies slots 1..Cols and row 2 occupies slots
% Cols+1..2*Cols.
%
% No extra marker is drawn for Review_PriorityData.TriggerFrame -- it's
% the same fusion/unbound frame Plot_Current_Trace already marks in green,
% so once the zoom pane is cropped to it, that existing marker is enough.
%
% VideoTimeVector and ClipWidth are pure pass-throughs to both panes'
% Plot_Current_Trace calls (see Plot_Current_Trace.m for what each is
% used for there).
%
% DiagnosticOverlay (optional, default []) is forwarded unchanged to both
% panes' Plot_Current_Trace calls -- currently only used to carry
% Options.ShowFocusDotsAlways's global focus-dot overlay onto H1 traces,
% same as every other segment.
%
% ManualFocusSubtractIndices (optional, default []) is likewise forwarded
% unchanged to both panes -- see Plot_Current_Trace.m for what it does.

    if nargin < 9
        DiagnosticOverlay = [];
    end
    if nargin < 10
        ManualFocusSubtractIndices = [];
    end

    assert(Segment.Rows == 2, 'Plot_Trace_Pair_H1:UnsupportedGrid', ...
        'Plot_Trace_Pair_H1 assumes a fixed 2-row layout (row 1 = full, row 2 = zoom) -- got Rows=%d.', ...
        Segment.Rows);

    FullSlot = PlotCounter;
    ZoomSlot = Segment.Cols + PlotCounter;

    FullAxes = FigureHandles.SubHandles(FullSlot);
    ZoomAxes = FigureHandles.SubHandles(ZoomSlot);

    % Full pane: Plot_Current_Trace against a synthetic single-axes
    % FigureHandles -- PlotCounter passes through uncapped for a correct
    % title (only the axes lookup is clamped to the lone handle).
    FullFH.MasterWindow = FigureHandles.MasterWindow;
    FullFH.SubHandles = FullAxes;
    Plot_Current_Trace(FullFH, CurrentVirusData, VideoTimeVector, ClipWidth, PlotCounter, CurrentTraceNumber, Options, DiagnosticOverlay, ManualFocusSubtractIndices);

    % Zoom pane: same call, then crop to the first/last 100 clipped-
    % coordinate frames per EdgeLocation and relabel the title.
    ZoomFH.MasterWindow = FigureHandles.MasterWindow;
    ZoomFH.SubHandles = ZoomAxes;
    Plot_Current_Trace(ZoomFH, CurrentVirusData, VideoTimeVector, ClipWidth, PlotCounter, CurrentTraceNumber, Options, DiagnosticOverlay, ManualFocusSubtractIndices);

    [ClippedTrace, ~] = Clip_Trace_For_Review(CurrentVirusData);
    ClipLen = numel(ClippedTrace);
    EdgeLocation = CurrentVirusData.Review_PriorityData.EdgeLocation;
    if strcmp(EdgeLocation, 'start')
        ZoomRange = [1, min(100, ClipLen)];
        ZoomLabel = 'first 100 frames';
    else
        ZoomRange = [max(1, ClipLen - 99), ClipLen];
        ZoomLabel = 'last 100 frames';
    end
    xlim(ZoomAxes, ZoomRange);
    title(ZoomAxes, sprintf('%d - Zoom: %s', PlotCounter, ZoomLabel));
end

function [FigureHandles] = Create_Master_Window(Rows, Cols)
%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
% Create_Master_Window  Build the master grid-of-subplots window.
%
% Tier-agnostic: Rows/Cols come from the caller's ReviewQueue.Segments
% row (Options.Low/Medium/High), not from any notion of tier baked in
% here -- Low, Medium, and High all use this same function.

    FigureHandles.MasterWindow = figure(1);
    clf(FigureHandles.MasterWindow);   % remove any leftover axes from a
        % previously-drawn segment with a different Rows/Cols grid shape --
        % tight_subplot only ADDS new axes, it doesn't clear old ones, so
        % without this, switching between segments (forward or via the
        % cross-segment 'b' go-back) leaves the old grid's axes/titles
        % overlapping the new one.
    set(FigureHandles.MasterWindow, 'WindowStyle', 'normal');
    set(0, 'DefaultAxesFontSize',11)

    Gap = [.04,.01];
    MarginsHeight = [.04,.04];
    MarginsWidth = [.03,.02];

    FigureHandles.SubHandles = tight_subplot(Rows, Cols, Gap, MarginsHeight, MarginsWidth);

end

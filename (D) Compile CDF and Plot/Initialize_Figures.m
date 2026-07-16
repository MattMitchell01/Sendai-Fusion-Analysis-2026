function [FigureHandles] = Initialize_Figures(Options)

%set(0, 'DefaultAxesFontSize',20)

    FigureHandles.CDFWindow=figure(1);
    FigureHandles.EffUnbWindow = figure(2);

    set(FigureHandles.EffUnbWindow, 'Position', [745.8000 341.8000 560 420]);
    set(FigureHandles.CDFWindow, 'Position', [184.2000 341.8000 560 420]);

end
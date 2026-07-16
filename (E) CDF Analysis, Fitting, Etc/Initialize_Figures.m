function [FigureHandles] = Initialize_Figures(Options)

% set(0, 'DefaultAxesFontSize',20)

    CDFFigSize = [556.8000 370.4000];

    FigureHandles.CDFWindow=figure(1);
    FigureHandles.EfficiencyWindow = figure(4);

        set(FigureHandles.CDFWindow, 'Position', [0 407.4000 CDFFigSize(1) CDFFigSize(2)]);
        set(FigureHandles.EfficiencyWindow, 'Position', [0 17.8000 CDFFigSize(1) 307.2000]);

    if strcmp(Options.SimpleComparisonOnly,'n')
        
        % if strcmp(Options.CalcRandomParam,'y')
        %     FigureHandles.NMinWindow = figure(6);
        %     set(FigureHandles.NMinWindow,'Position',[-5 50 552.8000 307.2000])
        % end

        if strcmp(Options.RunFit_SingleFitMultipleCDFs,'y') || strcmp(Options.RunFit_MultipleFitsSingleCDF,'y') 
            FigureHandles.FitWindow=figure(2);
            set(FigureHandles.FitWindow,'Position',[CDFFigSize(1) 407.4000 CDFFigSize(1) CDFFigSize(2)])
            
            FigureHandles.ResidualsWindow=figure(3);
            set(FigureHandles.ResidualsWindow,'Position',[CDFFigSize(1) 17.8000 CDFFigSize(1) 307.2000])
        end

    %       OLD STUFF BELOW, STILL NEED TO DEAL WITH
    %     FigureHandles.ResidualsWindow=figure(3);
    %     FigureHandles.FitWindow=figure(2);
    % 
    %     FigureHandles.NotNormalizedFusionWindow = figure(7);
    % 
    %     if strcmp(Options.ShowBeginningIntensity,'y')
    %         FigureHandles.IntensityWindow = figure(5);
    %         set(FigureHandles.IntensityWindow,'Position',[500 50 400 400])
    %     end
    % 
    %     if strcmp(Options.RunBootstrap,'y')
    %         FigureHandles.BootstrapWindow = figure(6);
    %         set(FigureHandles.BootstrapWindow,'Position',[-5 50 560 420])
    %     end
    % 
    % %     Fuse2Wind=figure(4);
    % %     FigureHandles.HistogramWindow = figure(4);
    % 
    % %     set(FigureHandles.FitWindow, 'Position', [562   366   560   420]);
    % %     set(FigureHandles.ResidualsWindow, 'Position', [1124  368  560  420]);
    % %     set(FigureHandles.Fuse1Wind, 'Position', [-5   364   560   420]);
    % %     
    %     set(FigureHandles.FitWindow, 'Position', [-5   364   560   420]);
    %     set(FigureHandles.ResidualsWindow, 'Position', [562   366   560   420]);
    % 
    %     set(FigureHandles.NotNormalizedFusionWindow, 'Position', [1343          99         560         420]);
    end

end
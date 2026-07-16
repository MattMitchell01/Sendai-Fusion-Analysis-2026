function [FigureHandles] = Setup_Figures(Options)
    
        if isfield(Options,'QuickScanThreshold')
            if strcmp(Options.QuickScanThreshold,'y')
                    FigureHandles.MasterWindow = figure(1);
                    set(FigureHandles.MasterWindow, 'Position', [2 53 1278 652]);
                    set(0, 'DefaultAxesFontSize',11)

                    NumPlotsX = 3;
                    NumPlotsY = 2;

                    Gap = [.04,.01];
                    MarginsHeight = [.04,.04];
                    MarginsWidth = [.03,.02];

                    FigureHandles.SubHandles = tight_subplot(NumPlotsY, NumPlotsX, Gap, MarginsHeight, MarginsWidth);                
                return
            end
        end

        FigureHandles.BinaryImageWindow = figure(3);
        set(FigureHandles.BinaryImageWindow,'Position',[1 -50 450 341]);
        FigureHandles.BackgroundTraceWindow = figure(4);
        set(FigureHandles.BackgroundTraceWindow,'Position',[452 -130 450 341]);
        FigureHandles.ImageWindow = figure(1);
        set(FigureHandles.ImageWindow,'Position',[1   479   451   338]);
        FigureHandles.CurrentTraceWindow = figure(2);
        set(FigureHandles.CurrentTraceWindow,'Position',[472   476   450   341]);
        
        if Options.UserGrabExampleTrace == 'y'
            FigureHandles.UserExampleTrace = figure(6);
            set(FigureHandles.UserExampleTrace,'Position',[250 -50 700 550]);
        end
        
end
function []= Display_All_Particles(VirusDataToSave,FigureHandles)

NumberParticles = length(VirusDataToSave);
set(0,'CurrentFigure',FigureHandles.ImageWindow);
                
for n = 1:NumberParticles
    CurrentParticleBox = VirusDataToSave(n).BoxAroundVirus;
    IsParticleGood = VirusDataToSave(n).IsVirusGood;

    CVB = CurrentParticleBox;
    BoxToPlot = [CVB.Bottom,CVB.Left;CVB.Bottom,CVB.Right;CVB.Top,CVB.Right;CVB.Top,CVB.Left;CVB.Bottom,CVB.Left];
    
    if strcmp(IsParticleGood,'y')
        LineColor = 'g-';

    elseif strcmp(IsParticleGood,'n')
        LineColor = 'r-';
        
    elseif strcmp(IsParticleGood,'Ignore')
        LineColor = 'y-';
    end
            
    plot(BoxToPlot(:,2),BoxToPlot(:,1),LineColor)
    hold on            
end

    drawnow
end
function [] = User_Grab_Example_Trace(FigureHandles,VirusDataToSave,ImageStackMatrix,RoughBack_Med,NumFrames)

% I haven't annotated this much at all. If you decide to use this option,
% you should modify it according to what you need, as it is currently set 
% up to capture a specific trace from the last set of data that I was 
% looking at. However, I decided to keep this script here as a starting point for you.

dbstop in User_Grab_Example_Trace at 94
VirusNumbertoKeep = 1;

CropWindow = FigureHandles.ImageWindow;

    % Define the ROI of the pH indicator
    %Choose the region manually
%     [~, Coordinates] = imcrop(CropWindow);
%         X = round(Coordinates(1));
%         Y = round(Coordinates(2));
%         Width = round(Coordinates(3));
%         Height = round(Coordinates(4));

        X = 14;
        Y = 50;
        Width = 16;
        Height = 14;
%         
        CurrentVirusBox.Left = X;
        CurrentVirusBox.Right = X+ Width;
        CurrentVirusBox.Top = Y +Height;
        CurrentVirusBox.Bottom = Y;
            
        CVB = CurrentVirusBox;
        BoxToPlot = [CVB.Bottom,CVB.Left;CVB.Bottom,CVB.Right;CVB.Top,CVB.Right;CVB.Top,CVB.Left;CVB.Bottom,CVB.Left];

        set(0,'CurrentFigure',FigureHandles.ImageWindow);
                plot(BoxToPlot(:,2),BoxToPlot(:,1),'g-')
                hold on
                drawnow
            

    CurrentTraceArray = ImageStackMatrix(...
                    CurrentVirusBox.Bottom:CurrentVirusBox.Top,...
                    CurrentVirusBox.Left:CurrentVirusBox.Right,...
                    1:NumFrames);

                CurrentTraceSumArray = sum(sum((CurrentTraceArray)));

            %Because CurrentTraceSumArray is a 3D array, the summed data is
            %transferred to row vector (CurrentTraceIntegrated), so it can be
            %plotted.
                CurrentTraceIntegrated = shiftdim(CurrentTraceSumArray(1,1,:),1);
                %CurrentTraceIntegrated_SimpleBackSub = CurrentTraceIntegrated -...
                %    RoughBackground(1:NumFrames)'.*((Options.ROI_Radius*2)+1)^2; 

                CurrentTraceIntegrated_SimpleBackSub = CurrentTraceIntegrated -...
                    RoughBack_Med(1:NumFrames)'.*((Width+1)*(Height+1)); 
                

%  PHData = CurrentTraceIntegrated_SimpleBackSub/max(CurrentTraceIntegrated_SimpleBackSub);
 PHData = CurrentTraceIntegrated./((Width+1)*(Height+1));
 
 PHData_Zero = PHData - min(PHData);
 PHData = PHData_Zero;
 
 BoxAroundVirus = VirusDataToSave(VirusNumbertoKeep).BoxAroundVirus;
 WidthVirus = BoxAroundVirus.Right - BoxAroundVirus.Left;
 HeightVirus = BoxAroundVirus.Bottom - BoxAroundVirus.Top;
 %VirusData = VirusDataToSave(VirusNumbertoKeep).Trace_BackSub./((WidthVirus+1)*(HeightVirus+1));
 VirusData = VirusDataToSave(VirusNumbertoKeep).Trace_BackSub;
 
 Time =(1: length(VirusData))*.288 -.288;
 
set(0,'CurrentFigure',FigureHandles.UserExampleTrace)
cla
%plot(CurrentTraceIntegrated);
%hold on
[Axes,VirusHandle,PHHandle] =plotyy(Time,VirusData,Time,PHData);
% PHLimit = 75;
% [Axes,VirusHandle,PHHandle] =plotyy(Time,VirusData,Time(1:PHLimit),PHData(1:PHLimit));
xlabel(Axes(1),'Time (s)')
ylabel(Axes(1),'Virus Intensity (AU)')
ylabel(Axes(2),'PH Indicator Intensity (AU)')
% ylim(Axes(2),[3.5*10^4 5*10^4])
% Axes(2).YTickLabel  = {'0.9','1','1.1','1.2','1.3','1.4'} ;
VirusHandle.Color = 'r';
PHHandle.Color = 'g';
Axes(1).YColor = 'k';
Axes(2).YColor = 'k';
xlim(Axes(2),[0 150])
xlim(Axes(1),[0 150])
ylim(Axes(2),[0 100])
ylim(Axes(1),[-2000 10000])
%ylim(Axes(1),[-20 80])
Axes(2).YTick = [0:20:100];
hold on 
                
end
function [Offset] = Single_Particle_Tracking(VideoMatrix,BWVideoMatrix,...
            FigureHandles,Options)
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
    
    % Tell the user what is going on
    disp('   You selected "Options.DeterminePixelOffset = y"')
    disp('   So, please select a bright, isolated spot in Figure 1.')
    disp('     I will track it over time to determine the offset.')
    disp('     Pro tip: dont select a spot too close to the edge or you will get an error')
    disp('     Pro tip: if clicking does nothing, click once on the MATLAB app/window')
    disp('     itself first (e.g. its Dock icon) to give it focus, then click your spot.')

    % NOTE: at this point Figure 1 is still showing FindingImage (the same
    % image particle boxes were located on), so XStart/YStart below are in
    % FindingImage's coordinate space. Offset.x/y are computed relative to
    % that click, NOT relative to VideoMatrix frame 1 -- see the note at
    % the bottom of this function. If the particle drifted more than
    % Options.SearchRadius between the finding frame and VideoMatrix frame
    % 1 (the start frame), the very first search below can still fail to
    % locate it and hit ErrorCantFindVirus -- that is a pre-existing
    % limitation of the search-radius approach, not something this
    % reference-frame fix addresses.


    %Manually choose the spot you want to analyze
    % Use figure(...) rather than set(0,'CurrentFigure',...) here --
    % the latter only changes which figure MATLAB considers "current"
    % internally, it does NOT raise/focus the window on screen. Since a
    % different figure (the FFT/SPT shift-plot window) is typically the
    % most recently focused window at this point, skipping the actual
    % raise/focus means ginput's crosshair can render on hover but real
    % clicks won't register as data points.
    figure(FigureHandles.ImageWindow);
    drawnow;
    [Pos.x,Pos.y] = ginput(1);

    disp('   Thanks. Tracking particle now.')
    
    Search_Radius = Options.SearchRadius;

    XStart = Pos.x; YStart = Pos.y;

    XOld = XStart; YOld = YStart;

    NumberFramesToAnalyze = size(VideoMatrix,3);
%     if isnan(Options.FrameNumberLimit)
%         NumberFramesToAnalyze = size(VideoMatrix,3);
%     else
%         NumberFramesToAnalyze = Options.FrameNumberLimit - (Options.StartAnalysisFrameNumber - 1);
%     end

    %Track the spot over time
    for CurrFrame = 1:NumberFramesToAnalyze

        ImAreaToSearch = VideoMatrix(...
                        round(YOld)-Search_Radius:round(YOld)+Search_Radius,...
                        round(XOld)-Search_Radius:round(XOld)+Search_Radius,...
                        CurrFrame);
        BWAreaToSearch = BWVideoMatrix(...
                        round(YOld)-Search_Radius:round(YOld)+Search_Radius,...
                        round(XOld)-Search_Radius:round(XOld)+Search_Radius,...
                        CurrFrame);

        VirusesFound = bwconncomp(BWAreaToSearch,8);
        NumVirusesFound = VirusesFound.NumObjects;

        if (NumVirusesFound == 0)
            %Means that it moves off the viewing area (or fully fuses)
            disp('   ---------Uh oh!!---------')
            disp('   Error-cant find virus. Maybe moved off viewing area or fused?')

            ErrorCantFindVirus
        end

        PropsOfVirusesFound = regionprops(VirusesFound, ImAreaToSearch, 'WeightedCentroid');

        if (NumVirusesFound == 1)
            XNew = PropsOfVirusesFound(1).WeightedCentroid(1)-Search_Radius-1+XOld;
            YNew = PropsOfVirusesFound(1).WeightedCentroid(2)-Search_Radius-1+YOld;
        elseif (NumVirusesFound > 1) %If there is more than one Virus found, use the one that is closest to the last measurement.
            DistToPrevCent = zeros(1,NumVirusesFound);
            for i = 1:NumVirusesFound
                XTest = PropsOfVirusesFound(i).WeightedCentroid(1)-Search_Radius-1+XOld;
                YTest = PropsOfVirusesFound(i).WeightedCentroid(2)-Search_Radius-1+YOld;
                DistToPrevCent(i) = sqrt((YTest-YOld)^2+(XTest-XOld)^2);
            end

            IdxVirusToUse = find(DistToPrevCent==min(DistToPrevCent));
            XNew = PropsOfVirusesFound(IdxVirusToUse).WeightedCentroid(1)-Search_Radius-1+XOld;
            YNew = PropsOfVirusesFound(IdxVirusToUse).WeightedCentroid(2)-Search_Radius-1+YOld;

        end

        XOld = XNew; YOld = YNew;
        Pos.y(CurrFrame) = YNew;
        Pos.x(CurrFrame) = XNew;

        
        %Plot out the tracking for diagnostic purposes
            figure(33)
            ImToShow = VideoMatrix(...
                           round(YStart)-2.5*Search_Radius:round(YStart)+2.5*Search_Radius,...
                           round(XStart)-2.5*Search_Radius:round(XStart)+2.5*Search_Radius,...
                           CurrFrame);
            imshow(ImToShow, [Options.MinImageShow, Options.MaxImageShow], 'InitialMagnification', 'fit');
            hold on;
            plot(XNew+2.5*Search_Radius+1-XStart,YNew+2.5*Search_Radius+1-YStart,'go')
            drawnow;


    end

    % Measure relative to the original click on FindingImage (XStart/YStart),
    % NOT Pos.x(1)/Pos.y(1) -- those get overwritten by the loop above with
    % the tracked position on VideoMatrix frame 1, which would make the
    % first offset entry a self-referential zero and discard the real
    % finding-image-to-start-frame shift.
    Offset.x = Pos.x - XStart;
    Offset.y = Pos.y - YStart;

    disp('   All done tracking. Moving on to trace extraction.')
    disp('   ...')
    
end

function [GaussianQuantResults]  = ...
    Gaussian_Quantification(NumFrames,CurrentVirusCroppedStack,FrameNumToFindVirus,... 
    CurrentVirusBox,ImageStackMatrix,FigureHandles,Options)
                
CVB = CurrentVirusBox;
BoxToPlot = [CVB.Bottom,CVB.Left;CVB.Bottom,CVB.Right;CVB.Top,CVB.Right;CVB.Top,CVB.Left;CVB.Bottom,CVB.Left];
        
    %Now we go through the all of the frames and calculate the SUV
    %intensity using a Gaussian Fit to determine the background
        TraceBackSub = zeros(1,NumFrames);
        OffsetFrom2DFitList = zeros(1,NumFrames);
        NoiseList = zeros(1,NumFrames);
        DidFitWork = zeros(1,NumFrames);
        LocalBackgroundFromGaussFit = zeros(1,NumFrames);
        Zmin = min(min(min(double(CurrentVirusCroppedStack(:,:,:)))));
        Zmax = max(max(max(double(CurrentVirusCroppedStack(:,:,:)))));
        CurrentTraceIntegrated_AreaUnderGauss = zeros(1,NumFrames);
        CentroidList = zeros(NumFrames,2);
%         NumReCenters = 0;
%         ReDoThisFrameNum = 'n';

        if exist('FitParameters','var') == 1
            clearvars FitParameters;
        end

        CurrentFrameNumber = FrameNumToFindVirus;
    while CurrentFrameNumber <= NumFrames
       if Options.DisplayFitFigures == 'y'
            set(0,'CurrentFigure',FigureHandles.ImageWindow);
%             set(FigureHandles.ImageWindow, 'CurrentAxes',MainImageAxes);
            hold off
            imshow(ImageStackMatrix(:,:,CurrentFrameNumber), [Options.MinImageShow, Options.MaxImageShow], 'InitialMagnification', 'fit');
            hold on
        end

        CurrCroppedImage = ImageStackMatrix(...
            CurrentVirusBox.Top:CurrentVirusBox.Bottom,...
            CurrentVirusBox.Left:CurrentVirusBox.Right,...
            CurrentFrameNumber);

        %To speed up (and improve the accuracy) of the vesicle
        %gaussian fitting, we feed it the fit parameters from the
        %last 2D fit as an initial guess.  Otherwise, it performs
        %two 1D gaussian fits to feed the 2D fit algorithm.
        if CurrentFrameNumber == FrameNumToFindVirus
            InitGuesses = NaN;
        elseif CurrentFrameNumber ~= 1 && exist('FitParameters','var') == 1
            InitGuesses = FitParameters;
            if FitParameters(4) < 0 || FitParameters(3) <0
                Bob = 2;
            end
        else
            InitGuesses = NaN;
        end

        try
            [OffsetFrom2DFit, Noise, Gaussian2DFitImage, FitParameters, TrackingInfo] = ...
                Virus_Gaussian_Fit_Quant(CurrCroppedImage,InitGuesses);
            OffsetFrom2DFitList(CurrentFrameNumber) = OffsetFrom2DFit;
            CentroidList(CurrentFrameNumber,:) = [CurrentVirusBox.Left+FitParameters(1)-1,CurrentVirusBox.Top+FitParameters(2)-1];

            if Options.DisplayFitFigures == 'y'
                set(0,'CurrentFigure',FigureHandles.ImageWindow);
%                 set(FigureHandles.ImageWindow, 'CurrentAxes',MainImageAxes);
                hold on
                plot(BoxToPlot(:,2),BoxToPlot(:,1),'g-')
%                 plot(CurrentVirusBox.Left+FitParameters(1)-1,CurrentVirusBox.Top+FitParameters(2)-1,'rx')

                hold off
            end

            NoiseList(CurrentFrameNumber) = Noise;
            LocalBackgroundFromGaussFit(CurrentFrameNumber) = abs(OffsetFrom2DFit) + Noise;
            NormalizedLocalThresh = LocalBackgroundFromGaussFit(CurrentFrameNumber)/2^16;
            if (NormalizedLocalThresh > 0 && NormalizedLocalThresh < 1)
                DidFitWork(CurrentFrameNumber) = 1;
            else
                DidFitWork(CurrentFrameNumber) = 0;
            end
        catch
            disp(strcat('Gaussian fit failed during SUV quantification on frame ',...
                num2str(CurrentFrameNumber), '.  Will use previous intensity value.'))
            DidFitWork(CurrentFrameNumber) = 0;
            if CurrentFrameNumber > 1
                OffsetFrom2DFitList(CurrentFrameNumber) = OffsetFrom2DFitList(CurrentFrameNumber-1);
                NoiseList(CurrentFrameNumber) = NoiseList(CurrentFrameNumber-1);
                LocalBackgroundFromGaussFit(CurrentFrameNumber) = ...
                    abs(OffsetFrom2DFitList(CurrentFrameNumber)) + NoiseList(CurrentFrameNumber);
            elseif CurrentFrameNumber == FrameNumToFindSUVs
                disp('Gaussian fit failed on first frame of SUV quantification. How??')
                continue
            end
        end

        if Options.DisplayFitFigures == 'y'
            set(0,'CurrentFigure',FigureHandles.VirusWindow);
            mesh(double(CurrCroppedImage));
            hold on
            GaussColor = double(Options.MinImageShow)*ones(size(Gaussian2DFitImage));
            surf(Gaussian2DFitImage,GaussColor);
            alpha(0.4);
            zlim([Zmin,Zmax]);
            hold off

            %imshow(CurrentVirusCroppedStack(:,:,CurrentFrameNumber),[Options.MinImageShow, Options.MaxImageShow],'InitialMagnification','fit');
            SUVTitle = {... %strcat('SUV ', num2str(n),'/', num2str(NumOfVesicleRegions)),...
                strcat('Fr:',num2str(CurrentFrameNumber),'/',num2str(NumFrames))};
            title(SUVTitle);
            drawnow
            hold off
        end

        CurrentTraceIntegrated_AreaUnderGauss(CurrentFrameNumber) = sum(sum(double(Gaussian2DFitImage)...
            -abs(OffsetFrom2DFitList(CurrentFrameNumber))));

        %Do another thresholding inside the cropped image using the
        %background calculated from the gaussian fit
        CurrBWCropped_FromLocalBack = CurrCroppedImage > LocalBackgroundFromGaussFit(CurrentFrameNumber);
        CurrBWCropped_FromLocalBack = bwareaopen(CurrBWCropped_FromLocalBack, Options.MinSUVSize, 4);
        CroppedImageComponents = bwconncomp(CurrBWCropped_FromLocalBack,4);
        CroppedImageProps = regionprops(CroppedImageComponents, CurrCroppedImage,...
            'Eccentricity', 'PixelValues', 'Area','Centroid');
        NumOfSUVsInCroppedRegion = length(CroppedImageProps);

        if NumOfSUVsInCroppedRegion >1
            CenterTest = [];
            %Assume true SUV is the one in the middle
            MiddleX = size(CurrCroppedImage,2)/2;
            MiddleY = size(CurrCroppedImage,1)/2;
            for s = 1:NumOfSUVsInCroppedRegion
                XY_CurrSUVFound = CroppedImageProps(s).Centroid;
                CenterTest(s) = ((MiddleX-XY_CurrSUVFound(1))^2 + (MiddleY-XY_CurrSUVFound(2))^2)^.5;
            end
            SUVChosen = find(CenterTest == min(CenterTest),1);
            TraceBackSub(CurrentFrameNumber) = ...
                sum(CroppedImageProps(SUVChosen).PixelValues-LocalBackgroundFromGaussFit(CurrentFrameNumber));

        elseif NumOfSUVsInCroppedRegion == 0
            TraceBackSub(CurrentFrameNumber) = 0;
        elseif NumOfSUVsInCroppedRegion == 1
            TraceBackSub(CurrentFrameNumber) = ...
                sum(CroppedImageProps(1).PixelValues-LocalBackgroundFromGaussFit(CurrentFrameNumber));
        end
        
        % If the current intensity value is zero, don't use the same fit
        % parameters next time. I noticed a tendency for the program to
        % always collapse to zero once a fit was originally zero.
        if CurrentTraceIntegrated_AreaUnderGauss(CurrentFrameNumber) <= 0 
            clearvars FitParameters;
        end
        
        CurrentFrameNumber = CurrentFrameNumber + 1;
        
    end
    GaussianQuantResults.TraceBackSub = TraceBackSub;
    GaussianQuantResults.TraceGauss = CurrentTraceIntegrated_AreaUnderGauss;
    
end
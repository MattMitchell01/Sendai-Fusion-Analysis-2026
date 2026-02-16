function [Results,VirusDataToSave, OtherDataToSave,Options] =...
    Find_And_Analyze_Particles(VideoFilePath,VideoFilename, ...
    Options)
    
    % First, we pull out all of the data into a big matrix, find the threshold 
    % for each image in the video, and define the time vector associated with A) 
    % the video and B) the data we are actually analyzing (excluding any frames 
    % that are ignored in the trace extraction). 
    % The procedure is different depending on whether we are extracting times from the metadata or not
    if strcmp(Options.ExtractTimesFromMetaData, 'y')
        
        [ImageWidth,NumFrames,ImageHeight,BitDepth,VideoMatrix,...
        BWVideoMatrix,ThresholdToFindParticles,TotalVideoIntensity,AverageVideoIntensity,...
        RoughBackground,FigureHandles,TimeVector,VideoTimeVector,StandardBindTime,FindingImage,Options.FocusFrameNumbers] =...
            Create_Video_Matrix_Auto_Time_Vector(VideoFilePath,Options);
    
    else
        [ImageWidth,NumFrames,ImageHeight,BitDepth,VideoMatrix,...
        BWVideoMatrix,ThresholdToFindParticles,TotalVideoIntensity,AverageVideoIntensity,...
        RoughBackground,FigureHandles,FindingImage] =...
            Create_Video_Matrix_Manual_Time_Vector(VideoFilePath,Options);
   
        TimeVector = Options.TimeVector; % TimeVector is the time associated with the actual
        % trace being measured (so excluding any frames which are ignored in the trace extraction). 
        % This will be the same length or less than VideoTimeVector
        
        VideoTimeVector = Options.TimeVector; % VideoTimeVector is the time of the entire video, 
        % including any frames which are ignored in the trace analysis
        
        StandardBindTime = Options.BindingTime; % StandardBindTime is the time (relative to the time vector) in which binding is considered to occur
    end
   
    
if Options.GrabTotalIntensityOnly ~= 'y'
    % Plot a trace of the background intensity
    set(0,'CurrentFigure', FigureHandles.BackgroundTraceWindow);
    hold on
    plot(RoughBackground,'r-')
    title('Calculated Background Intensity Versus Frame Number')
        
    %Set up counters
    NumGoodParticles = 0;
    NumBadParticles = 0;
    
    % If we are determining the sample-wide XY offset, then do so now
    if strcmp(Options.DeterminePixelOffset, 'y')

        % If FFT method is prefered, perform the calculation
        if strcmp(Options.FFTOffset, 'y')
            [Offset] = Fast_Fourier_Transform(VideoMatrix, ImageWidth, ImageHeight, NumFrames, Options);
            figure;
            plot(Offset.x, 'b');
            hold on;
            plot(Offset.y, 'r');
            hold off;
    
            legend('X shift', 'Y shift');
            xlabel('Frame');
            ylabel('Pixels');
            title('XY Drift Over Time');
            
            Options.Offsets = Offset;
            
            % Prompt the user to review a plot of the X Y offset vectors
            disp(' ')
            disp('Review the XY shift plots.')
            disp('Good FFT shifts should be generally smooth and continuous over time.')
            disp('Bad shifts often show sudden large jumps or irregular behavior.')
            disp(' ')
            userResp = input('Do you approve of the XY shifting graph? (y = yes, anything else = no): ','s');
            
            if strcmpi(strtrim(userResp),'y')
                disp(' ')
                Options.SingleParticleTrackingOffset = 'n';     % In case somebody had both FFT and SPT methods enabled
                disp('FFT-based XY shifting approved. Proceeding...')
            else
                disp(' ')
                disp('Looks like FFT did not work this time.')
                disp('You will have to determine XY shifts using single particle tracking.')
        
                Options.FFTOffset= 'n';
                Options.SingleParticleTrackingOffset = 'y';
            end
        end
        
        % If FFT failed or we want to do SPT instead...
        if strcmp(Options.SingleParticleTrackingOffset,'y')
            [Offset] = Single_Particle_Tracking(VideoMatrix,BWVideoMatrix,...
                FigureHandles,Options);
            Options.Offsets = Offset;
            figure;
            plot(Offset.x, 'b');
            hold on;
            plot(Offset.y, 'r');
            hold off;
    
            legend('X shift', 'Y shift');
            xlabel('Frame');
            ylabel('Pixels');
            title('XY Drift Over Time');
            disp(' ')
            userResp = input('Do you approve of the XY shifting graph? (y = yes, anything else = no): ','s');
        end

    end

    
    % Now we find all of the particles in the finding image, decide whether 
    % they are "good" or "bad" particles, and grab the integrated intensity 
    % trace within the region of interest around the particle for the entire length 
    % of the video. All of this data is then saved.
    % NOTE: There is a for loop here that only runs for one iteration. The reason 
    % that there is a for loop is because this offers flexibility if you wish to 
    % analyze multiple frames (such as in the case if you are tracking the particles 
    % frame by frame). Obviously, in such a case you would need to considerably modify the script.
    for CurrFrameNum = Options.FrameNumToFindParticles:Options.FrameNumToFindParticles
        
        %Re-define the finding image (CurrImage). This image will be an average to boost signal-to-noise.
        %Also re-define the logical image of this averaged finding image.
        CurrImage = FindingImage;
        CurrImage = uint16(CurrImage);
        CurrentRoughBackground = mean(min(CurrImage));
        CurrentThreshold = (CurrentRoughBackground + Options.Threshold)/2^BitDepth;
        BinaryCurrImage = im2bw(CurrImage, CurrentThreshold);
        BinaryCurrImage = bwareaopen(BinaryCurrImage , Options.MinParticleSize, 8);

        %All of the isolated regions are "particles" and will
        %be analyzed.
            ParticleComponentArray = bwconncomp(BinaryCurrImage,8);
            ParticleProperties = regionprops(ParticleComponentArray, CurrImage, 'Centroid',...
                'Eccentricity', 'PixelValues', 'Area','PixelIdxList');
            NumberOfParticlesFound = length(ParticleProperties);
            
        %Re-plot the finding image
        if strcmp(Options.DisplayAllFigures,'y')
            set(0,'CurrentFigure',FigureHandles.ImageWindow);
            hold off
            imshow(CurrImage, [Options.MinImageShow, Options.MaxImageShow], 'InitialMagnification', 'fit','Border','tight');
            title('Finding Image');
            hold on
            
            set(0,'CurrentFigure',FigureHandles.BinaryImageWindow);
            imshow(BinaryCurrImage, 'InitialMagnification', 'fit','Border','tight');
            title('Finding Image, Thresholded');
            drawnow
        end
            
            
        %Analyze each particle region
        for n = 1:NumberOfParticlesFound
            CurrentParticleProperties = ParticleProperties(n);
            CurrVesX = round(ParticleProperties(n).Centroid(1)); 
            CurrVesY = round(ParticleProperties(n).Centroid(2));
                CurrentParticleProperties.Centroid = [CurrVesX, CurrVesY];
            
            %Apply many tests to see if Particle is good
            [IsParticleGood,  ~, CurrentParticleBox, ~, ~, ~, ~, ~,...
                ~, ~, ReasonParticleFailed] =...
                Simplified_Test_Goodness(CurrImage,CurrentParticleProperties,BitDepth,...
                CurrentThreshold, Options.MinParticleSize, Options.MaxEccentricity,ImageWidth,...
                ImageHeight,Options.MaxParticleSize,BinaryCurrImage,Options);
            
            if strcmp(IsParticleGood,'y')
                LineColor = 'g-';
                NumGoodParticles = NumGoodParticles + 1;
                
            elseif strcmp(IsParticleGood,'n')
                LineColor = 'r-';
                NumBadParticles = NumBadParticles + 1;
                if strcmp(Options.DisplayRejectionReasons,'y')
                    disp(ReasonParticleFailed)    
                end
            end
                    
            %Plot a box around the Particle. Green particles are "good" and 
            %red particles are "bad".
                CVB = CurrentParticleBox;
                BoxToPlot = [CVB.Bottom,CVB.Left;CVB.Bottom,CVB.Right;CVB.Top,CVB.Right;CVB.Top,CVB.Left;CVB.Bottom,CVB.Left];

            if strcmp(Options.DisplayAllFigures,'y')
                set(0,'CurrentFigure',FigureHandles.ImageWindow);
                plot(BoxToPlot(:,2),BoxToPlot(:,1),LineColor)
                hold on
                drawnow
            end
                         
            % The integrated intensity trace within the region of interest around 
            % the particle is then calculated. The background subtracted intensity 
            % trace is also calculated, where the background is calculated as an 
            % average across the entire image. If we are determining an
            % offset over time, then we apply that as well.
            if strcmp(Options.DeterminePixelOffset,'y')
                CurrentTraceIntensity = zeros(1,NumFrames);
                CurrentTraceIntensityBackSub = zeros(1,NumFrames);
                
                for Frame = 1:NumFrames
                    
                    CurrBox_Offset.Top = round(CurrentParticleBox.Top+Offset.y(Frame));
                    CurrBox_Offset.Bottom = round(CurrentParticleBox.Bottom+Offset.y(Frame));
                    CurrBox_Offset.Left = round(CurrentParticleBox.Left+Offset.x(Frame));
                    CurrBox_Offset.Right = round(CurrentParticleBox.Right+Offset.x(Frame));
                    
                    % Test if the offset box is off the viewing range. If
                    % it is, ignore and set trace to zero everywhere. If
                    % not, then quantify
                    if CurrBox_Offset.Top < 1 || CurrBox_Offset.Left < 1 ||...
                            CurrBox_Offset.Bottom > ImageHeight || CurrBox_Offset.Right > ImageWidth
                        
                        IsParticleGood = 'Ignore';
                        ReasonParticleFailed = 'Moved Off Viewing Area';
                        
                        if strcmp(Options.DisplayRejectionReasons,'y')
                            disp(ReasonParticleFailed)    
                        end
                        
                        % Now we plot this one as a yellow box
                        LineColor = 'y-';
                        
                        if strcmp(Options.DisplayAllFigures,'y')
                            set(0,'CurrentFigure',FigureHandles.ImageWindow);
                            plot(BoxToPlot(:,2),BoxToPlot(:,1),LineColor)
                            hold on
                            drawnow
                        end
                        
                        CurrentTraceIntensity = zeros(1,NumFrames);
                        CurrentTraceIntensityBackSub = zeros(1,NumFrames);
                        break
                    else
                        CurrentFrame_ParticleCropped = VideoMatrix(...
                            CurrBox_Offset.Top:CurrBox_Offset.Bottom,...
                            CurrBox_Offset.Left:CurrBox_Offset.Right,...
                            Frame);

                        CurrentTraceIntensity(Frame) = sum(sum((CurrentFrame_ParticleCropped)));

                        CurrentTraceIntensityBackSub(Frame) = CurrentTraceIntensity(Frame) -...
                            (RoughBackground(Frame).*(CurrentParticleBox.Bottom - CurrentParticleBox.Top + 1)^2)';
                    end
                end
                
            else
                CurrentParticleCroppedVideo = VideoMatrix(...
                    CurrentParticleBox.Top:CurrentParticleBox.Bottom,...
                    CurrentParticleBox.Left:CurrentParticleBox.Right,...
                    1:NumFrames);

                CurrentTraceSumArray = sum(sum((CurrentParticleCroppedVideo)));

                %Because CurrentTraceSumArray is a 3D array, the summed data is
                %transferred to row vector (CurrentTraceIntensity), so it can be
                %plotted.
                CurrentTraceIntensity = shiftdim(CurrentTraceSumArray(1,1,:),1);

                CurrentTraceIntensityBackSub = CurrentTraceIntensity -...
                    (RoughBackground(1:NumFrames).*(CurrentParticleBox.Bottom - CurrentParticleBox.Top + 1)^2)';                 
            end

            % ALTERNATE: Use gaussian fit to determine local background intensity
            % as an alternative method for background subtraction
            if strcmp(Options.UseGaussianIntensity,'y')
                [GaussianQuantResults]  = ...   
                    Gaussian_Quantification(NumFrames,CurrentParticleCroppedVideo,Options.FrameNumToFindParticles,... 
                    CurrentParticleBox,VideoMatrix,FigureHandles,Options);
            end

            %The FigureHandles.CurrentTraceWindow is set as the current figure.  The axes are 
            %cleared (in case there was a previous trace) and then the new trace is plotted.
            if strcmp(Options.DisplayAllFigures,'y')
                set(0,'CurrentFigure',FigureHandles.CurrentTraceWindow)
                cla
                plot(CurrentTraceIntensityBackSub,LineColor)
                hold on
                if strcmp(Options.UseGaussianIntensity,'y')
                    plot(GaussianQuantResults.TraceBackSub,LineColor)
                    plot(GaussianQuantResults.TraceGauss,'b-')
                end
                title(strcat('Particle :', num2str(n),'/', num2str(NumberOfParticlesFound)));
                drawnow
            else
                if rem(n,20)==0
                    set(0,'CurrentFigure',FigureHandles.CurrentTraceWindow)
                    title(strcat('Particle :', num2str(n),'/', num2str(NumberOfParticlesFound)));
                    drawnow
                end
            end

            
            %Copy the data for each particle in a structure
                VirusDataToSave(n).Trace = CurrentTraceIntensity;
                VirusDataToSave(n).Trace_BackSub = CurrentTraceIntensityBackSub;
                VirusDataToSave(n).FrameNumFound = CurrFrameNum;
                VirusDataToSave(n).Coordinates = ParticleProperties(n).Centroid;
                VirusDataToSave(n).Eccentricity = ParticleProperties(n).Eccentricity;
                VirusDataToSave(n).Area = ParticleProperties(n).Area;
                VirusDataToSave(n).FullFilePath = VideoFilePath;
                VirusDataToSave(n).StreamFilename = VideoFilename;
                VirusDataToSave(n).BoxAroundVirus = CurrentParticleBox;
                VirusDataToSave(n).TimeVector = TimeVector;
                VirusDataToSave(n).VirusIDNumber = n;
                
                VirusDataToSave(n).BindingTime = StandardBindTime;
               
                VirusDataToSave(n).FocusFrameNumbers_Shifted  = Options.FocusFrameNumbers - (Options.StartAnalysisFrameNumber-1);
                VirusDataToSave(n).IgnoreFrameNumbers_Shifted  = Options.IgnoreFrameNumbers - (Options.StartAnalysisFrameNumber-1);
                VirusDataToSave(n).IsVirusGood = IsParticleGood;
                VirusDataToSave(n).ReasonVirusFailed = ReasonParticleFailed;
                    
        end
       
    end
    
    % If chosen, plot a specific trace in a separate window, together 
    % with the intensity trace of a pH indicator ROI. If you decide to use 
    % this, you should modify User_Grab_Example_Trace.m as it is currently set up 
    % for the last trace that I grabbed.
    if Options.UserGrabExampleTrace == 'y'
        User_Grab_Example_Trace(FigureHandles,VirusDataToSave,VideoMatrix,RoughBackground,NumFrames);
    end
    
elseif Options.GrabTotalIntensityOnly == 'y'
    VirusDataToSave=[];
end

    % Record any other data that we wish to save
    OtherDataToSave.ThresholdsUsed = ThresholdToFindParticles;
    OtherDataToSave.RoughBackground = RoughBackground;
    OtherDataToSave.TotalVideoIntensity = TotalVideoIntensity; 
    OtherDataToSave.AverageVideoIntensity = AverageVideoIntensity; 
    OtherDataToSave.VideoTimeVector = VideoTimeVector;    
    OtherDataToSave.FocusFrameNumbers = Options.FocusFrameNumbers; 
    OtherDataToSave.Options = Options;
    OtherDataToSave.StandardBindTime = StandardBindTime;
    OtherDataToSave.FindingImage = FindingImage;
    OtherDataToSave.FrameNumToFindParticles = Options.FrameNumToFindParticles;
   
    % Plot the total video intensity
    set(0,'CurrentFigure',FigureHandles.CurrentTraceWindow)
    cla
    plot(AverageVideoIntensity)
    title('Average frame intensity versus frame number')
    
    % Record any results that we wish to display in the command prompt window
    Results.Filename = VideoFilename;
    
    if Options.GrabTotalIntensityOnly ~= 'y'
        Results.NumBadParticles = NumBadParticles;
        Results.NumGoodParticles = NumGoodParticles;
    end
    
    Display_All_Particles(VirusDataToSave,FigureHandles)
end
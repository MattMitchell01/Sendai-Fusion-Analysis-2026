function [ImageWidth,NumFramesActuallyAnalyzed,ImageHeight,BitDepth,VideoMatrix,...
    BWVideoMatrix,ThresholdToFindParticles,TotalVideoIntensity,AverageVideoIntensity,...
    RoughBackground,FigureHandles,TimeVector,VideoTimeVector,StandardBindTime,FindingImage, FocusIndices] =...
    Create_Video_Matrix_Auto_Time_Vector(VideoFilePath,...
    Options)

    % First, we use the bio formats open function to get the data. Then pull out the information we need
    % If you want to mess around with this, you will need to look at the bio formats help documentation 
    % for OME.TIF files
    data = bfopen(VideoFilePath);
    

    % Extract the OME-XML metadata object returned by bfopen
    % This contains image dimensions, pixel type, timing, Z positions, etc.
    metadata = data{1,4};

    % Read spatial dimensions of the image frames
    ImageWidth  = metadata.getPixelsSizeX(0).getValue();  % pixels
    ImageHeight = metadata.getPixelsSizeY(0).getValue();  % pixels
    
    % First frame to include in analysis
    StartFrameNumber = Options.StartAnalysisFrameNumber;
    
    % Total number of frames in the OME-TIFF.
    NumFramesInFullVideo = metadata.getPixelsSizeT(0).getValue();
    
    % Determine the pixel type/depth
    pixelsTypeStr = char(metadata.getPixelsType(0));   % e.g. 'uint16'
    BitDepth = str2double(regexprep(pixelsTypeStr,'\D',''));
    
    % Extract Frame Timestamps from the metadata
    for i = 1:NumFramesInFullVideo
        VideoTimeVector(i) = str2double(metadata.getPlaneDeltaT(0,i-1).value()); % in milliseconds (ms)
    end



    %
    % Now we want to extract the Focus Indices using the Zposition Data
    % from the Metadata
    %
    
    % Preallocate vector to store Z position for each frame.
    PositionZ = nan(1, NumFramesInFullVideo);
    
    for i = 1:NumFramesInFullVideo
        z = metadata.getPlanePositionZ(0, i-1);
        % Read the Z position for plane i (0-based indexing in OME).
    
        if ~isempty(z)
            PositionZ(i) = z.value().doubleValue();
        end
    end
    
    % The list of Frame Indices where focus events will occur
    FocusIndices = [];
    
    for i = StartFrameNumber:(NumFramesInFullVideo-1)
        % Compare consecutive Z positions starting from the analysis frame.
    
        if ~isnan(PositionZ(i)) && ~isnan(PositionZ(i+1)) && ...
           (PositionZ(i) ~= PositionZ(i+1))
            FocusIndices(end+1) = i + 1;
            % Focus event appears in the frame AFTER the Z value changes.
        end
    end

    
    % User will now check Focus Indices
    disp('Detected Focus Indices:');
    disp(FocusIndices);
    
    while true
        resp = input('Do you approve of these focus indices? (y/n): ', 's');
    
        if strcmpi(resp, 'y')
            disp('Great — moving on!');
            fprintf('\n')
            fprintf('\n')
            break
        elseif strcmpi(resp, 'n')
            disp('READ READ READ --> Then you will have to find the focus indices yourself and use the input .txt file.');
            error('Execution stopped by user.');
        else
            disp('Invalid input. Please enter y or n.');
        end
    end


       

    % Depending on the situation, we truncate the time vector. TimeVector is the time associated with the actual
    % trace being measured (so excluding any frames which are ignored in the trace extraction). 
    % This will be the same length or less than VideoTimeVector
        if strcmp(Options.DetermineBindingTimeFromTimestamp, 'y')
            % If we are determining the time zero from an image timestamp, then 
            % we do that here. We will assume that the standard bind time is 0, 
            % and make that so by subtracting the time zero frame
            StandardBindTime = 0;
            TimeVector = VideoTimeVector - VideoTimeVector(Options.TimeZeroFrameNumber);
            TimeVector = TimeVector(StartFrameNumber:end); 
                % exclude the first image(s) from the trace analysis, since it was 
                % only used as a time zero marker, or otherwise excluded
        else
            StandardBindTime = Options.BindingTime;
            TimeVector = VideoTimeVector - VideoTimeVector(1);
        end

        if isnan(Options.FrameNumberLimit) || strcmp(Options.CombineMultipleVideos,'y')
        else
            if Options.FrameNumberLimit > NumFramesInFullVideo
                disp('Oops! You manually chose a frame number limit larger than the number of frames in the video.  Double check your options.');
            end
            NumFramesInFullVideo = Options.FrameNumberLimit;
            TimeVector = TimeVector(1:NumFramesInFullVideo);
        end
  
        NumFramesActuallyAnalyzed = length(StartFrameNumber:NumFramesInFullVideo);
        
        % Convert time vector to minutes (default is ms)
        TimeVector = TimeVector/(1000*60);

    % Now we pre-allocate data to store our video matrix, and display our finding image
        VideoMatrix = zeros(ImageHeight, ImageWidth, NumFramesActuallyAnalyzed, 'uint16');
        
        %Create a logical matrix the same size as the video matrix.
        BWVideoMatrix = VideoMatrix > 0; 
        
        %Preallocate various vectors as well
        ThresholdToFindParticles = zeros(NumFramesActuallyAnalyzed,1);
        TotalVideoIntensity = zeros(NumFramesActuallyAnalyzed,1);
        AverageVideoIntensity = zeros(NumFramesActuallyAnalyzed,1);
        RoughBackground = zeros(NumFramesActuallyAnalyzed,1);
        
        %Set up figures
        [FigureHandles] = Setup_Figures(Options);
        
        % Display the finding image first thing. This will help you determine 
        % if you have set a good threshold.
            VideoInfo.Width =  ImageWidth; %in pixels
            VideoInfo.Height =  ImageHeight; %in pixels
            VideoInfo.BitDepth = BitDepth;
        [FindingImage] = Display_Finding_Image(VideoFilePath,Options,VideoInfo,FigureHandles,data,StartFrameNumber);
        
    %This for loop populates the VideoMatrix with the data from each image
    %in the video.  The 1st two dimensions are the x,y of the image plane and the 3rd 
    %dimension is the frame number. We adjust for the start frame to account for the times when we exclude 
    %the first frame(s) during automatic time extraction.
    
        NewFrameNum = 0;
    for b = StartFrameNumber:NumFramesActuallyAnalyzed + (StartFrameNumber-1)
        CurrentFrameImage = data{1, 1}{b, 1};
        
        NewFrameNum = NewFrameNum + 1;
        
        if strcmp(Options.TopHatBackgroundSubtraction, 'y')
        
            CurrentFrameImage = Top_Hat_Background_Subtraction(CurrentFrameImage, Options);

        end
            
       
        VideoMatrix(:,:,NewFrameNum) = CurrentFrameImage;
  
                            
        % For each frame, the background intensity, average intensity, and 
        % integrated intensity are calculated. The threshold for each image 
        % is also calculated (this would be used if particles were being found 
        % or tracked in each image, which is currently not being done).
        RoughBackground(NewFrameNum) = mean(median(CurrentFrameImage));
        TotalVideoIntensity(NewFrameNum) = sum(sum(CurrentFrameImage));
        AverageVideoIntensity(NewFrameNum) = mean(mean(CurrentFrameImage));
        ThresholdToFindParticles(NewFrameNum) = (RoughBackground(NewFrameNum) + Options.Threshold)/2^BitDepth;
        
        %We apply the threshold to create a big logical matrix
        CurrThresh = ThresholdToFindParticles(NewFrameNum);
        BWVideoMatrix(:,:,NewFrameNum) = im2bw(CurrentFrameImage, CurrThresh);
        BWVideoMatrix(:,:,NewFrameNum) = bwareaopen(BWVideoMatrix(:,:,NewFrameNum), Options.MinParticleSize, 8);
        
        %Display the progress of loading the frames
        if rem(b,20)==0
            set(0,'CurrentFigure',FigureHandles.CurrentTraceWindow);
            title(strcat('Compiling Frame :', num2str(NewFrameNum),'/', num2str(NumFramesActuallyAnalyzed)));
            drawnow
        end
        
    end
    
    if strcmp(Options.CombineMultipleVideos,'y')
        
        for i= 2:Options.NumberOfVideosToCombine
            [VideoMatrix,RoughBackground,...
                TotalVideoIntensity,AverageVideoIntensity,ThresholdToFindParticles,...
                BWVideoMatrix,NumFramesActuallyAnalyzed,TimeVector] = ...
            Stitch_Next_Video(i,Options,VideoMatrix,RoughBackground,...
                TotalVideoIntensity,AverageVideoIntensity,ThresholdToFindParticles,...
                BWVideoMatrix,NumFramesActuallyAnalyzed,TimeVector,FigureHandles);
        end
  
    end
    
 end

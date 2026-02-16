function [ImageWidth,NumFramesActuallyAnalyzed,ImageHeight,BitDepth,VideoMatrix,...
    BWVideoMatrix,ThresholdToFindParticles,TotalVideoIntensity,AverageVideoIntensity,...
    RoughBackground,FigureHandles,FindingImage] =...
    Create_Video_Matrix_Manual_Time_Vector(VideoFilePath,...
    Options)

    %The first image in the video is read and displayed. A 3D array
    %(VideoMatrix) is created which will contain the data for all images in
    %the video.  VideoMatrix is pre-allocated with zeros to make the 
    %for loop faster.
    
    VideoInfo = imfinfo(VideoFilePath);
    if isnan(Options.FrameNumberLimit)
        NumFrames = length(VideoInfo);
    else
        NumFrames = Options.FrameNumberLimit;
    end
    
    StartFrameNumber = Options.StartAnalysisFrameNumber;
    NumFramesActuallyAnalyzed = length(StartFrameNumber:NumFrames);
    
        ImageWidth = VideoInfo.Width; %in pixels
        ImageHeight = VideoInfo.Height; %in pixels
        BitDepth = VideoInfo.BitDepth;
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
        [FindingImage] = Display_Finding_Image(VideoFilePath,Options,VideoInfo,FigureHandles,NaN,Options.StartAnalysisFrameNumber);
        
    %This for loop populates the VideoMatrix with the data from each image
    %in the video.  The 1st two dimensions are the x,y of the image plane and the 3rd 
    %dimension is the frame number.
    for b = StartFrameNumber:NumFramesActuallyAnalyzed + (StartFrameNumber-1)
        CurrentFrameImage = imread(VideoFilePath,b);
        
        NewFrameNum = b - (StartFrameNumber-1);
        % We adjust for the start frame to account for the times when we exclude 
        % the first frame(s).
        
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
            title(strcat('Loading Frame :', num2str(b),'/', num2str(NumFramesActuallyAnalyzed)));
            drawnow
        end
        
    end
    
 end
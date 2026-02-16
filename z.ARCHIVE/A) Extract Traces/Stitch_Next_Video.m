function [VideoMatrix,RoughBackground,...
                TotalVideoIntensity,AverageVideoIntensity,ThresholdToFindParticles,...
                BWVideoMatrix,NumFramesActuallyAnalyzed,TimeVector] = ...
            Stitch_Next_Video(CurrentVideoNum,Options,VideoMatrix,RoughBackground,...
                TotalVideoIntensity,AverageVideoIntensity,ThresholdToFindParticles,...
                BWVideoMatrix,NumFramesActuallyAnalyzed,TimeVector,FigureHandles)
            
    % First, we use the bio formats open function to get the data, and then pull out 
    % the relevant metadata by creating the appropriate key. If you want to mess 
    % around with this, you will need to look at the bio formats help documentation 
    % for OME.TIF files
    
    VideoFilePath = strcat(Options.DefaultPathnamesToCombine{1,CurrentVideoNum},Options.VideoFilenamesToCombine{1,CurrentVideoNum});
    data = bfopen(VideoFilePath);

    metadata = data{1, 2};
    NumFramesInThisVideo = Options.NumFrames(CurrentVideoNum);
    
    if NumFramesInThisVideo < 101
        for i = 1:NumFramesInThisVideo
            if i < 11
                ThisVideoTimeVector(i) = str2double(metadata.get(strcat('Plane #','0', num2str(i-1),...
                ' ElapsedTime-ms')));
                % VideoTimeVector is the time of the entire video, 
                % including any frames which are ignored in the trace analysis
            else
                ThisVideoTimeVector(i) = str2double(metadata.get(strcat('Plane #', num2str(i-1),...
                ' ElapsedTime-ms')));
            end
        end
        BitDepth = str2double(metadata.get(strcat('Plane #00 BitDepth')));
        ImageWidth = str2double(metadata.get(strcat('Plane #00 Width'))); %in pixels
        ImageHeight = str2double(metadata.get(strcat('Plane #00 Height'))); %in pixels
        
    else
        for i = 1:NumFramesInThisVideo
            if i < 11
                ThisVideoTimeVector(i) = str2double(metadata.get(strcat('Plane #','00', num2str(i-1),...
                ' ElapsedTime-ms')));
            elseif i >= 11 && i < 101
                ThisVideoTimeVector(i) = str2double(metadata.get(strcat('Plane #','0', num2str(i-1),...
                ' ElapsedTime-ms')));
            else
                ThisVideoTimeVector(i) = str2double(metadata.get(strcat('Plane #',num2str(i-1),...
                ' ElapsedTime-ms')));
            end

        end

        BitDepth = str2double(metadata.get(strcat('Plane #000 BitDepth')));
        ImageWidth = str2double(metadata.get(strcat('Plane #000 Width'))); %in pixels
        ImageHeight = str2double(metadata.get(strcat('Plane #000 Height'))); %in pixels
    end
  
        
        % Convert time vector to minutes (default is ms)
        ThisVideoTimeVector = ThisVideoTimeVector/(100*60);
        
    % Now we pre-allocate data to store our video matrix, and display our finding image
        ThisVideoMatrix = zeros(ImageHeight, ImageWidth, NumFramesInThisVideo, 'uint16');
        
        %Create a logical matrix the same size as the video matrix.
        BWThisVideoMatrix = ThisVideoMatrix > 0; 
        
        %Preallocate various vectors as well
        ThresholdToFindParticles_ThisVideo = zeros(NumFramesInThisVideo,1);
        TotalVideoIntensity_ThisVideo = zeros(NumFramesInThisVideo,1);
        AverageVideoIntensity_ThisVideo = zeros(NumFramesInThisVideo,1);
        RoughBackground_ThisVideo = zeros(NumFramesInThisVideo,1);

        
    %This for loop populates the VideoMatrix with the data from each image
    %in the video.  The 1st two dimensions are the x,y of the image plane and the 3rd 
    %dimension is the frame number.
    for b = 1:NumFramesInThisVideo
        CurrentFrameImage = data{1, 1}{b, 1};
        
        ThisVideoMatrix(:,:,b) = CurrentFrameImage;
                            
        % For each frame, the background intensity, average intensity, and 
        % integrated intensity are calculated. The threshold for each image 
        % is also calculated (this would be used if particles were being found 
        % or tracked in each image, which is currently not being done).
        RoughBackground_ThisVideo(b) = mean(median(CurrentFrameImage));
        TotalVideoIntensity_ThisVideo(b) = sum(sum(CurrentFrameImage));
        AverageVideoIntensity_ThisVideo(b) = mean(mean(CurrentFrameImage));
        ThresholdToFindParticles_ThisVideo(b) = (RoughBackground_ThisVideo(b) + Options.Threshold)/2^BitDepth;
        
        %We apply the threshold to create a big logical matrix
        CurrThresh = ThresholdToFindParticles_ThisVideo(b);
        BWThisVideoMatrix(:,:,b) = im2bw(CurrentFrameImage, CurrThresh);
        BWThisVideoMatrix(:,:,b) = bwareaopen(BWThisVideoMatrix(:,:,b), Options.MinParticleSize, 8);
        
        %Display the progress of loading the frames
        if rem(b,20)==0
            set(0,'CurrentFigure',FigureHandles.CurrentTraceWindow);
            title(strcat('Stitching Frame :', num2str(b),'/', num2str(NumFramesInThisVideo),'; Movie #',num2str(CurrentVideoNum)));
            drawnow
        end
        
    end 
    
    
    % Now we combine the video matrix and relevant variables extracted for the current 
    % video with the video matrix and variables extracted from the previous video(s)
    VideoMatrix = cat(3,VideoMatrix,ThisVideoMatrix);
    RoughBackground = [RoughBackground; RoughBackground_ThisVideo];
    TotalVideoIntensity = [TotalVideoIntensity; TotalVideoIntensity_ThisVideo];
    AverageVideoIntensity = [AverageVideoIntensity; AverageVideoIntensity_ThisVideo];
    ThresholdToFindParticles = [ThresholdToFindParticles; ThresholdToFindParticles_ThisVideo];
    BWVideoMatrix = cat(3,BWVideoMatrix,BWThisVideoMatrix);

        % Now we tack on the time vector for the current video to the end of the time 
        % vector from the last video, after adding the time And taking it relative to the last video
        ThisVideoTimeVector = ThisVideoTimeVector + TimeVector(end) + Options.TimeDelayBetweenVideos(CurrentVideoNum-1);
        TimeVector = [TimeVector ThisVideoTimeVector];
        NumFramesActuallyAnalyzed = length(TimeVector);
    
end
function []= Quick_Scan_Threshold(FileOptions,NumberOfParameters,ZoomImage)

BasicOptions = FileOptions(1).Parameter(1).Options;

    VideoFileName = BasicOptions.VideoFilename;
    DefaultPathname = BasicOptions.DefaultPathname;
    VideoFilePath = strcat(DefaultPathname,VideoFileName);

    if strcmp(BasicOptions.ExtractAnalysisInputsFromConsole, 'y')
        [BasicOptions] = Console_Extract_Analysis_Inputs(BasicOptions);
    elseif strcmp(BasicOptions.ExtractInputsFromTextFile, 'y')
        [BasicOptions] = Extract_Analysis_Inputs(BasicOptions,strcat(DefaultPathname,BasicOptions.AnalysisTextFilename));
    end


    FrameRange = BasicOptions.FrameNumToFindParticles - BasicOptions.FindFramesToAverage:BasicOptions.FrameNumToFindParticles + BasicOptions.FindFramesToAverage;
    NumFrames = length(FrameRange);
    StartFrameNumber = BasicOptions.StartAnalysisFrameNumber;
    
if strcmp(BasicOptions.ExtractTimesFromMetaData, 'y')    
    
    data = bfopen(VideoFilePath);

    metadata = data{1, 2};
    NumFramesInFullVideo = size(data{1,1},1);
    
    if NumFramesInFullVideo < 101
        BitDepth = str2double(metadata.get(strcat('Plane #00 BitDepth')));
        ImageWidth = str2double(metadata.get(strcat('Plane #00 Width'))); %in pixels
        ImageHeight = str2double(metadata.get(strcat('Plane #00 Height'))); %in pixels
        
    elseif NumFramesInFullVideo < 1001
        BitDepth = str2double(metadata.get(strcat('Plane #000 BitDepth')));
        ImageWidth = str2double(metadata.get(strcat('Plane #000 Width'))); %in pixels
        ImageHeight = str2double(metadata.get(strcat('Plane #000 Height'))); %in pixels
    else
        BitDepth = str2double(metadata.get(strcat('Plane #0000 BitDepth')));
        ImageWidth = str2double(metadata.get(strcat('Plane #0000 Width'))); %in pixels
        ImageHeight = str2double(metadata.get(strcat('Plane #0000 Height'))); %in pixels
    end

        VideoMatrix = zeros(ImageHeight, ImageWidth, length(FrameRange), 'uint16');

        FrameCounter = 0;
    for b = FrameRange

        CurrentFrameImage = data{1, 1}{b, 1};

        FrameCounter = FrameCounter + 1;
        % We adjust for the start frame to account for the times when we exclude 
        % the first frame(s) during automatic time extraction.
        VideoMatrix(:,:,FrameCounter) = CurrentFrameImage;
    end
    
else
    VideoInfo = imfinfo(VideoFilePath);
    
        ImageWidth = VideoInfo.Width; %in pixels
        ImageHeight = VideoInfo.Height; %in pixels
        BitDepth = VideoInfo.BitDepth;   
        
        VideoMatrix = zeros(ImageHeight, ImageWidth, length(FrameRange), 'uint16');

    for b = 1:NumFrames
            CurrentRealFrameNumber = FrameRange(b);
            CurrentFrameImage = imread(VideoFilePath,CurrentRealFrameNumber);
            VideoMatrix(:,:,b) = CurrentFrameImage;
    end

end

        
    %Set up figures
    [FigureHandles] = Setup_Figures(BasicOptions);

    ZoomCoords = [0.9*ImageHeight/4,...
                3.1*ImageHeight/4,...
                0.5*ImageWidth/4,...
                3.5*ImageWidth/4];
    ZoomCoords = round(ZoomCoords);

%Display the finding image. 
    FindingImage = mean(VideoMatrix(:,:,:),3);
    FindingImage = uint16(FindingImage);
    CurrentRoughBackground = mean(min(FindingImage));
    
    set(0,'CurrentFigure',FigureHandles.MasterWindow);
    FindingImagePosition = 2;
    set(FigureHandles.MasterWindow,'CurrentAxes',FigureHandles.SubHandles(FindingImagePosition));
    cla
    if strcmp(ZoomImage,'y')
        imshow(FindingImage(ZoomCoords(1):ZoomCoords(2),ZoomCoords(3):ZoomCoords(4)),...
            [BasicOptions.MinImageShow, BasicOptions.MaxImageShow], 'InitialMagnification', 'fit','Border','tight');
        title(strcat('Center Zoom Image, Frame=',num2str(BasicOptions.FrameNumToFindParticles)))
    else 
        imshow(FindingImage, [BasicOptions.MinImageShow, BasicOptions.MaxImageShow], 'InitialMagnification', 'fit','Border','tight');
        title('Raw Image')
    end
    hold on

% Display each of the thresholded images
    ImageOrder = [1, 4, 5, 6, 3];
    
    
    for j=  1:NumberOfParameters
        
        if j > length(ImageOrder)
            disp('Error: You Chose Too Many Parameters');
            ThisWillCauseProgramToEnd
        end
        
        CurrentHandlePosition = ImageOrder(j);
        
        CurrentOptions = FileOptions(1).Parameter(j).Options;
        FindingThreshold(j) = (CurrentRoughBackground + CurrentOptions.Threshold)/2^BitDepth;
        BinaryFindingImages(:,:,j) = im2bw(FindingImage, FindingThreshold(j));
        BinaryFindingImages(:,:,j) = bwareaopen(BinaryFindingImages(:,:,j), CurrentOptions.MinParticleSize, 8);
        ParticleComponentArray = bwconncomp(BinaryFindingImages(:,:,j),8);
        ParticleProperties = regionprops(ParticleComponentArray, FindingImage, 'Centroid',...
            'Eccentricity', 'PixelValues', 'Area','PixelIdxList');
        NumberOfParticlesFound(j) = length(ParticleProperties);
        
        %Plot the images
            set(0,'CurrentFigure',FigureHandles.MasterWindow);
            set(FigureHandles.MasterWindow,'CurrentAxes',FigureHandles.SubHandles(CurrentHandlePosition));
            cla

            if strcmp(ZoomImage,'y')
                imshow(BinaryFindingImages(ZoomCoords(1):ZoomCoords(2),ZoomCoords(3):ZoomCoords(4),j),...
                    'InitialMagnification', 'fit','Border','tight');
            else 
                imshow(BinaryFindingImages(:,:,j), 'InitialMagnification', 'fit','Border','tight');
            end
            title(strcat('TH=',num2str(CurrentOptions.Threshold),'; Tot Num Particles=',num2str(NumberOfParticlesFound(j))))
            drawnow

    end

disp("   Thresholds scanned and images displayed.")
disp("   Hopefully you see a good threshold!")
disp('   Program terminating below...')
disp('====================================')
ThisWillCauseProgramToEnd
end
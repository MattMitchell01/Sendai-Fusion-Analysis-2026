function [FindingImage] = Display_Finding_Image(VideoFilePath,Options,VideoInfo,FigureHandles,data,StartFrameNumber)

FrameRange = Options.FrameNumToFindParticles - Options.FindFramesToAverage:...
    Options.FrameNumToFindParticles + Options.FindFramesToAverage;

ImageWidth = VideoInfo.Width; %in pixels
ImageHeight = VideoInfo.Height; %in pixels
BitDepth = VideoInfo.BitDepth;

VideoMatrix = zeros(ImageHeight, ImageWidth, length(FrameRange), 'uint16');

if strcmp(Options.ExtractTimesFromMetaData, 'y')
        
    FrameCounter = 0;
    for b = FrameRange

        CurrentFrameImage = data{1, 1}{b, 1};

        FrameCounter = FrameCounter + 1;
        VideoMatrix(:,:,FrameCounter) = CurrentFrameImage;
    end

else
    FrameCounter = 0;
    for b = FrameRange
        CurrentFrameImage = imread(VideoFilePath,b);
        
        FrameCounter = FrameCounter + 1;
        VideoMatrix(:,:,FrameCounter) = CurrentFrameImage;
    end
end  

    FindingImage = mean(VideoMatrix(:,:,:),3);

%Display the finding image, including logical image. This is done 
%before loading the rest of the frames, just in case the threshold 
%needs to be changed and you don't want to have to wait for it 
%to load up all the frames every time.
    FindingImage = uint16(FindingImage);
    CurrentRoughBackground = mean(min(FindingImage));
    FindingThreshold = (CurrentRoughBackground + Options.Threshold)/2^BitDepth;
    BinaryFindingImage = im2bw(FindingImage, FindingThreshold);
    BinaryFindingImage = bwareaopen(BinaryFindingImage , Options.MinParticleSize, 8);

%Plot the images
    set(0,'CurrentFigure',FigureHandles.ImageWindow);
    hold off
    imshow(FindingImage, [Options.MinImageShow, Options.MaxImageShow], 'InitialMagnification', 'fit','Border','tight');
    hold on

    set(0,'CurrentFigure',FigureHandles.BinaryImageWindow);
    imshow(BinaryFindingImage, 'InitialMagnification', 'fit','Border','tight');
    drawnow
end
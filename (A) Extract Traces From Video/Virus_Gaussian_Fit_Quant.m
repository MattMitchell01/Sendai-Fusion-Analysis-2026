function [OffsetFrom2DFit, Noise, Gaussian2DFit, FittedParameters2D, TrackingInfo] =...
    Vesicle_Gaussian_Fit_For_Quant(ImageToFit,InitGuesses)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Vesicle_Gaussian_Fit.m
%
%This program was adapted from one written by PMB.  It takes the cropped
%image of a vesicle and fits it to a 2D Gaussian.  From that fit, a local
%background is calculated which will be used to determine the absolute
%intensity of the vesicle being analyzed.
%
%The 2D Gaussian fit is actually an unconstrained minimization and is
%performed in several steps.  First, the center of mass of the image is
%calculated.  The values from that calculation are then used to perform a
%1D Gaussian fit in the x and y directions of the image, using the slices
%of the image that bisect the centroid.  Finally, the values from the two
%1D fits are used as guesses for the 2D Gaussian fit.
%
%Written April, 2011 by Bob Rawle
%Modified Aug, 2012 by Bob Rawle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The image to fit is converted to double precision numbers.
ImageToFit = double(ImageToFit);

%The tolerance will be used for the minimizations below
Tolerance = 1e-9;
%The Optimization options are set for the minimizations
OptimizationOptions = optimset('Display','off','TolFun',Tolerance,'LargeScale','off');

[ImageToFitSizeY ImageToFitSizeX] = size(ImageToFit);

%The weighted centroid (i.e. center of mass) of the cropped image is
%calculated, as well as the standard deviation of the centroid position.
%These values will be used as initial guesses for the X and Y 1D Gaussian
%fits.  Note: Center_Of_Mass is a nested function (below)
if isnan(InitGuesses)
    [CentroidXFromCoM,CentroidYFromCoM,StdevXFromCoM,StdevYFromCoM] =...
        Center_Of_Mass(ImageToFit, ImageToFitSizeY, ImageToFitSizeX);
    MaxIntensity = max(max(ImageToFit));
    Offset = ImageToFit(1,1);

    %Using the values obtained from the Center of Mass calculation for the
    %centroid X position, a 1D Gaussian fit (an unconstrained minimization
    %actually) is performed for the X direction.  This yields new guesses for the CenterX and
    %StdevX which will be used in the 2D fit.
        %First we choose the slice of the image along the x axis which runs
        %through the Y Centroid obtained from the CoM calculation above, set up
        %the indices and the initial guesses which will be fed to the fit.
        CenterXSliceOfImage = ImageToFit(round(CentroidYFromCoM),:);
        X1DIndices = 1:ImageToFitSizeX; 
        InitGuessesFor1DXFit = [CentroidXFromCoM,StdevXFromCoM,MaxIntensity,Offset]; 

        %Then we run the 1D X fit (i.e. minimization of a 1D gaussian over our
        %data (center slice) using the initial guesses specified).  Note:
        %Gaussian_1D is a nested function below.
        FittedParameters1DX = fminunc(@Gaussian_1D,InitGuessesFor1DXFit,OptimizationOptions,CenterXSliceOfImage,X1DIndices);

        %The outputs of the fit are assigned to appropriate variables.
        CentroidXFrom1DFit = FittedParameters1DX(1);
        StdevXFrom1DFit = FittedParameters1DX(2);
        MaxIntensityFrom1DXFit = FittedParameters1DX(3);
        OffsetFrom1DXFit = FittedParameters1DX(4);


    %Using the values obtained from the Center of Mass calculation for the
    %centroid y position, and the 1D X fit, another 1D Gaussian fit is performed
    %for the Y direction.  This yields new guesses for the CenterY and
    %StdevY which will be used in the 2D fit.
        %First we choose the slice of the image along the y axis which runs
        %through the X Centroid obtained from the 1D gaussian fit above, set up
        %the indices and the initial guesses which will be fed to the fit.
        CenterYSliceOfImage = ImageToFit(:,round(CentroidXFrom1DFit))';
        Y1DIndices = 1:ImageToFitSizeY; %A vector of the y indices along the slice
        InitGuessesFor1DYFit = [CentroidYFromCoM,StdevYFromCoM,MaxIntensityFrom1DXFit,OffsetFrom1DXFit]; %The initial parameters which will be put into the 1D fit

        %Then we run the 1D Y fit (i.e. minimization of a 1D gaussian over our
        %data (center slice) using the initial guesses specified). Note:
        %Gaussian_1D is a nested function below.
        FittedParameters1DY = fminunc(@Gaussian_1D,InitGuessesFor1DYFit,OptimizationOptions,CenterYSliceOfImage,Y1DIndices);

        %The outputs of the fit are assigned to appropriate variables.
        CentroidYFrom1DFit = FittedParameters1DY(1);
        StdevYFrom1DFit = FittedParameters1DY(2);
        MaxIntensityFrom1DFit = FittedParameters1DY(3);
        OffsetFrom1DFit = FittedParameters1DY(4);
        InitGuessesFor2DFit = [CentroidXFrom1DFit,CentroidYFrom1DFit,StdevXFrom1DFit,StdevYFrom1DFit,MaxIntensityFrom1DFit, OffsetFrom1DFit];
else
    %If initial guesses were provided by the input, then we use those
    %instead of doing the 1D Gaussian fits described above (i.e. we
    %probably fit this same vesicle in the previous frame and so we use the
    %fit params from that as the initial guesses).
    InitGuessesFor2DFit = InitGuesses;
end

%Finally, the 2D Gaussian fit (another unconstrained minimization) is
%performed using the whole cropped image and the guesses obtained from the
%1D fits.
    %The indices are defined, which will be sent to the 2D Gaussian nested
    %function.  The initial guesses are also specified.
    [X,Y] = meshgrid(1:ImageToFitSizeX,1:ImageToFitSizeY);
    
    %The 2D fit is performed (minimization of a 2D gaussian minus our data
    %using the initial guesses specified).  Note: Gaussian_2D nested
    %function below.
    FittedParameters2D = fminunc(@Gaussian_2D,InitGuessesFor2DFit,OptimizationOptions,ImageToFit,X,Y);
    
    %The outputs of the fit are assigned to appropriate variables.
    CentroidXFrom2DFit = FittedParameters2D(1);
    CentroidYFrom2DFit = FittedParameters2D(2);
    StdevXFrom2DFit = FittedParameters2D(3);
    StdevYFrom2DFit = FittedParameters2D(4);
    MaxIntensityFrom2DFit = FittedParameters2D(5);
    OffsetFrom2DFit = FittedParameters2D(6);
    TotalStdevFrom2DFit = sqrt(StdevXFrom2DFit.^2+StdevYFrom2DFit.^2);

%Capture tracking information
    TrackingInfo.YShift = CentroidYFrom2DFit-ImageToFitSizeY/2; 
    TrackingInfo.XShift = CentroidXFrom2DFit - ImageToFitSizeX/2; 
    
    DistMoved = ((TrackingInfo.YShift)^2 + (TrackingInfo.XShift)^2)^.5;
    
    %If SUV has moved too far, then flag it for tracking, otherwise, leave
    %it be.
    if (DistMoved > ImageToFitSizeX/5)
        TrackingInfo.ShouldITrack = 'y';
    else
        TrackingInfo.ShouldITrack = 'n';
    end
%Using the 2D Gaussian fit, the noise is calculated and that will be used to
%determine the local threshold (i.e. background) for the current vesicle 
%being analyzed.
[x,y] = meshgrid(1:ImageToFitSizeX,1:ImageToFitSizeY);
Gaussian2DFit = abs(MaxIntensityFrom2DFit)*(exp(-0.5*(x-CentroidXFrom2DFit).^2./...
        (StdevXFrom2DFit^2)-0.5*(y-CentroidYFrom2DFit).^2./...
        (StdevYFrom2DFit^2)))+OffsetFrom2DFit;
Noise = std2(Gaussian2DFit-ImageToFit); 


    
    %Now the nested functions defining 1D and 2D Gaussians are created
    function [z] = Gaussian_1D(p,ImageSliceVector,XIndices)
        %p is the vector containing the initial guesses for the 1D gaussian
        %function
        
        zx = p(3)*exp(-0.5*(XIndices-p(1)).^2./(p(2)^2))+p(4) - ImageSliceVector;

        z = sum(zx.^2);
    end

    function [z] = Gaussian_2D(p,ImageToFit,X,Y)
        %X and Y are indices defined by meshgrid and are the same size
        %as the ImageToFit.
        %p is the vector containing the initial guesses for the 2D gaussian
        %function.
        ztmp = p(5)*(exp(-0.5*(X-p(1)).^2./(p(3)^2)-0.5*(Y-p(2)).^2./(p(4)^2))) +p(6) - ImageToFit;

        z = sum(sum(ztmp.^2));
    end
    
    %The nested Center_Of_Mass function is defined
    function [WeightedCentroidX,WeightedCentroidY,StdevXFromCoM,StdevYFromCoM] = ...
            Center_Of_Mass(ImageToFit, ImageToFitSizeY, ImageToFitSizeX)

        %We sum across the x and y directions to yield vectors.
        VectorOfSumX = sum(ImageToFit);
        VectorOfSumY = sum(ImageToFit');
        
        %If an element in either of the two vectors somehow ended up less
        %than zero, we change it to zero by logical indexing.
        VectorOfSumX = VectorOfSumX.*(VectorOfSumX>0);
        VectorOfSumY = VectorOfSumY.*(VectorOfSumY>0);
        
        %X and Y indices the size of the cropped image are created.
        XIndices = [1:ImageToFitSizeX];
        YIndices = [1:ImageToFitSizeY];
        
        %The weighted centroid X and Y are calculated (i.e. the center of
        %mass)
        WeightedCentroidX = sum(VectorOfSumX.*XIndices)/sum(VectorOfSumX);
        WeightedCentroidY = sum(VectorOfSumY.*YIndices)/sum(VectorOfSumY);
        
        %The stdev of the weighted centroids are calculated.  Note that we
        %multiply by VectorOfSum and then divide by the total sum of that
        %same vector because we are in essence applying a stdev calculation
        %where each of the elements has been weighted by a different
        %probability.
        StdevXFromCoM = sqrt(sum(VectorOfSumX.*(abs(XIndices-WeightedCentroidX).^2))/sum(VectorOfSumX));
        StdevYFromCoM = sqrt(sum(VectorOfSumY.*(abs(YIndices-WeightedCentroidY).^2))/sum(VectorOfSumY));

    end



end

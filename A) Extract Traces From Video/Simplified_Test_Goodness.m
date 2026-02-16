function [GoodParticle, CroppedImageProps, CurrParticleBox, OffsetFrom2DFit, Noise,...
    DidFitWork, CroppedParticleImageThresholded, SizeOfSquareAroundCurrParticle,...
    CurrParticleEccentricity, NewArea, WhyParticleFailed] =...
    Simplified_Test_Goodness(PrimaryImage,PropsOfCurrParticle,DataBits,...
    NormalizedPrimaryImageThresh, MinParticleSize, MaxEccentricityAllowed,ImageWidth,...
    ImageHeight, MaxParticleSize,BinaryPrimaryImage,Options)

%--------------------------------------------------------------------------
% Simplified_Test_Goodness.m
% 
% This function tests whether or not a particle is "good" and should be
% quantified.  The various tests are described below.  As part of the test,
% a gaussian fit is performed (this is done to ensure that the eccentricity
% and pixel values calculated for small particles are not incorrectly
% represented in samples where there are very big particles which dominate
% the original thresholding).
%
%--------------------------------------------------------------------------

%The outputs of the function are originally defined as NaN so that if the particle
%doesn't pass the test, the function doesn't crash because its outputs are
%not defined.  These values will not be used if the particle doesn't pass
%the tests below.
CroppedImageProps.Exists = 'n';
CurrParticleBox.Exists = 'n';
OffsetFrom2DFit = NaN;
Noise = NaN;
CroppedParticleImageThresholded = NaN;
DidFitWork = NaN;
CurrParticleEccentricity = NaN;
NewArea = NaN;
WhyParticleFailed = 'N/A';
GoodParticle = 'n';

%Values for the current particle are extracted from the PropsOfCurrParticle
%CurrParticleEccentricity = PropsOfCurrParticle.Eccentricity;
CurrParticleCentroid = (PropsOfCurrParticle.Centroid);
CurrParticleArea = PropsOfCurrParticle.Area; 
CurrParticlePixelValues = PropsOfCurrParticle.PixelValues;

%From those values, variables are defined which will be used below
CurrParticleXCenter = CurrParticleCentroid(1); 
CurrParticleYCenter = CurrParticleCentroid(2);
AreaForROI = max([PropsOfCurrParticle.Area, (Options.MinROISize)^2]);
    if ~isnan(Options.MaxROISize)
        if AreaForROI > (Options.MaxROISize)^2
            AreaForROI = (Options.MaxROISize)^2;
        end
    end
    
SizeOfSquareAroundCurrParticle = (sqrt(AreaForROI));  % was formerly (sqrt(AreaForROI)*2)
SaturationTest = length(find(2^DataBits - CurrParticlePixelValues <= 1));
    %The SaturationTest will be non-zero if any of the pixels are
    %saturated.

    %Coordinates are set up to define the area around the particle
    %of interest (i.e. the ROI)  The trigger "Exists" is changed to 'y'.
    CurrParticleBox.Left = max([round(CurrParticleXCenter) - round(SizeOfSquareAroundCurrParticle/2),...
        1]);
    CurrParticleBox.Right = min([round(CurrParticleXCenter) + round(SizeOfSquareAroundCurrParticle/2),...
        ImageWidth]);
    CurrParticleBox.Top = max([round(CurrParticleYCenter) - round(SizeOfSquareAroundCurrParticle/2),...
        1]);
    CurrParticleBox.Bottom = min([round(CurrParticleYCenter) + round(SizeOfSquareAroundCurrParticle/2),...
        ImageHeight]);
    
    CurrParticleBox.Exists = 'y';


%For the current particle, a border around the entire image is defined
%to make sure that the particle isn't too close to the side.
MinXBorder = SizeOfSquareAroundCurrParticle/2 + 1;
MaxXBorder = ImageWidth - SizeOfSquareAroundCurrParticle/2 - 1;
MinYBorder = SizeOfSquareAroundCurrParticle/2 + 1;
MaxYBorder = ImageHeight - SizeOfSquareAroundCurrParticle/2 - 1;



%------Test if the current particle being analyzed is "good"----------------
%There are a series of tests and unfortunately, they need to be applied
%sequentially.  First, we test if the particle isn't too close to the
%edge.  We also make sure that none of the pixels are saturated.
if (CurrParticleXCenter > MinXBorder &&...
    CurrParticleXCenter < MaxXBorder &&...
    CurrParticleYCenter > MinYBorder &&...
    CurrParticleYCenter < MaxYBorder &&...
    SaturationTest == 0 &&...
    CurrParticleArea < MaxParticleSize)

    %Then there are four more tests that need to be performed before the
    %particle can be considered "good".  1) determine that no
    %other particles lie within the region that will be used for the
    %gaussian fit.  2) make sure that the particle is not near the edge
    %of an area that was excluded (i.e. that pixels w/ value = zero will not be included in
    %the gaussian fit).  3) check the particle eccentricity 
    %To perform these tests, a few calculations need to be done. 
    % Note: this second test will never be a problem if
    % a patch mask was not defined.


    %A cropped image around the particle is created which will be used
    %in the gaussian fit.
    CurrParticleCroppedImage = PrimaryImage(CurrParticleBox.Top:CurrParticleBox.Bottom,...
        CurrParticleBox.Left:CurrParticleBox.Right);
    
    
    %A binary image is created to check if the region of analysis
    %around the current particle happens to overlap with more than 3
    %pixels of another particle (i.e. 3 connected pixels above the
    %threshold).
%     BinaryCroppedImage = im2bw(CurrParticleCroppedImage, NormalizedPrimaryImageThresh);
%     BinaryCroppedImage = bwareaopen(BinaryCroppedImage, 2, 4);
    BinaryCroppedImage = BinaryPrimaryImage(CurrParticleBox.Top:CurrParticleBox.Bottom,...
        CurrParticleBox.Left:CurrParticleBox.Right);
    CroppedImageComponents = bwconncomp(BinaryCroppedImage,4);
    CroppedImageProps = regionprops(CroppedImageComponents, 'Centroid', 'Eccentricity');
    NumOfParticlesInCroppedRegion = length(CroppedImageProps);

    %Pixels w/ value = zero in the ROI are detected (indicating that
    %the particle is near the edge of the patch).
    NumOfZeroPixels = length(find(CurrParticleCroppedImage == 0));

    %Now we test that there are no zero pixels and that there is only one
    %particle in the cropped region.
    if NumOfZeroPixels == 0 && NumOfParticlesInCroppedRegion == 1

        %Now we do a trick to capture both bright, big particles and small
        %dim particles on the same sample.  If the area of the particle after
        %the initial thresholding is big enough, then we examine its
        %eccentricity directly and determine the local background using a
        %gaussian fit.  If the fit works, then the particle is "good".
        if CurrParticleArea > 12.5
            CurrParticleEccentricity = PropsOfCurrParticle.Eccentricity;
            if CurrParticleEccentricity <= MaxEccentricityAllowed
                try
                    [OffsetFrom2DFit, Noise] = Particle_Gaussian_Fit(CurrParticleCroppedImage);
                    LocalBackgroundFromGaussFit = abs(OffsetFrom2DFit) + Noise;
                    NormalizedLocalThresh = LocalBackgroundFromGaussFit/2^DataBits;
                    if (NormalizedLocalThresh > 0 && NormalizedLocalThresh < 1)
                        DidFitWork = 1;
                    else
                        DidFitWork = 0;
                    end
                catch
                    %disp('Gaussian fit failed.  Particle Ignored.')
                    DidFitWork = 0;
                end
                
                if DidFitWork == 1
                    %The local background is applied to the cropped image
                    %as a threshold so that it will be quantified
                    %correctly (i.e. accounts for times when the initial
                    %threshold set happens to be higher than the local
                    %threshold which is found).
                    BinaryCroppedImage = im2bw(CurrParticleCroppedImage, NormalizedLocalThresh);
                    BinaryCroppedImage = bwareaopen(BinaryCroppedImage, MinParticleSize, 4);
                    CroppedImageComponents = bwconncomp(BinaryCroppedImage,4);
                    CroppedImageProps = regionprops(CroppedImageComponents, CurrParticleCroppedImage,...
                        'Eccentricity', 'PixelValues', 'Area');
                    NumOfParticlesInCroppedRegion = length(CroppedImageProps);
                    %If, after this local fit, there is only one particle in
                    %the cropped image, then analyze it.
                    if NumOfParticlesInCroppedRegion == 1
                        CurrParticleEccentricity = CroppedImageProps.Eccentricity;
                        NewArea = CroppedImageProps.Area;
                        CroppedParticleImageThresholded = CurrParticleCroppedImage;
                        CroppedParticleImageThresholded(BinaryCroppedImage == 0) = 0;
                        GoodParticle = 'y';
                    else
                        GoodParticle = 'n';
                        WhyParticleFailed = 'Too Many Regions, Normal Gauss';
                    end

                else
                    GoodParticle = 'n';
                    WhyParticleFailed = 'Fit failed';
                end
            else
                GoodParticle = 'n';
                WhyParticleFailed = 'Eccentricity Big Particle';
            end
        %If the particle is too small, then we perform a gaussian fit
        %directly on the raw data and use the local background determined
        %from the fit to calculate the eccentricity (and also to include
        %more pixels in the intensity calculation).
        else
            try
                [OffsetFrom2DFit, Noise] = Particle_Gaussian_Fit(CurrParticleCroppedImage);
                LocalBackgroundFromGaussFit = abs(OffsetFrom2DFit) + Noise;
                NormalizedLocalThresh = LocalBackgroundFromGaussFit/2^DataBits;
                DidFitWork = 1;
            catch
                if strcmp(Options.DisplayRejectionReasons,'y')
                    disp('Gaussian fit failed.  Particle Ignored.')
                end
                DidFitWork = 0;
            end



            if DidFitWork == 1
                %A binary image is created to determine the eccentricity of
                %the particle, now that a local threshold has been applied
                %(this is done to account for images that have bright
                %particles and dim particles on the same sample).
                BinaryCroppedImage = im2bw(CurrParticleCroppedImage, NormalizedLocalThresh);
                BinaryCroppedImage = bwareaopen(BinaryCroppedImage, MinParticleSize, 4);
                CroppedImageComponents = bwconncomp(BinaryCroppedImage,4);
                CroppedImageProps = regionprops(CroppedImageComponents, CurrParticleCroppedImage,...
                    'Eccentricity', 'PixelValues', 'Area');
                NumOfParticlesInCroppedRegion = length(CroppedImageProps);

                %If there is only one particle in the region (after the new thresholding), test the
                %eccentricity and determine particle goodness.
                if NumOfParticlesInCroppedRegion == 1 
                                   
                    CurrParticleEccentricity = CroppedImageProps.Eccentricity;
                    NewArea = CroppedImageProps.Area;
                    CroppedParticleImageThresholded = CurrParticleCroppedImage;
                    CroppedParticleImageThresholded(BinaryCroppedImage == 0) = 0;
                    
                    if CurrParticleEccentricity <= MaxEccentricityAllowed
                        GoodParticle = 'y';
                    else

                        GoodParticle = 'n';
                        WhyParticleFailed = 'Eccentricity Sm Particle';
                    end

                else
                    GoodParticle = 'n';
                    WhyParticleFailed = 'Too Many Regions, Small Gauss';
                end
            else
                GoodParticle = 'n';
                WhyParticleFailed = 'Fit failed';
            end
        end
    else
        GoodParticle = 'n';
        if NumOfZeroPixels ~= 0
            WhyParticleFailed = 'Zero Pixels';
        elseif NumOfParticlesInCroppedRegion ~= 1
            WhyParticleFailed = 'Too Many Regions, First Pass';
        end
    end

else
    GoodParticle = 'n';
    
    if SaturationTest ~= 0
        WhyParticleFailed = 'Saturation';
    elseif CurrParticleArea >= MaxParticleSize
        WhyParticleFailed = '>= Max Size';
    else
        WhyParticleFailed = 'Edge';
    end
end
%--------------------------------------------------------------------------
    
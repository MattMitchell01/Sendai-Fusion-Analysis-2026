function [GoodVesicle, CroppedImageProps, CurrVesBox, OffsetFrom2DFit, Noise,...
    DidFitWork, CroppedVesImageThresholded, SizeOfSquareAroundCurrVesicle,...
    CurrVesicleEccentricity, NewArea, WhyVesicleFailed] =...
    Test_Vesicle_Goodness(PrimaryImage,PropsOfCurrVes,DataBits,...
    NormalizedPrimaryImageThresh, MinVesSize, MaxEccentricityAllowed)

%--------------------------------------------------------------------------
% Test_Vesicle_Goodness.m
% Called from Quantify_Vesicles_N_Colors.m
% 
% This function tests whether or not a vesicle is "good" and should be
% quantified.  The various tests are described below.  As part of the test,
% a gaussian fit is performed (this is done to ensure that the eccentricity
% and pixel values calculated for small vesicles are not incorrectly
% represented in samples where there are very big vesicles which dominate
% the original thresholding).  Since the Gaussian fit is performed here, we
% send the output data from that (Offset and Noise) back to
% Quantify_Vesicles_N_Colors to be used in the quantification should the
% vesicle be found to be "good".
%
%--------------------------------------------------------------------------

%The outputs of the function are originally defined as NaN so that if the vesicle
%doesn't pass the test, the function doesn't crash because its outputs are
%not defined.  These values will not be used if the vesicle doesn't pass
%the tests below.
CroppedImageProps.Exists = 'n';
CurrVesBox.Exists = 'n';
OffsetFrom2DFit = NaN;
Noise = NaN;
CroppedVesImageThresholded = NaN;
DidFitWork = NaN;
CurrVesicleEccentricity = NaN;
NewArea = NaN;
WhyVesicleFailed = NaN;


%Values for the current vesicle are extracted from the PropsOfCurrVes
%CurrVesicleEccentricity = PropsOfCurrVes.Eccentricity;
CurrVesicleCentroid = (PropsOfCurrVes.Centroid);
CurrVesicleArea = PropsOfCurrVes.Area; 
CurrVesiclePixelValues = PropsOfCurrVes.PixelValues;

%From those values, variables are defined which will be used below
CurrVesicleXCenter = CurrVesicleCentroid(1); 
CurrVesicleYCenter = CurrVesicleCentroid(2);
AreaForROI = max([PropsOfCurrVes.Area, 12.25]); %The min square size is 7x7 pixels (even if the vesicle is super small)
if AreaForROI > 49 
    AreaForROI = 49; %The maximum square size will be 15x15 pixels.
end
SizeOfSquareAroundCurrVesicle = (sqrt(AreaForROI)*2);
SaturationTest = length(find(2^DataBits - CurrVesiclePixelValues <= 1));
    %The SaturationTest will be non-zero if any of the pixels are
    %saturated.


%For the current vesicle, a border around the entire image is defined
%to make sure that the vesicle isn't too close to the side.
MinXBorder = SizeOfSquareAroundCurrVesicle/2 + 1;
MaxXBorder = size(PrimaryImage,2) - SizeOfSquareAroundCurrVesicle/2 - 1;
MinYBorder = SizeOfSquareAroundCurrVesicle/2 + 1;
MaxYBorder = size(PrimaryImage,1) - SizeOfSquareAroundCurrVesicle/2 - 1;



%------Test if the current vesicle being analyzed is "good"----------------
%There are a series of tests and unfortunately, they need to be applied
%sequentially.  First, we test if the vesicle isn't too close to the
%edge.  We also make sure that none of the pixels are saturated.
if CurrVesicleXCenter > MinXBorder &&...
    CurrVesicleXCenter < MaxXBorder &&...
    CurrVesicleYCenter > MinYBorder &&...
    CurrVesicleYCenter < MaxYBorder &&...
    SaturationTest == 0

    %Then there are three more tests that need to be performed before the
    %vesicle can be considered "good".  1) determine that no
    %other vesicles lie within the region that will be used for the
    %gaussian fit.  2) make sure that the vesicle is not near the edge
    %of the patch (i.e. that pixels w/ value = zero will be included in
    %the gaussian fit).  3) check the vesicle eccentricity 
    %To perform these tests, a few calculations need to be done. 
    % Note: this second test will never be a problem if
    % a patch mask was not defined.

    %Coordinates are set up to define the area around the vesicle
    %of interest (i.e. the ROI)  The trigger "Exists" is changed to 'y'.
    CurrVesBox.Left = round(CurrVesicleXCenter) - round(SizeOfSquareAroundCurrVesicle/2);
    CurrVesBox.Right = round(CurrVesicleXCenter) + round(SizeOfSquareAroundCurrVesicle/2);
    CurrVesBox.Top = round(CurrVesicleYCenter) - round(SizeOfSquareAroundCurrVesicle/2);
    CurrVesBox.Bottom = round(CurrVesicleYCenter) + round(SizeOfSquareAroundCurrVesicle/2);
    CurrVesBox.Exists = 'y';

    %A cropped image around the vesicle is created which will be used
    %in the gaussian fit.
    CurrVesicleCroppedImage = PrimaryImage(CurrVesBox.Top:CurrVesBox.Bottom,...
        CurrVesBox.Left:CurrVesBox.Right);
    %figure(51)
    %imshow(CurrVesicleCroppedImage, [130, 400],'InitialMagnification', 'fit');
    %pause(0.02);
    
    %A binary image is created to check if the region of analysis
    %around the current vesicle happens to overlap with more than 3
    %pixels of another vesicle (i.e. 3 connected pixels above the
    %threshold).
    BinaryCroppedImage = im2bw(CurrVesicleCroppedImage, NormalizedPrimaryImageThresh);
    BinaryCroppedImage = bwareaopen(BinaryCroppedImage, 2, 4);
    CroppedImageComponents = bwconncomp(BinaryCroppedImage,4);
    CroppedImageProps = regionprops(CroppedImageComponents, 'Centroid', 'Eccentricity');
    NumOfVesiclesInCroppedRegion = length(CroppedImageProps);

    %Pixels w/ value = zero in the ROI are detected (indicating that
    %the vesicle is near the edge of the patch).
    NumOfZeroPixels = length(find(CurrVesicleCroppedImage == 0));

    %Now we test that there are no zero pixels and that there is only one
    %vesicle in the cropped region.
    if NumOfZeroPixels == 0 && NumOfVesiclesInCroppedRegion == 1

        %Now we do a trick to capture both bright, big vesicles and small
        %dim vesicles on the same sample.  If the area of the vesicle after
        %the initial thresholding is big enough, then we examine its
        %eccentricity directly and determine the local background using a
        %gaussian fit.  If the fit works, then the vesicle is "good".
        if CurrVesicleArea > 12.5
            CurrVesicleEccentricity = PropsOfCurrVes.Eccentricity;
            if CurrVesicleEccentricity <= MaxEccentricityAllowed
                try
                    [OffsetFrom2DFit, Noise] = Vesicle_Gaussian_Fit(CurrVesicleCroppedImage);
                    LocalBackgroundFromGaussFit = abs(OffsetFrom2DFit) + Noise;
                    NormalizedLocalThresh = LocalBackgroundFromGaussFit/2^DataBits;
                    if (NormalizedLocalThresh > 0 && NormalizedLocalThresh < 1)
                        DidFitWork = 1;
                    else
                        DidFitWork = 0;
                    end
                catch
                    disp('Gaussian fit failed.  Vesicle Ignored.')
                    DidFitWork = 0;
                end
                
                if DidFitWork == 1
                    %The local background is applied to the cropped image
                    %as a threshold so that it will be quantified
                    %correctly (i.e. accounts for times when the initial
                    %threshold set happens to be higher than the local
                    %threshold which is found).
                    BinaryCroppedImage = im2bw(CurrVesicleCroppedImage, NormalizedLocalThresh);
                    BinaryCroppedImage = bwareaopen(BinaryCroppedImage, MinVesSize, 4);
                    CroppedImageComponents = bwconncomp(BinaryCroppedImage,4);
                    CroppedImageProps = regionprops(CroppedImageComponents, CurrVesicleCroppedImage,...
                        'Eccentricity', 'PixelValues', 'Area');
                    NumOfVesiclesInCroppedRegion = length(CroppedImageProps);
                    CurrVesicleEccentricity = CroppedImageProps.Eccentricity;
                    NewArea = CroppedImageProps.Area;
                    CroppedVesImageThresholded = CurrVesicleCroppedImage;
                    CroppedVesImageThresholded(BinaryCroppedImage == 0) = 0;
                    
                    %If, after this local fit, there is only one vesicle in
                    %the cropped image, then analyze it.
                    if NumOfVesiclesInCroppedRegion == 1
                        GoodVesicle = 'y';
                    else
                        GoodVesicle = 'n';
                        WhyVesicleFailed = 'Too Many Regions';
                    end

                else
                    GoodVesicle = 'n';
                    WhyVesicleFailed = 'Fit failed';
                end
            else
                GoodVesicle = 'n';
                WhyVesicleFailed = 'Eccentricity Big Ves';
            end
        %If the vesicle is too small, then we perform a gaussian fit
        %directly on the raw data and use the local background determined
        %from the fit to calculate the eccentricity (and also to include
        %more pixels in the intensity calculation).
        else
            try
                [OffsetFrom2DFit, Noise] = Vesicle_Gaussian_Fit(CurrVesicleCroppedImage);
                LocalBackgroundFromGaussFit = abs(OffsetFrom2DFit) + Noise;
                NormalizedLocalThresh = LocalBackgroundFromGaussFit/2^DataBits;
                DidFitWork = 1;
            catch
                disp('Gaussian fit failed.  Vesicle Ignored.')
                DidFitWork = 0;
            end



            if DidFitWork == 1
                %A binary image is created to determine the eccentricity of
                %the vesicle, now that a local threshold has been applied
                %(this is done to account for images that have bright
                %vesicles and dim vesicles on the same sample).
                BinaryCroppedImage = im2bw(CurrVesicleCroppedImage, NormalizedLocalThresh);
                BinaryCroppedImage = bwareaopen(BinaryCroppedImage, MinVesSize, 4);
                CroppedImageComponents = bwconncomp(BinaryCroppedImage,4);
                CroppedImageProps = regionprops(CroppedImageComponents, CurrVesicleCroppedImage,...
                    'Eccentricity', 'PixelValues', 'Area');
                NumOfVesiclesInCroppedRegion = length(CroppedImageProps);

                %If there is only one vesicle in the region (after the new thresholding), test the
                %eccentricity and determine vesicle goodness.
                if NumOfVesiclesInCroppedRegion == 1 
                                   
                    CurrVesicleEccentricity = CroppedImageProps.Eccentricity;
                    NewArea = CroppedImageProps.Area;
                    CroppedVesImageThresholded = CurrVesicleCroppedImage;
                    CroppedVesImageThresholded(BinaryCroppedImage == 0) = 0;
                    
                    if CurrVesicleEccentricity <= MaxEccentricityAllowed
                        GoodVesicle = 'y';
                    else

                        GoodVesicle = 'n';
                        WhyVesicleFailed = 'Eccentricity Sm Ves';
                    end

                else
                    GoodVesicle = 'n';
                    WhyVesicleFailed = 'Too Many Regions';
                end
            else
                GoodVesicle = 'n';
                WhyVesicleFailed = 'Fit failed';
            end
        end
    else
        GoodVesicle = 'n';
        if NumOfZeroPixels ~= 0
            WhyVesicleFailed = 'Zero Pixels';
        elseif NumOfVesiclesInCroppedRegion ~= 1
            WhyVesicleFailed = 'Too Many Regions';
        end
    end

else
    GoodVesicle = 'n';
    if SaturationTest ~= 0
        WhyVesicleFailed = 'Saturation';
    else
        WhyVesicleFailed = 'Edge';
    end
end
%--------------------------------------------------------------------------
    
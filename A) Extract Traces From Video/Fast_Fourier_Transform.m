function [Offsets] = Fast_Fourier_Transform(VideoMatrix, ImageWidth, ImageHeight, NumFrames, Options)

%====================================================================================
% FAST-FOURIER-TRANSFORM XY SHIFT ESTIMATION (PHASE CROSS-CORRELATION)
%====================================================================================
%
% ALGORITHM FLOW (PER FRAME)
% --------------------------
% 1. Extract the same region of interest (ROI) from:
%       - the reference image (frame 1 of VideoMatrix = the user-defined start frame)
%       - the current image
%
% 2. Z-score normalize each ROI:
%       - subtract the mean   → removes overall brightness differences
%       - divide by std dev   → removes contrast / variance differences
%    This ensures the correlation depends only on spatial structure, not
%    illumination power or exposure.
%
% 3. Compute the FFT of both ROIs.
%
% 4. Perform phase-only cross-correlation in Fourier space:
%       R = Ref_ROI_FFT_Matrix .* conj(Curr_ROI_FFT_Matrix) ./ |Ref_ROI_FFT_Matrix .* conj(Curr_ROI_FFT_Matrix)|
%    Keeping only phase makes the method robust to intensity scaling.
%
% 5. Inverse FFT → correlation surface.
%    The location of the peak gives the relative X–Y shift.
%
% 6. Refine the peak location to subpixel accuracy using a local quadratic
%    (parabolic) fit in X and Y.
% 
% 7. Multiply the FFT output by -1. Why --> FFT calculated the corrections
%     that make the two frames align. Logic for shifting the ROI's in 
%     Find_And_Analyze_Particles needs the physical drift number (not correction).  
%
%
% OUTPUT
% ------
% Offsets.x and Offsets.y contain the X and Y shifts of each frame relative
% to the reference. The first entry is [0, 0] by definition.
%

    % By Matthew D. Mitchell, Rawle Lab, Williams College (Jan 2026)
%====================================================================================


    disp('   Determining XY shift using FFT cross-correlation...')

    % ROI defined as percentages of image size
    Region_Percents = Options.FFTRegion;

    Region_Left_X  = round(ImageWidth  * Region_Percents(1));
    Region_Right_X = round(ImageWidth  * Region_Percents(2));
    Region_Top_Y   = round(ImageHeight * Region_Percents(1));
    Region_BOT_Y   = round(ImageHeight * Region_Percents(2));

    % Extract reference ROI (frame 1 of VideoMatrix is the user "starting image")
    Ref_ROI = double(VideoMatrix( ...
        Region_Top_Y:Region_BOT_Y, ...
        Region_Left_X:Region_Right_X, ...
        1));

    % Z-score Normalization 
    Ref_ROI  = Ref_ROI  - mean(Ref_ROI(:));
    Ref_ROI  = Ref_ROI  / std(Ref_ROI(:));
    
    % FFT of reference
    Ref_ROI_FFT_Matrix = fft2(Ref_ROI);

    % Initialize output index
    k = 1;

    for frame = 1:NumFrames

        if frame == 1
            Offsets.x(k) = 0;
            Offsets.y(k) = 0;
            k = k + 1;
            continue
        end

        % Extract current ROI
        Curr_ROI = double(VideoMatrix( ...
            Region_Top_Y:Region_BOT_Y, ...
            Region_Left_X:Region_Right_X, ...
            frame));

        % Z-score normalization
        Curr_ROI = Curr_ROI - mean(Curr_ROI(:));
        Curr_ROI = Curr_ROI / std(Curr_ROI(:));

        % FFT of current frame ROI
        Curr_ROI_FFT_Matrix = fft2(Curr_ROI);

        CrossPowerSpectrumMatrix = Ref_ROI_FFT_Matrix .* conj(Curr_ROI_FFT_Matrix);
        CrossPowerSpectrumMatrix = CrossPowerSpectrumMatrix ./ (abs(CrossPowerSpectrumMatrix) + eps);

        % Correlation surface in spatial domain
        CorrelationSurfaceMatrix = real(ifft2(CrossPowerSpectrumMatrix));
        CorrelationSurfaceMatrix = fftshift(CorrelationSurfaceMatrix);

        % Find integer peak of correlation
        [~, PeakLinearIndex_int] = max(CorrelationSurfaceMatrix(:));
        [peakY, peakX] = ind2sub(size(CorrelationSurfaceMatrix), PeakLinearIndex_int);

        % Convert peak location to signed pixel shift (integer part)
        centerX = floor(size(CorrelationSurfaceMatrix,2)/2) + 1;
        centerY = floor(size(CorrelationSurfaceMatrix,1)/2) + 1;

        IntegerShiftX = peakX - centerX;
        IntegerShiftY = peakY - centerY;

        % -------------------------------
        % SUBPIXEL REFINEMENT (quadratic fit)
        % -------------------------------
        % Subpixel refinement via local quadratic interpolation:
        % Use the correlation peak and its immediate neighbors in X and Y
        % to estimate the true peak position between pixels.
        %
        % delta = (f_left - f_right) / (2*(f_left - 2*f_center + f_right))
        %
        % Applied only when the peak is not on the edge of the correlation map.
        SubpixelShiftX = 0;
        SubpixelShiftY = 0;

        if peakX > 1 && peakX < size(CorrelationSurfaceMatrix,2) && ...
           peakY > 1 && peakY < size(CorrelationSurfaceMatrix,1)

            PeakValue = CorrelationSurfaceMatrix(peakY, peakX);

            PeakXMinus1 = CorrelationSurfaceMatrix(peakY, peakX-1);
            PeakXPlus1  = CorrelationSurfaceMatrix(peakY, peakX+1);

            PeakYMinus1 = CorrelationSurfaceMatrix(peakY-1, peakX);
            PeakYPlus1  = CorrelationSurfaceMatrix(peakY+1, peakX);

            DenomX = (PeakXMinus1 - 2*PeakValue + PeakXPlus1);
            DenomY = (PeakYMinus1 - 2*PeakValue + PeakYPlus1);

            if abs(DenomX) > eps
                SubpixelShiftX = 0.5 * (PeakXMinus1 - PeakXPlus1) / DenomX;
            end
            if abs(DenomY) > eps
                SubpixelShiftY = 0.5 * (PeakYMinus1 - PeakYPlus1) / DenomY;
            end
        end

        % FFT calculated the correction necessary to make the images align.
        % But logic later in the pipeline uses the physical drift so we
        % must multiply the Offsets by -1. 
        Offsets.x(k) = -1*(IntegerShiftX + SubpixelShiftX); 
        Offsets.y(k) = -1*(IntegerShiftY + SubpixelShiftY);

        k = k + 1;
    end

    disp('   Done determining XY shift using FFT.')
end
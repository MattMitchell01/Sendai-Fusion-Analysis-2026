function smoothed = Calculate_Sliding_Window(vector, windowWidth, method, raiseToZero, clipWidth)
%
% ------------------------------------------------------------------------
% Written by Matthew D. Mitchell, Rawle Lab, Williams College, 2026.
% New script developed as part of the 2026 overhaul and expansion of the
% Sendai fusion analysis pipeline originally created by Bob Rawle
% (Kasson Lab, University of Virginia, 2016; Rawle et al., Disentangling
% Viral Membrane Fusion from Receptor Binding Using Synthetic DNA-Lipid
% Conjugates, Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% pipeline subsequently updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
% Calculate_Sliding_Window  Sliding window smoothing including the center point.
%                Window is always symmetric: if fewer points exist on one
%                side, the other side is shrunk to match.
%
%                The first/last clipWidth frames of `vector` are used as
%                real smoothing context (so the window can grow toward
%                windowWidth using genuine data instead of an artificially
%                shrunk boundary window) but are dropped from the returned
%                `smoothed` output — they are never meant to be scored.
%
% Inputs:
%   vector      - Input trace (numeric vector)
%   windowWidth - Full window width (integer); uses floor(width/2) on each side
%   method      - 'median' or 'mean'
%   raiseToZero - 1 = shift output up so minimum equals zero; 0 = no shift
%   clipWidth   - Number of frames to drop from each end of the OUTPUT after
%                 smoothing. Their real values are still used as window
%                 context for the frames that remain.
%
% Output:
%   smoothed    - Smoothed trace, length numel(vector) - 2*clipWidth

halfWidth = floor(windowWidth / 2);
n         = numel(vector);

if n - 2*clipWidth < 1
    error('Calculate_Sliding_Window: clipWidth (%d) leaves no frames after trimming a trace of length %d.', clipWidth, n);
end

smoothed  = zeros(size(vector));

for i = 1:n
    hw        = min([halfWidth, i-1, n-i]);   % symmetric: shrink both sides equally
    neighbors = vector(i-hw : i+hw);          % includes center point i

    switch method
        case 'median'
            smoothed(i) = median(neighbors, 'omitnan');
        case 'mean'
            smoothed(i) = mean(neighbors, 'omitnan');
        otherwise
            error('Calculate_Sliding_Window: method must be ''median'' or ''mean''.');
    end
end

smoothed = smoothed(clipWidth+1 : n-clipWidth);

if raiseToZero
    minVal = min(smoothed);
    if minVal < 0
        smoothed = smoothed + abs(minVal);
    end
end

end

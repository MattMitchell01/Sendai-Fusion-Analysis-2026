function [CurrTrace_Corrected] = Run_Med(CurrTrace_Corrected,Options)
%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
% Run_Med  Sliding-window running median, mirroring Part B's own
% Calculate_Sliding_Window.m method: the window shrinks SYMMETRICALLY as
% either edge is approached (hw = min(halfWidth, n-1, TraceLength-n)), not
% a fixed-width window that expands asymmetrically toward the far edge.
%
% At index 1 (and the last index) hw=0, so the window is just that one
% point -- no smoothing at all there. Moving one step in from either edge
% grows the window by one frame on each side (back AND forward, kept
% equal) until RunMedHalfLength is reached on both sides in the interior.

    RunMedHalfLength = Options.RunMedHalfLength;
    TraceLength = length(CurrTrace_Corrected);
    OldTrace = CurrTrace_Corrected;
    NewTrace = zeros(size(CurrTrace_Corrected));

    for n = 1:TraceLength
        hw = min([RunMedHalfLength, n-1, TraceLength-n]);
        NewTrace(n) = median(OldTrace(n-hw : n+hw), 'omitnan');
    end

    CurrTrace_Corrected = NewTrace;

end





function Display_Reviewed_Viruses(FindingImage, CombinedAnalyzedTraceData, MinImageShow, MaxImageShow, FrameNumToFindParticles)
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
% Display_Reviewed_Viruses  Plots the finding image with a box drawn
%                            around every trace, colored by final
%                            (reviewer-corrected) Designation.
%
% Adapted from Part B's Display_All_Analyzed_Viruses.m
% ('../Part (B)/_Dev_B) Lipid Mixing Trace Analysis/Display_All_Analyzed_Viruses.m'),
% same signature/logic, run at the end of a completed Part C review
% session against DataToSave.CombinedAnalyzedTraceData (post-correction),
% not the raw algorithm output. Renamed (not just reused) so it's never
% confused with the Part B original when both are on a MATLAB path at once.
%
% Inputs:
%   FindingImage              - Image matrix, DataToSave.OtherDataToSave.FindingImage
%   CombinedAnalyzedTraceData - Struct array, DataToSave.CombinedAnalyzedTraceData
%   MinImageShow              - Lower intensity display bound
%   MaxImageShow              - Upper intensity display bound
%   FrameNumToFindParticles   - Global frame number of the finding image, for the title

figure('Name', 'Reviewed Viruses', 'NumberTitle', 'off');
imshow(FindingImage, [MinImageShow, MaxImageShow], 'InitialMagnification', 'fit', 'Border', 'tight');
title(sprintf('Finding Image, Fr = %d  (%d traces, reviewed)', FrameNumToFindParticles, numel(CombinedAnalyzedTraceData)));
hold on;

for i = 1:numel(CombinedAnalyzedTraceData)

    switch CombinedAnalyzedTraceData(i).Designation
        case '1 Fuse'
            lineColor = 'g-';
        case '2 Fuse'
            lineColor = 'b-';
        case 'Unbound'
            lineColor = 'y-';
        case 'No Fusion'
            lineColor = 'r-';
        otherwise   % 'Other'
            lineColor = 'm-';
    end

    if ~isfield(CombinedAnalyzedTraceData, 'BoxAroundVirus') || isempty(CombinedAnalyzedTraceData(i).BoxAroundVirus)
        continue;
    end

    box  = CombinedAnalyzedTraceData(i).BoxAroundVirus;
    boxX = [box.Left, box.Right, box.Right, box.Left,  box.Left];
    boxY = [box.Top,  box.Top,   box.Bottom, box.Bottom, box.Top];
    plot(boxX, boxY, lineColor);

end

hold off;

end

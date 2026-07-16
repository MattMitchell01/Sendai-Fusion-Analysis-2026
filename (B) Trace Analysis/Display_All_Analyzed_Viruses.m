function Display_All_Analyzed_Viruses(FindingImage, AnalyzedTraceData, MinImageShow, MaxImageShow, FrameNumToFindParticles)
%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
% Display_All_Analyzed_Viruses  Plots the finding image with a box drawn
%                                around every analyzed trace, colored by
%                                final Designation.
%
% Adapted from Bob's original of the same name
% ('Fusion Analysis (26.0)/B) Lipid Mixing Trace Analysis/Display_All_Analyzed_Viruses.m'),
% which drove a live GUI figure via FigureHandles.ImageWindow and read each
% trace's box from AnalyzedTraceData(i).BoxCoords (an Nx2 closed polygon).
% This pipeline has no GUI figure-handle plumbing and does not produce
% BoxCoords, so this version opens its own figure and builds the polygon
% from AnalyzedTraceData(i).BoxAroundVirus (.Top/.Bottom/.Left/.Right --
% Part A's Find_And_Analyze_Particles.m field, passed through unchanged).
%
% Inputs:
%   FindingImage          - Image matrix, raw.OtherDataToSave.FindingImage
%   AnalyzedTraceData     - Struct array (DataToSave.CombinedAnalyzedTraceData)
%   MinImageShow          - Lower intensity display bound (raw.OtherDataToSave.Options.MinImageShow)
%   MaxImageShow          - Upper intensity display bound (raw.OtherDataToSave.Options.MaxImageShow)
%   FrameNumToFindParticles - Global frame number of the finding image, for the title

figure('Name', 'Analyzed Viruses', 'NumberTitle', 'off');
imshow(FindingImage, [MinImageShow, MaxImageShow], 'InitialMagnification', 'fit', 'Border', 'tight');
title(sprintf('Finding Image, Fr = %d  (%d traces)', FrameNumToFindParticles, numel(AnalyzedTraceData)));
hold on;

for i = 1:numel(AnalyzedTraceData)

    switch AnalyzedTraceData(i).Designation
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

    if ~isfield(AnalyzedTraceData, 'BoxAroundVirus') || isempty(AnalyzedTraceData(i).BoxAroundVirus)
        continue;
    end

    box  = AnalyzedTraceData(i).BoxAroundVirus;
    boxX = [box.Left, box.Right, box.Right, box.Left,  box.Left];
    boxY = [box.Top,  box.Top,   box.Bottom, box.Bottom, box.Top];
    plot(boxX, boxY, lineColor);

end

hold off;

end

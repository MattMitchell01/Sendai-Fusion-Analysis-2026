function Label = Get_H4_Failed_Test_Label(CurrentVirusData)
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
% Get_H4_Failed_Test_Label  Human-readable "which test(s) failed" label for
% an H4 trace (No Fusion, exactly one scoring test narrowly failed).
%
% Reads Review_PriorityData.SeekFusionH4/.SeekUnboundsH4 -- Part B's
% Assign_Review_Priority.m sets .Rule='H4' and these two sub-structs
% together, unconditionally, in the same code block, so both are
% guaranteed present whenever Rule=='H4'. Each carries a .Triggered flag
% plus per-test frame lists (empty if that specific test didn't fail);
% exactly which of the six ScoreNNFrames fields is non-empty identifies
% the failed test -- see the Part B Seek_Fusion/Seek_Unbounds
% Quick Reference tables for the score-to-test mapping.

    PriorityData = CurrentVirusData.Review_PriorityData;
    Reasons = {};

    if isfield(PriorityData, 'SeekFusionH4') && PriorityData.SeekFusionH4.Triggered
        sf = PriorityData.SeekFusionH4;
        if ~isempty(sf.Score29Frames)
            Reasons{end+1} = 'Fusion T2 (fractional rise)';
        end
        if ~isempty(sf.Score27Frames)
            Reasons{end+1} = 'Fusion T3 (sustained elevation)';
        end
        if ~isempty(sf.Score23Frames)
            Reasons{end+1} = 'Fusion T4 (pre-rise level)';
        end
    end

    if isfield(PriorityData, 'SeekUnboundsH4') && PriorityData.SeekUnboundsH4.Triggered
        ub = PriorityData.SeekUnboundsH4;
        if ~isempty(ub.Score13Frames)
            Reasons{end+1} = 'Unbound T2 (fractional drop)';
        end
        if ~isempty(ub.Score11Frames)
            Reasons{end+1} = 'Unbound T3 (sustained floor)';
        end
        if ~isempty(ub.Score7Frames)
            Reasons{end+1} = 'Unbound T4 (pre-drop elevation)';
        end
    end

    if isempty(Reasons)
        % Defensive fallback -- shouldn't happen for a real H4 trace, since
        % Rule=='H4' implies at least one of SeekFusionH4/SeekUnboundsH4
        % triggered with at least one non-empty ScoreNNFrames field.
        Label = 'Failed test: unknown';
    else
        Label = ['Failed: ' strjoin(Reasons, ', ')];
    end
end

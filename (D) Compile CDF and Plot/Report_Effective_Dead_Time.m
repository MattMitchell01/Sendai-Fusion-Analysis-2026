function Report_Effective_Dead_Time(EffectiveDeadTimes, SetDeadTime)
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
% Report_Effective_Dead_Time  Compares each input file's true "effective
% dead time" (OtherDataToSave.EffectiveDeadTime -- the real-world time, in
% minutes, before which no frame was actually scored by Part B) against
% Options.DeadTime, the lag time this analysis assumes. Warns if any file's
% true effective dead time is longer than what's assumed, since wait times
% in that window are not reliable.
%
% EffectiveDeadTimes - vector, one entry per input file. NaN for files
%                      missing OtherDataToSave or OtherDataToSave.EffectiveDeadTime
%                      (old-format Part A/B output).
% SetDeadTime        - Options.DeadTime, minutes.

    NumberOfFiles = numel(EffectiveDeadTimes);

    disp('   Effective dead time check:')

    if all(isnan(EffectiveDeadTimes))
        disp('      !!!! Could not check -- none of the selected file(s) have an OtherDataToSave.EffectiveDeadTime field (old-format file(s)).')
        return
    end

    if any(isnan(EffectiveDeadTimes))
        NumberMissing = sum(isnan(EffectiveDeadTimes));
        disp(strcat("      ***", num2str(NumberMissing), " of ", num2str(NumberOfFiles), ...
            " file(s) are missing OtherDataToSave.EffectiveDeadTime (old-format) and were excluded from this check."))
    end

    MaxEffectiveDeadTime = max(EffectiveDeadTimes, [], 'omitnan');

    disp(strcat("      Longest effective dead time across ", num2str(NumberOfFiles), ...
        " file(s) = ", num2str(MaxEffectiveDeadTime), " min"))
    disp(strcat("      Options.DeadTime (set by user)               = ", num2str(SetDeadTime), " min"))

    if MaxEffectiveDeadTime > SetDeadTime
        disp('      !!!! WARNING: The longest effective dead time exceeds Options.DeadTime.')
        disp('      !!!! Some events assumed to fall after the dead time may not actually have been scored yet.')
        disp('      !!!! Reconsider Options.DeadTime / Options.TimeCutoffLow in Setup_Options_Compile_CDF.m.')
    end

end

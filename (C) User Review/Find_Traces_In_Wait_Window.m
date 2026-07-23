function FilteredIndices = Find_Traces_In_Wait_Window(CombinedAnalyzedTraceData, LowMinutes, HighMinutes)
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
% Find_Traces_In_Wait_Window  Global indices of every '1 Fuse' trace whose
% BindtoFusionTime falls in [LowMinutes, HighMinutes] (inclusive) AND is
% currently trustworthy enough to actually be plotted by Part D's CDF.
%
% "Trustworthy" mirrors Extract_Data.m's own filter EXACTLY --
% (D) Compile CDF and Plot/Extract_Data.m, lines ~43-45 -- so a trace only
% comes back here if it's one of the ones genuinely responsible for a
% jump/artifact at this point in the real CDF, not a trace merely flagged
% '1 Fuse' but not yet corrected (ChangedByUser=='Incorrect Designation-Not
% Changed'), which doesn't contribute a wait time to the CDF at all.
% Keep this filter in sync with Extract_Data.m if that one ever changes.

    NumTraces = numel(CombinedAnalyzedTraceData);
    FilteredIndices = [];

    for k = 1:NumTraces
        Trace = CombinedAnalyzedTraceData(k);

        if ~strcmp(Trace.Designation, '1 Fuse')
            continue
        end
        if ~(strcmp(Trace.ChangedByUser, 'Not analyzed') || ...
                strcmp(Trace.ChangedByUser, 'Reviewed, Fuse frame chosen by user'))
            continue
        end

        WaitTime = Trace.FusionData.BindtoFusionTime(1);
        if WaitTime >= LowMinutes && WaitTime <= HighMinutes
            FilteredIndices(end+1) = k; %#ok<AGROW>
        end
    end
end

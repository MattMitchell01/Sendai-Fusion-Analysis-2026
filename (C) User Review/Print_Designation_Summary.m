function Print_Designation_Summary(CombinedAnalyzedTraceData)
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
% Print_Designation_Summary  Console tally of final (reviewer-corrected)
% Designation counts/percentages, printed at the end of a completed Part C
% review session. Order matches the user's requested category order:
% 1 Fuse, 2 Fuse, Unbound, No Fusion, Other.

Total = numel(CombinedAnalyzedTraceData);
Designations = {CombinedAnalyzedTraceData.Designation};

Labels = {'1 Fuse', '2 Fuse', 'Unbound', 'No Fusion', 'Other'};

for k = 1:numel(Labels)
    Count = sum(strcmp(Designations, Labels{k}));
    if Total > 0
        Pct = 100 * Count / Total;
    else
        Pct = 0;
    end
    fprintf('   %-10s %5d  (%5.1f%%)\n', Labels{k}, Count, Pct);
end

fprintf('   %-10s %5d\n', 'Total', Total);

end

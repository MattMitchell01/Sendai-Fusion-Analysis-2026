function GridConfig = Resolve_Grid_Config(Tier, SubgroupName, Options)
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
% Resolve_Grid_Config  Which Options grid struct applies to a given
% Tier/SubgroupName -- the single source of truth for this mapping, shared
% by Build_Review_Queue.m (first-run preprocessing) and
% Refresh_Segment_Grids (Start_User_Review.m, run on every resume) so the
% two can never drift out of sync with each other.

    switch Tier
        case 'Low'
            GridConfig = Options.Low;
        case 'Medium'
            GridConfig = Options.Medium;
        case 'High'
            if isfield(Options, 'HighOverrides') && isfield(Options.HighOverrides, SubgroupName)
                GridConfig = Options.HighOverrides.(SubgroupName);
            else
                GridConfig = Options.High;
            end
        otherwise
            GridConfig = Options.Low;   % 'Uncategorized' -- matches Build_Review_Queue.m's own fallback
    end
end

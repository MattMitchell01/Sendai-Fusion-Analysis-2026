function designation = Classify_Trace(ubClusterCount, sfClusterCount, landingData)
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
% Classify_Trace  Assigns a final designation from Seek_Unbounds and
%                 Seek_Fusion scores and Detect_Landing output.
%
% Inputs:
%   ubClusterCount - integer; number of unbound clusters (score >= 15)
%   sfClusterCount - integer; number of fusion clusters (score == 31)
%   landingData    - struct from Detect_Landing (may be empty struct)
%
% Output:
%   designation - one of: '1 Fuse', '2 Fuse', 'Unbound', 'No Fusion', 'Other'
%
% Classification rules (evaluated in priority order):
%   1. Landing detected                              → 'Other'
%   2. 2+ unbound clusters                          → 'Other'
%   3. Fusion cluster(s) present AND unbound         → 'Other'
%   4. Unbound only (no fusion clusters)             → 'Unbound'
%   5. Exactly 1 fusion cluster (no unbound)         → '1 Fuse'
%   6. Exactly 2 fusion clusters (no unbound)        → '2 Fuse'
%   7. 3+ fusion clusters (no unbound)               → 'Other'
%   8. Nothing found                                 → 'No Fusion'

landingDetected = isfield(landingData, 'LandingDetected') && landingData.LandingDetected;

if landingDetected
    designation = 'Other';
elseif ubClusterCount > 1
    designation = 'Other';
elseif sfClusterCount >= 1 && ubClusterCount >= 1
    designation = 'Other';
elseif ubClusterCount == 1
    designation = 'Unbound';
elseif sfClusterCount == 1
    designation = '1 Fuse';
elseif sfClusterCount == 2
    designation = '2 Fuse';
elseif sfClusterCount > 2
    designation = 'Other';
else
    designation = 'No Fusion';
end

end

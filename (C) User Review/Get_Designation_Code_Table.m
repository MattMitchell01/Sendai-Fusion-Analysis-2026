function CodeTable = Get_Designation_Code_Table()
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
% Get_Designation_Code_Table  The shared designation-code table used by
% Correct_Designations.m (top-level PlotNumber.Code prompt, keyed by the
% padded MilliCode) and Apply_Designation_Code.m / Parse_Designation_Code.m
% (mid-picker redirect prompt, ALSO keyed by MilliCode -- computed from a
% decimal-only input like ".11", the exact same fractional convention as
% the main prompt, just without a PlotNumber in front). DecimalCode is a
% pure display convenience (the digits a reviewer types after the leading
% dot) -- e.g. DecimalCode 11 == MilliCode 110 == typed ".11" either place,
% so there's only one code convention in the whole tool, never two to mix
% up. The repick code for each designation is just its own digit doubled
% (.1 -> .11, .2 -> .22, .3 -> .33) -- deliberately, so the repick/fresh
% pairing is easy to remember instead of each needing its own arbitrary
% suffix. See "Reference Guide For Prompt Codes.rtf":
%   MilliCode  DecimalCode
%   0          0     No Fusion
%   100        1     1 Fuse             (1 click: fuse frame)
%   110        11    1 Fuse, fix pick   (1 click: fuse frame)
%   200        2     2 Fuse             (2 clicks: fuse frames)
%   220        22    2 Fuse, fix picks  (2 clicks: fuse frames)
%   300        3     Unbound            (1 click: unbind frame -- mandatory, not optional)
%   330        33    Unbound, fix pick  (1 click: unbind frame)
%   900        9     Other

    CodeTable = struct('MilliCode', {0,100,110,200,220,300,330,900}, ...
        'DecimalCode', {0,1,11,2,22,3,33,9}, ...
        'Designation', {'No Fusion','1 Fuse','1 Fuse','2 Fuse','2 Fuse','Unbound','Unbound','Other'}, ...
        'Mode', {'new','new','repick','new','repick','new','repick','new'}, ...
        'NumClicks', {0,1,1,2,2,1,1,0}, ...
        'ClickTarget', {'none','Fuse','Fuse','Fuse','Fuse','Unbound','Unbound','none'});
end

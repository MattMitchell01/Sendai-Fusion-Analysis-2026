function [MilliCode, Valid] = Parse_Designation_Code(InputStr)
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
% Parse_Designation_Code  Parses a reviewer-typed decimal-only code (e.g.
% ".11", no integer part) into the padded MilliCode Apply_Designation_Code.m
% expects. Same fractional convention as the main PlotNumber.Code round
% prompt (Correct_Designations.m): the input's fractional part times 1000,
% rounded -- so ".1" is MilliCode 100, ".11" is MilliCode 110, etc. This is
% deliberately decimal-only (a leading "." is required, bare integers like
% the old "11" are rejected) so the mid-picker redirect prompt never mixes
% a second code convention in with the main prompt's.
%
% Valid=false for blank/non-numeric/no-decimal-point/negative/unrecognized
% input -- callers should check Valid before trusting MilliCode (NaN when
% Valid is false).

    InputStr = strtrim(InputStr);

    if isempty(InputStr) || ~contains(InputStr, '.')
        MilliCode = NaN;
        Valid = false;
        return
    end

    Value = str2double(InputStr);
    if isnan(Value) || Value < 0
        MilliCode = NaN;
        Valid = false;
        return
    end

    MilliCode = round(rem(Value, 1) * 1000);
    CodeTable = Get_Designation_Code_Table();

    if isempty(find([CodeTable.MilliCode] == MilliCode, 1))
        MilliCode = NaN;
        Valid = false;
        return
    end

    Valid = true;
end

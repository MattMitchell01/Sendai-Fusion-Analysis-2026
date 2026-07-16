function [CorrectedAnalysisData, Success] = Apply_Designation_Code(MilliCode, CorrectedAnalysisData, ...
    TraceNumberIndex, VideoTimeVector, Options, FigureHandles, ConfirmZeroClick, ClipWidth, ManualFocusSubtractIndices)
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
% Apply_Designation_Code  Looks up MilliCode in Get_Designation_Code_Table,
% applies the 'new'-vs-'repick' mistake guard (reading the CURRENT live
% Designation, so it reflects whatever a prior step of a redesignation
% chain just set), and dispatches to the click picker if one is needed.
%
% Shared by Correct_Designations.m (top-level PlotNumber.Code prompt) and
% the mid-picker redirect prompts inside Fix_Fusion_Wait_Time.m /
% Fix_Unbind_Wait_Time.m -- both apply a code through this single function
% so the guard/assignment/dispatch logic exists in exactly one place.
%
% ConfirmZeroClick (default false via nargin<7): when a MilliCode needs
% zero clicks (No Fusion / Other), the top-level round prompt applies it
% silently as it always has (ConfirmZeroClick=false). But when reached via
% a mid-picker redirect, the reviewer may be correcting the designation
% itself (not just a picked frame) -- so ConfirmZeroClick=true still shows
% an explicit "are you sure" confirmation, which itself can redirect again
% (recursion, same as the click-based pickers).
%
% Success=false (struct returned untouched) when MilliCode is unrecognized,
% or the mistake guard fires (e.g. re-flagging a trace as something it
% already is). Callers should treat Success=false as "nothing changed,
% report why and let the caller decide what to do next" -- never as an
% error to throw.

    if nargin < 7
        ConfirmZeroClick = false;
    end
    if nargin < 8
        ClipWidth = [];
    end
    if nargin < 9
        ManualFocusSubtractIndices = [];
    end

    Success = true;
    CodeTable = Get_Designation_Code_Table();
    MatchIdx = find([CodeTable.MilliCode] == MilliCode, 1);

    if isempty(MatchIdx)
        fprintf('Unrecognized designation code .%d.\n', MilliCode);
        Success = false;
        return
    end

    Entry = CodeTable(MatchIdx);
    PreviousDesignation = CorrectedAnalysisData(TraceNumberIndex).Designation;

    if strcmp(Entry.Mode,'repick')
        if ~strcmp(PreviousDesignation, Entry.Designation)
            fprintf('Trace is currently "%s", not "%s" -- repick code .%d only applies when the designation is already correct.\n', ...
                PreviousDesignation, Entry.Designation, Entry.DecimalCode);
            Success = false;
            return
        end
        CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Correct Designation, Incorrect Wait Time';
    else
        if strcmp(PreviousDesignation, Entry.Designation)
            fprintf('Trace is already "%s" -- nothing to change.\n', Entry.Designation);
            Success = false;
            return
        end
        CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Incorrect Designation-Changed';
    end

    CorrectedAnalysisData(TraceNumberIndex).Designation = Entry.Designation;
    CorrectedAnalysisData(TraceNumberIndex).FusionData.Designation = Entry.Designation;

    if Entry.NumClicks > 0 && strcmp(Options.FixWaitTime,"y")
        if strcmp(Entry.ClickTarget,'Fuse')
            CorrectedAnalysisData = Fix_Fusion_Wait_Time(CorrectedAnalysisData,VideoTimeVector,FigureHandles,TraceNumberIndex,Options,Entry.NumClicks,ClipWidth,ManualFocusSubtractIndices);
        elseif strcmp(Entry.ClickTarget,'Unbound')
            CorrectedAnalysisData = Fix_Unbind_Wait_Time(CorrectedAnalysisData,VideoTimeVector,FigureHandles,TraceNumberIndex,Options,ClipWidth,ManualFocusSubtractIndices);
        end
    elseif ConfirmZeroClick
        Confirmed = false;
        while ~Confirmed
            fprintf('Re-designated to "%s". We all good here? (Enter/y = yes, or type a new code to change designation)\n', Entry.Designation);
            s = strtrim(lower(input('> ', 's')));

            if isempty(s) || strcmp(s,'y')
                Confirmed = true;
            elseif strcmp(s,'n')
                fprintf('No pick to redo for this designation -- type a new code if you want to change it, or press Enter/y to confirm.\n');
            else
                [NewMilliCode, Valid] = Parse_Designation_Code(s);
                if ~Valid
                    fprintf('Unrecognized input -- enter blank/y to confirm, or a designation code (.0,.1,.11,.2,.22,.3,.33,.9).\n');
                else
                    [CorrectedAnalysisData, RecurseSuccess] = Apply_Designation_Code(NewMilliCode, ...
                        CorrectedAnalysisData, TraceNumberIndex, VideoTimeVector, Options, FigureHandles, true, ClipWidth, ManualFocusSubtractIndices);
                    if RecurseSuccess
                        Confirmed = true;   % nested call (or its own further chain) already confirmed
                    end
                    % RecurseSuccess==false: reason already printed, reprompt the same question
                end
            end
        end
    end
end

function [RerunThisRound, CorrectedAnalysisData, ErrorCounter] = Correct_Designations(IncorrectPlotIndices,...
    ~,CurrentTraceRange,CorrectedAnalysisData, ErrorCounter,Options,VideoTimeVector,FigureHandles,ClipWidth,ManualFocusSubtractIndices)
%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
% Correct_Designations  Apply reviewer-typed PlotNumber.Code corrections.
%
% Tier-agnostic -- Get_Designation_Code_Table.m has no notion of which
% tier a trace came from; it just maps a milli-code to a designation
% change, so Low/Medium/High review loops can all call this unchanged.
%
% Code = PlotNumber.DesignationCode. The fractional part is parsed as an
% integer 0-999 (round(rem(Code,1)*1000)) rather than compared with
% float-range checks like the old .09<code<.11 style, so e.g. .11 and
% .115 can never be confused. See "Reference Guide For Prompt Codes.rtf"
% and Get_Designation_Code_Table.m for the full table.
%
% The actual guard/assignment/click-dispatch logic lives in
% Apply_Designation_Code.m, shared with the mid-picker redirect prompts in
% Fix_Fusion_Wait_Time.m/Fix_Unbind_Wait_Time.m -- called here with
% ConfirmZeroClick=false so a 0-click code (.0/.9) still applies silently
% at this top-level round prompt, exactly as before. (2nd argument,
% formerly PreviousAnalysisData, is no longer read: Apply_Designation_Code
% checks the mistake guard against CorrectedAnalysisData's own live
% Designation instead, which also means two corrections for the same trace
% within one round's list -- e.g. a redesignation followed by a repick --
% now compose correctly instead of the repick spuriously failing against a
% stale round-start snapshot.)
%
% A bad entry within a multi-entry list (e.g. "5.100, 9.100, 12.777, 3.100")
% is resolved on the spot, in place, rather than aborting the whole batch:
% whichever half of that one entry was bad (plot number or designation
% code) gets a targeted re-prompt for just that half, and the loop then
% carries on to the remaining entries. This replaced an earlier version
% that RerunThisRound='y'+break on the first bad entry -- which silently
% dropped every entry after the bad one (never applied, not even reported)
% while forcing the reviewer to retype the WHOLE list including entries
% that had already been applied for real a moment earlier (and would then
% likely fail again on retry, since Apply_Designation_Code's mistake guard
% rejects re-flagging a trace as what it already is).

    RerunThisRound = 'n';
    NumberIncorrect = length(IncorrectPlotIndices);

    for j = 1:NumberIncorrect
        CurrentIndex = floor(IncorrectPlotIndices(j));
        MilliCode = round(rem(IncorrectPlotIndices(j), 1) * 1000);
        MaxRoundPlots = length(CurrentTraceRange);

        % Bad plot number for THIS entry -- re-ask for a full
        % PlotNumber.Code for just this one, without touching any other
        % entry in the list (already-applied ones stay applied; later
        % entries are handled in their own loop iteration below).
        while isnan(CurrentIndex) || CurrentIndex < 1 || CurrentIndex > MaxRoundPlots
            fprintf('Invalid plot number %g (entry %d of %d) -- valid plot numbers are 1 to %d.\n', ...
                IncorrectPlotIndices(j), j, NumberIncorrect, MaxRoundPlots);
            RawFix = strtrim(input('Re-enter as PlotNumber.Code for just this entry (blank to skip it): ','s'));
            if isempty(RawFix) || ~contains(RawFix,'.') || isnan(str2double(RawFix))
                if ~isempty(RawFix)
                    fprintf('Still not a valid PlotNumber.Code -- skipping this entry.\n');
                end
                CurrentIndex = NaN;
                break
            end
            FixVal = str2double(RawFix);
            CurrentIndex = floor(FixVal);
            MilliCode = round(rem(FixVal, 1) * 1000);
        end

        if isnan(CurrentIndex)
            continue   % skipped -- move on to the next entry, if any
        end

        TraceNumberIndex = CurrentTraceRange(CurrentIndex);

        if strcmp(Options.QuickModeNoCorrection,'y')
            CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Incorrect Designation-Not Changed';
            continue
        end

        % Bad/unrecognized designation code for THIS entry -- the plot
        % number half was already fine, so re-ask for just the decimal
        % code (Parse_Designation_Code's decimal-only form, same convention
        % as the mid-picker redirect prompt), not the whole PlotNumber.Code
        % again.
        while true
            [CorrectedAnalysisData, Success] = Apply_Designation_Code(MilliCode, CorrectedAnalysisData, ...
                TraceNumberIndex, VideoTimeVector, Options, FigureHandles, false, ClipWidth, ManualFocusSubtractIndices);
            if Success
                break
            end
            RawFix = strtrim(input(sprintf(['Enter a corrected code for plot %d -- .0,.1,.11,.2,.22,.3,.33,.9 ' ...
                '(blank to skip this trace for now): '], CurrentIndex),'s'));
            if isempty(RawFix)
                break   % leave this trace uncorrected, move on
            end
            [NewMilliCode, Valid] = Parse_Designation_Code(RawFix);
            if ~Valid
                fprintf('Unrecognized code -- skipping this entry.\n');
                break
            end
            MilliCode = NewMilliCode;
        end
    end

    ErrorCounter = ErrorCounter + NumberIncorrect;
end

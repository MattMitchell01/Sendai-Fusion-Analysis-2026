function [CorrectedAnalysisData] = Fix_Fusion_Wait_Time(CorrectedAnalysisData, VideoTimeVector, FigureHandles, TraceNumberIndex, Options, NumClicksNeeded, ClipWidth, ManualFocusSubtractIndices)
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
% Fix_Fusion_Wait_Time  Blown-up single-trace figure for click-picking fuse
% frame(s). Tier-agnostic -- works the same regardless of which tier
% TraceNumberIndex's trace came from. Named to match Fix_Unbind_Wait_Time.m,
% its Unbound-side counterpart.
%
% NumClicksNeeded: 1 for a 1-Fuse trace (single fusion event), 2 for a
% 2-Fuse trace (two sequential ginput clicks, one per fusion event,
% sorted ascending so FuseFrameNumbers(1) is always the earlier one).
%
% Fuse frame(s) are stored back into FusionData.FuseFrameNumbers in
% GLOBAL coordinates (Part B's convention -- see Analyze_Trace_Data.m:
% 220-221 in Part B) even though picking happens against the CLIPPED
% trace, per Clip_Trace_For_Review.m.
%
% Only ever overwrites FuseFrameNumbers/BindtoFusionTime -- never
% SeekFusionData (FuseFramesGlobal/FuseFramesClipped/
% BindtoFusionNumFrames/BindtoFusionTime), Part B's own frozen record of
% its original call (Analyze_Trace_Data.m), guaranteed present on every
% trace, same untouchable-by-Part-C invariant as the top-level
% AlgoDesignation field. FusionData itself has no BindtoFusionNumFrames/
% StandardBindFrameNum fields at all (removed from the schema -- see Part
% B's design notes on output fields) -- BindFrameNum below is only
% ever a local variable here, used for the StandardBindTime/
% BindtoFusionTime calculation, never written back to the struct.
%
% The post-click "We all good here?" prompt also accepts a new designation
% code (decimal-only, e.g. ".3" for Unbound -- see Get_Designation_Code_Table.m),
% redirecting via Apply_Designation_Code without the reviewer needing to
% back out to the main round prompt. FusionData is held in a LOCAL copy and
% only written back to CorrectedAnalysisData at the very end, guarded by
% ~Redirected -- so if the reviewer redirects away from this designation
% (the click(s) just made turn out to be for the wrong event entirely),
% that abandoned click never overwrites whatever the redirect correctly set.

    TraceStruct   = CorrectedAnalysisData(TraceNumberIndex);
    FusionData    = TraceStruct.FusionData;
    FrameNumFound = TraceStruct.FrameNumFound;
    % StandardBindFrameNum == FrameNumFound always (Part B convention --
    % see PART_B_HANDOFF_NOTES.md section 1a, UniversalData no longer
    % exists).
    BindFrameNum  = FrameNumFound;
    StandardBindTime = VideoTimeVector(BindFrameNum) / (1000*60);   % ms -> minutes

    Redirected = false;
    AskUserAgain = 'y';
    while strcmp(AskUserAgain,'y')

        set(0,'CurrentFigure',FigureHandles.FixWaitPlot)
        [~, ClippedFrameNums] = Plot_Trace_With_Focus_Markers(TraceStruct, ...
            strcat("Trace ID = ", num2str(TraceNumberIndex)), ClipWidth, Options, ManualFocusSubtractIndices);

        LineToPlot = ylim;

        if NumClicksNeeded == 1
            fprintf('Select a Frame (If designation is wrong you can enter a new code in a second).\n');
        else
            fprintf('Select %d Frames, one at a time (If designation is wrong you can enter a new code in a second).\n', NumClicksNeeded);
        end

        % First click = bright green, second click (2-Fuse only) = dark
        % green -- not magenta/red, since the trace itself is already
        % plotted red for a 2-Fuse trace and a matching marker would
        % compete with it. Matches Plot_Current_Trace.m's master-grid
        % marker colors for the same two fuse frames.
        ClickColors = {[0 1 0], [0 0.5 0]};

        PickedGlobalFrames = zeros(1,NumClicksNeeded);
        for c = 1:NumClicksNeeded
            [FrameNum, ~] = ginput(1);

            DistToFrames = (ClippedFrameNums - FrameNum).^2;
            IdxOfMinValue = find(DistToFrames == min(DistToFrames),1);
            ClippedPick = ClippedFrameNums(IdxOfMinValue);

            PickedGlobalFrames(c) = ClippedPick + FrameNumFound;   % clip(1) == global frame FrameNumFound+1

            plot([ClippedPick, ClippedPick], LineToPlot, '--', 'Color', ClickColors{c})
            drawnow
        end

        PickedGlobalFrames = sort(PickedGlobalFrames);
        FusionData.FuseFrameNumbers = PickedGlobalFrames;

        FusionData.BindtoFusionTime = VideoTimeVector(PickedGlobalFrames) / (1000*60) - StandardBindTime;
            % Vectorized over ALL picked frames, matching Part B's own
            % formula exactly (Analyze_Trace_Data.m, current
            % bindToFusionTime computation) -- a 2-Fuse trace gets one value
            % per fuse frame (bind-to-first-fuse AND bind-to-second-fuse),
            % not just the first.
            %
            % Indexes VideoTimeVector (raw video time in ms, indexed by true
            % GLOBAL frame number), NOT TraceStruct.TimeVector -- the
            % latter is a per-trace array that Part A slices to start at
            % StartAnalysisFrameNumber (index 1 == global frame
            % FrameNumFound+1) and pre-converts to minutes
            % (Create_Video_Matrix_Auto_Time_Vector.m), so indexing it with
            % a GLOBAL frame number like PickedGlobalFrames silently reads
            % the wrong frame's time (off by FrameNumFound) or throws an
            % out-of-bounds error for any event near the end of a long
            % trace -- exactly the bug PART_B_HANDOFF_NOTES.md section 4
            % describes Part B having just fixed on its own side by
            % switching to VideoTimeVector. This mirrors that fix exactly.

        title(strcat("Bind to Fuse wait time(s) = ", num2str(FusionData.BindtoFusionTime), " min"))

        ConfirmLoop = true;
        while ConfirmLoop
            fprintf('We all good here? (Enter/y = yes, n = redo click(s), or a new code: .0,.1,.11,.2,.22,.3,.33,.9)\n');
            s = strtrim(lower(input('> ', 's')));

            if strcmp(s,'n')
                disp('Lets try again then')
                ConfirmLoop = false;   % outer while redoes the click(s)
            elseif isempty(s) || strcmp(s,'y')
                AskUserAgain = 'n';
                ConfirmLoop = false;
                CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Reviewed, Fuse frame chosen by user';
            else
                [MilliCode, Valid] = Parse_Designation_Code(s);
                if ~Valid
                    fprintf('Unrecognized input -- enter blank/y to confirm, n to redo, or a designation code (.0,.1,.11,.2,.22,.3,.33,.9).\n');
                    % ConfirmLoop stays true -- reprompt this same question, no replot/reclick
                else
                    [CorrectedAnalysisData, Success] = Apply_Designation_Code(MilliCode, CorrectedAnalysisData, ...
                        TraceNumberIndex, VideoTimeVector, Options, FigureHandles, true, ClipWidth, ManualFocusSubtractIndices);
                    if Success
                        Redirected = true;
                        AskUserAgain = 'n';
                        ConfirmLoop = false;   % nested call already ran its own confirm loop -- don't touch ChangedByUser here
                    end
                    % Success==false: Apply_Designation_Code already printed why -- reprompt this same question
                end
            end
        end
    end

    if ~Redirected
        CorrectedAnalysisData(TraceNumberIndex).FusionData = FusionData;
    end

end

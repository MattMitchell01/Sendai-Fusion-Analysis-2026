function [CorrectedAnalysisData] = Fix_Unbind_Wait_Time(CorrectedAnalysisData, VideoTimeVector, FigureHandles, TraceNumberIndex, Options, ClipWidth, ManualFocusSubtractIndices, ClearedIndexRanges)
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
% Fix_Unbind_Wait_Time  Blown-up single-trace figure for click-picking the
% unbind frame on a trace designated Unbound. Tier-agnostic, mirrors
% Fix_Fusion_Wait_Time.m's interaction pattern (single click) but writes
% into the top-level UnboundData sub-struct (mirrors FusionData's shape,
% kept separate so Unbound-specific data never mixes into FusionData)
% instead of FusionData itself. Named to match Fix_Fusion_Wait_Time.m, its
% Fuse-side counterpart.
%
% Frame is stored in GLOBAL coordinates (Part B's convention -- see
% Analyze_Trace_Data.m in Part B) even though picking happens against the
% CLIPPED trace, per Clip_Trace_For_Review.m.
%
% BindtoUnboundTime is computed from VideoTimeVector (raw video time in
% ms, indexed by true GLOBAL frame number), converted to minutes here
% (/(1000*60)) -- NOT from TraceStruct.TimeVector, which Part A slices to
% start at StartAnalysisFrameNumber and pre-converts to minutes
% (Create_Video_Matrix_Auto_Time_Vector.m), so it can only safely be
% indexed by a CLIPPED-coordinate frame, never a global one. This matches
% Part B's own current BindtoFusionTime/BindtoUnboundTime formula exactly
% (Analyze_Trace_Data.m) -- see PART_B_HANDOFF_NOTES.md section 4.
%
% Only ever overwrites UnboundFrameNumbers/BindtoUnboundTime -- never
% SeekUnboundsData (UnboundFramesGlobal/UnboundFramesClipped/
% BindtoUnboundNumFrames/BindtoUnboundTime), Part B's own frozen record of
% its original call (Analyze_Trace_Data.m): its own Find_Unbound_Frame
% call already populates SeekUnboundsData with a real value whenever it
% detected an unbind event (empty only if it detected none), same
% untouchable-by-Part-C invariant as the top-level AlgoDesignation field.
% A reviewer who agrees with Part B's designation but thinks its picked
% frame is wrong should use the repick code (.33), not the fresh code
% (.3), for exactly this reason. UnboundData itself has no
% BindtoUnboundNumFrames/StandardBindFrameNum fields at all (removed from
% the schema -- see the Part B design notes on output fields) --
% BindFrameNum below is only ever a local variable here, used for the
% StandardBindTime/BindtoUnboundTime calculation, never written back to
% the struct.
%
% The post-click "We all good here?" prompt also accepts a new designation
% code (decimal-only, e.g. ".0" for No Fusion -- see Get_Designation_Code_Table.m),
% redirecting via Apply_Designation_Code without the reviewer needing to
% back out to the main round prompt. UnboundData is held in a LOCAL copy
% and only written back to CorrectedAnalysisData at the very end, guarded
% by ~Redirected -- so if the reviewer redirects away from this
% designation (the click just made turns out to be for the wrong event
% entirely), that abandoned click never overwrites whatever the redirect
% correctly set.

    TraceStruct   = CorrectedAnalysisData(TraceNumberIndex);
    UnboundData   = TraceStruct.UnboundData;
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
            strcat("Trace ID = ", num2str(TraceStruct.VirusIDNumber), " -- pick UNBIND frame"), ClipWidth, Options, ManualFocusSubtractIndices, ClearedIndexRanges);

        LineToPlot = ylim;

        fprintf('Select a Frame (If designation is wrong you can enter a new code in a second).\n');
        [FrameNum, ~] = ginput(1);

        DistToFrames = (ClippedFrameNums - FrameNum).^2;
        IdxOfMinValue = find(DistToFrames == min(DistToFrames),1);
        ClippedPick = ClippedFrameNums(IdxOfMinValue);

        PickedGlobalFrame = ClippedPick + FrameNumFound;   % clip(1) == global frame FrameNumFound+1

        plot([ClippedPick, ClippedPick], LineToPlot, 'g--')
        drawnow

        UnboundData.UnboundFrameNumbers    = PickedGlobalFrame;
        % Indexes VideoTimeVector (raw video time in ms, GLOBAL-frame-
        % indexed), NOT TraceStruct.TimeVector -- see the matching comment
        % in Fix_Fusion_Wait_Time.m for why indexing the latter with a
        % global frame number is wrong (silently reads the wrong frame's
        % time, or throws out-of-bounds near the end of a long trace).
        UnboundData.BindtoUnboundTime      = VideoTimeVector(PickedGlobalFrame) / (1000*60) - StandardBindTime;

        title(strcat("Bind to Unbind wait time = ", num2str(UnboundData.BindtoUnboundTime), " min"))

        ConfirmLoop = true;
        while ConfirmLoop
            fprintf('We all good here? (Enter/y = yes, n = redo click, or a new code: .0,.1,.11,.2,.22,.3,.33,.9)\n');
            s = strtrim(lower(input('> ', 's')));

            if strcmp(s,'n')
                disp('Lets try again then')
                ConfirmLoop = false;   % outer while redoes the click
            elseif isempty(s) || strcmp(s,'y')
                AskUserAgain = 'n';
                ConfirmLoop = false;
                CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Reviewed, Unbind frame chosen by user';
            else
                [MilliCode, Valid] = Parse_Designation_Code(s);
                if ~Valid
                    fprintf('Unrecognized input -- enter blank/y to confirm, n to redo, or a designation code (.0,.1,.11,.2,.22,.3,.33,.9).\n');
                    % ConfirmLoop stays true -- reprompt this same question, no replot/reclick
                else
                    [CorrectedAnalysisData, Success] = Apply_Designation_Code(MilliCode, CorrectedAnalysisData, ...
                        TraceNumberIndex, VideoTimeVector, Options, FigureHandles, true, ClipWidth, ManualFocusSubtractIndices, ClearedIndexRanges);
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
        CorrectedAnalysisData(TraceNumberIndex).UnboundData = UnboundData;
    end

end

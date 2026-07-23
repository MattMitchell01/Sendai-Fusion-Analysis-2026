function Plot_Focus_Dot_Markers(ClippedTrace, ClippedFocusFrames)
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
% Plot_Focus_Dot_Markers  Draws the red-before/black-on dot pair at each
% focus event onto whatever axes/figure is already current (caller has
% already done hold on), against the already-plotted ClippedTrace.
%
% Shared by Plot_Trace_With_Focus_Markers.m (the blown-up fuse/unbind
% frame picker) and Plot_Current_Trace.m (the master grid, only for
% High-tier subgroups whose DiagnosticOverlay requests it -- H2a/H2b,
% where seeing every focus event next to the fuse/no-fusion call is the
% whole point of the review) so this exact look only needs to be
% maintained in one place.
%
% No blue "after" dot -- removed per the user, who doesn't care about the
% frame after a focus event, only the frame before (red) and the focus
% frame itself (black).

    % Draw in two passes (red, then black) rather than one pass per focus
    % event -- guarantees black always ends up on top wherever two focus
    % events' dots land on the same clipped index (e.g. two focus frames one
    % frame apart). Priority: black > red. A single interleaved loop (draw
    % both dots for event 1, then both for event 2, ...) let a later event's
    % red dot overpaint an earlier event's black dot whenever focus events
    % were close together -- this reordering fixes that.
    for f = 1:numel(ClippedFocusFrames)
        FocusIdx = round(ClippedFocusFrames(f));
        if FocusIdx-1 >= 1
            plot(FocusIdx-1, ClippedTrace(FocusIdx-1), 'ro', 'MarkerFaceColor','r')
        end
    end
    for f = 1:numel(ClippedFocusFrames)
        FocusIdx = round(ClippedFocusFrames(f));
        if FocusIdx >= 1 && FocusIdx <= numel(ClippedTrace)
            plot(FocusIdx, ClippedTrace(FocusIdx), 'ko', 'MarkerFaceColor','k')
        end
    end
end

function [Stats_Of_Failures,ShouldISkipToNextSUV,Dock_FrameNums,Fuse_FrameNums,DtoF_Time] =...
    Refine_Docking_Estimate(CurrTrace,Dock_FrameNums,Fuse_FrameNums,DtoF_NumFrames,...
    Stats_Of_Failures,TimeInterval)

ShouldISkipToNextSUV = 'n';

%Select range to look at from raw trace (not running median)
    HalfRange = 30;
    OldDockFrameNum = Dock_FrameNums(1);
    OldHFuseFrameNum = Fuse_FrameNums(1);
    %Make sure that the selection range doesn't go outside the
    %data
        while (OldDockFrameNum - HalfRange) <= 0
            HalfRange = HalfRange - 1;
        end

if HalfRange < 5
    Reason_Failed = 'Fast fuse event at beginning of stream' ;
    Cross_Out_Plot(TraceWindow,Reason_Failed)
        Stats_Of_Failures.FastFuseBeg = Stats_Of_Failures.FastFuseBeg + 1;
    ShouldISkipToNextSUV = 'y';
end
%Collapse data into baseline and docked
    %Calc baseline
        Baseline = median(CurrTrace(OldDockFrameNum-HalfRange:OldDockFrameNum-1));
    %Calc int of docked SUV
            NumFramesToMed = min(HalfRange,DtoF_NumFrames(1));
        Docked_Int = median(CurrTrace(OldDockFrameNum:OldDockFrameNum+NumFramesToMed));

    %Collapsed docking test.  The first frame which is much closer to the
    %docked intensity than to the baseline is considered the
    %docking event.
        FramesToSelect_ForDocking = OldDockFrameNum-HalfRange:OldDockFrameNum+NumFramesToMed;
        TraceSelected_Dock = CurrTrace(FramesToSelect_ForDocking);
            DistToBaseline = (TraceSelected_Dock - Baseline).^2;
            DistToDocked = (TraceSelected_Dock - Docked_Int).^2;
        Collapsed_Trace_Dock = DistToDocked < 0.3*DistToBaseline;
        Collapsed_Trace_Dock_FrameNums = FramesToSelect_ForDocking(Collapsed_Trace_Dock);
        if isempty(Collapsed_Trace_Dock_FrameNums)
            NewDockFrameNum = OldDockFrameNum;
        else
            NewDockFrameNum = Collapsed_Trace_Dock_FrameNums(1);
        end
%Check to make sure that something screwy hasn't happened.
%If the new docking frame is off by more than 3 frames from
%the old one, assume that the old one is correct.
    if abs(NewDockFrameNum - OldDockFrameNum) > 3
        NewDockFrameNum = OldDockFrameNum;
    end

%Overwrite old values with the refined ones
    Dock_FrameNums(1) = NewDockFrameNum;
    DtoF_NumFrames = Fuse_FrameNums - Dock_FrameNums(1);
    DtoF_Time = TimeInterval*(DtoF_NumFrames); 
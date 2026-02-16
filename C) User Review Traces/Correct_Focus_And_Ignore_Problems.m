function [CurrTrace,CurrTimeVector,VectorOfFrameNumbers] = Correct_Focus_And_Ignore_Problems(CurrTrace,CurrTimeVector,UniversalData)
% WARNING: THE ALGORITHM BELOW EXACTLY MATCHES WHAT IS IN THE TRACE
% ANALYSIS PROGRAM. IF YOU CHANGE ANYTHING BELOW, MAKE SURE YOU CHANGE IT
% THERE AS WELL TO BE CONSISTENT.

% Focus problems are corrected by determining the difference in intensity 
% before and after the focus frame number (taking the median over some 
% number of frames before and after the focus frame number). This 
% difference is then subtracted from the intensity trace for all frame 
% numbers after the focus problem (essentially erasing the focus problem 
% from the trace which will be used to perform the analysis). The actual 
% frame number where the focus problem occurred is set to the median 
% intensity value before the focus event.

% NOTE: the numbers here are the shifted numbers (shifted to 
% account for time zero or other frames ignored at the beginning). They 
% will correspond to the frames as they appear in the traces within this 
% program, and not necessarily in the original data set in FIJI or whatever

        VectorOfFrameNumbers = 1:length(CurrTrace);

% Ignore frame numbers are simply deleted from the trace, together with their 
% corresponding values in the time vector
    if strcmp(UniversalData.IgnoreProblems,'y')

        for m=1:length(UniversalData.IgnoreFrameNumbers)
            CurrentIgnoreFrameNumber = UniversalData.IgnoreFrameNumbers(m);
            LogicalIndexToKeep = VectorOfFrameNumbers ~= CurrentIgnoreFrameNumber;

            VectorOfFrameNumbers = VectorOfFrameNumbers(LogicalIndexToKeep);
            CurrTrace = CurrTrace(LogicalIndexToKeep);
            CurrTimeVector = CurrTimeVector(LogicalIndexToKeep);
        end
    end

%Now we take care of the focus frame numbers (we do it after the ignore
%frames so those don't get included in how we correct for focus frames.

    if strcmp(UniversalData.FocusProblems,'y')
            for m=1:length(UniversalData.FocusFrameNumbers)
                currentfocusframenumber = UniversalData.FocusFrameNumbers(m);
                currfocusindex = find(VectorOfFrameNumbers == currentfocusframenumber);

                widthtoaverage = 3;
                offset = 1;
                numberpointseithersidetoreplace = 0;
                if currfocusindex+ widthtoaverage + offset < length(CurrTrace)
                    intafterfocus = median(CurrTrace(currfocusindex+offset:currfocusindex+offset+widthtoaverage));
                    intbeforefocus = median(CurrTrace(max(currfocusindex-offset-widthtoaverage,1):currfocusindex-offset));
                    diffromfocus = intafterfocus - intbeforefocus;

                    for b = currfocusindex-numberpointseithersidetoreplace:currfocusindex+numberpointseithersidetoreplace
                        CurrTrace(b) = intbeforefocus;
                    end
                    CurrTrace(b+1:end) = CurrTrace(b+1:end) - diffromfocus;
                end
            end
    end

end
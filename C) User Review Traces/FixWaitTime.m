function [CorrectedAnalysisData] = FixWaitTime(CorrectedAnalysisData,UniversalData,FigureHandles,TraceNumberIndex,Options)
            
         
            FusionData = CorrectedAnalysisData(TraceNumberIndex).FusionData;
            CurrTrace = CorrectedAnalysisData(TraceNumberIndex).Trace_BackSub;
            CurrTimeVector = CorrectedAnalysisData(TraceNumberIndex).TimeVector;
            BindFrameNum = UniversalData.StandardBindFrameNum;

            %Correct focus and ignore frames for the current trace
            [CurrTrace_Corrected,CurrTimeVector_Corrected,CurrFrameNums_Corrected] = Correct_Focus_And_Ignore_Problems(CurrTrace,CurrTimeVector,UniversalData);

            if strcmp(Options.UseRunMed,"y")
                [CurrTrace_Corrected] = Run_Med(CurrTrace_Corrected,Options);
            end

            %Set up while loop to determine new wait time. Ask user if it
            %looks good. If not, choose again.
            AskUserAgain = 'y';
            while AskUserAgain =='y'

                set(0,'CurrentFigure',FigureHandles.FixWaitPlot)
                cla
    
                plot(CurrTrace_Corrected,'b-')
                    % 240711 NOTE TO FUTURE BOB: we DON'T plot the corrected frame nums on x axis so
                    % that we will calculate the time to fusion correctly
                    % (since we go by the corrected time vector as done below). 
                    % This means that ignored frames will just be deleted
                    % and won't appear in the plotted trace. This is
                    % different than the flavi analysis, which is done
                    % from the difference in frame numbers.

                hold on
                title(strcat("Trace ID = ", num2str(TraceNumberIndex)))
    
                LineToPlot = ylim;
                XToPlot = [BindFrameNum, BindFrameNum];
                plot(XToPlot,LineToPlot,'k--')
    
                [FrameNum,Yval] = ginput(1);
    
                
                DistToFrames = (CurrFrameNums_Corrected - FrameNum).^2;

                IdxOfMinValue = find(DistToFrames == min(DistToFrames),1);
                
                Fuse1FrameNum = CurrFrameNums_Corrected(IdxOfMinValue);

     
                FusionData.FuseFrameNumbers = Fuse1FrameNum;
                FusionData.BindtoFusionNumFrames = Fuse1FrameNum-BindFrameNum;

                FuseFrameIdx = find(CurrFrameNums_Corrected == Fuse1FrameNum);
                FusionData.BindtoFusionTime = CurrTimeVector_Corrected(FuseFrameIdx)...
                    - UniversalData.StandardBindTime; 
                    %This is how it is calculated in the trace analysis program 
                    % Note: there was an update on July 16,2024, before that
                    % point it was calculated as
                    % TimeVector(FusionData.BindtoFusionNumFrames)-UniversalData.StandardBindTime;, which
                    % typically would have put us off by 1 frame (so ~30
                    % sec error most likely). This shouldn't really change
                    % much (small error relative to timescale of fusion)
                    % but is something to be aware of.
                  
                set(0,'CurrentFigure',FigureHandles.FixWaitPlot)
                hold on
                LineToPlot = ylim;
                XToPlot = [FusionData.FuseFrameNumbers(1), FusionData.FuseFrameNumbers(1)];
                plot(XToPlot,LineToPlot,'g--')
                drawnow
    
                title(strcat("Bind to Fuse wait time = ", num2str(FusionData.BindtoFusionTime), " sec"))


                % Ask user via console (y / n)
                fprintf('We all good here? (y/n)\n');
                s = strtrim(lower(input('> ', 's')));

                if strcmp(s,'n')
                    disp('Lets try again then')
                    AskUserAgain = 'y';
                else
                    % Default behavior: treat anything other than 'n' as yes
                    AskUserAgain = 'n';   

                    % Change this flag so that the wait time will be
                    % included in CDF program.
                    CorrectedAnalysisData(TraceNumberIndex).ChangedByUser = 'Reviewed, Fuse frame chosen by user';
                end
            end

            
            CorrectedAnalysisData(TraceNumberIndex).FusionData = FusionData;

            
end
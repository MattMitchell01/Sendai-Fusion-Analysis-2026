function Start_User_Review(varargin)

disp('====================================')
disp('Initiating User Review Trace Program.......')

close all

% Identify the path to the data you wish to analyze
    [DataFilename,DefaultPathname] = Load_Data(varargin);
        % Nested function

    DataFilePath = strcat(DefaultPathname,DataFilename);

% Setup options. Ask if they look ok.
    [Options] = Setup_Options_User_Review();

    Options
        AreOptionsOk = input('   Options look good? (y to proceed) :','s');
            
            if strcmp(AreOptionsOk,'y')
                disp('   Excellent. Lets proceed with the analysis.')
                disp('   ...')
            else
                disp('   Options no good? Then change em and run the program again!')
                disp('   Program terminated.')
                disp('====================================')
                return
            end
    
% Extract data that we will need
    InputData = open(DataFilePath);
    DataToSave = InputData.DataToSave;
    UniversalData = InputData.DataToSave.OtherDataToSave.UniversalData;
    PreviousAnalysisData = InputData.DataToSave.CombinedAnalyzedTraceData;
    NumTraces = length(PreviousAnalysisData);
    
% Define variables we will need
    CorrectedAnalysisData = PreviousAnalysisData;
    NumTracesToReview = NumTraces - Options.StartingTraceNumber - 1;
    NumReviewRounds = ceil(NumTracesToReview/Options.TotalNumPlots);
    TraceCounter = Options.StartingTraceNumber;
    ErrorCounter = 0;

    DataCounters.CurrentTraceNumber = Options.StartingTraceNumber;
    DataCounters.CurrentErrorRate = 0;

% Create master window with subplots
    [FigureHandles] = Create_Master_Window(Options);
    if strcmp(Options.FixWaitTime,"y")
        FigureHandles.FixWaitPlot = figure(2);
    end
    if strcmp(Options.GrabExampleTrace,'y')
        FigureHandles.ExTrace = figure(3);
    end

% Review plots round by round
    for b = 1:NumReviewRounds

        disp(strcat('Round-', num2str(b),'-of-', num2str(NumReviewRounds)))
        
        if b ~= NumReviewRounds
            CurrentTraceRange = TraceCounter:TraceCounter + Options.TotalNumPlots - 1;
            TraceCounter = max(CurrentTraceRange) +1;
        else
            CurrentTraceRange = TraceCounter:NumTraces;
            for d = length(CurrentTraceRange) + 1: Options.TotalNumPlots
                set(FigureHandles.MasterWindow,'CurrentAxes',FigureHandles.SubHandles(d));
                cla
            end
        end
        
        % Load and plot the current round
        PlotCounter = 1;
        for i = CurrentTraceRange
            
            CurrentTraceNumber = i;
            CurrentVirusData = PreviousAnalysisData(i);
            CurrTrace = CurrentVirusData.Trace_BackSub;
            CurrTimeVector = CurrentVirusData.TimeVector;
            
            % Correct focus and ignore problems
            [CurrTrace_Corrected,CurrTimeVector_Corrected,~] = Correct_Focus_And_Ignore_Problems(CurrTrace,CurrTimeVector,UniversalData);
            
            if strcmp(Options.UseRunMed,"y")
                [CurrTrace_Corrected] = Run_Med(CurrTrace_Corrected,Options);
            end
            
            
            [FigureHandles] = Plot_Current_Trace(FigureHandles,CurrentVirusData,UniversalData,...
                CurrTrace_Corrected,PlotCounter,CurrentTraceNumber,Options);
            
            

            if strcmp(Options.GrabExampleTrace,'y') && PlotCounter == 1

                %Plot in separate figure to grab example trace
                set(0,'CurrentFigure',FigureHandles.ExTrace)
                plot(CurrTimeVector_Corrected,CurrTrace_Corrected,'b-')
                xlabel("Time (min)")
            end

            CorrectedAnalysisData(i).ChangedByUser = 'Reviewed By User';
            PlotCounter = PlotCounter +1;

        end

        if ~strcmp(Options.GrabExampleTrace,'y')
            
            % Ask user if we need to change any of the designations on the current round of plots
            RerunThisRound = 'y';
            

            %==============================================================================
            % New Version of Prompting the User: Through the Console
            %==============================================================================
            % old way used Matlab's modal dialog box which would bug out
            
            while RerunThisRound == 'y'

                fprintf('%d/%d; Enter IncorrectPlotIndices\n', b, NumReviewRounds);
                fprintf('Format: 5.3, 9.3\n');
                fprintf('Press Enter = No Correction Needed | q = quit\n');
            
                s = strtrim(input('> ', 's'));
            
                % Enter → no correction
                if isempty(s)
                    RerunThisRound = 'n';
                    continue
                end
            
                % Quit
                if strcmpi(s, 'q')
                    disp('   You Chose To Quit')
                    disp('   Program terminated.')
                    disp('====================================')
                    return
                end
            
                % Replace commas with spaces so str2num works
                s = strrep(s, ',', ' ');
            
                % Convert to numeric vector
                IncorrectPlotIndices = str2num(s); 
            
                % Validate
                if isempty(IncorrectPlotIndices) || ~isvector(IncorrectPlotIndices)
                    disp('Invalid input. Use format: number.code, number.code   e.g. 5.3, 9.3')
                    RerunThisRound = 'y';
                    continue
                end 
            
                % Call Correct_Designations
                [RerunThisRound, CorrectedAnalysisData, ErrorCounter] = Correct_Designations( ...
                    IncorrectPlotIndices, PreviousAnalysisData, CurrentTraceRange, ...
                    CorrectedAnalysisData, ErrorCounter, Options, UniversalData, FigureHandles);
            
            end
            
            DataCounters.CurrentErrorCount = ErrorCounter;
            DataCounters.CurrentTraceNumber = CurrentTraceNumber;
            DataCounters.CurrentErrorRate = ErrorCounter/CurrentTraceNumber;
            
            DataCounters
            
            if  strcmp(Options.SaveAtEachStep,'y')
                DataToSave.DataCounters = DataCounters;
                Save_Data_At_Each_Step(DataFilename,DefaultPathname,DataToSave,CorrectedAnalysisData,Options)
            end
    
        else
            break %get out of for loop
        end
    end
      
    if ~strcmp(Options.GrabExampleTrace,'y')
        disp("   You have reviewed all the traces!")
        disp("   Give a fist pump. You're done!")
    
        if ~strcmp(Options.SaveAtEachStep,'y')
            DataToSave.DataCounters = DataCounters;
            Save_Data_At_Each_Step(DataFilename,DefaultPathname,DataToSave,CorrectedAnalysisData,Options)
        end
        
        
        disp('====================================')

    
    else
        disp('===Example Trace Has Been Generated. Youre done!====')
    end

end

function [DataFilenames,DefaultPathname] = Load_Data(varargin)

    disp('   Please select the file of analyzed traces.')
    disp('   It should be the output of the Trace Analysis or User Review Program (if re-reviewing).')
    disp('   It should be a .mat file in the save location you specified...')

    if length(varargin) == 1
        [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed',...
            char(varargin{1}),'Multiselect', 'on');
    else
        [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed', 'Multiselect', 'on');
    end

    disp("   Awesome - let's continue!")
end

function Save_Data_At_Each_Step(DataFilename,DefaultPathname,DataToSave,CorrectedAnalysisData,Options)

SaveDataFolder = DefaultPathname;

SaveDataFolder = strcat(SaveDataFolder,'/AnalysisReviewed/');
if exist(SaveDataFolder,'dir') == 0
    mkdir(SaveDataFolder);
end

DataFilenameWOExt = DataFilename(1:end-4);
    if ~isempty(CorrectedAnalysisData)
        DataToSave.ReviewOptions = Options;
        DataToSave.CombinedAnalyzedTraceData = CorrectedAnalysisData;

        CurrSaveFilePath = strcat(SaveDataFolder,DataFilenameWOExt,Options.Label,'.mat');
        save(CurrSaveFilePath,'DataToSave');

        disp('   ...')
        disp("   Output file saved to : " + CurrSaveFilePath)
    end
end
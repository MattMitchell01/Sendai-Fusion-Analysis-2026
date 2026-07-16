function Start_CDFAnalysis(varargin)

disp('====================================')
disp('Initiating Fit CDF Program.......')

    close all
    
   % Navigate to files
        [DataFilenames,DefaultPathname] = Select_Files(varargin);
        
        if iscell(DataFilenames)
            NumberDataFiles = length(DataFilenames);
        else 
            NumberDataFiles = 1;
        end    
    
    % Load options and verify
        [Options] = Setup_Options_CDFAnalysis();
    
        Options
            % AreOptionsOk = input('   Options look good? (y to proceed) :','s');
            % 
            %     if strcmp(AreOptionsOk,'y')
            %         disp('   Excellent. Lets proceed.')
            %         disp('   ...')
            %     else
            %         disp('   Options no good? Then change em and run the program again!')
            %         disp('   Program terminated.')
            %         disp('====================================')
            %         return
            %     end
    
    %Set up empty result structure and initialize figures

    AllResults = [];
    CDFReports = [];
    [FigureHandles] = Initialize_Figures(Options);

    % Extract data for each file
    for FileNumber = 1:NumberDataFiles
        
            if iscell(DataFilenames)
                CurrentFilename = DataFilenames{1,FileNumber};
            else
                CurrentFilename = DataFilenames;
            end 
            CurrDataFilePath = strcat(DefaultPathname,CurrentFilename);
            InputData = open(CurrDataFilePath);
            CurrCDFData = InputData.CDFData;
    
            % Set up data label for current file
                if strcmp(Options.DataLabel,'Filename')
                    IdxOfDot = find(CurrentFilename=='.');
                    FilenameWOExtension = CurrentFilename(1:IdxOfDot-1);
                    CurrCDFData.DataLabelForPlot = FilenameWOExtension;
        
                elseif strcmp(Options.DataLabel,'CDF Name')
                    CurrCDFData.DataLabelForPlot = CurrCDFData.Name;
                else
                    disp("!!!Error: Options.DataLabel is incorrect")
                    StopProgramNow
                end

            AllDataToAnalyze(FileNumber).CDFData = CurrCDFData;
    end

    % Set up results report (basic info, will add as we go along)
        [CDFReports] = CDF_Report_Setup(AllDataToAnalyze);
    
    % Plot all CDFs in CDF window (we always do this).
        [FigureHandles] = Plot_CDFs(AllDataToAnalyze,Options,FigureHandles);

    % Plot efficiencies (depending on options this gets done in different
    % ways (e.g. bootstrap error, etc))
        [FigureHandles] = Plot_Efficiencies(AllDataToAnalyze,Options,FigureHandles);
        
    % We skip over everything else if SimpleComparisonOnly = 'y'
    if strcmp(Options.SimpleComparisonOnly,'n') 
        % Run different fits and statistical analyses, depending on what
        % was selected in the Options.

        if strcmp(Options.CalcRandomParam,'y')
            [AllDataToAnalyze,CDFReports] = Randomness_Parameter_Analysis(AllDataToAnalyze,CDFReports, Options);
        end

        if strcmp(Options.RunFit_SingleFitMultipleCDFs,'y') || strcmp(Options.RunFit_MultipleFitsSingleCDF,'y')
            [AllDataToAnalyze,CDFReports] = Fit_CDFs_And_Plot_Fits(AllDataToAnalyze,CDFReports, FigureHandles, Options);
        end
    
    
 
    
        disp("   ...")
        disp("   Fitting and statistics reports complete!")

    elseif strcmp(Options.SimpleComparisonOnly,'y')
        disp("   Simple comparison analysis completed.")
    end

    % Display Results Report
        CDF_Report_Display(CDFReports);
    
    disp("   I think my work here is done.")
    disp('====================================')

end
    
function [DataFilenames,DefaultPathname] = Select_Files(varargin)
%First, we load the .mat data files.
    disp('   Please select the file(s) of the compiled CDF.')
    disp('   They should be the output of the D) Compile CDF program')
    disp('   They should be .mat files in the same folder...')

        if length(varargin) == 1
            [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed',...
                char(varargin{1}),'Multiselect', 'on');
        elseif length(varargin) == 2
            DefaultPathname = varargin{1,1}; DataFilenames = varargin{1,2};
        else
            [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed', 'Multiselect', 'on');
        end

    disp("   Awesome - let's continue!")
end
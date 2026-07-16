function Start_Compile_CDF(varargin)
%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%

disp('====================================')
disp('Initiating CDF Program.......')

    close all
    
    % Navigate to files
        [DataFilenames,DefaultPathname] = Select_Files(varargin);
            
        if iscell(DataFilenames)
            NumberDataFiles = length(DataFilenames);
        else 
            NumberDataFiles = 1;
        end    
    
    % Setup options
        [Options] = Setup_Options_Compile_CDF();
    
        Options
            AreOptionsOk = input('   Options look good? (y to proceed) :','s');

            if strcmp(AreOptionsOk,'y')
                disp('   Excellent. Lets proceed.')
                disp('   ...')
            else
                disp('   Options no good? Then change em and run the program again!')
                disp('   Program terminated.')
                disp('====================================')
                return
            end

    %Set up empty result structure and initialize figures
        [FigureHandles] = Initialize_Figures(Options);

    % Deal with combining files and correcting data if chosen
        if strcmp(Options.CombineMultipleFiles,'y') && NumberDataFiles > 1

            % Combine multiple files first
            [DataToSave, EffectiveDeadTimes] = Combine_Trace_Analysis_Data(DefaultPathname, DataFilenames, NumberDataFiles);
            DataFileName = strcat(Options.ComboFileName,'.mat');
            DataFilePath = strcat(DefaultPathname,DataFileName);

        elseif strcmp(Options.CombineMultipleFiles,'n') && NumberDataFiles == 1
            % If not combining data, just load the data set
                Options.ComboFileName = "Ignored";
                DataFileName = DataFilenames;
                DataFilePath = strcat(DefaultPathname,DataFileName);

                InputData = open(DataFilePath);
                DataToSave = InputData.DataToSave;

                if isfield(DataToSave,'OtherDataToSave') && isfield(DataToSave.OtherDataToSave,'EffectiveDeadTime')
                    EffectiveDeadTimes = DataToSave.OtherDataToSave.EffectiveDeadTime;
                else
                    EffectiveDeadTimes = NaN;
                end

            % If correcting data, do that, then continue on
            if Options.ExcludeOutlierPoints == 'y'
                [DataToSave] = Correct_Trace_Analysis_Data(DataToSave, Options);
            end

        else
            disp('   Error: You selected more than one file but didnt choose CombineMultipleFiles option')
            return
        end

    % Check each file's effective dead time against the user-set dead time
        Report_Effective_Dead_Time(EffectiveDeadTimes, Options.DeadTime);

    % Compile the CDF and plot
        [ErrorFlag,DataToSave,CDFData] = Compile_CDF_and_Plot(DataFileName,...
            DataFilePath,FigureHandles,Options,DataToSave);

    % Quit if error
        if ErrorFlag == true
            return
        end

    % Save data

        SaveFilename = CDFData.Name + ".mat";
            % Name is defined in Compile_CDF_and_Plot
        SaveFilePath = strcat(DefaultPathname,SaveFilename);
        %save(SaveFilePath,'DataToSave','CDFData','Options');
        save(SaveFilePath, 'DataToSave', 'CDFData', 'Options', '-v7.3'); %Using .MAT v7.3 for files >2GB


    % Finish
        disp("   ...")
        disp("   CDF compilation done.")
        % disp("   Output file name : " + Filename)
        disp("   Output file saved to : " + SaveFilePath)
        disp("   I think my work here is done.")
        disp('====================================')

end
    

function [DataFilenames,DefaultPathname] = Select_Files(varargin)
    disp('   Please select the file(s) of analyzed traces.')
    disp('   They should be the output of the B) Trace Analysis or C) User Review Program or (if you want to re-do a compilation) this program D).')
    disp('   They should be .mat files in the same folder...')

        if length(varargin) == 1
            [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed',...
                char(varargin{1}),'Multiselect', 'on');
        elseif length(varargin) == 2
            DefaultPathname = varargin{1,1}; DataFilenames = varargin{1,2};
        else
            [DataFilenames, DefaultPathname] = uigetfile('*.mat','Select .mat files to be analyzed', 'Multiselect', 'on');
        end    %[SortedBindtoFList,CurrentColor.DataPoints]=Generate_Test_Data();

    disp("   Awesome - let's continue!")
end

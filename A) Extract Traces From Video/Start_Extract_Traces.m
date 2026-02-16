function [] = Start_Extract_Traces(varargin)

% - - - - - - - - - - - - - - - - - - - - -

% Input:
% Run_Me_To_Start(), in this case the user navigates to the image video
%       stacks and also chooses the parent folder where the output analysis files will be saved
%   OR
% Run_Me_To_Start(DefaultPath), where DefaultPath is the directory to which 
%       the user will be automatically directed to find the image video
%       stacks. After choosing the video stacks, the user then chooses the 
%       parent folder where the output analysis files will be saved.
%   OR
% Run_Me_To_Start(DefaultPath,SavePath), where DefaultPath is as above, and 
%       SavePath is the parent folder where the output analysis files will be saved

% Output:
% A .mat file is created which saves all of the variables in the current 
% workspace. This file will be the input for the Lipid Mixing Trace Analysis program. Within 
% this file, the intensity traces for each viral particle which has been 
% found, together with additional data for each virus, will be in the 
% VirusDataToSave structure, as defined in Find_And_Analyze_Particles.m

% Note: This program has been designed to process many videos sequentially,
% but it has been tested with individual video streams, so keep that 
% in mind if you choose to process many videos at once.

% By Bob Rawle, Kasson Lab, University of Virginia, 2016
% Published online in conjunction with:
% Rawle et al., Disentangling Viral Membrane Fusion from Receptor Binding 
% Using Synthetic DNA-Lipid Conjugates, Biophysical Journal (2016) 
% http://dx.doi.org/10.1016/j.bpj.2016.05.048

% Updated by Bob Rawle, Williams College, 2023 and 2024

% Updated by Matthew D. Mitchell, Rawle Lab, Williams College (Jan 2026)
% - - - - - - - - - - - - - - - - - - - - -

disp('====================================')
disp('Initiating Extract Traces Program.......')

%Define which options will be used
[Options] = Setup_Options_Extract_Traces();


close all

%Load the image files, chosen by the user
[NumberOfFiles,FileOptions,NumberOfParameters] = Load_Video_Files(Options,varargin);

% Analyze each video stream one by one
for i = 1:NumberOfFiles
    for j = 1:NumberOfParameters
        CurrentOptions = FileOptions(i).Parameter(j).Options;
        CurrentFilename = CurrentOptions.VideoFilename;
        DefaultPathname = CurrentOptions.DefaultPathname;
        SaveParentFolder = CurrentOptions.SaveParentFolder;
        
        CurrStackFilePath = strcat(DefaultPathname,CurrentFilename);

        % Extract focus frame numbers, frame to find
        % the viruses, etc. from the corresponding text file if it is there
            
            if strcmp(CurrentOptions.ExtractAnalysisInputsFromConsole,'y')
                [CurrentOptions] = Console_Extract_Analysis_Inputs(Options);
                CurrentOptions
            
            elseif strcmp(CurrentOptions.ExtractInputsFromTextFile,'y')
                CurrentAnalysisTextFilePath = strcat(DefaultPathname,Options.AnalysisTextFilename);
                [CurrentOptions] = Extract_Analysis_Inputs(CurrentOptions,CurrentAnalysisTextFilePath);
                CurrentOptions  
                disp("   Confirmed: " + Options.AnalysisTextFilename + " was extracted and relevant options overwritten")
            else
                CurrentOptions  
                disp("   No analysis input text file used. Original options kept.")
            end
        % Print out options to commandline. Ask if they look ok.
            AreOptionsOk = input('   We good with the options above, boss? (y to proceed) :','s');
            
            if strcmp(AreOptionsOk,'y')
                disp('   Excellent. Lets proceed with the analysis.')
                disp('   ...')
            else
                disp('   Options no good? Then change em and run the program again!')
                disp('   Program terminated.')
                disp('====================================')
                return
            end

        % Now we call the function to find the virus particles and extract
        % their fluorescence intensity traces
        [Results,VirusDataToSave, OtherDataToSave,CurrentOptions] =...
            Find_And_Analyze_Particles(CurrStackFilePath,CurrentFilename, ...
                CurrentOptions);


        % Results are displayed in the command prompt window
        disp('   ...')
        disp('   All traces extracted successfully.')
        disp('   ...')

        disp( strcat('---Results: File_', num2str(i),'_of_',num2str(NumberOfFiles),'--- Param_',...
            num2str(j),'_of_',num2str(NumberOfParameters),'---'))
        Results



        % If selected, info is automatically grabbed from the data filenames and/or pathnames to make more 
        % informative save folder directory and output analysis filenames. The save 
        % folder is then created inside the parent directory.
        if strcmp(CurrentOptions.AutoCreateLabels,'y')
            [DataFileLabel,SaveDataPathname] = Create_Save_Folder_And_Grab_Data_Labels(DefaultPathname,...
                SaveParentFolder,CurrentOptions);

            % Add on extra label
            DataFileLabel = strcat(DataFileLabel,CurrentOptions.ExtraLabel);
        else
            % Otherwise, the label and save folder are defined as below.
            DataFileLabel = 'TestLabel';
            SaveDataPathname = SaveParentFolder;
            mkdir(SaveDataPathname);
        end

        % Analysis output file is saved to the save folder. All variables are saved.
    %     save(strcat(SaveDataPathname,DataFileLabel,'-Traces','.mat'));
        CurrSaveFilePath = strcat(SaveDataPathname,DataFileLabel,'','.mat');

        save(CurrSaveFilePath);
        disp('   ...')
        disp("   Output file saved to : " + CurrSaveFilePath)
        disp('------------------------------')


    end
end

disp("   Looks like we're all done here, boss!")
disp('   Thank you.  Come again.')
disp('====================================')
end
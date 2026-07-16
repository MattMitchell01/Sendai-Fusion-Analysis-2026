function [NumberOfFiles,FileOptions,NumberOfParameters] = Load_Video_Files(Options,varargin)
%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%

% Output:
% SaveFolderDir = directory where the analysis files will be saved
% StackFilenames = file names of the image video stacks (should be a stack of .tif)
% DefaultPathname = directory where the image video stacks are located

FileOptions = [];
InputPaths = varargin{1,1};
saveFile = fullfile(fileparts(mfilename('fullpath')), 'last_console_inputs.mat');

if ~strcmp(Options.CombineMultipleVideos,'y')

    % --- File selection (with memory) ---
    useSaved = false;
    if isfile(saveFile)
        try
            savedState = load(saveFile);
            pathOk = isfield(savedState, 'LastVideoDirectory') && ...
                     isfield(savedState, 'LastVideoFilename') && ...
                     ischar(savedState.LastVideoDirectory) && ...
                     ischar(savedState.LastVideoFilename) && ...
                     ~isempty(savedState.LastVideoDirectory) && ...
                     ~isempty(savedState.LastVideoFilename);
            if pathOk
                fprintf('\n');
                disp('  Last used metadata file:')
                fprintf('    %s%s\n\n', savedState.LastVideoDirectory, savedState.LastVideoFilename);
                resp = strtrim(input('  Use this file again? (y/n): ', 's'));
                if strcmpi(resp, 'y')
                    DefaultPathname = savedState.LastVideoDirectory;
                    StackFilenames  = savedState.LastVideoFilename;
                    useSaved = true;
                end
            end
        catch
        end
    end

    if ~useSaved
        disp('   All righty. Please select the metadata file for the data you wish to analyze.')
        disp('   It should be in the same folder as the data itself...')
        if ~isempty(InputPaths)
            startDir = InputPaths{1,1};
        else
            startDir = pwd;
        end
        [StackFilenames, DefaultPathname] = uigetfile('*.*','Select the metadata file', startDir, 'Multiselect', 'on');

        if isequal(StackFilenames, 0)
            error('No metadata file selected. Program terminated.')
        end

        % Save for next time
        if isfile(saveFile)
            savedState = load(saveFile);
        else
            savedState = struct();
        end
        savedState.LastVideoDirectory = DefaultPathname;
        savedState.LastVideoFilename  = StackFilenames;
        save(saveFile, '-struct', 'savedState');
    end

    % --- Save folder selection ---
    if length(InputPaths) >= 2
        SaveFolderDir = InputPaths{1,2};
    elseif length(InputPaths) == 1
        disp("   Great. Now select the location of the save folder...")
        SaveFolderDir = uigetdir(InputPaths{1},'Choose the directory where data folder will be saved');
    else
        disp("    Great. Now select the location of the save folder...")
        SaveFolderDir = uigetdir(DefaultPathname,'Choose the directory where data folder will be saved');
    end

    disp("   Awesome - let's continue!")

    %Determine the number of files selected by the user
    if iscell(StackFilenames)
        NumberOfFiles = length(StackFilenames);
    else
        NumberOfFiles = 1;
    end

else
    SaveFolderDir = InputPaths{1,2};
    for i= 1:Options.NumberOfVideosToCombine
        disp("   Select the metadata file for video #"+num2str(i))
        [StackFilenames{1,i}, DefaultPathnamesToCombine{1,i}] = uigetfile('*.*','Select the metadata file',...
            InputPaths{1,1},'Multiselect', 'on');
    end
    DefaultPathname = DefaultPathnamesToCombine{1,1};
    NumberOfFiles = 1;
    % We call this 1 file because that is what we will ultimately be combining it into
end



    for i= 1:NumberOfFiles
        CurrentOptions = Options;
        CurrentOptions.DefaultPathname = DefaultPathname;
        CurrentOptions.SaveParentFolder = SaveFolderDir;
        if ~strcmp(Options.CombineMultipleVideos,'y')
            if NumberOfFiles == 1
                CurrentOptions.VideoFilename = StackFilenames;
            else
                CurrentOptions.VideoFilename = StackFilenames{1,i};
            end
        else
            CurrentOptions.VideoFilename = StackFilenames{1,1};
            % We set the VideoFile name as referring to the first video.
            % That way, if we display the finding image in a later program
            % (e.g. trace analysis) it won't cause any problems

            CurrentOptions.VideoFilenamesToCombine = StackFilenames;
            CurrentOptions.DefaultPathnamesToCombine = DefaultPathnamesToCombine;
        end

        if strcmp(Options.ScanParameters,'y')
            [FileOptions,NumberOfParameters] = Setup_Parameter_Scan(FileOptions,i,CurrentOptions);
        else
            FileOptions(i).Parameter(1).Options = CurrentOptions;
            NumberOfParameters = 1;
        end
    end
end

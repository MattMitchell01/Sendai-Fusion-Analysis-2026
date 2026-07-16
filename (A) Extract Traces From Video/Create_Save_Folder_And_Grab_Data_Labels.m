function [DataFileLabel,SaveDataPathname] = Create_Save_Folder_And_Grab_Data_Labels(DefaultPathname,...
    SaveParentFolder,Options)

% Information is automatically grabbed from the filenames and/or pathnames to make more 
% informative save folder directory and output analysis filenames.

% WARNING: This script is highly specialized to the way that I format my folder 
% names (which contain much of the information about my experiment). So if you 
% want to use this script, you should modify the appropriate sections below 
% to match the way you format your data.

    Datalabel=  DefaultPathname;
            IndexofSlash = find(Datalabel=='/');
            DataLabelForSaveFolder = Datalabel(IndexofSlash(end- 3) : IndexofSlash(end- 2));
            InfoLabel = Datalabel(IndexofSlash(end- 3) : end);
            
            % Grab information which will be used to label the output analysis file.
            FileFolderInfo = Datalabel(IndexofSlash(end-1):end);
            IndextoStop = find(FileFolderInfo == '-');
            
            if numel(IndextoStop)==0
                DataFileLabel = 'NoLabel';
            else
                DataFileLabel = FileFolderInfo(2:IndextoStop(1)-1);
            end
            
            
            if strcmp(Options.GrabTotalIntensityOnly,'y')
                SaveFolderName = strcat(DataLabelForSaveFolder,'TotalIntensityOnly','/');
                DataFileLabel = '';
            else 
                SaveFolderName = strcat(DataLabelForSaveFolder,'Traces','/');
            end
    
% The folder where the output analysis files will be saved is created.
    SaveDataPathname = strcat(SaveParentFolder,SaveFolderName);
    mkdir(SaveDataPathname);

% Another folder, which copies the information in the parent folder which contains
% the data files, is also created in the same parent directory as the save folder.
% This is just for convenience in being able to access the information which is
% contained in the parent folder containing the data files. 
    InfoDataPathname = strcat(SaveParentFolder, InfoLabel);
    mkdir( InfoDataPathname);
end
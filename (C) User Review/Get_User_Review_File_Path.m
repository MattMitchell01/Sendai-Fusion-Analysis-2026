function UserReviewFilePath = Get_User_Review_File_Path(SelectedFilePath, Options)
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
% Get_User_Review_File_Path  Deterministically compute the one User Review
% File path for a given selected input file.
%
% The "User Review File" is the single .mat file per dataset that holds
% your corrected trace data plus CurrentQueuePosition (where you left
% off) -- written by Save_User_Review_File.m and loaded on resume.
%
% If SelectedFilePath already lives inside the review folder
% (Options.ReviewFolderName), it IS the User Review File -- return it
% as-is (resume case). Otherwise the path is computed fresh, one level
% below the selected file's parent folder.
%
% Pure function: never reads or writes anything, so it can be called
% safely on every run regardless of what gets (re)selected. This is what
% fixes the old nested-folder bug, where the save folder was recomputed
% from whatever file happened to be selected, so re-selecting an
% already-reviewed file stacked another folder layer each time.

    [ParentDir, BaseName, Ext] = fileparts(SelectedFilePath);
    [~, ParentFolderName]      = fileparts(ParentDir);

    if strcmp(ParentFolderName, Options.ReviewFolderName)
        UserReviewFilePath = SelectedFilePath;
    else
        UserReviewFilePath = fullfile(ParentDir, Options.ReviewFolderName, [BaseName Options.Label Ext]);
    end
end

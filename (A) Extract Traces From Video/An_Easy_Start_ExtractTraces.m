%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
% REMEMBER: Always check Setup_Options_Extract_Traces before you begin!!!

disp('====================================')
disp('   Hey there friend. Lets get started, shall we?')

    % Matt recommends pointing this at your master data folder - the one that
    % holds all of your individual assay data folders. That way, when the
    % metadata file picker opens, it opens close to what you actually want
    % instead of some random default folder, making it quicker to select.
    DataLocation = '/Users/littlem/Research/data';

    % If you'd rather manually choose the data location every time instead,
    % comment out the line above and uncomment the two lines below.
        % disp('   Choose the folder where your data lives...')
        % DataLocation = uigetdir();

    % I recommend pointing this at your master analysis folder - the parent
    % directory where all of your save folders live.
    ParentSaveFolderLocation = '/Users/littlem/Research/analysis';

    % If you'd rather manually choose the save folder location every time
    % instead, comment out the line above and uncomment the two lines below.
        % disp('   Choose the directory where the data folder will be saved...')
        % ParentSaveFolderLocation = uigetdir();

disp('   Awesome. Sending that to the start program now...')
% If you have already chosen the data and save locations you can simply copy and paste the
% line below for quicker workflow.
Start_Extract_Traces(DataLocation,ParentSaveFolderLocation);

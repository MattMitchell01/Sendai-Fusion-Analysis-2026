% REMEMBER: Always check Setup_Options_Extract_Traces before you begin!!!




    disp('====================================')
    disp('   Hey there friend. Lets get started, shall we?')
    disp('   First, choose the folder where you want to pick the data...')
        DataLocation = uigetdir();
        %DataLocation = '';

    
            
    disp('   Okay, now choose the directory where data folder will be saved...')
    disp('   This will probably be the parent save folder called: RBC analysis data')
        ParentSaveFolderLocation = uigetdir();
        %ParentSaveFolderLocation = '';

disp('   Awesome. Sending that to the start program now...')
% If you have already chosen the data and save locations you can simply copy and paste the
% line below for quicker workflow.
Start_Extract_Traces(DataLocation,ParentSaveFolderLocation);
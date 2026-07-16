
% REMEMBER: Always check Setup_Options_Fit before you begin!!!


disp('====================================')
disp('   Lets do this!')
disp('   First, choose the folder where you want to pick the data...')
    DataLocation_FitCDF = uigetdir();


disp('   Awesome. Sending that to the start program now...')
% If you have already chosen the data location you can simply copy and paste the
% line below for quicker workflow.
Start_CDFAnalysis(DataLocation_FitCDF);
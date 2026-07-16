%
% ------------------------------------------------------------------------
% Modified by Matthew D. Mitchell, Rawle Lab, Williams College, 2026
% (parameter / configuration updates to Bob Rawle's original script;
% original pipeline: Rawle et al., Biophysical Journal 2016,
% doi:10.1016/j.bpj.2016.05.048).
% ------------------------------------------------------------------------
%
% REMEMBER: Always check Setup_Options_Compile_CDF before you begin!!!


disp('====================================')
disp('   Lets do this!')
disp('   First, choose the folder where you want to pick the data...')
    %DataLocation_CDF = uigetdir();
    DataLocation_CDF = '/Users/littlem/Research/analysis';


disp('   Awesome. Sending that to the start program now...')
% If you have already chosen the data location you can simply copy and paste the
% line below for quicker workflow.
Start_Compile_CDF(DataLocation_CDF);

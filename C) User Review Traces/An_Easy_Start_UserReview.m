% REMEMBER: Always check Setup_Options_User_Review before you begin!!!
% You may also want to open: Reference Guide For Prompt Codes.txt


disp('====================================')
disp('   Howdy, partner. Lets get started, shall we?')
disp('   First, choose the folder where you want to pick the data.')
disp('   This should be the folder with the analyzed traces: output of Program B)...')
    DataLocation_UserReview = uigetdir();
    


disp('   Awesome. Sending that to the start program now...')
% If you have already chosen the data location you can simply copy and paste the
% line below for quicker workflow.
Start_User_Review(DataLocation_UserReview);
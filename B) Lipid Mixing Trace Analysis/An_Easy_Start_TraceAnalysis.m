% REMEMBER: Always check Setup_Options_Trace_Analysis.m before you begin!!!


disp('====================================')
disp('   Hey there friend. Lets get started, shall we?')
disp('   First, choose the folder where you want to pick the data...')
disp('   (This should be your extracted traces: the output of program A.)')
    DataLocation_TraceAnalysis = uigetdir();
    %DataLocation_TraceAnalysis = '';


disp('   Awesome. Sending that to the start program now...')
% If you have already chosen the data location you can simply copy and paste the
% line below for quicker workflow.
Start_Trace_Analysis(DataLocation_TraceAnalysis);
function Start_Trace_Analysis(DataLocation)
%
% ------------------------------------------------------------------------
% Substantially revised and expanded by Matthew D. Mitchell,
% Rawle Lab, Williams College, 2026.
% Original script by Bob Rawle (Kasson Lab, University of Virginia, 2016;
% Rawle et al., Biophysical Journal 2016, doi:10.1016/j.bpj.2016.05.048;
% updated by Prof. Bob Rawle, Williams College, 2024).
% ------------------------------------------------------------------------
%
% Start_Trace_Analysis  User-facing entry point for the new Part B pipeline.
%
% Input:
%   DataLocation - folder containing the Program A output .mat file.
%                  A file picker opens directed to this folder.
%
% Output:
%   A .mat file saved to <DataLocation>/Analysis/<filename>-AnalyzedTraces.mat
%   containing DataToSave.CombinedAnalyzedTraceData — compatible with Part C, D, E.

addpath(genpath(fileparts(mfilename('fullpath'))));

disp('====================================')
disp('   Please select the file of extracted traces from the raw video.')
disp('   It should be the output of the Extract Traces program.')
disp('   It should be a .mat file in the save location you specified...')

    [DataFilename, DataPathname] = uigetfile('*.mat', ...
        'Select extracted traces .mat file', DataLocation);

    if isequal(DataFilename, 0)
        disp('   No file selected. Program terminated.')
        disp('====================================')
        return
    end

    disp("   Awesome - let's continue!")

dataFilePath = fullfile(DataPathname, DataFilename);

% Analyze_Trace_Data derives its own save path from dataFilePath (one level
% up from the input file's folder, e.g. ExptFolder/Traces/file.mat ->
% ExptFolder/Analysis/label.mat) -- it takes no saveFilePath argument.
Analyze_Trace_Data(dataFilePath);

disp('   Thank you.  Come again.')
disp('====================================')

end

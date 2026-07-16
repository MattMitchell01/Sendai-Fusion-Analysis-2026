function [Options] = TextFile_Extract_Analysis_Inputs(Options,CurrentAnalysisTextFilePath)
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
% The script will automatically extract relevant numbers as specified below.
% For example, a typical text might be 'find 2 tzero 1 start 3'.
% Default sequence is 1, 2, 3: tzero=1 (time-zero frame), find=2 (finding
% image frame), start=3 (start analysis frame).
%
% Focus frame numbers are NOT read from analysis_inputs.txt. They are
% detected automatically from the OME-TIFF Z-position metadata in
% Create_Video_Matrix_Auto_Time_Vector.m and override anything the text
% file could have provided, so a 'foc' entry here would be ignored anyway.
% Do not add one.
%
% Do not add explanatory text into analysis_inputs.txt itself — the parser
% below does a plain substring search for each keyword, so stray occurrences
% of 'find', 'start', or 'tzero' inside comment text (e.g. "finding",
% "starting") can silently corrupt or break the extraction.

AnalysisInputsText = extractFileText(CurrentAnalysisTextFilePath);
AnalysisInputsText = char(AnalysisInputsText);

        Key = 'tzero';
        IndexOfKey = strfind(AnalysisInputsText, Key);
        if ~isempty(IndexOfKey)
           Options.TimeZeroFrameNumber = sscanf(AnalysisInputsText(IndexOfKey+length(Key)+1:end), '%i');
        end

        Key = 'start';
        IndexOfKey = strfind(AnalysisInputsText, Key);
        if ~isempty(IndexOfKey)
           Options.StartAnalysisFrameNumber = sscanf(AnalysisInputsText(IndexOfKey+length(Key)+1:end), '%i');
        end
        
        Key = 'find';
        IndexOfKey = strfind(AnalysisInputsText, Key);
        if ~isempty(IndexOfKey)
           Options.FrameNumToFindParticles = sscanf(AnalysisInputsText(IndexOfKey+length(Key)+1:end), '%i');
        end

        if isfield(Options,'ManuallyCorrectFind')
            if strcmp(Options.ManuallyCorrectFind,'y')
                Options.FrameNumToFindParticles  = Options.CorrectFindNumber;
            end
        end
        
        % If there are multiple videos to combine, then extract their
        % per-video frame counts as well.
        if strcmp(Options.CombineMultipleVideos,'y')
            
            Key = 'NumFrames';
            IndexOfKey = strfind(AnalysisInputsText, Key);
            if ~isempty(IndexOfKey)
               Options.NumFrames = sscanf(AnalysisInputsText(IndexOfKey+length(Key)+1:end), '%i');
            end
            
            for i= 2:Options.NumberOfVideosToCombine
                CurrentAnalysisTextFilePath = strcat(Options.DefaultPathnamesToCombine{1,i},Options.AnalysisTextFilename);
                AnalysisInputsText = extractFileText(CurrentAnalysisTextFilePath);
                AnalysisInputsText = char(AnalysisInputsText);
                
                Key = 'NumFrames';
                IndexOfKey = strfind(AnalysisInputsText, Key);
                if ~isempty(IndexOfKey)
                   NewNumFrames = sscanf(AnalysisInputsText(IndexOfKey+length(Key)+1:end), '%i');
                   Options.NumFrames = [Options.NumFrames, NewNumFrames];
                end
                
%                 WOULD HAVE TO FIX THIS BEFORE USING:
%                 Key = 'TimeDelay';
%                 IndexOfKey = strfind(AnalysisInputsText, Key);
%                 if ~isempty(IndexOfKey)
%                    NewTimeDelay = sscanf(AnalysisInputsText(IndexOfKey+length(Key)+1:end), '%i');
%                    Options.TimeDelayBetweenVideos = [Options.TimeDelayBetweenVideos NewTimeDelay];
%                 end
                
            end
        end

    % ---- Show the extracted analysis inputs so the user can eyeball them ----
    % Mirrors the clear per-field display the console entry point
    % (Console_Extract_Analysis_Inputs.m) uses, so the text-file path gives the
    % same at-a-glance check. The caller (Start_Extract_Traces.m) still asks for
    % final confirmation before proceeding. num2str([]) prints blank, which
    % flags any input whose keyword was missing from analysis_inputs.txt.
    fprintf('\n');
    disp('   Analysis inputs extracted from analysis_inputs.txt:')
    fprintf('      1. Time-zero frame     (tzero) :  %s\n', num2str(Options.TimeZeroFrameNumber));
    fprintf('      2. Finding image frame  (find)  :  %s\n', num2str(Options.FrameNumToFindParticles));
    fprintf('      3. Start analysis frame (start) :  %s\n', num2str(Options.StartAnalysisFrameNumber));
    fprintf('\n');

end

function [Options] = Setup_Options_Compile_CDF()
%
% ------------------------------------------------------------------------
% Modified by Matthew D. Mitchell, Rawle Lab, Williams College, 2026
% (parameter / configuration updates to Bob Rawle's original script;
% original pipeline: Rawle et al., Biophysical Journal 2016,
% doi:10.1016/j.bpj.2016.05.048).
% ------------------------------------------------------------------------
%

% -----General Options-----    
    
    Options.CombineMultipleFiles = 'y'; % 'y' if you want to combine multiple
            % analyzed files into a single CDF (same experimental
            % condition). 
            % Files to be combined need to be in the same folder so you can select
            % them together.
        Options.ComboFileName = 'ComboTest'; 
            % if 'y' above, this is the name of the output combo file. 
            % Otherwise output file will be extracted from input filename
            % "-CDF" and other labels will be attached as well depending on
            % options chosen below.

    % NOTES:
    % 1) If your CDF has 0 or 1 fusion events, then you will likely get an
    % error as well - the program wasn't designed for those situtations.

%------Efficiency and Unbound Calc Options-------
    % Determine what is included in the efficiency calculation. Basically, 
    % when we divide by the total number of viruses analyzed, do we exclude
    % any or not (such as excluding all of the viruses that unbind during the analysis window)?
        Options.ToIncludeInEfficiencyCalculation = 'Fuse1 and No Fusion';
        % Choices are:
            % 'Fuse1 and No Fusion'
            % 'Fuse1, Fuse2, and No Fusion'
            % 'Fuse1, Fuse2, Unbound, and No Fusion'
            % 'All'

    % Choose y if you want to calc the fraction unbound
    % This will show NumberUnbound./(NumberFuse1 + NumberNoFuse +
    % NumberUnbound); See Extract_Data.m if you want to change calculation.
        Options.ShowFractionUnbound = 'y';

%--------Data editing/trimming options: Don't change these options unless you ask Bob first ---------
    % NOTES: 
    % 1) ALL DATA EDITING AND TRIMMING NEEDS TO BE DONE AT THIS STAGE.
    % YOU CAN'T DO IT IN PROGRAM E).
    % 2) IF YOU TRIM/EDIT, YOU SHOULD NOTE IT IN YOUR LAB NOTEBOOK!

    % Exclude outlier points in the middle:
        % All fusion events in between the exclusion range will
        % not be included in the analysis.
        % NOTE: can't do this when combining multiple data sets. If you
        % want to correct and combine, you should first correct the
        % erroneous data set and then combine the corrected set with the
        % other data sets
        Options.ExcludeOutlierPoints = 'n';
            Options.LowEndOfExclusion = 110; % in whatever units your waiting time is in
            Options.HiEndOfExclusion = 116; % in whatever units your waiting time is in
            if Options.ExcludeOutlierPoints == 'y'
                Options.CombineMultipleFiles = 'n';
            end

    % Trim ends of data:
        % All fusion events with waiting times smaller or higher than the time cut off will
        % not be included in the analysis.
        
            Options.TimeCutoffLow = 8; % in whatever units your waiting time is in
            Options.TimeCutoffHigh = NaN;
                % NaN means no cut off on the high end

        % Define dead time (in whatever units your waiting time is in)
            Options.DeadTime = Options.TimeCutoffLow; % SeV people do this
            % Options.DeadTime = 0; % DENV people do this
                % This is the lag time at the beginning of the experiment before we
                % can measure fusion. It will be used to adjust many aspects of the
                % analysis below, so chat with Bob before changing.
                % **Right now only relevant for SeV fusion. DENV people should set it to 0.**

%-------------Other options: Don't change these options unless you ask Bob first -----------------

    Options.FusionTrigger = 'Binding';
        % 'Binding' or 'pH'
        % The data is extracted differently depending on the trigger

end

function [Options] = Setup_Options_CDFAnalysis()
%
% ------------------------------------------------------------------------
% Modified by Matthew D. Mitchell, Rawle Lab, Williams College, 2026
% (parameter / configuration updates to Bob Rawle's original script;
% original pipeline: Rawle et al., Biophysical Journal 2016,
% doi:10.1016/j.bpj.2016.05.048).
% ------------------------------------------------------------------------
%

% -----General Options-----
       
    % Decide on data label for plots
        Options.DataLabel = 'Filename';
            % 'Filename' = use the filename (w/o extension) to label data
            % 'CDF Name' = use the name given to the CDF when it was
            % compiled (this would also usually be the filename unless you
            % changed it)
   
    % Simple Comparison/Plotting Only
        Options.SimpleComparisonOnly = 'n';
            % This won't do any fitting or stats analysis. Just a quick
            % plotting of CDFs and Efficiencies. All Optional Analysis
            % Modules will be ignored.

    % Show comparison CDFs as normalized
        Options.NormalizeCDFsInComparisonPlot = 'y';
            % This won't apply to the fitting (for now)

% -----Optional Analysis Modules-----
% HELPFUL NOTE: Make sure you set Options.SimpleComparisonOnly = 'n' if you want to
% run any of these.

    Options.AccountForDeadTime = 'y';
        % Account for dead time in rand param calc and fitting modules
        % Dead time was defined in options for program D) Compile CDF, so refer
        % there if you want to change the number used for the dead time.
        % **Right now only relevant for SeV fusion. DENV people should set it to 'n'.**

    % ----Randomness Param and Nmin Calcs Module----
    % Will also calculate bootstrap confidence intervals for each statistic.
        Options.CalcRandomParam = 'y'; %'y' or 'n'
        Options.ConfInt_RandParam = 95;
        Options.NumberBootstraps_RandParam = 1000;

    % ------------Fit CDF Module----------------
        Options.RunFit_SingleFitMultipleCDFs = 'y'; %'y' or 'n'
            % If choose this option, only the first fit below will be run
            % on all data sets
        Options.RunFit_MultipleFitsSingleCDF = 'n'; %'y' or 'n'
            % If choose this option, all fits below will be run
            % on the first data set only
        Options.FitToPerform(1).FitType = '1 Exp';
        %Options.FitToPerform(1).FitType = 'Gamma';
        %Options.FitToPerform(2).FitType = 'Gamma';
        % Current FitType Options
            % 'Gamma'
            % '1 Exp'
    

% ------LEGACY OPTIONS BELOW: WORK IN PROGRESS----------------------------------------------
% 
% % All fusion events with waiting times smaller or higher than the time cut off will
% % not be included in the analysis
%     Options.TimeCutoffLow = 0; % in whatever units your waiting time is in
%     Options.TimeCutoffHigh = 90;
%         % NaN means no cut off on the high end
%         % Note: If TimeCutoffHigh = a number, the efficiency will only
%         % include fusion events with a wait time (would only impact old
%         % data where we sometimes had user corrected fusion events with
%         % no wait time)
% 
%  % NOTES:
%     % 1) If you only run one CDF, you may get an error when trying to plot the
%     % bar graph of the efficiencies with the error. This is normal. You can
%     % get around it by running multiple CDFs.
%     % 2) If your CDF has 0 or 1 fusion events, then you will likely get an
%     % error as well - the program wasn't designed for those situtations.
% 
% 

% 
% % ------Statistics Options--------
% 
%     % Choose y if you want to run a two-sided KS test on each CDF, comparing it
%     % to the next CDF in the for loop
%         Options.RunKSTest = 'n';
% 
%     % Choose y if you want to visualize and run the confidence intervals from the
%     % bootstrap calculation for the CDF
%         Options.RunBootstrap = 'y';
%         Options.ConfidenceInterval = 95;
%         Options.RunBootstrapMedian = 'n';
%     %     Options.ShowBootstrap = 'y';
%         Options.NumberBootstraps = 10000;
% 
    % Choose y if you want to calculate bootstrap errors for the fusion
    % efficiencies. Otherwise, errors will be set to zero.
        % Options.BootstrapEfficiencies = 'y';
        % Options.EfficiencyConfidenceInterval = 95;
        % Options.EfficiencyNumberBootstraps = 10000;
% 
% 
% %----- Options You are Unlikely to Change------
%     % Correct time extraction error (should usually keep as 'n')
%     Options.CorrectTimeVector = 'n';
% 
%     % This will show a plot of the median intensity of the viruses at the
%     % beginning of the trace (� some percentile)
%         Options.ShowBeginningIntensity = 'n';
% 
%     % Add data from nontraditional source (like manually entered data from
%     % another paper), and define how many extra data sets will be included
%         Options.AddExtraData = 'n';
%         Options.NumberExtra = 2;
% 
%     % Choose y if you want to scale the CDF according to the intensity jump
%     % upon fusion (this will only work if the intensity jump was calculated,
%     % and really only makes sense for fusion to tethered vesicles)
%         Options.ShowIntensityComparison= 'n';
% 


end

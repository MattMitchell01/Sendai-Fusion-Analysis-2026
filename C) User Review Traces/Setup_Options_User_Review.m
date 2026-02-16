function [Options] = Setup_Options_User_Review()
    
    % ------Easy reference guide for prompt codes-----
        % To enter at prompt: PlotNumber.DesignationCode
        % DesignationCode as follows:
        % .0 = No Fusion
        % .1 = 1 Fuse
        % .12 = 1 Fuse designation is already correct, but wait time is wrong
        % .2 = Abnormal (e.g. slow) fusion
        % .3 = Unbound event (see sharp decrease)
        % .9 = Hard To Classify, Ignore This One

    % ------Options to double check-------
        Options.Label = '-Revd';
        Options.ApplyFilter = 'n';
            % 'y' or 'n'. 
            % WARNING: If 'y', look in Plot_Current_Trace to make sure that
            % you are applying the filter you want. It will make some
            % traces be plotted in yellow according to your filter
            % condition.

    % -------Options you are less likely to change regularly-------
        Options.StartingTraceNumber = input("   Please enter starting trace number (1 is typical): ");

        Options.NumPlotsX = 7;
        Options.NumPlotsY = 4;
        Options.TotalNumPlots = Options.NumPlotsX*Options.NumPlotsY;
        Options.SaveAtEachStep = 'y';
        Options.QuickModeNoCorrection = 'n';
    
        Options.FixWaitTime = 'y'; %allow user to manually fix wait time when they choose 1 fuse.

        Options.UseRunMed = 'y'; % show running median instead of raw trace
        Options.RunMedHalfLength = 1; % num of data points on either side to include in running median

        Options.ShowBindFrame = 'n'; %Will draw a black dashed line
        Options.ExpandXAxis = 'y'; %will expand the x-axis so the beginning/end of the trace aren't obscured by y-axis
        
        Options.GrabExampleTrace = 'n'; %Use to grab example trace. 
            % The rest of user review will be ignored. It will grab the trace of the starting trace number.


        Options.AddPresetOptions = 'n';
        if strcmp(Options.AddPresetOptions, 'y')
            PresetOptionsDir = '/Users/bobrawle/Matlab/Virus Fusion Scripts/Preset Options/User Review Traces';
            [PresetOptionsFile, PresetOptionsDir] = uigetfile('*.m','Select pre-set options .m file',...
                char(PresetOptionsDir),'Multiselect', 'off');
            run(strcat(PresetOptionsDir,PresetOptionsFile));
        end
end
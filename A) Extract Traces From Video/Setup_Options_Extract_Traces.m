function [Options] = Setup_Options_Extract_Traces()

% ==============Options you'll probably change regularly==================
          
    Options.ScanParameters = 'n';
    Options.ScanRange = [100, 300];
        % y or n. Choose y if you are going to be extracting traces from 
        % the same file with different parameters. 
        % If so, you should make sure that everything is set up correctly 
        % in Setup_Parameter_Scan.m
        % Usually you will use this to determine the best threshold. The
        % program will show you the output of the threshold range you chose
        % and then terminate.
        
    Options.Threshold = 300; %Rem: Should not change much but dye [ ] in virus can be very different even for == samples
        % This is the number of counts above background which will be used
        % to detect virus particles. You will need to optimize this number 
        % for each set of imaging conditions and/or each data set. An optimal 
        % threshold value has been reached when you are able to see that the
        % program is accurately finding the particles in each image. To avoid
        % bias introduced by the optimization process, you should make sure 
        % that changing the optimal threshold value somewhat doesn't 
        % significantly affect your data. Assuming similar particle 
        % densities/intensities between data sets, and that the same imaging 
        % conditions are used, you shouldn't need to change the threshold 
        % value much if at all between data sets.

    Options.ExtraLabel =  strcat('_ExtractedTracestest, th = ', num2str(Options.Threshold));
        % This extra label will be added on to the end of the resulting filename 
        % for the output .mat file. Leave blank if you don't want anything added on.

    % ---------Image Visualization Options---------
    
    Options.MinImageShow = 400; % Normal = 350
    Options.MaxImageShow = 800; % Normal = 600/1600
        % These determine the minimum and maximum intensity counts that will 
        % be used to set the contrast for the grayscale images that are displayed.
        % The minimum value will be displayed as black and the maximum value 
        % will be displayed as white. 

%================= Options You Might Change While Optimizing For a Particular Data Type (after consulting Bob) ================
        
        Options.DeterminePixelOffset = 'y';
        % y or n. Do we need to account for XY Sample-Wide Drift

            Options.FFTOffset = 'y';
            % y or n. Uses Fast-Fourier-Transform to determine pixel offset over time.
            
                Options.FFTRegion = [0.30,0.70];    % Normal = [0.30,0.70] (Middle 40% of the image)


            Options.SingleParticleTrackingOffset = 'n';
            % y or n. Choose y if you want to determine pixel offset over time using single-particle tracking.
            % this will be done by selecting a stationary particle you wish to track and
            % then using that to determine the offset.  
            
                Options.SearchRadius = 20;
                % Sets the maximum box size in which the single particle will be tracked
      

        Options.TopHatBackgroundSubtraction = 'n';
            Options.DiskRadius = 10;
        % Top Hat background subtraction algorithm to remove the background
        % Disk Radius = radius of the disk in pixels that will be used to
        % determine 
    
    % ---------Inputs specific to each data set--------   
        Options.ExtractTimesFromMetaData = 'y';
            % y or n. Choose y if you want to automatically extract the time values from the metadata. This requires that
            % your data be in the OME.TIF format, and that you have the metadata.txt file. 
            % Otherwise, you will be manually defining the time vector below, and 
            % it will be assumed that your data is in the TIF format.
        Options.DetermineBindingTimeFromTimestamp = 'y';
            % y or n. Choose y if the timestamp from the first image is when the
            % viruses were added. Then the program will pull that time
            % automatically. It excludes the first image(s) from the trace analysis, 
            % since it was only used as a time zero marker. Note that y only makes 
            % sense if you are extracting times from the metadata.    

        Options.ExtractAnalysisInputsFromConsole = 'y';
            % y or n. Choose y if you want to enter the tzero, start, find,
            % ignore frames in via the console.
            % Otherwise you will have to use the .txt file

        Options.ExtractInputsFromTextFile = 'n';
            % 'y' OR 'n'
            % Choose 'y' to automatically extract the following inputs from 
            % an associated text file: TimeZeroFrameNumber, StartAnalysisFrameNumber, FrameNumToFindParticles,
            % FocusFrameNumbers, etc. In that case, the numbers put in for those 
            % values above will be ignored. Also, you should examine the notes 
            % in Extract_Analysis_Inputs.m for the proper formatting of the text file.
            
            Options.AnalysisTextFilename = "analysis_inputs.txt";
                % This is the file name from which the analysis inputs will automatically
                % be extracted. This file should be located in the same location as the 
                % data you are analyzing
        
        % -----------Manually defined inputs -----------
        % WARNING: These are only relevant if you don't extract inputs from
        % external .txt file, which is ununusual. So usually ignore these.
    
            Options.BindingTime = -3.5; % in min, doesn't have to be integer
                % This is the time of virus binding if that is known a priori.
                % Often this will happen if the viruses bound during a short time
                % window, were rinsed, and then imaging began. In that case, the
                % time of binding is set at the time when the viruses were added.
                % NOTE: This number will be ignored if you choose to extract inputs
                % from filename, or if the program determines this time
                % automatically (if Options.DetermineBindingTimeFromTimestamp = 'y').
            Options.FrameNumToFindParticles = 1;
                % This is the frame number that will be used as the finding image, 
                % to find viral particles. 
                % NOTE: This number will be ignored if you choose to extract inputs from associated text file.
                % NOTE: If you are averaging find frames, make sure to give yourself an offset
            Options.FindFramesToAverage = 0;
                % Number of frames before/after to average for the finding image
                % 1 means that you will average 3 frames, 2 means 5 frames, and so on
            Options.FocusFrameNumbers = 'NaN';
                % These are frame numbers where a focus event occurred and needs to be 
                % corrected in the intensity traces for each viral particle. This 
                % correction will be performed in the Trace Analysis program. 
                % 'NaN' for no events. [FrameNumber1,FrameNumber2, etc.] for multiple 
                % frames. NaN if not
                % NOTE: This number will be ignored if you choose to extract inputs from associated text file.
            Options.IgnoreFrameNumbers = 'NaN';
                % These are frame numbers that will be ignored in the trace analysis, 
                % but will not be corrected before and after as they are for focus events
                % 'NaN' for no events. [FrameNumber1,FrameNumber2, etc.] for multiple 
                % frames. 
                % NOTE: This number will be ignored if you choose to extract inputs from associated text file.
        
        %------Inputs that would be unusual to change for Sendai analysis------
            Options.TimeVector = [];
                    t = [0:5:95];
                    to_add = [0, 0.05, 0.1, 0.15, 0.2];
                for i = 1:length(t)
                    Options.TimeVector = [Options.TimeVector, t(i)+to_add];
                end
                % Input time vector used for all data. This will be
                % ignored if times extracted from metadata.
    
            if strcmp(Options.DetermineBindingTimeFromTimestamp, 'y')
                Options.TimeZeroFrameNumber = 1;
                % This is the frame number that will be used to determine time zero. Default is 1
                % NOTE: This number will be ignored if you choose to extract
                % inputs from associated text file.
                Options.StartAnalysisFrameNumber = Options.TimeZeroFrameNumber + 1;
                % This is the frame number where the viruses will begin to be analyzed. 
                % Anything before this will be ignored (although it may be used to 
                % determine time zero depending on what you have chosen above).
                % Default is TimeZeroFrameNumber + 1.
                % NOTE: This number will be ignored if you choose to extract inputs from associated text file.
            else
                Options.TimeZeroFrameNumber = NaN;
                Options.StartAnalysisFrameNumber = 1;
            end
    
    % ---------Parameters Used To Find Particles/Assess Particle 'Goodness'---------
        Options.MinParticleSize = 4;
            % This is the minimum particle size (defined as the number of connected pixels 
            % above the threshold) in order for a particle to be found. 
            % Particles smaller than this size will not be found (i.e. the program 
            % will assume that they are noise).
        Options.MaxParticleSize = 100; 
            % This is the maximum particle size (defined as the number of connected pixels 
            % above the threshold) in order for a particle to be considered "good". 
            % Particles larger than this size will be designated as "bad".
        Options.MaxEccentricity = 0.8;
            % This is the maximum eccentricity that a particle can have in order to still 
            % be considered "good". 0 = perfect circle, 1 = straight line. If the 
            % eccentricity is too high, that may indicate that the particle being 
            % analyzed is actually two diffraction limited particles close together.

% ===========Options you are unlikely to ever change===============

      
    % ---------Workflow Options---------
        Options.AutoCreateLabels = 'y';
             % 'y' OR 'n'
            % Choose 'y' to automatically use information in the pathname and/or 
            % filename to create labels for the output data file and/or save folder. 
            % If so, you should modify Create_Save_Folder_And_Grab_Data_Labels.m so it 
            % extracts the information that you want.
            
    % --------- Combine Multiple Videos Options---------
        Options.CombineMultipleVideos = 'n';
            % 'y' OR 'n'
            % Choose 'y'if you want to sequentially combine multiple videos (such 
            % as if your time-lapse got interrupted and you had to re-start a second 
            % set). If you do this, you will also need to list the options below
            % NOTE: this is set up assuming that you are extracting inputs from 
            % the text file, and that you are extracting times from the metadata. 
            % It also assumes that you are providing the DefaultPath and SavePath as inputs to Run_Me_To_Start
            % If those situations don't apply, then you'll need to double check that 
            % this works properly
            
         Options.NumberOfVideosToCombine = NaN;
            % This is the number of videos that you will be combining
            
         Options.TimeDelayBetweenVideos = 1.1; % Default units = minutes
            % This is the time that elapsed between the end of the previous video 
            % and the start of the next one. If you have more than 2 videos, then 
            % this will be a vector rather than a single value.
    
    % --------- Legacy Options You are Unlikely to Change---------
        Options.ManuallyCorrectFind = 'n';
            Options.CorrectFindNumber = NaN;  
            Options.DisplayAllFigures = 'y';
        Options.DisplayRejectionReasons = 'y';
    
        Options.UserGrabExampleTrace = 'n';
            % 'y' OR 'n'
            % Choose 'y' if you want to plot an example trace in a separate window. 
            % If you choose 'y', you should modify User_Grab_Example_Trace.m with 
            % the specific details of the trace that you want to plot and how you 
            % want to plot it.
    
        Options.MinROISize = 4; 
        Options.MaxROISize = NaN;
            % These determine the minimum and maximum size allowed for the region 
            % of interest around each particle. In between those values, the ROI 
            % scales with the size of the particle itself (i.e. bigger particles 
            % have larger regions of interest). The ROI is a square, and the 
            % values indicate the number of pixels along one side of the square 
            % (so 5 means a 5x5 pixel ROI). The minimum size must always be 
            % specified. 'NaN' can be used to indicate no maximum size.
    
        Options.FrameNumberLimit = NaN;
            % Determines the number of frames which will be loaded and analyzed. 
            % Use 'NaN' to indicate no limit (i.e. all frames will be included).        
            % You would likely only use this option if you wanted to quickly 
            % assess how the program is running without having 
            % to wait for an entire video to load. Alternatively, you could 
            % use it to exclude frames at the end of the video from your analysis.
    
            % ---------Types Of Analysis Options---------
            Options.UseGaussianIntensity = 'n';
                % 'y' OR 'n'
                % Choose 'y' if you want to use a Gaussian fit to determine the local 
                % background around a viral particle and use that to determine the 
                % background subtracted intensity of the particle. Choose 'n' if 
                % you want to determine the background as an average across the entire 
                % image (i.e. the background value will be the same for all viral 
                % particles within a given frame). Default is 'n'.
            Options.GrabTotalIntensityOnly = 'n';
                % 'y' OR 'n'
                % Choose 'y' if you don't want to analyze each particle individually and 
                % only want to grab the intensity trace of the entire video. Not typical.

end
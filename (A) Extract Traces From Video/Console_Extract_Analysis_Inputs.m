function [Options] = Console_Extract_Analysis_Inputs(Options)

% By Matthew D. Mitchell, Rawle Lab, Williams College (Jan 2026)

% Extract_Analysis_Inputs
% Console-prompt version (no text file).
% Remembers the previous run's inputs and offers to reuse them.
% Outputs are standardized:
%   - TimeZeroFrameNumber, FrameNumToFindParticles, StartAnalysisFrameNumber: scalar double

    saveFile = fullfile(fileparts(mfilename('fullpath')), 'last_console_inputs.mat');

    % ---- Try to load previous inputs ----
    hasPrev = false;
    if isfile(saveFile)
        try
            prev = load(saveFile);
            hasPrev = isfield(prev, 'TimeZeroFrameNumber') && ...
                      isfield(prev, 'FrameNumToFindParticles') && ...
                      isfield(prev, 'StartAnalysisFrameNumber');
        catch
            hasPrev = false;
        end
    end

    if hasPrev
        % Display previous values and offer to reuse them
        fprintf('\n');
        disp('  Previous analysis inputs:')
        fprintf('    1. Time-zero frame     (tzero) :  %s\n', num2str(prev.TimeZeroFrameNumber));
        fprintf('    2. Finding image frame  (find)  :  %s\n', num2str(prev.FrameNumToFindParticles));
        fprintf('    3. Start analysis frame (start) :  %s\n', num2str(prev.StartAnalysisFrameNumber));
        fprintf('\n');

        resp = strtrim(input('  All good with these values? (y to proceed, n to change): ', 's'));

        % Seed current values from previous run
        tzero  = prev.TimeZeroFrameNumber;
        find_  = prev.FrameNumToFindParticles;
        start_ = prev.StartAnalysisFrameNumber;

        if ~strcmpi(resp, 'y')
            confirmed = false;
            while ~confirmed
                % Let user change whichever fields they want
                while true
                    numStr = strtrim(input('  Which number to change? (1-3, or 0 when done): ', 's'));
                    num = str2double(numStr);
                    if num == 0
                        break
                    elseif num == 1
                        tzero  = promptScalarDouble('    New time-zero frame (tzero): ');
                    elseif num == 2
                        find_  = promptScalarDouble('    New finding image frame (find): ');
                    elseif num == 3
                        start_ = promptScalarDouble('    New start analysis frame (start): ');
                    else
                        disp('    Enter 1, 2, 3, or 0 to finish.');
                    end
                end

                % Show the updated values and ask for confirmation
                fprintf('\n');
                disp('====================================')
                disp('  Updated analysis inputs:')
                fprintf('    1. Time-zero frame     (tzero) :  %s\n', num2str(tzero));
                fprintf('    2. Finding image frame  (find)  :  %s\n', num2str(find_));
                fprintf('    3. Start analysis frame (start) :  %s\n', num2str(start_));
                disp('====================================')
                confirm = strtrim(input('  All good? (y to proceed): ', 's'));
                if strcmpi(confirm, 'y')
                    confirmed = true;
                end
            end
        end

    else
        % First run — prompt for all values, offering the standard
        % default sequence (1, 2, 3): tzero=1, find=2, start=3.
        fprintf('\n\n');
        disp('  Please enter the following analysis inputs:')
        tzero  = promptScalarDouble('    Time-zero frame     (tzero) [default: 1]: ');
        find_  = promptScalarDouble('    Finding image frame  (find)  [default: 2]: ');
        start_ = promptScalarDouble('    Start analysis frame (start) [default: 3]: ');
    end

    % ---- Write back to Options ----
    Options.TimeZeroFrameNumber      = tzero;
    Options.FrameNumToFindParticles  = find_;
    Options.StartAnalysisFrameNumber = start_;

    % ---- Save for next run (merge so video path fields are not overwritten) ----
    if isfile(saveFile)
        fileState = load(saveFile);
    else
        fileState = struct();
    end
    fileState.TimeZeroFrameNumber      = tzero;
    fileState.FrameNumToFindParticles  = find_;
    fileState.StartAnalysisFrameNumber = start_;
    save(saveFile, '-struct', 'fileState');

end

% ===================== Local Functions =====================

function val = promptScalarDouble(promptText)
    while true
        s = strtrim(input(promptText, 's'));
        if isempty(s)
            disp('    Enter a single number.');
            continue
        end
        val = str2double(s);
        if isfinite(val) && isscalar(val)
            val = double(val);
            return
        end
        disp('    Invalid input. Enter one numeric value (e.g., 3).');
    end
end

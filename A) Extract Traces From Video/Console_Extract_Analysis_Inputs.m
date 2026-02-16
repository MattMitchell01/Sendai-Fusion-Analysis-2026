function [Options] = Console_Extract_Analysis_Inputs(Options)

% By Matthew D. Mitchell, Rawle Lab, Williams College (Jan 2026)

% Extract_Analysis_Inputs
% Console-prompt version (no text file).
% Outputs are standardized:
%   - TimeZeroFrameNumber, FrameNumToFindParticles, StartAnalysisFrameNumber: scalar double
%   - IgnoreFrameNumbers: row vector double (empty allowed)

    % ---- Required scalar inputs ---
    fprintf('\n\n');
    disp("Please Enter the Following Analysis Inputs into the Console")
    Options.TimeZeroFrameNumber      = promptScalarDouble('Time-zero frame (tzero): ');
    Options.FrameNumToFindParticles  = promptScalarDouble('Finding image frame (find): ');
    Options.StartAnalysisFrameNumber = promptScalarDouble('Start analysis frame (start): ');

    % ---- Ignore list ----
    % User enters: "NaN" or "" for none, OR a space/comma-separated list like "5 6 7" or "5,6,7"
    ignoreStr = input('Ignore frames (space/comma-seperated list, or NaN for none): ', 's');
    Options.IgnoreFrameNumbers = parseNumberListToRowDouble(ignoreStr);
    
end

% ===================== Local Functions =====================

function val = promptScalarDouble(promptText)
    while true
        s = strtrim(input(promptText, 's'));
        if isempty(s)
            disp('  Enter a single number.');
            continue
        end
        val = str2double(s);
        if isfinite(val) && isscalar(val)
            val = double(val);
            return
        end
        disp('  Invalid input. Enter one numeric value (e.g., 3).');
    end
end

function vec = parseNumberListToRowDouble(s)
    s = strtrim(s);

    if isempty(s) || strcmpi(s,'nan')
        vec = double([]);   % no ignored frames
        return
    end

    % Allow commas or spaces
    s = regexprep(s, '[,;]+', ' ');
    raw = sscanf(s, '%f');

    if isempty(raw)
        vec = double([]);
    else
        vec = double(raw(:).'); % row vector
    end
end
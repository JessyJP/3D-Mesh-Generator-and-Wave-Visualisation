function printProgressBar(iteration, total, varargin)
    % printProgressBar(iteration, total, varargin)
    %
    % This function prints a progress bar in the MATLAB command window that
    % overwrites itself with each update, providing a visual indication of
    % progress for long-running operations.
    %
    % Inputs:
    %   iteration - The current iteration number (should not exceed 'total').
    %   total - The total number of iterations.
    %   varargin (optional) - Additional named parameters:
    %       'barLength' - The length of the progress bar in characters (default: 50).
    %       'prefixText' - A prefix text to display before the progress bar (default: 'Progress:').
    %       'refreshRate' - Minimum time in seconds between updates (default: 0.5).
    %
    % Example usage:
    %   N = 100;
    %   for i = 1:N
    %       printProgressBar(i, N, 'barLength', 40, 'prefixText', 'Loading:', 'refreshRate', 0.1);
    %       pause(0.05); % Simulate some work
    %   end

    % Persistent variables to maintain state between function calls
    persistent lastUpdate; % Timestamp of the last progress bar update
    persistent lastLength; % Length of the last printed progress bar string
    
    % Initialize persistent variables on the first function call
    if isempty(lastUpdate)
        lastUpdate = tic; % Start a timer
        lastLength = 0;    % Initialize the length of the last printed string
    end
    
    % Parse optional input arguments using inputParser
    p = inputParser;
    addOptional(p, 'barLength', 50, @isnumeric); % Length of the progress bar
    addOptional(p, 'prefixText', 'Progress:', @ischar); % Prefix text
    addOptional(p, 'refreshRate', 0.5, @isnumeric); % Refresh rate in seconds
    parse(p, varargin{:});
    
    % Extract the values of optional parameters
    barLength = p.Results.barLength;
    prefixText = p.Results.prefixText;
    refreshRate = p.Results.refreshRate;
    
    % Update the progress bar only if the specified refresh time has elapsed or on the last iteration
    if toc(lastUpdate) > refreshRate || iteration == total
        percentComplete = floor(iteration / total * 100); % Calculate percent completion
        barsComplete = floor(iteration / total * barLength); % Calculate the number of '=' to print
        barString = [repmat('=', 1, barsComplete), repmat(' ', 1, barLength - barsComplete)]; % Construct the progress bar string
        
        % Clear the previous output by printing backspaces
        if lastLength > 0
            fprintf(repmat('\b', 1, lastLength));
        end
        
        % Generate and print the new progress bar string
        progressStr = sprintf('\r%s [%s] %3d%%', prefixText, barString, percentComplete);
        fprintf('%s', progressStr);
        
        % Update the state for the next call
        lastLength = length(progressStr); % Remember the length of the printed string
        lastUpdate = tic; % Reset the timer
    end
    
    % At the completion of the process, move to the next line and reset persistent variables
    if iteration >= total
        fprintf('\n');
        clear lastUpdate lastLength;
    end
end

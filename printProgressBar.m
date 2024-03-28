function printProgressBar(iteration, total, varargin)
    % printProgressBar(iteration, total, varargin)
    %
    % This function prints a progress bar in the MATLAB command window that
    % overwrites itself with each update, providing a visual indication of
    % progress for long-running operations. It supports multiple progress bars
    % by using a unique identifier for each bar. If no identifier is provided,
    % a default ID of 0 is used.
    %
    % Inputs:
    %   iteration - The current iteration number (should not exceed 'total').
    %   total - The total number of iterations for the process.
    %   varargin - Additional optional parameters:
    %       'id' - A unique identifier for the progress bar. This can be any value
    %              that can be converted to a string, such as a number or string,
    %              allowing for multiple progress bars to be managed independently.
    %              The default value is 0.
    %       'barLength' - The length of the progress bar in characters (default: 50).
    %       'prefixText' - Text to display before the progress bar (default: 'Progress:').
    %       'refreshRate' - Minimum time in seconds between updates to prevent
    %                       excessive output (default: 0.5 seconds).
    %
    % Example usage:
    %   N = 100;
    %   for i = 1:N
    %       printProgressBar(i, N, 'barLength', 40, 'prefixText', 'Loading:', 'refreshRate', 0.1);
    %       pause(0.05); % Simulating work
    %   end

    % Parse optional input arguments
    p = inputParser;
    addParameter(p, 'id', 0);
    addParameter(p, 'barLength', 50, @isnumeric);
    addParameter(p, 'prefixText', 'Progress:', @ischar);
    addParameter(p, 'refreshRate', 0.5, @isnumeric);
    parse(p, varargin{:});

    id = p.Results.id;
    barLength = p.Results.barLength;
    prefixText = p.Results.prefixText;
    refreshRate = p.Results.refreshRate;

    % Convert id to string to use as a key in the map
    idStr = num2str(id);

    % Persistent variable to store state of each progress bar
    persistent progressStates;

    % Initialize the progressStates map on the first function call
    if isempty(progressStates)
        progressStates = containers.Map('KeyType', 'char', 'ValueType', 'any');
    end

    % Initialize state for new progress bar
    if ~isKey(progressStates, idStr)
        progressStates(idStr) = struct('lastUpdate', tic, 'lastLength', 0);
    end

    % Retrieve the state for the current progress bar
    state = progressStates(idStr);

    % Update the progress bar only if the specified refresh time has elapsed or on the last iteration
    if toc(state.lastUpdate) > refreshRate || iteration == total
        percentComplete = floor(iteration / total * 100);
        barsComplete = floor(iteration / total * barLength);
        barString = [repmat('=', 1, barsComplete), repmat(' ', 1, barLength - barsComplete)];

        % Clear the previous output by printing backspaces
        if state.lastLength > 0
            fprintf(repmat('\b', 1, state.lastLength));
        end

        % Generate and print the new progress bar string
        progressStr = sprintf('\r%s [%s] %3d%%', prefixText, barString, percentComplete);
        fprintf('%s', progressStr);

        % Update the state for the next call
        state.lastLength = length(progressStr);
        state.lastUpdate = tic;
        progressStates(idStr) = state; % Save the updated state back to the map
    end

    % At the completion of the process, move to the next line and reset persistent variables
    if iteration >= total
        fprintf('\n');
        progressStates.remove(idStr); % Remove the completed progress bar state
        if isempty(progressStates.Count)
            clear progressStates; % Clear the map if it's empty to reset state
        end
    end
end

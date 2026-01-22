function [macIndex, macValue] = modal_assurance_criterion(modes, modeRef, varargin)
    % Calculate the Modal Assurance Criterion (MAC) between a set of modes
    % and a reference mode.
    % Inputs:
    %   modes: the set of modes.
    %   modeRef: the reference mode.
    %   doPrint: flag to print the MAC value and index (optional, default
    %       is false).
    % Outputs:
    %   macIndex: the index of the mode with the maximum MAC value.
    %   macValue: the maximum MAC value.

    % Parse inputs
    p = inputParser;
    addOptional(p, 'doPrint', false);
    parse(p, varargin{:});
    doPrint = p.Results.doPrint;

    % Initialize variables
    modeDot = dot(modeRef, modeRef);
    macValues = zeros(size(modes, 2), 1);

    % Loop over modes
    for ii = 1:size(modes, 2)
        thisMode = modes(:, ii);
        macValues(ii) = dot(thisMode, modeRef).^2 / (modeDot * dot(thisMode, thisMode));
    end

    % Find max mac
    [macValue, macIndex] = max(macValues);

    % Print
    if doPrint
        fprintf("\nIdx %d has MAC %.4f\n", macIndex, macValue)
    end
end

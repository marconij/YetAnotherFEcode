function err = error_measure(obj, rho, varargin)
    % Calculate the error measure from the invariance equation.
    % See Pozzi 2024 (10.1007/s11071-024-09881-5) for more details.
    % Inputs:
    %   rho: the reduced amplitude.
    %   style1stOrder: the style of the 1st order system matrices (optional, default is 'L1').
    %       Valid values are 'L1' and 'L2'.
    %       See Jain 2022 (10.1007/s11071-021-06957-4) for more details.
    % Outputs:
    %   err: the error measure.

    % Parse the inputs
    p = inputParser;
    addOptional(p, 'style1stOrder', 'L1', @(x) ischar(x) && ismember(x, {'L1', 'L2'}));
    parse(p, varargin{:});
    style1stOrder = p.Results.style1stOrder;

    % Number of theta values
    nTheta = 30;
    theta = linspace(0, 2*pi, 30);

    % Loop over theta values
    err = 0;
    for ii = 1:nTheta
        % Compute the reduced coordinates
        p = rho * [exp(1i*theta(ii)); exp(-1i*theta(ii))];

        % Compute the residual error of the invariant equation
        [resNorm, norms] = obj.compute_invariance_residual(p, style1stOrder);

        % Update the error measure
        err = err + resNorm / nTheta / max(norms);
    end
end

function [domega, domegaDamped] = adjoint_sensitivity_omega(obj)
    % Compute the adjoint sensitivity of the eigenfrequency omega.
    % Outputs:
    %   domega: the adjoint sensitivity of omega.
    %   domegaDamped: the adjoint sensitivity of the damped omega.

    % Check if the partial derivatives are available
    if isempty(obj.pM)
        error(['The derivatives of the structural matrices are not available. ', ...
            'Please set the `computePartialDerivatives` flag to true when building the system.']);
    end
    if isempty(obj.pLambda_omega)
        error(['The partial deivatives of the modal analysis are not available. ', ...
            'Please set the `computePartialDerivatives` flag to true for the modal analysis.']);
    end

    % Initialize
    domega = zeros(obj.ndv, 1);
    domegaDamped = zeros(obj.ndv, 1);

    % Loop over the design variables
    for dv = 1:obj.ndv
        % Sensitivity of the frequency
        domega(dv) = obj.phi.' * (-obj.omega^2 * obj.pM(dv).coeffs + obj.pK(dv).coeffs) * obj.phi / (2 * obj.omega);

        % Sensitivity of the damped frequency
        domegaDamped(dv) = (sqrt(1 - obj.csi^2) - obj.omega * obj.csi / sqrt(1 - obj.csi^2) * obj.pcsi_omega) * domega(dv);
    end
end

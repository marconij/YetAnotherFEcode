function [domega, domegaDamped] = adjoint_sensitivity_omega_topopt(obj)
    % Compute the adjoint sensitivity of the eigenfrequency omega for a
    % topology optimization problem.
    % Outputs:
    %   domega: the adjoint sensitivity of omega.
    %   domegaDamped: the adjoint sensitivity of the damped omega.

    % Check if the partial derivatives are available
    if isempty(obj.Me)
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

    % Convert the mode shape into its unconstrained form
    phiUnc = obj.yafecAssembly.unconstrain_vector(obj.phi);

    % Loop over the design variables
    for dv = 1:obj.ndv
        % Extract element quantities
        thisElement = obj.yafecAssembly.Mesh.Elements(dv).Object;
        phiE = phiUnc(thisElement.iDOFs);

        % Sensitivity of the frequency
        domega(dv) = phiE.' * (-obj.omega^2 * obj.Me + obj.Ke) * phiE / (2 * obj.omega);

        % Sensitivity of the damped frequency
        domegaDamped(dv) = (sqrt(1 - obj.csi^2) - obj.omega * obj.csi / sqrt(1 - obj.csi^2) * obj.pcsi_omega) * domega(dv);
    end
end

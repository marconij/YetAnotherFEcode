function modal_analysis_sensitivity(obj)
    % Compute the direct sensitivity of the eigenvalues and mode shapes
    % with respect to the system parameters.
    % This method must be called before computing the direct sensitivity
    % of the SSM.

    % Check if the partial derivatives are available
    if isempty(obj.pM)
        error(['The derivatives of the structural matrices are not available. ', ...
            'Please set the `computePartialDerivatives` flag to true when building the system.']);
    end
    if isempty(obj.pLambda_omega)
        error(['The partial deivatives of the modal analysis are not available. ', ...
            'Please set the `computePartialDerivatives` flag to true for the modal analysis.']);
    end

    % Initialize sensitivity
    obj.dphi = repmat(struct('coeffs', []), 1, obj.ndv);
    obj.dPhi = repmat(struct('coeffs', []), 1, obj.ndv);
    obj.domega = zeros(obj.ndv, 1);
    obj.dcsi = zeros(obj.ndv, 1);
    obj.ddelta = zeros(obj.ndv, 1);
    obj.domegaDamped = zeros(obj.ndv, 1);
    obj.dlambda = zeros(obj.ndv, 1);
    obj.dLambda = zeros(obj.ndv, 2);
    obj.dC = repmat(struct('coeffs', []), 1, obj.ndv);

    % Build matrix
    A = [obj.K - obj.omega^2*obj.M, -2*obj.omega*obj.M*obj.phi;
         -2*obj.omega*obj.phi.'*obj.M, 0];

    % Build right-hand sides
    b = zeros(obj.n + 1, obj.ndv);
    for dv = 1:obj.ndv
        b(:, dv) = [-(obj.pK(dv).coeffs - obj.omega^2*obj.pM(dv).coeffs)*obj.phi;
                    obj.omega*obj.phi.'*obj.pM(dv).coeffs*obj.phi];
    end

    % Solve linear system
    X = A\b;

    % Extract sensitivity
    for dv = 1:obj.ndv
        % Mode shape
        obj.dphi(dv).coeffs = X(1:end-1, dv);
        obj.dPhi(dv).coeffs = [obj.dphi(dv).coeffs, obj.dphi(dv).coeffs];

        % Eigenfrequency
        obj.domega(dv) = X(end, dv);
        obj.dcsi(dv) = obj.pcsi_omega * obj.domega(dv);
        obj.ddelta(dv) = 2 * obj.omega * obj.dcsi(dv) + 2 * obj.csi * obj.domega(dv);
        obj.domegaDamped(dv) = sqrt(1 - obj.csi^2) * obj.domega(dv) - ...
            obj.omega * obj.csi / sqrt(1 - obj.csi^2) * obj.dcsi(dv);

        % Eigenvalue
        obj.dlambda(dv) = obj.plambda_omega * obj.domega(dv);
        obj.dLambda(dv, :) = [obj.dlambda(dv), conj(obj.dlambda(dv))];

        % Total derivative of the damping matrix
        obj.dC(dv).coeffs = obj.pC(dv).coeffs + obj.pC_omega * obj.domega(dv);
    end
end

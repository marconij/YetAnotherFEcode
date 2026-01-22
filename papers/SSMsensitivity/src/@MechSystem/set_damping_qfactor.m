function set_damping_qfactor(obj, Qfactor, varargin)
    % Set damping using Q-factor.
    % Inputs:
    %   Qfactor: Q-factor value for the damping.
    %   computePartialDerivatives: flag to enable the computation of the
    %       partial derivatives (optional, default is false).
    %   rayleigh: choose between "alpha" (mass-proportional damping) and
    %       "beta" (stiffness-proportional damping). Optional,
    %       default is "alpha".
    %   topOptProblem: flag to enable the topology optimization (optional,
    %       default is false).

    % Parse the inputs
    p = inputParser;
    addOptional(p, 'computePartialDerivatives', false);
    addOptional(p, 'rayleigh', 'alpha');
    addOptional(p, 'topOptProblem', false);
    parse(p, varargin{:});
    computePartialDerivatives = p.Results.computePartialDerivatives;
    rayleigh = p.Results.rayleigh;
    topOptProblem = p.Results.topOptProblem;

    % Check if the modal analysis has been performed
    if isempty(obj.omega)
        error('The modal analysis must be performed before setting the damping.');
    end

    % Set damping type
    obj.damping_type = 2;
    obj.Qfactor = Qfactor;

    % Damping ratio and modal damping
    obj.csi = 1 / (2 * Qfactor);
    obj.delta = 2 * obj.omega * obj.csi;

    % Compute the corresponding Rayleigh coefficients
    switch rayleigh
        case 'alpha'
            obj.rayleigh(1) = obj.omega / Qfactor;
            obj.rayleigh(2) = 0;
            % Damping matrix 
            obj.C = obj.omega / Qfactor * obj.M;
        case 'beta'
            obj.rayleigh(1) = 0;
            obj.rayleigh(2) = 1 / (obj.omega * Qfactor);
            % Damping matrix 
            obj.C = obj.rayleigh(2) * obj.K;
        otherwise
            error('choose either "alpha" or "beta"!')
    end

    % Damped eigenfrequency
    obj.omegaDamped = obj.omega * sqrt(1 - obj.csi^2);

    % Complex eigenvalue
    obj.lambda = -obj.omega * obj.csi + 1i * obj.omega * sqrt(1 - obj.csi^2);

    % Mass-normalized left mode shape
    obj.tau = 1 / (1i * 2 * obj.omega * sqrt(1 - obj.csi^2)); % mass-normalization coefficient
    obj.psi = conj(obj.tau) * obj.phi;

    % Master subspace
    obj.Lambda = [obj.lambda; conj(obj.lambda)]; % (2 x 1) complex eigenvalues
    obj.Phi = [obj.phi, obj.phi]; % (n x 2) mass-normalized right mode shapes
    obj.Psi = [obj.psi, conj(obj.psi)]; % (n x 2) mass-normalized left mode shapes

    % Partial derivatives
    if computePartialDerivatives
        % Partial derivatives of the damping matrix
        switch rayleigh
            case 'alpha'
                obj.pC_omega = obj.M / Qfactor;
            case 'beta'
                obj.pC_omega = - obj.K / (Qfactor * obj.omega^2);
        end
        if topOptProblem
            obj.Ce = obj.rayleigh(1) * obj.Me + obj.rayleigh(2) * obj.Ke;
        else
            obj.pC = repmat(struct('coeffs', []), obj.ndv, 1);
            for dv = 1:obj.ndv
                obj.pC(dv).coeffs = obj.rayleigh(1) * obj.pM(dv).coeffs + obj.rayleigh(2) * obj.pK(dv).coeffs;
            end
        end

        % Partial derivatives of the eigenvalues
        obj.pcsi_omega = 0;
        obj.pdelta_omega = 2 * obj.csi;
        obj.plambda_omega = -obj.csi + 1i * sqrt(1 - obj.csi^2);
        obj.pLambda_omega = [obj.plambda_omega; conj(obj.plambda_omega)];
    end
end

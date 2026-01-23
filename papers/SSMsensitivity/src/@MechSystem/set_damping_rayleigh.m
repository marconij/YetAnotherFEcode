function set_damping_rayleigh(obj, rayleigh, varargin)
    % Set the Rayleigh damping coefficients.
    % Inputs:
    %   rayleigh: vector of Rayleigh damping coefficients [alpha, beta].
    %   computePartialDerivatives: flag to enable the computation of the
    %       partial derivatives (optional, default is false).
    %   topOptProblem: flag to enable the topology optimization (optional,
    %       default is false).

    % Parse the inputs
    p = inputParser;
    addOptional(p, 'computePartialDerivatives', false);
    addOptional(p, 'topOptProblem', false);
    parse(p, varargin{:});
    computePartialDerivatives = p.Results.computePartialDerivatives;
    topOptProblem = p.Results.topOptProblem;

    % Check if the modal analysis has been performed
    if isempty(obj.omega)
        error('The modal analysis must be performed before setting the damping.');
    end

    % Set damping type
    obj.damping_type = 1;
    obj.rayleigh = rayleigh;

    % Damping ratio and modal damping
    obj.csi = (rayleigh(1) + rayleigh(2) * obj.omega^2) / (2 * obj.omega);
    obj.delta = 2 * obj.omega * obj.csi;
    obj.Qfactor = 1 / (2 * obj.csi);

    % Damping matrix
    obj.C = obj.rayleigh(1) * obj.M + obj.rayleigh(2) * obj.K;

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
        obj.pC_omega = sparse(obj.n, obj.n);
        if topOptProblem
            obj.Ce = obj.rayleigh(1) * obj.Me + obj.rayleigh(2) * obj.Ke;
        else
            obj.pC = repmat(struct('coeffs', []), obj.ndv, 1);
            for dv = 1:obj.ndv
                obj.pC(dv).coeffs = obj.rayleigh(1) * obj.pM(dv).coeffs + obj.rayleigh(2) * obj.pK(dv).coeffs;
            end
        end

        % Partial derivatives of the eigenvalues
        obj.pcsi_omega = (obj.rayleigh(2) * obj.omega^2 - obj.rayleigh(1)) / (2 * obj.omega^2);
        obj.pdelta_omega = 2 * obj.csi + 2 * obj.omega * obj.pcsi_omega;
        obj.plambda_omega = -obj.csi + 1i * sqrt(1 - obj.csi^2) + ...
            (-obj.omega - 1i * obj.csi * obj.omega / sqrt(1 - obj.csi^2)) * obj.pcsi_omega;
        obj.pLambda_omega = [obj.plambda_omega; conj(obj.plambda_omega)];
    end
end

function [D, V] = modal_analysis(obj, varargin)
    % Perform the linear eigenvalue analysis of the mechanical system.
    % Inputs:
    %   nModes: number of modes to compute (optional, default is 10).
    %   modeIndex: index of the mode to extract (optional, default is 1).
    %       This value is overwritten if the reference mode shape is
    %       provided.
    %   modeRef: reference mode shape for the Modal Assurance Criterion
    %       (optional, default is []). If the reference mode shape is
    %       provided, the modeIndex is computed as the mode that has the
    %       highest MAC value the reference mode shape.
    %   computePartialDerivatives: flag to enable the computation of the
    %       partial derivatives (optional, default is false).
    % Outputs:
    %   D: the diagonal matrix with the eigenvalues (nModes x nModes).
    %   V: the matrix with the eigenvectors (n x nModes).

    % Parse the inputs
    p = inputParser;
    addOptional(p, 'nModes', 10, @isnumeric);
    addOptional(p, 'modeIndex', 1, @isnumeric);
    addOptional(p, 'modeRef', []);
    addOptional(p, 'computePartialDerivatives', false);
    parse(p, varargin{:});
    nModes = p.Results.nModes;
    modeIndex = p.Results.modeIndex;
    modeRef = p.Results.modeRef;
    computePartialDerivatives = p.Results.computePartialDerivatives;

    % Linear eigenvalue problem
    [V, D] = eigs(obj.K, obj.M, nModes, 'sm');
    omega2 = diag(D);

    % Check if the reference mode shape is provided
    if ~isempty(modeRef)
        % Compute the MAC
        modeIndex = modal_assurance_criterion(V, modeRef);
    end

    % Extract the mode
    obj.omega = sqrt(omega2(modeIndex));
    obj.phi = V(:, modeIndex);
    obj.phi = obj.phi ./ sqrt(obj.phi.' * obj.M * obj.phi); % mass-normalized mode shape

    % No damping condition (default)
    % Set the damping with the set_damping method
    obj.csi = 0;
    obj.Qfactor = Inf;
    obj.rayleigh = [0, 0];
    obj.delta = 0;
    % obj.csi = (obj.rayleigh(1) + obj.rayleigh(2) * obj.omega^2) / (2 * obj.omega);
    % obj.delta = 2 * obj.omega * obj.csi;
    % pcsi_omega = (obj.rayleigh(2) * obj.omega^2 - obj.rayleigh(1)) / (2 * obj.omega^2);

    % Complex eigenvalue
    % obj.lambda = -obj.omega * obj.csi + 1i * obj.omega * sqrt(1 - obj.csi^2);
    obj.lambda = 1i * obj.omega;

    % Mass-normalized left mode shape
    % obj.tau = 1 / (1i * 2 * obj.omega * sqrt(1 - obj.csi^2)); % mass-normalization coefficient
    obj.tau = 1 / (2i * obj.omega);
    obj.psi = conj(obj.tau) * obj.phi;

    % Damped eigenfrequency
    % obj.omegaDamped = obj.omega * sqrt(1 - obj.csi^2);
    obj.omegaDamped = obj.omega;

    % Master subspace
    obj.Lambda = [obj.lambda; conj(obj.lambda)]; % (2 x 1) complex eigenvalues
    obj.Phi = [obj.phi, obj.phi]; % (n x 2) mass-normalized right mode shapes
    obj.Psi = [obj.psi, conj(obj.psi)]; % (n x 2) mass-normalized left mode shapes

    % Partial derivatives of the eigenvalues
    if computePartialDerivatives
        % obj.pcsi_omega = pcsi_omega;
        % obj.pdelta_omega = 2 * obj.csi + 2 * obj.omega * obj.pcsi_omega;
        % obj.plambda_omega = -obj.csi + 1i * sqrt(1 - obj.csi^2) + ...
        %     (-obj.omega - 1i * obj.csi * obj.omega / sqrt(1 - obj.csi^2)) * obj.pcsi_omega;
        % obj.pLambda_omega = [obj.plambda_omega; conj(obj.plambda_omega)];
        obj.pC_omega = sparse(obj.n, obj.n);
        obj.pcsi_omega = 0;
        obj.pdelta_omega = 0;
        obj.plambda_omega = 1i;
        obj.pLambda_omega = [obj.plambda_omega; conj(obj.plambda_omega)];
    end
end

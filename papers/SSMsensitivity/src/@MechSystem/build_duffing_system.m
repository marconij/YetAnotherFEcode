function build_duffing_system(obj, m, k, k2, k3, varargin)
    % Build the Duffing oscillator system.
    % Inputs:
    %   m: mass of each mass.
    %   k: stiffness coefficient.
    %   k2: quadratic stiffness coefficient.
    %   k3: cubic stiffness coefficient.
    %   computePartialDerivatives: flag to enable the computation of the
    %       partial derivatives (optional, default is false).

    % Parse the inputs
    p = inputParser;
    addOptional(p, 'computePartialDerivatives', false);
    parse(p, varargin{:});
    computePartialDerivatives = p.Results.computePartialDerivatives;

    % Set the number of DOFs
    obj.n = 1;

    % Create mass matrix
    obj.M = m * speye(obj.n);

    % Create the sparse matrices and tensors
    Ktilde = sparse(1, 1, 1, obj.n, obj.n);
    T2tilde = sptensor(ones(1, 3), 1, [obj.n, obj.n, obj.n]);
    T3tilde = sptensor(ones(1, 4), 1, [obj.n, obj.n, obj.n, obj.n]);

    % Create the stiffness matrix and tensors
    obj.K = k * Ktilde;
    obj.T2 = k2 * T2tilde;
    obj.T3 = k3 * T3tilde;

    % Damping matrix
    obj.C = sparse(obj.n, obj.n);
    
    % Number of design variables
    obj.ndv = 4;

    % Compute the partial derivatives
    if computePartialDerivatives
        % Mass matrix partial derivatives
        obj.pM = repmat(struct('coeffs', []), obj.ndv, 1);
        obj.pM(1).coeffs = speye(obj.n);
        obj.pM(2).coeffs = sparse(obj.n, obj.n);
        obj.pM(3).coeffs = sparse(obj.n, obj.n);
        obj.pM(4).coeffs = sparse(obj.n, obj.n);

        % Damping matrix partial derivatives
        obj.pC = repmat(struct('coeffs', sparse(obj.n, obj.n)), obj.ndv, 1);

        % Stiffness matrix partial derivatives
        obj.pK = repmat(struct('coeffs', []), obj.ndv, 1);
        obj.pK(1).coeffs = sparse(obj.n, obj.n);
        obj.pK(2).coeffs = Ktilde;
        obj.pK(3).coeffs = sparse(obj.n, obj.n);
        obj.pK(4).coeffs = sparse(obj.n, obj.n);

        % Quadratic stiffness tensor partial derivatives
        obj.pT2 = repmat(struct('coeffs', []), obj.ndv, 1);
        obj.pT2(1).coeffs = sptensor([obj.n, obj.n, obj.n]);
        obj.pT2(2).coeffs = sptensor([obj.n, obj.n, obj.n]);
        obj.pT2(3).coeffs = T2tilde;
        obj.pT2(4).coeffs = sptensor([obj.n, obj.n, obj.n]);

        % Cubic stiffness tensor partial derivatives
        obj.pT3 = repmat(struct('coeffs', []), obj.ndv, 1);
        obj.pT3(1).coeffs = sptensor([obj.n, obj.n, obj.n, obj.n]);
        obj.pT3(2).coeffs = sptensor([obj.n, obj.n, obj.n, obj.n]);
        obj.pT3(3).coeffs = sptensor([obj.n, obj.n, obj.n, obj.n]);
        obj.pT3(4).coeffs = T3tilde;
    end
end

function build_shaw_pierre_system(obj, n, m, k, k2, k3, idxLumped, varargin)
    % This function builds the modified Shaw-Pierre system described in
    % Ponsioen et al. 2019. https://doi.org/10.1007/s11071-019-05023-4
    % Inputs:
    %   n: number of masses.
    %   m: mass of each mass.
    %   k: stiffness coefficient.
    %   k2: quadratic stiffness coefficient.
    %   k3: cubic stiffness coefficient.
    %   idxLumped: index of the mass with the lumped nonlinearities.
    %   computePartialDerivatives: flag to enable the computation of the
    %       partial derivatives (optional, default is false).
    %   freeEnd: flag to leave a free end (optional, default is false).

    % Parse the inputs
    p = inputParser;
    addOptional(p, 'computePartialDerivatives', false);
    addOptional(p, 'freeEnd', false);
    parse(p, varargin{:});
    computePartialDerivatives = p.Results.computePartialDerivatives;
    freeEnd = p.Results.freeEnd;

    % Set the number of DOFs
    if freeEnd
        obj.n = n+1;
    else
        obj.n = n;
    end

    % Create mass matrix
    obj.M = m*speye(obj.n);

    % Set the coefficients at element level in a oscillator chain
    
    % Linear
    subs = [1, 1;
            1, 2;
            2, 1;
            2, 2];
    vals = [1; -1; -1; 1];

    % % Quadratic (x2 - x1)^2
    % subs2 = [1, 1, 1;
    %          1, 1, 2;
    %          1, 2, 1
    %          1, 2, 2;
    %          2, 1, 1;
    %          2, 1, 2;
    %          2, 2, 1;
    %          2, 2, 2];
    % vals2 = [-1; 1; 1;-1; 1; -1; -1; 1];

    % % Cubic (x2 - x1)^3
    % subs3 = [ones(8,1), subs2; 2*ones(8,1), subs2];
    % vals3 = [-vals2; vals2]; 

    % Initialize the sparse matrices and tensors
    n_el = n+1; % number of elements (springs)
    N = n+2; % total number of DOFs including Dirichlet DOFs
    SUBS = zeros(n_el*4, 2);
    VALS = zeros(n_el*4, 1);
    % SUBS2 = zeros(n_el*8, 3);
    % VALS2 = zeros(n_el*8, 1);
    % SUBS3 = zeros(n_el*16, 4);
    % VALS3 = zeros(n_el*16, 1);

    % Perform sparse assembly
    for j = 1:n_el
        % second order matrices
        index = j-1; % zero index for the j-th element
        SUBS(4*(j-1)+1:4*j, :) = repmat(index,[4,2]) + subs;
        VALS(4*(j-1)+1:4*j, :) = vals;
        % % second-order nonlinearity     
        % SUBS2(8*(j-1)+1:8*j, :) = repmat(index,[8,3]) + subs2;
        % VALS2(8*(j-1)+1:8*j, :) = vals2;
        % % third-order nonlinearity
        % SUBS3(16*(j-1)+1:16*j, :) = repmat(index,[16,4]) + subs3;
        % VALS3(16*(j-1)+1:16*j, :) = vals3;
    end

    % Create the sparse matrices and tensors
    Ktilde = sparse(SUBS(:, 1), SUBS(:, 2), VALS, N, N);
    T2tilde = sptensor([idxLumped, idxLumped, idxLumped], 1, [n, n, n]);
    T3tilde = sptensor([idxLumped, idxLumped, idxLumped, idxLumped], 1, [n, n, n, n]);
    % T2tilde = sptensor(SUBS2, VALS2, [N, N, N]);
    % T3tilde = sptensor(SUBS3, VALS3, [N, N, N, N]);

    % Apply the boundary conditions (first and last mass are fixed)
    if freeEnd
        Ktilde = Ktilde(2:end, 2:end);
        % T2tilde = T2tilde(2:end, 2:end, 2:end);
        % T3tilde = T3tilde(2:end, 2:end, 2:end, 2:end);
    else
        Ktilde = Ktilde(2:end-1, 2:end-1);
        % T2tilde = T2tilde(2:end-1, 2:end-1, 2:end-1);
        % T3tilde = T3tilde(2:end-1, 2:end-1, 2:end-1, 2:end-1);
    end

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

        % Stiffness matrix partial derivatives
        obj.pK = repmat(struct('coeffs', []), obj.ndv, 1);
        obj.pK(1).coeffs = sparse(obj.n, obj.n);
        obj.pK(2).coeffs = Ktilde;
        obj.pK(3).coeffs = sparse(obj.n, obj.n);
        obj.pK(4).coeffs = sparse(obj.n, obj.n);

        % Damping matrix partial derivatives
        obj.pC = repmat(struct('coeffs', sparse(obj.n, obj.n)), obj.ndv, 1);

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

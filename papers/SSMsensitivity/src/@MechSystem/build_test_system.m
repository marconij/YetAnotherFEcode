function build_test_system(obj, n, mVec, kVec, k2Vec, k3Vec, varargin)
    % Build the test system. This system is a uniform spring-mass chain
    % with properties that depend on multiple parameters:
    % - m = m_1 + ... + m_n
    % - k = k_1 + ... + k_n
    % - k2 = k2_1 + ... + k2_n
    % - k3 = k3_1 + ... + k3_n
    % Inputs:
    %   n: number of masses.
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
    obj.n = n;

    % Sum the masses and stiffnesses
    m = sum(mVec);
    k = sum(kVec);
    k2 = sum(k2Vec);
    k3 = sum(k3Vec);

    % Create mass matrix
    obj.M = m*speye(n);

    % Set the coefficients at element level in a oscillator chain
    
    % Linear
    subs = [1, 1;
            1, 2;
            2, 1;
            2, 2];
    vals = [1; -1; -1; 1];

    % Quadratic (x2 - x1)^2
    subs2 = [1, 1, 1;
             1, 1, 2;
             1, 2, 1
             1, 2, 2;
             2, 1, 1;
             2, 1, 2;
             2, 2, 1;
             2, 2, 2];
    vals2 = [-1; 1; 1;-1; 1; -1; -1; 1];

    % Cubic (x2 - x1)^3
    subs3 = [ones(8,1), subs2; 2*ones(8,1), subs2];
    vals3 = [-vals2; vals2]; 

    % Initialize the sparse matrices and tensors
    n_el = n+1; % number of elements
    N = n+2; % total number of DOFs including Dirichlet DOFs
    SUBS = zeros(n_el*4, 2);
    VALS = zeros(n_el*4, 1);
    SUBS2 = zeros(n_el*8, 3);
    VALS2 = zeros(n_el*8, 1);
    SUBS3 = zeros(n_el*16, 4);
    VALS3 = zeros(n_el*16, 1);

    % Perform sparse assembly
    for j = 1:n_el
        % second order matrices
        index = j-1; % zero index for the j-th element
        SUBS(4*(j-1)+1:4*j, :) = repmat(index,[4,2]) + subs;
        VALS(4*(j-1)+1:4*j, :) = vals;
        % second-order nonlinearity     
        SUBS2(8*(j-1)+1:8*j, :) = repmat(index,[8,3]) + subs2;
        VALS2(8*(j-1)+1:8*j, :) = vals2;
        % third-order nonlinearity
        SUBS3(16*(j-1)+1:16*j, :) = repmat(index,[16,4]) + subs3;
        VALS3(16*(j-1)+1:16*j, :) = vals3;
    end

    % Create the sparse matrices and tensors
    Ktilde = sparse(SUBS(:, 1), SUBS(:, 2), VALS, N, N);
    T2tilde = sptensor(SUBS2, VALS2, [N, N, N]);
    T3tilde = sptensor(SUBS3, VALS3, [N, N, N, N]);

    % Apply the boundary conditions (first and last mass are fixed)
    Ktilde = Ktilde(2:end-1, 2:end-1);
    T2tilde = T2tilde(2:end-1, 2:end-1, 2:end-1);
    T3tilde = T3tilde(2:end-1, 2:end-1, 2:end-1, 2:end-1);

    % Create the stiffness matrix and tensors
    obj.K = k * Ktilde;
    obj.T2 = k2 * T2tilde;
    obj.T3 = k3 * T3tilde;

    % Damping matrix
    obj.C = sparse(obj.n, obj.n);
    
    % Number of design variables
    obj.ndv = length(mVec) + length(kVec) + length(k2Vec) + length(k3Vec);

    % Compute the partial derivatives
    if computePartialDerivatives
        % Initialize the partial derivatives
        obj.pM = repmat(struct('coeffs', sparse(n, n)), obj.ndv, 1);
        obj.pC = repmat(struct('coeffs', sparse(obj.n, obj.n)), obj.ndv, 1);
        obj.pK = repmat(struct('coeffs', sparse(n, n)), obj.ndv, 1);
        obj.pT2 = repmat(struct('coeffs', sptensor([n, n, n])), obj.ndv, 1);
        obj.pT3 = repmat(struct('coeffs', sptensor([n, n, n, n])), obj.ndv, 1);

        % Contribution of m
        for i = 1:length(mVec)
            obj.pM(i).coeffs = speye(n);
        end

        % Contribution of k
        for i = length(mVec) + 1:length(mVec) + length(kVec)
            obj.pK(i).coeffs = Ktilde;
        end

        % Contribution of k2
        for i = length(mVec) + length(kVec) + 1:length(mVec) + length(kVec) + length(k2Vec)
            obj.pT2(i).coeffs = T2tilde;
        end

        % Contribution of k3
        for i = length(mVec) + length(kVec) + length(k2Vec) + 1:obj.ndv
            obj.pT3(i).coeffs = T3tilde;
        end
    end
end

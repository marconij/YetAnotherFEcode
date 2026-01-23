function build_spring_mass_chain(obj, n, mu, varargin)
    % Build the spring-mass chain system.
    % Inputs:
    %   n: number of masses.
    %   mu: parameters of the system. The first n elements are the masses.
    %       The next n+1 elements are the linear stiffness coefficients.
    %       The next n+1 elements are the quadratic stiffness coefficients.
    %       The next n+1 elements are the cubic stiffness coefficients.
    %   computePartialDerivatives: flag to enable the computation of the
    %       partial derivatives (optional, default is false).

    % Parse the inputs
    p = inputParser;
    addOptional(p, 'computePartialDerivatives', false);
    parse(p, varargin{:});
    computePartialDerivatives = p.Results.computePartialDerivatives;

    % Set the number of DOFs
    obj.n = n;

    % Extract the parameters
    m = mu(1:n);
    k = mu(n+1:2*n+1);
    k2 = mu(2*n+2:3*n+2);
    k3 = mu(3*n+3:4*n+3);

    % Create mass matrix
    obj.M = sparse(diag(m));

    % Set the coefficients at element level in a oscillator chain

    % Linear stiffness subs
    subs = [1, 1;
            1, 2;
            2, 1;
            2, 2];
    vals = [1; -1; -1; 1];

    % Quadratic (x2 - x1)^2
    subs2 = [1, 1, 1;
             1, 1, 2;
             1, 2, 1;
             1, 2, 2;
             2, 1, 1;
             2, 1, 2;
             2, 2, 1;
             2, 2, 2];
    vals2 = [-1; 1; 1; -1; 1; -1; -1; 1];

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
        VALS(4*(j-1)+1:4*j, :) = k(j) * vals;
        % second-order nonlinearity     
        SUBS2(8*(j-1)+1:8*j, :) = repmat(index,[8,3]) + subs2;
        VALS2(8*(j-1)+1:8*j, :) = k2(j) * vals2;
        % third-order nonlinearity
        SUBS3(16*(j-1)+1:16*j, :) = repmat(index,[16,4]) + subs3;
        VALS3(16*(j-1)+1:16*j, :) = k3(j) * vals3;
    end

    % Create the sparse matrices and tensors
    obj.K = sparse(SUBS(:, 1), SUBS(:, 2), VALS, N, N);
    obj.T2 = sptensor(SUBS2, VALS2, [N, N, N]);
    obj.T3 = sptensor(SUBS3, VALS3, [N, N, N, N]);

    %{
    % Test the nonlinearity
    
    % Random displacement
    w = rand(N, 1);
    % w(1) = 0; w(end) = 0;

    % Compute the quadratic nonlinearity
    f2 = zeros(size(w));
    f2(1) = -k2(1) * (w(2) - w(1))^2;
    for i = 2:length(w)-1
        f2(i) = k2(i-1) * (w(i) - w(i-1))^2 - k2(i) * (w(i+1) - w(i))^2;
    end
    f2(end) = k2(end) * (w(end) - w(end-1))^2;

    % Compute the cubic nonlinearity
    f3 = zeros(size(w));
    f3(1) = -k3(1) * (w(2) - w(1))^3;
    for i = 2:length(w)-1
        f3(i) = k3(i-1) * (w(i) - w(i-1))^3 - k3(i) * (w(i+1) - w(i))^3;
    end
    f3(end) = k3(end) * (w(end) - w(end-1))^3;

    % Compute the tensor contraction
    f2t = ttv(obj.T2, {w, w}, -1);
    f3t = ttv(obj.T3, {w, w, w}, -1);

    % Errors
    norm(f2 - f2t.data)
    norm(f3 - f3t.data)
    %}

    % Apply the boundary conditions (first and last mass are fixed)
    obj.K  = obj.K(2:end-1, 2:end-1);
    obj.T2 = obj.T2(2:end-1, 2:end-1, 2:end-1);
    obj.T3 = obj.T3(2:end-1, 2:end-1, 2:end-1, 2:end-1);

    % Damping matrix
    obj.C = sparse(obj.n, obj.n);
    
    % Number of design variables
    obj.ndv = length(mu);

    % Compute the partial derivatives
    if computePartialDerivatives
        % Initialize partial derivatives
        obj.pM = repmat(struct('coeffs', sparse(n, n)), 1, obj.ndv);
        obj.pK = repmat(struct('coeffs', sparse(n, n)), 1, obj.ndv);
        obj.pC = repmat(struct('coeffs', sparse(n, n)), 1, obj.ndv);
        obj.pT2 = repmat(struct('coeffs', sptensor([n, n, n])), 1, obj.ndv);
        obj.pT3 = repmat(struct('coeffs', sptensor([n, n, n, n])), 1, obj.ndv);

        % Loop over masses
        for dv = 1:n
            % Mass matrix
            obj.pM(dv).coeffs = sparse(dv, dv, 1, n, n);
        end

        % Loop over springs
        for i = 1:n_el
            % Zero index for the i-th element
            index = i-1;

            % Stiffness matrix
            dvK = n + i;
            pSUBS = repmat(index,[4,2]) + subs;
            pVALS = vals;
            obj.pK(dvK).coeffs = sparse(pSUBS(:, 1), pSUBS(:, 2), pVALS, N, N);
            obj.pK(dvK).coeffs = obj.pK(dvK).coeffs(2:end-1, 2:end-1);

            % Quadratic nonlinearity
            dvT2 = n + n_el + i;
            pSUBS2 = repmat(index,[8,3]) + subs2;
            pVALS2 = vals2;
            obj.pT2(dvT2).coeffs = sptensor(pSUBS2, pVALS2, [N, N, N]);
            obj.pT2(dvT2).coeffs = obj.pT2(dvT2).coeffs(2:end-1, 2:end-1, 2:end-1);

            % Cubic nonlinearity
            dvT3 = n + 2*n_el + i;
            pSUBS3 = repmat(index,[16,4]) + subs3;
            pVALS3 = vals3;
            obj.pT3(dvT3).coeffs = sptensor(pSUBS3, pVALS3, [N, N, N, N]);
            obj.pT3(dvT3).coeffs = obj.pT3(dvT3).coeffs(2:end-1, 2:end-1, 2:end-1, 2:end-1);
        end

        % % Finite differences
        % h = 1e-6;
        % for dv = 1:length(mu)
        %     muPert = mu;
        %     muPert(dv) = muPert(dv) + h;
        %     sysPert = MechSystem();
        %     sysPert = sysPert.build_spring_mass_chain(n, muPert);
        % 
        %     % Finite differences
        %     obj.pM(dv).coeffs = (sysPert.M - obj.M) / h;
        %     obj.pK(dv).coeffs = (sysPert.K - obj.K) / h;
        %     obj.pT2(dv).coeffs = (sysPert.T2 - obj.T2) / h;
        %     obj.pT3(dv).coeffs = (sysPert.T3 - obj.T3) / h;
        % end
    end
end

function build_yafec_mpc_model(obj, mu, yafec_mpc_constructor, varargin)
    % Build the yafec model representing the mechanical system.
    % Inputs:
    %   mu: the parameter vector.
    %   yafec_mpc_constructor: function that produces the yafec mpc
    %       assembly having as input a parameter vector mu.
    %   computePartialDerivatives: flag to enable the computation of the
    %       partial derivatives (optional, default is false).

    % Parse the inputs
    p = inputParser;
    addOptional(p, 'computePartialDerivatives', false);
    parse(p, varargin{:});
    computePartialDerivatives = p.Results.computePartialDerivatives;

    % Create the assembly
    obj.yafecAssembly = feval(yafec_mpc_constructor, mu);

    % Assmeble mass matrix
    obj.M = obj.yafecAssembly.mass_matrix();

    % Assemble stiffness matrix
    nDOFs = obj.yafecAssembly.Mesh.nDOFs;
    obj.K = obj.yafecAssembly.tangent_stiffness_and_force(zeros(nDOFs, 1));

    % Assemble stiffness tensors
    [obj.T2, obj.T3] = obj.yafecAssembly.mpc_tensors();

    % Number of constrained DOFs
    obj.n = size(obj.M, 1);

    % Damping matrix
    obj.C = sparse(obj.n, obj.n);

    % Number of design variables
    obj.ndv = length(mu);

    % Compute the partial derivatives
    if computePartialDerivatives
        % Initialize the sensitivity matrices
        obj.pM = repmat(struct('coeffs', []), obj.ndv, 1);
        obj.pK = repmat(struct('coeffs', []), obj.ndv, 1);
        obj.pC = repmat(struct('coeffs', sparse(obj.n, obj.n)), obj.ndv, 1);
        obj.pT2 = repmat(struct('coeffs', []), obj.ndv, 1);
        obj.pT3 = repmat(struct('coeffs', []), obj.ndv, 1);

        % Finite differences parameters
        mu0 = mu;
        dMu = 1e-4;

        % Use finite differences
        for dv = 1:obj.ndv
            % Perturb the parameter vector
            muPert = mu0;
            muPert(dv) = muPert(dv) + dMu;

            % Assemble the model
            yafecAssemblyPert = feval(yafec_mpc_constructor, muPert);

            % Mass matrix
            Mpert = yafecAssemblyPert.mass_matrix();
            obj.pM(dv).coeffs = (Mpert - obj.M) / dMu;

            % Stiffness matrix
            Kpert = yafecAssemblyPert.tangent_stiffness_and_force(zeros(nDOFs, 1));
            obj.pK(dv).coeffs = (Kpert - obj.K) / dMu;

            % Stiffness tensors 2 and 3
            [T2pert, T3pert] = yafecAssemblyPert.mpc_tensors();
            obj.pT2(dv).coeffs = (T2pert - obj.T2) / dMu;
            obj.pT3(dv).coeffs = (T3pert - obj.T3) / dMu;
        end
    end
end

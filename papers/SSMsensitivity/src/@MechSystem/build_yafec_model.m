function build_yafec_model(obj, mu, yafec_assembly_constructor, varargin)
    % Build the yafec model representing the mechanical system.
    % Inputs:
    %   mu: the parameter vector.
    %   yafec_assembly_constructor: function that produces the yafec
    %       assembly having as input a parameter vector mu.
    %   computePartialDerivatives: flag to enable the computation of the
    %       partial derivatives (optional, default is false).

    % Parse the inputs
    p = inputParser;
    addOptional(p, 'computePartialDerivatives', false);
    parse(p, varargin{:});
    computePartialDerivatives = p.Results.computePartialDerivatives;

    % Create the assembly
    obj.yafecAssembly = feval(yafec_assembly_constructor, mu);

    % Number of unconstrained and constrained DOFs
    nDOFs = obj.yafecAssembly.Mesh.nDOFs;
    obj.n = length(obj.yafecAssembly.Mesh.EBC.unconstrainedDOFs);

    % Assmeble mass matrix
    obj.M = obj.yafecAssembly.mass_matrix();
    obj.M = obj.yafecAssembly.constrain_matrix(obj.M);

    % Assemble stiffness matrix
    obj.K = obj.yafecAssembly.tangent_stiffness_and_force(zeros(nDOFs, 1));
    obj.K = obj.yafecAssembly.constrain_matrix(obj.K);

    % Assemble stiffness tensors
    obj.T2 = obj.yafecAssembly.tensor('T2', [nDOFs nDOFs nDOFs], [2, 3]);
    obj.T2 = obj.yafecAssembly.constrain_tensor(obj.T2);
    obj.T3 = obj.yafecAssembly.tensor('T3', [nDOFs nDOFs nDOFs nDOFs], [2, 3, 4]);
    obj.T3 = obj.yafecAssembly.constrain_tensor(obj.T3);

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
            yafecAssemblyPert = feval(yafec_assembly_constructor, muPert);

            % Mass matrix
            Mpert = yafecAssemblyPert.mass_matrix();
            Mpert = yafecAssemblyPert.constrain_matrix(Mpert);
            obj.pM(dv).coeffs = (Mpert - obj.M) / dMu;

            % Stiffness matrix
            Kpert = yafecAssemblyPert.tangent_stiffness_and_force(zeros(nDOFs, 1));
            Kpert = yafecAssemblyPert.constrain_matrix(Kpert);
            obj.pK(dv).coeffs = (Kpert - obj.K) / dMu;

            % Stiffness tensor 2
            T2pert = yafecAssemblyPert.tensor('T2', [nDOFs nDOFs nDOFs], [2, 3]);
            T2pert = yafecAssemblyPert.constrain_tensor(T2pert);
            obj.pT2(dv).coeffs = (T2pert - obj.T2) / dMu;

            % Stiffness tensor 3
            T3pert = yafecAssemblyPert.tensor('T3', [nDOFs nDOFs nDOFs nDOFs], [2, 3, 4]);
            T3pert = yafecAssemblyPert.constrain_tensor(T3pert);
            obj.pT3(dv).coeffs = (T3pert - obj.T3) / dMu;
        end
    end
end

function build_topopt_system(obj, yafecAssembly, varargin)
    % Build the yafec model representing the mechanical system for topology
    % optimization problems. In this case, the design variables are the
    % densities of the elements.
    % Inputs:
    %   yafecAssembly: the yafec Assembly object.
    %   mu: the design variables (optional, default is []).
    %   computePartialDerivatives: flag to enable the computation of the
    %       partial derivatives (optional, default is false).

    % Parse the inputs
    p = inputParser;
    addOptional(p, 'mu', []);
    addOptional(p, 'computePartialDerivatives', false);
    parse(p, varargin{:});
    computePartialDerivatives = p.Results.computePartialDerivatives;
    mu = p.Results.mu;

    % Number of unconstrained and constrained DOFs
    obj.yafecAssembly = yafecAssembly;
    obj.n = length(yafecAssembly.Mesh.EBC.unconstrainedDOFs);

    % Number of design variables (same as the number of elements)
    obj.ndv = yafecAssembly.Mesh.nElements;

    % Check if the design variables are provided
    if ~isempty(mu)
        % Update the yafec model
        obj.update_topopt_system(mu);
    end

    % Compute the partial derivatives of matrices and tensors with respect
    % to the design variables.
    % For topology optimization problems, the partial derivatives are
    % constant and can be computed once. Moreover, the partial derivatives
    % are equal for all the elements. This allows the computation of the
    % partial derivatives at the element level, thus avoiding the
    % full assembly of the matrices and tensors.
    if computePartialDerivatives
        % Compute the partial derivatives at element level
        obj.Me = yafecAssembly.Mesh.Elements(1).Object.mass_matrix();
        obj.Ke = yafecAssembly.Mesh.Elements(1).Object.tangent_stiffness_and_force(zeros(8,1));
        obj.Ce = zeros(size(obj.Ke));
        obj.T2e = sptensor(yafecAssembly.Mesh.Elements(1).Object.T2());
        obj.T3e = sptensor(yafecAssembly.Mesh.Elements(1).Object.T3());

        % % Initialize partial derivatives
        % warning('Add the computation of the partial derivatives at element level.');
        % obj.pM = repmat(struct('coeffs', []), 1, obj.ndv);
        % obj.pK = repmat(struct('coeffs', []), 1, obj.ndv);
        % obj.pC = repmat(struct('coeffs', sparse(obj.n, obj.n)), 1, obj.ndv);
        % obj.pT2 = repmat(struct('coeffs', []), 1, obj.ndv);
        % obj.pT3 = repmat(struct('coeffs', []), 1, obj.ndv);
        % 
        % % Sparse assembly
        % nDOFs = yafecAssembly.Mesh.nDOFs;
        % for dv = 1:obj.ndv
        %     % Extract the DOFs relative to the current element
        %     thisElement = yafecAssembly.Mesh.Elements(dv).Object;
        %     iDOFs = thisElement.iDOFs;
        %     I = kron(ones(size(iDOFs)), iDOFs);
        %     J = kron(iDOFs, ones(size(iDOFs)));
        % 
        %     % Mass matrix
        %     pM = sparse(I, J, obj.Me(:), nDOFs, nDOFs);
        %     obj.pM(dv).coeffs = yafecAssembly.constrain_matrix(pM);
        % 
        %     % Stiffness matrix
        %     pK = sparse(I, J, obj.Ke(:), nDOFs, nDOFs);
        %     obj.pK(dv).coeffs = yafecAssembly.constrain_matrix(pK);
        % 
        %     % Quadratic stiffness tensor
        %     SIZE = [nDOFs, nDOFs, nDOFs];
        %     pT2 = sptensor(iDOFs(obj.T2e.subs), obj.T2e.vals, SIZE);
        %     obj.pT2(dv).coeffs = yafecAssembly.constrain_tensor(pT2);
        % 
        %     % Cubic stiffness tensor
        %     SIZE = [nDOFs, nDOFs, nDOFs, nDOFs];
        %     pT3 = sptensor(iDOFs(obj.T3e.subs), obj.T3e.vals, SIZE);
        %     obj.pT3(dv).coeffs = yafecAssembly.constrain_tensor(pT3);
        % end
    end
end

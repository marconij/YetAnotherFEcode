function update_topopt_system(obj, mu)
    % Update the yafec model representing the mechanical system for
    % topology optimization problems. In this case, the design variables
    % are the densities of the elements.
    % Inputs:
    %   mu: the design variables.

    % Number of unconstrained and constrained DOFs
    nDOFs = obj.yafecAssembly.Mesh.nDOFs;

    % Assmeble mass matrix
    obj.M = obj.yafecAssembly.mass_matrix_uniform('weights', mu);
    obj.M = obj.yafecAssembly.constrain_matrix(obj.M);

    % Assemble stiffness matrix
    obj.K = obj.yafecAssembly.tangent_stiffness_and_force_uniform(zeros(nDOFs, 1), 'weights', mu);
    obj.K = obj.yafecAssembly.constrain_matrix(obj.K);

    % Assemble stiffness tensors
    obj.T2 = obj.yafecAssembly.tensor_uniform('T2', [nDOFs nDOFs nDOFs], [2, 3], 'weights', mu);
    obj.T2 = obj.yafecAssembly.constrain_tensor(obj.T2);
    obj.T3 = obj.yafecAssembly.tensor_uniform('T3', [nDOFs nDOFs nDOFs nDOFs], [2, 3, 4], 'weights', mu);
    obj.T3 = obj.yafecAssembly.constrain_tensor(obj.T3);

    % Damping matrix
    obj.C = sparse(obj.n, obj.n);
end

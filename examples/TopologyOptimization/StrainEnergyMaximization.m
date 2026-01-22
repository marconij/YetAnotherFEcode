clear; clc; close all;

%% Problem settings
% Material
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio 
thickness = .1;     % [m] beam's out-of-plane thickness
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
myElementConstructor = @()Quad4Element(thickness, myMaterial);

% Mesh
lx = 1; ly = 1;
nelx = 100; nely = 100;
[nodes, elements, nset] = mesh_2Drectangle(lx, ly, nelx, nely, 'QUAD4');
myMesh = Mesh(nodes);
myMesh.create_elements_table(elements, myElementConstructor);

% Fixed boundary conditions
fixedNodes = find(all(abs(nodes - [0, lx]) < [0.41*lx, 1e-12], 2));
myMesh.set_essential_boundary_condition(fixedNodes, 1:2, 0);

% Assembly
myAssembly = Assembly(myMesh);

% External action
F0 = -1e3*0;
uC0 = -0.1*ly;
F = zeros(myMesh.nDOFs, 1);
boundaryNodes = find(all(abs(nodes - [lx, 0.35*ly]) < [1e-12, 0.051*ly], 2));
if uC0 ~= 0
    myAssembly.Mesh.set_essential_boundary_condition(boundaryNodes, 2, uC0);
elseif F0 ~= 0
    force_dofs = myMesh.get_DOF_from_nodeIDs(boundaryNodes);
    F(force_dofs(:, 2), 1) = F0;
end
Fu = myAssembly.constrain_vector(F);
uC = myMesh.EBC.constrainedDOFs(:, 2);

% Elements centroid
coord = zeros(myMesh.nElements, 2);
for ii = 1:myMesh.nElements
    coord(ii, :) = mean(myMesh.Elements(ii).Object.nodes);
end

% Element stiffness matrix
Ke = myMesh.Elements(1).Object.stiffness_matrix();

% Area
Ae = myMesh.Elements(1).Object.area;
Atot = Ae * myMesh.nElements;

%% Initialize topology optimization
radius = 2;
beta = 10; eta = 0.5;
dMinSimp = 1e-6; p = 3;

% Initialize object
to = TopologyOptimization([nelx, nely], coord, radius, beta, eta, dMinSimp, p);

% Initial layout
to.initialize_density(0.5);
to.set_density_box([lx, ly], 0.6 * [lx, ly], 0);
to.set_density_box([lx, 0.35*ly], [0.1 * lx, 0.05 * ly], 1);

% Initial layout
figure();
plot_layout(to.nel, to.d, to.mapFea2To);
title('Initial Layout', 'Interpreter', 'latex');
drawnow;

%% Initialize optimizer
m = 1;
move = 0.01;
mma = MMA(m, move, to);

% Iterations
maxIter = 200;

% History file
history = NaN(m + 1, maxIter);
densHistory = NaN(to.nElements, maxIter);

% Initialize figure
figure(); drawnow;

%% Main loop

% Header
fprintf("\nIteration - Objective - Constraints\n");

% Start timer
tStart = tic;

% Loop
for iter = 1:maxIter
    % Apply filtering and projection stages
    to.filter();
    to.projection();
    to.simp();

    % Current area
    A = Ae * sum(to.d_proj);

    % Assemble matrix
    K = myAssembly.matrix_uniform('stiffness_matrix', 'weights', to.d_simp);
    Kuu = myAssembly.constrain_matrix(K);
    Kuc = K(myAssembly.Mesh.EBC.unconstrainedDOFs, myAssembly.Mesh.EBC.constrainedDOFs(:, 1));
    Kcu = K(myAssembly.Mesh.EBC.constrainedDOFs(:, 1), myAssembly.Mesh.EBC.unconstrainedDOFs);
    Kcc = K(myAssembly.Mesh.EBC.constrainedDOFs(:, 1), myAssembly.Mesh.EBC.constrainedDOFs(:, 1));

    % Solve stationary problem
    uU = Kuu \ (Fu - Kuc * uC);
    u = myAssembly.unconstrain_vector(uU);
    J = 0.5 * dot(u, K*u) - dot(Fu, uU);

    % Store initial value
    if iter == 1
        J0 = abs(J);
    end

    % Physical density sensitivity
    sensPh = to.simp_sensitivity();

    % Compute physical sensitivity
    dJdd = SensitivityLibrary.strain_energy_nonhomogeneous(myAssembly, u, Ke, 'sensPh', sensPh);
    dAdd = Ae * ones(myMesh.nElements, 1);

    % Compute filter sensitivity
    dJ = to.filter_sensitivity(dJdd);
    dA = to.filter_sensitivity(dAdd);
    
    % Print current iteration
    fprintf("\n%4d %16.4e %16.4e", iter, J, A / Atot);

    % Plot current layout
    plot_layout(to.nel, to.d_proj, to.mapFea2To); drawnow;

    % Design variables (n x 1)
    xval  = to.d;

    % Objective function and sensitivity (n x 1)
    f0val = -J / J0;
    df0dx = -dJ(:) / J0;

    % Constraints (m x 1) and sensitivity (m x n)
    fval  = [A / Atot - 0.4];
    dfdx  = [dA(:).' / Atot];

    % Save current iteration
    history(:, iter) = [f0val; fval];
    densHistory(:, iter) = to.d_proj;
   
    % Convergence criterion
    if iter > 5 % check convergence after 5 iterations
        fval_tol = 1e-3;
        if all(fval < fval_tol) % all the constraints are satisfied
            err = abs(1 - history(1, iter-3:iter-1) / history(1, iter));
            err_tol = 1e-3;
            if all(err < err_tol) % check relative variation
                break
            end
        end
    end

    % MMA step
    to.d = mma.optimize(iter, xval, f0val, df0dx, fval, dfdx);
end

% Stop timer and display elapsed time
tElapsed = toc(tStart);
fprintf('\n\nEnd of the optimization.\n');
fprintf('Elapsed time is %f seconds.\n', tElapsed)

%% Optimal results

% History
figure();
plot_history(history);

% Optimal layout
figure();
plot_layout(to.nel, to.d_proj, to.mapFea2To);
title('Optimal Layout', 'Interpreter', 'latex');

% Create gif of the density evolution
% create_gif(to.nel, densHistory, 'mapFea2To', to.mapFea2To, 'fileName', 'ElectrostaticOptimization');

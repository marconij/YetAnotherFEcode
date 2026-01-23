clear; clc; close all;

%% Problem settings
% Material
E = 1; % Young's modulus [Pa]
rho = 1; % density [kg/m^3]
nu = 0.3; % Poisson's ratio 
thickness = 1; % [m] beam's out-of-plane thickness
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
myElementConstructor = @()Quad4Element(thickness, myMaterial);

% Mesh
lx = 2; ly = 0.5;
nelx = 200; nely = 50;
[nodes, elements, nset] = mesh_2Drectangle(lx, ly, nelx, nely, 'QUAD4');
myMesh = Mesh(nodes);
myMesh.create_elements_table(elements, myElementConstructor);

% Boundary conditions
myMesh.set_essential_boundary_condition(nset{1}, 1:2, 0);

% Assembly
myAssembly = Assembly(myMesh);

% Nodal force
F = zeros(myMesh.nDOFs, 1);
loadNode = find_node(lx, ly/2, [], nodes);
loadDofs = get_index(loadNode, myMesh.nDOFPerNode);
F0 = 1;
F(loadDofs(1)) = -F0;
Fc = myAssembly.constrain_vector(F);

% Elements centroid
coord = zeros(myMesh.nElements, 2);
for ii = 1:myMesh.nElements
    coord(ii, :) = mean(myMesh.Elements(ii).Object.nodes);
end

% Element stiffness matrix
Ke = myMesh.Elements(1).Object.tangent_stiffness_and_force(zeros(8,1));

% Area
Ae = myMesh.Elements(1).Object.area;
Atot = Ae * myMesh.nElements;

% Null displacement vector
% u0 = zeros(myMesh.nDOFs, 1);

%% Initialize topology optimization
radius = 4;
beta = 10; eta = 0.5;
dMinSimp = 1e-6; p = 3;

% Initialize object
to = TopologyOptimization([nelx, nely], coord, radius, beta, eta, dMinSimp, p);

% Initial layout
to.initialize_density(0.5);

% Initial layout
figure();
plot_layout(to.nel, to.d, to.mapFea2To);
title('Initial Layout', 'Interpreter', 'latex');
drawnow

%% Initialize optimizer
m = 2; % 1 for compliance minimization, 2 for buckling maximization
move = 0.02;
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

    % Assemble matrices
    K = myAssembly.stiffness_matrix_uniform('weights', to.d_simp);
    Kc = myAssembly.constrain_matrix(K);

    % Solve stationary problem
    uc = Kc \ Fc;
    u = myAssembly.unconstrain_vector(uc);

    % Evaluate compliance
    C = dot(Fc, uc);
    if iter == 1
        C0 = C;
    end

    % Evaluate buckling
    if m > 1
        % Geometric stiffness matrix
        G = myAssembly.geometric_stiffness_matrix(u, 'weights', to.d_simp);
        Gc = myAssembly.constrain_matrix(G);

        % Solve buckling problem
        nModes = 10;
        [V, D] = eigs(Gc, Kc, nModes, 'largestabs');
        P = -1 ./ diag(D);

        % Extract mode
        modeIndex = 1;
        Pi = P(modeIndex);
        vc = V(:, modeIndex);
        vc = vc / sqrt(vc.' * (-Gc) * vc);
        v = myAssembly.unconstrain_vector(vc);

        % Store initial value
        if iter == 1
            P0 = Pi;
        end
    else
        Pi = 0;
    end

    % Physical density sensitivity
    sensPh = to.simp_sensitivity();

    % Compute physical sensitivity
    dCdd = SensitivityLibrary.compliance(myMesh, u, Ke, 'sensPh', sensPh);
    dAdd = Ae * ones(myMesh.nElements, 1);

    % Compute filter sensitivity
    dC = to.filter_sensitivity(dCdd);
    dA = to.filter_sensitivity(dAdd);

    % Compute buckling sensitivity
    if m > 1
        if iter == 1
            pG_u = SensitivityLibrary.initialize_pG_u(myAssembly);
        end
        dPdd = SensitivityLibrary.buckling(myAssembly, Kc, Pi, u, v, Ke, pG_u, to.d_simp, ...
            'sensPhK', sensPh, 'sensPhG', sensPh);
        dP = to.filter_sensitivity(dPdd);
    end
    
    % Print current iteration
    fprintf("\n%4d %16.4e %16.4e %16.4e", iter, C, Pi, A / Atot);

    % Plot current layout
    plot_layout(to.nel, to.d_proj, to.mapFea2To); drawnow;

    % Design variables (n x 1)
    xval  = to.d;

    % Optimization problem
    if m == 1 % compliance minimization
        % Objective function and sensitivity (n x 1)
        f0val = C / C0;
        df0dx = dC(:) / C0;

        % Constraints (m x 1) and sensitivity (m x n)
        fval  = [A / Atot - 0.6];
        dfdx  = [dA(:).' / Atot];
    else % buckling maximization with compliance constraint
        % Objective function and sensitivity (n x 1)
        f0val = -Pi / P0;
        df0dx = -dP(:)/ P0;

        % Constraints (m x 1) and sensitivity (m x n)
        Clim = 10;
        fval = [A / Atot - 0.6; (C - Clim) / C0];
        dfdx = [dA(:).' / Atot; dC(:).' / C0];
    end

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

%% Post processing

% Solve stationary problem
uc = Kc \ Fc;
u = myAssembly.unconstrain_vector(uc);

% Evaluate compliance
C = dot(Fc, uc);

% Geometric stiffness matrix
G = myAssembly.geometric_stiffness_matrix(u, 'weights', to.d_simp);
Gc = myAssembly.constrain_matrix(G);

% Solve buckling problem
nModes = 10;
[V, D] = eigs(Gc, Kc, nModes, "largestabs");
P = -1 ./ diag(D);

% Extract mode
modeIdx = 1;
Pi = P(modeIdx);
vc = V(:, modeIdx);
vc = vc ./ sqrt(vc.' * (-Gc) * vc); % mass normalization
v = myAssembly.unconstrain_vector(vc);

% Print
fprintf("\nCompliance: %.4e\nBuckling load: %.4e\n\n", C, Pi)

% Buckling mode
figure
plot_layout_deformed(to.nel, to.d_proj, nodes, v, to.mapFea2To)

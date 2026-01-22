clear; clc; close all;

%% Problem settings
% Material
E = 148e9; % Young's modulus [Pa]
rho = 2330e-6; % density [ng/um^3]
nu = 0.23; % Poisson's ratio
thickness = 24; % [um] beam's out-of-plane thickness
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
myElementConstructor = @()Quad4Element(thickness, myMaterial);

% Mesh
lx = 500; ly = 100; % [um]
nelx = 100; nely = 20;
[nodes, elements, nset] = mesh_2Drectangle(lx, ly, nelx, nely, 'QUAD4');
myMesh = Mesh(nodes);
myMesh.create_elements_table(elements, myElementConstructor);

% Boundary conditions
myMesh.set_essential_boundary_condition(nset{1}, 1:2, 0);
myMesh.set_essential_boundary_condition(nset{3}, 1, 0);

% Assembly
myAssembly = Assembly(myMesh);
nDOFs = myMesh.nDOFs;

% Elements centroid
coord = zeros(myMesh.nElements, 2);
for ii = 1:myMesh.nElements
    coord(ii, :) = mean(myMesh.Elements(ii).Object.nodes);
end

% Element quantities
Ke = myMesh.Elements(1).Object.tangent_stiffness_and_force(zeros(8,1));
Me = myMesh.Elements(1).Object.mass_matrix();

% Area
Ae = myMesh.Elements(1).Object.area;
Atot = Ae * myMesh.nElements;

% Null displacement vector
u0 = zeros(myMesh.nDOFs, 1);

%% SSM settings

% Damping
rayleigh = [0.01, 0];

% Expansion order
maxOrder = 5;

% Target point
outDofs = myMesh.get_DOF_from_location([lx, ly/2]);
outDof = myAssembly.free2constrained_index(outDofs(2));
xT = 0.2 * ly;

% Create system
sys = MechSystem();
sys.build_topopt_system(myAssembly, ...
    'computePartialDerivatives', true)

%% Initialize topology optimization
radius = 4;
beta = 10; eta = 0.5;
dMinSimp = 1e-6; p = 1;
dp = 0.25; dpIter = 5; pMax = 8;

% Initialize object
to = TopologyOptimization([nelx, nely], coord, radius, beta, eta, dMinSimp, p);

% Initial layout
to.initialize_density(0.5);

% Modify regions of the domain
to.set_density_box([lx, ly], [0.2 * lx, 1e10], 1);

% Initial layout
figure();
plot_layout(to.nel, to.d, to.mapFea2To);
title('Initial Layout', 'Interpreter', 'latex');
drawnow

%% Initialize optimizer
m = 2;
move = 0.01;
mma = MMA(m, move, to);

% Iterations
maxIter = 200;

% History file
history = NaN(m + 1, maxIter);
densHistory = NaN(to.nElements, maxIter);
timeHistory = NaN(1, maxIter);
convHistory = zeros(1, maxIter);
pHistory = NaN(1, maxIter);

% Initialize figure
figure(); drawnow;

%% Main loop

% Header
fprintf("\nIteration - Objective - Constraints\n");

% Start timer
tStart = tic;

% Loop
for iter = 1:maxIter
    tIter = tic;
    % Update p
    if mod(iter, dpIter) == 0 && to.p < pMax
        to.p = to.p + dp;
        fprintf("\nUpdate p: %.4e", to.p);
    end
    pHistory(iter) = to.p;

    % Apply filtering and projection stages
    to.filter();
    to.projection();
    to.simp();

    % Current area
    A = Ae * sum(to.d_proj);

    % Update the topology optimization system
    sys.update_topopt_system(to.d_simp);

    % Solve eigenvalue problem
    % At the first iteration, we specify the index of the taregt mode.
    % Then, we store the mode shape and use it in the next iterations as
    % the reference one for the MAC.
    if iter == 1
        sys.modal_analysis('nModes', 10, 'modeIndex', 1, ...
            'computePartialDerivatives', true);

        % Store initial value and reference mode shape
        modeRef = sys.phi;
    else
        sys.modal_analysis('nModes', 10, 'modeRef', modeRef, ...
            'computePartialDerivatives', true);
    end

    % Damping
    sys.set_damping_rayleigh(rayleigh, ...
        'computePartialDerivatives', true, 'topOptProblem', true);
    omega = sys.omegaDamped;

    % Compute SSM
    ssm = MechSSM(sys);
    ssm.compute_ssm(maxOrder, 'storeCoefficients', true);
    
    % Target point
    rhoT = ssm.compute_rho_from_x(xT, outDof);
    OmegaT = ssm.compute_omega(rhoT);
    [xT, xkT] = ssm.compute_x(rhoT, outDof);

    % Store initial values
    if iter == 1
        omega0 = omega;
        OmegaT0 = OmegaT;
    end

    % Physical density sensitivity
    sensPh = to.simp_sensitivity();

    % Compute adjoint sensitivity
    [domegaNat, domegaDamped] = sys.adjoint_sensitivity_omega_topopt;
    domegadd = sensPh .* domegaDamped;
    % domegadd = sensPh .* domegaNat;

    % Compute adjoint sensitivity
    sens = SensitivityBackbone(ssm);
    dOmegaT = sens.adjoint_sensitivity_Omega_topopt(rhoT, xkT, xT, outDof);
    dOmegaTdd = sensPh .* dOmegaT;

    % Compute area sensitivity
    dAdd = Ae * ones(myMesh.nElements, 1);

    % Compute filter sensitivity
    domega = to.filter_sensitivity(domegadd);
    dOmegaT = to.filter_sensitivity(dOmegaTdd);
    dA = to.filter_sensitivity(dAdd);

    % Print current iteration
    fprintf("\n%4d %16.4e %16.4e %16.4e", iter, omega/2/pi, OmegaT/2/pi, A / Atot);

    % Plot current layout
    plot_layout(to.nel, to.d_proj, to.mapFea2To); drawnow;

    % Design variables (n x 1)
    xval  = to.d;

    % Objective function and sensitivity (n x 1)
    f0val = -omega / omega0;
    df0dx = -domega(:) / omega0;

    % Constraints (m x 1) and sensitivity (m x n)
    fval  = [A / Atot - 0.5; (OmegaT - 550*2*pi) / OmegaT0];
    dfdx  = [dA(:).' / Atot; dOmegaT(:).' / OmegaT0];

    % Save current iteration
    history(:, iter) = [f0val; fval];
    densHistory(:, iter) = to.d_proj;

    % MMA step
    to.d = mma.optimize(iter, xval, f0val, df0dx, fval, dfdx);

    % Record time
    timeHistory(iter) = toc(tIter);
    fprintf('\tIteration time is %f seconds', timeHistory(iter))
    
    % Convergence criterion
    if iter > 5 % check convergence after 5 iterations
        fval_tol = 1e-3;
        if all(fval < fval_tol) % all the constraints are satisfied
            nConvIter = 2;
            err = abs(1 - history(1, iter - (nConvIter:-1:1)) / history(1, iter));
            err_tol = 1e-3;
            if all(err < err_tol) % check relative variation
                fprintf("\nConvergence has been reached at iteration %d", iter)
                convHistory(iter) = 1;
            end
        end
    end
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

% Computation time
figure();
hold on; grid on; box on; axis tight;
plot(timeHistory, 'LineWidth', 2);
xlabel('Iterations', 'Interpreter', 'latex');
ylabel('Time [s]', 'Interpreter', 'latex');
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');

% Convergence
figure();
hold on; grid on; box on; axis tight;
plot(convHistory, 'LineWidth', 2);
xlabel('Iterations', 'Interpreter', 'latex');
ylabel('Conv flag', 'Interpreter', 'latex');
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');
return
%% Post-processing
convIter = 98;
d_proj = densHistory(:, convIter);
d_simp = dMinSimp + (1 - dMinSimp) * d_proj .^ pHistory(convIter);

% Layout
figure();
plot_layout(to.nel, d_proj, to.mapFea2To);

% Update the topology optimization system
sys.update_topopt_system(d_simp);

% Solve eigenvalue problem
sys.modal_analysis('nModes', 10, 'modeRef', modeRef);

% Damping
sys.set_damping_rayleigh(rayleigh, ...
    'computePartialDerivatives', true, 'topOptProblem', true);
omega = sys.omegaDamped;

% Compute SSM
ssm = MechSSM(sys);
ssm.compute_ssm(maxOrder);

% Target point
rhoT = ssm.compute_rho_from_x(xT, outDof);
OmegaT = ssm.compute_omega(rhoT);

% Physical amplitude
rho = linspace(0, 1.1*rhoT, 101);
Omega = ssm.compute_omega(rho);
x = ssm.compute_x(rho, outDof);

% Print
disp('Iterations:')
disp(convIter)
s = seconds(sum(timeHistory(1:convIter))); s.Format = 'hh:mm:ss';
disp('Total time [hh:mm:ss]:')
disp(s)
disp('Time per iteration [s]:')
disp(mean(timeHistory(1:convIter)))
disp('Eigenfrequency [kHz]:')
disp(sys.omegaDamped/2/pi)
disp('Frequency [kHz]:')
disp(OmegaT/2/pi)

% Plot
figure
hold on; grid on; box on; axis square tight;
plot(Omega / (2*pi), x, 'k', 'LineWidth', 2)
plot(OmegaT / (2*pi), xT, '.r', 'MarkerSize', 20)
xlabel('$\Omega$ [kHz]', 'Interpreter', 'latex')
ylabel('$x$ [$\mu$m]', 'Interpreter', 'latex')
legend('Backbone', 'Target Point', 'Interpreter', 'latex', 'Location', 'southwest')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
xlim([610, 700])

% Computation time
timeHist = timeHistory(1:convIter);
meanTime = mean(timeHist);
meanTime2 = mean(timeHist(timeHist < 150));
figure('Position', [488   338   1000   420]);
hold on; grid on; box on; axis tight;
plot(timeHist, 'LineWidth', 2);
plot([0, convIter], meanTime * [1, 1], '--', 'LineWidth', 2);
plot([0, convIter], meanTime2 * [1, 1], '--', 'LineWidth', 2);
xlabel('Iterations', 'Interpreter', 'latex');
ylabel('Time per iteration [s]', 'Interpreter', 'latex');
legend('', 'Mean time', 'Mean time w/o outliers', 'Location', 'northeast', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex');
xlim([0, convIter])

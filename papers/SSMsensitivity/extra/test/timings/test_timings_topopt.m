% Time the computation of the SSM and its sensitivity
clear; clc; close all;
addpath(genpath('../../src'));

%% Problem settings
% Material
E         = 70e9;   % Young's modulus [Pa]
rho       = 2700;   % density [kg/m^3]
nu        = 0.33;   % Poisson's ratio
thickness = 0.1;    % [m] beam's out-of-plane thickness
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
myElementConstructor = @()Quad4Element(thickness, myMaterial);

% Dimensions
lx = 2.0; ly = 1.0;

%% Timings

% Number of elements
nxVec = 20:10:100;
nyVec = nxVec / 2;
nVec = nxVec .* nyVec;

% Expansion order
maxOrderVec = [3, 5, 7];

% Initialize
timeSSM = zeros(length(nVec), length(maxOrderVec));
timeSens = zeros(length(nVec), length(maxOrderVec));

% Loop over the number of dofs
for i = 1:length(nVec)
    % Number of elements
    nx = nxVec(i);
    ny = nyVec(i);
    n = nVec(i);

    % Create density
    d = ones(n, 1);

    % Create mesh
    [nodes, elements, nset] = mesh_2Drectangle(lx, ly, nx, ny, 'QUAD4');
    myMesh = Mesh(nodes);
    myMesh.create_elements_table(elements, myElementConstructor);
    
    % Boundary conditions
    myMesh.set_essential_boundary_condition(nset{1}, 1:2, 0);
    myMesh.set_essential_boundary_condition(nset{3}, 1, 0);
    
    % Assembly
    myAssembly = Assembly(myMesh);

    % SSM settings
    rayleigh = struct('alphaR', 0.01, 'betaR', 0);
    outDofs = myMesh.get_DOF_from_location([lx, ly/2]);
    outDof = myAssembly.free2constrained_index(outDofs(2));
    xTarget = 0.2 * ly;

    % Create system
    sys = MechSystem();
    sys.build_topopt_system(myAssembly, rayleigh, ...
        'computePartialDerivatives', true);
    sys.update_topopt_system(d);

    % Modal analysis
    sys.modal_analysis('nModes', 10, 'modeIndex', 1, ...
        'computePartialDerivatives', true);
    
    % Loop over the expansion order
    for j = 1:length(maxOrderVec)
        % Expansion order
        maxOrder = maxOrderVec(j);

        % Compute SSM
        tStart = tic;
        ssm = MechSSM(sys);
        ssm.compute_ssm(maxOrder, 'storeCoefficients', true);
        rhoTarget = ssm.compute_rho_from_x(xTarget, outDof);
        [xTarget, xkTarget] = ssm.compute_x(rhoTarget, outDof);
        Omega = ssm.compute_omega(rhoTarget);
        timeSSM(i, j) = toc(tStart);
        fprintf('n = %5d, maxOrder = %3d, SSM = %5.2e\n', n, maxOrder, timeSSM(i, j));

        % Compute sensitivity
        tStart = tic;
        sens = SensitivityBackbone(ssm);
        dOmega = sens.adjoint_sensitivity_Omega_topopt(rhoTarget, xkTarget, xTarget, outDof);
        timeSens(i, j) = toc(tStart);
        fprintf('n = %5d, maxOrder = %3d, Sens = %5.2e\n', n, maxOrder, timeSens(i, j));
    end
end

%% Plot
fig = figure;
tiledlayout('flow', 'TileSpacing', 'tight')

% Create grid
[xx, yy] = meshgrid(nVec, maxOrderVec);
xx = xx.';
yy = yy.';

% SSM
nexttile
hold on; grid on; box on; axis tight; view(3);
mesh(xx, yy, timeSSM)
xlabel('Number of elements [-]', 'Interpreter', 'latex')
ylabel('Expansion order [-]', 'Interpreter', 'latex')
zlabel('Time [s]', 'Interpreter', 'latex')
title('SSM', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Sensitivity adjoint
nexttile
hold on; grid on; box on; axis tight; view(3);
mesh(xx, yy, timeSens)
xlabel('Number of elements [-]', 'Interpreter', 'latex')
ylabel('Expansion order [-]', 'Interpreter', 'latex')
zlabel('Time [s]', 'Interpreter', 'latex')
title('Sensitivity Adjoint', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Save
saveas(fig, 'sens_topopt.fig')

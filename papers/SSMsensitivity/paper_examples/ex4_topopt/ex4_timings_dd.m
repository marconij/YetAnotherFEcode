% Time the computation of the SSM and its sensitivity
% Uncomment from line 74 in build_topopt_system
clear; clc; close all;
addpath(genpath('../../src'));

%% Problem settings
% Material
E         = 148e9;  % Young's modulus [Pa]
rho       = 2330;   % density [kg/m^3]
nu        = 0.23;   % Poisson's ratio
thickness = 24;     % [m] beam's out-of-plane thickness
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myMaterial.PLANE_STRESS = true;	% set "false" for plane_strain
myElementConstructor = @()Quad4Element(thickness, myMaterial);

% Dimensions
lx = 500; ly = 100;

% Number of elements
nx = 50;
ny = nx/5;
n = nx*ny;

% Expansion order
maxOrder = 3;

%% Timings

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
rayleigh = [0.01, 0];
outDofs = myMesh.get_DOF_from_location([lx, ly/2]);
outDof = myAssembly.free2constrained_index(outDofs(2));
xTarget = 0.2 * ly;

% Start timer
tStart = tic;

% Create system
sys = MechSystem();
sys.build_topopt_system(myAssembly, ...
    'computePartialDerivatives', true);
sys.update_topopt_system(d);

% Modal analysis
sys.modal_analysis('nModes', 10, 'modeIndex', 1, ...
    'computePartialDerivatives', true);
sys.set_damping_rayleigh(rayleigh, 'computePartialDerivatives', true, ...
    'topOptProblem', false)
sys.modal_analysis_sensitivity;

% Compute SSM
ssm = MechSSM(sys);
ssm.compute_ssm(maxOrder, 'storeCoefficients', true, 'computeSensitivity', true);
rhoTarget = ssm.compute_rho_from_x(xTarget, outDof);
[xTarget, xkTarget] = ssm.compute_x(rhoTarget, outDof);
Omega = ssm.compute_omega(rhoTarget);

% Compute sensitivity
% sens = SensitivityBackbone(ssm);
% dOmega = sens.adjoint_sensitivity_Omega_topopt(rhoTarget, xkTarget, xTarget, outDof);

% Time
timeDd = toc(tStart);
fprintf('n = %5d, maxOrder = %3d, SensDd = %5.2e\n', n, maxOrder, timeDd);

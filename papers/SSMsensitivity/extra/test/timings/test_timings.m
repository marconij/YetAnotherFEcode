% Time the computation of the SSM and its sensitivity
clear; clc; close all;
addpath(genpath('../../src'));

%% Constant parameters

% Parameters
m = 1;
k = 1;
k2 = 0.3;
k3 = 0.5;

% Rayleigh damping
rayleigh = [0.01, 0.1];

% Target physical amplitude
xTarget = 3;

% Output dof
outDof = 5;

%% Timings

% Number of dofs
nVec = [5, 50, 500];

% Expansion order
maxOrderVec = [5, 7, 9, 11];

% Initialize
timeSSM = zeros(length(nVec), length(maxOrderVec));
timeSensDD = zeros(length(nVec), length(maxOrderVec));
timeSensAdj = zeros(length(nVec), length(maxOrderVec));

% Loop over the number of dofs
for i = 1:length(nVec)
    % Number of dofs
    n = nVec(i);
    
    % Loop over the expansion order
    for j = 1:length(maxOrderVec)
        % Expansion order
        maxOrder = maxOrderVec(j);

        % Create system
        sys = MechSystem();
        sys.build_uniform_spring_mass_chain(n, m, k, k2, k3, ...
            'computePartialDerivatives', true);

        % Modal analysis
        sys.modal_analysis('computePartialDerivatives', true);
        sys.set_damping_rayleigh(rayleigh, 'computePartialDerivatives', true);

        % Modal analysis sensitivity
        sys.modal_analysis_sensitivity;

        % Initialize the SSM
        ssm = MechSSM(sys);

        % Compute the SSM w/o sensitivity
        tStart = tic;
        ssm.compute_ssm(maxOrder, 'computeSensitivity', false, ...
            'storeCoefficients', true);
        rhoTarget = ssm.compute_rho_from_x(xTarget, outDof);
        xTarget = ssm.compute_x(rhoTarget, outDof);
        ssm.compute_omega(rhoTarget);
        timeSSM(i, j) = toc(tStart);
        fprintf('n = %5d, maxOrder = %3d, SSM = %5.2e\n', n, maxOrder, timeSSM(i, j));

        % Compute the SSM w/ sensitivity
        tStart = tic;
        ssm.compute_ssm(maxOrder, 'computeSensitivity', true, ...
            'storeCoefficients', true);
        rhoTarget = ssm.compute_rho_from_x(xTarget, outDof);
        [xTarget, xkTarget] = ssm.compute_x(rhoTarget, outDof);
        OmegaTarget = ssm.compute_omega(rhoTarget);
        drho = ssm.compute_rho_derivative(rhoTarget, outDof, xkTarget);
        dOmegaDD = ssm.compute_omega_derivative(rhoTarget, drho);
        timeSensDD(i, j) = toc(tStart) - timeSSM(i, j);
        fprintf('n = %5d, maxOrder = %3d, SensDD = %5.2e\n', n, maxOrder, timeSensDD(i, j));

        % Compute adjoint sensitivity
        tStart = tic;
        sens = SensitivityBackbone(ssm);
        dOmegaAdj = sens.adjoint_sensitivity_Omega(rhoTarget, xkTarget, xTarget, outDof);
        timeSensAdj(i, j) = toc(tStart);
        fprintf('n = %5d, maxOrder = %3d, SensAdj = %5.2e\n', n, maxOrder, timeSensAdj(i, j));

        % Check
        sensDD = dOmegaDD;
        sensAdj = dOmegaAdj;
        absErr = norm(sensDD - sensAdj);
        fprintf('dOmegaDD = %5.2e, dOmegaAdj = %5.2e, absErr = %5.2e\n', ...
            norm(sensDD), norm(sensAdj), absErr);
    end
end

%% Plot
figure
tiledlayout('flow', 'TileSpacing', 'compact')

% Create grid
[xx, yy] = meshgrid(nVec, maxOrderVec);
xx = xx.';
yy = yy.';

% SSM
% nexttile
% hold on; grid on; box on; axis tight; view(3);
% mesh(xx, yy, timeSSM)
% xlabel('Number of dofs [-]', 'Interpreter', 'latex')
% ylabel('Expansion order [-]', 'Interpreter', 'latex')
% zlabel('Time [s]', 'Interpreter', 'latex')
% title('SSM', 'Interpreter', 'latex')
% set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Sensitivity direct differentiation
nexttile
hold on; grid on; box on; axis tight; view(3);
mesh(xx, yy, timeSensDD)
xlabel('Number of dofs [-]', 'Interpreter', 'latex')
ylabel('Expansion order [-]', 'Interpreter', 'latex')
zlabel('Time [s]', 'Interpreter', 'latex')
title('Sensitivity Direct Differentiation', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Sensitivity adjoint
nexttile
hold on; grid on; box on; axis tight; view(3);
mesh(xx, yy, timeSensAdj)
xlabel('Number of dofs [-]', 'Interpreter', 'latex')
ylabel('Expansion order [-]', 'Interpreter', 'latex')
zlabel('Time [s]', 'Interpreter', 'latex')
title('Sensitivity Adjoint', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Save
saveas(gcf, 'sensitivity.fig')

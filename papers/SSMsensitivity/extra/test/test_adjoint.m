clear; clc; close all;
addpath(genpath('../../src'));

%% System
% Parameters
n = 10;
m = 1;
k = 1;
k2 = 0.3;
k3 = 0.5;

% Rayleigh damping
alphaR = 0.01;
betaR = 0.1;
Qfactor = 100;

% Create system
sys = MechSystem();
sys.build_uniform_spring_mass_chain(n, m, k, k2, k3, 'computePartialDerivatives', true);

% Index of the master mode
modeIdx = 1;

% Modal analysis
sys.modal_analysis('modeIndex', modeIdx, 'computePartialDerivatives', true);
sys.set_damping_rayleigh([alphaR, betaR], 'computePartialDerivatives', true);
% sys.set_damping_qfactor(Qfactor, 'computePartialDerivatives', true);

% Modal analysis sensitivity
sys.modal_analysis_sensitivity;

% Adjoint sensitivity of of omega
[domega, domegaDamped] = sys.adjoint_sensitivity_omega;

%% Spectral Submanifold
% Maximum expansion order
maxOrder = 7;

% Initialize the SSM
ssm = MechSSM(sys);

% Compute the SSM
ssm.compute_ssm(maxOrder, 'storeCoefficients', true, 'computeSensitivity', true);

% Target physical amplitude
xTarget = 4;

% Output dof
outDof = 5;

% Compute the corresponding reduced amplitude
rhoTarget = ssm.compute_rho_from_x(xTarget, outDof);

% Reduced amplitudes
rho = linspace(0, 1.1*rhoTarget, 201);
Omega = ssm.compute_omega(rho);

% Physical amplitudes
[x, xk] = ssm.compute_x(rho, outDof);

%% Plot

% % Plot reduced response
% figure
% hold on; grid on; box on; axis square tight;
% plot(Omega, rho, 'k', 'LineWidth', 2)
% xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
% ylabel('$\rho$ [-]', 'Interpreter', 'latex')
% title('Reduced amplitude', 'Interpreter', 'latex')
% set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
% 
% % Plot manifold
% figure
% hold on; grid on; box on; axis square tight; view(3);
% surf(rho.' .* cos(ssm.theta), rho.' .* sin(ssm.theta), xk, ...
%     'FaceColor', [0.0745 0.6235 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7)
% [~, idx] = min(abs(rho-rhoTarget));
% plot3(rho(idx) * cos(ssm.theta), rho(idx) * sin(ssm.theta), xk(idx, :), 'r', 'LineWidth', 2);
% L = light;
% L.Position = [0 0 1];
% lighting gouraud
% xlabel('$\rho cos(\theta)$', 'Interpreter', 'latex')
% ylabel('$\rho sin(\theta)$', 'Interpreter', 'latex')
% zlabel('$x$ [m]', 'Interpreter', 'latex')
% title('Manifold', 'Interpreter', 'latex')
% set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
% 
% % Plot physical response
% figure
% hold on; grid on; box on; axis square tight;
% plot(Omega, x, 'k', 'LineWidth', 2)
% xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
% ylabel('$x$ [m]', 'Interpreter', 'latex')
% title('Physical amplitude', 'Interpreter', 'latex')
% set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
% if exist('test/testSSM.mat','file') > 0
%     plot(testSSM.Omega, testSSM.ZoutNorm, '--', 'LineWidth', 2) % SSM validation
% end

%% Target point

% Compute target points
[xTarget, xkTarget] = ssm.compute_x(rhoTarget, outDof);
OmegaTarget = ssm.compute_omega(rhoTarget);

% Compute the derivative of rho
drhoDd = ssm.compute_rho_derivative(rhoTarget, outDof, xkTarget);

% Compute the derivative of Omega
dOmegaDd = ssm.compute_omega_derivative(rhoTarget, drhoDd);

% Plot physical response
% plot(OmegaTarget, xTarget, '.r', 'MarkerSize', 20)

% Initialize the adjint sensitivity object
sens = SensitivityBackbone(ssm);

%% Partial derivatives with respect to rho

% Partial derivative of Omega
pOmega_rho = sens.partial_derivative_Omega_rho(rhoTarget);

% Partial derivative of x
px_rho = sens.partial_derivative_x_rho(rhoTarget, xkTarget, xTarget, outDof);

% Adjoint equation
rhoAdj = sens.adjoint_equation_rho(rhoTarget, xkTarget, xTarget, outDof);

%% Partial derivatives with respect to w

% Partial derivative of x
px_w = sens.partial_derivative_x_w(rhoTarget, xkTarget, xTarget, outDof);
% fprintf('px_w:\n')
% for i = 2:length(px_w)
%     fprintf('\t%.4e\n', ...
%         norm(px_w(i).coeffs - conj(px_w(i).coeffs(:, end:-1:1))))
% end

% Adjoint equation
wAdj = sens.adjoint_equation_w(rhoTarget, xkTarget, xTarget, outDof, rhoAdj);
% fprintf('wAdj:\n')
% for i = 2:length(wAdj)
%     fprintf('\t%.4e\n', ...
%         norm(wAdj(i).coeffs - conj(wAdj(i).coeffs(:, end:-1:1))))
% end

%% Partial derivatives with respect to phi

% Partial derivative of Omega
[pOmega_phi, pCohom_phi] = sens.partial_derivative_Omega_phi(rhoTarget, wAdj);

% Partial derivative of x
px_phi = sens.partial_derivative_x_phi(rhoTarget, xkTarget, xTarget, outDof);

%% Partial derivatives with respect to omega

% Partial derivative of Omega
[pOmega_omega, pCohom_omega] = sens.partial_derivative_Omega_omega(rhoTarget, wAdj);

% Partial derivative of x
px_omega = 0;

% Adjoint equation
[omegaAdj, phiAdj] = sens.adjoint_equation_omega_phi(rhoTarget, xkTarget, xTarget, outDof, rhoAdj, wAdj);

%% Adjoint sensitivity

% Partial derivative of Omega
[pOmega_mu, pCohom_mu] = sens.partial_derivative_Omega_mu(rhoTarget, wAdj);

% Adjoint sensitivity of Omega
dOmega = sens.adjoint_sensitivity_Omega(rhoTarget, xkTarget, xTarget, outDof);

%% Sensitivity check

% Perturbation step
dMu = 1e-7;

% Collect design variables
mu0 = [m, k, k2, k3];

% Header
fprintf('\nAnalytical, finite difference, absolute error, relative error\n')

% Loop over design variables
for dv = 1:sys.ndv
    fprintf('\nDesign variable %d\n', dv)

    % Perturb design variable
    mu = mu0;
    mu(dv) = mu(dv) + dMu;

    % Extract parameters
    mPert = mu(1);
    kPert = mu(2);
    k2Pert = mu(3);
    k3Pert = mu(4);

    % Create system
    sysPert = MechSystem();
    sysPert.build_uniform_spring_mass_chain(n, mPert, kPert, k2Pert, k3Pert);

    % Modal analysis
    sysPert.modal_analysis('modeIndex', modeIdx);
    sysPert.set_damping_rayleigh([alphaR, betaR])
    % sysPert.set_damping_qfactor(Qfactor);

    % SSM
    ssmPert = MechSSM(sysPert);
    ssmPert.compute_ssm(maxOrder);
    rhoPert = ssmPert.compute_rho_from_x(xTarget, outDof);
    OmegaPert = ssmPert.compute_omega(rhoPert);

    % Modal analysis sensitivity
    sensAdj = domega(dv);
    sensFD = (sysPert.omega - sys.omega) / dMu;
    test_sensitivity_value_function('omega', sensAdj, sensFD);
    sensAdj = sys.pcsi_omega * domega(dv);
    sensFD = (sysPert.csi - sys.csi) / dMu;
    test_sensitivity_value_function('csi', sensAdj, sensFD);
    sensAdj = domegaDamped(dv);
    sensFD = (sysPert.omegaDamped - sys.omegaDamped) / dMu;
    test_sensitivity_value_function('omegaD', sensAdj, sensFD);

    % Sensitivity of the backbone
    sensAdj = dOmega(dv);
    sensFD = (OmegaPert - OmegaTarget) / dMu;
    test_sensitivity_value_function('Omega', sensAdj, sensFD);
end

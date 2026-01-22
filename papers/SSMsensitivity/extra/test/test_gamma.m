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
Qfactor = 20;

% Create system
sys = MechSystem();
sys.build_uniform_spring_mass_chain(n, m, k, k2, k3, 'computePartialDerivatives', true);

% Index of the master mode
modeIdx = 1;

% Modal analysis
sys.modal_analysis('modeIndex', modeIdx, 'computePartialDerivatives', true);
sys.set_damping_rayleigh([alphaR, betaR], 'computePartialDerivatives', true);
% sys.set_damping_qfactor(Qfactor, 'computePartialDerivatives', true);
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

% Reduced backbone coefficient
gammaOrder = 7;
gamma = imag(ssm.R(gammaOrder).coeffs(1));

%% Sensitivity

% Extract direct differentiation sensitivity
dgammaDd = zeros(sys.ndv, 1);
for dv = 1:sys.ndv
    dgammaDd(dv) = 1i/2 * (ssm.R(gammaOrder).dv(dv).coeffs(2) - ssm.R(gammaOrder).dv(dv).coeffs(1));
end

% Compute adjoint sensitivity
sens = SensitivityGamma(ssm, gammaOrder);
dgammaAdj = sens.adjoint_sensitivity;

% Compare sensitivity
errAbs = norm(dgammaAdj - dgammaDd);
errRel = norm(dgammaAdj ./ dgammaDd - 1);

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
    sysPert.set_damping_rayleigh([alphaR, betaR]);
    % sysPert.set_damping_qfactor(Qfactor);

    % SSM
    ssmPert = MechSSM(sysPert);
    ssmPert.compute_ssm(maxOrder);
    gammaPert = imag(ssmPert.R(gammaOrder).coeffs(1));

    % Sensitivity of the backbone
    sensAdj = dgammaAdj(dv);
    sensFD = (gammaPert - gamma) / dMu;
    test_sensitivity_value_function('Omega', sensAdj, sensFD);
end

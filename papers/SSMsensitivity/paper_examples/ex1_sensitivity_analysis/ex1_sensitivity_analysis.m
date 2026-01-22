clear; clc; close all;
addpath(genpath('../../src'))

mblue = '#0072BD';
mbluelight = '#9CCDEF';
mred = '#D95319';
mredlight = '#F4C6B2';
myellow = '#EDB120';

%% Settings

% Parameters
m = 1;
k = 1;
k2 = 0.5;
k3 = 0.2;

% System
n = 1;
rayleigh = [0, 0.1];
modeIndex = 1;

% SSM
outDof = n+1;
maxOrder = 5;
xT = 1.5;

%% Analysis

% Create system
sys = MechSystem();
sys.build_uniform_spring_mass_chain(n, m, k, k2, k3, ...
    'freeEnd', true, 'computePartialDerivatives', true);

% Modal analysis
sys.modal_analysis('modeIndex', modeIndex, ...
    'computePartialDerivatives', true);
sys.set_damping_rayleigh(rayleigh, 'computePartialDerivatives', true)

% Create SSM
ssm = MechSSM(sys);
ssm.compute_ssm(maxOrder, 'storeCoefficients', true);

% Compute target points
nRhoSens = 11;
rhoT = ssm.compute_rho_from_x(xT, outDof);
rho = linspace(0, rhoT, nRhoSens);
[x0, xk] = ssm.compute_x(rho, outDof);
Omega0 = ssm.compute_omega(rho);

% Sensitivity
sens = SensitivityBackbone(ssm);
dOmega = zeros(sys.ndv, length(x0));
for i = 1:length(x0)
    if x0(i) == 0
        [~, domegaDamped] = sys.adjoint_sensitivity_omega;
        dOmega(:, i) = domegaDamped;
    else
        dOmega(:, i) = sens.adjoint_sensitivity_Omega(rho(i), xk(i, :), x0(i), outDof);
    end
end

% Compute backbone
nRho = 51;
rhoT = ssm.compute_rho_from_x(xT, outDof);
rho = linspace(0, rhoT, nRho);
[x, xk] = ssm.compute_x(rho, outDof);
Omega = ssm.compute_omega(rho);

%% Plot

% Create figure
fig = figure;
tl = tiledlayout('flow', 'TileSpacing', 'compact');

% Parameters
mu0 = [m, k, k2, k3];
idxPert = [1, 2, 3, 4]; % we do not perturb k2
pertStep = [0.01 * m, 0.01 * k, 0.03 * k2, 0.03 * k3];
title_str = {'$m \pm \delta m$', '$k \pm \delta k$', '$k_2 \pm \delta k_2$', '$k_3 \pm \delta k_3$'};

% Loop over parameters
for i = 1:length(idxPert)
    % Plot nominal backbone
    nexttile(tl);
    hold on; grid on; box on; axis tight;
    plot(Omega / sys.omegaDamped, x, 'k', 'LineWidth', 2)

    % Perturbation
    mu = mu0;
    mu(idxPert(i)) = mu(idxPert(i)) - pertStep(1, i);
    [OmegaPertL, xPertL] = compute_backbone(mu, n, rayleigh, modeIndex, outDof, maxOrder, xT, nRho);
    plot(OmegaPertL / sys.omegaDamped, xPertL, '-', 'Color', mredlight, 'LineWidth', 2)

    mu = mu0;
    mu(idxPert(i)) = mu(idxPert(i)) + pertStep(1, i);
    [OmegaPertU, xPertU] = compute_backbone(mu, n, rayleigh, modeIndex, outDof, maxOrder, xT, nRho);
    plot(OmegaPertU / sys.omegaDamped, xPertU, '-', 'Color', mredlight, 'LineWidth', 2)

    OmegaPert = Omega0 + dOmega(idxPert(i), :) * pertStep(1, i) .* [-1; 1];
    plot(OmegaPert / sys.omegaDamped, x0, 'x', 'Color', mred, 'LineWidth', 2, 'MarkerSize', 10)
    
    % Decorations
    xlabel('$\Omega / \omega_0$', 'Interpreter', 'latex')
    ylabel('$x$', 'Interpreter', 'latex')
    title(title_str{i}, 'Interpreter', 'latex')
    set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
    xlim(1 + 0.02*[-1, 1])
end

fig.Position(1:2) = 0.5*fig.Position(1:2);
fig.Position(3:4) = [560*2 420*1.5];

%% Functions
function [Omega, x] = compute_backbone(mu, n, rayleigh, modeIndex, outDof, maxOrder, xT, nRho)
    % Extract parameters
    m = mu(1);
    k = mu(2);
    k2 = mu(3);
    k3 = mu(4);

    % Create system
    sys = MechSystem();
    sys.build_uniform_spring_mass_chain(n, m, k, k2, k3, rayleigh, ...
        'freeEnd', true);
    
    % Modal analysis
    sys.modal_analysis('modeIndex', modeIndex);
    sys.set_damping_rayleigh(rayleigh)
    
    % Create SSM
    ssm = MechSSM(sys);
    ssm.compute_ssm(maxOrder);
    
    % Compute target points
    rhoT = ssm.compute_rho_from_x(xT, outDof);
    rho = linspace(0, rhoT, nRho);
    x = ssm.compute_x(rho, outDof);
    Omega = ssm.compute_omega(rho);
end

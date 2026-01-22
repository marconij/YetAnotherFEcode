clear; clc; close all;
addpath(genpath('../../src'));

%% Settings
m = 1;
k = 1;
k2 = 0.01;
k3 = 0.01;

% Rayleigh damping
rayleigh = [0, 0.1];

% SSM
order = 3;
outDof = 1;
inDof = 1;
xT = 8;

% Forcing levels
epsilon = [0.1, 0.4, 0.6];

%% System

% Build system
sys = MechSystem;
sys.build_duffing_system(m, k, k2, k3, ...
    'computePartialDerivatives', true);

% Modal analysis
sys.modal_analysis('nModes', 1, 'modeIndex', 1, ...
    'computePartialDerivatives', true);
sys.set_damping_rayleigh(rayleigh, 'computePartialDerivatives', true);
% sys.set_damping_qfactor(10, 'computePartialDerivatives', true);
sys.modal_analysis_sensitivity;

%% SSM

% Build SSM
ssm = MechFRC(sys);
ssm.compute_ssm(order, 'computeSensitivity', true);

% Non-autonomous reduced dynamics
ssm.compute_S(inDof);

%% FRCs

% Initialize vectors
nEps = length(epsilon);
nRho = 51;
omegaPk = zeros(1, nEps);
xPk = zeros(1, nEps);
rhoPk = zeros(1, nEps);
omegaFr = zeros(2*nRho, nEps);
xFr = zeros(2*nRho, nEps);
psiFr = zeros(2*nRho, nEps);
omegaBb = zeros(nRho, nEps);
xBb = zeros(nRho, nEps);

% Loop over forcing levels
rho0 = 0.1;
for i = 1:length(epsilon)
    % FR peak
    [omegaPk(i), xPk(i), rhoPk(i)] = ssm.compute_frc_peak(outDof, epsilon(i), rho0);
    rho0 = rhoPk(i);
    
    % Frequency response
    [omegaFr(:, i), xFr(:, i), psiFr(:, i), omegaBb(:, i), xBb(:, i)] = ssm.compute_frc(rhoPk(i), outDof, epsilon(i), 'nRho', nRho);
end

%% Plot

% Figure
fig = figure;
tl = tiledlayout(2, 1, 'TileSpacing', 'compact');
nexttile(tl, 1)
hold on; grid on; box on; axis tight;
nexttile(tl, 2)
hold on; grid on; box on; axis tight;

% Loop over forcing levels
rho0 = 0.1;
for i = 1:length(epsilon)
    % Amplitude
    nexttile(tl, 1)
    plot(omegaFr(:, i) / sys.omega, xFr(:, i), '.-', 'LineWidth', 2, 'MarkerSize', 20, 'DisplayName', num2str(epsilon(i), '$\\epsilon = %.1f$'))

    % Phase
    nexttile(tl, 2)
    plot(omegaFr(:, i) / sys.omega, psiFr(:, i) / pi, '.-', 'LineWidth', 2, 'MarkerSize', 20)
end

% Decorations
nexttile(tl, 1)
plot(omegaBb(:, end) / sys.omega, xBb(:, end), '--k', 'LineWidth', 2, 'DisplayName', 'Backbone')
xlabel('$\Omega / \omega_0$ [-]')
ylabel('$|x_{out}|$ [m]')
legend('Location', 'northwest')
xlim([0.5, 1.5])

nexttile(tl, 2)
yline(-1/2, '--k', 'LineWidth', 2)
xlabel('$\Omega / \omega_0$ [-]')
ylabel('$\psi / \pi$ [-]')
xlim([0.5, 1.5])

% Size
fig.Position(2) = 0.1*fig.Position(2);
fig.Position(4) = 2*fig.Position(4);
% exportgraphics(fig, 'fr.png', 'Resolution', 300)

%% Sensitivity analysis

% Compute sensitivity
mu0 = [m, k, k2, k3];
eps0 = 0.5;
[omegaPk0, xPk0, rhoPk0, domegaPk, dxPk, drhoPk] = fun(mu0, rayleigh, order, inDof, outDof, eps0, true);

% Concatenate
J0 = [omegaPk0, xPk0, rhoPk0];
dJ = [domegaPk, dxPk, drhoPk];

% Finite differences
delta = [1e-3, 1e-4, 1e-5, 1e-6];
dJ_fd = zeros(sys.ndv, length(delta), length(J0));
errAbs = zeros(sys.ndv, length(delta), length(J0));
errRel = zeros(sys.ndv, length(delta), length(J0));
for i = 1:sys.ndv
    % Loop over perturbations
    for j = 1:length(delta)
        % Perturb
        muP = mu0;
        muP(i) = muP(i) + delta(j);

        % Compute
        [omegaPkP, xPkP, rhoPkP] = fun(muP, rayleigh, order, inDof, outDof, eps0, false);
        JP = [omegaPkP, xPkP, rhoPkP];

        % Finite differences
        for k = 1:length(J0)
            dJ_fd(i, j, k) = (JP(k) - J0(k)) / delta(j);
            errAbs(i, j, k) = abs(dJ_fd(i, j, k) - dJ(i, k));
            errRel(i, j, k) = errAbs(i, j, k) / abs(dJ(i, k));
        end
    end

    % Display
    fprintf('Variable %d: %12.4e %12.4e %12.4e\n', i, dJ(i, :))
    for k = 1:length(J0)
        fprintf('\tJ%d:\n', k)
        for j = 1:length(delta)
            fprintf('\t\tStep %.2e: %12.4e %12.4e %12.4e\n', delta(j), dJ_fd(i, j, k), errAbs(i, j, k), errRel(i, j, k))
        end
    end
end

% Plot
title_str = {'$\Omega$', '$x$', '$\rho$'};
fig = figure;
tiledlayout(1, length(J0), 'TileSpacing', 'compact')
for k = 1:length(J0)
    nexttile
    semilogy(1:sys.ndv, errRel(:, 1, k), '.-', 'LineWidth', 2, 'MarkerSize', 20, 'DisplayName', num2str(delta(1), '%.0e'))
    hold on; grid on; box on; axis tight;
    for j = 2:length(delta)
        semilogy(1:sys.ndv, errRel(:, j, k), '.-', 'LineWidth', 2, 'MarkerSize', 20, 'DisplayName', num2str(delta(j), '%.0e'))
    end
    xlabel('Design variable', 'Interpreter', 'latex')
    ylabel('Relative error', 'Interpreter', 'latex')
    title(['Sensitivity of ', title_str{k}], 'Interpreter', 'latex')
    legend show
    set(gca, 'FontSize', 18, 'TickLabelInterpreter', 'latex')
end
fig.Position(2) = 0.1*fig.Position(2);
fig.Position(3) = 3*fig.Position(3);
fig.Position(4) = 1.5*fig.Position(4);
% exportgraphics(fig, 'fr_sensitivity.png', 'Resolution', 300)




















%% Functions

function [omegaPk, xPk, rhoPk, domegaPk, dxPk, drhoPk] = fun(mu, rayleigh, order, inDof, outDof, epsilon, computeSensitivity)
    % Extract parameters
    m = mu(1);
    k = mu(2);
    k2 = mu(3);
    k3 = mu(4);

    % Create system
    sys = MechSystem;
    sys.build_duffing_system(m, k, k2, k3, ...
        'computePartialDerivatives', computeSensitivity);
    
    % Modal analysis
    sys.modal_analysis('nModes', 1, 'modeIndex', 1, ...
        'computePartialDerivatives', computeSensitivity);
    sys.set_damping_rayleigh(rayleigh, 'computePartialDerivatives', computeSensitivity);
    % sys.set_damping_qfactor(10, 'computePartialDerivatives', computeSensitivity);
    if computeSensitivity
        sys.modal_analysis_sensitivity;
    end

    % Build SSM
    ssm = MechFRC(sys);
    ssm.compute_ssm(order, 'computeSensitivity', computeSensitivity);

    % Non-autonomous reduced dynamics
    ssm.compute_S(inDof);

    % Compute peak value
    rho0 = 0.1;
    [omegaPk, xPk, rhoPk, xiPk] = ssm.compute_frc_peak(outDof, epsilon, rho0);
    if rhoPk < 0
        disp('rhoPk < 0')
    end
    
    % Compute sensitivity
    if computeSensitivity
        [domegaPk, dxPk, drhoPk] = ssm.compute_frc_peak_derivative(outDof, epsilon, rhoPk, xPk, xiPk);
    else
        domegaPk = [];
        dxPk = [];
        drhoPk = [];
    end
end

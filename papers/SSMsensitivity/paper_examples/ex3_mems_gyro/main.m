% Optimization of a MEMS gyroscope
clear; close all; clc;

addpath(genpath('../../src/'))

%% Settings and initial conditions
mu(1) =  20; % delta length (drive beams)
mu(2) = -10; % delta length (sense beams)
mu(3) =   1; % delta width (drive beams)
mu(4) =   1; % delta width (sense beams)
mu(5) =   0; % delta width (sense beam connections)

mu(6:10)  = [10 10 0 0 0]; % first 5 harmonics for the drive beams (top)
mu(11:15) = [10 10 0 0 0]; % first 5 harmonics for the drive beams (bottom)

% Rayleigh damping
rayleigh = [0.01, 0];

% Model constructor
yafec_mpc_constructor = @(x) build_model(x);

% Optimization arguments
args.yafec_mpc_constructor = yafec_mpc_constructor;
args.rayleigh = rayleigh;
args.order = 3;
args.maxOrder = 9;
args.errTol = 1e-1;

%% Preliminary analysis @ mu = 0

% Build system
sys = MechSystem;
sys.build_yafec_mpc_model(mu*0, yafec_mpc_constructor);

% Plot geometry
figure
plot_geometry(sys)
title('Initial geometry with $p = 0$', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Modal analysis
[D, V] = sys.modal_analysis('nModes', 10, 'modeIndex', 1);
sys.set_damping_rayleigh(rayleigh);

% Drive mode
phiDrive = sys.phi;
fDrive = sys.omega / (2*pi) / 1e3;

% Sense mode
phiSense = V(:, 2);
phiSense = phiSense / sqrt(phiSense.' * sys.M * phiSense); % mass normalization
fSense = sqrt(D(2, 2)) / (2*pi) / 1e3;

% Plot mode shapes
scale = 10;
figure
tiledlayout('flow', 'TileSpacing', 'compact')
nexttile
plot_geometry(sys, 'u', phiDrive, 'scale', scale)
axis equal tight;
title(sprintf('Drive mode, f = %.2f kHz', fDrive), 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
nexttile
plot_geometry(sys, 'u', phiSense, 'scale', scale)
axis equal tight;
title(sprintf('Sense mode, f = %.2f kHz', fSense), 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% SSM
ssm = MechSSM(sys);
ssm.compute_ssm(args.order);

% Target displacement
args.outDof = sys.n - 1; % y displacement of the drive frame
args.xT = 3; % [um]

% Backbone
rhoT = ssm.compute_rho_from_x(args.xT, args.outDof);
rho = linspace(0, 1.1 * rhoT, 51);
OmegaInit = ssm.compute_omega(rho);
xInit = ssm.compute_x(rho, args.outDof);
    
% Plot backbone
figure
hold on; grid on; box on; axis square tight;
plot(OmegaInit / 2 / pi / 1000, xInit, 'k', 'LineWidth', 2)
xlabel('$\Omega$ [kHz]', 'Interpreter', 'latex')
ylabel('$x$ [$\mu$m]', 'Interpreter', 'latex')
title('Initial backbone', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

%% Optimization settings

% Initial conditions
mu0 = mu*0;

% Bounds
muL = -[25, 25, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5];
muU =  [25, 25, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5];

% Arguments
args.modeRefDrive = phiDrive;
args.modeRefSense = phiSense;
args.omegaDriveT = 28 * 1e3 * 2*pi; % [rad/s]
args.omegaSenseT = 26 * 1e3 * 2*pi; % [rad/s]
args.OmegaT = 26 * 1e3 * 2*pi; % [rad/s]

% Options
options = optimoptions('fmincon', 'Display', 'none');
options.Algorithm = "interior-point";
options.EnableFeasibilityMode = true;
options.SpecifyObjectiveGradient = true;
options.SpecifyConstraintGradient = true;
% options.ConstraintTolerance = 1e-3;
% options.StepTolerance = 1e-3;
options.Display = "final";

% Optimize
t0 = tic;
% [muOpt, orderOpt, orderHistory] = run_fmincon(mu0, [], [], [], [], muL, muU, args, options);
[muOpt, orderOpt, orderHistory] = run_fmincon_no_obj(mu0, [], [], [], [], muL, muU, args, options);
topt = toc(t0);
fprintf('Optimization time: %d min %.0f s\n\n', floor(topt/60), rem(topt,60))

%% Optimization results

% Error and order
fig = figure;
% fig.Position(3) = 560*2;
yyaxis left
semilogy(orderHistory.iter, orderHistory.err, '.-', 'LineWidth', 2, 'MarkerSize', 20)
hold on; grid on; box on; axis tight;
yline(args.errTol .* [1, 0.01], '--', 'LineWidth', 2)
ylabel('Error measure', 'Interpreter', 'latex')
ylim([1e-4, 1e0])
yyaxis right
plot(orderHistory.iter, orderHistory.order, '.-', 'LineWidth', 2, 'MarkerSize', 20)
ylabel('Order', 'Interpreter', 'latex')
ylim([2, args.maxOrder + 1])
xlabel('Iterations', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Build system
sys = MechSystem;
sys.build_yafec_mpc_model(muOpt, yafec_mpc_constructor);

% Plot geometry
figure
plot_geometry(sys)
title('Optimal geometry', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Modal analysis
[D, V] = sys.modal_analysis('nModes', 10, 'modeRef', args.modeRefDrive);
sys.set_damping_rayleigh(rayleigh);

% Drive mode
phiDrive = sys.phi;
fDrive = sys.omega / (2*pi) / 1e3;

% Sense mode
senseIndex = modal_assurance_criterion(V, args.modeRefSense);
phiSense = V(:, senseIndex);
phiSense = phiSense / sqrt(phiSense.' * sys.M * phiSense); % mass normalization
fSense = sqrt(D(senseIndex)) / (2*pi) / 1e3;

% Plot mode shapes
scale = 10;
figure
tiledlayout('flow', 'TileSpacing', 'compact')
nexttile
plot_geometry(sys, 'u', phiDrive, 'scale', scale)
axis equal tight;
title(sprintf('Drive mode, f = %.2f kHz', fDrive), 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
nexttile
plot_geometry(sys, 'u', phiSense, 'scale', scale)
axis equal tight;
title(sprintf('Sense mode, f = %.2f kHz', fSense), 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% SSM
ssm = MechSSM(sys);
ssm.compute_ssm(orderOpt);

% Backbone
rhoT = ssm.compute_rho_from_x(max(args.xT), args.outDof);
rho = linspace(0, 1.1 * rhoT, 51);
Omega = ssm.compute_omega(rho);
[x, xk] = ssm.compute_x(rho, args.outDof);

% Plot backbone
figure
hold on; grid on; box on; axis square tight;
plot(Omega / 2 / pi / 1000, x, 'k', 'LineWidth', 2)
plot(args.OmegaT / 2 / pi / 1000, args.xT, '.r', 'MarkerSize', 20)
xlabel('$\Omega$ [kHz]', 'Interpreter', 'latex')
ylabel('$x$ [$\mu$m]', 'Interpreter', 'latex')
legend('Backbone', 'Target Points', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Plot manifold
figure
hold on; grid on; box on; axis tight;
view(3)
surf(rho.' * cos(ssm.theta), rho.' * sin(ssm.theta), xk, ...
    'FaceColor', [0.0745 0.6235 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7)
L = light;
L.Position = [0 0 1];
lighting gouraud
axis square tight
xlabel('$\rho \cos(\theta)$', 'Interpreter', 'latex')
ylabel('$\rho \sin(\theta)$', 'Interpreter', 'latex')
zlabel('$x$', 'Interpreter', 'latex')
title('Spectral Submanifold', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
view(45,15)

%% Convergenc analysis

% Initialize figure
figure
hold on; grid on; box on; axis square tight;
    
% Plot initial backbone
plot(OmegaInit / 2 / pi / 1000, xInit, 'k--', 'LineWidth', 2, 'DisplayName', 'Initial')

% Loop over orders
colorsLines = lines(7);
orderMin = max(orderOpt, 3);
orders = (0:2:4) + orderMin;
for i = 1:length(orders)
    % Compute SSM
    order = orders(i);
    ssm.compute_ssm(order);

    % Backbone
    rhoT = ssm.compute_rho_from_x(max(args.xT), args.outDof);
    rho = linspace(0, 1.1 * rhoT, 51);
    Omega = ssm.compute_omega(rho);
    [x, xk] = ssm.compute_x(rho, args.outDof);
    
    % Plot backbone
    plot(Omega / 2 / pi / 1000, x, 'LineWidth', 2, 'Color', colorsLines(i, :), ...
        'DisplayName', num2str(order, '$\\mathcal{O}$(%d)'))
end

% Decorations
plot(args.OmegaT / 2 / pi / 1000, args.xT, '.r', 'MarkerSize', 20, ...
    'DisplayName', 'Target point')
xlabel('$\Omega$ [kHz]', 'Interpreter', 'latex')
ylabel('$x$ [$\mu$m]', 'Interpreter', 'latex')
title('Convergence analysis', 'Interpreter', 'latex')
legend('Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

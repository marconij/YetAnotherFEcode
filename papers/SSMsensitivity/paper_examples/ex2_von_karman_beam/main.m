% Optimization of a von Karman beam
clear; clc; close all;
addpath(genpath('../../src'));

%% Settings
% Parameters
h = 10e-3;
L = 1;
mu = [0, 0, h, L];
yafec_assembly_constructor = @(x) define_beam_model(x);

% System
sys = MechSystem;
sys.build_yafec_model(mu, yafec_assembly_constructor);

% Create geometry
nodes = sys.yafecAssembly.Mesh.nodes;
nx = nodes(:, 1) * 1e3;
ny = nodes(:, 2) * 1e3;
sx = [nx; flip(nx)];
sy = [ny - mu(3)/2*1e3; flip(ny) + mu(3)/2*1e3];

% Plot geometry
figure
hold on; grid on; box on; axis tight;
patch(sx, sy, [156, 205, 239] / 255, 'FaceAlpha', 1)
plot(nx, ny, 'k.-', 'LineWidth', 2, 'MarkerSize', 20)
xlabel('$x$ [mm]', 'Interpreter', 'latex')
ylabel('$y$ [mm]', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% System settings
args.yafec_assembly_constructor = yafec_assembly_constructor;

% Initial eigenvalue analysis to identify the shape of the first mode
sys.modal_analysis('modeIndex', 1);
args.modeRef = sys.phi;
args.omega0 = sys.omegaDamped;

% Output DOF: vertical displacement at the middle of the beam
args.outDof = (sys.yafecAssembly.Mesh.nElements / 2) * 3 + 2;

% SSM settings
args.order = 3;
args.maxOrder = 11;
args.tol = 1e-1;

%% Optimization problem
% Define the target points as frequency - amplitude pairs.
args.OmegaTarget = [1, 0.95] * args.omega0;
args.xTarget = [0.2, 0.4] * h;

% Specify if the frequency - amplitude pairs are treated as:
% equality constraint (0)
% upper inequality constraint (1)
% lower inequality constraint (-1).
args.coeffTarget = [0, 0];

% Define the initial conditions for the design variables.
mu0 = mu;

% Define the lower and upper bounds for the design variables.
muL = [  0,   0, h/10, 0.5];
muU = [2*h, 2*h, 10*h, 1.5];

% Define the linear equality and inequality constraints (if any).
A = []; b = [];
Aeq = []; beq = [];

%% Optimizer options
options = optimoptions('fmincon', 'Display', 'iter');
options.Algorithm = "interior-point";
options.EnableFeasibilityMode = true;
options.SpecifyObjectiveGradient = true;
options.SpecifyConstraintGradient = true;
options.ConstraintTolerance = 1e-3;
options.StepTolerance = 1e-3;
options.Display = "final";

%% Solve the optimization problem
% Start timer
tStart = tic;

% Solve
[mu, orderNew, history] = run_fmincon(mu0, A, b, Aeq, beq, muL, muU, args, options);

% Stop timer and print results
tElapsed = toc(tStart);
fprintf('Elapsed time is %.2f seconds.\n', tElapsed)

%% Post processing

% System
sysOpt = MechSystem;
sysOpt.build_yafec_model(mu, @(x) define_beam_model(x));

% Create geometry
nodes = sysOpt.yafecAssembly.Mesh.nodes;
nx = nodes(:, 1) * 1e3;
ny = nodes(:, 2) * 1e3;
sx = [nx; flip(nx)];
sy = [ny - mu(3)/2*1e3; flip(ny) + mu(3)/2*1e3];

% Plot geometry
figure
hold on; grid on; box on; axis tight;
patch(sx, sy, [156, 205, 239] / 255, 'FaceAlpha', 1)
plot(nx, ny, 'k.-', 'LineWidth', 2, 'MarkerSize', 20)
xlabel('$x$ [mm]', 'Interpreter', 'latex')
ylabel('$y$ [mm]', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Modal analysis
sysOpt.modal_analysis('modeRef', args.modeRef);

% Create SSM
ssmOpt = MechSSM(sysOpt);
ssmOpt.compute_ssm(orderNew);

% Create backbone
xMax = max(args.xTarget);
rhoMax = ssmOpt.compute_rho_from_x(xMax, args.outDof);
rho = linspace(0, rhoMax, 200);
Omega = ssmOpt.compute_omega(rho);
x = ssmOpt.compute_x(rho, args.outDof);

% Plot
figure
hold on; grid on; box on; axis square tight;
plot(Omega, x*1e3, 'k', 'LineWidth', 2)
plot(args.OmegaTarget, args.xTarget*1e3, '.r', 'MarkerSize', 20)
xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
ylabel('$x$ [mm]', 'Interpreter', 'latex')
title('Physical amplitude', 'Interpreter', 'latex')
legend('Backbone', 'Target Points', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

%% Create animation

% Backbone evolution
fig_bb = figure('Color', 'white');
filename = 'backbone.gif';
delay_time = 0.5;
for i = 1:length(history)
    % Create plot
    cla
    hold on; grid on; box on; axis square tight;
    plot(history(i).Omega, history(i).x*1e3, 'k', 'LineWidth', 2)
    plot(args.OmegaTarget, args.xTarget*1e3, '.r', 'MarkerSize', 20)
    xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
    ylabel('$x$ [mm]', 'Interpreter', 'latex')
    title(num2str(i - 1, 'Iteration %d'), 'Interpreter', 'latex')
    set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
    xlim([200, 250])
    hold off;
    drawnow

    % Save gif
    frame = getframe(fig_bb);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time, 'BackgroundColor', 0);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
end

% Geometry evolution
fig_beam = figure('Color', 'white');
filename = 'beam.gif';
delay_time = 0.5;
for i = 1:length(history)
    % System
    sys = MechSystem;
    sys.build_yafec_model(history(i).mu, @(x) define_beam_model(x));
    
    % Create geometry
    nodes = sys.yafecAssembly.Mesh.nodes;
    nx = nodes(:, 1) * 1e3;
    ny = nodes(:, 2) * 1e3;
    sx = [nx; flip(nx)];
    sy = [ny - history(i).mu(3)/2*1e3; flip(ny) + history(i).mu(3)/2*1e3];
    
    % Plot geometry
    cla
    hold on; grid on; box on; axis tight;
    patch(sx, sy, [156, 205, 239] / 255, 'FaceAlpha', 1)
    plot(nx, ny, 'k.-', 'LineWidth', 2, 'MarkerSize', 20)
    xlabel('$x$ [mm]', 'Interpreter', 'latex')
    ylabel('$y$ [mm]', 'Interpreter', 'latex')
    title(num2str(i - 1, 'Iteration %d'), 'Interpreter', 'latex')
    set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
    xlim([0, 1200])
    ylim([-22, 22])

    hold off;
    drawnow

    % Save gif
    frame = getframe(fig_beam);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time, 'BackgroundColor', 0);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
end

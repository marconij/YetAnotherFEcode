clear; clc; close all;
addpath(genpath('../../src'))

% Plot settings
pretty_plot;

%% Settings

% Parameters
n = 2;
m = 1; % kg
k = 3; % N/m
c1 = 0.03; % N s/m
c2 = sqrt(3) * 0.03; % N s/m
kappa = 0.4; % N/m^3
alpha = -0.6; % N (s/m)^3

% Rayleigh dampin
rayleigh = [(c1 - c2)/m, c2/k];

% Forcing
inDof = 1;
outDof = 1;

% SSM settings
maxOrder = 3;

%% SSM analysis

% Create system
sys = MechSystem();
sys.build_shaw_pierre_system(n, m, k, 0, kappa, inDof);

% Modal analysis
sys.modal_analysis('modeIndex', 1);
sys.set_damping_rayleigh(rayleigh);
sys.set_damping_nonlinear(0, alpha/kappa);

% Create SSM
ssm = MechFRC(sys);
ssm.compute_ssm(maxOrder);
ssm.compute_S(inDof);

% Extract coefficients
gamma1 = ssm.R(1).coeffs(1);
gamma3 = ssm.R(3).coeffs(1);
S0 = ssm.S0;

% Merging value
epsLim = 1 / abs(S0) * sqrt(-4*real(gamma1)^3 / (27 * real(gamma3)));

% Forcings
epsT = epsLim * (0.5:0.1:1.5);
nEps = length(epsT);
epsPlotIdx = [1, find(epsT / epsLim == 1), nEps];

%% Compute frequency response

% Initialize vectors
nRho = 101;
omegaPk = zeros(1, nEps);
xPk = zeros(1, nEps);
rhoPk = zeros(1, nEps);
omegaFr = zeros(2*nRho, nEps);
xFr = zeros(2*nRho, nEps);
psiFr = zeros(2*nRho, nEps);
rhoFr = zeros(2*nRho, nEps);
omegaBb = zeros(nRho, nEps);
xBb = zeros(nRho, nEps);
rhoBb = zeros(nRho, nEps);
isStable = false(2*nRho, nEps);

% Loop over forcing levels
rho0 = 0.1;
for i = nEps:-1:1
    % FR peak
    [omegaPk(i), xPk(i), rhoPk(i)] = ssm.compute_frc_peak(outDof, epsT(i), rho0);
    rho0 = rhoPk(i);
    
    % Frequency response
    [omegaFr(:, i), xFr(:, i), psiFr(:, i), omegaBb(:, i), xBb(:, i), rhoBb(:, i)] = ssm.compute_frc(rhoPk(i), outDof, epsT(i), 'nRho', nRho);
    rhoFr(:, i) = [rhoBb(:, i); flip(rhoBb(:, i))];

    % Stability
    for j = 1:length(omegaFr(:, i))
        [~, isStable(j, i)] = ssm.compute_jacobian(rhoFr(j, i), real(omegaFr(j, i)));
    end
end

%% Frequency-phase-amplitude plot

% Limits on omega axis
omegaLim = (1 + [-1, 1]*0.02) * sys.omegaDamped;
psiLim = [-2*pi, 0];
rhoLim = [0, 1.1*max(rhoPk)];

% Plot frequency response
ssm.plot_frequency_response(epsT(epsPlotIdx), omegaLim, rhoLim, 'plotPhase', true);

% Loop over forcing levels
for epsIdx = epsPlotIdx
    ssm.plot_frequency_phase_amplitude(epsT(epsIdx), omegaLim, psiLim, [0, 1.1*rhoPk(epsIdx)]);
    % ssm.plot_frequency_phase_amplitude(epsT(epsIdx), omegaLim, psiLim, [0, 1.1*rhoPk(epsIdx)], ...
    %     'frcData', {[omegaFr(:, epsIdx), psiFr(:, epsIdx), rhoFr(:, epsIdx)]}, 'plot2D', true);
end

%% Stability plot

% Initialize figure
figure;
tiledlayout('flow', 'TileSpacing', 'compact')

% Loop over forcing levels
for i = epsPlotIdx
    % Extract this vectors
    thisOmega = omegaFr(:, i) / sys.omegaDamped;
    thisRho = rhoFr(:, i);
    thisPsi = psiFr(:, i);
    thisX = xFr(:, i);
    thisStable = isStable(:, i);

    % Saddle points
    sp = find(abs(diff(thisStable)) > 0);

    % Stable and unstable branches
    omegaStb = thisOmega; omegaStb(~thisStable) = NaN;
    rhoStb = thisRho; rhoStb(~thisStable) = NaN;
    omegaUnst = thisOmega; omegaUnst(thisStable) = NaN;
    rhoUnst = thisRho; rhoUnst(thisStable) = NaN;

    % Initialize plot
    nexttile
    hold on; grid on; box on; axis tight;
    plot(omegaStb, rhoStb, '-', 'LineWidth', 2, 'Color', mblue);
    plot(omegaUnst, rhoUnst, '-', 'LineWidth', 2, 'Color', mred);
    plot(thisOmega(sp), thisRho(sp), '.', 'MarkerSize', 20, 'MarkerFaceColor', myellow);

    % Decorations
    xlabel('$\Omega / \omega_0$ [-]')
    ylabel('$\rho$ [-]')
    xlim(omegaLim / sys.omegaDamped)
end

%% Frequency-phase-amplitude plot


%% Frequency-force-amplitude plot

% Define vectors
nPts = 101;
omegaV = linspace(omegaLim(1), omegaLim(2), nPts);
epsV = linspace(epsT(1), epsT(end), nPts);
rhoV = linspace(0, 1.1 * max(rhoPk, [], 'all'), nPts + 1); rhoV(1) = [];

% 3D grid
[omegaM, epsM, rhoM] = meshgrid(omegaV, epsV, rhoV);

% Define function
M = (ssm.compute_a(rhoM)).^2 + (ssm.compute_b(rhoM) - omegaM.*rhoM).^2 - (epsM * abs(ssm.S0)).^2;

% Plot
figure
hold on; grid on; box on; axis square tight;
patch(isosurface(omegaM / sys.omegaDamped, epsM / epsLim, rhoM, M, 0), ...
    'EdgeColor', 'none', 'FaceColor', '#0072BD', 'FaceAlpha', 0.5);
xlabel('$\Omega / \omega_0$ [-]')
ylabel('$\epsilon / \epsilon_m$ [-]')
zlabel('$\rho$ [-]')
xlim(minmax(omegaV) / sys.omegaDamped)
ylim(minmax(epsV) / epsLim)
zlim(minmax(rhoV))
view(3)

% Lighting and view
camlight;
lighting gouraud;

%% Animation

% Create grids
[omega3, eps3, rho3] = meshgrid(omegaV, epsV, rhoV);
[omega2, rho2] = meshgrid(omegaV, rhoV);

% Create level set functions
M2 = (ssm.compute_a(rho2)).^2 + (ssm.compute_b(rho2) - omega2.*rho2).^2;
M3 = (ssm.compute_a(rho3)).^2 + (ssm.compute_b(rho3) - omega3.*rho3).^2 - (eps3 * abs(ssm.S0)).^2;

% Create backbone
rhoBbP = rhoV;
omegaBbP = ssm.compute_b(rhoBbP) ./ rhoBbP;

% Create figure
figAnimation = figure('Color', 'w');
figAnimation.Position(3) = 2*figAnimation.Position(3);
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');

% Initialize frequency-force-amplitude plot
nexttile(tl, 1);
hold on; grid on; box on; axis square tight; view(3);

% Initialize frequency-amplitude plot
nexttile(tl, 2);
hold on; grid on; box on; axis square tight;

% Gif settings
fileName = 'isola_animation.gif';
delayTime = 0.1;
isFirstFrame = true;

% Loop over forcings
for i = 1:4:length(epsV)
    % Set 3D plot
    nexttile(tl, 1);
    cla;
    xlabel('$\Omega / \omega_0$ [-]')
    ylabel('$\epsilon / \epsilon_m$ [-]')
    zlabel('$\rho$ [-]')
    xlim(minmax(omegaV) / sys.omegaDamped)
    ylim(minmax(epsV) / epsLim)
    zlim(minmax(rhoV))
    
    % Lighting and view
    camlight;
    lighting gouraud;

    % Plot frequency-force-amplitude surface
    patch(isosurface(omega3 / sys.omegaDamped, eps3 / epsLim, rho3, M3, 0), ...
        'EdgeColor', 'none', 'FaceColor', mblue, 'FaceAlpha', 0.5);

    % Plot force plane
    surf(omega2 / sys.omegaDamped, epsV(i) / epsLim * ones(size(omega2)), rho2, ...
        'EdgeColor', 'none', 'FaceColor', myellow, 'FaceAlpha', 0.5)

    % Plot frequency response
    s = contourslice(omega3 / sys.omegaDamped, eps3 / epsLim, rho3, M3, [], epsV(i) / epsLim, [], [0, 0]);
    for j = 1:length(s)
        s(j).EdgeColor = mred;
        s(j).LineWidth = 2;
    end

    % Plot backbone
    % plot3(omegaBbP / sys.omegaDamped, epsV(i) / epsLim * ones(size(omegaBbP)), rhoBbP, '--', 'Color', 'k', 'LineWidth', 2)

    % Set 2D plot
    nexttile(tl, 2);
    cla;
    xlabel('$\Omega / \omega_0$ [-]')
    ylabel('$\rho$ [-]')
    xlim(minmax(omegaV) / sys.omegaDamped)
    ylim(minmax(rhoV))

    % Plot frequency response
    level2 = (epsV(i) * abs(ssm.S0)).^2;
    contour(omega2 / sys.omegaDamped, rho2, M2, [level2, level2], ...
        'LineWidth', 2, 'EdgeColor', mred, 'FaceColor', 'none')

    % Plot backbone
    plot(omegaBbP / sys.omegaDamped, rhoBbP, '--', 'Color', 'k', 'LineWidth', 2)

    % Save gif
    drawnow
    frame = getframe(figAnimation);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if isFirstFrame == true
        isFirstFrame = false;
        imwrite(imind, cm, fileName, 'gif', 'Loopcount', inf, 'DelayTime', delayTime, 'BackgroundColor', 0);
    else
        imwrite(imind, cm, fileName, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end




















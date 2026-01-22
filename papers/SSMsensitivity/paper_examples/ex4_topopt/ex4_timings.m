% Time the computation of the SSM and its sensitivity
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

%% Timings

% Number of elements
nxVec = [50, 80, 100, 120, 150, 180, 200];
nyVec = nxVec / 5;
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
    rayleigh = [0.01, 0];
    outDofs = myMesh.get_DOF_from_location([lx, ly/2]);
    outDof = myAssembly.free2constrained_index(outDofs(2));
    xTarget = 0.2 * ly;

    % Create system
    sys = MechSystem();
    sys.build_topopt_system(myAssembly, ...
        'computePartialDerivatives', true);
    sys.update_topopt_system(d);

    % Modal analysis
    sys.modal_analysis('nModes', 10, 'modeIndex', 1, ...
        'computePartialDerivatives', true);
    sys.set_damping_rayleigh(rayleigh, 'computePartialDerivatives', true, ...
        'topOptProblem', true)
    
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

% Time complexity
tt = linspace(0, 10000, 101);
% idx = 1;
% Ox = @(t) t(idx, end) / nVec(idx) * tt;
% Oxlogx = @(t) t(idx, end) / (nVec(idx) * log2(nVec(idx))^2) * tt .* log2(tt).^2;
% Oxsqrtx = @(t) t(idx, end) / (nVec(idx).^(3/2)) * tt.^(3/2);
% Ox2 = @(t) t(idx, end) / nVec(idx)^2 * tt.^2;

% Fitting
% nVec = nVec(3:end);
% timeSSM = timeSSM(3:end, :);
% timeSens = timeSens(3:end, :);
[OxSSM, gofOxSSM] = fit(nVec.', timeSSM(:, end), 'poly1');
[Ox, gofOx] = fit(nVec.', timeSens(:, end), 'poly1');
[Oxlogx, gofOxlogx] = fit(nVec.', timeSens(:, end), fittype('a*x*log2(x) + b'));
[Oxsqrtx, gofOxsqrtx] = fit(nVec.', timeSens(:, end), fittype('a*x*sqrt(x) + b'));
[Ox2, gofOx2] = fit(nVec.', timeSens(:, end), 'poly2');
disp(['Ox: ', num2str(gofOx.rmse), ', Oxlogx: ', num2str(gofOxlogx.rmse), ...
      ', Oxsqrtx: ', num2str(gofOxsqrtx.rmse), ', Ox2: ', num2str(gofOx2.rmse)])
disp(['Ox: ', num2str(gofOx.rsquare), ', Oxlogx: ', num2str(gofOxlogx.rsquare), ...
      ', Oxsqrtx: ', num2str(gofOxsqrtx.rsquare), ', Ox2: ', num2str(gofOx2.rsquare)])

% Create figure
% figure('Units','normalized','Position',[.1 .2 .8 .6],'Color','w')
% tiledlayout(1,2,"TileSpacing","compact","Padding","compact")
fig2D = figure;
fig2D.Position(2) = 0.5*fig2D.Position(2);
fig2D.Position(3) = 1300;
fig2D.Position(3:4) = fig2D.Position(3:4)*1.5;
tiledlayout(1, 2, 'TileSpacing', 'tight')

% Font size and colors
fs = 20;
newcolors = parula(7);
colororder(newcolors(2:2:end-1, :))
% m = colormap("parula"); close
% m = m(round(linspace(size(m,1)/10, size(m,1)*9/10, length(maxOrderVec)+2)), :);

% Plot SSM timings
nexttile
hold on; grid on; box on; axis tight;
bSSM = bar(nVec, timeSSM);
xLim = xlim;
yLim = ylim;
plot(tt, OxSSM(tt), 'b--', 'LineWidth', 2)
% plot(tt, Oxlogx(tt), 'g--', 'LineWidth', 2)
% plot(tt, Oxsqrtx(tt), 'r--', 'LineWidth', 2)
% plot(tt, Ox2(tt), 'c--', 'LineWidth', 2)
xlim(xLim)
ylim(yLim)
xlabel('Design variables [-]','Interpreter','latex')
ylabel('Time [s]','Interpreter','latex')
title('Backbone','Interpreter','latex')
% l=legend('O(3)','O(5)','O(7)', 'Interpreter','latex', 'location', 'northeastoutside');
set(gca, 'FontSize', fs, 'TickLabelInterpreter', 'latex')
xticks(nVec)

% Plot sensitivity timings
nexttile
hold on; grid on; box on; axis tight;
bSens = bar(nVec, timeSens);
xLim = xlim;
yLim = ylim;
plot(tt, Ox(tt), 'b--', 'LineWidth', 2)
plot(tt, Oxlogx(tt), 'g--', 'LineWidth', 2)
plot(tt, Oxsqrtx(tt), 'r--', 'LineWidth', 2)
plot(tt, Ox2(tt), 'y--', 'LineWidth', 2)
xlim(xLim)
ylim(yLim)
xlabel('Design variables [-]','Interpreter','latex')
ylabel('Time [s]','Interpreter','latex')
title('Adjoint Sensitivity', 'Interpreter','latex')
% legend('$\mathcal{O}(3)$','$\mathcal{O}(5)$','$\mathcal{O}(7)$', '$O(n)$', '$O(n \log n)$', '$O(n \sqrt{n})$', '$O(n^2)$', 'Interpreter','latex', 'location', 'northeastoutside');
legend('$\mathcal{O}(3)$','$\mathcal{O}(5)$','$\mathcal{O}(7)$', 'Interpreter','latex', 'location', 'northeastoutside');
set(gca, 'FontSize', fs, 'TickLabelInterpreter', 'latex')
xticks(nVec)

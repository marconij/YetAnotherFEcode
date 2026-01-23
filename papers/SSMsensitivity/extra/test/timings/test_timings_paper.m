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

% Number of dofs
nDOFs = 501;

% Output dof
outDof = (nDOFs + 1)/2;

%% Timings

% Number of design variables
muVec = 4 * (1:5);

% Expansion order
maxOrderVec = 3:7;

% Initialize
timeSdd = zeros(length(muVec), length(maxOrderVec));
timeAdj = zeros(length(muVec), length(maxOrderVec));

% Loop over the number of dofs
for i = 1:length(muVec)
    % Number of design variables
    ndv = muVec(i);
    fprintf('ndv = %5d\n', ndv);
    
    % Loop over the expansion order
    for j = 1:length(maxOrderVec)
        % Expansion order
        maxOrder = maxOrderVec(j);
        fprintf('\tmaxOrder = %3d\n', maxOrder);

        % Create system
        mVec = m / (ndv/4) * ones(1, ndv/4);
        kVec = k / (ndv/4) * ones(1, ndv/4);
        k2Vec = k2 / (ndv/4) * ones(1, ndv/4);
        k3Vec = k3 / (ndv/4) * ones(1, ndv/4);
        sys = MechSystem();
        sys.build_test_system(nDOFs, mVec, kVec, k2Vec, k3Vec, ...
            'computePartialDerivatives', true);

        % Modal analysis
        sys.modal_analysis('computePartialDerivatives', true);
        sys.set_damping_rayleigh(rayleigh, 'computePartialDerivatives', true);

        % Direct differentiation sensitivity
        tStart = tic;
        sys.modal_analysis_sensitivity;
        ssm = MechSSM(sys);
        ssm.compute_ssm(maxOrder, 'computeSensitivity', true, ...
            'storeCoefficients', false);
        rhoTarget = ssm.compute_rho_from_x(xTarget, outDof);
        [xTarget, xkTarget] = ssm.compute_x(rhoTarget, outDof);
        OmegaTarget = ssm.compute_omega(rhoTarget);
        drho = ssm.compute_rho_derivative(rhoTarget, outDof, xkTarget);
        dOmegaDD = ssm.compute_omega_derivative(rhoTarget, drho);
        timeSdd(i, j) = toc(tStart);
        fprintf('\t\tTime Sdd = %5.2e\n', timeSdd(i, j));

        % Adjoint sensitivity
        tStart = tic;
        ssm = MechSSM(sys);
        ssm.compute_ssm(maxOrder, 'computeSensitivity', false, ...
            'storeCoefficients', true);
        rhoTarget = ssm.compute_rho_from_x(xTarget, outDof);
        [xTarget, xkTarget] = ssm.compute_x(rhoTarget, outDof);
        OmegaTarget = ssm.compute_omega(rhoTarget);
        sens = SensitivityBackbone(ssm);
        dOmegaAdj = sens.adjoint_sensitivity_Omega(rhoTarget, xkTarget, xTarget, outDof);
        timeAdj(i, j) = toc(tStart);
        fprintf('\t\tTime Adj = %5.2e\n', timeAdj(i, j));

        % Check
        errAbs = norm(dOmegaDD - dOmegaAdj);
        errRel = errAbs / norm(dOmegaDD);
        fprintf('\t\tdOmegaDD = %5.2e, dOmegaAdj = %5.2e, errAbs = %5.2e, errRel = %5.2e\n', ...
            norm(dOmegaDD), norm(dOmegaAdj), errAbs, errRel);
    end
end

%% Post-processing

ySdd = timeSdd(:, end);
rsqSdd = R2(muVec(:), ySdd);

yAdj = timeAdj(:, end);
rsqAdj = R2(muVec(:), yAdj);

%% Plot

% Create grid
[xx, yy] = meshgrid(muVec, maxOrderVec);
xx = xx.';
yy = yy.';

% Create figure
fig = figure;
tiledlayout(2, 1, 'TileSpacing', 'tight')

% Sensitivity direct differentiation
% f_sdd = figure;
nexttile
hold on; grid on; box on; axis tight; view(3);
mesh(xx, yy, timeSdd)
xlabel('Design variables [-]', 'Interpreter', 'latex')
ylabel('Expansion order [-]', 'Interpreter', 'latex')
zlabel('Time [s]', 'Interpreter', 'latex')
title('Sensitivity Direct Differentiation', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Sensitivity adjoint
% f_adj = figure;
nexttile
hold on; grid on; box on; axis tight; view(3);
mesh(xx, yy, timeAdj)
xlabel('Design variables [-]', 'Interpreter', 'latex')
ylabel('Expansion order [-]', 'Interpreter', 'latex')
zlabel('Time [s]', 'Interpreter', 'latex')
title('Sensitivity Adjoint', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Save
% saveas(f_sdd, 'sensitivity_smc_sdd.fig')
% saveas(f_adj, 'sensitivity_smc_adj.fig')
saveas(fig, 'time_sens.fig')

%%
function rsq = R2(x, y)
    yfit = polyval(polyfit(x, y, 1), x);
    yresid = y - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1) * var(y);
    rsq = 1 - SSresid/SStotal;
end

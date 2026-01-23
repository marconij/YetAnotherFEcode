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
alphaR = 0;
betaR = 0.1;
Qfactor = 100;

% Create system
sys = MechSystem();
sys.build_uniform_spring_mass_chain(n, m, k, k2, k3, 'computePartialDerivatives', true);

% Modal analysis
sys.modal_analysis('computePartialDerivatives', true);
sys.set_damping_rayleigh([alphaR, betaR], 'computePartialDerivatives', true)
% sys.set_damping_qfactor(Qfactor, 'computePartialDerivatives', true)

% Modal analysis sensitivity
sys.modal_analysis_sensitivity;

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

%% Validation with SSMTool
if exist('test/testSSM.mat','file') > 0
    % Load
    testSSM = load("test/testSSM.mat");

    % Check the SSM
    fprintf('Check the SSM\n')
    fprintf('R\n')
    for mOrder = 1:2:maxOrder
        err = norm(ssm.R(mOrder).coeffs(1) - testSSM.R0(mOrder).coeffs(1, 2));
        fprintf('\tOrder %d: %.2e\n', mOrder, err)
    end
    fprintf('w\n')
    for mOrder = 1:maxOrder
        err = norm(ssm.w(mOrder).coeffs - testSSM.W0(mOrder).coeffs(1:end/2, end:-1:1));
        fprintf('\tOrder %d: %.2e\n', mOrder, err)
    end
    fprintf('dw\n')
    for mOrder = 1:maxOrder
        err = norm(ssm.dw(mOrder).coeffs - testSSM.W0(mOrder).coeffs(end/2+1:end, end:-1:1));
        fprintf('\tOrder %d: %.2e\n', mOrder, err)
    end
end

%% Plot

% Plot reduced response
figure
hold on; grid on; box on; axis square tight;
plot(Omega, rho, 'k', 'LineWidth', 2)
xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
ylabel('$\rho$ [-]', 'Interpreter', 'latex')
title('Reduced amplitude', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Plot manifold
figure
hold on; grid on; box on; axis square tight; view(3);
surf(rho.' .* cos(ssm.theta), rho.' .* sin(ssm.theta), xk, ...
    'FaceColor', [0.0745 0.6235 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7)
[~, idx] = min(abs(rho-rhoTarget));
plot3(rho(idx) * cos(ssm.theta), rho(idx) * sin(ssm.theta), xk(idx, :), 'r', 'LineWidth', 2);
L = light;
L.Position = [0 0 1];
lighting gouraud
xlabel('$\rho cos(\theta)$', 'Interpreter', 'latex')
ylabel('$\rho sin(\theta)$', 'Interpreter', 'latex')
zlabel('$x$ [m]', 'Interpreter', 'latex')
title('Manifold', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')

% Plot physical response
figure
hold on; grid on; box on; axis square tight;
plot(Omega, x, 'k', 'LineWidth', 2)
xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
ylabel('$x$ [m]', 'Interpreter', 'latex')
title('Physical amplitude', 'Interpreter', 'latex')
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
if exist('test/testSSM.mat','file') > 0
    plot(testSSM.Omega, testSSM.ZoutNorm, '--', 'LineWidth', 2) % SSM validation
end

%% Sensitivity

% Compute target points
[xTarget, xkTarget] = ssm.compute_x(rhoTarget, outDof);
OmegaTarget = ssm.compute_omega(rhoTarget);

% Plot physical response
plot(OmegaTarget, xTarget, '.r', 'MarkerSize', 20)

% Compute the derivative of rho
drho = ssm.compute_rho_derivative(rhoTarget, outDof, xkTarget);

% Compute the derivative of Omega
dOmega = ssm.compute_omega_derivative(rhoTarget, drho);

%% Sensitivity check

% Perturbation step
dMu = 1e-7;

% Collect design variables
mu0 = [m, k, k2, k3];

% Header
fprintf('Analytical, finite difference, absolute error, relative error\n')

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
    sysPert.modal_analysis;
    sysPert.set_damping_rayleigh([alphaR, betaR]);
    % sysPert.set_damping_qfactor(Qfactor);

    % SSM
    ssmPert = MechSSM(sysPert);
    ssmPert.compute_ssm(maxOrder, 'storeCoefficients', true);
    rhoPert = ssmPert.compute_rho_from_x(xTarget, outDof);
    OmegaPert = ssmPert.compute_omega(rhoPert);

    % Modal analysis sensitivity
    fprintf('\n\tOrder 1\n')
    sens = sys.domega(dv);
    sensFD = (sysPert.omega - sys.omega) / dMu;
    test_sensitivity_value_function('omega', sens, sensFD);
    sens = sys.dcsi(dv);
    sensFD = (sysPert.csi - sys.csi) / dMu;
    test_sensitivity_value_function('csi', sens, sensFD);
    sens = sys.domegaDamped(dv);
    sensFD = (sysPert.omegaDamped - sys.omegaDamped) / dMu;
    test_sensitivity_value_function('omegaD', sens, sensFD);
    sens = sys.dphi(dv).coeffs;
    sensFD = (sysPert.phi - sys.phi) / dMu;
    test_sensitivity_value_function('phi', sens, sensFD);
    sens = sys.dlambda(dv);
    sensFD = (sysPert.lambda - sys.lambda) / dMu;
    test_sensitivity_value_function('lam', sens, sensFD);

    % Loop over orders
    for mOrder = 2:ssm.maxOrder
        fprintf('\n\tOrder %d\n', mOrder)

        % Loop over multi-indices
        for mIdx = 1:ssm.MIs(mOrder).n
            % Current multi-index
            MI = ssm.MIs(mOrder).coeffs(:, mIdx);
            % fprintf('\nMI %d%d\n', MI(1), MI(2))

            % Check if the auxiliary coefficients have been stored
            if ssm.storeCoefficients && ssmPert.storeCoefficients
                % Sensitivity of LambdaM
                LambdaMName = ['LaM', num2str(MI(1)), num2str(MI(2))];
                LambdaMm0 = ssm.LambdaM(mOrder).coeffs(mIdx);
                LambdaMm = ssmPert.LambdaM(mOrder).coeffs(mIdx);
                dLambdaMm = ssm.LambdaM(mOrder).dv(dv).coeffs(mIdx);
                dLambdaMmFD = (LambdaMm - LambdaMm0) / dMu;
                test_sensitivity_value_function(LambdaMName, dLambdaMm, dLambdaMmFD);

                % Sensitivity of f
                fName = ['f_', num2str(MI(1)), num2str(MI(2))];
                fm0 = ssm.f(mOrder).coeffs(:, mIdx);
                fm = ssmPert.f(mOrder).coeffs(:, mIdx);
                dfm = ssm.f(mOrder).dv(dv).coeffs(:, mIdx);
                dfmFD = (fm - fm0) / dMu;
                test_sensitivity_value_function(fName, dfm, dfmFD);

                % Sensitivity of V
                VName = ['V_', num2str(MI(1)), num2str(MI(2))];
                Vm0 = ssm.V(mOrder).coeffs(:, mIdx);
                Vm = ssmPert.V(mOrder).coeffs(:, mIdx);
                dVm = ssm.V(mOrder).dv(dv).coeffs(:, mIdx);
                dVmFD = (Vm - Vm0) / dMu;
                test_sensitivity_value_function(VName, dVm, dVmFD);

                % Sensitivity of Y
                YName = ['Y_', num2str(MI(1)), num2str(MI(2))];
                Ym0 = ssm.Y(mOrder).coeffs(:, mIdx);
                Ym = ssmPert.Y(mOrder).coeffs(:, mIdx);
                dYm = ssm.Y(mOrder).dv(dv).coeffs(:, mIdx);
                dYmFD = (Ym - Ym0) / dMu;
                test_sensitivity_value_function(YName, dYm, dYmFD);

                % Sensitivity of C
                CName = ['C_', num2str(MI(1)), num2str(MI(2))];
                Cm0 = ssm.C(mOrder).coeffs(:, mIdx);
                Cm = ssmPert.C(mOrder).coeffs(:, mIdx);
                dCm = ssm.C(mOrder).dv(dv).coeffs(:, mIdx);
                dCmFD = (Cm - Cm0) / dMu;
                test_sensitivity_value_function(CName, dCm, dCmFD);

                % Sensitivity of h
                hName = ['h_', num2str(MI(1)), num2str(MI(2))];
                hm0 = ssm.h(mOrder).coeffs(mIdx);
                hm = ssmPert.h(mOrder).coeffs(mIdx);
                dhm = ssm.h(mOrder).dv(dv).coeffs(mIdx);
                dhmFD = (hm - hm0) / dMu;
                test_sensitivity_value_function(hName, dhm, dhmFD);
            end

            % Sensitivity of R only if |MI(1) - MI(2)| = 1
            if abs(MI(1) - MI(2)) == 1
                if MI(1) - MI(2) == 1
                    RIdx = 1;
                else
                    RIdx = 2;
                end
                RName = ['R_', num2str(MI(1)), num2str(MI(2))];
                Rm0 = ssm.R(mOrder).coeffs(RIdx);
                Rm = ssmPert.R(mOrder).coeffs(RIdx);
                dRm = ssm.R(mOrder).dv(dv).coeffs(RIdx);
                dRmFD = (Rm - Rm0) / dMu;
                test_sensitivity_value_function(RName, dRm, dRmFD);
            end

            % Sensitivity of w
            wName = ['w_', num2str(MI(1)), num2str(MI(2))];
            wm0 = ssm.w(mOrder).coeffs(:, mIdx);
            wm = ssmPert.w(mOrder).coeffs(:, mIdx);
            dwm = ssm.w(mOrder).dv(dv).coeffs(:, mIdx);
            dwmFD = (wm - wm0) / dMu;
            test_sensitivity_value_function(wName, dwm, dwmFD);

            % Sensitivity of dw
            dwName = ['dw_', num2str(MI(1)), num2str(MI(2))];
            dwm0 = ssm.dw(mOrder).coeffs(:, mIdx);
            dwm = ssmPert.dw(mOrder).coeffs(:, mIdx);
            ddwm = ssm.dw(mOrder).dv(dv).coeffs(:, mIdx);
            ddwmFD = (dwm - dwm0) / dMu;
            test_sensitivity_value_function(dwName, ddwm, ddwmFD);
        end
    end

    % Sensitivity of the backbone (rho and Omega)
    fprintf('\n\tBackbone\n')
    sens = drho(dv);
    sensFD = (rhoPert - rhoTarget) / dMu;
    test_sensitivity_value_function('rho', sens, sensFD);
    sens = dOmega(dv);
    sensFD = (OmegaPert - OmegaTarget) / dMu;
    test_sensitivity_value_function('Omega', sens, sensFD);
end

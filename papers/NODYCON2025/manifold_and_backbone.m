function [Omega, theta, wManifold, wBB_rms, wBB_infNorm] = manifold_and_backbone(M, K, T2, T3, omega, phi, rho, idx, doPlot)
    % Compute the 3rd order manifold and the relative backbone curve.
    % While gamma_value and gamma_sensitivity carry out the LSM computation up to
    % the 2nd order until gamma is obtained, this function computes the also 3rd 
    % order contributions to compute the manifold w to the 3rd order.
    % Inputs:
    %   M: mass matrix
    %   K: stiffness matrix
    %   T2: quadratic stiffness tensor
    %   T3: cubic stiffness tensor
    %   omega: eigenfrequency
    %   phi: mode shape
    %   rhoBB: amplitude in the reduced space
    %   idx: index of the output degrees of freedom
    %   doPlot: flag to plot the results (optional, default is false)
    % Outputs:
    %   Omega: backbone curve frequency
    %   theta: polar coordinates for the manifold
    %   wManifold: manifold in the physical space
    %   wBB_rms: backbone curve in the physical space (RMS value for the output dof)
    %   wBB_infNorm: backbone curve in the physical space (inf-norm value for the output dof)

    % Check inputs
    if nargin < 9
        doPlot = false;
    end

    % Order 2: index 20
    Lambda20 = 2i * omega;
    f20 = double(ttv(T2, {phi, phi}, [3, 2]));
    L20 = K + Lambda20^2 * M;
    w20 = lsqminnorm(L20, -f20);
    
    % Order 2: index 11
    f11 = 2 * f20;
    L11 = K;
    w11 = lsqminnorm(L11, -f11);
    
    % Order 2: index 02
    w02 = conj(w20);

    % Order 3: index 30
    Lambda30 = 3i * omega;
    f30 = double(ttv(T3, {phi, phi, phi}, [4, 3, 2]) + ...
                 ttv(T2, {w20, phi}, [3, 2]) + ...
                 ttv(T2, {phi, w02}, [3, 2]));
    L30 = K + Lambda30^2 * M;
    w30 = lsqminnorm(L30, -f30);

    % Order 3: index 21
    Lambda21 = 1i * omega;
    f21 = double(ttv(T2, {w20 + w11, phi}, [3, 2]) + ...
                 ttv(T2, {w20 + w11, phi}, [2, 3]) + ...
             3 * ttv(T3, {phi, phi, phi}, [4, 3, 2]));
    L21 = K + Lambda21^2 * M;
    w21 = lsqminnorm(L21, -f21);
    
    % Order 2: index 02
    w12 = conj(w21);
    
    % Order 2: index 03
    w03 = conj(w30);
    
    % Compute gamma
    gamma = 1/(2*omega) * phi.' * f21;

    % Backbone in the reduced space
    Omega = omega + gamma * rho.^2;

    % Theta points
    nTheta = 2^7;
    theta = linspace(0, 2*pi, nTheta);

    % Compute manifold
    wManifold = zeros(length(rho), length(theta));
    for ii = 1:length(rho)
        for kk = 1:length(theta)
            p = rho(ii) * exp([1i; -1i] * theta(kk));
            wManifold(ii, kk) = phi(idx) * p(1) + phi(idx) * p(2) + ...
                w20(idx) * p(1)^2 + w11(idx) * p(1) * p(2) + w02(idx) * p(2)^2 + ...
                w30(idx) * p(1)^3 + w21(idx) * p(1)^2 * p(2) + ...
                w12(idx) * p(1) * p(2)^2 + w03(idx) * p(2)^3;
        end
    end
    wManifold = real(wManifold);
    
    % inf-norm of the backbone (used for comparison with SSMtool)
    for ii = 1:length(rho)
        wBB_infNorm(ii) = norm(wManifold(ii,:),'inf');
    end

    % Compute backbone in the physical space
    wBB_rms = sqrt(sum(wManifold.^2, 2) / nTheta);

    % Plot the results
    if doPlot
        fprintf(' dof: %d \n', idx)
        fprintf(' gamma = %g\n', gamma)

        % Plot manifold
        figure
        hold on; grid on; box on;
        surf(rho.' * cos(theta), rho.' * sin(theta), wManifold, ...
            'FaceColor', [0.0745 0.6235 1], 'EdgeColor', 'none', 'FaceAlpha', 0.7)
        L = light;
        L.Position = [0 0 1];
        lighting gouraud
        axis square tight
        xlabel('$\rho \cos(\theta)$', 'Interpreter', 'latex')
        ylabel('$\rho \sin(\theta)$', 'Interpreter', 'latex')
        zlabel(sprintf('$w_{%d}$',idx), 'Interpreter', 'latex')
        title('Lyapunov Subcenter Manifold','Interpreter','latex')
        set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex')
        view(45,15)

        % Plot backbone
        figure
        plot(Omega/2/pi/1000, wBB_rms, 'LineWidth', 2)
        xlabel('$f$ [kHz]', 'Interpreter', 'latex')
        ylabel('RMS($w$)', 'Interpreter', 'latex')
        set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex')
        axis square tight
        grid on

        % Plot backbone (compare vs SSMtool output)
        figure
        plot(Omega, wBB_infNorm, 'LineWidth', 2)
        xlabel('$\Omega$', 'Interpreter', 'latex')
        ylabel(sprintf('$||w_{%d}||_\\infty$',idx), 'Interpreter', 'latex')
        set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex')
        axis square tight
        grid on
    end
end

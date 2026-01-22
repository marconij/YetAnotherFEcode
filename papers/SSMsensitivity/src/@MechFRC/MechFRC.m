classdef MechFRC < MechSSM
    % A class for computing FRCs associated with 2D autonomous SSMs for
    % lightly-damped mechanical systems.
    % Properties:
    %   Omega0: the external frequency.
    %   dOmega0: the sensitivity of the external frequency.
    %   S0: the non-autonomous reduced dynamics coefficient.
    %   dS0: the sensitivity of the non-autonomous reduced dynamics
    %       coefficient.
    %   Lambda0: the non-autonomous eigenvalue coefficient.
    %   dLambda0: the sensitivity of the non-autonomous eigenvalue
    %       coefficient.
    % Methods:
    %   MechFRC: the constructor.
    %   compute_S: compute the non-autonomous reduced dynamics coefficient.
    %   compute_a: compute the coefficient a.
    %   compute_b: compute the coefficient b.
    %   compute_phase: compute the phase angle psi.
    %   state_equation_fr_peak: compute the state equation for the FR peak.
    %   compute_frequency: compute the response frequency omega.
    %   compute_amplitude: compute the physical amplitude x.
    %   compute_frc_peak: compute the FR peak given the amplitude of the
    %       external force.
    %   compute_frc_peak_derivative: compute the sensitivity of the FR
    %       peak given the amplitude of the external force.
    %   compute_frc: compute the frequency response curve given the peak
    %       reduced amplitude. The backbone curve is also computed.
    %   compute_jacobian: compute the Jacobian of the frequency response
    %       curve. The stability of the periodic response is determined by
    %       the eigenvalues of the Jacobian matrix.
    % Methods (static):
    %   wrapPhase: wrap the phase angle x to the interval [-pi, 0].

    properties
        inDof
        Omega0, dOmega0
        Lambda0, dLambda0
        S0, dS0
    end

    methods
        function obj = MechFRC(sys, varargin)
            % Construct an instance of this class.
            % Inputs:
            %   sys: the MechSystem object.
            %   nTheta: number of time instants (optional, default is 2^7).
            %   rmsStyle: flag for the RMS computation style (optional,
            %       default is 'SSMTool').
            % Outputs:
            %   obj: the MechFRC object.

            % Call the constructor of the superclass
            obj = obj@MechSSM(sys, varargin{:});

            % Set the external frequency
            obj.Omega0 = obj.sys.omegaDamped;
            obj.dOmega0 = obj.sys.domegaDamped;
        end

        function compute_S(obj, inDof)
            % Compute the non-autonomous reduced dynamics coefficient.
            % Inputs:
            %   inDof: the input DOF where the load is applied.

            % Non-autonomous eigenvalue coefficient
            obj.Lambda0 = 1i * obj.Omega0;
            if obj.computeSensitivity
                obj.dLambda0 = zeros(obj.sys.ndv, 1);
                for dv = 1:obj.sys.ndv
                    obj.dLambda0(dv) = 1i * obj.dOmega0(dv);
                end
            else
                obj.dLambda0 = [];
            end

            % Compute the non-autonomous reduced dynamics coefficient
            den = obj.sys.delta + obj.sys.lambda + obj.Lambda0;
            obj.S0 = -1/2 * obj.sys.phi(inDof) / den;

            % Compute the derivative of S
            if obj.computeSensitivity
                obj.dS0 = zeros(obj.sys.ndv, 1);
                for dv = 1:obj.sys.ndv
                    dnum = -1/2 * obj.sys.dphi(dv).coeffs(inDof);
                    dden = obj.sys.ddelta(dv) + obj.sys.dlambda(dv) + obj.dLambda0(dv);
                    obj.dS0(dv) = (dnum - obj.S0 * dden) / den;
                end
            else
                obj.dS0 = [];
            end

            % Store the input DOF
            obj.inDof = inDof;
        end

        function a = compute_a(obj, rho)
            % Compute the coefficient a.
            % Inputs:
            %   rho: the reduced amplitude.
            % Outputs:
            %   a: the coefficient a.

            a = 0;
            for mOrder = 1:length(obj.R)
                a = a + real(obj.R(mOrder).coeffs(1)) * rho.^mOrder;
            end
        end

        function b = compute_b(obj, rho)
            % Compute the coefficient b.
            % Inputs:
            %   rho: the reduced amplitude.
            % Outputs:
            %   b: the coefficient b.

            b = 0;
            for mOrder = 1:length(obj.R)
                b = b + imag(obj.R(mOrder).coeffs(1)) * rho.^mOrder;
            end
        end

        function psi = compute_phase(obj, rho, omega, epsilon)
            % Compute the phase angle psi.
            % Inputs:
            %   rho: the reduced amplitude.
            %   omega: the response frequency.
            %   epsilon: the amplitude of the external force.
            % Outputs:
            %   psi: the phase angle psi.

            a = obj.compute_a(rho);
            b = obj.compute_b(rho) - real(omega) .* rho;
            f = epsilon * obj.S0;

            psi = atan2((b * real(f) - a * imag(f)), (-a * real(f) - b * imag(f))); % [-pi; pi]
            psi = obj.wrapPhase(psi); % wrap to [-2*pi; 0]
        end

        function r = state_equation_fr_peak(obj, rho, epsilon)
            % Compute the state equation for the FR peak.
            % Inputs:
            %   rho: the reduced amplitude.
            %   epsilon: the amplitude of the external force.
            % Outputs:
            %   r: the residual.
        
            a = obj.compute_a(rho);
            f = epsilon * obj.S0;
            r = abs(f).^2 - a.^2;
        end

        function omega = compute_frequency(obj, rho)
            % Compute the response frequency omega.
            % Inputs:
            %   rho: the reduced amplitude.
            % Outputs:
            %   omega: the response frequency.
            omega = obj.compute_b(rho) ./ rho;
        end

        function [x, xi] = compute_amplitude(obj, rho, psi, outDof)
            % Compute the physical amplitude x.
            % Inputs:
            %   rho: the reduced amplitude.
            %   psi: the phase angle psi.
            %   outDof: the output DOF.
            % Outputs:
            %   x: the physical amplitude.
            %   xi: the physical amplitude at all time instants.

            % Loop over time instants
            xi = zeros(1, obj.nTheta);
            for i = 1:obj.nTheta
                % Adjust the phase angle
                theta = obj.theta(i) + pi/2 + psi; % the +pi/2 is needed for compatibility with the SSMTool (ask Matteo)

                % Compute the reduced polar coordinates
                p = rho * exp([1i * theta; -1i * theta]);

                % Loop over orders
                for mOrder = 1:obj.maxOrder
                    % Loop over multi-indices
                    for mIdx = 1:obj.MIs(mOrder).n
                        % Current multi-index
                        m = obj.MIs(mOrder).coeffs(:, mIdx);

                        % Add the contribution
                        xi(i) = xi(i) + obj.w(mOrder).coeffs(outDof, mIdx) * p(1).^m(1) * p(2).^m(2);
                    end
                end
            end

            % Make sure the response is a real number
            xi = real(xi);

            % Compute the RMS of the physical amplitude
            x = sqrt(sum(xi.^2) / obj.nTimeSteps);
        end

        function [omegaPk, xPk, rhoPk, xiPk] = compute_frc_peak(obj, outDof, epsilon, varargin)
            % Compute the FR peak given the amplitude of the external
            % force.
            % Inputs:
            %   outDof: the output DOF.
            %   epsilon: the amplitude of the external force.
            %   rho0: the initial reduced amplitude (optional, default is
            %       0.1).
            % Outputs:
            %   omegaPk: the peak frequency.
            %   xPk: the peak amplitude.
            %   rhoPk: the peak reduced amplitude.
            %   xiPk: the peak amplitude at all time instants.

            % Parse the inputs
            p = inputParser;
            addOptional(p, 'rho0', 0.1, @isnumeric);
            parse(p, varargin{:});
            rho0 = p.Results.rho0;

            % Find peak rho
            [rhoPk, ~, exitflag] = fzero(@(x) obj.state_equation_fr_peak(x, epsilon), rho0);
            if exitflag <= 0
                fprintf('Equation not solved (%d) for epsilon %.4f\n', exitflag, epsilon)
            end

            % Check if the solution is negative
            if rhoPk < 0
                % Try to find a solution in the positive domain
                rho0 = -rhoPk;
                [rhoPk, ~, exitflag] = fzero(@(x) obj.state_equation_fr_peak(x, epsilon), rho0);
                if exitflag <= 0
                    fprintf('Equation not solved (%d) for epsilon %.4f\n', exitflag, epsilon)
                end
            end

            % The phase shift at resonance is -pi/2
            psiPk = -pi/2;

            % Compute peak frequency
            omegaPk = obj.compute_frequency(rhoPk);

            % Compute peak amplitude
            [xPk, xiPk] = obj.compute_amplitude(rhoPk, psiPk, outDof);
        end

        function [domegaPk, dxPk, drhoPk] = compute_frc_peak_derivative(obj, outDof, epsilon, rhoPk, xPk, xiPk)
            % Compute the sensitivity of the FR peak.

            % Sensitivity of a
            a = obj.compute_a(rhoPk);
            da = epsilon^2 * (conj(obj.S0) * obj.dS0 + conj(obj.dS0) * obj.S0) / (2*a);

            % Sensitivity of reduced amplitude
            num = da;
            den = 0;
            for mOrder = 1:2:obj.maxOrder
                den = den + real(obj.R(mOrder).coeffs(1)) * mOrder * rhoPk^(mOrder-1);
                for dv = 1:obj.sys.ndv
                    num(dv) = num(dv) - real(obj.R(mOrder).dv(dv).coeffs(1)) * rhoPk^mOrder;
                end
            end
            drhoPk = num / den;

            % Sensitivity of frequency
            domegaPk = imag(obj.sys.dlambda);
            pomega_rho = 0;
            for mOrder = 3:2:obj.maxOrder
                pomega_rho = pomega_rho + (mOrder - 1) * imag(obj.R(mOrder).coeffs(1)) * rhoPk^(mOrder - 2);
                for dv = 1:obj.sys.ndv
                    domegaPk(dv) = domegaPk(dv) + imag(obj.R(mOrder).dv(dv).coeffs(1)) * rhoPk^(mOrder - 1);
                end
            end
            domegaPk = domegaPk + drhoPk * pomega_rho;

            % Sensitivity of physical amplitude
            psiPk = -pi/2;
            pxi_rho = zeros(1, obj.nTheta);
            dxPk = 0;
            for k = 1:obj.nTheta
                % Adjust the phase angle
                theta = obj.theta(k) + pi/2 + psiPk; % the +pi/2 is needed for compatibility with the SSMTool (ask Matteo)
    
                % Compute the reduced polar coordinates
                pTilde = exp([1i * theta; -1i * theta]);
                p = rhoPk * pTilde;
    
                % Loop over orders
                dxikPk = zeros(obj.sys.ndv, 1);
                for mOrder = 1:obj.maxOrder
                    % Loop over multi-indices
                    for mIdx = 1:obj.MIs(mOrder).n
                        % Current multi-index
                        m = obj.MIs(mOrder).coeffs(:, mIdx);
    
                        % Add contribution
                        pxi_rho(k) = pxi_rho(k) + obj.w(mOrder).coeffs(outDof, mIdx) * prod(pTilde.^m) * mOrder * rhoPk^(mOrder - 1);
                        for dv = 1:obj.sys.ndv
                            dxikPk(dv) = dxikPk(dv) + obj.w(mOrder).dv(dv).coeffs(outDof, mIdx) * prod(p.^m);
                        end
                    end
                end
    
                % Assemble derivative
                dxikPk = dxikPk + drhoPk * pxi_rho(k);
                dxPk = dxPk + xiPk(k) * dxikPk;
            end
            dxPk = real(dxPk) / (obj.nTimeSteps * xPk);
        end

        function [omegaFr, xFr, psiFr, omegaBb, xBb, rhoBb] = compute_frc(obj, rhoPk, outDof, epsilon, varargin)
            % Compute the frequency response curve given the peak reduced
            % amplitude. The backbone curve is also computed.
            % Inputs:
            %   rho: the peak reduced amplitude.
            %   outDof: the output DOF.
            %   epsilon: the amplitude of the external force.
            %   nRho: the number of reduced amplitudes (optional, default
            %       is 51).
            % Outputs:
            %   omegaFr: the response frequency.
            %   xFr: the physical amplitude.
            %   psiFr: the phase angle.
            %   omegaBb: the backbone frequency.
            %   xBb: the backbone amplitude.
            %   rhoBb: the backbone reduced amplitude.

            % Parse the inputs
            p = inputParser;
            addOptional(p, 'nRho', 51);
            parse(p, varargin{:});
            nRho = p.Results.nRho;

            % Create vector of reduced amplitudes
            rhoBb = linspace(0, rhoPk, nRho + 1);
            rhoBb(1) = []; % remove the zero amplitude

            % Initialize the output vectors
            omegaBb = zeros(1, nRho);
            xBb = zeros(1, nRho);
            omegaFr = zeros(1, 2*nRho);
            xFr = zeros(1, 2*nRho);
            psiFr = zeros(1, 2*nRho);

            % Loop over reduced amplitudes
            for i = 1:nRho
                % Compute backbone
                omegaBb(i) = obj.compute_frequency(rhoBb(i));
                xBb(i) = obj.compute_amplitude(rhoBb(i), -pi/2, outDof);

                % Compute response frequency
                a = obj.compute_a(rhoBb(i));
                f = epsilon * obj.S0;
                delta = abs(f).^2 - a.^2;
                if abs(delta) < eps
                    delta = 0; % make sure to avoid negative zero
                end
                omegaFr(i) = omegaBb(i) - sqrt(delta) ./ rhoBb(i);
                omegaFr(end - i + 1) = omegaBb(i) + sqrt(delta) ./ rhoBb(i);

                % Compute phase angle
                psiFr(i) = obj.compute_phase(rhoBb(i), omegaFr(i), epsilon);
                psiFr(end - i + 1) = obj.compute_phase(rhoBb(i), omegaFr(end - i + 1), epsilon);

                % Compute physical amplitude
                xFr(i) = obj.compute_amplitude(rhoBb(i), psiFr(i), outDof);
                xFr(end - i + 1) = obj.compute_amplitude(rhoBb(i), psiFr(end - i + 1), outDof);
            end
        end

        function [J, isStable] = compute_jacobian(obj, rho, Omega)
            % Compute the Jacobian of the frequency response curve.
            % The stability of the periodic response is determined by the 
            % eigenvalues of the Jacobian matrix. In particular, the
            % periodic response is stable if tr(J) < 0 and det(J) > 0.
            % Inputs:
            %   rho: the reduced amplitude.
            %   Omega: the response frequency.
            % Outputs:
            %   J: the Jacobian matrix.
            %   isStable: a flag indicating whether the periodic response
            %       is stable (1) or not (0).

            % Coefficients a and b
            a = 0;
            b = 0;
            for mOrder = 1:length(obj.R)
                a = a + real(obj.R(mOrder).coeffs(1)) * rho.^mOrder;
                b = b + imag(obj.R(mOrder).coeffs(1)) * rho.^mOrder;
            end

            % Partial derivatives
            pa_rho = real(obj.R(1).coeffs(1));
            pb_rho = imag(obj.R(1).coeffs(1));
            for mOrder = 2:length(obj.R)
                pa_rho = pa_rho + real(obj.R(mOrder).coeffs(1)) * mOrder * rho.^(mOrder - 1);
                pb_rho = pb_rho + imag(obj.R(mOrder).coeffs(1)) * mOrder * rho.^(mOrder - 1);
            end

            % Jacobian matrix
            J = [pa_rho, rho * Omega - b;
                (pb_rho - Omega) / rho, a / rho];
            
            % Check stability
            % E = eig(J);
            % isStable = all(real(E) < 0);
            isStable = trace(J) < 0 && det(J) > 0;
        end

        function fig = plot_frequency_response(obj, epsV, omegaV, rhoV, varargin)
            % Plot the frequency response curve.
            % Inputs:
            %   epsV: the vector of amplitudes of the external force.
            %   omegaV: the vector of frequency values.
            %   rhoV: the vector of amplitude values.
            %   nPts: the number of points for the meshgrid (optional,
            %       default is 101).
            %   nMaxEps: the maximum number of epsilons to be plotted
            %       (optional, default is 7).
            %   plotPhase: flag to also plot the phase (optional, default is
            %       false).
            % Outputs:
            %   fig: the figure handle.

            % Parse inputs
            p = inputParser;
            addOptional(p, 'nPts', 101, @isnumeric);
            addOptional(p, 'nMaxEps', 7, @isnumeric);
            addOptional(p, 'plotPhase', false, @islogical);
            parse(p, varargin{:});
            nPts = p.Results.nPts;
            nMaxEps = p.Results.nMaxEps;
            plotPhase = p.Results.plotPhase;

            % Check if the number of epsilons is too large
            if length(epsV) > nMaxEps
                warning('Too many epsilons, only the first %d will be plotted', nMaxEps);
                epsV = epsV(1:nMaxEps);
            end

            % Check if input vectors are valid
            if length(omegaV) == 2
                omegaLim = omegaV;
                omegaV = linspace(omegaLim(1), omegaLim(2), nPts);
            end
            if length(rhoV) == 2
                rhoLim = rhoV;
                rhoV = linspace(rhoLim(1), rhoLim(2), nPts + 1);
                rhoV(1) = []; % remove zero amplitude
            end
            if plotPhase
                psiV = linspace(-2*pi, 0, nPts);
            end

            % Compute backbone
            omegaBb = obj.compute_b(rhoV) ./ rhoV;
            psiBb = -pi/2 * ones(size(rhoV));
            rhoBb = rhoV;

            % Compute frequency response level set function
            [omega2, rho2] = meshgrid(omegaV, rhoV);
            ls = obj.compute_a(rho2).^2 + (obj.compute_b(rho2) - omega2.*rho2).^2;

            % Initialize figure
            fig = figure;
            tiledlayout('flow', 'TileSpacing', 'compact')
            nexttile
            hold on; grid on; box on; axis tight;

            % Plot frequency response curves for different epsilons
            colorLines = lines(length(epsV));
            for i = 1:length(epsV)
                contour(omega2 / obj.sys.omegaDamped, rho2, ls, (epsV(i) * abs(obj.S0))^2 * [1, 1], ...
                    'LineColor', colorLines(i, :), 'LineWidth', 2, 'DisplayName', num2str(epsV(i), '$\\epsilon = %.2e$'));
            end

            % Plot backbone
            plot(omegaBb / obj.sys.omegaDamped, rhoBb, '--', 'Color', 'k', 'LineWidth', 2, 'DisplayName', 'Backbone');

            % Decorate plot
            xlabel('$\Omega / \omega_0$ [-]')
            ylabel('$\rho$ [-]')
            xlim(minmax(omegaV) / obj.sys.omegaDamped)
            ylim(minmax(rhoV))
            legend('Location', 'best');

            % Check if phase plot is requested
            if plotPhase
                % Initialize phase plot
                nexttile
                hold on; grid on; box on; axis tight;

                % Plot phase curves for different epsilons
                for i = 1:length(epsV)
                    % Compute boundaries from contour matrix
                    c = contourc(omegaV, rhoV, ls, (epsV(i) * abs(obj.S0))^2 * [1, 1]);
                    boundaries = obj.extract_boundaries(c);

                    % Loop over boundaries
                    for j = 1:length(boundaries)
                        % Extract boundary points
                        omegaFr = boundaries{j}(1, :);
                        rhoFr = boundaries{j}(2, :);

                        % Compute phase
                        psiFr = obj.compute_phase(rhoFr, omegaFr, epsV(i));

                        % Plot phase
                        plot(omegaFr / obj.sys.omegaDamped, psiFr / pi, ...
                            'Color', colorLines(i, :), 'LineWidth', 2);
                    end
                end

                % Plot -pi/2 line
                plot(omegaV / obj.sys.omegaDamped, psiBb / pi, '--', 'Color', 'k', 'LineWidth', 2);

                % Decorate plot
                xlabel('$\Omega / \omega_0$ [-]')
                ylabel('$\psi / \pi$ [-]')
                xlim(minmax(omegaV) / obj.sys.omegaDamped)
                ylim(minmax(psiV) / pi)
            end
        end

        function fig = plot_frequency_phase_amplitude(obj, epsT, omegaV, psiV, rhoV, varargin)
            % Plot the zero level sets of the state equations in the
            % frequency-phase-amplitude space.
            % Inputs:
            %   epsT: the amplitude of the external force.
            %   omegaV: the vector of frequency values.
            %   psiV: the vector of phase values.
            %   rhoV: the vector of amplitude values.
            %   nPts: the number of points for the meshgrid (optional,
            %       default is 101).
            %   frcData: the frequency response data to be overlaid on the
            %       3D plot (optional, default is empty, in which case the
            %       data is computed).
            %   plot2D: flag to also plot the 2D projections (optional,
            %       default is false).
            % Outputs:
            %   fig: the figure handle.

            % Parse inputs
            p = inputParser;
            addOptional(p, 'nPts', 101, @isnumeric);
            addOptional(p, 'frcData', [], @iscell);
            addOptional(p, 'plot2D', false, @islogical);
            parse(p, varargin{:});
            nPts = p.Results.nPts;
            frcData = p.Results.frcData;
            plot2D = p.Results.plot2D;

            % Check if input vectors are valid
            if length(omegaV) == 2
                omegaLim = omegaV;
                omegaV = linspace(omegaLim(1), omegaLim(2), nPts);
            end
            if length(psiV) == 2
                psiLim = psiV;
                psiV = linspace(psiLim(1), psiLim(2), nPts);
            end
            if length(rhoV) == 2
                rhoLim = rhoV;
                rhoV = linspace(rhoLim(1), rhoLim(2), nPts + 1);
                rhoV(1) = []; % remove zero amplitude
            end

            % Compute backbone
            omegaBb = obj.compute_b(rhoV) ./ rhoV;
            psiBb = -pi/2 * ones(size(rhoV));
            rhoBb = rhoV;

            % Compute frequency response data if not provided
            if isempty(frcData)
                % 2D grid
                [omega2, rho2] = meshgrid(omegaV, rhoV);

                % Level set function
                ls = obj.compute_a(rho2).^2 + (obj.compute_b(rho2) - omega2.*rho2).^2 - (epsT * abs(obj.S0))^2;

                % Compute boundaries from contour matrix
                c = contourc(omegaV, rhoV, ls, [0, 0]);
                boundaries = obj.extract_boundaries(c);

                % Create frcData
                frcData = cell(size(boundaries));
                for j = 1:length(boundaries)
                    % Extract boundary points
                    omegaFr = boundaries{j}(1, :);
                    rhoFr = boundaries{j}(2, :);

                    % Compute phase
                    psiFr = obj.compute_phase(rhoFr, omegaFr, epsT);

                    % Store data
                    frcData{j} = [omegaFr; psiFr; rhoFr]';
                end
            end

            % 3D grid
            [omegaM, psiM, rhoM] = meshgrid(omegaV, psiV, rhoV);

            % Define functions
            M1 = obj.compute_a(rhoM) + epsT * (real(obj.S0) * cos(psiM) + imag(obj.S0) * sin(psiM));
            M2 = obj.compute_b(rhoM) - omegaM.*rhoM + epsT * (imag(obj.S0) * cos(psiM) - real(obj.S0) * sin(psiM));

            % Initialize figure
            fig = figure;
            tiledlayout('flow', 'TileSpacing', 'compact')

            % Frequency-phase-amplitude plot
            nexttile
            hold on; grid on; box on; axis square tight;
            patch(isosurface(omegaM / obj.sys.omegaDamped, psiM / pi, rhoM, M1, 0), ...
                'EdgeColor', 'none', 'FaceColor', '#0072BD', 'FaceAlpha', 0.5);
            patch(isosurface(omegaM / obj.sys.omegaDamped, psiM / pi, rhoM, M2, 0), ...
                'EdgeColor', 'none', 'FaceColor', '#EDB120', 'FaceAlpha', 0.5);
            for j = 1:length(frcData)
                omegaFr = frcData{j}(:, 1);
                psiFr = frcData{j}(:, 2);
                rhoFr = frcData{j}(:, 3);
                plot3(omegaFr / obj.sys.omegaDamped, psiFr / pi, rhoFr, '-', 'Color', '#D95319', 'LineWidth', 2);
            end
            plot3(omegaBb / obj.sys.omegaDamped, psiBb / pi, rhoBb, '--', 'Color', 'k', 'LineWidth', 2);
            xlabel('$\Omega / \omega_0$ [-]')
            ylabel('$\psi / \pi$ [-]')
            zlabel('$\rho$ [-]')
            legend(['$\mathcal{M}_1$', '$\mathcal{M}_2$', repmat({''}, 1, length(frcData) + 1)]);
            xlim(minmax(omegaV) / obj.sys.omegaDamped)
            ylim(minmax(psiV) / pi)
            zlim(minmax(rhoV))
            view(3)

            % Check if 2D plots are requested
            if plot2D
                % Frequency-amplitude
                nexttile
                hold on; grid on; box on; axis square tight;
                patch(isosurface(omegaM / obj.sys.omegaDamped, psiM / pi, rhoM, M1, 0), ...
                    'EdgeColor', 'none', 'FaceColor', '#0072BD', 'FaceAlpha', 0.5);
                patch(isosurface(omegaM / obj.sys.omegaDamped, psiM / pi, rhoM, M2, 0), ...
                    'EdgeColor', 'none', 'FaceColor', '#EDB120', 'FaceAlpha', 0.5);
                for j = 1:length(frcData)
                    omegaFr = frcData{j}(:, 1);
                    psiFr = frcData{j}(:, 2);
                    rhoFr = frcData{j}(:, 3);
                    plot3(omegaFr / obj.sys.omegaDamped, psiFr / pi, rhoFr, '-', 'Color', '#D95319', 'LineWidth', 2);
                end
                plot3(omegaBb / obj.sys.omegaDamped, psiBb / pi, rhoBb, '--', 'Color', 'k', 'LineWidth', 2);
                xlabel('$\Omega / \omega_0$ [-]')
                ylabel('$\psi / \pi$ [-]')
                zlabel('$\rho$ [-]')
                xlim(minmax(omegaV) / obj.sys.omegaDamped)
                ylim(minmax(psiV) / pi)
                zlim(minmax(rhoV))
                view([0, 0])

                % Frequency-phase
                nexttile
                hold on; grid on; box on; axis square tight;
                patch(isosurface(omegaM / obj.sys.omegaDamped, psiM / pi, rhoM, M1, 0), ...
                    'EdgeColor', 'none', 'FaceColor', '#0072BD', 'FaceAlpha', 0.5);
                patch(isosurface(omegaM / obj.sys.omegaDamped, psiM / pi, rhoM, M2, 0), ...
                    'EdgeColor', 'none', 'FaceColor', '#EDB120', 'FaceAlpha', 0.5);
                for j = 1:length(frcData)
                    omegaFr = frcData{j}(:, 1);
                    psiFr = frcData{j}(:, 2);
                    rhoFr = frcData{j}(:, 3);
                    plot3(omegaFr / obj.sys.omegaDamped, psiFr / pi, rhoFr, '-', 'Color', '#D95319', 'LineWidth', 2);
                end
                plot3(omegaBb / obj.sys.omegaDamped, psiBb / pi, rhoBb, '--', 'Color', 'k', 'LineWidth', 2);
                xlabel('$\Omega / \omega_0$ [-]')
                ylabel('$\psi / \pi$ [-]')
                zlabel('$\rho$ [-]')
                xlim(minmax(omegaV) / obj.sys.omegaDamped)
                ylim(minmax(psiV) / pi)
                zlim(minmax(rhoV))
                view([0, 90])
            end
        end
    end

    methods(Static)
        function x = wrapPhase(x)
            % Wrap the phase angle x to the interval [-2*pi, 0].
            % Inputs:
            %   x: the phase angle x.
            % Outputs:
            %   x: the wrapped phase angle x.

            x(x > 0) = x(x > 0) - 2*pi;
        end

        function boundaries = extract_boundaries(C)
            % Extract the boundaries of the contour matrix C.
            % Inputs:
            %   C: the contour matrix.
            % Outputs:
            %   boundaries: a cell array containing the boundaries.

            % Initialize output
            boundaries = {};

            % Initialize counters
            idx = 1;
            bIdx = 1;

            % Loop over contour matrix
            while idx < size(C, 2)
                nPts = C(2, idx);
                boundaries{bIdx} = C(:, idx + 1:idx + nPts);
                idx = idx + nPts + 1;
                bIdx = bIdx + 1;
            end
        end
    end
end

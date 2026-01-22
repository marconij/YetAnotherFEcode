function [mu, orderNew, orderHistory] = run_fmincon_no_obj(mu0, A, b, Aeq, beq, xl, xu, args, options, varargin)
    % Run the optimization with 'fmincon'.
    % This allows to use global variables to update the expansion order.
    % Inputs:
    %   mu0: the initial design variables.
    %   A: the linear inequality constraint matrix.
    %   b: the linear inequality constraint vector.
    %   Aeq: the linear equality constraint matrix.
    %   beq: the linear equality constraint vector.
    %   xl: the lower bound vector.
    %   xu: the upper bound vector.
    %   args: the structure with the system and optimization arguments. It
    %       is used as a global variable to update the expansion order.
    %   options: the optimization options.
    %   doPlot: plot the backbone evolution during the optimization
    %       (optional, default is false).
    % Outputs:
    %   mu: the optimal design variables.
    %   orderNew: the new expansion order.

    % Parse inputs
    p = inputParser;
    p.addOptional('doPlot', false);
    parse(p, varargin{:});
    doPlot = p.Results.doPlot;

    % Define output function
    options.OutputFcn = @(x, optimValues, state) output_fun(x, optimValues, state, doPlot);

    % Initialize order history
    orderHistory = struct('iter', [], 'order', [], 'err', []);

    % Solve with fmincon
    mu = fmincon(@(x) fun(x), mu0, A, b, Aeq, beq, xl, xu, @(x) nonlcon(x), options);

    % Update expansion order
    orderNew = args.order;

    %% Objective function and sensitivities
    function [fval, fgrad] = fun(mu)
        % Objective function for optimization.
        % Inputs:
        %   mu: the design variables.
        % Outputs:
        %   fval: the objective function value.
        %   fgrad: the objective function gradient.

        fval = 0;
        fgrad = zeros(1, length(mu));
    end
    
    %% Nonlinear constraint function and sensitivities
    function [c, ceq, dc, dceq] = nonlcon(mu)
        % Nonlinear constraints function for optimization. This function
        % computes the equality and inequality constraints and their
        % sensitivities.
        % Inputs:
        %   mu: the design variables.
        % Outputs:
        %   c: the inequality constraints.
        %   ceq: the equality constraints.
        %   dc: the inequality constraints sensitivities.
        %   dceq: the equality constraints sensitivities.
        
        % Build system
        sys = MechSystem;
        sys.build_yafec_mpc_model(mu, args.yafec_mpc_constructor, ...
            'computePartialDerivatives', true);
        
        % Modal analysis
        [D, V] = sys.modal_analysis('modeRef', args.modeRefDrive, ...
            'computePartialDerivatives', true);
        sys.set_damping_rayleigh(args.rayleigh, ...
            'computePartialDerivatives', true);

        % Drive mode
        omegaDrive = sys.omega;

        % Sensitivity of the drive frequency
        dOmegaDrive = sys.adjoint_sensitivity_omega;

        % Sense mode
        senseIndex = modal_assurance_criterion(V, args.modeRefSense);
        phiSense = V(:, senseIndex);
        phiSense = phiSense / sqrt(phiSense.' * sys.M * phiSense); % mass normalization
        omegaSense = sqrt(D(senseIndex, senseIndex));

        % Sensitivity of the sense frequency
        dOmegaSense = zeros(size(dOmegaDrive));
        for i = 1:sys.ndv
            dOmegaSense(i) = phiSense.' * (-omegaSense^2 * sys.pM(i).coeffs + sys.pK(i).coeffs) * phiSense / (2 * omegaSense);
        end
    
        % Create SSM
        ssm = MechSSM(sys);
        ssm.compute_ssm(args.order, 'storeCoefficients', true);

        % Target point
        rhoT = ssm.compute_rho_from_x(args.xT, args.outDof);
        [~, xkT] = ssm.compute_x(rhoT, args.outDof);
        OmegaT = ssm.compute_omega(rhoT);

        % Sensitivity object
        sens = SensitivityBackbone(ssm);
        dOmegaT = sens.adjoint_sensitivity_Omega(rhoT, xkT, args.xT, args.outDof);
  
        % Equality constraints (1 x nCeq)
        ceq = [omegaDrive - args.omegaDriveT, omegaSense - args.omegaSenseT, OmegaT - args.OmegaT];
        dceq = [dOmegaDrive, dOmegaSense, dOmegaT];

        % Inequality constraints
        c = [];
        dc = [];
    end
    
    %% Output function for fmincon
    function stop = output_fun(mu, optimValues, state, doPlot)
        % Output function for the optimization with 'fmincon'.
        % Inputs:
        %   mu: the design variables.
        %   optimValues: the optimization values.
        %   state: the optimization state (init, iter, done).
        %   doPlot: plot the backbone evolution during the optimization.
        % Outputs:
        %   stop: the stop flag.
        
        % Output
        stop = false;
        
        % Iteration number
        iter = optimValues.iteration;
        
        % Build system
        sys = MechSystem;
        sys.build_yafec_mpc_model(mu, args.yafec_mpc_constructor);
        
        % Modal analysis
        [D, V] = sys.modal_analysis('modeRef', args.modeRefDrive);
        sys.set_damping_rayleigh(args.rayleigh);

        % Reference mode indices
        driveIndex = modal_assurance_criterion(V, args.modeRefDrive);
        senseIndex = modal_assurance_criterion(V, args.modeRefSense);

        % Referecne frequencies
        omegaDrive = sqrt(D(driveIndex, driveIndex));
        omegaSense = sqrt(D(senseIndex, senseIndex));
        
        % Create SSM
        ssm = MechSSM(sys);
        ssm.compute_ssm(args.order);

        % Compute omega
        rhoT = ssm.compute_rho_from_x(args.xT, args.outDof);
        OmegaT = ssm.compute_omega(rhoT);
        
        % Switch state
        switch state
            case 'init'
                % Display initial message
                disp('Starting optimization...')

                % Plot backbone
                if doPlot
                    figure();
                    hold on; grid on; box on; axis square tight;
                    plot_backbone(iter, ssm, args.outDof, rhoT, args.OmegaT, args.xT)
                end
            case 'iter'
                % Evaluate the residual error of the invariance equation
                err = ssm.error_measure(rhoT, 'style1stOrder', 'L2');

                % Update order history
                orderHistory.iter(end+1) = iter;
                orderHistory.order(end+1) = args.order;
                orderHistory.err(end+1) = err;
        
                % Display iteration and design variables
                fprintf('Iter %4d: fD = %6.2f kHz, fS = %6.2f kHz, fT = %6.2f kHz, err = %6.2e\n', ...
                    iter, omegaDrive / (2*pi) / 1e3, omegaSense / (2*pi) / 1e3, OmegaT / (2*pi) / 1e3, err)
        
                % Check the residual error. The order is increased if the
                % residual is too high and decreased if the residual is too
                % low. The order is not allowed to be less than 3 or
                % greater than args.maxOrder.
                warning('off', 'backtrace')
                if err > args.errTol && args.order < args.maxOrder
                    args.order = args.order + 2;
                    warning('Residual is %.2e. Order increased (%d->%d). Continuing.', err, args.order - 2, args.order)
                elseif err > args.errTol && args.order >= args.maxOrder
                    warning('Residual is %.2e but max order has been reached. Continuing.', err)
                elseif err < args.errTol / 100 && args.order >=5
                    args.order = args.order - 2;
                    warning('Residual is %.2e. Order reduced (%d->%d). Continuing.', err, args.order + 2, args.order)    
                end
        
                % Plot backbone
                if doPlot
                    plot_backbone(iter, ssm, args.outDof, rhoT, args.OmegaT, args.xT)
                end
            case 'done'
                % Display final message
                disp('Done!');
        
                % Plot backbone
                if doPlot
                    plot_backbone(iter, ssm, args.outDof, rhoT, args.OmegaT, args.xT)
                end
            otherwise
                % Nothing to do
        end
    end
    
    %% Function to plot the backbone
    function plot_backbone(iter, ssm, outDof, rhoT, OmegaT, xT)

        % Compute backbone
        rho = linspace(0, 1.1 * rhoT, 51);
        Omega = ssm.compute_omega(rho);
        x = ssm.compute_x(rho, outDof);

        % Plot backbone
        plot(Omega / (2*pi) / 1e3, x, 'k', 'LineWidth', 2)
        plot(OmegaT / (2*pi) / 1e3, xT, '.r', 'MarkerSize', 20)
        xlabel('$\Omega$ [kHz]')
        ylabel('$x$ [$\mu$m]')
        title(sprintf('Iteration %d', iter))
        legend('Backbone', 'Target Points', 'Interpreter', 'latex')
        drawnow;
    end
end

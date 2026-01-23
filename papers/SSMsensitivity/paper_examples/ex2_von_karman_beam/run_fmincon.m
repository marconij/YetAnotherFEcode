function [mu, orderNew, history] = run_fmincon(mu0, A, b, Aeq, beq, xl, xu, args, options, varargin)
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
    %   args: the arguments structure. This is a global variable.
    %   options: the optimization options.
    %   doPlot: plot the backbone evolution during the optimization (optional, default is true).
    % Outputs:
    %   mu: the optimal design variables.
    %   orderNew: the new expansion order.

    % Parse inputs
    p = inputParser;
    p.addOptional('doPlot', true);
    parse(p, varargin{:});
    doPlot = p.Results.doPlot;

    % Define output function
    options.OutputFcn = @(x, optimValues, state) output_fun(x, optimValues, state, doPlot);

    % Build system
    sys0 = MechSystem;
    sys0.build_yafec_model(mu0, args.yafec_assembly_constructor);
    sys0.modal_analysis('modeRef', args.modeRef);
    ssm0 = MechSSM(sys0);
    ssm0.compute_ssm(args.order);
    rho0 = ssm0.compute_rho_from_x(max(args.xTarget), args.outDof);
    rhoVec = linspace(0, rho0, 51);
    x0 = ssm0.compute_x(rhoVec, args.outDof);
    Omega0 = ssm0.compute_omega(rhoVec);

    % Create history
    history = struct('mu', mu0, 'Omega', Omega0, 'x', x0);

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

        fval = mu(2) * mu(4);
        fgrad = [0, mu(4), 0, mu(2)];
    end
    
    %% Nonlinear constraint function and sensitivities
    function [c, ceq, dc, dceq] = nonlcon(mu)
        % Nonlinear constraints function for optimization. This function computes
        % the equality and inequality constraints and their sensitivities.
        % Inputs:
        %   mu: the design variables.
        % Outputs:
        %   c: the inequality constraints.
        %   ceq: the equality constraints.
        %   dc: the inequality constraints sensitivities.
        %   dceq: the equality constraints sensitivities.
        
        % Build system
        sys = MechSystem;
        sys.build_yafec_model(mu, args.yafec_assembly_constructor, ...
            'computePartialDerivatives', true);
        
        % Modal analysis
        sys.modal_analysis('modeRef', args.modeRef, ...
            'computePartialDerivatives', true);
        sys.modal_analysis_sensitivity;
        
        % Create SSM
        ssm = MechSSM(sys);
        ssm.compute_ssm(args.order, 'computeSensitivity', true);

        % Exctract nonlinear constraints
        outDof = args.outDof;
        xTarget = args.xTarget;
        OmegaTarget = args.OmegaTarget;
        coeffTarget = args.coeffTarget;
    
        % Find number of equality and inequality constraints
        iCeq = find(coeffTarget == 0);
        iC = find(abs(coeffTarget) == 1);
    
        % Resize vectors
        ceq = zeros(1, length(iCeq));
        dceq = zeros(length(mu), length(iCeq));
        c = zeros(1, length(iC));
        dc = zeros(length(mu), length(iC));

        % Loop over equality constraints
        for ii = 1:length(iCeq)
            % Chek if the target point is equal to 0
            if xTarget(iCeq(ii)) == 0
                % Compute target points and sensitivities
                fval = sys.omegaDamped;
                fgrad = sys.domegaDamped;
            else
                % Compute target points and sensitivities
                rho = ssm.compute_rho_from_x(xTarget(iCeq(ii)), outDof);
                [~, xk] = ssm.compute_x(rho, outDof);
                fval = ssm.compute_omega(rho);
                drho = ssm.compute_rho_derivative(rho, outDof, xk);
                fgrad = ssm.compute_omega_derivative(rho, drho);
            end
    
            % Store equality constraints and sensitivities
            ceq(ii) = fval - OmegaTarget(iCeq(ii));
            dceq(:, ii) = fgrad;
        end
    
        % Loop over inequality constraints
        for ii = 1:length(iC)
            % Chek if the target point is equal to 0
            if xTarget(iC(ii)) == 0
                % Compute target points and sensitivities
                fval = sys.omegaDamped;
                fgrad = sys.domegaDamped;
            else
                % Compute target points and sensitivities
                rho = ssm.compute_rho_from_x(xTarget(iC(ii)), outDof);
                [~, xk] = ssm.compute_x(rho, outDof);
                fval = ssm.compute_omega(rho);
                drho = ssm.compute_rho_derivative(rho, outDof, xk);
                fgrad = ssm.compute_omega_derivative(rho, drho);
            end
    
            % Store inequality constraints and sensitivities
            c(ii) = (fval - OmegaTarget(iC(ii))) * coeffTarget(iC(ii));
            dc(:, ii) = fgrad * coeffTarget(iC(ii));
        end
    end
    
    %% Output function for fmincon
    function stop = output_fun(mu, optimValues, state, doPlot)
        % Output function for the optimization with 'fmincon'.
        % Inputs:
        %   mu: the design variables.
        %   optimValues: the optimization values.
        %   state: the optimization state (init, iter, done).
        %   doPlot: plot the backbone backbone evolution during the optimization.
        % Outputs:
        %   stop: the stop flag.
        
        % Output
        stop = false;
        
        % Iteration number
        iter = optimValues.iteration;
        
        % Build system
        sys = MechSystem;
        sys.build_yafec_model(mu, args.yafec_assembly_constructor);
        
        % Modal analysis
        sys.modal_analysis('modeRef', args.modeRef);
        
        % Create SSM
        ssm = MechSSM(sys);
        ssm.compute_ssm(args.order);
        
        % Compute omega
        xMax = max(args.xTarget);
        rhoMax = ssm.compute_rho_from_x(xMax, args.outDof);
        rho = linspace(0, rhoMax, 51);
        x = ssm.compute_x(rho, args.outDof);
        Omega = ssm.compute_omega(rho);

        % Store history
        history(iter + 2).mu = mu;
        history(iter + 2).Omega = Omega;
        history(iter + 2).x = x;
        
        % Switch state
        switch state
            case 'init'
                % Initialize backbone plot
                if doPlot
                    figure();
                    hold on; grid on; box on; axis square tight;
                    plot(Omega, x, 'k', 'LineWidth', 2)
                    plot(args.OmegaTarget, args.xTarget, '.r', 'MarkerSize', 20)
                    xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
                    ylabel('$x$ [m]', 'Interpreter', 'latex')
                    title('Physical amplitude', 'Interpreter', 'latex')
                    legend('Backbone', 'Target Points', 'Interpreter', 'latex')
                    set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
                    drawnow;
                end
            case 'iter'
                % Evaluate the residual error of the invariance equation
                err = ssm.error_measure(rhoMax, 'style1stOrder', 'L2');
        
                % Tolerance on the residual error (adjust this value if needed)
                % tol = 1e-1;
        
                % Check the residual error. The order is increased if the residual is too high
                % and decreased if the residual is too low. The order is not allowed to be less
                % than 3 or greater than args.maxOrder.
                warning('off', 'backtrace')
                if err > args.tol && args.order < args.maxOrder
                    args.order = args.order + 2;
                    warning('Residual is %.2e. Order increased (%d->%d). Continuing.', err, args.order - 2, args.order)
                elseif err > args.tol && args.order >= args.maxOrder
                    warning('Residual is %.2e but max order has been reached. Continuing.', err)
                elseif err < args.tol / 100 && args.order >=5
                    args.order = args.order - 2;
                    warning('Residual is %.2e. Order reduced (%d->%d). Continuing.', err, args.order + 1, args.order)    
                end
        
                % Display iteration and design variables
                fprintf('Iter %d: [%.2e, %.2e, %.2e, %.2e], res: %.2e\n', iter, mu, err)
        
                % Plot backbone
                if doPlot
                    plot(Omega, x, 'k', 'LineWidth', 2)
                    plot(args.OmegaTarget, args.xTarget, '.r', 'MarkerSize', 20)
                    xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
                    ylabel('$x$ [m]', 'Interpreter', 'latex')
                    title('Physical amplitude', 'Interpreter', 'latex')
                    legend('Backbone', 'Target Points', 'Interpreter', 'latex')
                    set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
                    drawnow;
                end
            case 'done'
                % Display final iteration and design variables
                disp(['Done: ', num2str(mu)]);
        
                % Plot backbone
                if doPlot
                    plot(Omega, x, 'k', 'LineWidth', 2)
                    plot(args.OmegaTarget, args.xTarget, '.r', 'MarkerSize', 20)
                    xlabel('$\Omega$ [rad/s]', 'Interpreter', 'latex')
                    ylabel('$x$ [m]', 'Interpreter', 'latex')
                    title('Physical amplitude', 'Interpreter', 'latex')
                    legend('Backbone', 'Target Points', 'Interpreter', 'latex')
                    set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
                    drawnow;
                end
            otherwise
                % Nothing to do
        end
    end
end

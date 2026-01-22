function [xsol,fval,history] = run_fmincon(x0, lb, ub, args)
% This function is a wrapper for fmincon with a custom output function
% to track the optimization history and with a mode tracking algorithm
% which requires the global variable args.imod to be updated.
% [for documentation see: Output Functions for Optimization Toolbox
% (https://www.mathworks.com/help/optim/ug/output-functions.html)]
%
% INPUTS:
%   x0: initial design variables
%   lb: lower bound on design variables
%   ub: upper bound on design variables
%   args: structure containing the following fields
%       .sys_fun: function handle to the system function
%       .imod: target mode indices
%       .f_target: target frequencies for the mode shapes
%       .imod_gamma: index of the mode used to compute gamma
%       .gamma_target: treshold value for the gamma parameter
%       .ConstraintTolerance: fmincon constraint tolerance value
% OUTPUTS:
%   xsol: optimal design variables
%   fval: objective function value at the optimal design variables
%   history: structure containing the optimization history

% Set up shared variables with outfun
history.x = [];     % iterates
history.V = {};     % eigenvectors
history.f0 = [];    % eigenfrequencies    
history.imod = [];  % mode indices

% Options for fmincon
options = optimoptions('fmincon', 'Display', 'iter-detailed');
options.Algorithm = "interior-point";
options.EnableFeasibilityMode = true;
options.SpecifyObjectiveGradient = true;
options.SpecifyConstraintGradient = true;
options.CheckGradients = false;
options.FiniteDifferenceType = 'central';
options.OutputFcn = @(x, optimValues, state) output_fun(x, optimValues, state);

% options.ObjectiveLimit = 1e-3; 
options.MaxIterations = 200;
options.ConstraintTolerance = args.ConstraintTolerance;

args.options = options;
args.dp = 1e-4; % finite difference step for dK and dM

% Linear constraints
A = []; b = [];
Aeq = []; beq = [];

[xsol, fval] = fmincon(@(x) obj_fun(x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlcon(x), options);

    function stop = output_fun(x, optimValues, state)
        % Output function for fmincon. This function is called at each
        % iteration of the optimization algorithm. It updates the history
        % structure and plots the optimization progress.
        % INPUTS:
        %   x: design variables
        %   optimValues: structure containing information about the current
        %                iteration
        %   state: current state of the optimization algorithm
        % OUTPUTS:
        %   stop: flag to stop the optimization algorithm

        stop = false;
        iter = optimValues.iteration + 1;
        history.x(iter,:) = x;
        history.V{iter} = args.Vold;
        history.f0(iter,:) = args.f0(args.imod);
        history.imod(iter,:) = args.imod;
        history.imod_gamma(iter) = args.imod_gamma;
        history.gamma(iter) = args.gamma;
        
        switch state
            case 'init'
                figure
                hold on
                xlabel('iter','fontsize',16)
                yyaxis left
                ylabel('f_0 [kHz]','fontsize',16)
                grid on; grid minor;
                yyaxis right
                ylabel('\gamma')
            case 'iter'
                yyaxis left
                cla
                plot(0:iter-1, history.f0/1000)
                yyaxis right
                cla
                semilogy(0:iter-1, history.gamma,'.-')
            case 'done'
                axis tight
                legend('Drive','Sense')
                fclose('all');
                
                options.SpecifyConstraintGradient = false;
                [c, ceq] = nonlcon(x);
                disp(' ')
                for ii = 1:length(c)
                    fprintf('Inequality constraint %d: %g\n', ii, c(ii))
                end
                for ii = 1:length(ceq)
                    fprintf('Equality constraint %d: %g\n', ii, ceq(ii))
                end
                disp(optimValues)
            otherwise
        end
        drawnow
    end

    function [fval, fgrad] = obj_fun(x)
        % Objective function for fmincon. Since the target frequencies are
        % enforced as equality constraints, and that the gamma parameter
        % is enforced as an inequality constraint, the objective function
        % is set to zero. The gradient of the objective function is also
        % zero.
        % INPUTS:
        %   x: design variables
        fval = 0;
        fgrad = zeros(1,15);
    end

    function [c, ceq, dc, dceq] = nonlcon(x)
        % Nonlinear constraint function for fmincon. The equality
        % constraints enforce the target frequencies. The gradient of the
        % equality constraints is computed if the option is enabled. The
        % inequality constraint enforces the gamma parameter to be less
        % than the target value. The gradient of the inequality constraint
        % is computed if the option is enabled.
        % INPUTS:
        %   x: design variables
        % OUTPUTS:
        %   c: inequality constraints
        %   ceq: equality constraints
        %   dc: gradient of inequality constraints
        %   dceq: gradient of equality constraints        

        sys = args.sys_fun(x);
        M = sys.mass_matrix;
        K1 = sys.tangent_stiffness_and_force(zeros(sys.nDOFs.all_structural,1));

        % Eigenfrequencies and Eigenvectors
        [V,D] = eigs(M\K1, max([3, args.imod+3]), 'smallestabs');
        om0 = sqrt(diag(D));
        [om0, ind] = sort(om0);
        V = V(:,ind);
        for ii = 1 : size(V,2)
            V(:,ii) = V(:,ii)/sqrt(V(:,ii)'*M*V(:,ii));
        end
        
        % MAC to track target modes
        if isfield(args,'Vold')
            track_modes(V, args.Vold); % updates args.imod
        end
        args.Vold = V;
        args.f0 = om0/2/pi;

        % Finite Differences (derivatives of K1, K2, K3 and M)
        [K2, K3] = sys.mpc_tensors;
        [dK1dp, dK2dp, dK3dp, dMdp] = finite_differences(x, args, K1, K2, K3, M);
        
        % Equality Constraints
        ceq = zeros(length(args.imod), 1);
        dceq = zeros(length(x), length(args.imod));
        for mm = 1 : length(args.imod)
            m = args.imod(mm);
            q = V(:,m)/sqrt(V(:,m)'*M*V(:,m));

            omT = 2*pi*args.f_target(mm);
            omi = om0(m);
            ceq(mm) = (omi/omT-1)^2;

            % Sensitivities
            if args.options.SpecifyConstraintGradient
                dOmega = zeros(1,length(x));
                for ii = 1 : length(x)
                    % eigenfrequency sensitivities
                    dOmega(ii) = 1/(2*omi)*q'*(dK1dp{ii} - omi^2*dMdp{ii})*q;
                end
                dceq(:, mm) = 2*(omi/omT-1)/omT*dOmega;
            else
                dceq = [];
            end
        end

        % inequality constraint on gamma
        if args.options.SpecifyConstraintGradient
            [gamma, dgamma] = gamma_sensitivity(M, K1, K2, K3, ...
                dMdp, dK1dp, dK2dp, dK3dp, ...
                om0(args.imod_gamma), V(:,args.imod_gamma));
            dc = dgamma;
        else
            % for finite differences
            gamma = gamma_value(M, K1, K2, K3, ...
                om0(args.imod_gamma), V(:,args.imod_gamma));
            dc = 0;
        end
        args.gamma = gamma;
        c = gamma/abs(args.gamma_target) - sign(args.gamma_target);
        dc = dc/abs(args.gamma_target);
    end

    function [dK1dp, dK2dp, dK3dp, dMdp] = finite_differences(x, args, K10, K20, K30, M0)
        % Finite differences for the derivatives of the mass and stiffness
        % matrices with respect to the design variables. The derivatives
        % are computed using the forward difference method.
        % INPUTS:
        %   x: design variables
        %   args: structure containing the system function handle and
        %         options (SpecifyConstraintGradient=true)
        %   K10: initial linear stiffness matrix
        %   K20: initial quadratic stiffness tensor
        %   K30: initial cubic stiffness tensor
        %   M0: initial mass matrix
        % OUTPUTS:
        %   dK1dp: derivatives of the linear stiffness matrix (cell array)
        %   dK2dp: derivatives of the quadratic stiffness tensor (cell array)
        %   dK3dp: derivatives of the cubic stiffness tensor (cell array)
        %   dMdp: derivatives of the mass matrix (cell array)
        if args.options.SpecifyConstraintGradient
            dK1dp = cell(length(x),1);
            dK2dp = cell(length(x),1);
            dK3dp = cell(length(x),1);
            dMdp = cell(length(x),1);
            for ii = 1 : length(x)
                xp = x;
                xp(ii) = x(ii)+args.dp;
                sys = args.sys_fun(xp);
                Mp = sys.mass_matrix;
                Kp = sys.tangent_stiffness_and_force(zeros(sys.nDOFs.all_structural,1));
                [K2p, K3p] = sys.mpc_tensors;
                dK1dp{ii} = (Kp-K10)/args.dp;
                dMdp{ii} = (Mp-M0)/args.dp;
                dK2dp{ii} = (K2p-K20)/args.dp;
                dK3dp{ii} = (K3p-K30)/args.dp;
            end
        else
            dK1dp = []; dK2dp = []; dK3dp = []; dMdp = [];
        end
    end

    function track_modes(Vnew, Vold)
        % Mode tracking algorithm based on the modal assurance criterion
        % (MAC) between the current and previous eigenvectors. The mode
        % indices are updated in the global variable args.imod. 
        % args.imod_gamma is also updated accordingly.
        % INPUTS:
        %   Vnew: current eigenvectors
        %   Vold: previous eigenvectors
        MAC = zeros(size(Vnew,2),1);
        imod_old = args.imod;
        for mm = 1 : length(args.imod)
            jj = args.imod(mm);
            for ii = 1 : size(Vnew,2)
                MAC(ii) = norm(Vold(:,jj)'*Vnew(:,ii))/((Vold(:,jj)'*Vold(:,jj))*(Vnew(:,ii)'*Vnew(:,ii)));
            end
            [~,ind] = max(MAC);
            if ind~=args.imod(mm)
                warning('mode veering %d-->%d', args.imod(mm), ind)
                args.imod(mm) = ind;
            end
        end
        args.imod_gamma = imod_old(args.imod==args.imod_gamma);
    end

end
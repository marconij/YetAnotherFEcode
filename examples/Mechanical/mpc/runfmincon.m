function [xsol,fval,history] = runfmincon(x0, A, b, Aeq, beq, lb, ub, args)
% for documentation see: Output Functions for Optimization Toolbox
% This function is a wrapper for fmincon with a custom output function
% to track the optimization history and with a mode tracking algorithm
% which requires the global variable args.imod to be updated.
% INPUTS:
%   x0: initial design variables
%   A: linear inequality constraint matrix
%   b: linear inequality constraint vector
%   Aeq: linear equality constraint matrix
%   beq: linear equality constraint vector
%   lb: lower bound on design variables
%   ub: upper bound on design variables
%   args: structure containing the following fields
%       .sys_fun: function handle to the system function
%       .f_target: target frequencies
%       .desired_percent_error: desired percent error in the eigenfrequencies
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
options.MaxIterations = 50;
options.ConstraintTolerance = (args.desired_percent_error/100)^2;

args.options = options;
args.dp = 1e-4; % finite difference step for dK and dM

[xsol, fval] = fmincon(@(x) obj_fun(x), x0, A, b, Aeq, beq, lb, ub, @(x) nonlcon(x), options);

    function stop = output_fun(x, optimValues, state)
        % Output function for fmincon. This function is called at each
        % iteration of the optimization algorithm. It updates the history
        % structure and the mode indices. The mode indices are updated
        % based on the modal assurance criterion (MAC) between the current
        % and previous eigenvectors. f0 ordering follows the mode indices.
        % The function also plots the eigenfrequencies at each iteration.
        % INPUTS:
        %   x: design variables
        %   optimValues: structure containing information about the current
        %                iteration
        %   state: current state of the optimization algorithm
        stop = false;
        iter = optimValues.iteration + 1;
        history.x(iter,:) = x;
        history.V{iter} = args.Vold;
        history.f0(iter,:) = args.f0(args.imod);
        history.imod(iter,:) = args.imod;
        
        switch state
            case 'init'
                figure
                hold on
                xlabel('iter')
                ylabel('f_0 [kHz]')
                grid on; grid minor;
            case 'iter'
                cla
                plot(history.f0/1000,'-o')
            case 'done'
                axis tight
                legend('Drive','Sense')
            otherwise
        end
        drawnow
    end

    function [fval, fgrad] = obj_fun(x)
        % Objective function for fmincon. Since the target frequencies are
        % enforced as equality constraints, the objective function is set to
        % zero. The gradient of the objective function is also zero.
        % INPUTS:
        %   x: design variables
        fval = 0;
        fgrad = zeros(1,15);
    end

    function [c, ceq, dc, dceq] = nonlcon(x)
        % Nonlinear constraint function for fmincon. The equality
        % constraints enforce the target frequencies. The gradient of the
        % equality constraints is computed if the option is enabled.
        % INPUTS:
        %   x: design variables
        % OUTPUTS:
        %   c: inequality constraints
        %   ceq: equality constraints
        %   dc: gradient of inequality constraints
        %   dceq: gradient of equality constraints
        c = 0;
        dc = zeros(length(x), 1);

        sys = args.sys_fun(x);
        M0 = sys.mass_matrix;
        K0 = sys.tangent_stiffness_and_force(zeros(sys.nDOFs.all_structural,1));

        % Eigenfrequencies and Eigenvectors
        [V,D] = eigs(M0\K0, max([3, args.imod+1]), 'smallestabs');
        om0 = sqrt(diag(D));
        [om0, ind] = sort(om0);
        V = V(:,ind);
        
        % MAC to track target modes
        if isfield(args,'Vold')
            track_modes(V, args.Vold); % updates args.imod
        end
        args.Vold = V;
        args.f0 = om0/2/pi;

        % Finite Differences (derivatives of M and K)
        [dKdp, dMdp] = finite_differences(x, args, K0, M0);
        
        % Equality Constraints
        ceq = zeros(length(args.imod), 1);
        dceq = zeros(length(x), length(args.imod));
        for mm = 1 : length(args.imod)
            m = args.imod(mm);
            q = V(:,m)/sqrt(V(:,m)'*M0*V(:,m));

            omT = 2*pi*args.f_target(mm);
            omi = om0(m);
            ceq(mm) = (omi/omT-1)^2;

            % Ssensitivities
            if args.options.SpecifyConstraintGradient
                dOmega = zeros(1,length(x));
                for ii = 1 : length(x)
                    % eigenfrequency sensitivities
                    dOmega(ii) = 1/(2*omi)*q'*(dKdp{ii} - omi^2*dMdp{ii})*q;
                end
                dceq(:, mm) = 2*(omi/omT-1)/omT*dOmega;
            else
                dceq = [];
            end
        end
        % Print constraints at each iteration
        s = "\b\t ceq:";
        for mm = 1 : length(args.imod)
            s = s+sprintf(' %g', ceq(mm));
        end
        fprintf(s+"\n")
    end

    function [dKdp, dMdp] = finite_differences(x, args, K0, M0)
        % Finite differences for the derivatives of the mass and stiffness
        % matrices with respect to the design variables. The derivatives
        % are computed using the forward difference method.
        % INPUTS:
        %   x: design variables
        %   args: structure containing the system function handle and
        %         options (SpecifyConstraintGradient=true)
        %   K0: initial stiffness matrix
        %   M0: initial mass matrix
        % OUTPUTS:
        %   dKdp: derivatives of the stiffness matrix (cell array)
        %   dMdp: derivatives of the mass matrix (cell array)
        if args.options.SpecifyConstraintGradient
            dKdp = cell(length(x));
            dMdp = cell(length(x));
            for ii = 1 : length(x)
                xp = x;
                xp(ii) = x(ii)+args.dp;
                sys = args.sys_fun(xp);
                Mp = sys.mass_matrix;
                Kp = sys.tangent_stiffness_and_force(zeros(sys.nDOFs.all_structural,1));
                dKdp{ii} = (Kp-K0)/args.dp;
                dMdp{ii} = (Mp-M0)/args.dp;
            end
        end
    end

    function track_modes(Vnew, Vold)
        % Mode tracking algorithm based on the modal assurance criterion
        % (MAC) between the current and previous eigenvectors. The mode
        % indices are updated in the global variable args.imod.
        % INPUTS:
        %   Vnew: current eigenvectors
        %   Vold: previous eigenvectors
        MAC = zeros(size(Vnew,2),1);
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
    end

end
function rho = compute_rho_from_x(obj, xTarget, outDof)
    % Compute the reduced amplitude rho from the RMS of the physical
    % amplitude x of the output dof.
    % Inputs:
    %   xTarget: RMS of the target physical amplitude of the output dof.
    %   outDof: the index of the output dof.
    % Outputs:
    %   rho: reduced amplitude.

    % Define the function to be solved
    % The computation of the partial derivatives is temporarily disabled
    fun = @(rho) obj.compute_x(rho, outDof);

    % Initialize the reduced amplitude
    % The initial guess is computed at the leading order, for which
    % an analytical expression is available.
    rho = obj.compute_rho_from_x_leading_order(xTarget, outDof);

    % Solve for each target
    for i = 1:length(xTarget)
        if xTarget(i) == 0
            rho(i) = 0;
        else
            % rho(i) = fsolve(@(x) fun(x) - xTarget(i), rho(i)); % use fsolve
            rho(i) = fzero(@(x) fun(x) - xTarget(i), rho(i)); % use fzero
        end
    end
end

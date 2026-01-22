function [resNorm, norms] = compute_invariance_residual(obj, p, style1stOrder)
    % Compute the residual error of the invariance equation.
    % Inputs:
    %   p: the reduced coordinates (2 x 1).
    %   style1stOrder: the style of the 1st order system matrices.
    %       Valid values are 'L1' and 'L2'.
    % Outputs:
    %   resNorm: the norm of the residual error.
    %   norms: the norms of the terms of the invariance equation. Note that
    %       the norms of the terms depend on the style of the 1st order
    %       system matrices.

    % Compute 1st order order coefficients W = [w; dw] and DW = [Dw; Ddw]
    W = obj.compute_W(p);
    DW = obj.compute_dW_t(p);

    % Extract 2nd order coefficients
    w = W(1:end/2, 1);
    dw = W(end/2+1:end, 1);
    Dw = DW(1:end/2, 1);
    Ddw = DW(end/2+1:end, 1);

    % Compute the 2nd order nonlinear force vector
    f = obj.sys.compute_nonlinear_force(w, dw);

    % Compute the terms of the 1st order invariance equation according to
    % the 1st order system style
    if strcmp(style1stOrder, 'L1') % SSM sensitivity paper
        BDW = [obj.sys.C * Dw + obj.sys.M * Ddw; obj.sys.M * Dw];
        AW = [-obj.sys.K * w; obj.sys.M * dw];
        FW = [-f; zeros(size(f))];

        % % For debugging
        % B = [obj.sys.C, obj.sys.M; 
        %     obj.sys.M, sparse(obj.sys.n, obj.sys.n)];
        % A = blkdiag(-obj.sys.K, obj.sys.M);
        % errBDW = norm(BDW - B * DW);
        % errAW = norm(AW - A * W);
    elseif strcmp(style1stOrder, 'L2') % LSM optimization paper
        BDW = [-obj.sys.K * Dw; obj.sys.M * Ddw];
        AW = [-obj.sys.K * dw; -obj.sys.K * w - obj.sys.C * dw];
        FW = [zeros(size(f)); -f];

        % % For debugging
        % B = blkdiag(-obj.sys.K, obj.sys.M);
        % A = [sparse(obj.sys.n, obj.sys.n), -obj.sys.K; 
        %     -obj.sys.K, -obj.sys.C];
        % errBDW = norm(BDW - B * DW);
        % errAW = norm(AW - A * W);
    else
        error('Invalid style1stOrder value. Valid values are: ''L1'', ''L2''.');
    end

    % Compute the residual error of the invariance equation
    resNorm = norm(BDW - AW - FW);
    norms = [norm(BDW), norm(AW), norm(FW)];
end

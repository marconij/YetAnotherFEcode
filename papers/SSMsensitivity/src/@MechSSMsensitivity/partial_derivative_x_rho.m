function px_rho = partial_derivative_x_rho(obj, rho, xk, x, outDof)
    % Compute the partial derivative of the physical amplitude x with
    % respect to the reduced amplitude rho.
    % Inputs:
    %   rho: the reduced amplitude.
    %   xk: the physical amplitudes at all the time steps.
    %   x: the rms of the physical amplitude.
    %   outDof: the index of the output dof.
    % Outputs:
    %   px_rho: the partial derivative of x with respect to rho.

    % Reduced coordinates
    p1Tilde = exp( 1i * obj.ssm.theta);
    p2Tilde = exp(-1i * obj.ssm.theta);
        
    % Derivative of xk with respect to rho
    pxk_rho = obj.ssm.w(1).coeffs(outDof, 1) * (p1Tilde + p2Tilde);

    % Loop over orders
    for mOrder = 2:obj.ssm.maxOrder
        % Loop over the multi-indices at the current order
        for i = 1:obj.ssm.MIs(mOrder).n
            % Current multi-index
            m = obj.ssm.MIs(mOrder).coeffs(:, i);

            % Partial derivative of xk with respect to rho
            pxk_rho = pxk_rho + obj.ssm.w(mOrder).coeffs(outDof, i) * (p1Tilde.^m(1) .* p2Tilde.^m(2)) * mOrder * rho^(mOrder - 1);
        end
    end

    % Complete the partial derivative
    px_rho = sum(real(pxk_rho) .* xk) / (obj.ssm.nTimeSteps * x);
end

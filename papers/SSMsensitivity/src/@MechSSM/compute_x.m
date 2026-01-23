function [x, xk] = compute_x(obj, rho, outDof)
    % Compute the RMS of the physical amplitude x(rho_i, theta_k) of the
    % output dof.
    % The RMS is performed along the phase dimension. The result is a
    % vector of RMS values, one for each reduced amplitude rho.
    % Inputs:
    %   rho: the vector of reduced amplitudes (nRho x 1).
    %   outDof: the index of the output dof.
    % Outputs:
    %   x: RMS of the physical amplitude of the output dof.
    %   xk: physical amplitudes of the output dof.

    % Make sure rho is a column vector and theta is a row vector
    rho = rho(:);

    % Normalized reduced coordinates
    p1Tilde = exp(1i * obj.theta);
    p2Tilde = exp(-1i * obj.theta);

    % Initialize the response with the first order term
    xk = (obj.w(1).coeffs(outDof, 1) * p1Tilde + obj.w(1).coeffs(outDof, 2) * p2Tilde) .* rho;

    % Loop over orders
    for mOrder = 2:obj.maxOrder
        % Loop over the multi-indices at the current order
        for i = 1:obj.MIs(mOrder).n
            % Current multi-index
            m = obj.MIs(mOrder).coeffs(:, i);

            % Add the contribution of the current multi-index
            xk = xk + obj.w(mOrder).coeffs(outDof, i) * p1Tilde.^m(1) .* p2Tilde.^m(2) .* rho.^mOrder;
        end
    end

    % Make sure the response is a real number
    xk = real(xk);

    % Compute the norm of the output dof
    x = sqrt(sum(xk.^2, 2) / obj.nTimeSteps);
end

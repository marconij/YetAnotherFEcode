function drho = compute_rho_derivative(obj, rho, outDof, xk)
    % Compute the derivative of the reduced amplitude at the output dof.
    % Inputs:
    %   rho: the reduced amplitude.
    %   outDof: the index of the output dof.
    %   xk: the physical amplitudes of the output dof.
    % Outputs:
    %   drho: the derivative of the reduced amplitude at the output dof.

    % Number of time steps
    nTheta = obj.nTheta;

    % Normalized reduced coordinatres
    p1Tilde = exp(1i * obj.theta);
    p2Tilde = exp(-1i * obj.theta);

    % Initialize the numerator and denominator
    num = zeros(obj.sys.ndv, nTheta);
    for dv = 1:obj.sys.ndv
        num(dv, :) = (obj.w(1).dv(dv).coeffs(outDof, 1) * p1Tilde + obj.w(1).dv(dv).coeffs(outDof, 2) * p2Tilde) .* rho;
    end
    den = obj.w(1).coeffs(outDof, 1) * p1Tilde + obj.w(1).coeffs(outDof, 2) * p2Tilde;

    % Loop over orders
    for mOrder = 2:obj.maxOrder
        % Loop over the multi-indices at the current order
        for mIdx = 1:obj.MIs(mOrder).n
            % Current multi-index
            m = obj.MIs(mOrder).coeffs(:, mIdx);

            % Update the denominator
            wm = obj.w(mOrder).coeffs(outDof, mIdx);
            den = den + wm * p1Tilde.^m(1) .* p2Tilde.^m(2) .* mOrder .* rho.^(mOrder - 1);

            % Loop over design variables
            for dv = 1:obj.sys.ndv
                % Extract derivative of w
                dwm = obj.w(mOrder).dv(dv).coeffs(outDof, mIdx);

                % Update the numerator
                num(dv, :) = num(dv, :) + dwm * p1Tilde.^m(1) .* p2Tilde.^m(2) .* rho.^mOrder;
            end
        end
    end

    % Complete the denominator
    den = sum(xk .* real(den), 2); % sum over theta

    % Complete the numerator and assemble the sensitivity
    % The result should be a column vector (nRho x 1)
    drho = zeros(obj.sys.ndv, 1);
    for dv = 1:obj.sys.ndv
        drho(dv) = -sum(xk .* real(num(dv, :)), 2) ./ den;
    end
end

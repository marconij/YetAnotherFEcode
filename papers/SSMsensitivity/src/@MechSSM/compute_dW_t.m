function DW = compute_dW_t(obj, p)
    % Compute the time derivative of the manifold parametrization in the
    % 1st order form.
    % Inputs:
    %   p: the reduced coordinates.
    % Outputs:
    %   DW: the time derivative of the manifold parametrization in the 1st
    %       order form.

    % Compute the time derivative of the reduced coordinates
    Dp = obj.compute_dp_t(p);

    % Loop over orders
    dWdp = zeros(2*obj.sys.n, 2);
    for mOrder = 1 : obj.maxOrder
        % Loop over multi-indices
        for mIdx = 1 : obj.MIs(mOrder).n
            % Current multi-index
            m = obj.MIs(mOrder).coeffs(:, mIdx);

            % Assemble manifold parametrization (1st order)
            Wm = [obj.w(mOrder).coeffs(:, mIdx); obj.dw(mOrder).coeffs(:, mIdx)];

            % Compute the partial derivatives with respect to the reduced coordinates.
            % If the associated multi-index component is zero, the derivative is zero.
            if m(1) > 0
                dWdp(:, 1) = dWdp(:, 1) + m(1) * Wm * p(1)^(m(1)-1) * p(2)^m(2);
            end
            if m(2) > 0
                dWdp(:, 2) = dWdp(:, 2) + m(2) * Wm * p(1)^m(1) * p(2)^(m(2)-1);
            end
        end
    end

    % Compute the time derivative of the manifold parametrization (1st order)
    % Also, ensure real values
    DW = real(dWdp * Dp);
end

function W = compute_W(obj, p)
    % Compute the 1st order form manifold parametrization: W = [w; dw].
    % Inputs:
    %   p: the reduced coordinates.
    % Outputs:
    %   W: the 1st order form manifold parametrization.

    % Loop over orders
    W = zeros(2*obj.sys.n, 1);
    for mOrder = 1 : obj.maxOrder
        % Loop over multi-indices
        for mIdx = 1 : obj.MIs(mOrder).n
            % Current multi-index
            m = obj.MIs(mOrder).coeffs(:, mIdx);

            % Assemble manifold parametrization (1st order)
            Wm = [obj.w(mOrder).coeffs(:, mIdx); obj.dw(mOrder).coeffs(:, mIdx)];

            % Sum contribution
            W = W + Wm * p(1)^m(1) * p(2)^m(2);
        end
    end

    % Ensure real values
    W = real(W);
end

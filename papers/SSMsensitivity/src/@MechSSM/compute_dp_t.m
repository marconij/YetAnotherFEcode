function Dp = compute_dp_t(obj, p)
    % Compute the time derivative of the reduced coordinates: Dp = R(p).
    % Inputs:
    %   p: the reduced coordinates.
    % Outputs:
    %   R: the time derivative of the reduced coordinates.

    % Loop over orders
    Dp = zeros(2, 1);
    for mOrder = 1:2:obj.maxOrder
        % Loop over multi-indices
        for mIdx = 1 : obj.MIs(mOrder).n
            % Current multi-index
            m = obj.MIs(mOrder).coeffs(:, mIdx);

            % Compute the time derivative of the reduced coordinates.
            % The only non-zero terms are the ones with m(1) - m(2) = 1 or -1.
            % If m(1) - m(2) = 1, the term is associated with the first reduced coordinate.
            % If m(1) - m(2) = -1, the term is associated with the second reduced coordinate.
            if m(1) - m(2) == 1
                Dp(1) = Dp(1) + obj.R(mOrder).coeffs(1) * p(1)^m(1) * p(2)^m(2);
            elseif m(1) - m(2) == -1
                Dp(2) = Dp(2) + obj.R(mOrder).coeffs(2) * p(1)^m(1) * p(2)^m(2);
            end
        end
    end
end

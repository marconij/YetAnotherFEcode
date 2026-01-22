function Omega = compute_omega(obj, rho)
    % Compute the response frequency Omega at the given reduced amplitudes
    % rho.
    % Inputs:
    %   rho: the reduced amplitude or vector of reduced amplitudes.
    % Outputs:
    %   Omega: the response frequency.

    % Initialize the response frequency with the first order term
    % The 1st order reduced dynamics coefficients are the eigenvalues of the system
    % Omega = 1i/2 * sum(obj.R(1).coeffs(2, :) - obj.R(1).coeffs(1, :)) * ones(size(rho));
    Omega = 1i/2 * (conj(obj.sys.lambda) - obj.sys.lambda) * ones(size(rho));

    % Loop over odd orders (the coefficients of even orders are zero)
    for mOrder = 3:2:obj.maxOrder
        R1 = obj.R(mOrder).coeffs(1); % 1st non-zero coefficient (m(1) - m(2) = 1)
        % R2 = obj.R(mOrder).coeffs(2); % 2nd non-zero coefficient (m(1) - m(2) = -1)
        R2 = conj(R1); % the 2nd non-zero coefficient is the complex conjugate of the 1st one
        Omega = Omega + 1i/2 * (R2 - R1) * rho.^(mOrder - 1);
    end
end

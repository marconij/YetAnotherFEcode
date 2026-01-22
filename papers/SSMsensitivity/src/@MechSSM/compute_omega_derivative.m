function dOmega = compute_omega_derivative(obj, rho, drho)
    % Compute the derivative of the response frequency Omega at the given 
    % reduced amplitude rho.
    % Inputs:
    %   rho: the reduced amplitude.
    %   drho: the derivative of the reduced amplitude.
    % Outputs:
    %   dOmega: the derivative of the response frequency.

    % Initialize the derivative of the response frequency
    dOmega = zeros(obj.sys.ndv, 1);

    % Compute the partial derivative of Omega wrt rho
    pOmega_rho = 0;
    for mOrder = 3:2:obj.maxOrder
        pOmega_rho = pOmega_rho + 1i/2 * (obj.R(mOrder).coeffs(2) - obj.R(mOrder).coeffs(1)) * (mOrder - 1) * rho.^(mOrder - 2);
    end

    % Loop over design variables
    for dv = 1:obj.sys.ndv
        % Initialize the derivative of the response frequency with the first order term
        dOmega(dv) = 1i/2 * (conj(obj.sys.dlambda(dv)) - obj.sys.dlambda(dv)) * ones(size(rho));

        % Loop over odd orders
        for mOrder = 3:2:obj.maxOrder
            dOmega(dv) = dOmega(dv) + 1i/2 * (obj.R(mOrder).dv(dv).coeffs(2) - obj.R(mOrder).dv(dv).coeffs(1)) * rho.^(mOrder - 1);
        end

        % Add the contribution of the partial derivative of Omega wrt rho
        dOmega(dv) = dOmega(dv) + pOmega_rho .* drho(dv);
    end
end

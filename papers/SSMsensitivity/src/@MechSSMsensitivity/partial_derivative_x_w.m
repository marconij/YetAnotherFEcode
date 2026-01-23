function px_w = partial_derivative_x_w(obj, rho, xk, x, outDof)
    % Compute the partial derivatives of the physical amplitude x with
    % respect to the manifold coefficients w.
    % Inputs:
    %   rho: the reduced amplitude.
    %   xk: the physical amplitudes at all the time steps.
    %   x: the rms of the physical amplitude.
    %   outDof: the index of the output dof.
    % Outputs:
    %   px_w: the partial derivatives of x with respect to w.

    % Initialize partial derivatives
    px_w = repmat(struct('coeffs', []), 1, obj.ssm.maxOrder);

    % % Loop over orders
    % for mOrder = 1:obj.ssm.maxOrder
    %     % Number of multi-indices at the current order
    %     nMIs = obj.ssm.MIs(mOrder).n;
    % 
    %     % Initialize the partial derivative of x with respect to wm
    %     px_w(mOrder).coeffs = zeros(obj.ssm.sys.n, nMIs);
    % 
    %     % Loop over the multi-indices at the current order
    %     for mIdx = 1:nMIs
    %         % Current multi-index
    %         m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);
    % 
    %         % Loop over time steps
    %         for k = 1:obj.ssm.nTheta
    %             % Reduced coordinates
    %             p1Tilde = exp(1i * obj.ssm.theta(k));
    %             p2Tilde = exp(-1i * obj.ssm.theta(k));
    % 
    %             % Partial derivative of xk with respect to wm
    %             px_w(mOrder).coeffs(outDof, mIdx) = px_w(mOrder).coeffs(outDof, mIdx) ...
    %                 + real(p1Tilde.^m(1) .* p2Tilde.^m(2)) * rho.^mOrder * xk(1, k);
    %         end
    % 
    %         % Divide
    %         px_w(mOrder).coeffs(outDof, mIdx) = px_w(mOrder).coeffs(outDof, mIdx) / (obj.ssm.nTimeSteps * x);
    %     end
    % end

    % Loop over orders
    for mOrder = 2:obj.ssm.maxOrder
        % Number of multi-indices at the current order
        nMIs = obj.ssm.MIs(mOrder).n;

        % Initialize the partial derivative of x with respect to wm
        px_w(mOrder).coeffs = zeros(obj.ssm.sys.n, nMIs);

        % Loop over the multi-indices at the current order
        for mIdx = 1:nMIs
            % Current multi-index
            m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);

            % Reduced coordinates at all the time steps
            p1 = rho * exp( 1i * obj.ssm.theta);
            p2 = rho * exp(-1i * obj.ssm.theta);

            % Multi-index exponent
            pm = p1.^m(1) .* p2.^m(2);

            % Compute the partial derivative by summing the contributions from all the time steps
            px_w(mOrder).coeffs(outDof, mIdx) = sum(pm .* xk) / (obj.ssm.nTimeSteps * x);
        end
    end
end

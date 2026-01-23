function v_pdwm_phi = partial_derivative_dw_phi(obj, mOrder, mIdx, pR_phi, v)
    % Compute the partial derivative of the manifold velocity coefficients
    % dwm with respect to the mode shape phi.
    % The derivative is premultiplied by the vector v to avoid storing the
    % full matrix.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the index of the multi-index in the multi-index array.
    %   pR_phi: the partial derivative of the reduced dynamics coefficients
    %       with respect to phi.
    %   v: the vector to premultiply the derivative by.
    % Outputs:
    %   v_pdw_phi: the partial derivative of dwm with respect to phi
    %       premultiplied by v.

    % Initialize with the partial derivative of the V vector
    v_pdwm_phi = obj.partial_derivative_V_phi(mOrder, mIdx, pR_phi, v);

    % Current multi-index
    m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);

    % Add the contribution of the reduced dynamics
    if m(1) - m(2) == 1
        v_pdwm_phi = v_pdwm_phi + v.' * obj.ssm.R(mOrder).coeffs(1) + ...
            (v.' * obj.ssm.sys.phi) * pR_phi(mOrder).coeffs(1, :);
    elseif m(1) - m(2) == -1
        v_pdwm_phi = v_pdwm_phi + v.' * obj.ssm.R(mOrder).coeffs(2) + ...
            (v.' * obj.ssm.sys.phi) * pR_phi(mOrder).coeffs(2, :);
    end
end

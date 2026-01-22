function v_pdw_w = partial_derivative_dw_w(obj, qOrder, qIdx, mOrder, mIdx, pR_w, v)
    % Compute the partial derivative of the manifold velocity coefficients
    % dw at the multi-index q
    % with respect to the manifold coefficients w at the multi-index m.
    % The derivative is premultiplied by the vector v to avoid storing the
    % full matrix.
    % Inputs:
    %   qOrder: the order of the multi-index q.
    %   qIdx: the index of the multi-index q in the multi-index array.
    %   mOrder: the order of the multi-index m.
    %   mIdx: the index of the multi-index m in the multi-index array.
    %   pR_w: the partial derivative of the reduced dynamics coefficients
    %       with respect to wm.
    %   v: the vector to premultiply the derivative by.
    % Outputs:
    %   v_pdw_phi: the partial derivative of dwm with respect to phi
    %       premultiplied by v.

    % Initialize with the partial derivative of the V vector
    v_pdw_w = obj.partial_derivative_V_w(qOrder, qIdx, mOrder, mIdx, pR_w, v);

    % Current multi-indices
    q = obj.ssm.MIs(qOrder).coeffs(:, qIdx);
    m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);

    % Contribution of the first term (LambdaQ * wq)
    if all(q == m)
        v_pdw_w = v_pdw_w + v.' * obj.ssm.LambdaM(qOrder).coeffs(qIdx);
    end

    % Add the contribution of the reduced dynamics
    if q(1) - q(2) == 1
        v_pdw_w = v_pdw_w + (v.' * obj.ssm.sys.phi) * pR_w(qOrder).coeffs(1, :);
    elseif q(1) - q(2) == -1
        v_pdw_w = v_pdw_w + (v.' * obj.ssm.sys.phi) * pR_w(qOrder).coeffs(2, :);
    end
end

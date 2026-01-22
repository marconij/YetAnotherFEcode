function v_pY_w = partial_derivative_Y_w(obj, qOrder, qIdx, mOrder, mIdx, pR_w, v)
    % Compute the partial derivative of the Y vector at the multi-index q
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
    %   v_pY_w: the partial derivative of Yq with respect to wm
    %       premultiplied by v.

    % Initialize
    v_pY_w = zeros(1, obj.ssm.sys.n);

    % Extract sets of multi-indices
    coeffsVm = obj.ssm.MIs(qOrder).coeffsV(qIdx);
    jSet = coeffsVm.j;
    uSet = coeffsVm.u;
    kSet = coeffsVm.k;

    % Compute the V vector at the multi-index m
    for i = 1:coeffsVm.n
        % Extract the multi-indices
        j = jSet(i);
        u = uSet(:, i);
        k = kSet(:, i);

        % Compute the order of the multi-indices
        uOrder = sum(u);
        kOrder = sum(k);

        % Compute the index positions
        uIdx = from_MI_to_position_idx(u);
        % kIdx = from_MI_to_position_idx(k);

        % Extract dw coefficients and partial derivatives
        dwu = obj.ssm.dw(uOrder).coeffs(:, uIdx);
        v_pdwu_wm = obj.partial_derivative_dw_w(uOrder, uIdx, mOrder, mIdx, pR_w, v);

        % Extract R coefficients and partial derivatives
        Rkj = obj.ssm.R(kOrder).coeffs(j);
        pRkj_wm = pR_w(kOrder).coeffs(j, :);

        % Compute the V vector at the multi-index m
        v_pY_w = v_pY_w + u(j) * ((v.' * dwu) * pRkj_wm + Rkj * v_pdwu_wm);
    end
end

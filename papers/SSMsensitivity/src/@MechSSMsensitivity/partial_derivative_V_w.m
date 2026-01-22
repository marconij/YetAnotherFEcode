function v_pV_w = partial_derivative_V_w(obj, qOrder, qIdx, mOrder, mIdx, pR_w, v)
    % Compute the partial derivative of the V vector at the multi-index q
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
    %   v_pV_w: the partial derivative of Vq with respect to wm
    %       premultiplied by v.

    % Initialize
    v_pV_w = zeros(1, obj.ssm.sys.n);

    % Multi-index m
    m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);

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

        % Multi-index u
        u = obj.ssm.MIs(uOrder).coeffs(:, uIdx);

        % Extract w coefficients
        wu = obj.ssm.w(uOrder).coeffs(:, uIdx);

        % Extract R coefficients
        Rkj = obj.ssm.R(kOrder).coeffs(j);
        pRkj_wm = pR_w(kOrder).coeffs(j, :);

        % Compute the V vector at the multi-index m
        v_pV_w = v_pV_w + u(j) * (v.' * wu) * pRkj_wm;

        % Add contribution if u is equal to m (w_u = w_m)
        if all(u == m)
            v_pV_w = v_pV_w + u(j) * Rkj * v.';
        end
    end
end

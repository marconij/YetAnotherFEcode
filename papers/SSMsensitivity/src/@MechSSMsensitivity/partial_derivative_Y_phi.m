function v_pYm_phi = partial_derivative_Y_phi(obj, mOrder, mIdx, pR_phi, v)
    % Compute the partial derivative of the V vector at the multi-index m
    % with respect to the mode shape phi.
    % The derivative is premultiplied by the vector v to avoid storing the
    % full matrix.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the position index of the multi-index in the multi-index
    %       array.
    %   pR_phi: the partial derivative of the reduced dynamics coefficients
    %       with respect to phi.
    % Outputs:
    %   v_pYm_phi: the partial derivative of Ym with respect to phi
    %       premultiplied by v.

    % Initialize
    v_pYm_phi = zeros(1, obj.ssm.sys.n);

    % Extract sets of multi-indices
    coeffsVm = obj.ssm.MIs(mOrder).coeffsV(mIdx);
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
        v_pdwu_phi = obj.partial_derivative_dw_phi(uOrder, uIdx, pR_phi, v);

        % Extract R coefficients and partial derivatives
        Rkj = obj.ssm.R(kOrder).coeffs(j);
        pRkj_phi = pR_phi(kOrder).coeffs(j, :);

        % Compute the V vector at the multi-index m
        v_pYm_phi = v_pYm_phi + u(j) * ((v.' * dwu) * pRkj_phi + Rkj * v_pdwu_phi);
    end
end

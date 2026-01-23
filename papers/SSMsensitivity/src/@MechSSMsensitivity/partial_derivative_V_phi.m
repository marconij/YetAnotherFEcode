function v_pVm_phi = partial_derivative_V_phi(obj, mOrder, mIdx, pR_phi, v)
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
    %   v_pVm_phi: the partial derivative of Vm with respect to phi
    %       premultiplied by v.

    % Initialize
    v_pVm_phi = zeros(1, obj.ssm.sys.n);

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

        % Extract w coefficients
        wu = obj.ssm.w(uOrder).coeffs(:, uIdx);

        % Extract R coefficients
        pRkj_phi = pR_phi(kOrder).coeffs(j, :);

        % Compute the V vector at the multi-index m
        v_pVm_phi = v_pVm_phi + u(j) * (v.' * wu) * pRkj_phi;
    end
end

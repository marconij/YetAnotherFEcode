function pVm_omega = partial_derivative_V_omega(obj, mOrder, mIdx, pR_omega)
    % Compute the partial derivative of the V vector at the multi-index m
    % with respect to the natural frequency omega.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the position index of the multi-index in the multi-index
    %       array.
    %   pR_omega: the partial derivative of the reduced dynamics
    %       coefficients with respect to omega.
    % Outputs:
    %   pVm_omega: the partial derivative of Vm with respect to omega.

    % Initialize
    pVm_omega = zeros(obj.ssm.sys.n, 1);

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
        pRkj_omega = pR_omega(kOrder).coeffs(j);

        % Compute the V vector at the multi-index m
        pVm_omega = pVm_omega + pRkj_omega * u(j) * wu;
    end
end

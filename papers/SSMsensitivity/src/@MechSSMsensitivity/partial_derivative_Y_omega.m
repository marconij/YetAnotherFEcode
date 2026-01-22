function pYm_omega = partial_derivative_Y_omega(obj, mOrder, mIdx, pR_omega, pdw_omega)
    % Compute the partial derivative of the Y vector at the multi-index m
    % with respect to the natural frequency omega.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the position index of the multi-index in the multi-index
    %       array.
    %   pR_omega: the partial derivative of the reduced dynamics
    %       coefficients with respect to omega.
    %   pdw_omega: the partial derivative of the manifold velocity
    %       coefficients with respect to omega.
    % Outputs:
    %   pYm_omega: the partial derivative of Ym with respect to omega.

    % Initialize
    pYm_omega = zeros(obj.ssm.sys.n, 1);

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
        pdwu_omega = pdw_omega(uOrder).coeffs(:, uIdx);

        % Extract R coefficients and partial derivatives
        Rkj = obj.ssm.R(kOrder).coeffs(j);
        pRkj_omega = pR_omega(kOrder).coeffs(j);

        % Compute the V vector at the multi-index m
        pYm_omega = pYm_omega + pRkj_omega * u(j) * dwu + Rkj * u(j) * pdwu_omega;
    end
end

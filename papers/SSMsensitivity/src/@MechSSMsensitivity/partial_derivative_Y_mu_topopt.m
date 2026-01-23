function pYm_mu = partial_derivative_Y_mu_topopt(obj, mOrder, mIdx, pR_mu, pdw_mu)
    % Compute the partial derivative of the Y vector at the multi-index m
    % with respect to the design variables mu for a topology optimization
    % problem.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the position index of the multi-index in the multi-index
    %       array.
    %   pR_mu: the partial derivative of the reduced dynamics coefficients
    %       with respect to mu.
    %   pdw_mu: the partial derivative of the manifold velocity
    %       coefficients with respect to mu.
    % Outputs:
    %   pYm_mu: the partial derivative of Ym with respect to mu.

    % Initialize
    pYm_mu = zeros(obj.ssm.sys.n, 1);

    % Extract sets of multi-indices
    coeffsVm = obj.ssm.MIs(mOrder).coeffsV(mIdx);
    jSet = coeffsVm.j;
    uSet = coeffsVm.u;
    kSet = coeffsVm.k;

    % Compute the Y vector at the multi-index m
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

        % Extract coefficients
        dwu = obj.ssm.dw(uOrder).coeffs(:, uIdx);
        Rkj = obj.ssm.R(kOrder).coeffs(j);

        % Extract partial derivatives of the coefficients
        pdwu_mu = pdw_mu(uOrder).pdv(:, uIdx);
        pRkj_mu = pR_mu(kOrder).pdv(j);

        % Compute the Y vector at the multi-index m
        pYm_mu = pYm_mu + u(j) * (pRkj_mu * dwu + Rkj * pdwu_mu);
    end
end

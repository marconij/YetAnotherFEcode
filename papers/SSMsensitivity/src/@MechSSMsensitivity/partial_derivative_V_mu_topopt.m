function pVm_mu = partial_derivative_V_mu_topopt(obj, mOrder, mIdx, pR_mu)
    % Compute the partial derivative of the V vector at the multi-index m
    % with respect to the design variables mu for a topology optimization
    % problem.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the position index of the multi-index in the multi-index
    %       array.
    %   pR_mu: the partial derivative of the reduced dynamics coefficients
    %       with respect to mu.
    % Outputs:
    %   pVm_mu: the partial derivative of Vm with respect to mu.

    % Initialize
    pVm_mu = zeros(obj.ssm.sys.n, 1);

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
        pRkj_mu = pR_mu(kOrder).pdv(j);

        % Compute the V vector at the multi-index m
        pVm_mu = pVm_mu + pRkj_mu * u(j) * wu;
    end
end

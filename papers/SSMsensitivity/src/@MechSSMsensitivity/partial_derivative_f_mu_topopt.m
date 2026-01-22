function pfm_mu = partial_derivative_f_mu_topopt(obj, mOrder, mIdx, iDOFs)
    % Compute the partial derivative of the nonlinear force contribution fm
    % at the multi-index m with respect to the design variables mu for a
    % topology optimization problem.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the index of the multi-index in the multi-index array.
    %   iDOFs: the index of the DOFs of the current design variable.
    % Outputs:
    %   pfm_mu: the partial derivative of fm with respect to mu.

    % Initialize
    pfm_mu = tensor(zeros(size(iDOFs)), length(iDOFs));

    % Check if the quadratic tensor has at least 1 non-zero element
    if nnz(obj.ssm.sys.T2e) > 0
        % Extract multi-indices
        coeffsNlForce2m = obj.ssm.MIs(mOrder).coeffsNlForce2(mIdx);
        uSet = coeffsNlForce2m.u;
        kSet = coeffsNlForce2m.k;

        % Compute the quadratic force term
        for i = 1:coeffsNlForce2m.n
            % Extract the multi-indices
            u = uSet(:, i);
            k = kSet(:, i);

            % Compute the order of the multi-indices
            uOrder = sum(u);
            kOrder = sum(k);

            % Compute the index positions
            uIdx = from_MI_to_position_idx(u);
            kIdx = from_MI_to_position_idx(k);

            % Extract w coefficients
            wuE = obj.wUnc(uOrder).coeffs(iDOFs, uIdx);
            wkE = obj.wUnc(kOrder).coeffs(iDOFs, kIdx);

            % Compute the quadratic force term
            pfm_mu = pfm_mu + ttv(obj.ssm.sys.T2e, {wuE, wkE}, [2, 3]);
        end
    end

    % Check if the cubic tensor has at least 1 non-zero element
    if nnz(obj.ssm.sys.T3e) > 0
        % Compute the cubic force term
        coeffsNlForce3m = obj.ssm.MIs(mOrder).coeffsNlForce3(mIdx);
        uSet = coeffsNlForce3m.u;
        kSet = coeffsNlForce3m.k;
        lSet = coeffsNlForce3m.l;
        
        for i = 1:coeffsNlForce3m.n
            % Extract the multi-indices
            u = uSet(:, i);
            k = kSet(:, i);
            l = lSet(:, i);

            % Compute the order of the multi-indices
            uOrder = sum(u);
            kOrder = sum(k);
            lOrder = sum(l);

            % Compute the index positions
            uIdx = from_MI_to_position_idx(u);
            kIdx = from_MI_to_position_idx(k);
            lIdx = from_MI_to_position_idx(l);

            % Extract w coefficients
            wuE = obj.wUnc(uOrder).coeffs(iDOFs, uIdx);
            wkE = obj.wUnc(kOrder).coeffs(iDOFs, kIdx);
            wlE = obj.wUnc(lOrder).coeffs(iDOFs, lIdx);

            % Compute the cubic force term
            pfm_mu = pfm_mu + ttv(obj.ssm.sys.T3e, {wuE, wkE, wlE}, [2, 3, 4]);
        end
    end

    % Extract vector from tensor
    pfm_mu = pfm_mu.data;
end

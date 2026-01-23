function pfm_mu = partial_derivative_f_mu(obj, mOrder, mIdx)
    % Compute the partial derivative of the nonlinear force contribution fm
    % at the multi-index m
    % with respect to the design variables mu.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the index of the multi-index in the multi-index array.
    % Outputs:
    %   pfm_mu: the partial derivative of fm with respect to mu.

    % Initialize
    pf0_mu = tensor(zeros(obj.ssm.sys.n, 1), obj.ssm.sys.n);
    pfm_mu = repmat(struct('coeffs', pf0_mu), 1, obj.ssm.sys.ndv);

    % Loop over design variables
    for dv = 1:obj.ssm.sys.ndv
        % Check if the quadratic tensor has at least 1 non-zero element
        if nnz(obj.ssm.sys.pT2(dv).coeffs) > 0
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
                wu = obj.ssm.w(uOrder).coeffs(:, uIdx);
                wk = obj.ssm.w(kOrder).coeffs(:, kIdx);

                % Compute the quadratic force term
                pfm_mu(dv).coeffs = pfm_mu(dv).coeffs + ...
                    ttv(obj.ssm.sys.pT2(dv).coeffs, {wu, wk}, [2, 3]);
            end
        end

        % Check if the cubic tensor has at least 1 non-zero element
        if nnz(obj.ssm.sys.pT3(dv).coeffs) > 0
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
                wu = obj.ssm.w(uOrder).coeffs(:, uIdx);
                wk = obj.ssm.w(kOrder).coeffs(:, kIdx);
                wl = obj.ssm.w(lOrder).coeffs(:, lIdx);

                % Compute the cubic force term
                pfm_mu(dv).coeffs = pfm_mu(dv).coeffs + ...
                    ttv(obj.ssm.sys.pT3(dv).coeffs, {wu, wk, wl}, [2, 3, 4]);
            end
        end

        % Extract vector from tensor
        pfm_mu(dv).coeffs = pfm_mu(dv).coeffs.data;
    end
end

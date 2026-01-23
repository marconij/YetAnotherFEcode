function v_pf_w = partial_derivative_f_w(obj, qOrder, qIdx, mOrder, mIdx, v)
    % Compute the partial derivative of the nonlinear force contribution fq
    % at the multi-index q with respect to the manifold coefficient wm at
    % the multi-index m.
    % The derivative is premultiplied by the vector v to avoid storing the
    % full matrix.
    % Inputs:
    %   qOrder: the order of the multi-index q.
    %   qIdx: the index of the multi-index q in the multi-index array.
    %   mOrder: the order of the multi-index m.
    %   mIdx: the index of the multi-index m in the multi-index array.
    %   v: the vector to premultiply the derivative by.
    % Outputs:
    %   v_pf_w: the partial derivative of fq with respect to wm
    %       premultiplied by v.

    % Initialize
    v_pf_w = tensor(zeros(obj.ssm.sys.n, 1), obj.ssm.sys.n);

    % Multi-index m
    m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);

    % Check if the quadratic tensor has at least 1 non-zero element
    if nnz(obj.ssm.sys.T2) > 0
        % Extract multi-indices
        coeffsNlForce2 = obj.ssm.MIs(qOrder).coeffsNlForce2(qIdx);
        uSet = coeffsNlForce2.u;
        kSet = coeffsNlForce2.k;

        % Compute the quadratic force term
        for i = 1:coeffsNlForce2.n
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

            % Check if u is equal to m (w_u = w_m)
            if all(u == m)
                v_pf_w = v_pf_w + ttv(obj.ssm.sys.T2, {v, wk}, [1, 3]);
            end

            % Check if k is equal to m (w_k = w_m)
            if all(k == m)
                v_pf_w = v_pf_w + ttv(obj.ssm.sys.T2, {v, wu}, [1, 2]);
            end
        end
    end

    % Check if the cubic tensor has at least 1 non-zero element
    if nnz(obj.ssm.sys.T3) > 0
        % Compute the cubic force term
        coeffsNlForce3 = obj.ssm.MIs(qOrder).coeffsNlForce3(qIdx);
        uSet = coeffsNlForce3.u;
        kSet = coeffsNlForce3.k;
        lSet = coeffsNlForce3.l;
        
        for i = 1:coeffsNlForce3.n
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

            % Check if u is equal to m (w_u = w_m)
            if all(u == m)
                v_pf_w = v_pf_w + ttv(obj.ssm.sys.T3, {v, wk, wl}, [1, 3, 4]);
            end

            % Check if k is equal to m (w_k = w_m)
            if all(k == m)
                v_pf_w = v_pf_w + ttv(obj.ssm.sys.T3, {v, wu, wl}, [1, 2, 4]);
            end

            % Check if l is equal to m (w_l = w_m)
            if all(l == m)
                v_pf_w = v_pf_w + ttv(obj.ssm.sys.T3, {v, wu, wk}, [1, 2, 3]);
            end
        end
    end

    % Extract vector from tensor and transpose the result
    v_pf_w = v_pf_w.data.';
end

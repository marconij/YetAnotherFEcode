function v_pfm_phi = partial_derivative_f_phi(obj, mOrder, mIdx, v)
    % Compute the partial derivative of the nonlinear force contribution fm
    % at the multi-index m with respect to the mode shape phi.
    % The derivative is premultiplied by the vector v to avoid storing the
    % full matrix.
    % The mode shape phi is equal to the 1st order manifold coefficients
    % (w_10 = w_01 = phi).
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the index of the multi-index in the multi-index array.
    %   v: the vector to premultiply the derivative by.
    % Outputs:
    %   v_pfm_phi: the partial derivative of fm with respect to phi
    %       premultiplied by v.

    % Initialize as a tensor
    v_pfm_phi = tensor(zeros(obj.ssm.sys.n, 1), obj.ssm.sys.n);

    % Check if the quadratic tensor has at least 1 non-zero element
    if nnz(obj.ssm.sys.T2) > 0
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

            % Check if u is equal to 10 or 01 (w_u = phi)
            if all(u == [1; 0]) || all(u == [0; 1])
                v_pfm_phi = v_pfm_phi + ttv(obj.ssm.sys.T2, {v, wk}, [1, 3]);
            end

            % Check if k is equal to 10 or 01 (w_k = phi)
            if all(k == [1; 0]) || all(k == [0; 1])
                v_pfm_phi = v_pfm_phi + ttv(obj.ssm.sys.T2, {v, wu}, [1, 2]);
            end
        end
    end

    % Check if the cubic tensor has at least 1 non-zero element
    if nnz(obj.ssm.sys.T3) > 0
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

            % Check if u is equal to 10 or 01 (w_u = phi)
            if all(u == [1; 0]) || all(u == [0; 1])
                v_pfm_phi = v_pfm_phi + ttv(obj.ssm.sys.T3, {v, wk, wl}, [1, 3, 4]);
            end

            % Check if k is equal to 10 or 01 (w_k = phi)
            if all(k == [1; 0]) || all(k == [0; 1])
                v_pfm_phi = v_pfm_phi + ttv(obj.ssm.sys.T3, {v, wu, wl}, [1, 2, 4]);
            end

            % Check if l is equal to 10 or 01 (w_l = phi)
            if all(l == [1; 0]) || all(l == [0; 1])
                v_pfm_phi = v_pfm_phi + ttv(obj.ssm.sys.T3, {v, wu, wk}, [1, 2, 3]);
            end
        end
    end

    % Extract vector from tensor and transpose the result
    v_pfm_phi = v_pfm_phi.data.';
end

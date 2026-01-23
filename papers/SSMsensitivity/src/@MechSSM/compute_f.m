function [fm, dfm] = compute_f(obj, mOrder, mIdx)
    % Compute the nonlinear force contribution at the multi-index m
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the position index of the multi-index in the multi-index
    %       array.
    % Outputs:
    %   fm: the nonlinear force contribution at the multi-index m.
    %   dfm: the derivative of the nonlinear force contribution at the
    %       multi-index m.

    % Initialize
    fm = tensor(zeros(obj.sys.n, 1), obj.sys.n);
    if obj.computeSensitivity
        dfm = repmat(struct('coeffs', fm), 1, obj.sys.ndv);
    else
        dfm = [];
    end

    % Quadratic contributions
    coeffsNlForce2m = obj.MIs(mOrder).coeffsNlForce2(mIdx);
    uSet = coeffsNlForce2m.u;
    kSet = coeffsNlForce2m.k;
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
        wu = obj.w(uOrder).coeffs(:, uIdx);
        wk = obj.w(kOrder).coeffs(:, kIdx);

        % Compute the quadratic displacement term
        if nnz(obj.sys.T2) > 0
            fm = fm + ttv(obj.sys.T2, {wu, wk}, [2, 3]);
        end

        % Extract dw coefficients
        dwu = obj.dw(uOrder).coeffs(:, uIdx);
        dwk = obj.dw(kOrder).coeffs(:, kIdx);

        % Compute the quadratic velocity term
        if nnz(obj.sys.D2) > 0
            fm = fm + ttv(obj.sys.D2, {dwu, dwk}, [2, 3]);
        end
    end

    % Cubic contributions
    coeffsNlForce3m = obj.MIs(mOrder).coeffsNlForce3(mIdx);
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
        wu = obj.w(uOrder).coeffs(:, uIdx);
        wk = obj.w(kOrder).coeffs(:, kIdx);
        wl = obj.w(lOrder).coeffs(:, lIdx);

        % Compute the cubic displacement term
        if nnz(obj.sys.T3) > 0
            fm = fm + ttv(obj.sys.T3, {wu, wk, wl}, [2, 3, 4]);
        end

        % Extract dw coefficients
        dwu = obj.dw(uOrder).coeffs(:, uIdx);
        dwk = obj.dw(kOrder).coeffs(:, kIdx);
        dwl = obj.dw(lOrder).coeffs(:, lIdx);

        % Compute the cubic velocity term
        if nnz(obj.sys.D3) > 0
            fm = fm + ttv(obj.sys.D3, {dwu, dwk, dwl}, [2, 3, 4]);
        end
    end

    % Extract vector from tensor
    fm = fm.data;

    % Compute derivative
    if obj.computeSensitivity
        % Quadratic contributions
        coeffsNlForce2m = obj.MIs(mOrder).coeffsNlForce2(mIdx);
        uSet = coeffsNlForce2m.u;
        kSet = coeffsNlForce2m.k;
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
            wu = obj.w(uOrder).coeffs(:, uIdx);
            wk = obj.w(kOrder).coeffs(:, kIdx);

            % Extract dw coefficients
            dwu = obj.dw(uOrder).coeffs(:, uIdx);
            dwk = obj.dw(kOrder).coeffs(:, kIdx);

            % Loop over design variables
            for dv = 1:obj.sys.ndv
                % Extract the derivative of the w coefficients
                Dwu = obj.w(uOrder).dv(dv).coeffs(:, uIdx);
                Dwk = obj.w(kOrder).dv(dv).coeffs(:, kIdx);

                % Extract the derivative of the dw coefficients
                Ddwu = obj.dw(uOrder).dv(dv).coeffs(:, uIdx);
                Ddwk = obj.dw(kOrder).dv(dv).coeffs(:, kIdx);

                % Check if the derivative of the quadratic displacement tensor has at least 1 non-zero element
                if nnz(obj.sys.pT2(dv).coeffs) > 0
                    dfm(dv).coeffs = dfm(dv).coeffs + ...
                        ttv(obj.sys.pT2(dv).coeffs, {wu, wk}, [2, 3]);
                end
                
                % Check if the quadratic displacement tensor has at least 1 non-zero element
                if nnz(obj.sys.T2) > 0
                    dfm(dv).coeffs = dfm(dv).coeffs + ...
                        ttv(obj.sys.T2, {Dwu, wk}, [2, 3]) + ...
                        ttv(obj.sys.T2, {wu, Dwk}, [2, 3]);
                end

                % Check if the derivative of the quadratic velocity tensor has at least 1 non-zero element
                % TODO: we assume that this tensor is constant, so this term is zero
                % if nnz(obj.sys.pD2(dv).coeffs) > 0
                %     dfm(dv).coeffs = dfm(dv).coeffs + ...
                %         ttv(obj.sys.pD2(dv).coeffs, {dwu, dwk}, [2, 3]);
                % end

                % Check if the quadratic velocity tensor has at least 1 non-zero element
                if nnz(obj.sys.D2) > 0
                    dfm(dv).coeffs = dfm(dv).coeffs + ...
                        ttv(obj.sys.D2, {Ddwu, dwk}, [2, 3]) + ...
                        ttv(obj.sys.D2, {dwu, Ddwk}, [2, 3]);
                end
            end
        end

        % Cubic contributions
        coeffsNlForce3m = obj.MIs(mOrder).coeffsNlForce3(mIdx);
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
            wu = obj.w(uOrder).coeffs(:, uIdx);
            wk = obj.w(kOrder).coeffs(:, kIdx);
            wl = obj.w(lOrder).coeffs(:, lIdx);

            % Extract dw coefficients
            dwu = obj.dw(uOrder).coeffs(:, uIdx);
            dwk = obj.dw(kOrder).coeffs(:, kIdx);
            dwl = obj.dw(lOrder).coeffs(:, lIdx);

            % Loop over design variables
            for dv = 1:obj.sys.ndv
                % Extract the derivative of the w coefficients
                Dwu = obj.w(uOrder).dv(dv).coeffs(:, uIdx);
                Dwk = obj.w(kOrder).dv(dv).coeffs(:, kIdx);
                Dwl = obj.w(lOrder).dv(dv).coeffs(:, lIdx);

                % Extract the derivative of the dw coefficients
                Ddwu = obj.dw(uOrder).dv(dv).coeffs(:, uIdx);
                Ddwk = obj.dw(kOrder).dv(dv).coeffs(:, kIdx);
                Ddwl = obj.dw(lOrder).dv(dv).coeffs(:, lIdx);

                % Check if the derivative of the cubic displacement tensor has at least 1 non-zero element
                if nnz(obj.sys.pT3(dv).coeffs) > 0
                    dfm(dv).coeffs = dfm(dv).coeffs + ...
                        ttv(obj.sys.pT3(dv).coeffs, {wu, wk, wl}, [2, 3, 4]);
                end

                % Check if the cubic displacement tensor has at least 1 non-zero element
                if nnz(obj.sys.T3) > 0
                    dfm(dv).coeffs = dfm(dv).coeffs + ...
                        ttv(obj.sys.T3, {Dwu, wk, wl}, [2, 3, 4]) + ...
                        ttv(obj.sys.T3, {wu, Dwk, wl}, [2, 3, 4]) + ...
                        ttv(obj.sys.T3, {wu, wk, Dwl}, [2, 3, 4]);
                end

                % Check if the derivative of the cubic velocity tensor has at least 1 non-zero element
                % TODO: we assume that this tensor is constant, so this term is zero
                % if nnz(obj.sys.pD3(dv).coeffs) > 0
                %     dfm(dv).coeffs = dfm(dv).coeffs + ...
                %         ttv(obj.sys.pD3(dv).coeffs, {dwu, dwk, dwl}, [2, 3, 4]);
                % end

                % Check if the cubic velocity tensor has at least 1 non-zero element
                if nnz(obj.sys.D3) > 0
                    dfm(dv).coeffs = dfm(dv).coeffs + ...
                        ttv(obj.sys.D3, {Ddwu, dwk, dwl}, [2, 3, 4]) + ...
                        ttv(obj.sys.D3, {dwu, Ddwk, dwl}, [2, 3, 4]) + ...
                        ttv(obj.sys.D3, {dwu, dwk, Ddwl}, [2, 3, 4]);
                end
            end
        end

        % Loop over design variables
        for dv = 1:obj.sys.ndv
            % Extract vector from tensor
            dfm(dv).coeffs = dfm(dv).coeffs.data;
        end
    end
end

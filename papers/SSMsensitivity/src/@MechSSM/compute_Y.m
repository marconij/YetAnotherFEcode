function [Ym, dYm] = compute_Y(obj, mOrder, mIdx)
    % Compute the Y vector at the multi-index m.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the position index of the multi-index in the multi-index
    %       array.
    % Outputs:
    %   Ym: the Y vector at the multi-index m.
    %   dYm: the derivative of the Y vector at the multi-index m.

    % Initialize
    Ym = zeros(obj.sys.n, 1);
    if obj.computeSensitivity
        dYm = repmat(struct('coeffs', zeros(obj.sys.n, 1)), 1, obj.sys.ndv);
    else
        dYm = [];
    end

    % Extract sets of multi-indices
    coeffsVm = obj.MIs(mOrder).coeffsV(mIdx);
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

        % Extract dw coefficients
        dwu = obj.dw(uOrder).coeffs(:, uIdx);

        % Extract R coefficients
        Rkj = obj.R(kOrder).coeffs(j);

        % Compute the Y vector at the multi-index m
        Ym = Ym + Rkj * u(j) * dwu;

        % Loop over design variables
        if obj.computeSensitivity
            for dv = 1:obj.sys.ndv
                % Extract the derivatives
                ddwu = obj.dw(uOrder).dv(dv).coeffs(:, uIdx);
                dRkj = obj.R(kOrder).dv(dv).coeffs(j);

                % Compute the derivative of the Y vector at the multi-index m
                dYm(dv).coeffs = dYm(dv).coeffs + dRkj * u(j) * dwu + Rkj * u(j) * ddwu;
            end
        end
    end
end

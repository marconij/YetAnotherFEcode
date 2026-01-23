function [Vm, dVm] = compute_V(obj, mOrder, mIdx)
    % Compute the V vector at the multi-index m.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the position index of the multi-index in the multi-index
    %       array.
    % Outputs:
    %   Vm: the V vector at the multi-index m.
    %   dVm: the derivative of the V vector at the multi-index m.

    % Initialize
    Vm = zeros(obj.sys.n, 1);
    if obj.computeSensitivity
        dVm = repmat(struct('coeffs', zeros(obj.sys.n, 1)), 1, obj.sys.ndv);
    else
        dVm = [];
    end

    % Extract sets of multi-indices
    coeffsVm = obj.MIs(mOrder).coeffsV(mIdx);
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
        wu = obj.w(uOrder).coeffs(:, uIdx);

        % Extract R coefficients
        Rkj = obj.R(kOrder).coeffs(j);

        % Compute the V vector at the multi-index m
        Vm = Vm + Rkj * u(j) * wu;

        % Loop over design variables
        if obj.computeSensitivity
            for dv = 1:obj.sys.ndv
                % Extract the derivatives
                dwu = obj.w(uOrder).dv(dv).coeffs(:, uIdx);
                dRkj = obj.R(kOrder).dv(dv).coeffs(j);

                % Compute the derivative of the V vector at the multi-index m
                dVm(dv).coeffs = dVm(dv).coeffs + dRkj * u(j) * wu + Rkj * u(j) * dwu;
            end
        end
    end
end

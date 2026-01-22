function [dwm, ddwm] = compute_dw(obj, mOrder, mIdx, LambdaM, Vm, dLambdaM, dVm)
    % Compute the manifold velocity coefficients at the multi-index m.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the position index of the multi-index in the multi-index array.
    %   LambdaM: the eigenvalue coefficient at the multi-index m.
    %   Vm: the V vector at the multi-index m.
    %   dLambdaM: the derivative of the eigenvalue coefficient at the
    %       multi-index m.
    %   dVm: the derivative of the V vector at the multi-index m.
    % Outputs:
    %   dwm: the manifold velocity coefficients at the multi-index m.
    %   ddwm: the derivative of the manifold velocity coefficients at the
    %       multi-index m.

    % Multi-index
    m = obj.MIs(mOrder).coeffs(:, mIdx);

    % Compute the manifold velocity coefficients at the multi-index m
    dwm = LambdaM * obj.w(mOrder).coeffs(:, mIdx) + Vm;
    if m(1) - m(2) == 1
        dwm = dwm + obj.sys.phi * obj.R(mOrder).coeffs(1);
    end

    % Compute the derivative of the manifold velocity coefficients at the multi-index m
    if obj.computeSensitivity
        ddwm = repmat(struct('coeffs', zeros(obj.sys.n, 1)), 1, obj.sys.ndv);
        for dv = 1:obj.sys.ndv
            ddwm(dv).coeffs = dLambdaM(dv).coeffs * obj.w(mOrder).coeffs(:, mIdx) + ...
                LambdaM * obj.w(mOrder).dv(dv).coeffs(:, mIdx) + dVm(dv).coeffs;
            if m(1) - m(2) == 1
                ddwm(dv).coeffs = ddwm(dv).coeffs + obj.sys.dphi(dv).coeffs * obj.R(mOrder).coeffs(1) + ...
                    obj.sys.phi * obj.R(mOrder).dv(dv).coeffs(1);
            end
        end
    else
        ddwm = [];
    end
end

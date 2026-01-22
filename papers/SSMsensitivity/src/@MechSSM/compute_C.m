function [Cm, dCm] = compute_C(obj, mOrder, mIdx, LambdaM, Vm, Ym, fm, dLambdaM, dVm, dYm, dfm)
    % Compute the C vector at the multi-index m.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the position index of the multi-index in the multi-index
    %       array.
    %   LambdaM: the eigenvalue coefficient at the multi-index m.
    %   Vm: the V vector at the multi-index m.
    %   Ym: the Y vector at the multi-index m.
    %   fm: the nonlinear force contribution at the multi-index m.
    %   dLambdaM: the derivative of the eigenvalue coefficient at the
    %       multi-index m.
    %   dVm: the derivative of the V vector at the multi-index m.
    %   dYm: the derivative of the Y vector at the multi-index m.
    %   dfm: the derivative of the nonlinear force contribution at the
    %       multi-index m.
    % Outputs:
    %   Cm: the C vector at the multi-index m.
    %   dCm: the derivative of the C vector at the multi-index m.

    % Compute the C vector at the multi-index m
    Cm = -(LambdaM * obj.sys.M + obj.sys.C) * Vm - obj.sys.M * Ym - fm;

    % Compute derivative
    if obj.computeSensitivity
        % Initialize
        dCm = repmat(struct('coeffs', zeros(obj.sys.n, 1)), 1, obj.sys.ndv);

        % Loop over design variables
        for dv = 1:obj.sys.ndv
            % Compute the derivative of the C vector at the multi-index m
            dCm(dv).coeffs = -(dLambdaM(dv).coeffs * obj.sys.M + LambdaM * obj.sys.pM(dv).coeffs + obj.sys.dC(dv).coeffs) * Vm - ...
                (LambdaM * obj.sys.M + obj.sys.C) * dVm(dv).coeffs - ...
                obj.sys.pM(dv).coeffs * Ym - obj.sys.M * dYm(dv).coeffs - dfm(dv).coeffs;
        end
    else
        dCm = [];
    end
end

function [Rm, dRm] = compute_R(obj, mOrder, mIdx, LambdaM, Cm, dLambdaM, dCm)
    % Compute the reduced dynamics coefficients at the multi-index m.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the position index of the multi-index in the multi-index
    %       array.
    %   LambdaM: the eigenvalue coefficient at the multi-index m.
    %   Cm: the C vector at the multi-index m.
    %   dLambdaM: the derivative of the eigenvalue coefficient at the
    %       multi-index m.
    %   dCm: the derivative of the C vector at the multi-index m.
    % Outputs:
    %   Rm: the reduced dynamics coefficients at the multi-index m.
    %   dRm: the derivative of the reduced dynamics coefficients at the
    %       multi-index m.

    % Compute Rm for the multi-index m for which m(1) - m(2) = 1
    num = dot(obj.sys.phi, Cm);
    den = LambdaM + obj.sys.lambda + obj.sys.delta;
    Rm = num / den;

    % Compute the derivative of Rm
    if obj.computeSensitivity
        dRm = repmat(struct('coeffs', 0), 1, obj.sys.ndv);
        for dv = 1:obj.sys.ndv
            dnum = dot(obj.sys.dphi(dv).coeffs, Cm) + dot(obj.sys.phi, dCm(dv).coeffs);
            dden = dLambdaM(dv).coeffs + obj.sys.dlambda(dv) + obj.sys.ddelta(dv);
            dRm(dv).coeffs = (dnum * den - num * dden) / den^2;
        end
    else
        dRm = [];
    end
end

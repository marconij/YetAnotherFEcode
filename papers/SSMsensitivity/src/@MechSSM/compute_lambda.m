function [LambdaM, dLambdaM] = compute_lambda(obj, mOrder, mIdx)
    % Compute the eigenvalue coefficient at the multi-index m.
    % Inputs:
    %   mOrder: the order of the multi-index.
    %   mIdx: the position index of the multi-index in the multi-index
    %       array.
    % Outputs:
    %   LambdaM: the eigenvalue coefficient at the multi-index m.
    %   dLambdaM: the derivative of the eigenvalue coefficient at the 
    %       multi-index m.

    % Extract current multi-index
    m = obj.MIs(mOrder).coeffs(:, mIdx);
    
    % Compute eigenvalue coefficient
    LambdaM = dot(m, obj.sys.Lambda);

    % Compute derivative
    if obj.computeSensitivity
        dLambdaM = repmat(struct('coeffs', []), 1, obj.sys.ndv);
        for dv = 1:obj.sys.ndv
            dLambdaM(dv).coeffs = dot(m, obj.sys.dLambda(dv, :));
        end
    else
        dLambdaM = [];
    end
end

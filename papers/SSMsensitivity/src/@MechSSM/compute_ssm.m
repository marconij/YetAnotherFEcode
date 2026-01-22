function compute_ssm(obj, maxOrder, varargin)
    % Compute the SSM parametrization coefficients.
    % Additionally, the SSM sensitivity can be computed using the direct
    % differentiation by setting the `computeSensitivity` flag to true.
    % Inputs:
    %   maxOrder: the maximum expansion order.
    %   storeCoefficients: flag to store the SSM coefficients (optional,
    %       default is false).
    %   computeSensitivity: flag to computate the SSM sensitivity using the
    %       direct differentiation (optional, default is false).

    % Parse the inputs
    p = inputParser;
    addOptional(p, 'storeCoefficients', false);
    addOptional(p, 'computeSensitivity', false);
    parse(p, varargin{:});
    obj.computeSensitivity = p.Results.computeSensitivity;
    obj.storeCoefficients = p.Results.storeCoefficients;

    % If the `computeSensitivity` flag is set to true, check if the modal
    % analysis sensitivity is available
    if obj.computeSensitivity
        if isempty(obj.sys.dLambda)
            error(['The modal analysis sensitivity is not available. ', ...
                'Call the `modal_analysis_sensitivity` method of the `MechSystem` object.']);
        end
    end

    % Set the expansion order
    obj.maxOrder = maxOrder;

    % Initialize the multi-index set
    obj.MIs = initialize_multi_index(maxOrder, 'MIs', obj.MIs);

    % Initialize coefficients
    obj.R = repmat(struct('coeffs', [], 'dv', []), 1, maxOrder);
    obj.w = repmat(struct('coeffs', [], 'dv', []), 1, maxOrder);
    obj.dw = repmat(struct('coeffs', [], 'dv', []), 1, maxOrder);

    % Initialize auxiliary coefficients
    if obj.storeCoefficients
        obj.LambdaM = repmat(struct('coeffs', [], 'dv', []), 1, maxOrder);
        obj.f = repmat(struct('coeffs', [], 'dv', []), 1, maxOrder);
        obj.V = repmat(struct('coeffs', [], 'dv', []), 1, maxOrder);
        obj.Y = repmat(struct('coeffs', [], 'dv', []), 1, maxOrder);
        obj.C = repmat(struct('coeffs', [], 'dv', []), 1, maxOrder);
        obj.h = repmat(struct('coeffs', [], 'dv', []), 1, maxOrder);
    end

    % 1st order coefficients
    obj.R(1).coeffs = obj.sys.Lambda(:).';
    obj.w(1).coeffs = obj.sys.Phi;
    obj.dw(1).coeffs = obj.sys.Lambda.' .* obj.sys.Phi;

    % 1st order sensitivity
    if obj.computeSensitivity
        % Initialize
        obj.R(1).dv = repmat(struct('coeffs', []), 1, obj.sys.ndv);
        obj.w(1).dv = repmat(struct('coeffs', []), 1, obj.sys.ndv);
        obj.dw(1).dv = repmat(struct('coeffs', []), 1, obj.sys.ndv);

        % Loop over design variables
        for dv = 1:obj.sys.ndv
            obj.R(1).dv(dv).coeffs = obj.sys.dLambda(dv, :);
            obj.w(1).dv(dv).coeffs = obj.sys.dPhi(dv).coeffs;
            obj.dw(1).dv(dv).coeffs = obj.sys.dLambda(dv, :) .* obj.sys.Phi + ...
                obj.sys.Lambda.' .* obj.sys.dPhi(dv).coeffs;
        end
    end

    % Loop over orders
    for mOrder = 2:obj.maxOrder
        % Number of multi-indices at this order
        nMIs = size(obj.MIs(mOrder).coeffs, 2);

        % Initialize coefficients
        obj.R(mOrder).coeffs = zeros(1, 2);
        obj.w(mOrder).coeffs = zeros(obj.sys.n, nMIs);
        obj.dw(mOrder).coeffs = zeros(obj.sys.n, nMIs);

        % Initialize auxiliary coefficients
        if obj.storeCoefficients
            obj.LambdaM(mOrder).coeffs = zeros(1, nMIs);
            obj.f(mOrder).coeffs = zeros(obj.sys.n, nMIs);
            obj.V(mOrder).coeffs = zeros(obj.sys.n, nMIs);
            obj.Y(mOrder).coeffs = zeros(obj.sys.n, nMIs);
            obj.C(mOrder).coeffs = zeros(obj.sys.n, nMIs);
            obj.h(mOrder).coeffs = zeros(obj.sys.n, nMIs);
        end

        % Initialize sensitivity
        if obj.computeSensitivity
            obj.R(mOrder).dv = repmat(struct('coeffs', zeros(1, 2)), 1, obj.sys.ndv);
            obj.w(mOrder).dv = repmat(struct('coeffs', zeros(obj.sys.n, nMIs)), 1, obj.sys.ndv);
            obj.dw(mOrder).dv = repmat(struct('coeffs', zeros(obj.sys.n, nMIs)), 1, obj.sys.ndv);

            % Initialize auxiliary sensitivity
            if obj.storeCoefficients
                obj.LambdaM(mOrder).dv = repmat(struct('coeffs', zeros(1, nMIs)), 1, obj.sys.ndv);
                obj.f(mOrder).dv = repmat(struct('coeffs', zeros(obj.sys.n, nMIs)), 1, obj.sys.ndv);
                obj.V(mOrder).dv = repmat(struct('coeffs', zeros(obj.sys.n, nMIs)), 1, obj.sys.ndv);
                obj.Y(mOrder).dv = repmat(struct('coeffs', zeros(obj.sys.n, nMIs)), 1, obj.sys.ndv);
                obj.C(mOrder).dv = repmat(struct('coeffs', zeros(obj.sys.n, nMIs)), 1, obj.sys.ndv);
                obj.h(mOrder).dv = repmat(struct('coeffs', zeros(obj.sys.n, nMIs)), 1, obj.sys.ndv);
            end
        end

        % For real systems, symmetric multi-indices have complex
        % conjugate coefficients. Therefore, we only need to solve
        % the cohomological equation for half of the multi-indices
        % and then take the complex conjugate of the coefficients
        % for the other half.
        nUniqueMIs = ceil(nMIs / 2);

        % Loop over the first half of the multi-indices
        for mIdx = 1:nUniqueMIs
            % Current multi-index
            m = obj.MIs(mOrder).coeffs(:, mIdx);

            % Eigenvalue coefficient
            [LambdaM, dLambdaM] = obj.compute_lambda(mOrder, mIdx);
            if obj.storeCoefficients
                obj.LambdaM(mOrder).coeffs(mIdx) = LambdaM;
                if obj.computeSensitivity
                    for dv = 1:obj.sys.ndv
                        obj.LambdaM(mOrder).dv(dv).coeffs(mIdx) = dLambdaM(dv).coeffs;
                    end
                end
            end

            % Compute the nonlinear force contribution at the multi-index m
            [fm, dfm] = obj.compute_f(mOrder, mIdx);
            if obj.storeCoefficients
                obj.f(mOrder).coeffs(:, mIdx) = fm;
                if obj.computeSensitivity
                    for dv = 1:obj.sys.ndv
                        obj.f(mOrder).dv(dv).coeffs(:, mIdx) = dfm(dv).coeffs;
                    end
                end
            end

            % Compute the V vector at the multi-index m
            [Vm, dVm] = obj.compute_V(mOrder, mIdx);
            if obj.storeCoefficients
                obj.V(mOrder).coeffs(:, mIdx) = Vm;
                if obj.computeSensitivity
                    for dv = 1:obj.sys.ndv
                        obj.V(mOrder).dv(dv).coeffs(:, mIdx) = dVm(dv).coeffs;
                    end
                end
            end

            % Compute the Y vector at the multi-index m
            [Ym, dYm] = obj.compute_Y(mOrder, mIdx);
            if obj.storeCoefficients
                obj.Y(mOrder).coeffs(:, mIdx) = Ym;
                if obj.computeSensitivity
                    for dv = 1:obj.sys.ndv
                        obj.Y(mOrder).dv(dv).coeffs(:, mIdx) = dYm(dv).coeffs;
                    end
                end
            end
            
            % Compute the C vector at the multi-index m
            [Cm, dCm] = obj.compute_C(mOrder, mIdx, LambdaM, Vm, Ym, fm, dLambdaM, dVm, dYm, dfm);
            if obj.storeCoefficients
                obj.C(mOrder).coeffs(:, mIdx) = Cm;
                if obj.computeSensitivity
                    for dv = 1:obj.sys.ndv
                        obj.C(mOrder).dv(dv).coeffs(:, mIdx) = dCm(dv).coeffs;
                    end
                end
            end

            % Compute the normal-form style parametrization coefficients.
            % The coefficients are non-zero only if the near-resonance
            % condition is satisfied, i.e., m(1) - m(2) = 1 or m(1) - m(2) = -1.
            % If m(1) - m(2) = 1, then only the first element of the
            % coefficients is non-zero. If m(1) - m(2) = -1, then only
            % the second element of the coefficients is non-zero.
            % Here, only the first case is considered.
            if m(1) - m(2) == 1
                [Rm, dRm] = obj.compute_R(mOrder, mIdx, LambdaM, Cm, dLambdaM, dCm);
                obj.R(mOrder).coeffs(1) = Rm;
                if obj.computeSensitivity
                    for dv = 1:obj.sys.ndv
                        obj.R(mOrder).dv(dv).coeffs(1) = dRm(dv).coeffs;
                    end
                end
            end

            % Compute hm = Dm * Rm + Cm
            % The coefficients for R is non-zero only for the multi-indices
            % such that m(1) - m(2) = 1 or m(1) - m(2) = -1.
            % Here, only the first case is considered.
            hm = Cm;
            dhm = dCm;
            if m(1) - m(2) == 1
                Dm = -((LambdaM + obj.sys.lambda) * obj.sys.M + obj.sys.C) * obj.sys.phi;
                hm = hm + Dm * obj.R(mOrder).coeffs(1);

                % Sensitivity
                if obj.computeSensitivity
                    % Loop over design variables
                    for dv = 1:obj.sys.ndv
                        dDm_dv = -((dLambdaM(dv).coeffs + obj.sys.dlambda(dv)) * obj.sys.M + ...
                            (LambdaM + obj.sys.lambda) * obj.sys.pM(dv).coeffs + obj.sys.dC(dv).coeffs) * obj.sys.phi - ...
                            ((LambdaM + obj.sys.lambda) * obj.sys.M + obj.sys.C) * obj.sys.dphi(dv).coeffs;
                        dhm(dv).coeffs = dhm(dv).coeffs + dDm_dv * obj.R(mOrder).coeffs(1) + ...
                            Dm * obj.R(mOrder).dv(dv).coeffs(1);
                    end
                end
            end
            if obj.storeCoefficients
                obj.h(mOrder).coeffs(:, mIdx) = hm;
                if obj.computeSensitivity
                    for dv = 1:obj.sys.ndv
                        obj.h(mOrder).dv(dv).coeffs(:, mIdx) = dhm(dv).coeffs;
                    end
                end
            end

            % Cohomological matrix
            Lm = obj.sys.K + LambdaM * obj.sys.C + LambdaM^2 * obj.sys.M;

            % Solve the cohomological equation
            obj.w(mOrder).coeffs(:, mIdx) = lsqminnorm(Lm, hm);
            if obj.computeSensitivity
                for dv = 1:obj.sys.ndv
                    rhs = dhm(dv).coeffs - (obj.sys.pK(dv).coeffs + dLambdaM(dv).coeffs * obj.sys.C + ...
                        LambdaM * obj.sys.dC(dv).coeffs + 2 * LambdaM * dLambdaM(dv).coeffs * obj.sys.M + ...
                        LambdaM^2 * obj.sys.pM(dv).coeffs) * obj.w(mOrder).coeffs(:, mIdx);
                    obj.w(mOrder).dv(dv).coeffs(:, mIdx) = lsqminnorm(Lm, rhs);
                end
            end

            % Compute the velocity coefficients
            [dwm, ddwm] = obj.compute_dw(mOrder, mIdx, LambdaM, Vm, dLambdaM, dVm);
            obj.dw(mOrder).coeffs(:, mIdx) = dwm;
            if obj.computeSensitivity
                for dv = 1:obj.sys.ndv
                    obj.dw(mOrder).dv(dv).coeffs(:, mIdx) = ddwm(dv).coeffs;
                end
            end

            % Compute the symmetric multi-index and its position index
            mConj = [m(2); m(1)];
            mIdxConj = from_MI_to_position_idx(mConj);

            % Compute the coefficients for the symmetric multi-index
            if mIdx ~= mIdxConj
                % Reduced dynamics coefficients
                if m(1) - m(2) == 1
                    obj.R(mOrder).coeffs(2) = conj(obj.R(mOrder).coeffs(1));

                    % Loop over design variables
                    if obj.computeSensitivity
                        for dv = 1:obj.sys.ndv
                            obj.R(mOrder).dv(dv).coeffs(2) = conj(obj.R(mOrder).dv(dv).coeffs(1));
                        end
                    end
                end

                % Manifold and velocity coefficients
                obj.w(mOrder).coeffs(:, mIdxConj) = conj(obj.w(mOrder).coeffs(:, mIdx));
                obj.dw(mOrder).coeffs(:, mIdxConj) = conj(obj.dw(mOrder).coeffs(:, mIdx));

                % Store the auxiliary coefficients
                if obj.storeCoefficients
                    obj.LambdaM(mOrder).coeffs(mIdxConj) = conj(obj.LambdaM(mOrder).coeffs(mIdx));
                    obj.f(mOrder).coeffs(:, mIdxConj) = conj(obj.f(mOrder).coeffs(:, mIdx));
                    obj.V(mOrder).coeffs(:, mIdxConj) = conj(obj.V(mOrder).coeffs(:, mIdx));
                    obj.Y(mOrder).coeffs(:, mIdxConj) = conj(obj.Y(mOrder).coeffs(:, mIdx));
                    obj.C(mOrder).coeffs(:, mIdxConj) = conj(obj.C(mOrder).coeffs(:, mIdx));
                    obj.h(mOrder).coeffs(:, mIdxConj) = conj(obj.h(mOrder).coeffs(:, mIdx));

                    % Store the auxiliary sensitivity coefficients
                    if obj.computeSensitivity
                        for dv = 1:obj.sys.ndv
                            obj.LambdaM(mOrder).dv(dv).coeffs(mIdxConj) = conj(obj.LambdaM(mOrder).dv(dv).coeffs(mIdx));
                            obj.f(mOrder).dv(dv).coeffs(:, mIdxConj) = conj(obj.f(mOrder).dv(dv).coeffs(:, mIdx));
                            obj.V(mOrder).dv(dv).coeffs(:, mIdxConj) = conj(obj.V(mOrder).dv(dv).coeffs(:, mIdx));
                            obj.Y(mOrder).dv(dv).coeffs(:, mIdxConj) = conj(obj.Y(mOrder).dv(dv).coeffs(:, mIdx));
                            obj.C(mOrder).dv(dv).coeffs(:, mIdxConj) = conj(obj.C(mOrder).dv(dv).coeffs(:, mIdx));
                            obj.h(mOrder).dv(dv).coeffs(:, mIdxConj) = conj(obj.h(mOrder).dv(dv).coeffs(:, mIdx));
                        end
                    end
                end

                % Loop over design variables
                if obj.computeSensitivity
                    for dv = 1:obj.sys.ndv
                        obj.w(mOrder).dv(dv).coeffs(:, mIdxConj) = conj(obj.w(mOrder).dv(dv).coeffs(:, mIdx));
                        obj.dw(mOrder).dv(dv).coeffs(:, mIdxConj) = conj(obj.dw(mOrder).dv(dv).coeffs(:, mIdx));
                    end
                end
            end
        end
    end
end

classdef SensitivityPeakFrequency < MechSSMsensitivity
    % A class to compute the adjoint sensitivity of the peak frequency.
    % Methods:
    %   SensitivityPeakFrequency: the constructor.
    %   adjoint_sensitivity: compute the adjoint sensitivity.
    %   adjoint_equation_rho: solve the adjoint equation for the reduced
    %       amplitude.
    %   adjoint_equation_w: solve the adjoint equation for the manifold
    %       coefficients.
    %   adjoint_equation_omega_phi: solve the adjoint equations for the
    %       natural frequency and the mode shape.
    %   partial_derivative_rho: compute the partial derivative with respect
    %       to the reduced amplitude.
    %   partial_derivative_w: compute the partial derivative with respect
    %       to the manifold coefficients.
    %   partial_derivative_omega: compute the partial derivative with
    %       respect to the natural frequency.
    %   partial_derivative_phi: compute the partial derivative with respect
    %       to the mode shape.
    %   partial_derivative_mu: compute the partial derivative with respect
    %       to the design variables.

    methods
        function obj = SensitivityPeakFrequency(ssm)
            % Construct an instance of this class.
            % Inputs:
            %   ssm: the MechSSM object.
            % Outputs:
            %   obj: an instance of this class.

            % Call the constructor of the superclass
            obj = obj@MechSSMsensitivity(ssm);
        end

        function dJ = adjoint_sensitivity(obj, rho, a, epsilon)
            % Compute the adjoint sensitivity of the peak frequency.
            % Inputs:
            %   rho: the reduced amplitude.
            %   a: the coefficient a.
            %   epsilon: the force amplitude.
            % Outputs:
            %   dJ: the sensitivity vector of the peak frequency.

            % Initialize the sensitivity vector
            dJ = zeros(obj.ssm.sys.ndv, 1);

            % Adjoint equation for the reduced amplitude
            rhoAdj = obj.adjoint_equation_rho(rho, a);

            % Adjoint equation for the manifold coefficients
            wAdj = obj.adjoint_equation_w(rho, a, rhoAdj);

            % Adjoint equation for the natural frequency and mode shape
            [omegaAdj, phiAdj] = obj.adjoint_equation_omega_phi(rho, a, epsilon, rhoAdj, wAdj);

            % Partial derivative of the frequency and the cohomological equation
            [pJ_mu, pCohom_mu, pg_mu] = obj.partial_derivative_mu(rho, a, wAdj);

            % Loop over design variables
            for dv = 1:obj.ssm.sys.ndv
                dJ(dv) = pJ_mu(dv) ...
                    + phiAdj.' * (obj.ssm.sys.pK(dv).coeffs - obj.ssm.sys.omega^2 * obj.ssm.sys.pM(dv).coeffs) * obj.ssm.sys.phi ...
                    + omegaAdj * obj.ssm.sys.phi.' * obj.ssm.sys.pM(dv).coeffs * obj.ssm.sys.phi ...
                    + pCohom_mu(dv) ...
                    + rhoAdj * pg_mu(dv);
            end
        end

        function rhoAdj = adjoint_equation_rho(obj, rho, a)
            % Solve the adjoint equation for the reduced amplitude.
            % Inputs:
            %   rho: the reduced amplitude.
            %   a: the coefficient a.
            % Outputs:
            %   rhoAdj: the adjoint variable for the reduced amplitude.
            
            % Partial derivatives with respect to the reduced amplitude
            [pJ_rho, pg_rho] = obj.partial_derivative_rho(rho, a);

            % Solve adjoint equation for the reduced amplitude
            rhoAdj = -pJ_rho / pg_rho;
        end

        function wAdj = adjoint_equation_w(obj, rho, a, rhoAdj)
            % Solve the adjoint equation for the manifold coefficients.
            % Inputs:
            %   rho: the reduced amplitude.
            %   a: the coefficient a.
            %   rhoAdj: the adjoint variable for the reduced amplitude.
            % Outputs:
            %   wAdj: the adjoint variables for the manifold coefficients.
            
            % Initialize the adjoint variables
            wAdj = repmat(struct('coeffs', []) , 1, obj.ssm.maxOrder);
        
            % Loop over orders (backward)
            for mOrder = obj.ssm.maxOrder:-1:2
                % Number of multi-indices at the current order
                nMIs = obj.ssm.MIs(mOrder).n;
        
                % Initialize the adjoint variables at the current order
                wAdj(mOrder).coeffs = zeros(obj.ssm.sys.n, nMIs);
        
                % Number of unique multi-indices
                nUniqueMIs = ceil(nMIs / 2);
        
                % Loop over the multi-indices at the current order
                for mIdx = 1:nUniqueMIs
                    % Current multi-index
                    m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);
        
                    % Partial derivatives with respect to the manifold coefficients
                    [pJ_wm, pCohom_wm, pg_wm] = obj.partial_derivative_w(mOrder, mIdx, rho, a, wAdj);
        
                    % Cohomological matrix
                    LambdaM = obj.ssm.LambdaM(mOrder).coeffs(mIdx);
                    Lm = obj.ssm.sys.K + LambdaM * obj.ssm.sys.C + LambdaM^2 * obj.ssm.sys.M;
        
                    % Solve the adjoint equation for the manifold coefficients
                    rhs = -pCohom_wm.' - pJ_wm.' - rhoAdj * pg_wm.';
                    wAdj(mOrder).coeffs(:, mIdx) = lsqminnorm(Lm, rhs);
        
                    % Compute the symmetric multi-index and its position index
                    mConj = [m(2); m(1)];
                    mIdxConj = from_MI_to_position_idx(mConj);
        
                    % Compute the coefficients for the symmetric multi-index
                    if mIdx ~= mIdxConj
                        wAdj(mOrder).coeffs(:, mIdxConj) = conj(wAdj(mOrder).coeffs(:, mIdx));
                    end
                end
            end
        end

        function [omegaAdj, phiAdj] = adjoint_equation_omega_phi(obj, rho, a, epsilon, rhoAdj, wAdj)
            % Solve together the adjoint equations for the natural
            % frequency and the mode shape.
            % Inputs:
            %   rho: the reduced amplitude.
            %   a: the coefficient a.
            %   epsilon: the force amplitude.
            %   rhoAdj: the adjoint variable for the reduced amplitude.
            %   wAdj: the adjoint variables for the manifold coefficients.
            % Outputs:
            %   omegaAdj: the adjoint variable for the natural frequency.
            %   phiAdj: the adjoint variable for the mode shape.

            % Create matrix A
            A = [obj.ssm.sys.K - obj.ssm.sys.omega^2 * obj.ssm.sys.M, 2 * obj.ssm.sys.M * obj.ssm.sys.phi; ...
                2 * obj.ssm.sys.phi.' * obj.ssm.sys.M, 0];
            
            % Partial derivatives with respect to the natural frequency
            [pJ_omega, pCohom_omega, pg_omega] = obj.partial_derivative_omega(rho, a, epsilon, wAdj);

            % Partial derivatives with respect to the mode shape
            [pJ_phi, wAdj_pCohom_phi, pg_phi] = obj.partial_derivative_phi(rho, a, epsilon, wAdj);

            % Right-hand side vector
            b = [-pJ_phi.' - wAdj_pCohom_phi.' - rhoAdj * pg_phi.'; ...
                 (pJ_omega + pCohom_omega + rhoAdj * pg_omega) / obj.ssm.sys.omega];
    
            % Solve the system
            x = A \ b;
        
            % Extract the adjoint variables
            phiAdj = x(1:end-1);
            omegaAdj = x(end);
        end
        
        function [pJ_rho, pg_rho] = partial_derivative_rho(obj, rho, a)
            % Compute the partial derivative of the frequency with respect
            % to the reduced amplitude.
            % Inputs:
            %   rho: the reduced amplitude.
            %   a: the coefficient a.
            % Outputs:
            %   pJ_rho: the partial derivative of the frequency with
            %       respect to the reduced amplitude.
            %   pg_rho: the partial derivative of the force residual with
            %       respect to the reduced amplitude.

            % Initialize the partial derivative
            pJ_rho = 0;
            pa_rho = real(obj.ssm.R(1).coeffs(1));

            % Loop over orders
            for mOrder = 3:2:obj.ssm.maxOrder
                % Add contribution to the frequency
                pJ_rho = pJ_rho + imag(obj.ssm.R(mOrder).coeffs(1)) * (mOrder-1) * rho.^(mOrder-2);

                % Add contribution to the coefficient a
                pa_rho = pa_rho + real(obj.ssm.R(mOrder).coeffs(1)) * mOrder * rho.^(mOrder-1);
            end

            % Partial derivative of the force residual
            pg_rho = 2*a*pa_rho;
        end

        function [pJ_wm, pCohom_wm, pg_wm] = partial_derivative_w(obj, mOrder, mIdx, rho, a, wAdj)
            % Compute the partial derivative of the frequency with respect
            % to the manifold coefficients.
            % Inputs:
            %   mOrder: the order of the multi-index.
            %   mIdx: the index of the multi-index.
            %   rho: the reduced amplitude.
            %   a: the coefficient a.
            %   wAdj: the adjoint variable for the manifold coefficients.
            % Outputs:
            %   pJ_wm: the partial derivative of the frequency with
            %       respect to the manifold coefficients.
            %   pCohom_wm: the partial derivative of the cohomological
            %       equation with respect to the manifold coefficients.
            %   pg_wm: the partial derivative of the force residual with
            %       respect to the manifold coefficients.
        
            % Mode shape
            phi = obj.ssm.sys.phi;
            Mphi = obj.ssm.sys.M * phi;
            Cphi = obj.ssm.sys.C * phi;
            
            % Initialize partial derivatives
            pR_wm = repmat(struct('coeffs', zeros(2, obj.ssm.sys.n)), 1, obj.ssm.maxOrder);
            pJ_wm = zeros(1, obj.ssm.sys.n);
            pCohom_wm = zeros(1, obj.ssm.sys.n);
            pa_wm = zeros(1, obj.ssm.sys.n);
        
            % Loop over orders that are greater than the current order
            for qOrder = mOrder+1:obj.ssm.maxOrder
                % Number of multi-indices at the current order
                nMIs = obj.ssm.MIs(qOrder).n;
        
                % Loop over the multi-indices at the current order
                for qIdx = 1:nMIs
                    % Current multi-index
                    q = obj.ssm.MIs(qOrder).coeffs(:, qIdx);
        
                    % Adjoint variable
                    wqAdj = wAdj(qOrder).coeffs(:, qIdx);
        
                    % Eigenvalue coefficient
                    LambdaQ = obj.ssm.LambdaM(qOrder).coeffs(qIdx);
        
                    % Auxiliary variables
                    LMCphi = LambdaQ * Mphi + Cphi;
                    LMCwqAdj = (LambdaQ * obj.ssm.sys.M + obj.ssm.sys.C) * wqAdj;
                    MwqAdj = obj.ssm.sys.M * wqAdj;
        
                    % Compute the nonlinear force contribution
                    phi_pfq_wm = obj.partial_derivative_f_w(qOrder, qIdx, mOrder, mIdx, phi);
                    wqAdj_pfq_wm = obj.partial_derivative_f_w(qOrder, qIdx, mOrder, mIdx, wqAdj);
        
                    % Compute the V vector
                    LMCphi_pVq_wm = obj.partial_derivative_V_w(qOrder, qIdx, mOrder, mIdx, pR_wm, LMCphi);
                    LMCwqAdj_pVq_wm = obj.partial_derivative_V_w(qOrder, qIdx, mOrder, mIdx, pR_wm, LMCwqAdj);
        
                    % Compute the Y vector
                    Mphi_pYq_wm = obj.partial_derivative_Y_w(qOrder, qIdx, mOrder, mIdx, pR_wm, Mphi);
                    MwqAdj_pYq_wm = obj.partial_derivative_Y_w(qOrder, qIdx, mOrder, mIdx, pR_wm, MwqAdj);
                    
                    % Compute the C vector
                    phi_pCq_wm = -Mphi_pYq_wm - LMCphi_pVq_wm - phi_pfq_wm;
                    wqAdj_pCq_wm = -MwqAdj_pYq_wm - LMCwqAdj_pVq_wm - wqAdj_pfq_wm;
        
                    % Compute the reduced dynamics coefficients
                    if q(1) - q(2) == 1
                        % Compute the partial derivative of R
                        pRq_wm = phi_pCq_wm / ...
                            (LambdaQ + obj.ssm.sys.Lambda(1) + obj.ssm.sys.delta);
                        pR_wm(qOrder).coeffs(1, :) = pRq_wm;
        
                        % Add the contribution to the response frequency
                        pJ_wm = pJ_wm - 1i/2 * pRq_wm * rho^(qOrder - 1);

                        % Add the contribution to the coefficient a
                        pa_wm = pa_wm + 1/2 * pRq_wm * rho^qOrder;
                    elseif q(1) - q(2) == -1
                        % Compute the partial derivative of R
                        pRq_wm = phi_pCq_wm / ...
                            (LambdaQ + obj.ssm.sys.Lambda(2) + obj.ssm.sys.delta);
                        pR_wm(qOrder).coeffs(2, :) = pRq_wm;
        
                        % Add the contribution to the response frequency
                        pJ_wm = pJ_wm + 1i/2 * pRq_wm * rho^(qOrder - 1);

                        % Add the contribution to the coefficient a
                        pa_wm = pa_wm + 1/2 * pRq_wm * rho^qOrder;
                    end
                    
                    % Compute the vector h
                    wqAdj_phq_wm = wqAdj_pCq_wm;
                    if q(1) - q(2) == 1
                        % Compute the D vector
                        Dq = -((LambdaQ + obj.ssm.sys.Lambda(1)) * obj.ssm.sys.M + obj.ssm.sys.C) * obj.ssm.sys.phi;
        
                        % Add the contribution to the cohomological equation
                        wqAdj_phq_wm = wqAdj_phq_wm + (wqAdj.' * Dq) * pR_wm(qOrder).coeffs(1, :);
                    elseif q(1) - q(2) == -1
                        % Compute the D vector
                        Dq = -((LambdaQ + obj.ssm.sys.Lambda(2)) * obj.ssm.sys.M + obj.ssm.sys.C) * obj.ssm.sys.phi;
        
                        % Add the contribution to the cohomological equation
                        wqAdj_phq_wm = wqAdj_phq_wm + (wqAdj.' * Dq) * pR_wm(qOrder).coeffs(2, :);
                    end
        
                    % Compute the cohomological equation
                    pCohom_wm = pCohom_wm - wqAdj_phq_wm;
                end
            end

            % Partial derivative of the force residual
            pg_wm = 2*a*pa_wm;
        end

        function [pJ_omega, pCohom_omega, pg_omega] = partial_derivative_omega(obj, rho, a, epsilon, wAdj)
            % Compute the partial derivative of the frequency with respect
            % to the natural frequency.
            % Inputs:
            %   rho: the reduced amplitude.
            %   a: the coefficient a.
            %   epsilon: the force amplitude.
            %   wAdj: the adjoint variables for the manifold coefficients.
            % Outputs:
            %   pJ_omega: the partial derivative of the frequency with
            %       respect to the natural frequency.
            %   pCohom_omega: the partial derivative of the cohomological
            %       equation with respect to the natural frequency.
            %   pg_omega: the partial derivative of the force residual with
            %       respect to the natural frequency.
        
            % Initialize partial derivatives
            pR_omega = repmat(struct('coeffs', []), 1, obj.ssm.maxOrder);
            pdw_omega = repmat(struct('coeffs', []), 1, obj.ssm.maxOrder);
        
            % Initialize partial derivatives
            pJ_omega = imag(obj.ssm.sys.pLambda_omega(1));
            pCohom_omega = 0;
            pa_omega = real(obj.ssm.sys.pLambda_omega(1)) * rho;
        
            % Loop over orders
            for mOrder = 2:obj.ssm.maxOrder
                % Number of multi-indices at this order
                nMIs = obj.ssm.MIs(mOrder).n;
        
                % Initialize coefficients
                pR_omega(mOrder).coeffs = zeros(1, 2);
                pdw_omega(mOrder).coeffs = zeros(obj.ssm.sys.n, nMIs);
        
                % Number of unique multi-indices
                nUniqueMIs = ceil(nMIs / 2);
        
                % Loop over unique multi-indices
                for mIdx = 1:nUniqueMIs
                    % Current multi-index
                    m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);
        
                    % Eigenvalue coefficient
                    LambdaM = obj.ssm.LambdaM(mOrder).coeffs(mIdx);
                    pLambdaM_omega = dot(m, obj.ssm.sys.pLambda_omega);
        
                    % Compute the V vector
                    Vm = obj.ssm.V(mOrder).coeffs(:, mIdx);
                    pVm_omega = obj.partial_derivative_V_omega(mOrder, mIdx, pR_omega);
        
                    % Compute the Y vector
                    pYm_omega = obj.partial_derivative_Y_omega(mOrder, mIdx, pR_omega, pdw_omega);
        
                    % Compute the C vector
                    pCm_omega = -(pLambdaM_omega * obj.ssm.sys.M + obj.ssm.sys.pC_omega) * Vm - ...
                        (LambdaM * obj.ssm.sys.M + obj.ssm.sys.C) * pVm_omega - ...
                        obj.ssm.sys.M * pYm_omega;
        
                    % Compute the reduced dynamics coefficients
                    if m(1) - m(2) == 1
                        % Extract the R coefficient
                        Rm = obj.ssm.R(mOrder).coeffs(1);
        
                        % Compute the partial derivative of R
                        pNum_omega = obj.ssm.sys.phi.' * pCm_omega;
                        den = LambdaM + obj.ssm.sys.Lambda(1) + obj.ssm.sys.delta;
                        pDen_omega = pLambdaM_omega + obj.ssm.sys.pLambda_omega(1) + obj.ssm.sys.pdelta_omega;
                        pRm_omega = (pNum_omega - Rm * pDen_omega) / den;
        
                        % Store the partial derivative of R
                        pR_omega(mOrder).coeffs(1) = pRm_omega;
        
                        % Add the contribution to the response frequency
                        pJ_omega = pJ_omega - 1i/2 * pRm_omega * rho^(mOrder - 1);

                        % Add the contribution to the coefficient a
                        pa_omega = pa_omega + 1/2 * pRm_omega * rho^mOrder;
                    end
                    
                    % Compute the vector h
                    phm_omega = pCm_omega;
                    if m(1) - m(2) == 1
                        % Compute vector D
                        Dm = -((LambdaM + obj.ssm.sys.Lambda(1)) * obj.ssm.sys.M + obj.ssm.sys.C) * obj.ssm.sys.phi;
        
                        % Partial derivative of D
                        pDm_omega = -((pLambdaM_omega + obj.ssm.sys.pLambda_omega(1)) * obj.ssm.sys.M + obj.ssm.sys.pC_omega) * obj.ssm.sys.phi;
        
                        % Compute the partial derivative
                        phm_omega = phm_omega + pRm_omega * Dm + Rm * pDm_omega;
                    end
        
                    % Partial derivative of L
                    pLm_omega = pLambdaM_omega * obj.ssm.sys.C + LambdaM * obj.ssm.sys.pC_omega ...
                        + 2 * LambdaM * pLambdaM_omega * obj.ssm.sys.M;
        
                    % Extract the coefficients
                    wm = obj.ssm.w(mOrder).coeffs(:, mIdx);
                    wmAdj = wAdj(mOrder).coeffs(:, mIdx);
        
                    % Partial derivative of the cohomological equation
                    cohomAdj = wmAdj.' * (pLm_omega * wm - phm_omega);
                    pCohom_omega = pCohom_omega + cohomAdj;
        
                    % Compute the velocity coefficients
                    pdwm_omega = pLambdaM_omega * wm + pVm_omega;
                    if m(1) - m(2) == 1
                        pdwm_omega = pdwm_omega + obj.ssm.sys.phi * pR_omega(mOrder).coeffs(1);
                    end
                    pdw_omega(mOrder).coeffs(:, mIdx) = pdwm_omega;
        
                    % Compute the symmetric multi-index and its position index
                    mConj = [m(2); m(1)];
                    mIdxConj = from_MI_to_position_idx(mConj);
        
                    % Compute the coefficients for the symmetric multi-index
                    if mIdx ~= mIdxConj
                        if m(1) - m(2) == 1
                            % Compute the partial derivative of R
                            pR_omega(mOrder).coeffs(2) = conj(pRm_omega);
            
                            % Add the contribution to the response frequency
                            pJ_omega = pJ_omega + 1i/2 * conj(pRm_omega) * rho^(mOrder - 1);

                            % Add the contribution to the coefficient a
                            pa_omega = pa_omega + 1/2 * conj(pRm_omega) * rho^mOrder;
                        end
        
                        % Compute the partial derivative of the cohomological equation
                        pCohom_omega = pCohom_omega + conj(cohomAdj);
        
                        % Store the partial derivative of the velocity coefficients
                        pdw_omega(mOrder).pdv(:, mIdxConj) = conj(pdwm_omega);
                    end
                end
            end
        
            % Ensure real values
            pJ_omega = real(pJ_omega);
            pCohom_omega = real(pCohom_omega);
            pa_omega = real(pa_omega);

            % Partial derivative of the non-autonomous coefficient
            den = obj.ssm.sys.delta + obj.ssm.sys.lambda + obj.ssm.Lambda0;
            pden_omega = obj.ssm.sys.pdelta_omega + obj.ssm.sys.plambda_omega + 1i * imag(obj.ssm.sys.plambda_omega);
            pS0_omega = -obj.ssm.S0 * pden_omega / den;

            % Partial derivative of the force residual
            pg_omega = 2*a*pa_omega - epsilon^2 * 2 * real(conj(obj.ssm.S0) * pS0_omega);
        end

        function [pJ_phi, wAdj_pCohom_phi, pg_phi] = partial_derivative_phi(obj, rho, a, epsilon, wAdj)
            % Comptue the partial derivative of the frequency with respect
            % to the mode shape.
            % Inputs:
            %   rho: the reduced amplitude.
            %   a: the coefficient a.
            %   epsilon: the force amplitude.
            %   wAdj: the adjoint variables for the manifold coefficients.
            % Outputs:
            %   pJ_phi: the partial derivative of the frequency with
            %       respect to the mode shape.
            %   wAdj_pCohom_phi: the partial derivative of the
            %       cohomological equation with respect to the mode shape
            %       premultiplied by the adjoint variable.
            %   pg_phi: the partial derivative of the force residual with
            %       respect to the mode shape.
        
            % Mode shape
            phi = obj.ssm.sys.phi;
            Mphi = obj.ssm.sys.M * phi;
            Cphi = obj.ssm.sys.C * phi;
        
            % Initialize
            pR_phi = repmat(struct('coeffs', zeros(2, obj.ssm.sys.n)), 1, obj.ssm.maxOrder);
            pJ_phi = zeros(1, obj.ssm.sys.n);
            wAdj_pCohom_phi = zeros(1, obj.ssm.sys.n);
            pa_phi = zeros(1, obj.ssm.sys.n);
        
            % Loop over orders
            for mOrder = 2:obj.ssm.maxOrder
                % Number of multi-indices at this order
                nMIs = obj.ssm.MIs(mOrder).n;
        
                % Number of unique multi-indices
                nUniqueMIs = ceil(nMIs / 2);
        
                % Loop over multi-indices
                for mIdx = 1:nUniqueMIs
                    % Current multi-index
                    m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);
        
                    % Adjoint variable
                    wmAdj = wAdj(mOrder).coeffs(:, mIdx);
        
                    % Eigenvalue coefficient
                    LambdaM = obj.ssm.LambdaM(mOrder).coeffs(mIdx);
        
                    % Auxiliary variables
                    LMCphi = LambdaM * Mphi + Cphi;
                    LMCwmAdj = (LambdaM * obj.ssm.sys.M + obj.ssm.sys.C) * wmAdj;
                    MwmAdj = obj.ssm.sys.M * wmAdj;
        
                    % Compute the nonlinear force contribution
                    phi_pfm_phi = obj.partial_derivative_f_phi(mOrder, mIdx, phi);
                    wmAdj_pfm_phi = obj.partial_derivative_f_phi(mOrder, mIdx, wmAdj);
        
                    % Compute the V vector
                    LMCphi_pVm_phi = obj.partial_derivative_V_phi(mOrder, mIdx, pR_phi, LMCphi);
                    LMCwmAdj_pVm_phi = obj.partial_derivative_V_phi(mOrder, mIdx, pR_phi, LMCwmAdj);
        
                    % Compute the Y vector
                    Mphi_pYm_phi = obj.partial_derivative_Y_phi(mOrder, mIdx, pR_phi, Mphi);
                    MwmAdj_pYm_phi = obj.partial_derivative_Y_phi(mOrder, mIdx, pR_phi, MwmAdj);
        
                    % Compute the C vector
                    Cm = obj.ssm.C(mOrder).coeffs(:, mIdx);
                    phi_pCm_phi = -Mphi_pYm_phi - LMCphi_pVm_phi - phi_pfm_phi;
                    wmAdj_pCm_phi = -MwmAdj_pYm_phi - LMCwmAdj_pVm_phi - wmAdj_pfm_phi;
        
                    % Compute the reduced dynamics coefficients
                    if m(1) - m(2) == 1
                        pRm_phi = (Cm.' + phi_pCm_phi) / ...
                            (LambdaM + obj.ssm.sys.Lambda(1) + obj.ssm.sys.delta);
                        pR_phi(mOrder).coeffs(1, :) = pRm_phi;
        
                        % Add the contribution to the response frequency
                        pJ_phi = pJ_phi - 1i/2 * pRm_phi * rho^(mOrder - 1);

                        % Add the contribution to the coefficient a
                        pa_phi = pa_phi + 1/2 * pRm_phi * rho^mOrder;
                    end
                    
                    % Compute the vector h
                    wmAdj_phm_phi = wmAdj_pCm_phi;
                    if m(1) - m(2) == 1
                        % Partial derivative of D
                        Dm = -((LambdaM + obj.ssm.sys.Lambda(1)) * obj.ssm.sys.M + obj.ssm.sys.C) * obj.ssm.sys.phi;
                        wmAdj_pDm_phi = -wmAdj.' * ((LambdaM + obj.ssm.sys.Lambda(1)) * obj.ssm.sys.M + obj.ssm.sys.C);
        
                        % Compute the partial derivative
                        wmAdj_phm_phi = wmAdj_phm_phi + (wmAdj.' * Dm) * pR_phi(mOrder).coeffs(1, :) + ...
                            obj.ssm.R(mOrder).coeffs(1) * wmAdj_pDm_phi;
                    end
        
                    % Compute the cohomological equation
                    wAdj_pCohom_phi = wAdj_pCohom_phi - wmAdj_phm_phi;
        
                    % Compute the symmetric multi-index and its position index
                    mConj = [m(2); m(1)];
                    mIdxConj = from_MI_to_position_idx(mConj);
        
                    % Compute the coefficients for the symmetric multi-index
                    if mIdx ~= mIdxConj
                        if m(1) - m(2) == 1
                            % Compute the partial derivative of R
                            pR_phi(mOrder).coeffs(2, :) = conj(pRm_phi);
            
                            % Add the contribution to the response frequency
                            pJ_phi = pJ_phi + 1i/2 * conj(pRm_phi) * rho^(mOrder - 1);

                            % Add the contribution to the coefficient a
                            pa_phi = pa_phi + 1/2 * conj(pRm_phi) * rho^mOrder;
                        end
        
                        % Compute the partial derivative of the cohomological equation
                        wAdj_pCohom_phi = wAdj_pCohom_phi - conj(wmAdj_phm_phi);
                    end
                end
            end
        
            % Ensure real values
            pJ_phi = real(pJ_phi);
            wAdj_pCohom_phi = real(wAdj_pCohom_phi);
            pa_phi = real(pa_phi);

            % Partial derivative of the non-autonomous coefficient
            num = -1/2;
            den = obj.ssm.sys.delta + obj.ssm.sys.lambda + obj.ssm.Lambda0;
            pS0_phi = zeros(1, obj.ssm.sys.n);
            pS0_phi(obj.ssm.inDof) = num / den;

            % Partial derivative of the force residual
            pg_phi = 2*a*pa_phi - epsilon^2 * 2 * real(conj(obj.ssm.S0) * pS0_phi);
        end

        function [pJ_mu, pCohom_mu, pg_mu] = partial_derivative_mu(obj, rho, a, wAdj)
            % Compute the partial derivative of the response frequency with
            % respect to the design variables.
            % Inputs:
            %   rho: the reduced amplitude.
            %   a: the coefficient a.
            %   wAdj: the adjoint variables for the manifold coefficients.
            % Outputs:
            %   pJ_mu: the partial derivative of the frequency with
            %       respect to the design variables.
            %   pCohom_mu: the partial derivative of the cohomological
            %       equation with respect to the design variables.
            %   pg_mu: the partial derivative of the coefficient a with
            %       respect to the design variables.
        
            % Initialize partial derivatives
            pR_mu = repmat(struct('pdv', []), 1, obj.ssm.maxOrder);
            pdw_mu = repmat(struct('pdv', []), 1, obj.ssm.maxOrder);
            pCohom_mu = zeros(obj.ssm.sys.ndv, 1);
            pJ_mu = zeros(obj.ssm.sys.ndv, 1);
            pa_mu = zeros(obj.ssm.sys.ndv, 1);
        
            % Loop over orders
            for mOrder = 2:obj.ssm.maxOrder
                % Number of multi-indices at this order
                nMIs = obj.ssm.MIs(mOrder).n;
        
                % Initialize coefficients
                pR_mu(mOrder).pdv = repmat(struct('coeffs', zeros(1, 2)), 1, obj.ssm.sys.ndv);
                pdw_mu(mOrder).pdv = repmat(struct('coeffs', zeros(obj.ssm.sys.n, nMIs)), 1, obj.ssm.sys.ndv);
        
                % Number of unique multi-indices
                nUniqueMIs = ceil(nMIs / 2);
        
                % Loop over multi-indices
                for mIdx = 1:nUniqueMIs
                    % Current multi-index
                    m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);
        
                    % Eigenvalue coefficient
                    LambdaM = obj.ssm.LambdaM(mOrder).coeffs(mIdx);
        
                    % Extract coefficients
                    wm = obj.ssm.w(mOrder).coeffs(:, mIdx);
                    wmAdj = wAdj(mOrder).coeffs(:, mIdx);
        
                    % Compute the nonlinear force contribution
                    pfm_mu = obj.partial_derivative_f_mu(mOrder, mIdx);
        
                    % Compute the V vector
                    Vm = obj.ssm.V(mOrder).coeffs(:, mIdx);
                    pVm_mu = obj.partial_derivative_V_mu(mOrder, mIdx, pR_mu);
        
                    % Compute the Y vector
                    Ym = obj.ssm.Y(mOrder).coeffs(:, mIdx);
                    pYm_mu = obj.partial_derivative_Y_mu(mOrder, mIdx, pR_mu, pdw_mu);
        
                    % Loop over design variables
                    pCm_mu = repmat(struct('coeffs', []), 1, obj.ssm.sys.ndv);
                    for dv = 1:obj.ssm.sys.ndv
                        % Compute the partial derivative of C
                        pCm_mu(dv).coeffs = -(LambdaM * obj.ssm.sys.pM(dv).coeffs + obj.ssm.sys.pC(dv).coeffs) * Vm - ...
                            (LambdaM * obj.ssm.sys.M + obj.ssm.sys.C) * pVm_mu(dv).coeffs - ...
                            obj.ssm.sys.pM(dv).coeffs * Ym - obj.ssm.sys.M * pYm_mu(dv).coeffs - ...
                            pfm_mu(dv).coeffs;
                    end
        
                    % Loop over design variables
                    for dv = 1:obj.ssm.sys.ndv
                        if m(1) - m(2) == 1
                            % Compute the partial derivative of R
                            pR_mu(mOrder).pdv(dv).coeffs(1) = obj.ssm.sys.phi.' * pCm_mu(dv).coeffs / ...
                                (LambdaM + obj.ssm.sys.Lambda(1) + obj.ssm.sys.delta);
                            
                            % Add the contribution to the response frequency
                            pJ_mu(dv) = pJ_mu(dv) - 1i/2 * pR_mu(mOrder).pdv(dv).coeffs(1) * rho^(mOrder - 1);

                            % Add the contribution to the coefficient a
                            pa_mu(dv) = pa_mu(dv) + 1/2 * pR_mu(mOrder).pdv(dv).coeffs(1) * rho^mOrder;
                        end
                    end
                    
                    % Compute the vector h
                    phm_mu = pCm_mu;
                    if m(1) - m(2) == 1
                        % Compute vector D
                        Dm = -((LambdaM + obj.ssm.sys.Lambda(1)) * obj.ssm.sys.M + obj.ssm.sys.C) * obj.ssm.sys.phi;
        
                        % Loop over design variables
                        for dv = 1:obj.ssm.sys.ndv
                            % Partial derivative of D
                            pDm_mu = -((LambdaM + obj.ssm.sys.Lambda(1)) * obj.ssm.sys.pM(dv).coeffs + ...
                                obj.ssm.sys.pC(dv).coeffs) * obj.ssm.sys.phi;
                            
                            % Compute the partial derivative
                            phm_mu(dv).coeffs = phm_mu(dv).coeffs + ...
                                pR_mu(mOrder).pdv(dv).coeffs(1) * Dm + obj.ssm.R(mOrder).coeffs(1) * pDm_mu;
                        end
                    end
        
                    % Loop over design variables
                    cohomAdj = repmat(struct('coeffs', []), 1, obj.ssm.sys.ndv);
                    for dv = 1:obj.ssm.sys.ndv
                        % Compute the partial derivative of L
                        pLm_mu = obj.ssm.sys.pK(dv).coeffs + LambdaM * obj.ssm.sys.pC(dv).coeffs + ...
                            LambdaM^2 * obj.ssm.sys.pM(dv).coeffs;
        
                        % Compute the partial derivative of the cohomological equation
                        cohomAdj(dv).coeffs = wmAdj.' * (pLm_mu * wm - phm_mu(dv).coeffs);
                        pCohom_mu(dv) = pCohom_mu(dv) + cohomAdj(dv).coeffs;
                    end
        
                    % Loop over design variables
                    for dv = 1:obj.ssm.sys.ndv
                        % Initialize with the partial derivative of the V vector
                        pdwm_mu = pVm_mu(dv).coeffs;
        
                        % Add the contribution from the reduced dynamics coefficients
                        if m(1) - m(2) == 1
                            pdwm_mu = pdwm_mu + ...
                                obj.ssm.sys.phi * pR_mu(mOrder).pdv(dv).coeffs(1);
                        end
        
                        % Store the partial derivative of the velocity coefficients
                        pdw_mu(mOrder).pdv(dv).coeffs(:, mIdx) = pdwm_mu;
                    end
        
                    % Compute the symmetric multi-index and its position index
                    mConj = [m(2); m(1)];
                    mIdxConj = from_MI_to_position_idx(mConj);
        
                    % Compute the coefficients for the symmetric multi-index
                    if mIdx ~= mIdxConj
                        % Loop over design variables
                        for dv = 1:obj.ssm.sys.ndv
                            if m(1) - m(2) == 1
                                % Compute the partial derivative of R
                                pR_mu(mOrder).pdv(dv).coeffs(2) = ...
                                    conj(pR_mu(mOrder).pdv(dv).coeffs(1));
        
                                % Add the contribution to the response frequency
                                pJ_mu(dv) = pJ_mu(dv) ...
                                    + 1i/2 * pR_mu(mOrder).pdv(dv).coeffs(2) * rho^(mOrder - 1);
                            
                                % Add the contribution to the coefficient a
                                pa_mu(dv) = pa_mu(dv) + 1/2 * pR_mu(mOrder).pdv(dv).coeffs(2) * rho^mOrder;
                            end
        
                            % Compute the partial derivative of the cohomological equation
                            pCohom_mu(dv) = pCohom_mu(dv) + conj(cohomAdj(dv).coeffs);
        
                            % Store the partial derivative of the velocity coefficients
                            pdw_mu(mOrder).pdv(dv).coeffs(:, mIdxConj) = ...
                                conj(pdw_mu(mOrder).pdv(dv).coeffs(:, mIdx));
                        end
                    end
                end
            end
        
            % Ensure real values
            pJ_mu = real(pJ_mu);
            pCohom_mu = real(pCohom_mu);
            pa_mu = real(pa_mu);

            % Partial derivative of the force residual
            pg_mu = 2*a*pa_mu;
        end
    end
end

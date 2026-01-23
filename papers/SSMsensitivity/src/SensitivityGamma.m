classdef SensitivityGamma < MechSSMsensitivity
    % A subclass of MechSSMsensitivity that is specifically designed to
    % handle the sensitivity analysis of the reduced backbone coefficients
    % gamma at a given order.
    % Properties:
    %   order: the order of the SSM coefficient gamma.
    % Methods:
    %   SensitivityGamma: the constructor.
    %   adjoint_sensitivity: compute the adjoint sensitivity of the reduced
    %       backbone coefficients gamma.
    %   adjoint_sensitivity_topopt: compute the adjoint sensitivity of the
    %       reduced backbone coefficients gamma for topology optimization.
    %   adjoint_equation_w: solve the adjoint equation of the manifold
    %       coefficients w.
    %   adjoint_equation_omega_phi: solve the adjoint equation of the
    %       natural frequency omega and the mode shape phi.
    %   partial_derivative_w: compute the partial derivative of the reduced
    %       backbone coefficients gamma with respect to the manifold
    %       coefficients w.
    %   partial_derivative_omega: compute the partial derivative of the
    %       reduced backbone coefficients gamma with respect to the natural
    %       frequency omega.
    %   partial_derivative_phi: compute the partial derivative of the
    %       reduced backbone coefficients gamma with respect to the mode
    %       shape phi.
    %   partial_derivative_mu: compute the partial derivative of the
    %       reduced backbone coefficients gamma with respect to the design
    %       variables mu.
    %   partial_derivative_mu_topopt: compute the partial derivative of the
    %       reduced backbone coefficients gamma with respect to the design
    %       variables mu for topology optimization.

    properties
        order % The order of the SSM coefficient gamma
    end

    methods
        function obj = SensitivityGamma(ssm, varargin)
            % Construct an instance of this class.
            % Inputs:
            %   ssm: the MechSSM object.
            %   order: the order of the SSM coefficient gamma.
            % Outputs:
            %   obj: an instance of this class.

            % Call the constructor of the superclass
            obj = obj@MechSSMsensitivity(ssm);

            % Parse the input arguments
            p = inputParser;
            addOptional(p, 'order', 0);
            parse(p, varargin{:});
            obj.order = p.Results.order;

            % Assign the default order if not specified
            if obj.order == 0
                obj.order = obj.ssm.maxOrder;
            end

            % Make sure the order is not greater than the SSM order
            if obj.order > obj.ssm.maxOrder
                error('The order of the SSM coefficient gamma cannot be greater than the order of the SSM.');
            end

            % Make sure the order is odd
            if mod(obj.order, 2) == 0
                error('The order of the SSM coefficient gamma must be odd.');
            end
        end

        function dJ = adjoint_sensitivity(obj)
            % Compute the adjoint sensitivity of the reduced backbone
            % coefficients gamma.
            % Outputs:
            %   dJ: the adjoint sensitivity of gamma.

            % Initialize the adjoint sensitivity of gamma
            dJ = zeros(obj.ssm.sys.ndv, 1);

            % Adjoint equation for the manifold coefficients w
            wAdj = obj.adjoint_equation_w();

            % Adjoint equation for the natural frequency omega and the mode
            % shape phi
            [omegaAdj, phiAdj] = obj.adjoint_equation_omega_phi(wAdj);

            % Partial derivative with respect to the design variables mu 
            [pJ_mu, pCohom_mu] = obj.partial_derivative_mu(wAdj);

            % Loop over the design variables
            for dv = 1:obj.ssm.sys.ndv
                % Add the contribution of gamma
                dJ(dv) = pJ_mu(dv) ...
                    + phiAdj.' * (obj.ssm.sys.pK(dv).coeffs - obj.ssm.sys.omega^2 * obj.ssm.sys.pM(dv).coeffs) * obj.ssm.sys.phi ...
                    + omegaAdj * obj.ssm.sys.phi.' * obj.ssm.sys.pM(dv).coeffs * obj.ssm.sys.phi ...
                    + pCohom_mu(dv);
            end
        end

        function dJ = adjoint_sensitivity_topopt(obj)
            % Compute the adjoint sensitivity of the reduced backbone
            % coefficients gamma for topology optimization.
            % Outputs:
            %   dJ: the adjoint sensitivity of gamma.

            % Extract FE assembly
            yafecAssembly = obj.ssm.sys.yafecAssembly;
            nDOFs = yafecAssembly.Mesh.nDOFs;

            % Initialize the adjoint sensitivity of gamma
            dJ = zeros(obj.ssm.sys.ndv, 1);

            % Adjoint equation for the manifold coefficients w
            wAdj = obj.adjoint_equation_w();

            % Adjoint equation for the natural frequency omega and the mode
            % shape phi
            [omegaAdj, phiAdj] = obj.adjoint_equation_omega_phi(wAdj);

            % Convert vectors into their unconstrained form
            obj.phiUnc = yafecAssembly.unconstrain_vector(obj.ssm.sys.phi);
            obj.phiAdjUnc = yafecAssembly.unconstrain_vector(phiAdj);
            
            obj.wUnc = repmat(struct('coeffs', []), 1, obj.order);
            obj.wAdjUnc = repmat(struct('coeffs', []), 1, obj.order);
            obj.wUnc(1).coeffs = [obj.phiUnc, obj.phiUnc];
            for mOrder = 2:(obj.order - 1)
                % Number of multi-indices at this order
                nMIs = obj.ssm.MIs(mOrder).n;
        
                % Initialize coefficients
                obj.wUnc(mOrder).coeffs = zeros(nDOFs, nMIs);
                obj.wAdjUnc(mOrder).coeffs = zeros(nDOFs, nMIs);
        
                % Loop over multi-indices
                for mIdx = 1:nMIs
                    obj.wUnc(mOrder).coeffs(:, mIdx) = yafecAssembly.unconstrain_vector(obj.ssm.w(mOrder).coeffs(:, mIdx));
                    obj.wAdjUnc(mOrder).coeffs(:, mIdx) = yafecAssembly.unconstrain_vector(wAdj(mOrder).coeffs(:, mIdx));
                end
            end

            obj.VUnc = repmat(struct('coeffs', []), 1, obj.order);
            obj.YUnc = repmat(struct('coeffs', []), 1, obj.order);
            for mOrder = 2:obj.order
                % Number of multi-indices at this order
                nMIs = obj.ssm.MIs(mOrder).n;
        
                % Initialize coefficients
                obj.VUnc(mOrder).coeffs = zeros(nDOFs, nMIs);
                obj.YUnc(mOrder).coeffs = zeros(nDOFs, nMIs);
        
                % Loop over multi-indices
                for mIdx = 1:nMIs
                    obj.VUnc(mOrder).coeffs(:, mIdx) = yafecAssembly.unconstrain_vector(obj.ssm.V(mOrder).coeffs(:, mIdx));
                    obj.YUnc(mOrder).coeffs(:, mIdx) = yafecAssembly.unconstrain_vector(obj.ssm.Y(mOrder).coeffs(:, mIdx));
                end
            end

            % Loop over the design variables (elements)
            for dv = 1:obj.ssm.sys.ndv
                % Extract the dofs of the current design variable
                thisElement = yafecAssembly.Mesh.Elements(dv).Object;
                iDOFs = thisElement.iDOFs;

                % Extract element quantities
                phiE = obj.phiUnc(iDOFs);
                phiAdjE = obj.phiAdjUnc(iDOFs);

                % Partial derivative with respect to the design variables mu 
                [pJ_mu, pCohom_mu] = obj.partial_derivative_mu_topopt(wAdj, iDOFs);

                % Add the contribution of gamma
                dJ(dv) = pJ_mu ...
                    + phiAdjE.' * (obj.ssm.sys.Ke - obj.ssm.sys.omega^2 * obj.ssm.sys.Me) * phiE ...
                    + omegaAdj * phiE.' * obj.ssm.sys.Me * phiE ...
                    + pCohom_mu;
            end
        end

        function wAdj = adjoint_equation_w(obj)
            % Solve the adjoint equation of the manifold coefficients w.
            % Outputs:
            %   wAdj: the adjoint variables for the manifold coefficients.

            % Initialize adjoint variables
            wAdj = repmat(struct('coeffs', []) , 1, obj.order);

            % Loop over orders (backward)
            for mOrder = obj.order-1:-1:2
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

                    % Partial derivatives of gamma and the cohomological equations
                    [pJ_wm, pCohom_wm] = obj.partial_derivative_w(mOrder, mIdx, wAdj);

                    % Cohomological matrix
                    LambdaM = obj.ssm.LambdaM(mOrder).coeffs(mIdx);
                    Lm = obj.ssm.sys.K + LambdaM * obj.ssm.sys.C + LambdaM^2 * obj.ssm.sys.M;

                    % Solve the adjoint equation for the manifold coefficients w
                    rhs = -pCohom_wm.' - pJ_wm.';
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

        function [omegaAdj, phiAdj] = adjoint_equation_omega_phi(obj, wAdj)
            % Solve together the adjoint equations for the natural
            % frequency omega and the mode shape phi.
            % Inputs:
            %   wAdj: the adjoint variables for the manifold coefficients.
            % Outputs:
            %   omegaAdj: the adjoint variable for omega.
            %   phiAdj: the adjoint variable for phi.

            % Create matrix A
            A = [obj.ssm.sys.K - obj.ssm.sys.omega^2 * obj.ssm.sys.M, 2 * obj.ssm.sys.M * obj.ssm.sys.phi; ...
                 2 * obj.ssm.sys.phi.' * obj.ssm.sys.M, 0];
            
            % Partial derivatives with respect to the mode shape phi
            [pJ_phi, pCohom_phi] = obj.partial_derivative_phi(wAdj);

            % Partial derivatives of J and the cohomological equations
            [pJ_omega, pCohom_omega] = obj.partial_derivative_omega(wAdj);

            % Create vector b
            b = [-pJ_phi.' - pCohom_phi.'; ...
                 (pJ_omega + pCohom_omega) / obj.ssm.sys.omega];
            
            % Solve the system
            x = A \ b;

            % Extract the adjoint variables
            phiAdj = x(1:end-1);
            omegaAdj = x(end);
        end

        function [pJ_wm, pCohom_wm] = partial_derivative_w(obj, mOrder, mIdx, wAdj)
            % Compute the partial derivative of the reduced backbone
            % coefficients gamma with respect to the manifold coefficient
            % wm.
            % Inputs:
            %   mOrder: the order of the multi-index.
            %   mIdx: the index of the multi-index in the multi-index
            %       array.
            %   wAdj: the adjoint variables for w.
            % Outputs:
            %   pJ_wm: the partial derivative of gamma with respect to wm.
            %   pCohom_wm: the partial derivative of the cohomological
            %       equation with respect to wm.

            % Mode shape
            phi = obj.ssm.sys.phi;
            Mphi = obj.ssm.sys.M * phi;
            Cphi = obj.ssm.sys.C * phi;

            % Initialize partial derivatives
            % pJ_wm = zeros(1, obj.ssm.sys.n);
            pCohom_wm = zeros(1, obj.ssm.sys.n);

            % Current multi-index
            pR_wm = repmat(struct('coeffs', zeros(2, obj.ssm.sys.n)), 1, obj.order);

            % Loop over orders that are greater than the current order
            % Stop at one order before the order of gamma
            for qOrder = mOrder+1:(obj.order-1)
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
                    elseif q(1) - q(2) == -1
                        % Compute the partial derivative of R
                        pRq_wm = phi_pCq_wm / ...
                            (LambdaQ + obj.ssm.sys.Lambda(2) + obj.ssm.sys.delta);
                        pR_wm(qOrder).coeffs(2, :) = pRq_wm;
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

            % Compute the last order
            qOrder = obj.order;
            nMIs = obj.ssm.MIs(qOrder).n;
            for qIdx = 1:nMIs
                % Current multi-index
                q = obj.ssm.MIs(qOrder).coeffs(:, qIdx);
        
                % Eigenvalue coefficient
                LambdaQ = obj.ssm.LambdaM(qOrder).coeffs(qIdx);

                % Auxiliary variables
                LMCphi = LambdaQ * Mphi + Cphi;

                % Compute the nonlinear force contribution
                phi_pfq_wm = obj.partial_derivative_f_w(qOrder, qIdx, mOrder, mIdx, phi);
    
                % Compute the V vector
                LMCphi_pVq_wm = obj.partial_derivative_V_w(qOrder, qIdx, mOrder, mIdx, pR_wm, LMCphi);
    
                % Compute the Y vector
                Mphi_pYq_wm = obj.partial_derivative_Y_w(qOrder, qIdx, mOrder, mIdx, pR_wm, Mphi);
                    
                % Compute the C vector
                phi_pCq_wm = -Mphi_pYq_wm - LMCphi_pVq_wm - phi_pfq_wm;

                % Compute the reduced dynamics coefficients
                if q(1) - q(2) == 1
                    % Compute the partial derivative of R
                    pRq_wm = phi_pCq_wm / ...
                        (LambdaQ + obj.ssm.sys.Lambda(1) + obj.ssm.sys.delta);
                    pR_wm(qOrder).coeffs(1, :) = pRq_wm;
                elseif q(1) - q(2) == -1
                    % Compute the partial derivative of R
                    pRq_wm = phi_pCq_wm / ...
                        (LambdaQ + obj.ssm.sys.Lambda(2) + obj.ssm.sys.delta);
                    pR_wm(qOrder).coeffs(2, :) = pRq_wm;
                end
            end

            % Compute the partial derivative of gamma
            pJ_wm = 1i/2 * (pR_wm(obj.order).coeffs(2, :) - pR_wm(obj.order).coeffs(1, :));
        end

        function [pJ_omega, pCohom_omega] = partial_derivative_omega(obj, wAdj)
            % Compute the partial derivative of the reduced backbone
            % coefficients gamma with respect to the natural frequency
            % omega.
            % Inputs:
            %   wAdj: the adjoint variables for the manifold coefficients.
            % Outputs:
            %   pJ_omega: the partial derivative of gamma with respect to
            %       omega.
            %   pCohom_omega: the partial derivative of the cohomological
            %       equation with respect to omega.

            % Initialize partial derivatives
            pR_omega = repmat(struct('coeffs', []), 1, obj.order);
            pdw_omega = repmat(struct('coeffs', []), 1, obj.order);
            pCohom_omega = 0;

            % Loop over orders
            for mOrder = 2:(obj.order - 1)
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
                        end
        
                        % Compute the partial derivative of the cohomological equation
                        pCohom_omega = pCohom_omega + conj(cohomAdj);
        
                        % Store the partial derivative of the velocity coefficients
                        pdw_omega(mOrder).pdv(:, mIdxConj) = conj(pdwm_omega);
                    end
                end
            end

            % Compute the last order
            mOrder = obj.order;
            nMIs = obj.ssm.MIs(mOrder).n;
            nUniqueMIs = ceil(nMIs / 2);
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
                    pR_omega(mOrder).coeffs(2) = conj(pRm_omega);
                end
            end

            % Compute the partial derivative of gamma
            pJ_omega = 1i/2 * (pR_omega(obj.order).coeffs(2) - pR_omega(obj.order).coeffs(1));

            % Ensure real values
            pJ_omega = real(pJ_omega);
            pCohom_omega = real(pCohom_omega);
        end

        function [pJ_phi, wAdj_pCohom_phi] = partial_derivative_phi(obj, wAdj)
            % Compute the partial derivative of the reduced backbone
            % coefficients gamma with respect to the mode shape phi.
            % Inputs:
            %   wAdj: the adjoint variables for the manifold coefficients.
            % Outputs:
            %   pJ_phi: the partial derivative of gamma with respect to
            %       phi.
            %   wAdj_pCohom_phi: the partial derivative of the
            %       cohomological equation with respect to phi
            %       pre-multiplied by wAdj.

            % Mode shape
            phi = obj.ssm.sys.phi;
            Mphi = obj.ssm.sys.M * phi;
            Cphi = obj.ssm.sys.C * phi;
        
            % Initialize partial derivatives
            wAdj_pCohom_phi = zeros(1, obj.ssm.sys.n);
            pR_phi = repmat(struct('coeffs', zeros(2, obj.ssm.sys.n)), 1, obj.order);
        
            % Loop over orders
            for mOrder = 2:(obj.order - 1)
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
                        end
        
                        % Compute the partial derivative of the cohomological equation
                        wAdj_pCohom_phi = wAdj_pCohom_phi - conj(wmAdj_phm_phi);
                    end
                end
            end

            % Compute the last order
            mOrder = obj.order;
            nMIs = obj.ssm.MIs(mOrder).n;
            nUniqueMIs = ceil(nMIs / 2);
            for mIdx = 1:nUniqueMIs
                % Current multi-index
                m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);
        
                % Eigenvalue coefficient
                LambdaM = obj.ssm.LambdaM(mOrder).coeffs(mIdx);
        
                % Auxiliary variables
                LMCphi = LambdaM * Mphi + Cphi;
        
                % Compute the nonlinear force contribution
                phi_pfm_phi = obj.partial_derivative_f_phi(mOrder, mIdx, phi);
        
                % Compute the V vector
                LMCphi_pVm_phi = obj.partial_derivative_V_phi(mOrder, mIdx, pR_phi, LMCphi);
        
                % Compute the Y vector
                Mphi_pYm_phi = obj.partial_derivative_Y_phi(mOrder, mIdx, pR_phi, Mphi);
        
                % Compute the C vector
                Cm = obj.ssm.C(mOrder).coeffs(:, mIdx);
                phi_pCm_phi = -Mphi_pYm_phi - LMCphi_pVm_phi - phi_pfm_phi;
        
                % Compute the reduced dynamics coefficients
                if m(1) - m(2) == 1
                    pRm_phi = (Cm.' + phi_pCm_phi) / ...
                        (LambdaM + obj.ssm.sys.Lambda(1) + obj.ssm.sys.delta);
                    pR_phi(mOrder).coeffs(1, :) = pRm_phi;
                    pR_phi(mOrder).coeffs(2, :) = conj(pRm_phi);
                end
            end

            % Compute the partial derivative of gamma
            pJ_phi = 1i/2 * (pR_phi(obj.order).coeffs(2, :) - pR_phi(obj.order).coeffs(1, :));
        
            % Ensure real values
            pJ_phi = real(pJ_phi);
            wAdj_pCohom_phi = real(wAdj_pCohom_phi);
        end

        function [pJ_mu, pCohom_mu] = partial_derivative_mu(obj, wAdj)
            % Compute the partial derivative of the reduced backbone
            % coefficients gamma with respect to the design variables mu.
            % Inputs:
            %   wAdj: the adjoint variables for the manifold coefficients.
            % Outputs:
            %   pJ_mu: the partial derivative of gamma with respect to mu.
            %   pCohom_mu: the partial derivative of the cohomological
            %       equation with respect to mu.

            % Initialize partial derivatives
            pR_mu = repmat(struct('pdv', []), 1, obj.order);
            pdw_mu = repmat(struct('pdv', []), 1, obj.order);
            pJ_mu = zeros(obj.ssm.sys.ndv, 1);
            pCohom_mu = zeros(obj.ssm.sys.ndv, 1);

            % Loop over orders
            for mOrder = 2:(obj.order - 1)
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

            % Compute the last order
            mOrder = obj.order;
            nMIs = obj.ssm.MIs(mOrder).n;
            nUniqueMIs = ceil(nMIs / 2);
            for mIdx = 1:nUniqueMIs
                % Current multi-index
                m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);

                % Eigenvalue coefficient
                LambdaM = obj.ssm.LambdaM(mOrder).coeffs(mIdx);

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
                        pR_mu(mOrder).pdv(dv).coeffs(2) = conj(pR_mu(mOrder).pdv(dv).coeffs(1));
                    end
                end
            end

            % Compute the partial derivative of gamma
            for dv = 1:obj.ssm.sys.ndv
                pJ_mu(dv) = 1i/2 * (pR_mu(obj.order).pdv(dv).coeffs(2) - pR_mu(obj.order).pdv(dv).coeffs(1));
            end

            % Ensure real values
            pJ_mu = real(pJ_mu);
            pCohom_mu = real(pCohom_mu);
        end

        function [pJ_mu, pCohom_mu] = partial_derivative_mu_topopt(obj, wAdj, iDOFs)
            % Compute the partial derivative of the reduced backbone
            % coefficients gamma with respect to the design variables mu
            % for topology optimization.
            % Inputs:
            %   wAdj: the adjoint variables for the manifold coefficients.
            %   iDOFs: the index of the DOFs of the current design
            %       variable.
            % Outputs:
            %   pJ_mu: the partial derivative of gamma with respect to mu.
            %   pCohom_mu: the partial derivative of the cohomological
            %       equation with respect to mu.

            % Initialize partial derivatives
            pR_mu = repmat(struct('pdv', []), 1, obj.order);
            pdw_mu = repmat(struct('pdv', []), 1, obj.order);
            % pJ_mu = 0;
            pCohom_mu = 0;

            % Extract element quantities
            phiE = obj.phiUnc(iDOFs);

            % Loop over orders
            for mOrder = 2:(obj.order - 1)
                % Number of multi-indices at this order
                nMIs = obj.ssm.MIs(mOrder).n;

                % Initialize coefficients
                pR_mu(mOrder).pdv = zeros(1, 2);
                pdw_mu(mOrder).pdv = zeros(obj.ssm.sys.n, nMIs);

                % Number of unique multi-indices
                nUniqueMIs = ceil(nMIs / 2);

                % Loop over multi-indices
                for mIdx = 1:nUniqueMIs
                    % Current multi-index
                    m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);

                    % Eigenvalue coefficient
                    LambdaM = obj.ssm.LambdaM(mOrder).coeffs(mIdx);

                    % Extract element quantities
                    wmE = obj.wUnc(mOrder).coeffs(iDOFs, mIdx);
                    wmAdjE = obj.wAdjUnc(mOrder).coeffs(iDOFs, mIdx);
        
                    % Extract coefficients
                    wmAdj = wAdj(mOrder).coeffs(:, mIdx);

                    % Compute the nonlinear force contribution
                    pfmE_mu = obj.partial_derivative_f_mu_topopt(mOrder, mIdx, iDOFs);

                    % Compute the V vector
                    VmE = obj.VUnc(mOrder).coeffs(iDOFs, mIdx);
                    pVm_mu = obj.partial_derivative_V_mu_topopt(mOrder, mIdx, pR_mu);

                    % Compute the Y vector
                    YmE = obj.YUnc(mOrder).coeffs(iDOFs, mIdx);
                    pYm_mu = obj.partial_derivative_Y_mu_topopt(mOrder, mIdx, pR_mu, pdw_mu);

                    % Compute the partial derivative of C
                    % The element and global parts are computed separately
                    pCm_mu = -(LambdaM * obj.ssm.sys.M + obj.ssm.sys.C) * pVm_mu - obj.ssm.sys.M * pYm_mu;
                    pCmE_mu = -(LambdaM * obj.ssm.sys.Me + obj.ssm.sys.Ce) * VmE - ...
                        obj.ssm.sys.Me * YmE - pfmE_mu;

                    % Compute the partial derivative of R
                    if m(1) - m(2) == 1
                        % Compute the partial derivative of R
                        pR_mu(mOrder).pdv(1) = (obj.ssm.sys.phi.' * pCm_mu + phiE.' * pCmE_mu) / ...
                            (LambdaM + obj.ssm.sys.Lambda(1) + obj.ssm.sys.delta);
                    end
                    
                    % Compute the vector h
                    phm_mu = pCm_mu;
                    phmE_mu = pCmE_mu;
                    if m(1) - m(2) == 1
                        % Compute vector D
                        Dm = -((LambdaM + obj.ssm.sys.Lambda(1)) * obj.ssm.sys.M + obj.ssm.sys.C) * obj.ssm.sys.phi;
        
                        % Partial derivative of D
                        pDmE_mu = -((LambdaM + obj.ssm.sys.Lambda(1)) * obj.ssm.sys.Me + ...
                            obj.ssm.sys.Ce) * phiE;
                        
                        % Compute the partial derivative
                        phm_mu = phm_mu + pR_mu(mOrder).pdv(1) * Dm;
                        phmE_mu = phmE_mu + obj.ssm.R(mOrder).coeffs(1) * pDmE_mu;
                    end
        
                    % Compute the partial derivative of L
                    pLmE_mu = obj.ssm.sys.Ke + LambdaM * obj.ssm.sys.Ce + ...
                        LambdaM^2 * obj.ssm.sys.Me;
        
                    % Compute the partial derivative of the cohomological equation
                    cohomAdj = wmAdjE.' * (pLmE_mu * wmE - phmE_mu) - wmAdj.' * phm_mu;
                    pCohom_mu = pCohom_mu + cohomAdj;
        
                    % Initialize with the partial derivative of the V vector
                    pdwm_mu = pVm_mu;
        
                    % Add the contribution from the reduced dynamics coefficients
                    if m(1) - m(2) == 1
                        pdwm_mu = pdwm_mu + ...
                            obj.ssm.sys.phi * pR_mu(mOrder).pdv(1);
                    end
        
                    % Store the partial derivative of the velocity coefficients
                    pdw_mu(mOrder).pdv(:, mIdx) = pdwm_mu;
        
                    % Compute the symmetric multi-index and its position index
                    mConj = [m(2); m(1)];
                    mIdxConj = from_MI_to_position_idx(mConj);
        
                    % Compute the coefficients for the symmetric multi-index
                    if mIdx ~= mIdxConj
                        if m(1) - m(2) == 1
                            % Compute the partial derivative of R
                            pR_mu(mOrder).pdv(2) = conj(pR_mu(mOrder).pdv(1));
                        end
        
                        % Compute the partial derivative of the cohomological equation
                        pCohom_mu = pCohom_mu + conj(cohomAdj);
        
                        % Store the partial derivative of the velocity coefficients
                        pdw_mu(mOrder).pdv(:, mIdxConj) = conj(pdwm_mu);
                    end
                end
            end

            % Compute the last order
            mOrder = obj.order;
            nMIs = obj.ssm.MIs(mOrder).n;
            nUniqueMIs = ceil(nMIs / 2);
            for mIdx = 1:nUniqueMIs
                % Current multi-index
                m = obj.ssm.MIs(mOrder).coeffs(:, mIdx);

                % Eigenvalue coefficient
                LambdaM = obj.ssm.LambdaM(mOrder).coeffs(mIdx);

                % Compute the nonlinear force contribution
                pfmE_mu = obj.partial_derivative_f_mu_topopt(mOrder, mIdx, iDOFs);

                % Compute the V vector
                VmE = obj.VUnc(mOrder).coeffs(iDOFs, mIdx);
                pVm_mu = obj.partial_derivative_V_mu_topopt(mOrder, mIdx, pR_mu);

                % Compute the Y vector
                YmE = obj.YUnc(mOrder).coeffs(iDOFs, mIdx);
                pYm_mu = obj.partial_derivative_Y_mu_topopt(mOrder, mIdx, pR_mu, pdw_mu);

                % Compute the partial derivative of C
                % The element and global parts are computed separately
                pCm_mu = -(LambdaM * obj.ssm.sys.M + obj.ssm.sys.C) * pVm_mu - obj.ssm.sys.M * pYm_mu;
                pCmE_mu = -(LambdaM * obj.ssm.sys.Me + obj.ssm.sys.Ce) * VmE - ...
                    obj.ssm.sys.Me * YmE - pfmE_mu;

                % Compute the partial derivative of R
                if m(1) - m(2) == 1
                    % Compute the partial derivative of R
                    pRm_mu = (obj.ssm.sys.phi.' * pCm_mu + phiE.' * pCmE_mu) / ...
                        (LambdaM + obj.ssm.sys.Lambda(1) + obj.ssm.sys.delta);
                    
                    % Store the partial derivative of R
                    pR_mu(mOrder).pdv(1) = pRm_mu;
                    pR_mu(mOrder).pdv(2) = conj(pRm_mu);
                end
            end

            % Compute the partial derivative of gamma
            pJ_mu = 1i/2 * (pR_mu(obj.order).pdv(2) - pR_mu(obj.order).pdv(1));

            % Ensure real values
            pJ_mu = real(pJ_mu);
            pCohom_mu = real(pCohom_mu);
        end
    end
end

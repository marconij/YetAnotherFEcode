classdef SensitivityBackbone < MechSSMsensitivity
    % SensitivityBackbone is a subclass of MechSSMsensitivity that is
    % specifically designed to handle the sensitivity analysis of the
    % backbone curve.
    % Methods:
    %   SensitivityBackbone: the constructor.
    %   adjoint_sensitivity_Omega: compute the adjoint sensitivity of the
    %       response frequency Omega.
    %   adjoint_sensitivity_Omega_topopt: compute the adjoint sensitivity
    %       of the response frequency Omega for topology optimization.
    %   adjoint_equation_rho: solve the adjoint equation of the reduced
    %       amplitude rho.
    %   adjoint_equation_w: solve the adjoint equation of the manifold
    %       coefficients w.
    %   adjoint_equation_omega_phi: solve the adjoint equation of the
    %       natural frequency omega and the mode shape phi.
    %   partial_derivative_Omega_rho: compute the partial derivative of the
    %       response frequency Omega with respect to the reduced amplitude
    %       rho.
    %   partial_derivative_Omega_w: compute the partial derivative of the
    %       response frequency Omega with respect to the manifold
    %       coefficients w.
    %   partial_derivative_Omega_omega: compute the partial derivative of
    %       the response frequency Omega with respect to the natural 
    %       requency omega.
    %   partial_derivative_Omega_phi: compute the partial derivative of the
    %       response frequency Omega with respect to the mode shape phi.
    %   partial_derivative_Omega_mu: compute the partial derivative of the
    %       response frequency Omega with respect to the design variables
    %       mu.
    %   partial_derivative_Omega_mu_topopt: compute the partial derivative
    %       of the response frequency Omega with respect to the design
    %       variables mu for topology optimization.

    methods
        function obj = SensitivityBackbone(ssm)
            % Construct an instance of this class.
            % Inputs:
            %   ssm: the MechSSM object.
            % Outputs:
            %   obj: an instance of this class.

            % Call the constructor of the superclass
            obj = obj@MechSSMsensitivity(ssm);
        end

        function dOmega = adjoint_sensitivity_Omega(obj, rho, xk, x, outDof)
            % Compute the adjoint sensitivity of the response frequency
            % Omega.
            % Inputs:
            %   rho: the reduced amplitude. It must be a scalar.
            %   xk: the physical amplitudes at all the time steps.
            %   x: the rms of the physical amplitude.
            %   outDof: the index of the output dof.
            % Outputs:
            %   dOmega: the adjoint sensitivity of Omega.
        
            % Initialize the adjoint sensitivity of Omega
            dOmega = zeros(obj.ssm.sys.ndv, 1);
        
            % Adjoint equation for the reduced amplitude rho
            rhoAdj = obj.adjoint_equation_rho(rho, xk, x, outDof);
        
            % Adjoint equation for the manifold coefficients w
            wAdj = obj.adjoint_equation_w(rho, xk, x, outDof, rhoAdj);
        
            % Adjoint equation for the natural frequency omega and the mode shape phi
            [omegaAdj, phiAdj] = obj.adjoint_equation_omega_phi(rho, xk, x, outDof, rhoAdj, wAdj);
        
            % Partial derivatives of Omega and the cohomological equation
            [pOmega_mu, pCohom_mu] = obj.partial_derivative_Omega_mu(rho, wAdj);
        
            % Loop over the design variables
            for dv = 1:obj.ssm.sys.ndv
                % Add the contribution of Omega
                dOmega(dv) = pOmega_mu(dv) ...
                    + phiAdj.' * (obj.ssm.sys.pK(dv).coeffs - obj.ssm.sys.omega^2 * obj.ssm.sys.pM(dv).coeffs) * obj.ssm.sys.phi ...
                    + omegaAdj * obj.ssm.sys.phi.' * obj.ssm.sys.pM(dv).coeffs * obj.ssm.sys.phi ...
                    + pCohom_mu(dv);
            end
        end

        function dOmega = adjoint_sensitivity_Omega_topopt(obj, rho, xk, x, outDof)
            % Compute the adjoint sensitivity of the response frequency
            % Omega for a topology optimization problem.
            % Inputs:
            %   rho: the reduced amplitude. It must be a scalar.
            %   xk: the physical amplitudes at all the time steps.
            %   x: the rms of the physical amplitude.
            %   outDof: the index of the output dof.
            % Outputs:
            %   dOmega: the adjoint sensitivity of Omega.
        
            % Extract FE assembly
            yafecAssembly = obj.ssm.sys.yafecAssembly;
            nDOFs = yafecAssembly.Mesh.nDOFs;
        
            % Initialize the adjoint sensitivity of Omega
            dOmega = zeros(obj.ssm.sys.ndv, 1);
        
            % Adjoint equation for the reduced amplitude rho
            rhoAdj = obj.adjoint_equation_rho(rho, xk, x, outDof);
        
            % Adjoint equation for the manifold coefficients w
            wAdj = obj.adjoint_equation_w(rho, xk, x, outDof, rhoAdj);
        
            % Adjoint equation for the natural frequency omega and the mode shape phi
            [omegaAdj, phiAdj] = obj.adjoint_equation_omega_phi(rho, xk, x, outDof, rhoAdj, wAdj);
        
            % Convert vectors into their unconstrained form
            obj.phiUnc = yafecAssembly.unconstrain_vector(obj.ssm.sys.phi);
            obj.phiAdjUnc = yafecAssembly.unconstrain_vector(phiAdj);
            
            obj.wUnc = repmat(struct('coeffs', []), 1, obj.ssm.maxOrder);
            obj.wAdjUnc = repmat(struct('coeffs', []), 1, obj.ssm.maxOrder);
            obj.VUnc = repmat(struct('coeffs', []), 1, obj.ssm.maxOrder);
            obj.YUnc = repmat(struct('coeffs', []), 1, obj.ssm.maxOrder);
        
            obj.wUnc(1).coeffs = [obj.phiUnc, obj.phiUnc];
            for mOrder = 2:obj.ssm.maxOrder
                % Number of multi-indices at this order
                nMIs = obj.ssm.MIs(mOrder).n;
        
                % Initialize coefficients
                obj.wUnc(mOrder).coeffs = zeros(nDOFs, nMIs);
                obj.wAdjUnc(mOrder).coeffs = zeros(nDOFs, nMIs);
                obj.VUnc(mOrder).coeffs = zeros(nDOFs, nMIs);
                obj.YUnc(mOrder).coeffs = zeros(nDOFs, nMIs);
        
                % Loop over multi-indices
                for mIdx = 1:nMIs
                    obj.wUnc(mOrder).coeffs(:, mIdx) = yafecAssembly.unconstrain_vector(obj.ssm.w(mOrder).coeffs(:, mIdx));
                    obj.wAdjUnc(mOrder).coeffs(:, mIdx) = yafecAssembly.unconstrain_vector(wAdj(mOrder).coeffs(:, mIdx));
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
        
                % Partial derivatives of Omega and the cohomological equation
                [pOmega_mu, pCohom_mu] = obj.partial_derivative_Omega_mu_topopt(rho, wAdj, iDOFs);
        
                % Assemble the Lagrangian
                dOmega(dv) = pOmega_mu ...
                    + phiAdjE.' * (obj.ssm.sys.Ke - obj.ssm.sys.omega^2 * obj.ssm.sys.Me) * phiE ...
                    + omegaAdj * phiE.' * obj.ssm.sys.Me * phiE ...
                    + pCohom_mu;
            end
        end

        function rhoAdj = adjoint_equation_rho(obj, rho, xk, x, outDof)
            % Solve the adjoint equation for the reduced amplitude rho.
            % Inputs:
            %   rho: the reduced amplitude. It must be a scalar.
            %   xk: the physical amplitudes at all the time steps.
            %   x: the rms of the physical amplitude.
            %   outDof: the index of the output dof.
            % Outputs:
            %   rhoAdj: the adjoint variable for the reduced amplitude.
        
            % Check that rho is a scalar
            if ~isscalar(rho)
                error('The computation of the partial derivatives is only supported for scalar reduced amplitudes.');
            end
        
            % Partial derivative of Omega with respect to rho
            pOmega_rho = obj.partial_derivative_Omega_rho(rho);
        
            % Partial derivative of x with respect to rho
            px_rho = obj.partial_derivative_x_rho(rho, xk, x, outDof);
        
            % Solve adjoint equation for rho
            rhoAdj = -pOmega_rho / px_rho;
        end

        function wAdj = adjoint_equation_w(obj, rho, xk, x, outDof, rhoAdj)
            % Solve the adjoint equation for the manifold coefficients w.
            % Inputs:
            %   rho: the reduced amplitude.
            %   xk: the physical amplitudes at all the time steps.
            %   x: the rms of the physical amplitude.
            %   outDof: the index of the output dof.
            %   rhoAdj: the adjoint variable for the reduced amplitude.
            % Outputs:
            %   wAdj: the adjoint variables for w.
        
            % Initialize adjoint variables
            wAdj = repmat(struct('coeffs', []) , 1, obj.ssm.maxOrder);
        
            % Compute partial derivative of x
            px_w = obj.partial_derivative_x_w(rho, xk, x, outDof);
        
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
        
                    % Partial derivatives of Omega and the cohomological equations
                    [pOmega_wm, pCohom_wm] = obj.partial_derivative_Omega_w(mOrder, mIdx, rho, wAdj);
        
                    % Cohomological matrix
                    LambdaM = obj.ssm.LambdaM(mOrder).coeffs(mIdx);
                    Lm = obj.ssm.sys.K + LambdaM * obj.ssm.sys.C + LambdaM^2 * obj.ssm.sys.M;
        
                    % Solve the adjoint equation for the manifold coefficients w
                    px_wm = px_w(mOrder).coeffs(:, mIdx);
                    rhs = -pCohom_wm.' - pOmega_wm.' - rhoAdj * px_wm;
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

        function [omegaAdj, phiAdj] = adjoint_equation_omega_phi(obj, rho, xk, x, outDof, rhoAdj, wAdj)
            % Solve together the adjoint equations for the natural
            % frequency omega and the mode shape phi.
            % Inputs:
            %   rho: the reduced amplitude.
            %   xk: the physical amplitudes at all the time steps.
            %   x: the rms of the physical amplitude.
            %   outDof: the index of the output dof.
            %   rhoAdj: the adjoint variable for the reduced amplitude.
            %   wAdj: the adjoint variables for the manifold coefficients.
            % Outputs:
            %   omegaAdj: the adjoint variable for omega.
            %   phiAdj: the adjoint variable for phi.
        
            % Create matrix A
            A = [obj.ssm.sys.K - obj.ssm.sys.omega^2 * obj.ssm.sys.M, 2 * obj.ssm.sys.M * obj.ssm.sys.phi; ...
                 2 * obj.ssm.sys.phi.' * obj.ssm.sys.M, 0];
        
            % Compute partial derivative of x
            px_phi = obj.partial_derivative_x_phi(rho, xk, x, outDof);
            
            % Partial derivatives of Omega and the cohomological equations
            [pOmega_omega, pCohom_omega] = obj.partial_derivative_Omega_omega(rho, wAdj);
            [pOmega_phi, wAdj_pCohom_phi] = obj.partial_derivative_Omega_phi(rho, wAdj);
        
            % Initialize vector b
            b = [-pOmega_phi.' - wAdj_pCohom_phi.' - rhoAdj * px_phi; ...
                 (pOmega_omega + pCohom_omega) / obj.ssm.sys.omega];
        
            % Solve the system
            x = A \ b;
        
            % Extract the adjoint variables
            phiAdj = x(1:end-1);
            omegaAdj = x(end);
        end
        
        function pOmega_rho = partial_derivative_Omega_rho(obj, rho)
            % Compute the partial derivative of the response frequency
            % Omega with respect to the reduced amplitude rho.
            % Inputs:
            %   rho: the reduced amplitude.
            % Outputs:
            %   pOmega_rho: the partial derivative of Omega with respect to
            %       rho.
        
            % Initialize partial derivatives
            pOmega_rho = 0;
        
            % Loop over odd orders
            for mOrder = 3:2:obj.ssm.maxOrder
                % Extract coefficients
                R1 = obj.ssm.R(mOrder).coeffs(1); % 1st non-zero coefficient (m(1) - m(2) =  1)
                R2 = obj.ssm.R(mOrder).coeffs(2); % 2nd non-zero coefficient (m(1) - m(2) = -1)
        
                % Partial derivative with respect to rho
                pOmega_rho = pOmega_rho + 1i/2 * (R2 - R1) * (mOrder - 1) * rho.^(mOrder - 2);
            end
        end

        function [pOmega_wm, pCohom_wm] = partial_derivative_Omega_w(obj, mOrder, mIdx, rho, wAdj)
            % Compute the partial derivative of the response frequency
            % Omega with respect to the manifold coefficient wm.
            % Inputs:
            %   mOrder: the order of the multi-index.
            %   mIdx: the index of the multi-index in the multi-index
            %       array.
            %   rho: the reduced amplitude.
            %   wAdj: the adjoint variables for w.
            % Outputs:
            %   pOmega_wm: the partial derivative of Omega with respect to
            %       wm.
            %   pCohom_wm: the partial derivative of the cohomological
            %       equation with respect to wm.
        
            % Mode shape
            phi = obj.ssm.sys.phi;
            Mphi = obj.ssm.sys.M * phi;
            Cphi = obj.ssm.sys.C * phi;
            
            % Initialize partial derivatives
            pOmega_wm = zeros(1, obj.ssm.sys.n);
            pCohom_wm = zeros(1, obj.ssm.sys.n);
        
            % Current multi-index
            pR_wm = repmat(struct('coeffs', zeros(2, obj.ssm.sys.n)), 1, obj.ssm.maxOrder);
        
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
                        pOmega_wm = pOmega_wm - 1i/2 * pRq_wm * rho^(qOrder - 1);
                    elseif q(1) - q(2) == -1
                        % Compute the partial derivative of R
                        pRq_wm = phi_pCq_wm / ...
                            (LambdaQ + obj.ssm.sys.Lambda(2) + obj.ssm.sys.delta);
                        pR_wm(qOrder).coeffs(2, :) = pRq_wm;
        
                        % Add the contribution to the response frequency
                        pOmega_wm = pOmega_wm + 1i/2 * pRq_wm * rho^(qOrder - 1);
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
        end

        function [pOmega_omega, pCohom_omega] = partial_derivative_Omega_omega(obj, rho, wAdj)
            % Compute the partial derivative of the natural frequency Omega
            % with respect to the reduced frequency omega.
            % Inputs:
            %   rho: the reduced amplitude.
            %   wAdj: the adjoint variables for the manifold coefficients.
            % Outputs:
            %   pOmega_omega: the partial derivative of Omega with respect
            %       to omega.
            %   pCohom_omega: the partial derivative of the cohomological
            %       equation with respect to omega.
        
            % Initialize partial derivatives
            pR_omega = repmat(struct('coeffs', []), 1, obj.ssm.maxOrder);
            pdw_omega = repmat(struct('coeffs', []), 1, obj.ssm.maxOrder);
            pCohom_omega = 0;
        
            % Initialize partial derivative of Omega
            pOmega_omega = 1i/2 * (obj.ssm.sys.pLambda_omega(2) - obj.ssm.sys.pLambda_omega(1));
        
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
        
                        % Add the contribution to the partial derivative of Omega
                        pOmega_omega = pOmega_omega - 1i/2 * pRm_omega * rho^(mOrder - 1);
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
            
                            % Add the contribution to the partial derivative of Omega
                            pOmega_omega = pOmega_omega + 1i/2 * conj(pRm_omega) * rho^(mOrder - 1);
                        end
        
                        % Compute the partial derivative of the cohomological equation
                        pCohom_omega = pCohom_omega + conj(cohomAdj);
        
                        % Store the partial derivative of the velocity coefficients
                        pdw_omega(mOrder).pdv(:, mIdxConj) = conj(pdwm_omega);
                    end
                end
            end
        
            % Ensure real values
            pOmega_omega = real(pOmega_omega);
            pCohom_omega = real(pCohom_omega);
        end

        function [pOmega_phi, wAdj_pCohom_phi] = partial_derivative_Omega_phi(obj, rho, wAdj)
            % Compute the partial derivative of the response frequency
            % Omega with respect to the mode shape phi.
            % Inputs:
            %   rho: the reduced amplitude.
            %   wAdj: the adjoint variables for the manifold coefficients.
            % Outputs:
            %   pOmega_phi: the partial derivative of Omega with respect to
            %       phi.
            %   wAdj_pCohom_phi: the partial derivative of the
            %       cohomological equations with respect to phi
            %       premultiplied by wAdj.
        
            % Mode shape
            phi = obj.ssm.sys.phi;
            Mphi = obj.ssm.sys.M * phi;
            Cphi = obj.ssm.sys.C * phi;
        
            % Initialize
            pOmega_phi = zeros(1, obj.ssm.sys.n);
            wAdj_pCohom_phi = zeros(1, obj.ssm.sys.n);
            pR_phi = repmat(struct('coeffs', zeros(2, obj.ssm.sys.n)), 1, obj.ssm.maxOrder);
        
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
                        pOmega_phi = pOmega_phi - 1i/2 * pRm_phi * rho^(mOrder - 1);
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
            
                            % Add the contribution to the partial derivative of Omega
                            pOmega_phi = pOmega_phi + 1i/2 * conj(pRm_phi) * rho^(mOrder - 1);
                        end
        
                        % Compute the partial derivative of the cohomological equation
                        wAdj_pCohom_phi = wAdj_pCohom_phi - conj(wmAdj_phm_phi);
                    end
                end
            end
        
            % Ensure real values
            pOmega_phi = real(pOmega_phi);
            wAdj_pCohom_phi = real(wAdj_pCohom_phi);
        end

        function [pOmega_mu, pCohom_mu] = partial_derivative_Omega_mu(obj, rho, wAdj)
            % Compute the partial derivative of the response frequency
            % Omega with respect to the design variables mu.
            % Inputs:
            %   rho: the reduced amplitude.
            %   wAdj: the adjoint variables for the manifold coefficients.
            % Outputs:
            %   pOmega_mu: the partial derivative of Omega with respect to
            %       mu.
            %   pCohom_mu: the partial derivative of the cohomological
            %       equation with respect to mu.
        
            % Initialize partial derivatives
            pR_mu = repmat(struct('pdv', []), 1, obj.ssm.maxOrder);
            pdw_mu = repmat(struct('pdv', []), 1, obj.ssm.maxOrder);
            pCohom_mu = zeros(obj.ssm.sys.ndv, 1);
            pOmega_mu = zeros(obj.ssm.sys.ndv, 1);
        
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
                            
                            % Add the contribution to the partial derivative of Omega
                            pOmega_mu(dv) = pOmega_mu(dv) - 1i/2 * pR_mu(mOrder).pdv(dv).coeffs(1) * rho^(mOrder - 1);
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
        
                                % Add the contribution to the partial derivative of Omega
                                pOmega_mu(dv) = pOmega_mu(dv) ...
                                    + 1i/2 * pR_mu(mOrder).pdv(dv).coeffs(2) * rho^(mOrder - 1);
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
            pOmega_mu = real(pOmega_mu);
            pCohom_mu = real(pCohom_mu);
        end

        function [pOmega_mu, pCohom_mu] = partial_derivative_Omega_mu_topopt(obj, rho, wAdj, iDOFs)
            % Compute the partial derivative of the response frequency
            % Omega with respect to the design variables mu for a topology
            % optimization problem.
            % Inputs:
            %   rho: the reduced amplitude.
            %   wAdj: the adjoint variables for the manifold coefficients.
            %   iDOFs: the index of the DOF of the current design variable.
            % Outputs:
            %   pOmega_mu: the partial derivative of Omega with respect to
            %       mu.
            %   pCohom_mu: the partial derivative of the cohomological
            %       equation with respect to mu.
        
            % Initialize partial derivatives
            pR_mu = repmat(struct('pdv', []), 1, obj.ssm.maxOrder);
            pdw_mu = repmat(struct('pdv', []), 1, obj.ssm.maxOrder);
            pCohom_mu = 0;
            pOmega_mu = 0;
        
            % Extract element quantities
            phiE = obj.phiUnc(iDOFs);
        
            % Loop over orders
            for mOrder = 2:obj.ssm.maxOrder
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
        
                        % Add the contribution to the partial derivative of Omega
                        pOmega_mu = pOmega_mu - 1i/2 * pR_mu(mOrder).pdv(1) * rho^(mOrder - 1);
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
            
                            % Add the contribution to the partial derivative of Omega
                            pOmega_mu = pOmega_mu + 1i/2 * pR_mu(mOrder).pdv(2) * rho^(mOrder - 1);
                        end
        
                        % Compute the partial derivative of the cohomological equation
                        pCohom_mu = pCohom_mu + conj(cohomAdj);
        
                        % Store the partial derivative of the velocity coefficients
                        pdw_mu(mOrder).pdv(:, mIdxConj) = conj(pdwm_mu);
                    end
                end
            end
        
            % Ensure real values
            pOmega_mu = real(pOmega_mu);
            pCohom_mu = real(pCohom_mu);
        end
    end
end

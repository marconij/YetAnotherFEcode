classdef MechSSMsensitivity < handle
    % A class to compute the adjoint sensitivity of the SSM coefficients.
    % Properties:
    %   ssm: the MechSSM object.
    %   phiUnc: the mode shape (unconstrained).
    %   phiAdjUnc: the adjoint of the mode shape (unconstrained).
    %   wUnc: the manifold coefficients (unconstrained).
    %   wAdjUnc: the adjoint of the manifold coefficients (unconstrained).
    %   VUnc: the vector V (unconstrained).
    %   YUnc: the vector Y (unconstrained).
    % Methods:
    %   MechSSMsensitivity: the constructor.
    %   partial_derivative_x_rho: compute the partial derivative of the
    %       physical amplitude x with respect to the reduced amplitude rho.
    %   partial_derivative_x_w: compute the partial derivative of the
    %       physical amplitude x with respect to the manifold coefficients
    %       w.
    %   partial_derivative_x_phi: compute the partial derivative of the
    %       physical amplitude x with respect to the mode shape phi.
    %   partial_derivative_f_w: compute the partial derivative of the 
    %       nonlinear force f with respect to the manifold coefficients w.
    %   partial_derivative_f_phi: compute the partial derivative of the 
    %       nonlinear force f with respect to the mode shape phi.
    %   partial_derivative_f_mu: compute the partial derivative of the 
    %       nonlinear force f with respect to the design variables mu.
    %   partial_derivative_V_w: compute the partial derivative of the 
    %       vector V with respect to the manifold coefficients w.
    %   partial_derivative_V_omega: compute the partial derivative of the 
    %       vector V with respect to the natural frequency omega.
    %   partial_derivative_V_phi: compute the partial derivative of the 
    %       vector V with respect to the mode shape phi.
    %   partial_derivative_V_mu: compute the partial derivative of the 
    %       vector V with respect to the design variables mu.
    %   partial_derivative_Y_w: compute the partial derivative of the 
    %       vector Y with respect to the manifold coefficients w.
    %   partial_derivative_Y_omega: compute the partial derivative of the 
    %       vector Y with respect to the natural frequency omega.
    %   partial_derivative_Y_phi: compute the partial derivative of the 
    %       vector Y with respect to the mode shape phi.
    %   partial_derivative_Y_mu: compute the partial derivative of the 
    %       vector Y with respect to the design variables mu.
    %   partial_derivative_dw_w: compute the partial derivative of the 
    %       velocity coefficients dw with respect to the manifold
    %       coefficients w.
    %   partial_derivative_dw_phi: compute the partial derivative of the 
    %       velocity coefficients dw with respect to the mode shape phi.
    %   partial_derivative_f_mu_topopt: compute the partial derivative of
    %       the nonlinear force f with respect to the design variables mu
    %       for topology optimization.
    %   partial_derivative_V_mu_topopt: compute the partial derivative of
    %       the vector V with respect to the design variables mu for
    %       topology optimization.
    %   partial_derivative_Y_mu_topopt: compute the partial derivative of
    %       the vector Y with respect to the design variables mu for
    %       topology optimization.
    
    properties
        ssm
        phiUnc, phiAdjUnc
        wUnc, wAdjUnc
        VUnc, YUnc
    end
    
    methods
        function obj = MechSSMsensitivity(ssm)
            % Construct an instance of this class.
            % Inputs:
            %   ssm: the MechSSM object.
            % Outputs:
            %   obj: the MechSSMsenstivity object.
            
            % Set the SSM object
            obj.ssm = ssm;

            % Check that the partial derivatives of the mechanical system have been computed
            % pM is used for parametric optimization problems, while Me is used for topology optimization problems
            if isempty(obj.ssm.sys.pM) && isempty(obj.ssm.sys.Me)
                error(['The partial derivatives of the mechanical system have not been computed. ' ...
                       'Please set the `computePartialDerivatives` flag to true when building the mechanical system.']);
            end

            % Check that the SSM coefficients have been stored
            if isempty(obj.ssm.LambdaM)
                error(['The SSM coefficients have not been stored. ', ...
                    'Please set the `storeCoefficients` flag to true when computing the SSM.']);
            end
        end
    end
end

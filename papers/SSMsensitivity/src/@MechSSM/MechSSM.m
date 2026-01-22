classdef MechSSM < handle
    % A class for computing 2D autonomous SSMs for lightly-damped
    % mechanical systems.
    % Properties:
    %   sys: the MechSystem object.
    %   maxOrder: the maximum order of the SSM.
    %   nTheta: number of time instants (default is 2^7).
    %   nTimeSteps: number of time steps for the RMS computation.
    %   theta: the time instants (row vectoer of phases [rad]).
    %   MIs: the structure of multi-indices.
    %   R: the structure of the reduced dynamics coefficients.
    %   w: the structure of the manifold coefficients.
    %   dw: the structure of the manifold velocity coefficients.
    %   LambdaM: the structure of the eigenvalue coefficients.
    %   f: the structure of the nonlinear force contributions.
    %   V: the structure of the V vectors.
    %   Y: the structure of the Y vectors.
    %   C: the structure of the C vectors.
    %   h: the structure of the h vectors (right-hand side of the
    %       cohomological equation).
    %   storeCoefficients: flag to store the SSM coefficients (default is
    %       false).
    %   computeSensitivity: flag to enable the sensitivity computation
    %       using the direct differentiation (default is false).
    % Methods:
    %   MechSSM: the constructor.
    %   compute_ssm: compute the SSM coefficients.
    %   compute_lambda: compute the eigenvalue coefficients at the
    %       multi-index m.
    %   compute_f: compute the nonlinear force contributions at the
    %       multi-index m.
    %   compute_V: compute the V vector at the multi-index m.
    %   compute_Y: compute the Y vector at the multi-index m.
    %   compute_R: compute the reduced dynamics coefficients at the
    %       multi-index m.
    %   compute_dw: compute the manifold velocity coefficients at the
    %       multi-index m.
    %   compute_x: compute the RMS of the physical amplitude over the time
    %       instants.
    %   compute_rho_from_x: compute the reduced amplitude from the RMS of
    %       the physical amplitude.
    %   compute_rho_from_x_leading_order: compute the reduced amplitude
    %       from the RMS of the physical amplitude at the leading order.
    %   compute_omega: compute the response frequency.
    %   compute_rho_derivative: compute the derivative of the reduced
    %       amplitude.
    %   compute_omega_derivative: compute the derivative of the response
    %       frequency.
    %   compute_dp_t: compute the time derivative of the reduced
    %       coordinates.
    %   compute_W: compute the 1st order form manifold parametrization:
    %       W = [w; dw].
    %   compute_dW_t: compute the time derivative of the 1st order form
    %       manifold parametrization: DW = [dw; ddw].
    %   compute_invariance_residual: compute the residual of the invariance
    %       equation.
    %   error_measure: compute the error measure from the invariance
    %       equation.

    properties
        sys
        maxOrder
        nTheta, nTimeSteps, theta
        MIs
        R, w, dw
        LambdaM, f, V, Y, C, h
        storeCoefficients
        computeSensitivity
    end

    methods
        function obj = MechSSM(sys, varargin)
            % Construct an instance of this class.
            % Inputs:
            %   sys: the MechSystem object.
            %   nTheta: number of time instants (optional, default is 2^7).
            %   rmsStyle: flag for the RMS computation style (optional,
            %       default is 'SSMTool').
            % Outputs:
            %   obj: the MechSSM object.

            % Parse the inputs
            p = inputParser;
            addOptional(p, 'nTheta', 2^7, @isnumeric);
            addOptional(p, 'rmsStyle', 'SSMTool', @(x) any(validatestring(x, {'RMS', 'SSMTool'})));
            parse(p, varargin{:});
            obj.nTheta = p.Results.nTheta;

            % Check the RMS style
            if strcmp(p.Results.rmsStyle, 'RMS')
                obj.nTimeSteps = obj.nTheta; % classic RMS style
            elseif strcmp(p.Results.rmsStyle, 'SSMTool')
                obj.nTimeSteps = obj.nTheta - 1; % SSMTool style
            else
                error('Invalid RMS style.');
            end

            % Initialize the mechanical system
            obj.sys = sys;

            % Initialize the time instants (row vector)
            obj.theta = linspace(0, 2*pi, obj.nTheta);
            % obj.theta(end) = []; % remove the last element (2*pi) since it is equal to the first element (0)
        end
    end
end

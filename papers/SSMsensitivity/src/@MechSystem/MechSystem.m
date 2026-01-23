classdef MechSystem < handle
    % Model of a mechanical system with quadratic and cubic geometric
    % nonlinearities.
    % Properties:
    %   n: number of constrained DOFs.
    %   yafecAssembly: yafec assembly.
    %   M: mass matrix.
    %   C: damping matrix.
    %   K: stiffness matrix.
    %   T2: quadratic geometric stiffness tensor.
    %   T3: cubic geometric stiffness tensor.
    %   D2: quadratic damping tensor.
    %   D3: cubic damping tensor.
    %   damping_type: type of damping (0: no damping, 1: Rayleigh damping,
    %       2: Q-factor damping).
    %   rayleigh: Rayleigh damping mass coefficients.
    %   Qfactor: Q-factor for the damping.
    %   omega: natural frequency.
    %   phi: mode shape.
    %   csi: damping ratio.
    %   delta: modal damping.
    %   omegaDamped: damped natural frequency.
    %   lambda: eigenvalue of the 1st order system.
    %   tau: mass normalization coefficient.
    %   psi: left mode shape.
    %   Lambda: eigenvalues of the 1st order system.
    %   Phi: mode shapes.
    %   Psi: left mode shapes.
    %   pcsi_omega: partial derivative of the damping ratio wrt the natural
    %       frequency.
    %   pdelta_omega: partial derivative of the modal damping wrt the
    %       natural frequency.
    %   plambda_omega: partial derivative of the eigenvalue wrt the natural
    %       frequency.
    %   pLambda_omega: partial derivative of the eigenvalues wrt the
    %       natural frequency.
    %   ndv: number of design variables.
    %   pM: partial derivative of the mass matrix.
    %   pC: partial derivative of the damping matrix.
    %   pK: partial derivative of the stiffness matrix.
    %   pT2: partial derivative of the quadratic geometric stiffness
    %       tensor.
    %   pT3: partial derivative of the cubic geometric stiffness tensor.
    %   pD2: partial derivative of the quadratic damping tensor.
    %   pD3: partial derivative of the cubic damping tensor.
    %   pC_omega: partial derivative of the damping matrix wrt the natural
    %       frequency.
    %   dC: total derivative of the damping matrix.
    %   Me: element mass matrix.
    %   Ce: element damping matrix.
    %   Ke: element stiffness matrix.
    %   T2e: element quadratic geometric stiffness tensor.
    %   T3e: element cubic geometric stiffness tensor.
    %   D2e: element quadratic damping tensor.
    %   D3e: element cubic damping tensor.
    %   domega: total derivative of the natural frequency.
    %   dphi: total derivative of the mode shape.
    %   dcsi: total derivative of the damping ratio.
    %   ddelta: total derivative of the modal damping.
    %   domegaDamped: total derivative of the damped natural frequency.
    %   dlambda: total derivative of the eigenvalue.
    %   dLambda: total derivative of the eigenvalues.
    %   dPhi: total derivative of the mode shapes.
    % Methods:
    %   MechSystem: constructor.
    %   compute_nonlinear_force: compute the nonlinear force vector.
    %   build_duffing_system: build the Duffing oscillator system.
    %   build_uniform_spring_mass_chain: build the uniform spring-mass
    %       chain system.
    %   build_spring_mass_chain: build the spring-mass chain system.
    %   build_yafec_model: build the yafec model.
    %   build_topopt_system: build the system for topology optimization.
    %   update_topopt_system: update the system for topology optimization.
    %   modal_analysis: perform the linear eigenvalue analysis.
    %   modal_assurance_criterion: compute the modal assurance criterion.
    %   modal_analysis_sensitivity: compute the sensitivity of the
    %       eigenvalues and mode shapes.
    %   set_damping_rayleigh: set the Rayleigh damping coefficients.
    %   set_damping_qfactor: set the Q-factor for the damping.
    %   adjoint_sensitivity_omega: compute the adjoint sensitivity of the
    %       natural frequency.
    %   adjoint_sensitivity_omega: compute the adjoint sensitivity of the
    %       natural frequency for a topology optimization problem.

    properties
        n % number of DOFs
        yafecAssembly % yafec assembly
        M, C, K, T2, T3, D2, D3 % system matrices
        damping_type % type of damping
        rayleigh % Rayleigh damping coefficients
        Qfactor % Q-factor value
        omega, phi, csi, delta, omegaDamped, lambda, tau, psi % 2nd order system eigenvalues and mode shapes
        Lambda, Phi, Psi % master subspace
        pcsi_omega, pdelta_omega, plambda_omega, pLambda_omega % partial derivatives wrt omega
        ndv % number of design variables
        pM, pC, pK, pT2, pT3, pD2, pD3 % partial derivatives of the system matrices
        pC_omega % partial derivative of the damping matrix wrt omega
        dC % total derivative of the damping matrix
        Me, Ce, Ke, T2e, T3e, D2e, D3e % element matrices and tensors
        domega, dphi, dcsi, ddelta, domegaDamped, dlambda, dLambda, dPhi % total derivatives
    end

    methods
        function obj = MechSystem()
            % The constructor.

            % Set default damping type
            obj.damping_type = 0; % 0: no damping, 1: Rayleigh damping, 2: Q-factor damping
        end
    end
end

function [gamma, dgamma] = gamma_sensitivity(M, K, K2, K3, ...
    dMdp, dK1dp, dK2dp, dK3dp, omega, phi)
    % Compute the 3rd order backbone coefficient (gamma) and its adjoint 
    % sensitivity (dgamma) for a parametric optimization problem. The multi-
    % index formulation is employed to compute the adjoint sensitivity.
    % The backbone curve can be written as:
    %   Omega = omega_0 + gamma * rho^2 + O(3)
    % where omega_0 is the linear frequency and rho is the reduced amplitude.
    % Inputs:
    %   M: mass matrix
    %   K: stiffness matrix
    %   K2: quadratic stiffness tensor
    %   K3: cubic stiffness tensor
    %   dMdp: sensitivity of the mass matrix
    %   dK1dp: sensitivity of the stiffness matrix
    %   dK2dp: sensitivity of the quadratic stiffness tensor
    %   dK3dp: sensitivity of the cubic stiffness tensor
    %   omega: eigenfrequency
    %   phi: mode shape
    % Outputs:
    %   gamma: 3rd order backbone coefficient
    %   dgamma: adjoint sensitivity of gamma

    % COMPUTE GAMMA
    % Order 2: index 20
    Lambda20 = 2i * omega;
    f20 = ttv(K2, {phi, phi}, [3, 2]).data;
    L20 = K + Lambda20^2 * M;
    w20 = lsqminnorm(L20, -f20);
    
    % Order 2: index 11
    f11 = 2 * f20;
    L11 = K;
    w11 = lsqminnorm(L11, -f11);
    
    % Compute gamma
    f21 = ttv(K2, {w20 + w11, phi}, [3, 2]).data + ...
          ttv(K2, {w20 + w11, phi}, [2, 3]).data + ...
          3 * ttv(K3, {phi, phi, phi}, [4, 3, 2]).data;
    gamma = 1/(2*omega) * phi.' * f21;

    % COMPUTE ADJOINT SENSITIVITY
    % Partial derivative wrt w20
    f21_w20 = ttv(K2, {phi}, 2) + ttv(K2, {phi}, 3);
    gamma_w20 = 1/(2*omega) * ttv(f21_w20, {phi}, 1).data;
    
    % Partial derivative wrt w11
    gamma_w11 = gamma_w20;
    
    % Partial derivative wrt omega
    gamma_omega = -gamma / omega;
    
    % Partial derivative wrt phi ( = phi = phi)
    f20_phi = f21_w20;
    w20_w11 = w20 + w11;
    f21_phi = ttv(K2, {w20_w11}, 2) + ttv(K2, {w20_w11}, 3) + ...
              3 * ttv(K3, {phi, phi}, [3, 2]) + ...
              3 * ttv(K3, {phi, phi}, [4, 2]) + ...
              3 * ttv(K3, {phi, phi}, [4, 3]);
    gamma_phi = 1/(2*omega) * (f21 + ttv(f21_phi, {phi}, 1).data);

    % Adjoint of w20
    lambda20 = lsqminnorm(L20.', -gamma_w20);
    
    % Adjoint of w11
    lambda11 = lsqminnorm(L11.', -gamma_w11);
    
    % Adjoint of phi and omega
    b_omega = gamma_omega - 8 * omega * lambda20.' * M * w20;
    b_phi = gamma_phi + ttv(f20_phi, lambda20 + 2 * lambda11, 1).data;
    b = [-b_phi; b_omega / omega];
    A = [(K - omega^2 * M), 2 * M * phi;
         2 * phi.' * M, 0];
    x = lsqminnorm(A, b); % A\x triggers warnings on accuracy
    lambda0 = x(1:end-1);
    lambda1 = x(end);    

    % Sensitivities
    np = length(dMdp);      % Number of parameters
    dgamma = zeros(np, 1);  % Initialize sensitivity
    for ii = 1 : np         % Loop over parameters

        dL20 = dK1dp{ii} + Lambda20^2 * dMdp{ii};
        dL11 = dK1dp{ii};

        dT2dT2 = dK2dp{ii} + permute(dK2dp{ii}, [1, 3, 2]);

        % f partial derivatives
        df20 = double( ttv(dK2dp{ii}, {phi, phi}, [3, 2]) );
        df11 = 2 * df20;
        df21 = double( ttv(dT2dT2, {w20 + w11, phi}, [3, 2]) + ...
               3 * ttv(dK3dp{ii}, {phi, phi, phi}, [4, 3, 2]) );
        
        % Compute sensitivity
        dgamma(ii, :) = 1/(2*omega) * phi.' * df21 + ...
                        lambda0.' * (dK1dp{ii} - omega^2 * dMdp{ii}) * phi + ...
                        lambda1 * phi.' * dMdp{ii} * phi + ...
                        lambda20.' * (dL20 * w20 + df20) + ...
                        lambda11.' * (dL11 * w11 + df11);
    end
end

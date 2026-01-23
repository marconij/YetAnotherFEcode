function gamma = gamma_value(M, K, T2, T3, omega, phi)
    % Compute the 3rd order backbone coefficient (gamma).
    % The backbone curve can be written as:
    %   Omega = omega_0 + gamma * rho^2 + O(3)
    % where omega_0 is the linear frequency and rho is the reduced amplitude.
    %
    % Note: this function is called only if the SpecifyObjectiveGradient
    % option is set to false in the optimization problem, otherwise the
    % gamma_sensitivity function is called instead.
    %
    % Inputs:
    %   myAssembly: instance of the Assembly class
    %   M: mass matrix
    %   K: stiffness matrix
    %   T2: 2nd order tensor
    %   T3: 3rd order tensor
    %   omega: frequency
    %   phi: mode shape
    % Outputs:
    %   gamma: 3rd order backbone coefficient

    % Order 2: index 20
    Lambda20 = 2i * omega;
    f20 = ttv(T2, {phi, phi}, [3, 2]).data;
    L20 = K + Lambda20^2 * M;
    w20 = lsqminnorm(L20, -f20);
    
    % Order 2: index 11
    f11 = 2 * f20;
    L11 = K;
    w11 = lsqminnorm(L11, -f11);
    
    % Compute gamma
    f21 = ttv(T2, {w20 + w11, phi}, [3, 2]).data + ...
          ttv(T2, {w20 + w11, phi}, [2, 3]).data + ...
          3 * ttv(T3, {phi, phi, phi}, [4, 3, 2]).data;
    gamma = 1/(2*omega) * phi.' * f21;
end

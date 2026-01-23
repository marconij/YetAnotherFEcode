function set_damping_nonlinear(obj, beta2, beta3)
    % Set the nonlinear damping coefficients for the mechanical system.
    % Inputs:
    %   beta2: quadratic damping coefficient.
    %   beta3: cubic damping coefficient.

    % Set the quadratic and cubic damping tensors
    obj.D2 = beta2 * obj.T2;
    obj.D3 = beta3 * obj.T3;
end

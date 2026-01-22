function px_phi = partial_derivative_x_phi(obj, rho, xk, x, outDof)
    % Compute the partial derivative of the physical amplitude x with
    % respect to the mode shape phi.
    % Inputs:
    %   rho: the reduced amplitude.
    %   xk: the physical amplitudes at all the time steps.
    %   x: the rms of the physical amplitude.
    %   outDof: the index of the output dof.
    % Outputs:
    %   px_phi: the partial derivative of x with respect to phi.

    % Initialize partial derivative
    px_phi = zeros(obj.ssm.sys.n, 1);

    % Reduced coordinates at all the time steps
    p1 = rho * exp( 1i * obj.ssm.theta);
    p2 = rho * exp(-1i * obj.ssm.theta);

    % Compute the partial derivative by summing the contributions from all the time steps
    px_phi(outDof, 1) = sum(xk .* real(p1 + p2)) / (obj.ssm.nTimeSteps * x);
end

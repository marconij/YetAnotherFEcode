function rho = compute_rho_from_x_leading_order(obj, xTarget, outDof)
    % Compute the reduced amplitude rho from the RMS of the physical
    % amplitude x of the output dof at the leading order.
    % At the leading order, the reduced amplitude can be computed
    % analytically.
    % Inputs:
    %   xTarget: RMS of the target physical amplitude of the output dof.
    %   outDof: the index of the output dof.
    % Outputs:
    %   rho: reduced amplitude.

    % Reduced coordinates
    p1Tilde = exp(1i * obj.theta);
    p2Tilde = exp(-1i * obj.theta);

    % Compute the RMS of the physical amplitude at the leading order with rho = 1
    xiTilde = obj.w(1).coeffs(outDof, 1) * p1Tilde + obj.w(1).coeffs(outDof, 2) * p2Tilde;
    xTilde = sqrt(sum(real(xiTilde).^2, 2) / (obj.nTheta - 1)); % SSMTool style

    % Compute the reduced amplitude
    rho = xTarget / xTilde;
end

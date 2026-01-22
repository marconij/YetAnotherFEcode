function f = compute_nonlinear_force(obj, x, dx)
    % Compute the nonlinear force vector.
    % Inputs:
    %   x: the displacement vector.
    %   dx: the velocity vector.
    % Outputs:
    %   f: the nonlinear force vector.

    % Initialize the nonlinear force vector
    f = tensor(zeros(obj.n, 1), obj.n);

    % Add quadratic displacement term
    if nnz(obj.T2) > 0
        f = f + ttv(obj.T2, {x, x}, -1);
    end

    % Add cubic displacement term
    if nnz(obj.T3) > 0
        f = f + ttv(obj.T3, {x, x, x}, -1);
    end

    % Add quadratic velocity term
    if nnz(obj.D2) > 0
        f = f + ttv(obj.D2, {dx, dx}, -1);
    end

    % Add cubic velocity term
    if nnz(obj.D3) > 0
        f = f + ttv(obj.D3, {dx, dx, dx}, -1);
    end

    % Extract the nonlinear force vector
    f = f.data;
end

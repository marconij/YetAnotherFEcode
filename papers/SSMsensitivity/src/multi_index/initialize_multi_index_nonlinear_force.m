function MIs = initialize_multi_index_nonlinear_force(MIs, maxOrder, varargin)
    % For each multi-index, initialize the multi-indices for the nonlinear
    % force contributions.
    % Inputs:
    %   MIs: the multi-index struct.
    %   maxOrder: the maximum expansion order.
    %   minOrder: the minimum expansion order (optional, default is 1).
    % Outputs:
    %   MIs: the multi-index struct with the nonlinear force multi-indices.

    % Parse the input arguments
    p = inputParser;
    addOptional(p, 'minOrder', 1);
    parse(p, varargin{:});
    minOrder = p.Results.minOrder;

    % Loop over the expansion orders
    for mOrder = minOrder:maxOrder
        % Number of multi-indices at the current order
        n = MIs(mOrder).n;

        % Initialize the struct of multi-indices
        MIs(mOrder).coeffsNlForce2 = repmat(struct('u', [], 'k', [], 'n', []), 1, n);
        MIs(mOrder).coeffsNlForce3 = repmat(struct('u', [], 'k', [], 'l', [], 'n', []), 1, n);

        % Loop over the multi-indices at the current order
        for mIdx = 1:n
            % Current multi-index
            m = MIs(mOrder).coeffs(:, mIdx);

            % Add all the combinations multi-indices u of order < mOrder
            % such that their sum is equal to m.
            % For the quadratic force part, the number of
            % multi-indices in the sum are 2, while they are 3 for
            % the cubic force part.
            [u, k] = generate_multi_index_quadratic_force(m);
            MIs(mOrder).coeffsNlForce2(mIdx).u = u;
            MIs(mOrder).coeffsNlForce2(mIdx).k = k;
            MIs(mOrder).coeffsNlForce2(mIdx).n = size(u, 2);

            [u, k, l] = generate_multi_index_cubic_force(m);
            MIs(mOrder).coeffsNlForce3(mIdx).u = u;
            MIs(mOrder).coeffsNlForce3(mIdx).k = k;
            MIs(mOrder).coeffsNlForce3(mIdx).l = l;
            MIs(mOrder).coeffsNlForce3(mIdx).n = size(u, 2);
        end
    end
end

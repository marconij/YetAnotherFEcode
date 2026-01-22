function MIs = initialize_multi_index_V(MIs, maxOrder, varargin)
    % For each multi-index, initialize the multi-indices for the assembly
    % of the V vector.
    % Inputs:
    %   MIs: the multi-index struct.
    %   maxOrder: the maximum expansion order.
    %   minOrder: the minimum expansion order (optional, default is 1).
    % Outputs:
    %   MIs: the multi-index struct with the multi-indices for the V vector.

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
        MIs(mOrder).coeffsV = repmat(struct('j', [], 'u', [], 'k', [], 'n', []), 1, n);

        % Loop over the multi-indices at the current order
        for mIdx = 1:n
            % Current multi-index
            m = MIs(mOrder).coeffs(:, mIdx);

            % Add all the combinations multi-indices u and k such that:
            %   - their order is between 2 and mOrder-1
            %   - their sum is equal to m + e1 or m + e2, where e1 = [1; 0] and e2 = [0; 1]
            %   - the element indicated by ei of the reduced dynamics coefficients R associated to the multi-index k is non-zero
            %     this happens for e1 if k(1) - k(2) = 1 and for e2 if k(1) - k(2) = -1
            [j, u, k] = generate_multi_index_V(m);
            MIs(mOrder).coeffsV(mIdx).j = j;
            MIs(mOrder).coeffsV(mIdx).u = u;
            MIs(mOrder).coeffsV(mIdx).k = k;
            MIs(mOrder).coeffsV(mIdx).n = size(u, 2);
        end
    end
end

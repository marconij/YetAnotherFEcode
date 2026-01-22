function MIs = initialize_multi_index(maxOrder, varargin)
    % Initialize the multi-index object up to the maximum given order.
    % If the object is already initialized, it will add the higher orders
    % if necessary.
    % Inputs:
    %   maxOrder: the maximum expansion order.
    %   MIs: the multi-index struct (optional, default is []).
    % Outputs:
    %   MIs: the multi-index struct with the multi-indices.

    % Parse the input arguments
    p = inputParser;
    addOptional(p, 'MIs', []);
    parse(p, varargin{:});
    MIs = p.Results.MIs;

    % Check if the multi-index object is already initialized
    if isempty(MIs)
        % Initialize multi-indices
        MIs = repmat(struct('coeffs', [], 'n', [], 'coeffsNlForce2', [], 'coeffsNlForce3', [], 'coeffsV', []), 1, maxOrder);
        for mOrder = 1:maxOrder
            MIs(mOrder).coeffs = generate_multi_index(mOrder);
            MIs(mOrder).n = mOrder + 1;
        end

        % Initialize the struct of multi-indices for the nonlinear force contributions
        MIs = initialize_multi_index_nonlinear_force(MIs, maxOrder);

        % Initialize the struct of multi-indices for the assembly of the V vector
        MIs = initialize_multi_index_V(MIs, maxOrder);
        
    elseif length(MIs) < maxOrder
        % Current order
        currentOrder = length(MIs);

        % Add the higher orders
        for mOrder = currentOrder + 1:maxOrder
            MIs(mOrder).coeffs = generate_multi_index(mOrder);
            MIs(mOrder).n = mOrder + 1;
        end

        % Initialize the struct of multi-indices for the nonlinear force contributions
        MIs = initialize_multi_index_nonlinear_force(MIs, maxOrder, 'minOrder', currentOrder + 1);

        % Initialize the struct of multi-indices for the assembly of the V vector
        MIs = initialize_multi_index_V(MIs, maxOrder, 'minOrder', currentOrder + 1);
    
    else
        % The multi-index object is already initialized
        % disp('No need to initialize the multi-index object.');
    end
end

function m = from_position_idx_to_MI(idx, order)
    % Convert the position index to the multi-index.
    % The position index is the index of the multi-index in the multi-index
    % array.
    % Inputs:
    %   idx: position index or array of position indices (1 x n).
    %   order: order of the multi-indices.
    % Outputs:
    %   m: multi-index or array of multi-indices (2 x n).
    val = idx(:).' - 1;
    m = [order - val; val];
end

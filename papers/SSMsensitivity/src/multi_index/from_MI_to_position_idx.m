function idx = from_MI_to_position_idx(m)
    % Convert the multi-index to the position index.
    % The position index is the index of the multi-index in the multi-index
    % array.
    % Inputs:
    %   m: multi-index or array of multi-indices (2 x n).
    % Outputs:
    %   idx: position index or array of position indices (1 x n).
    idx = m(2, :) + 1;
end

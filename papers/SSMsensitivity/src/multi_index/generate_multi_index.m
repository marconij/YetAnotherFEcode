function m = generate_multi_index(order)
    % Generate the multi-indices at the given order.
    % Each column represents a 2D multi-index.
    m = [order:-1:0; 0:1:order];
end
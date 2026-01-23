function [uSet, kSet] = generate_multi_index_quadratic_force(m)
    % Compute the set of all the 2 multi-indices for the composition of the
    % quadratic force term at the given multi-index m.
    % The two multi-indices are such that their sum is equal to m.
    % Inputs:
    %   m: multi-index.
    % Outputs:
    %   uSet: 1st set of multi-indices (2 x n).
    %   kSet: 2nd set of multi-indices (2 x n).

    % % Print for debugging
    % fprintf('MI m:\n')
    % disp(m);
    % fprintf('MIs u, k, u+k:\n')
    % n = 0;
    
    % Find current order
    mOrder = sum(m);

    % Initialize sets
    uSet = [];
    kSet = [];

    % Loop over the 1st MI u order
    for uOrder = 1:mOrder-1
        % Loop over the MIs u
        for uIdx = 1:uOrder+1
            % Current MI u
            u = [uOrder - (uIdx - 1); uIdx - 1];

            % Compute the 2nd MI k
            k = m - u;

            % Check and store
            if all(k >= 0) && sum(k) > 0 && sum(k) < mOrder
                % % Print for debugging
                % disp([u, k, u+k])
                % n = n + 1;

                % Store
                uSet(:, end+1) = u;
                kSet(:, end+1) = k;
            end
        end
    end

    % % Print for debugging
    % fprintf('There are %d MIs.\n\n', n);
    % fprintf('--------------------------------------------------\n\n')
end

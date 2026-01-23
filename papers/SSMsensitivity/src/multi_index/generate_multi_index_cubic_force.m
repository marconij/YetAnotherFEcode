function [uSet, kSet, lSet] = generate_multi_index_cubic_force(m)
    % Compute the set of 3 multi-indices for the composition of the cubic
    % force term at the given multi-index m.
    % The three multi-indices are such that their sum is equal to m.
    % Inputs:
    %   m: multi-index.
    % Outputs:
    %   uSet: 1st set of multi-indices (2 x n).
    %   kSet: 2nd set of multi-indices (2 x n).
    %   lSet: 3rd set of multi-indices (2 x n).

    % % Print for debugging
    % fprintf('MI m:\n')
    % disp(m);
    % fprintf('MIs u, k, l, u+k+l:\n')
    % n = 0;
    
    % Find current order
    mOrder = sum(m);

    % Initialize sets
    uSet = [];
    kSet = [];
    lSet = [];

    % Loop over the 1st MI order
    for uOrder = 1:mOrder-2
        % Loop over the MIs u
        for uIdx = 1:uOrder+1
            % Current MI u
            u = [uOrder - (uIdx - 1); uIdx - 1];

            % Loop over the 2nd MI order
            for kOrder = 1:mOrder-uOrder-1
                % Loop over the MIs k
                for kIdx = 1:kOrder+1
                    % Current MI k
                    k = [kOrder - (kIdx - 1); kIdx - 1];

                    % Compute 3rd MI l
                    l = m - k - u;

                    % Check and store
                    if all(l >= 0) && sum(l) > 0 && sum(l) < mOrder
                        % % Print for debugging
                        % disp([u, k, l, u+k+l])
                        % n = n + 1;

                        % Store
                        uSet(:, end+1) = u;
                        kSet(:, end+1) = k;
                        lSet(:, end+1) = l;
                    end
                end
            end
        end
    end

    % % Print for debugging
    % fprintf('There are %d MIs.\n\n', n);
    % fprintf('--------------------------------------------------\n\n')
end

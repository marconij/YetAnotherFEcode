function [jSet, uSet, kSet] = generate_multi_index_V(m)
    % Compute the set of multi-indices for the assembly of the V vector at
    % the given multi-index m.
    % Inputs:
    %   m: multi-index.
    % Outputs:
    %   jSet: subscripts of the set of unitary multi-indices e (1 x n).
    %   uSet: 1st set of multi-indices (2 x n).
    %   kSet: 2nd set of multi-indices (2 x n).

    % % Print for debugging
    % fprintf('MI m:\n')
    % disp(m);
    % fprintf('MIs e, u, k, u+k-e:\n')
    % n = 0;
    
    % Find current order
    mOrder = sum(m);

    % Initialize sets
    jSet = [];
    % eSet = [];
    uSet = [];
    kSet = [];

    % Loop over the 2nd MI order (odd orders only, starting from 3)
    for kOrder = 3:2:mOrder-1
        % Loop over unitary MI e
        for j = 1:2
            % Compute unitary MI e
            % e = zeros(2, 1);
            % e(j) = 1;

            % Compute MI k depending on e:
            %   - if j = 1, then e = [1; 0] and k = [(kOrder + 1) / 2; (kOrder - 1) / 2], such that k(1) - k(2) = 1
            %   - if j = 2, then e = [0; 1] and k = [(kOrder - 1) / 2; (kOrder + 1) / 2], such that k(1) - k(2) = -1
            k = (kOrder - 1) / 2 * ones(2, 1);
            k(j) = k(j) + 1;

            % Compute the MI u = m + e - k
            % u = m + e - k;
            u = m - k;
            u(j) = u(j) + 1;

            % Check
            uOrder = sum(u);
            if all(u >= 0) && uOrder > 1 && uOrder < mOrder
                % % Print for debugging
                % disp([e, u, k, u+k-e])
                % n = n+1;

                % Store
                jSet(1, end+1) = j;
                % eSet(:, end+1) = e;
                uSet(:, end+1) = u;
                kSet(:, end+1) = k;
            end
        end
    end

    % % Print for debugging
    % fprintf('There are %d MIs.\n\n', n);
    % fprintf('--------------------------------------------------\n\n')
end

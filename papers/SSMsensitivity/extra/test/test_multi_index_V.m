clear; clc; close all;

% Multi-index
m = randi(5, 2, 1) - 1;
% m = [1; 0];
mOrder = sum(m);
fprintf('MI m:\n')
disp(m);

% Initialize
fprintf('MIs e, u, k, u+k-e:\n')
n = 0;

% Loop over the 2nd MI order (odd only, starting from 3)
for kOrder = 3:2:mOrder-1
    % Loop over unitary MI e
    for j = 1:2
        % Current MI e
        e = zeros(2, 1);
        e(j) = 1;

        % Compute MI k depending on e
        k  = [(kOrder + (3 - 2*j)) / 2; (kOrder - (3 - 2*j)) / 2];

        % Compute the MI u
        u = m + e - k;

        % Check
        uOrder = sum(u);
        if all(u >= 0) && uOrder > 1 && uOrder < mOrder
            disp([e, u, k, u+k-e])
            n = n+1;
        end
    end
end
fprintf('There are %d MIs.\n\n', n);
fprintf('--------------------------------------------------\n\n')

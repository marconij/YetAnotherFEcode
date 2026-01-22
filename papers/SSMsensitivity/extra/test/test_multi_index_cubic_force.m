clear; clc; close all;

% Multi-index
% m = randi(5, 2, 1) - 1;
m = [3; 2];
mOrder = sum(m);
fprintf('MI m:\n')
disp(m);

% Initialize
fprintf('MIs u, k, u+k:\n')
n = 0;

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

                % Check
                if all(l >= 0) && sum(l) > 0 && sum(l) < mOrder
                    disp([u, k, l, u+k+l])
                    n = n + 1;
                end
            end
        end
    end
end
fprintf('There are %d MIs.\n\n', n);
fprintf('--------------------------------------------------\n\n')


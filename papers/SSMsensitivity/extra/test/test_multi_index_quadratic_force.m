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
for uOrder = 1:mOrder-1
    % Loop over the MIs u
    for uIdx = 1:uOrder+1
        % Current MI u
        u = [uOrder - (uIdx - 1); uIdx - 1];

        % Compute the 2nd MI k
        k = m - u;

        % Check
        if all(k >= 0) && sum(k) > 0 && sum(k) < mOrder
            disp([u, k, u+k])
            n = n + 1;
        end
    end
end
fprintf('There are %d MIs.\n\n', n);
fprintf('--------------------------------------------------\n\n')

%%

u1Max = min(m(1), mOrder - 1);
u1Min = 0;
for u1 = u1Max:-1:u1Min
    u2Max = min(m(2), mOrder - 1 - u1);
    u2Min = max(0, 1 - u1);
    for u2 = u2Max:-1:u2Min
        u = [u1; u2];
        k = m - u;
        disp([u, k])
    end
end

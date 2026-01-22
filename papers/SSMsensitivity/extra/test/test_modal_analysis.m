clear; clc; close all;

%% Parameters
m = 1;
k = 10;
c = 0.1;

% 2nd order matrices
M = m * eye(2);
K = k * [2, -1; -1, 2];
C = c * [2, -1; -1, 2];

% 1st order matrices (w/o damping)
A0 = blkdiag(-K, M);
B0 = [zeros(size(M)), M; M, zeros(size(M))];

% 1st order matrices (w/ damping)
A = blkdiag(-K, M);
B = [C, M; M, zeros(size(M))];

%% Undamped system
% 2nd order system
[PHI0, OMEGA02, PSI0] = eig(K, M);

% First mode
omega = sqrt(OMEGA02(1, 1)); % eigenfrequency [rad/s]
phi = PHI0(:, 1); % mode shape
phi = phi ./ sqrt(phi.'*M*phi); % mass normalization

% Build 1st order solution
lambda0 = 1i*omega;
LAMBDA0 = [lambda0; conj(lambda0)];
v0 = [phi; lambda0*phi];
V0 = [v0, conj(v0)];
tau0 = 1 / (2*lambda0); % mass normalization
psi0 = conj(tau0) * phi;
u0 = [psi0; conj(lambda0)*psi0];
U0 = [u0, conj(u0)];

% 1st order system
disp('1st order undamped eigenvalue problem (right eigenvectors):')
disp(sum(abs(A0 * V0 - B0 * V0 * diag(LAMBDA0)), 'all'))

disp('1st order undamped eigenvalue problem (left eigenvectors):')
disp(sum(abs(U0' * A0 - diag(LAMBDA0) * U0' * B0), 'all'))

disp('1st order undamped eigenvalue problem (mass normalization):')
disp(sum(abs(U0' * B0 * V0 - eye(2)), 'all'))
% [V01, D01, U01] = eig(A0, B0);
% LAMBDA01 = diag(D01);
% LAMBDA01 = LAMBDA01(3:4);
% V01 = V01(:, 3:4);
% U01 = U01(:, 3:4);

%% Damped system (assume mass-normalized mode shapes)

% Damping parameters
csi = (phi.' * C * phi) / (2*omega);
tau = -1i / (2*omega * sqrt(1 - csi^2));
psi = conj(tau) * phi; % mass normalization

% Build 1st order solution
lambda = -csi*omega + 1i*omega * sqrt(1 - csi^2);
LAMBDA = [lambda; conj(lambda)];
v = [phi; lambda*phi];
V = [v, conj(v)];
u = [psi; conj(lambda)*psi];
U = [u, conj(u)];

% Quadratic eigenvalue problem
[X, e, s] = polyeig(K, C, M);
disp('Quadratic eigenvalue problem:')
disp(sum(abs(M*[phi, phi]*diag(LAMBDA).^2 + C*[phi, phi]*diag(LAMBDA) + K*[phi, phi]), 'all'))

% 1st order system
disp('1st order damped eigenvalue problem (right eigenvectors):')
disp(sum(abs(A * V - B * V * diag(LAMBDA)), 'all'))

disp('1st order damped eigenvalue problem (left eigenvectors):')
disp(sum(abs(U' * A - diag(LAMBDA) * U' * B), 'all'))

disp('1st order damped eigenvalue problem (mass normalization):')
disp(sum(abs(U' * B * V - eye(2)), 'all'))
% [V1, D1, U1] = eig(A, B);
% LAMBDA1 = diag(D1);
% LAMBDA1 = LAMBDA1(3:4);
% V1 = V1(:, 3:4);
% U1 = U1(:, 3:4);









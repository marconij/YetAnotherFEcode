function FRC = FRC_using_SSMtool(obj,order,imod,Qfactor,epsilon,omegaRange_norm,SSMTool_version)
% Compute the Frequency Response Curve using SSMtool and plot
% Inputs:
%   obj: instance of the class MultiPointConstraints
%   order: order of the manifold to be computed
%   imod: mode index
%   Qfactor: Q factor vector
%   epsilon: amplitude of the external force
%   omegaRange_norm: frequency range normalized by the imod-th natural frequency
%   SSMTool_version (optional): either '2.5' or '2.6' (default)
% Outputs:
%   FRC: frequency response curve

if nargin < 7
    SSMTool_version = '2.6';
end

t0 = tic;

% add SSMtool path, but remove the YetAnotherFEcode included therein.
run    ..\..\..\SSMTool\install.m
rmpath(fullfile('..\..\..\SSMTool\ext\tensor_toolbox'))
rmpath(genpath('..\..\..\SSMTool\ext\YetAnotherFEcode\src'))
rmpath(genpath('..\..\..\SSMTool\ext\YetAnotherFEcode\examples'))
fprintf([' SSMTool (v' SSMTool_version ') added to path.\n\n'])

% Mass and stiffness matrices
u0 = zeros(obj.nDOFs.all_structural, 1);
[K,~] = obj.tangent_stiffness_and_force(u0);
M = obj.mass_matrix();

% RAYLEIGH: C = alpha*M + beta*K0 *****************************************
n_VMs = 3;
[~,omega2] = eigs(K,M,n_VMs,'SM');
omega = sort(sqrt(diag(omega2)));

[alpha, beta] = damping_rayleigh('MK', Qfactor, omega(imod+[0 1])/2/pi, 1);
C = alpha*M + beta*K;
csi = 1./(2*Qfactor);       % damping ratio
% *************************************************************************

% Nonlinear force (quadratic and cubic stiffness tensors)
[K2,K3] = obj.mpc_tensors;
fnl{1} = K2;
fnl{2} = K3;

outdof = size(K,1)-1; % y-displacement of the drive frame

% Dynamical system setup
if strcmpi(SSMTool_version,'2.6')
    DS = DynamicalSystem(2);
else
    DS = DynamicalSystem();
end
set(DS,'M',M,'C',C,'K',K);
set(DS,'fnl',fnl)
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')

% Linear Modal analysis and SSM setup
DS.linear_spectral_analysis();

% external force
fext = 1*sparse(outdof,ones(size(outdof)),1,obj.nDOFs.mpc,1);
% epsilon = 0.42;
kappas = [-1; 1];
coeffs = [fext fext]/2;
DS.add_forcing(coeffs, kappas, epsilon);
% DS.Fext.data(2).F_n_k.coeffs

% *Choose Master subspace (perform resonance analysis)*
S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
masterModes = 2*imod-[1,0]; 
S.choose_E(masterModes);

% setup options
set(S.Options, 'reltol', 1,'IRtol',0.02,'notation', 'multiindex','contribNonAuto',true,'COMPtype','second')
set(S.FRCOptions, 'nt', 2^7, 'nRho', 200, 'nPar', 200, 'nPsi', 200, 'rhoScale', 2)
set(S.FRCOptions, 'method','continuation ep')
set(S.FRCOptions, 'outdof',outdof)

set(S.contOptions,'MaxRes',100) % max residual for continuation
set(S.contOptions,'h_max',1e3,'h_min',1e-8) % max and min step size for continuation
set(S.contOptions,'h_fac_max',1e3,'h_fac_min',1e-3) % max and min step size factor for continuation
set(S.contOptions,'PtMX',1000) % max number of continuation steps
set(S.FRCOptions,'initialSolver','fsolve') % initial solver for the continuation problem

% FRC
% choose frequency range around the first natural frequency
omega0 = abs(imag(S.E.spectrum(2)));
omegaRange = omega0*omegaRange_norm;
disp('omega range:')
disp(omegaRange)

% extract forced response curve
% FRC = S.extract_FRC('freq',omegaRange,order);
FRC = S.SSM_isol2ep('contep',masterModes,order,1,'freq',omegaRange,outdof);

% Linear FRF
w = [linspace(omegaRange(1), omega0*sqrt(1-csi(imod)^2), 10000) ...
     linspace(omega0*sqrt(1-csi(imod)^2), omegaRange(2), 10000)];
u_out = w*0;
for ii = 1 : length(w)
    u = (-w(ii)^2*M + 1i*w(ii)*C + K)\(epsilon*fext);
    u_out(ii) = abs(u(outdof));
end

plot(w, u_out,'k-','LineWidth',1.5,'DisplayName','linear')
axis tight
xlabel('$\Omega$','Interpreter','latex')
set(gca,"FontSize",14)
grid on

toc(t0)

warning off
rmpath(genpath('..\..\..\..\SSMTool\src'))
rmpath(genpath('..\..\..\..\SSMTool\ext'))
warning on
fprintf([' SSMTool (v' SSMTool_version ') removed from path.\n\n'])
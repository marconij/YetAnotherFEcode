function bb = verify_using_SSMtool(obj,order,imod,SSMTool_version)
% Verify the backbone curve using SSMtool
% Inputs:
%   obj: instance of the class MultiPointConstraints
%   order: order of the manifold to be computed
%   imod: mode index
%   SSMTool_version (optional): either '2.5' or '2.6' (default)
% Outputs:
%   bb: backbone curve

if nargin < 4
    SSMTool_version = '2.6';
end

% add SSMtool path, but remove the YetAnotherFEcode included therein.
%  Note: the following commands assume that the "SSMTool" folder is at the
%  same level as the YetAnotherFEcode folder.
run    ..\..\..\SSMTool\install.m
rmpath(fullfile('..\..\..\SSMTool\ext\tensor_toolbox'))
rmpath(genpath('..\..\..\SSMTool\ext\YetAnotherFEcode\src'))
rmpath(genpath('..\..\..\SSMTool\ext\YetAnotherFEcode\examples'))
fprintf([' SSMTool (v' SSMTool_version ') added to path.\n\n'])

% Mass and stiffness matrices
u0 = zeros(obj.nDOFs.all_structural, 1);
[K,~] = obj.tangent_stiffness_and_force(u0);
M = obj.mass_matrix();

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
set(DS,'M',M,'C',K*0,'K',K);
set(DS,'fnl',fnl)
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')

% Linear Modal analysis and SSM setup
DS.linear_spectral_analysis();

% *Choose Master subspace (perform resonance analysis)*
S = SSM(DS);

% setup options
set(S.Options, 'reltol', 1, 'IRtol', 0.02, 'contribNonAuto', false,...
    'notation', 'multiindex', 'COMPtype', 'second')
set(S.Options,'solver','backslash')
set(S.FRCOptions, 'nt', 2^7, 'nRho', 200, 'nPar', 200, 'nPsi', 200, 'rhoScale', 2 )
set(S.FRCOptions, 'outdof',outdof)

% choose frequency range around the first natural frequency
omega0 = imag(DS.spectrum.Lambda(imod*2-1));
omegaRange = omega0*[0 inf];

% compute backbone curve
% to achieve (approximately) umax, we rescale rhomax based on the max 
% displacement of the linear eigenvector at the outdof:
umax = 5; % [um] max desired displacement
rho_max = abs(umax / sum(S.System.spectrum.V(outdof,imod*2+[-1 0])));

bb = S.extract_backbone(imod*2+[-1 0],omegaRange,order,rho_max);
drawnow

warning off
rmpath(genpath('..\..\..\..\SSMTool\src'))
rmpath(genpath('..\..\..\..\SSMTool\ext'))
warning on
fprintf([' SSMTool (v' SSMTool_version ') removed from path.\n\n'])
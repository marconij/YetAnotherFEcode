function [Omega, wBB_rms, wBB_infNorm] = analyze_system(obj)
    % Analyze the system computing eigenmodes, backbone curve and Lyapunov
    % Subcenter manifold of the drive mode.
    % Input:
    %   obj: instance of class MultiPointConstraints
    % Output:
    %   Omega: frequency of the backbone curve
    %   wBB_rms: RMS value of the backbone curve (computed at outdof)
    %   wBB_infNorm: infinity norm of the backbone curve (computed at outdof)

% linearized system matrices
u0   = zeros( obj.Mesh.nDOFs, 1);
Kmpc = obj.tangent_stiffness_and_force(u0);
Mmpc = obj.mass_matrix;

% Eigenvalue problem                                               
n_VMs = 7; % first n_VMs modes with lowest frequency calculated 
[Phi1,om2] = eigs(Kmpc, Mmpc, n_VMs, 'SM');
[om, ind] = sort(sqrt(diag(om2)));
f0 = om/2/pi;
Phi = Phi1(:,ind);

% Mass normalize vibration modes
for ii = 1 : size(Phi,2)
    Phi(:,ii) = Phi(:,ii)/sqrt(Phi(:,ii)'*Mmpc*Phi(:,ii));
end

% Plot the first two modes (sense and drive)
scale = 100;
for imod = 1 : 2
    figure
    obj.mpc_plot(Phi(:,imod),scale)
    u_mass = Phi(end-5:end-3,imod);
    [~,ind] = max(abs(u_mass));
    if ind==1
        t = 'Sense, ';
    elseif ind==2
        t = 'Drive, ';
        imod_gamma = imod;
    else
        t = '';
    end
    title(sprintf('%sMode %d, f=%.0f Hz', t, imod, f0(imod)))
end

% Manifold and Backbone
[T2,T3] = obj.mpc_tensors; % get stiffness tensors
outdof = size(Kmpc,1)-1;   % y-displacement of the drive frame
doPlot = 1;
rho_max = 2e-4; % note: here eigenmodes are mass-normalized
rho = linspace(0,rho_max,100); % vector of the reduced amplitude of the backbone curve
fprintf(' mode: %d\n',imod_gamma)
[Omega, theta, wManifold, wBB_rms, wBB_infNorm] = manifold_and_backbone ( ...
    Mmpc, Kmpc, T2, T3, om(imod_gamma), Phi(:,imod_gamma), rho, outdof, doPlot);

end
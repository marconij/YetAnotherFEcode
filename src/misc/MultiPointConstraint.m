classdef MultiPointConstraint < Assembly
    % This class implements the methods to construct a system using von 
    % Kàrman beams connected to rigid bodies using a Multi-Point constraint
    % approach. 

    properties
        mpc         % [#slave, #master]
        masterNodes % coordinates of new master nodes
        Tmpc        % mpc transformation matrix
        masterDOFs  % a 3-columns matrix mapping master DOFs between full and mpc vectors: [DOFs in "full" vector (structural free & constrained, slaves and all master DOFs), DOFs in the "mpc" vector (structural free & all master DOFs), is constrained?]
        Mm          % cell containing the master lumped mass matrices
        Km          % cell containing the master lumped stiffness matrices
        PLOT        % structure with:
                    %   .elem:      cell with element nodes for plotting
                    %   .mass_xy:   cell with master mass shapes for plotting
        nNodes      % info about number of nodes
        nDOFs       % info about number of DOFs
    end

    methods
        function self = MultiPointConstraint(Mesh, mpc, masterNodes, Mm, Km, par)
            % Constructor of the class.
            % INPUTS:
            %   Mesh:       mesh object
            %   mpc:        matrix with slave-master nodes pairs
            %   masterNodes:matrix with master nodes coordinates and constraints
            %   Mm:         cell with master lumped mass matrices
            %   Km:         cell with master lumped stiffness matrices
            %   par:        parallel computing flag
            % OUTPUTS:
            %   self:       MultiPointConstraint object

            if nargin < 4
                for ii = 1 : size(masterNodes,1)
                    Mm{ii} = zeros(3);
                    Km{ii} = zeros(3);
                end
                par = false;
            elseif nargin < 5
                for ii = 1 : size(masterNodes,1)
                    Km{ii} = zeros(3);
                end
                par = false;
            elseif nargin < 6
                par = false;
            end
            if sum(masterNodes(:,3:5)~=0 & masterNodes(:,3:5)~=1,'all')>0
                error('master node constrains can be specified as 1 (fixed) or 0 (free)')
            end

            self@Assembly(Mesh,par);
            
            if ~isequal(class(self.Mesh.Elements(1).Object),'BeamElement')
                error('This class works only with BeamElement!')
            end

            self.mpc = mpc;
            self.masterNodes = masterNodes;
            
            self.Mm = Mm;
            self.Km = Km;
            
            % nodes and DOFs information
            self.nNodes.structure = self.Mesh.nNodes;
            self.nNodes.slaves = size(mpc,1);
            self.nNodes.masters = size(unique(mpc(:,2)),1);
            
            self.nDOFs.all_structural = self.Mesh.nDOFs;
            self.nDOFs.all_slaves = self.nNodes.slaves * 3;
            self.nDOFs.all_masters = self.nNodes.masters * 3;
            self.nDOFs.constrained_structural = size(self.Mesh.EBC.constrainedDOFs,1);
            self.nDOFs.constrained_masters = sum(masterNodes(:,3:5),"all");
            self.nDOFs.free_structural = self.nDOFs.all_structural - self.nDOFs.constrained_structural;
            self.nDOFs.free_masters = sum(masterNodes(:,3:5)==0,"all");
            self.nDOFs.mpc = self.nDOFs.free_structural - self.nDOFs.all_slaves + self.nDOFs.free_masters;

            self.Tmpc = construct_Tmatrix(self);
        end

        function T = construct_Tmatrix(self)
            % This method constructs the transformation matrix T that maps 
            % the full displacement vector to the MPC displacement vector. 
            % The latter contains the free structural DOFs and all the 
            % master DOFs. For this reason, after transformation, the 
            % constrained master DOFs have to be removed from the MPC 
            % vector/matrices/tensors.

            constrained_dofs = self.Mesh.EBC.constrainedDOFs(:,1);
            constrained_nodes = unique(ceil(constrained_dofs/3));
            
            ncd = self.nDOFs.constrained_structural;

            ns = self.nNodes.slaves;
            nm = self.nNodes.masters;
            nn = self.nNodes.structure;

            T = zeros((nn+nm)*3-ncd, (nn-ns+nm)*3-ncd);
            ind = 1;
            for ii = 1 : nn % loop over all structure nodes

                node_dofs = self.Mesh.get_DOF_from_nodeIDs(ii);

                if ~isempty(find(constrained_nodes==ii, 1))
                    % CASE: CONSTRAINED NODE
                    % if only some, but not all DOFs of this constrained 
                    % node are fixed, include an identity matrix of the 
                    % proper size in the T matrix.
                    kk = 1;
                    unc_node_dofs = [];
                    I = [];
                    % Loop over the DOFs of the constrained node:
                    for jj = 1 : length(node_dofs) 
                        % check if the DOF is fixed
                        if ~sum(node_dofs(jj)==constrained_dofs)==1
                            % if the DOF is unconstrained, append a 1 to
                            % the identity matrix and store the DOF#
                            I(kk,kk) = 1; %#ok<AGROW>
                            unc_node_dofs(kk) = node_dofs(jj); %#ok<AGROW>
                            kk = kk+1;
                        end
                    end
                    % if there are unconstrained DOFs, include them in the
                    % T matrix, otherwise continue.
                    if ~isempty(unc_node_dofs)
                        node_dofs = self.free2constrained_index( unc_node_dofs );
                        T(node_dofs, ind : ind+size(I,1)-1) = I;
                        ind = ind + size(I,1);
                    end
                elseif sum(self.mpc(:,1)==ii)>0
                    % if the current node ii is a slave node, then compute
                    % the R matrix and include it in the transformation
                    % matrix T.
                    coord_slave = self.Mesh.nodes(ii,:);
                    mpc_ind = find(self.mpc(:,1)==ii);
                    master_node = self.mpc(mpc_ind, 2);
                    coord_master = self.masterNodes(master_node, 1:2);
                    delta = coord_slave - coord_master;
                    c = self.masterNodes(master_node, 3:5)==0;
                    R = [c(1) 0 -delta(2)*c(3);
                         0 c(2)  delta(1)*c(3);
                         0 0              c(3)];
                    node_dofs = self.free2constrained_index( node_dofs );
                    T(node_dofs, end-(nm-master_node+1)*3+1:end-(nm-master_node+1)*3+3) = R;
                else
                    % put identity
                    node_dofs = self.free2constrained_index( node_dofs );
                    T(node_dofs, ind : ind+2) = eye(3);
                    ind = ind+3;
                end
            end
            % map master DOFs from full to mpc vector and put an identity
            % matrix in the T matrix at the proper position:
            cont = 0;
            self.masterDOFs = [];
            for ii = 1 : nm
                % "full" x vector contains free dofs, slave dofs, and master dofs
                master_full_dofs = nn*3-ncd+1+cont : nn*3-ncd+3+cont;

                % "mpc" x vector contains free dofs and master dofs
                master_mpc_dofs = ((nn-ns+nm)*3-ncd)-(nm-ii+1)*3+1 : ((nn-ns+nm)*3-ncd)-(nm-ii+1)*3+3;

                T(master_full_dofs, master_mpc_dofs) = eye(3);
                self.masterDOFs = [self.masterDOFs; [master_full_dofs' master_mpc_dofs' self.masterNodes(ii,3:5)']];

                cont = cont+3;
            end
            T = sparse(T);
        end

        function [K_mpc,F_mpc] = tangent_stiffness_and_force(self,u,Km)
            % This method computes the tangent stiffness matrix and the
            % internal force vector for the system with Multi-Point
            % constraints. The method returns the tangent stiffness matrix
            % and the internal force vector with the constrained master 
            % DOFs removed.
            % INPUTS:
            %   u:  displacement vector
            %   Km: cell with master lumped stiffness matrices (optional)
            % OUTPUTS:
            %   K_mpc: tangent stiffness matrix
            %   F_mpc: internal force vector
            
            if nargin<3
                Km = self.Km;
            end
            Ktemp=[];
            for ii = 1 : length(Km)
                Ktemp = blkdiag(Ktemp, Km{ii});
            end
            Km = Ktemp;
            if size(Km,1) ~= self.nDOFs.all_masters
                error('Wrong size for the master nodes stiffness matrix')
            end

            [K,F] = tangent_stiffness_and_force@Assembly(self,u);
            Kc = self.constrain_matrix(K);
            Fc = self.constrain_vector(F);

            K_mpc = self.Tmpc' * blkdiag(Kc, Km) * self.Tmpc;
            F_mpc = self.Tmpc' * [Fc; zeros(self.nDOFs.all_masters,1)];

            c = self.masterDOFs(self.masterDOFs(:,3)==1,2); % constrained master dofs
            K_mpc(c,:) = [];
            K_mpc(:,c) = [];
            K_mpc = (K_mpc+K_mpc')/2;
            F_mpc(c) = [];
        end

        function M_mpc = mass_matrix(self,Mm)
            % This method computes the mass matrix for the system with
            % Multi-Point constraints. The method returns the mass matrix
            % with the constrained master DOFs removed.
            % INPUTS:
            %   Mm: cell with master lumped mass matrices (optional)
            % OUTPUTS:
            %   M_mpc: Multi-Point constrained mass matrix

            if nargin<2
                Mm = self.Mm;
            end
            Mtemp=[];
            for ii = 1 : length(Mm)
                Mtemp = blkdiag(Mtemp, Mm{ii});
            end
            Mm = Mtemp;
            if size(Mm,1) ~= self.nDOFs.all_masters
                error('Wrong size for the master nodes mass matrix')
            end

            M = mass_matrix@Assembly(self);
            Mc = self.constrain_matrix(M);

            M_mpc = self.Tmpc' * blkdiag(Mc, Mm) * self.Tmpc;

            c = self.masterDOFs(self.masterDOFs(:,3)==1,2); % constrained master dofs
            M_mpc(c,:) = [];
            M_mpc(:,c) = [];
            M_mpc = (M_mpc+M_mpc')/2;
        end

        function [K2_mpc, K3_mpc] = mpc_tensors(self)
            % This method computes the second and third order tensors for
            % the system with Multi-Point constraints. The method returns
            % the third and fourth order tensors with the constrained master
            % DOFs removed.
            % OUTPUTS:
            %   T2_mpc: quadratic stiffness (third order) tensor
            %   T3_mpc: cubic stiffness (fourth order) tensor

            % tic
            % fprintf(' Building tensors... \t\t')
            K2 = self.tensor('T2',[self.Mesh.nDOFs, self.Mesh.nDOFs, self.Mesh.nDOFs], [2,3]);
            K3 = self.tensor('T3',[self.Mesh.nDOFs, self.Mesh.nDOFs, self.Mesh.nDOFs, self.Mesh.nDOFs], [2,3,4]);

            T = sptensor(self.Tmpc);
            K2c = self.constrain_tensor(K2);
            K3c = self.constrain_tensor(K3);

            nm = self.nDOFs.all_masters;
            d  = self.nDOFs.free_structural + nm; % that is: size(self.Tmpc,1)
            K2_temp = sptensor([],[],[d,d,d]);
            K2_temp(1:end-nm, 1:end-nm, 1:end-nm) = K2c;
            K3_temp = sptensor([],[],[d,d,d,d]);
            K3_temp(1:end-nm, 1:end-nm, 1:end-nm, 1:end-nm) = K3c;
            % toc

            K2_mpc = ttt(ttt(ttt(T,K2_temp,1,1),T,3,1),T,2,1);
            K3_mpc = ttt(ttt(ttt(ttt(T,K3_temp,1,1),T,3,1),T,2,1),T,2,1);

            % APPLY MASTER NODE CONSTRAINTS
            % unfortunately, slices of a sparse tensor cannot be removed.
            % One option is to use a full matrix, partition it, and then
            % convert it back to a sparse tensor. This approach is however
            % memory hungry and very slow. For mildly large system, it is
            % not applicable. We tested other tensor packages, but each of
            % them is lacking some important features.
            % Here we deviced a strategy to circumvent the problem, working
            % with sparse tensors only. To this aim, The Tensor Toolbox for
            % MATLAB v3.6 is necessary (as it includes the "squash"
            % method). Notice that yafec 2.0 (current version) includes
            % the Tensor Toolbox v3.1, so "squash" has to be manually
            % added. For convenience, we pasted the function at the end of
            % this script.
            % https://www.tensortoolbox.org/index.html

            c = self.masterDOFs(self.masterDOFs(:,3)==1,2); % constrained master dofs
            
            if ~isempty(c)
                n = self.nDOFs.free_structural - self.nDOFs.all_slaves + self.nDOFs.all_masters; % that is: size(self.Tmpc, 2)
                offset = 7 + rand(1); % random offset to avoid zero entries
                o = ones(n,1);
                o(c) = 0; % set to zero the constrained master dofs
    
                I2 = sptensor([(1:n)' (1:n)' (1:n)'], o, [n,n,n]);
                % add the offset diagonal tensor I2 in order to preserve the zero 
                % entries in T2_mpc corresponding to free master nodes
                ti = K2_mpc + I2*offset; 
                % now squash the tensor to remove the zero-slices (corresponding
                % to the constrained master DOFs)
                ts = squash(ti);
                m = size(ts,1);
                % now remove the offset diagonal tensor
                K2_mpc = ts - sptensor([(1:m)' (1:m)' (1:m)'], ones(m,1)*offset);
                
                % same procedure for the third order tensor:
                I3 = sptensor([(1:n)' (1:n)' (1:n)' (1:n)'], o, [n,n,n,n]);
                ti = K3_mpc + I3*offset;
                ts = squash(ti);
                m = size(ts,1);
                K3_mpc = ts - sptensor([(1:m)' (1:m)' (1:m)' (1:m)'], ones(m,1)*offset);
            end

        end

        function [u_structure_unconstrained, u_masters] = unconstrain_mpc_vector(self,u)
            % This method removes the constraints from the MPC displacement
            % vector. The method returns the unconstrained structural 
            % displacement vector and the master displacement vector in full
            % form, i.e. with all structural (including constrained and slave
            % DOFs) and master DOFs.
            % INPUTS:
            %   u:  MPC displacement vector (can be a matrix)
            % OUTPUTS:
            %   u_structure_unconstrained: unconstrained structural displacement vector (all DOFs)
            %   u_masters: master displacement vector (all DOFs)

            nm_all_dof = size(self.masterDOFs,1);
            nm_free_dof = sum(self.masterDOFs(:,3)==0);

            u_structure = u(1:end-nm_free_dof,:); % master DOFs are removed
            u_master = u(end-nm_free_dof+1:end,:); % only master DOFs
            u_master2 = zeros(size(self.masterDOFs,1), size(u,2));
            u_master2(self.masterDOFs(:,3)==0, :) = u_master;

            u = [u_structure; u_master2]; % mpc "all-master-DOFs" vector (can be multiplied by Tmpc)

            u_all = self.Tmpc*u; % slave DOFs included

            u_structure_constrained = u_all(1 : end-nm_all_dof,:); % structural DOFs only
            u_masters = u_all(end-nm_all_dof+1 : end,:); % master DOFs only
            
            % unconstrain the structural displacement vector(s):
            u_structure_unconstrained = zeros(self.Mesh.nNodes*3, size(u,2));
            for ii = 1 : size(u,2)
                u_structure_unconstrained(:,ii) = self.unconstrain_vector(u_structure_constrained(:,ii));
            end
        end

        function mpc_plot(self,U0,scale)
            % This method plots the deformed structure with the master nodes
            % and the master mass shapes. If called with no arguments, it
            % plots the undeformed structure with the master nodes and the
            % master mass shapes.
            % INPUTS:
            %   U0:     displacement vector (optional)
            %   scale:  scale factor for the deformed structure (optional)
            % REQUIREMENTS:
            %   The method requires that the PLOT structure is defined in the
            %   MPC object. The PLOT structure is a field of the class that
            %   contains the element nodes and the master mass shapes for
            %   plotting.
            
            if nargin<2
                scale = 0;
            elseif nargin<3
                scale = 10;
            end

            if scale~=0
                [U, Um] = self.unconstrain_mpc_vector(U0);
                a = max(abs([U(1:3:end); U(2:3:end); Um(1:3:end); Um(2:3:end)]));
                U  =  U/a*scale;
                Um = Um/a*scale;
                u = U(1:3:end);
                v = U(2:3:end);
            end
            x = self.Mesh.nodes(:, 1);
            y = self.Mesh.nodes(:, 2);

            elements = self.PLOT.elem;

            % undeformed structure
            for ii = 1 : length(elements)
                e = elements{ii}';
                e = unique(e(:));
                plot(x(e),y(e), 'k.-')
                hold on
            end
            mm = self.PLOT.mass_xy;
            for ii = 1 : length(mm)
                patch(mm{ii}(:,1), mm{ii}(:,2),[0.0745 0.6235 1],'FaceAlpha',0)
            end

            % deformed structure
            if scale~=0
                for ii = 1 : length(elements)
                    e = elements{ii}';
                    e = unique(e(:));
                    m = sqrt(u(e).^2 + v(e).^2);
                    patch([x(e)+u(e); NaN], [y(e)+v(e); NaN], [m; NaN], ...
                        'EdgeColor','interp','Marker','.', ...
                        'MarkerFaceColor','flat','LineWidth',2)
                end
                for ii = 1 : length(mm)

                    th = -Um(ii*3);
                    R = [cos(th) sin(th); -sin(th) cos(th)];
                    xm = mm{ii}(:,1)';
                    ym = mm{ii}(:,2)';
                    xG = mean(xm);
                    yG = mean(ym);
                    Q = R*[xm-xG; ym-yG];
                    rx = Q(1,:);
                    ry = Q(2,:);

                    patch(xG+rx+Um(ii*3-2), yG+ry+Um(ii*3-1), [0.0745    0.6235    1.0000])
                    plot(xG+Um(ii*3-2),     yG+Um(ii*3-1), 'ko','MarkerFaceColor','r')
                end
                colormap jet; colorbar;
            end

            axis equal; grid on;
            xlabel('x')
            ylabel('y')
            
            % in case scale=0, just plot geometry with masters and slaves
            if scale==0
                for ii = 1 : length(mm) % loop on master nodes
                    xm = mm{ii}(:,1);
                    ym = mm{ii}(:,2);
                    xG = mean(unique(xm));
                    yG = mean(unique(ym));
                    
                    slave_IDs = self.mpc(self.mpc(:,2)==ii);
                    xs = x(slave_IDs);
                    ys = y(slave_IDs);
                    for jj = 1 : length(slave_IDs)
                        plot([xG xs(jj)],[yG ys(jj)],'r--') % plot master-slaves connection
                        plot(xs(jj), ys(jj), 'ko','MarkerFaceColor','y') % plot slave node
                    end
                    plot(xG, yG, 'ko','MarkerFaceColor','r') % plot master node
                    text(xG, yG, [' ' num2str(ii)],'HorizontalAlignment','left','VerticalAlignment','top')
                    c = self.masterNodes(ii,3:5)==1;
                    if sum(c)>0
                        for kk = 1:sum(c)
                            plot(xG, yG, 'bo','MarkerSize',6*(1+kk)) % plot constrains on master nodes
                        end
                    end
                end
            end
        end

    end
end

function [Y,S] = squash(X)
%SQUASH Remove empty slices from a sparse tensor.
%
%   Y = SQUASH(X) returns a sparse tensor Y with the same elements as
%   X but with all the empty slices removed.  The indices appearing in the
%   tensor X for each dimension n are remapped to the range [1:M_n] where M_n
%   is the number of unique indices for dimension n.
%
%   [Y,S] = squash(X) also returns a cell-array of length ndims(X) which
%   specifies which indices in X the indices in Y correspond to, i.e.,
%   X.subs(:,n) == S{n}(Y.subs(:,n)) for each n.
%
%   See also SPTENSOR.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

if ~isa(X,'sptensor')
    error('Input must be an sptensor!');
end

d = ndims(X);
subs = zeros(size(X.subs));
sz = zeros(1,d);
if nargout == 2
    S = cell(d,1);
end
for n=1:d
    [s,~,j] = unique(X.subs(:,n));
    subs(:,n) = j;
    sz(n) = length(s);
    if nargout == 2
        S{n} = s;
    end
end

Y = sptensor(subs, X.vals, sz);
end
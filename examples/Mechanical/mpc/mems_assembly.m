function obj = mems_assembly(model,plot_flag)

if nargin==1
    plot_flag=0;
end

thickness = model.thickness;
width     = model.width;
Nel       = model.nElements;
joints    = model.joints;

% create material
myMaterial  = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',model.E,'DENSITY',model.rho,'POISSONS_RATIO',model.nu);

nb = 0;
for ii = 1 : length(joints)
    nb = nb+size(joints{ii},1)-1;
end

if isfield(model,'harmonics')
    H = model.harmonics;
else
    H = zeros(nb,1);
end
nh = size(H,2);

if length(Nel) ~= nb || length(width) ~= nb
    error('Nel and width must have the same size as the number of beams')
end

%% create nodes and elements
nodes = [];
node_info = [];
elements = [];
EL = 1;
count = 1; % count beams
for ii = 1 : length(joints)

    for jj = 1 : size(joints{ii},1)-1 % loop over pairs of joints

        x0 = joints{ii}(jj,  1);
        xf = joints{ii}(jj+1,1);
        y0 = joints{ii}(jj,  2);
        yf = joints{ii}(jj+1,2);
        I0 = joints{ii}(jj,  3);
        If = joints{ii}(jj+1,3);
        
        nel = Nel(count);

        xl = linspace(x0, xf, nel+1)';
        yl = linspace(y0, yf, nel+1)';
        info = [I0; zeros(length(xl)-2,1); If];
        
        if sum(H(count,:))~=0
            % apply a change of shape to the beams as a sum of harmonic,
            % sinusoidal components
            Lb = sqrt(diff(joints{ii}(jj:jj+1,1))^2 + diff(joints{ii}(jj:jj+1,2))^2);
            th = atan2((joints{ii}(jj+1,2)-joints{ii}(jj,2)) , ...
                (joints{ii}(jj+1,1)-joints{ii}(jj,1)));
            R = [cos(th), -sin(th); sin(th), cos(th)];
    
            xx = linspace(0, Lb, nel+1);
            yy = 0;
            for kk = 1 : nh
                yy = yy + H(count,kk) * sin(kk*pi/Lb*xx);
            end
            rr = R*[xx; yy];
            xl = x0 + rr(1,:)';
            yl = y0 + rr(2,:)';
        end

        if jj>1
            xl(1) = [];
            yl(1) = [];
            info(1) = [];
        end

        nodes = [nodes; [xl yl]]; %#ok<AGROW>
        node_info = [node_info; info]; %#ok<AGROW>
        
        % define elements with cells and different constructors to be able
        % to specify different properties for each beam
        elements{count} = [(EL:(EL+nel-1))' (EL+1:(EL+nel))'];
        element_constructor{count} = @()BeamElement(thickness, width(count), myMaterial);
        
        xm = mean(nodes((EL:(EL+nel))',1));
        ym = mean(nodes((EL:(EL+nel))',2));
        beam_tag(count,:) = [xm, ym];
        
        count = count + 1;
        EL = EL + nel;
    end
    EL = EL+1;
end

%% MPC assembly

% Create mesh
myMesh = Mesh(nodes);
myMesh.create_elements_table(elements, element_constructor);

% Constraints
constrained_nodes = find(node_info == -1);
myMesh.set_essential_boundary_condition(constrained_nodes, 1:3, 0);

% MPC information
nodes_master = [model.xg model.yg];
nodes_master_constr = model.mass_constraints;

masterNodes = [nodes_master nodes_master_constr];
slave_node_IDs  = find(node_info ~= -1 & node_info ~= 0);
master_IDs = node_info(node_info ~= -1 & node_info ~= 0);
mpc = [slave_node_IDs, master_IDs];

% lumped mass and stiffness matrices
Mm = cell(length(model.A),1);
for ii = 1 : length(model.A)
    Mm{ii} = zeros(3);
    Mm{ii}(1,1) = myMaterial.DENSITY*thickness*model.A(ii);
    Mm{ii}(2,2) = myMaterial.DENSITY*thickness*model.A(ii);
    Mm{ii}(3,3) = myMaterial.DENSITY*thickness*(model.Ixx(ii) + model.Iyy(ii));
end

% Build MPC object
obj  = MultiPointConstraint(myMesh, mpc, masterNodes, Mm);
obj.PLOT.mass_xy = model.xy;
obj.PLOT.elem = elements;
obj.PLOT.xy_tag = beam_tag;

if plot_flag == 1
    figure
    obj.mpc_plot
    title('MPC model')
end

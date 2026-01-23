function plot_geometry(sys, varargin)
    % Plot the MEMS geometry. If a displacement field is provided, the
    % deformed geometry is plotted.
    % Inputs:
    %   sys: the MPC mechanical system.
    %   u: the displacement field (optional, default is []).
    %   scale: the scale factor for the displacement field (optional,
    %       default is 1).

    % Parse inputs
    p = inputParser;
    addOptional(p, 'u', [])
    addOptional(p, 'scale', 1)
    parse(p, varargin{:})
    u = p.Results.u;
    scale = p.Results.scale;

    % Colors
    driveBeamColor = [217, 83, 25] / 255; % matlab red
    senseBeamColor = [0, 114, 189] / 255; % matlab blue
    driveFrameColor = [244, 198, 178] / 255; % matlab light red
    senseFrameColor = [156, 205, 239] / 255; % matlab light blue

    % Extract yafec MPC assembly and coordinates
    yafecAssembly = sys.yafecAssembly;
    x = yafecAssembly.Mesh.nodes(:, 1);
    y = yafecAssembly.Mesh.nodes(:, 2);
    mm = yafecAssembly.PLOT.mass_xy;

    % Update nodes coordinates if a displacement field is provided
    if ~isempty(u)
        [U, Um] = yafecAssembly.unconstrain_mpc_vector(u);
        a = max(abs([U(1:3:end); U(2:3:end); Um(1:3:end); Um(2:3:end)]));
        U  =  U/a*scale;
        Um = Um/a*scale;
        u = U(1:3:end);
        v = U(2:3:end);

        % Update coordinates
        x = x + u;
        y = y + v;

        % Loop over masses
        for i = 1:length(mm)
            % Extract nominal coordinates
            xm = mm{i}(:,1).';
            ym = mm{i}(:,2).';

            % Center of mass
            xG = mean(xm);
            yG = mean(ym);

            % Rotation around the center of mass
            th = -Um(i*3);
            R = [cos(th) sin(th); -sin(th) cos(th)];
            Q = R*[xm-xG; ym-yG];
            rx = Q(1,:);
            ry = Q(2,:);

            % Update mass coordinates
            mm{i} = [xG + rx + Um(i*3-2); yG + ry + Um(i*3-1)].';
        end
    end

    % Initialize figure
    hold on; grid on; box on; axis equal tight;

    % Plot beams
    elements = yafecAssembly.PLOT.elem;
    for i = 1:length(elements)
        % Extract drive beam nodes
        e = elements{i}.';
        e = unique(e(:));
        xDrive = x(e);
        yDrive = y(e);

        % Set color
        if i < 5
            color = driveBeamColor;
        else
            color = senseBeamColor;
        end

        % Plot beam
        plot(xDrive, yDrive, '.-', 'LineWidth', 2, 'Color', color);
    end

    % Plot masses
    patch(mm{1}(:, 1), mm{1}(:, 2), senseFrameColor);
    patch(mm{2}(:, 1), mm{2}(:, 2), driveFrameColor);

    % Extract constrained nodes from constrained DOFs
    constrainedDOFs = yafecAssembly.Mesh.EBC.constrainedDOFs;
    constrainedNodes = yafecAssembly.Mesh.get_nodeIDs_from_DOF(constrainedDOFs(:, 1));

    % Plot fixed nodes
    xs = yafecAssembly.Mesh.get_location_from_nodeIDs(constrainedNodes);
    plot(xs(:, 1), xs(:, 2), 'ko', 'MarkerFaceColor', 'y');
end

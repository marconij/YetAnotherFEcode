classdef ComsolTO < handle
    % A class to import topology optimization into COMSOL Multiphysics
    
    properties
        boundaries = {};
        model;
        name = 'comsol_model';
        box_counter = 0;
        fix_counter = 0;
        roll_counter = 0;
        load_counter = 0;
        study_counter = 0;
        solution_counter = 0;
        plot_counter = 0;
    end
    
    methods
        function obj = ComsolTO(coord, d, varargin)
            % Constructor.
            % Inputs:
            %   coord: a matrix (n x nDim) with the elements coordinates.
            %   d: the density vector.
            %   doPlot: a flag to plot the data (optional, default is false).
            % Outputs:
            %   obj: the ComsolTO object.

            % Parse the inputs
            p = inputParser;
            addOptional(p, 'doPlot', false, @islogical);
            parse(p, varargin{:});
            doPlot = p.Results.doPlot;

            % Sort coordinates and density
            nDim = size(coord, 2); % geometric dimension
            [coord, idx] = sortrows(coord, nDim:-1:1);
            d = d(idx);

            % Unique coordinates of the element centers
            elx = uniquetol(coord(:, 1), 1e-6);
            ely = uniquetol(coord(:, 2), 1e-6);

            % Element size
            dx = elx(2) - elx(1);
            dy = ely(2) - ely(1);

            % Number of elements
            nx = length(elx);
            ny = length(ely);
            dd = reshape(d, [nx, ny]).';

            % Plot density
            if doPlot
                figure;
                hold on; box on; axis equal tight;
                surf(elx, ely, dd, 'LineStyle', 'none')
                view(2)
                colormap(gca, flip(gray))
                set(gca, 'FontSize', 20, 'xtick', [], 'ytick', []);
            end

            % Build level-set
            ls = zeros(size(dd) + 1);
            ls(1:end-1, 1:end-1) = ls(1:end-1, 1:end-1) + dd;
            ls(2:end, 1:end-1) = ls(2:end, 1:end-1) + dd;
            ls(1:end-1, 2:end) = ls(1:end-1, 2:end) + dd;
            ls(2:end, 2:end) = ls(2:end, 2:end) + dd;
            ls = (ls / 4) * 2 - 1;

            % Zero padding for visualization
            x = -1:(nx+1);
            y = -1:(ny+1);
            zz = -ones(ny + 3, nx + 3);
            zz(2:end-1, 2:end-1) = ls;

            % Plot level set
            if doPlot
                figure;
                hold on; box on; axis equal tight;
                contourf(x, y, zz, [0, 0], 'color', 'k', 'facecolor', 'k');
                set(gca, 'FontSize', 20, 'xtick', [], 'ytick', []);
            end

            % Extract boundary points
            level = 0; % level set value to extract the boundaries
            doSimplify = true; % flag to enable/disable simplification of the boundaries
            tol = 1e-3; % default value for the Ramer–Douglas–Peucker algorithm
            obj.boundaries = obj.extract_boundaries(x, y, zz, level, ...
                'doSimplify', doSimplify, 'tol', tol);

            % Plot boundary points
            if doPlot
                figure;
                hold on; box on; axis equal tight;
                contourf(x, y, zz, [0, 0], 'color', 'k', 'facecolor', 'k');
                set(gca, 'FontSize', 20, 'xtick', [], 'ytick', []);

                % Add boundaries
                for i = 1:length(obj.boundaries)
                    plot3(obj.boundaries{i}(:, 1), obj.boundaries{i}(:, 2), ones(size(obj.boundaries{i}(:, 1))), '.r', 'MarkerSize', 20)
                end
            end

            % Convert boundary points in physical coordinates
            for i = 1:length(obj.boundaries)
                obj.boundaries{i} = obj.boundaries{i} .* [dx, dy];
            end
        end

        function write_points_to_file(obj, filename, varargin)
            % Write the coordinates of the boundary points to a file.
            % Inputs:
            %   filename: the name of the file w/o extension.
            %   scale: scaling factor for the coordinates.
            %   origin: origin point for the coordinates.

            % Parse the inputs
            p = inputParser;
            addOptional(p, 'scale', 1, @isnumeric);
            addOptional(p, 'origin', [0, 0], @isnumeric);
            parse(p, varargin{:});
            scale = p.Results.scale;
            origin = p.Results.origin;

            % Loop over boundaries
            for i = 1:length(obj.boundaries)
                % Get boundary points
                bp = origin + obj.boundaries{i} * scale;
                
                % Open file
                fid = fopen([filename, '_', num2str(i, '%2d'), '.txt'], 'w');
                if fid == -1
                    error('Cannot open file %s', [filename, '_', num2str(i, '%2d'), '.txt']);
                end

                % Write points to file
                for j = 1:size(bp, 1)
                    fprintf(fid, '%.6f %.6f\n', bp(j, 1), bp(j, 2));
                end

                % Close file
                fclose(fid);
            end
        end
        
        function create_model(obj, varargin)
            % Create the comsol model.
            % Inputs:
            %   modelName: the name of the model (optional, default is
            %       'comsol_model').
            %   lengthUnit: the length unit for the model (optional,
            %       default is 'm').
            %   meshSize: the size of the mesh elements. It's a number
            %   between 1 (extremely fine) and 9 (extremely coarse), with a
            %   default of 5 (normal).

            % Parse the inputs
            p = inputParser;
            addOptional(p, 'modelName', 'comsol_model');
            addOptional(p, 'lengthUnit', 'm');
            addOptional(p, 'meshSize', 5);
            parse(p, varargin{:});
            obj.name = p.Results.modelName;
            lengthUnit = p.Results.lengthUnit;
            meshSize = p.Results.meshSize;

            % Import comsol model
            import com.comsol.model.*
            import com.comsol.model.util.*

            % Create model
            obj.model = ModelUtil.create('Model');
            obj.model.label([obj.name, '.mph']);
    
            % Find max and min coordinates
            x1 = Inf; x2 = -Inf;
            y1 = Inf; y2 = -Inf;
            for i = 1:length(obj.boundaries)
                bpMin = min(obj.boundaries{i});
                x1 = min(x1, bpMin(1));
                y1 = min(y1, bpMin(2));
            
                bpMax = max(obj.boundaries{i});
                x2 = max(x2, bpMax(1));
                y2 = max(y2, bpMax(2));
            end
    
            % Parameters
            obj.model.param.set('x1', num2str(x1, '%.e'));
            obj.model.param.set('x2', num2str(x2, '%.e'));
            obj.model.param.set('y1', num2str(y1, '%.e'));
            obj.model.param.set('y2', num2str(y2, '%.e'));
            obj.model.param.set('lx', 'x2-x1');
            obj.model.param.set('ly', 'y2-y1');
    
            % Create component
            obj.model.component.create('comp1', true);
            obj.model.component('comp1').geom.create('geom1', 2);
            obj.model.component('comp1').mesh.create('mesh1');

            % Set length unit
            obj.model.component('comp1').geom('geom1').lengthUnit(lengthUnit);

            % Set view
            obj.model.component('comp1').view('view1').axis.set('xmin', 'x1 - 0.2*lx');
            obj.model.component('comp1').view('view1').axis.set('xmax', 'x2 + 0.2*lx');
            obj.model.component('comp1').view('view1').axis.set('ymin', 'y1 - 0.2*ly');
            obj.model.component('comp1').view('view1').axis.set('ymax', 'y2 + 0.2*ly');

            % Geometry timer
            tStart = tic;
            fprintf('\nBuilding geometry ...\n')
    
            % Add polygons
            n_polygons = length(obj.boundaries);
            pol_names = cell(n_polygons, 1);
            for i = 1:n_polygons
            
                % Print
                fprintf('\tPolygon %d out of %d\n', i, n_polygons)
            
                % Polygon name
                pol_names{i} = strcat(['pol', num2str(i)]);
            
                % Create new polygon
                obj.model.component('comp1').geom('geom1').create(pol_names{i}, 'Polygon');
            
                % Table of points in the correct format
                points_cell = cellstr(string(obj.boundaries{i}));
            
                % Give source table
                obj.model.component('comp1').geom('geom1').feature(pol_names{i}).set('source', 'table');
                obj.model.component('comp1').geom('geom1').feature(pol_names{i}).set('table', points_cell);
                obj.model.component('comp1').geom('geom1').run(pol_names{i});
            end
            
            % Set difference
            if length(obj.boundaries) > 1
                fprintf('\tSetting difference\n')
                obj.model.component('comp1').geom('geom1').create('dif1', 'Difference');
                obj.model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'pol1'});
                obj.model.component('comp1').geom('geom1').feature('dif1').selection('input2').set(pol_names(2:end));
                obj.model.component('comp1').geom('geom1').run('dif1');
            end
            
            % Build geometry
            obj.model.component('comp1').geom('geom1').run;
            fprintf('done in %.2f seconds\n', toc(tStart))
            
            % Mesh
            tStart = tic;
            fprintf('\nCreating mesh ...\n')
            obj.model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
            obj.model.component('comp1').mesh('mesh1').feature('ftri1').create('size1', 'Size');
            obj.model.component('comp1').mesh('mesh1').feature('ftri1').feature('size1').set('hauto', meshSize);
            obj.model.component('comp1').mesh('mesh1').run;
            fprintf('done in %.2f seconds\n', toc(tStart))

            % Print mesh info
            meshInfo = mphmeshstats(obj.model, 'mesh1');
            fprintf('%d domain elements and %d boundary elements\n', meshInfo.numelem([2, 1]))
        end

        function add_material_properties(obj, density, young, poisson)
            % Add material properties.
            % Inputs:
            %   density: string with the density value (exaple: '1[kg/m^3]').
            %   young: string with the Young's modulus value (exaple: '1[MPa]').
            %   poisson: string with the Poisson's ratio value (exaple: '0.3').

            obj.model.component('comp1').material.create('mat1', 'Common');
            obj.model.component('comp1').material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
            obj.model.component('comp1').material('mat1').propertyGroup('def').set('density', density);
            obj.model.component('comp1').material('mat1').propertyGroup('Enu').set('E', young);
            obj.model.component('comp1').material('mat1').propertyGroup('Enu').set('nu', poisson);
        end

        function add_physics(obj, type2d, lz)
            % Add physics.
            % Inputs:
            %   type2d: string with the 2D type of physics.
            %       Options: 'PlaneStrain', 'PlaneStress', 'GeneralizedPlaneStrain'.
            %   lz: string with the thickness value (exaple: '1[m]').

            obj.model.component('comp1').physics.create('solid', 'SolidMechanics', 'geom1');
            obj.model.component('comp1').physics('solid').prop('Type2D').set('Type2D', type2d);
            obj.model.component('comp1').physics('solid').prop('d').set('d', lz);
            obj.model.component('comp1').physics('solid').feature('lemm1').set('geometricNonlinearity', 'linear');
        end

        function add_selection_box(obj, coord, tol, label)
            % Add selection box.
            % Inputs:
            %   coord: the coordinates of the center of the box.
            %   tol: the half sizes of the box.
            %   label: the label of the box.

            % Update the box counter
            obj.box_counter = obj.box_counter + 1;
            boxName = ['box_', label];

            % Find max and min coordinates
            x1 = coord(1) - tol(1);
            x2 = coord(1) + tol(1);
            y1 = coord(2) - tol(2);
            y2 = coord(2) + tol(2);

            % Create selection box
            obj.model.component('comp1').selection.create(boxName, 'Box');
            obj.model.component('comp1').selection(boxName).set('entitydim', 1);
            obj.model.component('comp1').selection(boxName).set('xmin', num2str(x1));
            obj.model.component('comp1').selection(boxName).set('xmax', num2str(x2));
            obj.model.component('comp1').selection(boxName).set('ymin', num2str(y1));
            obj.model.component('comp1').selection(boxName).set('ymax', num2str(y2));
            obj.model.component('comp1').selection(boxName).set('condition', 'inside');
        end

        function add_fixed_boundary(obj, boxLabel)
            % Add boundary conditions.
            % Inputs:
            %   boxLabel: the label of the box where the boundary condition is applied.

            % Update fix counter
            obj.fix_counter = obj.fix_counter + 1;
            fixName = ['fix', num2str(obj.fix_counter)];

            % Create fixed boundary condition
            obj.model.component('comp1').physics('solid').create(fixName, 'Fixed', 1);
            obj.model.component('comp1').physics('solid').feature(fixName).selection.named(boxLabel);
        end

        function add_roller(obj, boxlabel)
            % Add roller boundary condition.
            % Inputs:
            %   boxLabel: the label of the box where the boundary condition is applied.

            % Update roller counter
            obj.roll_counter = obj.roll_counter + 1;
            rollName = ['roller_', num2str(obj.roll_counter)];

            % Create roller boundary condition
            obj.model.component('comp1').physics('solid').create(rollName, 'Roller', 1);
            obj.model.component('comp1').physics('solid').feature(rollName).selection.named(boxlabel);
        end

        function add_boundary_load(obj, boxLabel, force)
            % Add boundary load.
            % Inputs:
            %   boxLabel: the label of the box where the boundary load is applied.
            %   force: cell string with the force values (exaple: {'1[N]', '0[N]', '0[N]'}).

            % Update load counter
            obj.load_counter = obj.load_counter + 1;
            loadName = ['bndl_', num2str(obj.load_counter)];

            % Create boundary load
            obj.model.component('comp1').physics('solid').create(loadName, 'BoundaryLoad', 1);
            obj.model.component('comp1').physics('solid').feature(loadName).selection.named(boxLabel);
            obj.model.component('comp1').physics('solid').feature(loadName).set('LoadType', 'TotalForce');
            obj.model.component('comp1').physics('solid').feature(loadName).set('Ftot', force);
        end

        function volFrac = compute_volume_fraction(obj)
            % Compute volume fraction.
            volFrac = mphint2(obj.model, '1/(lx*ly)', 2);
            volFrac = volFrac(1);
        end

        function C = stationary_study(obj)
            % Perform stationary study.
            % Outputs:
            %   C: the compliance value.

            % Update counters
            obj.study_counter = obj.study_counter + 1;
            obj.solution_counter = obj.solution_counter + 1;
            obj.plot_counter = obj.plot_counter + 1;

            % Names
            stdName = ['std', num2str(obj.study_counter)];
            solName = ['sol', num2str(obj.solution_counter)];
            plotName = ['pg', num2str(obj.plot_counter)];
            dataName = ['dset', num2str(obj.solution_counter)];

            % Start timer
            tStart = tic;
            fprintf('\nStationary study ...\n')

            % Create study
            obj.model.study.create(stdName);
            obj.model.study(stdName).create('stat', 'Stationary');

            % Create solution
            obj.model.sol.create(solName);
            obj.model.sol(solName).study(stdName);
            obj.model.sol(solName).attach(stdName);
            obj.model.sol(solName).create('st1', 'StudyStep');
            obj.model.sol(solName).create('v1', 'Variables');
            obj.model.sol(solName).create('s1', 'Stationary');
            obj.model.sol(solName).feature('s1').create('fc1', 'FullyCoupled');
            obj.model.sol(solName).feature('s1').feature.remove('fcDef');

            % Print number of dofs
            info = mphxmeshinfo(obj.model, 'soltag', solName, 'studysteptag', 'st1');
            fprintf('number of degrees of freedom: %d ...\n', info.ndofs)

            % Run solution
            obj.model.sol(solName).attach(stdName);
            obj.model.sol(solName).feature('st1').label('Compile Equations: Stationary');
            obj.model.sol(solName).feature('v1').label('Dependent Variables 1.1');
            obj.model.sol(solName).feature('s1').label('Stationary Solver 1.1');
            obj.model.sol(solName).feature('s1').feature('dDef').label('Direct 1');
            obj.model.sol(solName).feature('s1').feature('aDef').label('Advanced 1');
            obj.model.sol(solName).feature('s1').feature('aDef').set('cachepattern', true);
            obj.model.sol(solName).feature('s1').feature('fc1').label('Fully Coupled 1.1');
            obj.model.sol(solName).runAll;

            % Create plot
            obj.model.result.create(plotName, 'PlotGroup2D');
            obj.model.result(plotName).create('surf1', 'Surface');
            obj.model.result(plotName).feature('surf1').set('expr', 'solid.misesGp');
            obj.model.result(plotName).feature('surf1').create('def', 'Deform');

            % Plot results
            obj.model.result(plotName).label('Stress (solid)');
            obj.model.result(plotName).set('data', dataName);
            obj.model.result(plotName).set('frametype', 'spatial');
            obj.model.result(plotName).feature('surf1').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x-coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y-coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z-coordinate'});
            obj.model.result(plotName).feature('surf1').set('colortable', 'Prism');
            obj.model.result(plotName).feature('surf1').set('threshold', 'manual');
            obj.model.result(plotName).feature('surf1').set('thresholdvalue', 0.2);
            obj.model.result(plotName).feature('surf1').set('resolution', 'normal');
            obj.model.result(plotName).feature('surf1').feature('def').set('scale', 632925.2448293482);
            obj.model.result(plotName).feature('surf1').feature('def').set('scaleactive', false);

            % Evaluate strain energy
            Ws_tot = mphglobal(obj.model, 'solid.Ws_tot');
            C = 2 * Ws_tot;

            % Print
            fprintf('done in %.2f seconds\n', toc(tStart))
        end

        function f = eigenfrequency_study(obj, shift, varargin)
            % Perform eigenfrequency study.
            % Inputs:
            %   shift: string with the shift value (exaple: '1[kHz]').
            %   neigs: number of eigenfrequencies (optional, default is 6).
            %   eigunit: string with the eigenfrequency unit (optional, default is 'Hz').
            % Outputs:
            %   f: the eigenfrequencies.

            % Parse the inputs
            p = inputParser;
            addOptional(p, 'eigunit', 'Hz');
            addOptional(p, 'neigs', 6);
            parse(p, varargin{:});
            eigunit = p.Results.eigunit;
            neigs = p.Results.neigs;

            % Update counters
            obj.study_counter = obj.study_counter + 1;
            obj.solution_counter = obj.solution_counter + 1;
            obj.plot_counter = obj.plot_counter + 1;

            % Names
            stdName = ['std', num2str(obj.study_counter)];
            solName = ['sol', num2str(obj.solution_counter)];
            plotName = ['pg', num2str(obj.plot_counter)];
            dataName = ['dset', num2str(obj.solution_counter)];

            % Start timer
            tStart = tic;
            fprintf('\nEigenfrequency study ...\n')

            % Settings
            obj.model.result.evaluationGroup.create('std1EvgFrq', 'EvaluationGroup');
            obj.model.result.evaluationGroup.create('std1mpf1', 'EvaluationGroup');
            obj.model.result.evaluationGroup('std1EvgFrq').create('gev1', 'EvalGlobal');
            obj.model.result.evaluationGroup('std1mpf1').create('gev1', 'EvalGlobal');
            obj.model.component('comp1').common.create('mpf1', 'ParticipationFactors');

            % Create study
            obj.model.study.create(stdName);
            obj.model.study(stdName).create('eig', 'Eigenfrequency');

            % Study settings
            obj.model.study(stdName).feature('eig').set('neigs', neigs);
            obj.model.study(stdName).feature('eig').set('neigsactive', true);
            obj.model.study(stdName).feature('eig').set('eigunit', eigunit);
            obj.model.study(stdName).feature('eig').set('shift', shift);

            % Create solution
            obj.model.sol.create(solName);
            obj.model.sol(solName).study(stdName);
            obj.model.sol(solName).attach(stdName);
            obj.model.sol(solName).create('st1', 'StudyStep');
            obj.model.sol(solName).create('v1', 'Variables');
            obj.model.sol(solName).create('e1', 'Eigenvalue');

            % Print number of dofs
            info = mphxmeshinfo(obj.model, 'soltag', solName, 'studysteptag', 'st1');
            fprintf('number of degrees of freedom: %d ...\n', info.ndofs)

            % Run solution
            obj.model.sol(solName).attach(stdName);
            obj.model.sol(solName).feature('st1').label('Compile Equations: Eigenfrequency');
            obj.model.sol(solName).feature('v1').label('Dependent Variables 1.1');
            obj.model.sol(solName).feature('e1').label('Eigenvalue Solver 1.1');
            obj.model.sol(solName).feature('e1').set('transform', 'eigenfrequency');
            obj.model.sol(solName).feature('e1').set('neigs', 10);
            obj.model.sol(solName).feature('e1').set('shift', '1[kHz]');
            obj.model.sol(solName).feature('e1').set('eigvfunscale', 'maximum');
            obj.model.sol(solName).feature('e1').set('eigvfunscaleparam', 1.1199999999999999E-9);
            obj.model.sol(solName).feature('e1').feature('dDef').label('Direct 1');
            obj.model.sol(solName).feature('e1').feature('aDef').label('Advanced 1');
            obj.model.sol(solName).feature('e1').feature('aDef').set('cachepattern', true);
            obj.model.sol(solName).runAll;

            % Results
            obj.model.result.evaluationGroup('std1EvgFrq').label(num2str(obj.study_counter, 'Eigenfrequencies (Study %d)'));
            obj.model.result.evaluationGroup('std1EvgFrq').set('data', dataName);
            obj.model.result.evaluationGroup('std1EvgFrq').set('looplevelinput', {'all'});
            obj.model.result.evaluationGroup('std1EvgFrq').feature('gev1').set('expr', {'2*pi*freq' 'imag(freq)/abs(freq)' 'abs(freq)/imag(freq)/2'});
            obj.model.result.evaluationGroup('std1EvgFrq').feature('gev1').set('unit', {'rad/s' '1' '1'});
            obj.model.result.evaluationGroup('std1EvgFrq').feature('gev1').set('descr', {'Angular frequency' 'Damping ratio' 'Quality factor'});
            obj.model.result.evaluationGroup('std1EvgFrq').feature('gev1').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x-coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y-coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z-coordinate'});
            obj.model.result.evaluationGroup('std1mpf1').label(num2str(obj.study_counter, 'Participation Factors (Study %d)'));
            obj.model.result.evaluationGroup('std1mpf1').set('data', dataName);
            obj.model.result.evaluationGroup('std1mpf1').set('looplevelinput', {'all'});
            obj.model.result.evaluationGroup('std1mpf1').feature('gev1').set('expr', {'mpf1.pfLnormX' 'mpf1.pfLnormY' 'mpf1.pfLnormZ' 'mpf1.pfRnormX' 'mpf1.pfRnormY' 'mpf1.pfRnormZ' 'mpf1.mEffLX' 'mpf1.mEffLY' 'mpf1.mEffLZ' 'mpf1.mEffRX'  ...
            'mpf1.mEffRY' 'mpf1.mEffRZ'});
            obj.model.result.evaluationGroup('std1mpf1').feature('gev1').set('unit', {'1' '1' '1' '1' '1' '1' 'kg' 'kg' 'kg' 'kg*m^2'  ...
            'kg*m^2' 'kg*m^2'});
            obj.model.result.evaluationGroup('std1mpf1').feature('gev1').set('descr', {'Participation factor, normalized, X-translation' 'Participation factor, normalized, Y-translation' 'Participation factor, normalized, Z-translation' 'Participation factor, normalized, X-rotation' 'Participation factor, normalized, Y-rotation' 'Participation factor, normalized, Z-rotation' 'Effective modal mass, X-translation' 'Effective modal mass, Y-translation' 'Effective modal mass, Z-translation' 'Effective modal mass, X-rotation'  ...
            'Effective modal mass, Y-rotation' 'Effective modal mass, Z-rotation'});
            obj.model.result.evaluationGroup('std1mpf1').feature('gev1').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x-coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y-coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z-coordinate'});
            obj.model.result.evaluationGroup('std1EvgFrq').run;
            obj.model.result.evaluationGroup('std1mpf1').run;

            % Plot
            obj.model.result.create(plotName, 'PlotGroup2D');
            obj.model.result(plotName).set('data', dataName);
            obj.model.result(plotName).create('surf1', 'Surface');
            obj.model.result(plotName).feature('surf1').create('def', 'Deform');
            obj.model.result(plotName).label('Mode Shape (solid)');
            obj.model.result(plotName).set('showlegends', false);
            obj.model.result(plotName).feature('surf1').set('const', {'solid.refpntx' '0' 'Reference point for moment computation, x-coordinate'; 'solid.refpnty' '0' 'Reference point for moment computation, y-coordinate'; 'solid.refpntz' '0' 'Reference point for moment computation, z-coordinate'});
            obj.model.result(plotName).feature('surf1').set('colortable', 'AuroraBorealis');
            obj.model.result(plotName).feature('surf1').set('threshold', 'manual');
            obj.model.result(plotName).feature('surf1').set('thresholdvalue', 0.2);
            obj.model.result(plotName).feature('surf1').set('resolution', 'normal');
            obj.model.result(plotName).feature('surf1').feature('def').set('scale', 3.702346150931581E7);
            obj.model.result(plotName).feature('surf1').feature('def').set('scaleactive', false);

            % Get eigenfrequencies
            egv = obj.model.result.evaluationGroup('std1EvgFrq').get('real');
            f = egv(:, 1);

            % Print
            fprintf('done in %.2f seconds\n', toc(tStart))
        end
    end

    methods(Static, Access = private)
        function boundaries = extract_boundaries(x, y, zz, level, varargin)
            % Extract boundaries from a level set function using contourc.
            % Inputs:
            %   x: x-coordinates of the grid.
            %   y: y-coordinates of the grid.
            %   zz: level set function values.
            %   level: level set value to extract the boundaries.
            %   doSimplify: flag to simplify the boundaries (optional,
            %       default is true).
            %   tol: tolerance for the Ramer–Douglas–Peucker algorithm
            %       (optional, default is 1e-3).
            % Outputs:
            %   boundaries: cell array with the boundaries. Each cell
            %       contains a matrix with the x and y coordinates of the
            %       boundary points.

            % Parse inputs
            p = inputParser;
            addOptional(p, 'doSimplify', true, @islogical);
            addOptional(p, 'tol', 1e-3, @isnumeric);
            parse(p, varargin{:});
            doSimplify = p.Results.doSimplify;
            tol = p.Results.tol;
            
            % Perform contourc
            bp = contourc(x, y, zz, [level, level]);
            
            % Extract points and connectivity
            boundaries = {};
            bp_x = [];
            bp_y = [];
            n_points_left = 0;
            for ii = 1:size(bp, 2)
                
                if n_points_left > 0
                    bp_x = [bp_x, bp(1, ii)];
                    bp_y = [bp_y, bp(2, ii)];
                    n_points_left = n_points_left - 1;
            
                    if n_points_left == 0
                        bp_temp.x = bp_x;
                        bp_temp.y = bp_y;
                        boundaries = [boundaries, [bp_x.', bp_y.']];
                    end
                else
                    bp_x = [];
                    bp_y = [];
                    n_points_left = bp(2, ii);
                end
            end
            
            % Remove trivial boundaries (god only knows why those are there)
            for ii = length(boundaries):-1:1
                temp = boundaries{ii};
                if size(temp, 1) < 4
                    boundaries(ii) = [];
                end
            end
            
            % Remove unnecessary points
            if doSimplify
                % % Plot points
                % figure
                % hold on; grid on; box on; axis equal tight; legend show;
                % for i = 1:length(boundaries)
                %     plot(boundaries{i}(:, 1), boundaries{i}(:, 2), '.b', ...
                %         'MarkerSize', 20, 'DisplayName', 'v0')
                % end
            
                % Loop over boundaries to remove points
                nx = 1; ny = 1;
                for ii = 1:length(boundaries)
                    nx = max([boundaries{ii}(:,1); nx]);
                    ny = max([boundaries{ii}(:,2); ny]);
                end
                tolll = 1/4;
                for ii = 1:length(boundaries)
                    % Schiaffo sul contorno i punti che sono quasi sul contorno
                    dove = abs(boundaries{ii}(:,1)) < tolll;
                    boundaries{ii}(dove,1) = 0;
                    dove = abs(boundaries{ii}(:,1)-nx) < tolll;
                    boundaries{ii}(dove,1) = nx;
                    dove = abs(boundaries{ii}(:,2)) < tolll;
                    boundaries{ii}(dove,2) = 0;
                    dove = abs(boundaries{ii}(:,2)-ny) < tolll;
                    boundaries{ii}(dove,2) = ny;
            
            
                    boundaries{ii} = reducepoly(boundaries{ii}, tol);
                end
            
                % % Plot points
                % for i = 1:length(boundaries)
                %     plot(boundaries{i}(:, 1), boundaries{i}(:, 2), '.r', ...
                %         'MarkerSize', 20, 'DisplayName', 'v1')
                % end
            end
        end
    end
end

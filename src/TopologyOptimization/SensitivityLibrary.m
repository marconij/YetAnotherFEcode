classdef SensitivityLibrary
    % SensitivityLibrary: A collection of sensitivity analysis methods
    % for topology optimization problems.
    % Methods:
    %   compliance: compliance sensitivity.
    %   compliance_out: output compliance sensitivity.
    %   frequency: frequency sensitivity.
    %   initialize_pG_u: the partial derivative of the geometric
    %       stiffness matrix with respect to the displacement field.
    %   buckling: buckling sensitivity.
    %   strain_energy_nonhomogeneous: strain energy sensitivity for
    %   non-homogeneous boundary conditions.

    methods (Static)
        function sens = compliance(myMesh, u, Ke, varargin)
            % Compute the compliance sensitivity.
            % Inputs:
            %   myMesh: the Mesh object.
            %   u: vector of size n x 1, the displacement field.
            %   Ke: matrix of size nE x nE, the element stiffness matrix.
            %   sensPh: vector of size m x 1, the sensitivity of the
            %       physical density with respect to the projected density
            %       (optional, default is ones(m, 1)).
            % Outputs:
            %   sens: vector of size m x 1, the compliance sensitivity.

            % Parse inputs
            p = inputParser;
            addOptional(p, 'sensPh', ones(myMesh.nElements, 1));
            parse(p, varargin{:});
            sensPh = p.Results.sensPh;

            % Loop over the elements
            sens = zeros(myMesh.nElements, 1);
            for ii = 1:myMesh.nElements
                iDOFs = myMesh.Elements(ii).Object.iDOFs;
                ue = u(iDOFs, 1);
                dK = sensPh(ii) * Ke;
                sens(ii, :) = -ue.' * dK * ue;
            end
        end

        function sens = compliance_out(myAssembly, Kc, fOutC, u, Ke, varargin)
            % Compute the output compliance sensitivity.
            % Inputs:
            %   myAssembly: the Assembly object.
            %   uc: vector of size nC x 1, the displacement field.
            %   Kc: matrix of size nC x nC, the stiffness matrix (with the
            %       lumped spring).
            %   fOutC: vector of size nC x 1, the output force vector.
            %   u: vector of size n x 1, the displacement field.
            %   Ke: matrix of size nE x nE, the element stiffness matrix.
            %   sensPh: vector of size m x 1, the sensitivity of the
            %       physical density with respect to the projected density
            %       (optional, default is ones(m, 1)).
            % Outputs:
            %   sens: vector of size m x 1, the compliance sensitivity.

            % Parse inputs
            p = inputParser;
            addOptional(p, 'sensPh', ones(myAssembly.Mesh.nElements, 1));
            parse(p, varargin{:});
            sensPh = p.Results.sensPh;

            % Adjoint variable
            adjC = -Kc \ fOutC;
            adj = myAssembly.unconstrain_vector(adjC);

            % Loop over the elements
            sens = zeros(myAssembly.Mesh.nElements, 1);
            for ii = 1:myAssembly.Mesh.nElements
                iDOFs = myAssembly.Mesh.Elements(ii).Object.iDOFs;
                uE = u(iDOFs, 1);
                adjE = adj(iDOFs, 1);
                dK = sensPh(ii) * Ke;
                sens(ii, :) = adjE.' * dK * uE;
            end
        end

        function sens = frequency(myMesh, omega, u, Ke, Me, varargin)
            % Compute the frequency sensitivity.
            % Inputs:
            %   myMesh: the Mesh object.
            %   omega: scalar, the eigenfrequency.
            %   u: vector of size n x 1, the mode shape.
            %   Ke: matrix of size nE x nE, the element stiffness matrix.
            %   Me: matrix of size nE x nE, the element mass matrix.
            %   sensPhK: vector of size m x 1, the sensitivity of the
            %       physical density with respect to the projected one for
            %       the stiffness matrix (optional, default is ones(m, 1)).
            %   sensPhM: vector of size m x 1, the sensitivity of the
            %       physical density with respect to the projected one for
            %       the mass matrix (optional, default is ones(m, 1)).
            % Outputs:
            %   sens: vector of size m x 1, the frequency sensitivity.
            % Notes:
            %   Make sure the mode shape u is mass-normalized.

            % Parse inputs
            p = inputParser;
            addOptional(p, 'sensPhK', ones(myMesh.nElements, 1));
            addOptional(p, 'sensPhM', ones(myMesh.nElements, 1));
            parse(p, varargin{:});
            sensPhK = p.Results.sensPhK;
            sensPhM = p.Results.sensPhM;

            % Loop over the elements
            sens = zeros(myMesh.nElements, 1);
            for ii = 1:myMesh.nElements
                iDOFs = myMesh.Elements(ii).Object.iDOFs;
                ue = u(iDOFs, 1);
                dK = sensPhK(ii) * Ke;
                dM = sensPhM(ii) * Me;
                sens(ii, :) = ue.' * (dK - omega^2 * dM) * ue;
            end
            sens = sens / (2 * omega) / (2 * pi);
        end

        function pG_u = initialize_pG_u(myAssembly)
            % Compute the partial derivative of the geometric stiffness
            % matrix with respect to the displacement field.
            % Inputs:
            %   myAssembly: the Assembly object.
            % Outputs:
            %   pG_u: cell array of size nE x 1, where each cell contains
            %       a matrix of size nDOFs x nDOFs, the partial derivative
            %       of the geometric stiffness matrix with respect to the
            %       displacement field for each element.

            elIdx = 1;
            nE = myAssembly.Mesh.Elements(elIdx).Object.nelDOFs;
            iDOFs = myAssembly.Mesh.Elements(elIdx).Object.iDOFs;
            pG_u = cell(nE, 1);
            for ii = 1:nE
                u = zeros(myAssembly.Mesh.nDOFs, 1);
                u(iDOFs(ii)) = 1;  % Perturb the i-th DOF
                pG_u{ii} = myAssembly.Mesh.Elements(elIdx).Object.geometric_stiffness_matrix(u);
            end
        end

        function sens = buckling(myAssembly, Kc, P, u, v, Ke, pG_u, d, varargin)
            % Compute the buckling sensitivity.
            % Inputs:
            %   myAssembly: the Assembly object.
            %   Kc: matrix of size nC x nC, the stiffness matrix.
            %   P: scalar, the buckling load factor.
            %   u: vector of size n x 1, the displacement field.
            %   v: vector of size n x 1, the mode shape.
            %   Ke: matrix of size nE x nE, the element stiffness matrix.
            %   pG_u: cell array of size nE x 1, the partial derivative of
            %       the geometric stiffness matrix with respect to the
            %       displacement field.
            %   d: vector of size nElements x 1, the physical density.
            %   sensPhK: vector of size m x 1, the sensitivity of the
            %       physical density with respect to the projected one for
            %       the stiffness matrix (optional, default is ones(m, 1)).
            %   sensPhG: vector of size m x 1, the sensitivity of the
            %       physical density with respect to the projected one for
            %       the geometric stiffness matrix (optional, default is
            %       ones(m, 1)).
            % Outputs:
            %   sens: vector of size m x 1, the buckling sensitivity.
            % Notes:
            %   Make sure the mode shape v is mass-normalized.

            % Parse inputs
            p = inputParser;
            addOptional(p, 'sensPhK', ones(myAssembly.Mesh.nElements, 1));
            addOptional(p, 'sensPhG', ones(myAssembly.Mesh.nElements, 1));
            parse(p, varargin{:});
            sensPhK = p.Results.sensPhK;
            sensPhG = p.Results.sensPhG;

            % Create right-hand side of the adjoint equation
            rhs = zeros(myAssembly.Mesh.nDOFs, 1);
            for ii = 1:myAssembly.Mesh.nElements
                iDOFs = myAssembly.Mesh.Elements(ii).Object.iDOFs;
                vE = v(iDOFs, 1);
                for jj = 1:length(iDOFs)
                    rhs(iDOFs(jj)) = rhs(iDOFs(jj)) - d(ii) * vE.' * pG_u{jj} * vE;
                end
            end
            rhs = P * rhs;
            rhsC = myAssembly.constrain_vector(rhs);

            % Adjoint variable
            adjC = Kc \ rhsC;
            adj = myAssembly.unconstrain_vector(adjC);

            % Loop over the elements
            sens = zeros(myAssembly.Mesh.nElements, 1);
            for ii = 1:myAssembly.Mesh.nElements
                iDOFs = myAssembly.Mesh.Elements(ii).Object.iDOFs;
                uE = u(iDOFs, 1);
                vE = v(iDOFs, 1);
                adjE = adj(iDOFs, 1);
                dK = sensPhK(ii) * Ke;
                dG = sensPhG(ii) * myAssembly.Mesh.Elements(ii).Object.geometric_stiffness_matrix(u);
                sens(ii, :) = adjE.' * dK * uE + vE.' * (dK + P * dG) * vE;
            end
        end

        function sens = strain_energy_nonhomogeneous(myAssembly, u, Ke, varargin)
            % Compute the strain energy sensitivity for non-homogeneous
            % boundary conditions.
            % Inputs:
            %   myAssembly: the Assembly object.
            %   u: vector of size n x 1, the displacement field.
            %   Ke: matrix of size nE x nE, the element stiffness matrix.
            %   sensPh: vector of size m x 1, the sensitivity of the
            %       physical density with respect to the projected density
            %       (optional, default is ones(m, 1)).
            % Outputs:
            %   sens: vector of size m x 1, the strain energy sensitivity.

            % Parse inputs
            p = inputParser;
            addOptional(p, 'sensPh', ones(myAssembly.Mesh.nElements, 1));
            parse(p, varargin{:});
            sensPh = p.Results.sensPh;

            % If the totoal energy is defined as
            % J = 1/2 * u.' * K * u - fU.' * uU,
            % then the partial derivative of J with respect to the
            % displacement field u is pJ_u = 0.
            % In this case, the adjoint variable is zero, and the
            % sensitivity is computed directly as the partial derivative
            % of J with respect to the density field.

            % Loop over the elements
            sens = zeros(myAssembly.Mesh.nElements, 1);
            for ii = 1:myAssembly.Mesh.nElements
                iDOFs = myAssembly.Mesh.Elements(ii).Object.iDOFs;
                uE = u(iDOFs, 1);
                dK = sensPh(ii) * Ke;
                sens(ii, :) = 0.5 * uE.' * dK * uE;
            end
        end
    end
end

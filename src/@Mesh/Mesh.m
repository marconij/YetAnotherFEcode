classdef Mesh < handle
    properties
        nElements               % total number of elements in the mesh (including natural boundary elements)
        nNodes                  % total number of nodes in the mesh
        nodes = []              % Matrix containing nodal coordinates, row-wise for each nodes
        EBC = []                % EssentialBoundary class object
        nDOFPerNode = 0         % number of degrees of freedom associated to each node
        nDOFs = 0
        nDim = 0                % no. of dimensions in the physical domain    
        Elements = []           % A structure array containing pointer to 
                                % each element and an indicator which is true
                                % when element belongs to a natural
                                % boundary  and false when it belongs to
                                % the domain of the Mesh.
        elementsConnectivity    % Elements connectivity
        DATA                    % Miscellaneous data structure with 
                                % arbitrary, user-defined fields
    end
    
    methods
        function self = Mesh(varargin)
            narginchk(0,1)            
            if nargin==1
                self.nodes = varargin{1};
            end
        end
        
        function [nNodes] = get.nNodes(self)
            nNodes = size(self.nodes,1);
        end
        
        function nDim = get.nDim(self)
            nDim = size(self.nodes,2);
        end
        
        function nElements = get.nElements(self)
            nElements = length(self.Elements);
        end
        
        function [nDOFs] = get.nDOFs(self)
            nDOFs = self.nNodes * self.nDOFPerNode;
        end
        
        function dofs = get_element_indices(self,j)
            % this function computes the indices corresponding to the j-th
            %  element DOFs in the full vector of unknowns
            dofs = self.Elements(j).iDOFs;
        end
        
        function DOFs = get_DOF_from_nodeIDs(self,nodeIDs)
            nodeIDs = reshape(nodeIDs,[],1);
            DOFs = (nodeIDs-1)*self.nDOFPerNode + (1:self.nDOFPerNode);
        end
        
        function nodeIDs = get_nodeIDs_from_DOF(self,DOFs)
            % Returns the node IDs corresponding to the DOFs
            DOFs = reshape(DOFs,[],1);
            nodeIDs = ceil(DOFs/self.nDOFPerNode);
            nodeIDs = unique(nodeIDs); % remove duplicates
        end

        function DOFs = get_DOF_from_location(self,outcoord)
            DOFs = zeros(size(outcoord,1), self.nDOFPerNode);
            for ii = 1 : size(outcoord,1)
                dist = vecnorm(self.nodes - repmat(outcoord(ii,:),[self.nNodes,1]),2,2);
                [~,outnode] = min(dist);            
                DOFs(ii,:) = (outnode-1)*self.nDOFPerNode + (1:self.nDOFPerNode);
            end
        end

        function nodeIDs = get_nodeIDs_from_location(self,outcoord)
            nodeIDs = zeros(size(outcoord,1), 1);
            for ii = 1 : size(outcoord,1)
                dist = vecnorm(self.nodes - repmat(outcoord(ii,:),[self.nNodes,1]),2,2);
                [~,ind] = min(dist);            
                nodeIDs(ii,1) = ind;
            end
        end

        function nodeIDs = get_nodeIDs_inside_box(self,box)
            % box: 2 x nDim matrix defining the two opposite vertices of a
            % rectangle (in 2D) or a parallelpiped (in 3D). The IDs of the
            % nodes inside the box are returned.
            x1 = box(1,1); x2 = box(2,1);
            y1 = box(1,2); y2 = box(2,2);
            x = self.nodes(:,1); 
            y = self.nodes(:,2);
            if numel(box)==4
                nodeIDs = x>=x1 & x <= x2 & y >= y1 & y <= y2;
            elseif numel(box)==6
                z1 = box(1,3); z2 = box(2,3);
                z = self.nodes(:,3);
                nodeIDs = x>=x1 & x <= x2 & y >= y1 & y <= y2 & z >= z1 & z <= z2;
            else
                error('Wrong box dimension (must be 2 x nDim)')
            end
            nodeIDs = find(nodeIDs);
            if isempty(nodeIDs)
                warning('No nodes inside this box!')
            end
        end

        function nodeIDs = get_nodeIDs_from_elemIDs(self,elemIDs)
            nodeIDs = cell(length(elemIDs));
            for e = 1 : length(elemIDs)
                nodeIDs{e}(:,1) = self.Elements(e).Object.nodeIDs;
            end
            nodeIDs = vertcat(nodeIDs{:});
            nodeIDs = unique(nodeIDs);
        end

        function elemIDs = get_elemIDs_from_nodeIDs(self,nodeIDs)
            % input: vector containing node IDs
            % output: IDs of the elements whose all nodes are also in nodeIDs
            A = self.elementsConnectivity;
            B = nodeIDs;
            ic = zeros(self.nElements,1);
            for ii = 1 : self.nElements
                C = setdiff(A(ii,:),B);
                ic(ii) = isempty(C);
            end
            elemIDs = find(ic);
        end

        function [elemIDs,nodeIDs] = get_elemIDs_inside_box(self,box)
            nodeIDs = self.get_nodeIDs_inside_box(box);
            elemIDs = self.get_elemIDs_from_nodeIDs(nodeIDs);
        end

        function coordinates = get_location_from_nodeIDs(self,nodeIDs)
            nodeIDs = nodeIDs(:);
            coordinates = self.nodes(nodeIDs,:);
        end
        
        function create_elements_table(self, elementsConnectivity,elementConstructor,varargin)
            % creates element_table containing element-wise data
            % Elements(j).Nodes = list of node id for j-th element
            %                       (different number of nodes per element is allowed)
            % Elements(j).Object = Element object for j-th element
            % Elements(j).isBoundary = true if j-th Element is Boundary element, false otherwise
            
            narginchk(3,5)

            self.elementsConnectivity = elementsConnectivity;
            %% check if there are multiple element types
            multielem = 0;
            if iscell(elementsConnectivity)
                ntypes = length(elementsConnectivity);
                if ntypes > 1
                    multielem = 1;
                end
            else
                ntypes = 1;
                econn{1} = elementsConnectivity;
                econstr{1} = elementConstructor;
                elementsConnectivity = econn;
                elementConstructor = econstr;
            end
            %% check element compatibility with the mesh
            for ee = 1 : ntypes
                test = elementConstructor{ee}();
                if self.nDOFPerNode==0 % default value,                
                    self.nDOFPerNode = test.nDOFPerNode;
                else    % check nDOFPerNode compatibility between elementConstructor and the Mesh
                    if ~isequal(self.nDOFPerNode, test.nDOFPerNode)
                        error('nDOFPerNode error: The supplied element must have the same number of DOFs per node (nDOFPerNode) as the Mesh')
                    end
                end
            end
            %% parse inputs
            % check whether these are natural boundary elements or not           
            
            p = inputParser;
            addParameter(p,'isBoundary',false,@(x)validateattributes(x,{'logical'},{'nonempty'}))
            parse(p,varargin{:});
            isBoundary = p.Results.isBoundary;
            
            %% Add the supplied elements to given list
            for ee = 1 : ntypes
                nElements = size(elementsConnectivity{ee},1); %#ok<*PROPLC>
                if ~isempty(elementsConnectivity{ee})
                    % initialize elements as a table
                    thisElementsTable = struct('Object', cell(nElements,1), 'isBoundary', isBoundary);
                    % node IDs for each element
                    for j = 1 : nElements
                        nodeIDs = elementsConnectivity{ee}(j,:);
                        thisElementsTable(j).Object = elementConstructor{ee}();
                        thisElementsTable(j).Object.nodeIDs = nodeIDs;
                        thisElementsTable(j).Object.nodes = self.nodes(nodeIDs,:);
                    end
                    % add new elements to exisiting  (possibly empty) set of elements
                    self.Elements = [self.Elements; thisElementsTable];
                end
            end
        end
        
        function set_essential_boundary_condition(self,DirichletNodes,ConstrainedDOF,value)
            if isempty(self.EBC)
                self.EBC = EssentialBoundary(self.nDOFs,self.nDOFPerNode);
            end            
            self.EBC.apply_Dirichlet_BC(DirichletNodes,ConstrainedDOF,value);
        end
        
        function reset_boundary(self)
            self.EBC = [];
            disp('Removed Dirichlet Boundary DOFs')
            
            % find boundary elements
            I = find([self.Elements.isBoundary]);            
            if ~isempty(I)
                disp('Following boundary elements will be deleted')
                disp(I)
                self.Elements(I) = [];
            else
                disp('No boundary elements present')
            end
        end
        
    end
    
end

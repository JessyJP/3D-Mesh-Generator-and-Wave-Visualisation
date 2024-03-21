%
% Class-name:   MeshGenerator 
% Description: This script generates 2D/3D mesh and visualizes it.
% 
% Author: JessyJP (2022) % License: GPLv3 @ LICENCE.md
%

%% Mesh generator
classdef MeshGenerator < handle
    % MeshGenerator Class Implementation
    % This script generates 2D/3D mesh and visualizes it using CubeVoxel instances.
    
    % Internal properties
    properties (SetAccess=private)
        cubeTemplate;% This is a cube prototype
        elementAlignemntTypes;
        elements;% Store the elements
        xyz;% The dimension and setup of the mesh
    end

    % Public properties
    properties (Access=public)
        cube % Working cube modified
        dimRemoveID % 0=>OFF, 1=>x, 2=>y, 3=>z
        Pplot % Plot settings
        fig % Figure for plotting
    end
    
    methods
        % Constructor
        function [mesh] = MeshGenerator()
            % Settings
            mesh.dimRemoveID = 0;% 0=>OFF 1=>x 2=>y 3=>z
            % Labels
            xyz(1).Label='x';
            xyz(2).Label='y';
            xyz(3).Label='z';
            for i = 1:numel(xyz)
                % Constant IDs
                xyz(i).ID=i;
                % Mesh properties: Default Limits
                xyz(i).Lim = [0 1]; 
                % Mesh properties: Default Number of Elements
                xyz(i).N = 1;
                % Element Gap factor between elements
                xyz(i).gapP= 0.1;% Input 
                % Element side length
                xyz(i).ElementL = mesh.calcElSpacing(xyz(i))*(1-xyz(i).gapP); 
                % Element side length gap
                xyz(i).gap =  mesh.calcElSpacing(xyz(i))*(xyz(i).gapP);
                % Mesh grid coordinates
                xyz(i).Coords = {};
                % Mesh grid coordinates
                xyz(i).mesh =[];
            end

            % Element Alignment
            mesh.elementAlignemntTypes = ["Symetric";"Center"];
            % Default Plot view direction
            mesh.Pplot = struct( ...
                'View', 3, 'gridON', true, 'scaleFactor', 1, 'axisEqual', true, ...
                'limX', [0, 1], 'limY', [0, 1], 'limZ', [0, 1]);


            % Assign Mesh
            mesh.xyz = struct2table(xyz,'RowNames',{xyz.Label});


            % Initialize cube template with default CubeVoxel
            mesh.cubeTemplate = CubeVoxel(0, [], [], [], [], 1.0);
            mesh.cube = mesh.cubeTemplate; % Copy the cube template to the public cube property            
        end
        
        % Setup mesh
        function setupMesh(mesh,xyz_)% mesh obj
            % Update the mesh parameters 
            if nargin == 2 && istable(xyz_)
                mesh.xyz = xyz_;
%                 disp(mesh.xyz);
                disp('Table updated!');
            else
                disp('Table not updated!');
            end
            
            % Setup Properties and Update colours
            mesh.cubeTemplate.updateColourPropertiesFrom(mesh.cube);
            % Update the colours of the template cube
            mesh.cube = mesh.cubeTemplate;
            
            % Default viewing Limits 
            mesh.Pplot.limX = mesh.xyz{'x','Lim'}*mesh.Pplot.scaleFactor;
            mesh.Pplot.limY = mesh.xyz{'y','Lim'}*mesh.Pplot.scaleFactor;
            mesh.Pplot.limZ = mesh.xyz{'z','Lim'}*mesh.Pplot.scaleFactor;                       
            
            % Make coordinates
            mesh.xyz('x',:) = mesh.setupCoordinate(mesh.xyz('x',:));
            mesh.xyz('y',:) = mesh.setupCoordinate(mesh.xyz('y',:));
            mesh.xyz('z',:) = mesh.setupCoordinate(mesh.xyz('z',:));

            % Remove the dimension component 
            if mesh.dimRemoveID > 0 && mesh.dimRemoveID <= 3
                % Remove the cube nodes that correspond to the removed dimension   
                ind = not(mesh.cube.nodesR(:,mesh.dimRemoveID) > 0);
                mesh.cube.nodesR = mesh.cube.nodesR(ind,:);

                % Get the right square face
                % For removing For X
                if mesh.dimRemoveID == 1                    
                    mesh.cube.FaceConnectNodes = [1 2 4 3];
                end
                % For removing For Y
                if mesh.dimRemoveID == 2                    
                    mesh.cube.FaceConnectNodes = [1 2 3 4];
                end
                % For removing For Z
                if mesh.dimRemoveID == 3
                    mesh.cube.FaceConnectNodes = [1 2 4 3];
                end
            end
            
            % Scale the cube
            mesh.cube = mesh.scale(mesh.cube,...
                [mesh.calcElSpacing(mesh.xyz('x',:))  ...
                 mesh.calcElSpacing(mesh.xyz('y',:)) ...
                 mesh.calcElSpacing(mesh.xyz('z',:))].*(1-mesh.xyz.gapP'), ...
                 [0 , 0, 0]);

            % Cube Test preview
%             mesh.PplotElements(mesh.translate(mesh.cube,[-1 -1 -1]*0.5*0));
            
            [X,Y,Z] = meshgrid(mesh.xyz{"x","Coords"}{:},mesh.xyz{"y","Coords"}{:},mesh.xyz{"z","Coords"}{:});
            % Assign Mesh
            mesh.xyz{"x","mesh"} = {X};
            mesh.xyz{"y","mesh"} = {Y};
            mesh.xyz{"z","mesh"} = {Z};
            
            Elements_ = [];
            offsetNodeInd = 0;
            % Make mesh
            for i = 1:numel(X)
                % Make elements
                el = mesh.cube.deepCopy();
                el.id = i;
                el.offsetNodeInd  = offsetNodeInd ;                
                % Update the face node index                
                offsetNodeInd = offsetNodeInd  + size(el.nodesR,1);
                
                % Translate the 
                el = el.translate([X(i),Y(i),Z(i)]);
                % For fast plotting
                el.FaceConnectNodes_Fast = el.FaceConnectNodes + (i-1)*(8*(mesh.dimRemoveID==0)+4*(mesh.dimRemoveID~=0)); 
%                 El.colourF_Fast     = repmat(El.colourF,   size(El.FaceConnectNodes_Fast,1),1);
%                 El.FaceAlpha_Fast  = repmat(El.FaceAlpha,size(El.FaceConnectNodes_Fast,1),1);
                % Gather elements together
                Elements_ = [Elements_;el];
            end
            mesh.elements = Elements_;
            % Plot Preview
%             mesh.PplotElements(mesh.translate(Elements_,[0 0 0]*0.5*0));


        end
                
        % Setup Coordinate
        function [c] = setupCoordinate(mesh,c)%s - general coordinates, P - properties

            % Setup Mesh properties: Element length
            c.ElementL = mesh.calcElSpacing(c)*(1-c.gapP);
            % Setup Mesh properties: Element gap
            c.gap      = mesh.calcElSpacing(c)*(c.gapP);

            % Remove a dimension for the grid
            c.Lim = c.Lim*not(mesh.dimRemoveID==c.ID);

            % Make a grid origins i.e. the [0 0 0] vertex of the Element cube polygon
            Coords = linspace(c.Lim(1),c.Lim(2),(c.N+1));
            Coords(end) =[]; 
            c.Coords = {Coords};
        end
        
        function size_element = calcElSpacing(~,c)
        % function size_element = sizeElement(~,c)
            size_element = (c.Lim(2)-c.Lim(1))./c.N;
        end
            
        % Subfunction for plotting element
        function plotElements(mesh,elements, plotEdgeAlpha)
            % Get handles and do plotting
            try 
                if not(isvalid(mesh.fig))
                    error();
                end
            catch
                mesh.fig = figure();
                mesh.fig.WindowStyle = 'docked';
            end
            clf(mesh.fig);
            view(mesh.Pplot.View);
            ax = mesh.fig.Children(1);
                
            if nargin < 3
                plotEdgeAlpha = false;
            end

            % Get axes. Setup plot limits and labels
%             ax = mesh.ax;
            if mesh.Pplot.axisEqual
                axis(ax,'equal');
            end
            ax.XLim = mesh.Pplot.limX;
            ax.YLim = mesh.Pplot.limY;
            ax.ZLim = mesh.Pplot.limZ;
            ax.XLabel.String = 'x';
            ax.YLabel.String = 'y';
            ax.ZLabel.String = 'z';
            if mesh.Pplot.gridON
                grid(ax,'on')
            end
            mesh.Pplot.axP = [];
            
            N_elements = numel(elements);
            % Plot elements one by one
            for i = 1:N_elements
                E = elements(i);

                if (not(plotEdgeAlpha))
                    mesh.Pplot.axP(i) = patch(ax,'Faces', E.FaceConnectNodes, 'Vertices', E.nodesR, ...
                    'Facecolor', E.colourF,'FaceAlpha',E.faceAlpha,'EdgeColor',E.colourE); 
            %         drawnow
                else
                    % Plot faces
                    mesh.Pplot.axP(i) = patch(ax, 'Faces', E.FaceConnectNodes, 'Vertices', E.nodesR, ...
                          'Facecolor', E.colourF, 'FaceAlpha', E.faceAlpha, 'EdgeColor', 'none');
                    
                    % Plot edges with specified edge alpha
                    edgeAlpha = E.edgeAlpha; % Example edge alpha value
                    for j = 1:size(E.FaceConnectNodes, 1)
                        face = E.FaceConnectNodes(j, :);
                        for k = 1:length(face)
                            % Determine the next vertex in the face (or loop back to the start)
                            nextIdx = mod(k, length(face)) + 1;
                            % Extract the coordinates for the current and next vertex
                            v1 = E.nodesR(face(k), :);
                            v2 = E.nodesR(face(nextIdx), :);
                            % Plot the edge as a line with specified transparency
                            line(ax, [v1(1), v2(1)], [v1(2), v2(2)], [v1(3), v2(3)], ...
                                 'Color', [E.colourE, edgeAlpha], 'LineWidth', 1.5);
                        end
                    end                    
                end
                printProgressBar(i,N_elements,"barLength",20)
            end
            drawnow;
        end

        % Subfunction for plotting element
        function plotMeshFast(mesh,elements)
            % Get handles and do plotting
            try 
                if not(isvalid(mesh.fig))
                    error;
                end
            catch
                mesh.fig = figure();
                mesh.fig.WindowStyle = 'docked';
            end
            clf(mesh.fig);
            view(mesh.Pplot.View);
            ax = mesh.fig.Children(1);
            % Get axes and setup plot limits
%             ax = mesh.ax;
            if mesh.Pplot.axisEqual
                axis(ax,'equal');
            end
            ax.XLim = mesh.Pplot.limX;
            ax.YLim = mesh.Pplot.limY;
            ax.ZLim = mesh.Pplot.limZ;
            ax.XLabel.String = 'x';
            ax.YLabel.String = 'y';
            ax.ZLabel.String = 'z';
            if mesh.Pplot.gridON
                grid(ax,'on')
            end
            mesh.Pplot.axP = [];
            
            % Get properties ax matrix
            FaceConnectNodes_Fast = mesh.getMatrix(elements,"FaceConnectNodes_Fast");            
            nodesR = mesh.getMatrix(elements,"nodesR");
%             colour =  mesh.getMatrix(elements,"Colour_Fast");
%             faceAlpha = mesh.getMatrix(elements,"FaceAlpha_Fast");
            
            mesh.Pplot.axP = ...
            patch(ax,'Faces', FaceConnectNodes_Fast, 'Vertices', nodesR,...
            'FaceColor',mesh.cube.colourF,'FaceAlpha',mesh.cube.faceAlpha, 'EdgeColor', mesh.cube.colourE); 
            drawnow;
        end
        
        % Subfunction for translating element
        function [elements] = translate(~,elements,translationVector)
            for i = 1:numel(elements)
                elements(i) =  elements(i).translate(translationVector);
            end
        end

        % Subfunction for scaling element
        function [elements] = scale(~, elements, scaleFactors, centerPoint)
            % Applies scaling to each voxel in the mesh.
            % mesh: Array or cell array of CubeVoxel instances.
            % scaleFactors: A 3-element vector specifying scale factors along x, y, z axes.
            % centerPoint: Optional. A 3-element vector specifying the center point to scale around.
            %              If not provided, each voxel's center is used.
        
            if nargin < 4
                useVoxelCenter = true;
            else
                useVoxelCenter = false;
            end
        
            for i = 1:length(elements)
                if useVoxelCenter
                    elements(i) = elements(i).scale(scaleFactors);
                else
                    elements(i) = elements(i).scale(scaleFactors, centerPoint);
                end
            end
        end

        function [elements] = rotate(~, elements, angle, axis, pivotPoint)
            % Applies rotation to each voxel in the mesh.
            % mesh: Array or cell array of CubeVoxel instances.
            % angle: Rotation angle in degrees.
            % axis: 'x', 'y', or 'z' axis to rotate around.
            % pivotPoint: Optional. A 3-element vector specifying the pivot point.
            %             If not provided, each voxel's center is used.
        
            if nargin < 5
                useVoxelCenter = true;
            else
                useVoxelCenter = false;
            end
        
            for i = 1:length(elements)
                if useVoxelCenter
                    elements(i) = elements(i).rotate(angle, axis);
                else
                    elements(i) = elements(i).rotate(angle, axis, pivotPoint);
                end
            end
        end


        % Transformation functions
        function [elements] = transformElements(~,elements,transformationF)
        % function [elements] = transformElements(~,elements,transformationF)   
            for i = 1:numel(elements)
                % Translate vector        
                elements(i) = transformationF(elements(i));                
            end
        end
        
        function [elements] = transformNodesR(~,elements,transformationF)
        % function [elements] = transformNodesR(~,elements,transformationF)   
            for i = 1:numel(elements)
                % Transform Nodes
                elements(i).nodesR = transformationF(elements(i).nodesR);                
            end
        end
        
        function elementsCopy = copyElements(mesh)
            % Creates a deep copy of all elements in the mesh.
            % Returns an array of CubeVoxel instances that are deep copies of the original elements.
            
            % Initialize an empty array to hold the copies.
            % Depending on how elements are stored, you might need to use {} for a cell array.
            elementsCopy = repmat(CubeVoxel(0, [], [], [], [], 1.0), size(mesh.elements));
            
            % Iterate over each element and create a deep copy.
            for i = 1:length(mesh.elements)
                elementsCopy(i) = mesh.elements(i).deepCopy();
            end
        end
        
        % Get method for the max XYZ
        function xyz = get.xyz(mesh)
             xyz = mesh.xyz;
        end
        
        % Convert elements to a matrix of coordinates
        function M = getMatrix(~,elements,property)
        % function M = getMatrix(~,elements,property)
        % M = cell2mat({elements.(property)}');
            M = cell2mat({elements.(property)}');
        end
        
        function vertexMatrix = elementsToVertexMatrix(mesh)
            % Extract all nodes from each element and concatenate them into a single matrix
            allNodes = vertcat(mesh.elements.nodesR);
            
            % Assuming each vertex has a unique ID across the mesh, sort or organize
            % the vertices if necessary (this step may adjust based on how IDs are assigned and used)
            % For simplicity, this example assumes vertices are already in the correct order
            
            vertexMatrix = allNodes; % Each row is [x, y, z] for a vertex
        end

        function vertexMatrixToElements(mesh, vertexMatrix)
            % Starting index for vertices in the matrix
            startIndex = 1;
            
            for i = 1:length(mesh.elements)
                % Number of vertices for the current element
                numVertices = size(mesh.elements(i).nodesR, 1);
                
                % Extract the subset of the vertex matrix for the current element
                endIndex = startIndex + numVertices - 1;
                mesh.elements(i).nodesR = vertexMatrix(startIndex:endIndex, :);
                
                % Update the start index for the next element
                startIndex = endIndex + 1;
            end
        end


        function exportToFile(mesh, filename, format)
            % Validate the format
            validFormats = ["OBJ", "STL"];
            if ~ismember(upper(format), validFormats)
                error(['Unsupported format: ', format]);
            end
            
            % Ensure filename ends with the correct extension
            [path, name, ext] = fileparts(filename);
            if ~strcmpi(ext, ['.', lower(format)])
                filename = fullfile(path, join([name, '.', lower(format)],""));
            end
            
            % Open the file
            fileID = fopen(filename, 'w');
            if fileID == -1
                error(['Failed to open file: ', filename]);
            end
            
            % Write the file header based on the format
            switch upper(format)
                case "OBJ"
                    fprintf(fileID, '# OBJ file exported from MATLAB\n');
                case "STL"
                    fprintf(fileID, 'solid mesh\n');
            end
            
            % Write the mesh data
            for i = 1:numel(mesh.elements)
                switch upper(format)
                    case "OBJ"
                        fprintf(fileID, '%s', mesh.elements(i).toOBJString());
                    case "STL"
                        fprintf(fileID, '%s', mesh.elements(i).toSTLString());
                end
                printProgressBar(i,numel(mesh.elements),"barLength",20);
            end
            
            % Write the file footer based on the format
            if upper(format) == "STL"
                fprintf(fileID, 'endsolid mesh\n');
            end
            
            % Close the file
            fclose(fileID);
        end

    end
end



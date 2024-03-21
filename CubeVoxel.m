classdef CubeVoxel < handle
    % CubeVoxel Class Implementation
    % A class representing a voxel (3D pixel) in the form of a cube.
    % This class provides methods for translating, scaling, rotating,
    % and applying custom transformations to the cube voxel.
    
    properties(Access=public)
        id % Identifier for the voxel
        nodesR % List of node positions defining the voxel
        FaceConnectNodes % Connectivity list defining the faces of the voxel
        colourF = [1.0, 1.0, 1.0] % Face color of the voxel
        colourE = [0, 0, 0] % Edge color of the voxel
        faceAlpha = 1.0 % Transparency of the voxel faces
        edgeAlpha = 1.0 % Transparency of the voxel faces
        center % Center of the voxel
        FaceConnectNodes_Fast;
        offsetNodeInd = 0;
    end
    
    methods
        function obj = CubeVoxel(id, nodesR, FaceConnectNodes, colourF, colourE, faceAlpha)
            % Constructor for CubeVoxel class with default parameters
            
            % Default parameters
            defaultNodesR = [0, 0, 0; 0, 0, 1; 1, 0, 1; 1, 0, 0; ...
                             0, 1, 0; 0, 1, 1; 1, 1, 1; 1, 1, 0];
            defaultFaceConnectNodes = [1, 4, 3, 2; 5, 8, 7, 6; 1, 2, 6, 5; ...
                                       3, 4, 8, 7; 2, 3, 7, 6; 1, 5, 8, 4];
            defaultColourF = [204 204 0]/255;
            defaultColourE = [100 100 100]/255;
            defaultFaceAlpha = 1.0;
            
            % Assign provided parameters or use defaults
            obj.id = id;
            obj.nodesR = defaultNodesR;
            obj.FaceConnectNodes = defaultFaceConnectNodes;
            obj.colourF = defaultColourF;
            obj.colourE = defaultColourE;
            obj.faceAlpha = defaultFaceAlpha;
            
            if nargin >= 2 && ~isempty(nodesR)
                obj.nodesR = nodesR;
            end
            if nargin >= 3 && ~isempty(FaceConnectNodes)
                obj.FaceConnectNodes = FaceConnectNodes;
            end
            if nargin >= 4 && ~isempty(colourF)
                obj.colourF = colourF;
            end
            if nargin >= 5 && ~isempty(colourE)
                obj.colourE = colourE;
            end
            if nargin >= 6 && ~isempty(faceAlpha)
                obj.faceAlpha = faceAlpha;
            end
            
            % Calculate the center based on the nodes
            obj.center = mean(obj.nodesR, 1);
        end
        
        function obj = translate(obj, translation_vector)
            % Translates the voxel by a given vector
            obj.nodesR = bsxfun(@plus, obj.nodesR, translation_vector);
            obj.center = mean(obj.nodesR, 1);
        end
        
        function obj = scale(obj, scale_factors, centerPoint)
            % Scales the voxel by given scale factors along each axis
            % Optionally, scales around a specified center point.
            %
            % Parameters:
            %   scale_factors: A 3-element vector specifying scale factors along each axis (x, y, z).
            %   centerPoint (optional): A 3-element vector specifying the center point to scale around.
            %                           If not provided, the voxel's current center is used.
            
            if nargin < 3
                % If centerPoint is not provided, use the voxel's current center
                centerPoint = obj.center;
            end
            
            % Ensure centerPoint is a row vector for the subtraction and addition operations
            if iscolumn(centerPoint)
                centerPoint = centerPoint';
            end
            
            % Perform scaling around the specified center
            obj.nodesR = bsxfun(@minus, obj.nodesR, centerPoint); % Translate to origin based on centerPoint
            obj.nodesR = bsxfun(@times, obj.nodesR, scale_factors); % Scale
            obj.nodesR = bsxfun(@plus, obj.nodesR, centerPoint); % Translate back to original position based on centerPoint
            
            % Update the center based on the new nodesR
            obj.center = mean(obj.nodesR, 1);
        end
       
        function obj = rotate(obj, angle, axis, pivotPoint)
            % Rotates the voxel around a specified axis by a given angle
            % angle is in degrees
            % pivotPoint is an optional parameter specifying the pivot point for rotation
            
            if nargin < 4
                pivotPoint = obj.center; % Use the voxel's center as the default pivot point
            end
            
            angle_rad = deg2rad(angle); % Convert angle to radians
            rotation_matrix = eye(3); % Initialize rotation matrix
            
            % Define rotation matrices for each axis
            switch axis
                case 'x'
                    rotation_matrix = [1, 0, 0; 0, cos(angle_rad), -sin(angle_rad); 0, sin(angle_rad), cos(angle_rad)];
                case 'y'
                    rotation_matrix = [cos(angle_rad), 0, sin(angle_rad); 0, 1, 0; -sin(angle_rad), 0, cos(angle_rad)];
                case 'z'
                    rotation_matrix = [cos(angle_rad), -sin(angle_rad), 0; sin(angle_rad), cos(angle_rad), 0; 0, 0, 1];
                otherwise
                    error('Axis must be ''x'', ''y'', or ''z''.');
            end
            
            % Translate nodes to the origin based on pivotPoint, apply rotation, and translate back
            obj.nodesR = bsxfun(@minus, obj.nodesR, pivotPoint) * rotation_matrix; % Apply rotation
            obj.nodesR = bsxfun(@plus, obj.nodesR, pivotPoint); % Translate back
            
            % Update the center based on the new nodesR
            obj.center = mean(obj.nodesR, 1);
        end
        
        % Additional custom transformation methods can be added here

        function obj = updateColourPropertiesFrom(obj, otherVoxel)
            % Updates the colour properties from another CubeVoxel instance.
            %
            % Parameters:
            %   otherVoxel: Another instance of CubeVoxel from which to copy
            %   the colour properties.
            
            if nargin < 2 || ~isa(otherVoxel, 'CubeVoxel')
                error('The input must be an instance of CubeVoxel.');
            end
            
            % Update the colour properties
            obj.colourF = otherVoxel.colourF;
            obj.colourE = otherVoxel.colourE;
            obj.faceAlpha = otherVoxel.faceAlpha;
        end
        
        function newObj = deepCopy(obj)
            % Creates a deep copy of the CubeVoxel instance.
            newObj = CubeVoxel(0); % Use the default constructor
            
            % Manually copy all properties
            p = properties(obj);
            for i = 1:length(p)
                newObj.(p{i}) = obj.(p{i});
            end
        end
    end

    methods
        function objString = toOBJString(obj)
            % Generates an OBJ format string for the voxel.
            objString = "";
            for i = 1:size(obj.nodesR, 1)
                objString = objString + sprintf('v %f %f %f\n', obj.nodesR(i, :));
            end
            for i = 1:size(obj.FaceConnectNodes, 1)
                % OBJ indices are 1-based
                indices = obj.FaceConnectNodes(i, :);
                indices = indices + obj.offsetNodeInd;
                % Create a string for the face indices
                faceStr = sprintf('%d ', indices);
                % Trim the last space and add to objString
                objString = objString + sprintf('f %s\n', strtrim(faceStr));
            end
        end

        function stlString = toSTLString(obj)
            % Generates an STL format string for the voxel.
            stlString = "solid voxel\n";
            for i = 1:size(obj.FaceConnectNodes, 1)
                vertices = obj.nodesR(obj.FaceConnectNodes(i, :), :);
                % Assuming calculateNormal can handle more than 3 vertices
                % and returns a proper normal for the polygon
                normal = obj.calculateNormal(vertices);
                % Convert polygon face to triangles and add to STL string
                stlString = stlString + obj.polygonToTrianglesSTL(vertices, normal);
            end
            stlString = stlString + "endsolid voxel\n";
        end

        function trianglesSTL = polygonToTrianglesSTL(~, vertices, normal)
            % This function converts a polygon into a series of triangles
            % following the specified algorithm and returns STL string
            % representation of these triangles.
            trianglesSTL = "";
            if size(vertices, 1) < 3
                error('Not enough vertices to form a polygon.');
            end
            
            for i = 3:size(vertices, 1)
                trianglesSTL = trianglesSTL + sprintf('facet normal %f %f %f\n', normal);
                trianglesSTL = trianglesSTL + "    outer loop\n";
                % First vertex
                trianglesSTL = trianglesSTL + sprintf('        vertex %f %f %f\n', vertices(1, :));
                % (i-1)th vertex
                trianglesSTL = trianglesSTL + sprintf('        vertex %f %f %f\n', vertices(i-1, :));
                % ith vertex
                trianglesSTL = trianglesSTL + sprintf('        vertex %f %f %f\n', vertices(i, :));
                trianglesSTL = trianglesSTL + "    endloop\nendfacet\n";
            end
        end
            
        
       function normal = calculateNormal(~, vertices)
            % Calculates the normal vector for a face given its vertices.
            % Assumes the vertices are provided in a counterclockwise order.
            
            % Calculate vectors along two edges of the triangle
            edge1 = vertices(2, :) - vertices(1, :);
            edge2 = vertices(3, :) - vertices(1, :);
            
            % The normal vector is perpendicular to both edges
            normal = cross(edge1, edge2);
            
            % Normalize the vector to have a unit length
            normal = normal / norm(normal);
        end
    end
end

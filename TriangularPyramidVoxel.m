classdef TriangularPyramidVoxel < CubeVoxel
    % TriangularPyramidVoxel Class Implementation
    % A class representing a voxel in the form of a triangular pyramid,
    % inheriting from CubeVoxel and overriding the constructor to specify
    % the pyramid's geometry.
    
    methods
        function obj = TriangularPyramidVoxel(id,nodesR_, FaceConnectNodes_, colourF, colourE, faceAlpha)
            % Constructor for TriangularPyramidVoxel class with default parameters
            
            % Define the nodes for a triangular pyramid (tetrahedron)
            % Assuming the base lies on the XY plane and the apex points along the Z axis
            nodesR = [0, 0, 0;    % Base vertex 1
                      1, 0, 0;    % Base vertex 2
                      0.5, sqrt(3)/2, 0;    % Base vertex 3
                      0.5, sqrt(3)/6, sqrt(6)/3];   % Apex vertex
            
            % Define the faces connecting the nodes
            FaceConnectNodes = [1, 2, 3;  % Base face
                                1, 2, 4;  % Side face 1
                                2, 3, 4;  % Side face 2
                                3, 1, 4]; % Side face 3
            
            % Call the superclass constructor with the pyramid-specific geometry
            % Note: The CubeVoxel constructor is designed to accept geometry parameters,
            % so we pass the pyramid's geometry along with any color and transparency settings.
            obj@CubeVoxel(id, nodesR, FaceConnectNodes, colourF, colourE, faceAlpha);
            
            % Since the superclass constructor calculates the center, we don't need to recalculate it here
            % unless the pyramid's geometry requires a different method of calculation.
        end
        
        % If needed, override other methods here to customize behavior for the triangular pyramid
    end
end

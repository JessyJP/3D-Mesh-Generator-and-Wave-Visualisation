"""
CubeVoxel Class Implementation

Description:
    A class representing a voxel (3D pixel) in the form of a cube. This class provides methods
    for translating, scaling, rotating, and applying custom transformations to the cube voxel.
    Rotations can be performed around the voxel's center or an optional external center.

Author:
    JessyJP (2024)

License:
    GPLv3 - This code is released under the GNU General Public License v3.0.
    See LICENSE.md for more details.
"""

from collections import OrderedDict
import numpy as np
from copy import deepcopy

## TODO: it would be good to have unit testing for this class

class CubeVoxel:
    rotation_matrix_cache = OrderedDict()  # Static cache for rotation matrices, using OrderedDict
    cache_limit = 1000  # Limit for the cache size
    
    def __init__(self, id, nodesR=None, faceConnectNodes=None, colourF=[1.0, 1.0, 1.0], colourE=[0, 0, 0], faceAlpha=1.0):
        """
        Initializes a CubeVoxel instance with specified properties or defaults.

        Parameters:
        - id (int): Identifier for the voxel.
        - nodesR (list, optional): List of node positions defining the voxel. Defaults to a unit cube.
        - faceConnectNodes (list, optional): Connectivity list defining the faces of the voxel. Defaults to a unit cube.
        - colourF (list, optional): Face color of the voxel. Defaults to white.
        - colourE (list, optional): Edge color of the voxel. Defaults to black.
        - faceAlpha (float, optional): Transparency of the voxel faces. Defaults to 1.0 (opaque).
        """
        if nodesR is None:
            nodesR = [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0],
                      [0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]]
        if faceConnectNodes is None:
            faceConnectNodes = [[1, 4, 3, 2], [5, 8, 7, 6], [1, 2, 6, 5],
                                [3, 4, 8, 7], [2, 3, 7, 6], [1, 5, 8, 4]]
            
        self.id = id
        self.nodesR = np.array(nodesR, dtype=float)
        self.faceConnectNodes = np.array(faceConnectNodes) - 1  # Adjusting for Python's 0-indexing
        self.faceConnectNodes_Fast = self.faceConnectNodes # For fast render
        self.colourF = np.array(colourF)
        self.colourE = np.array(colourE)
        self.faceAlpha = faceAlpha
        self.center = np.mean(self.nodesR, axis=0)
    
    def copy(self):
        """
        Creates a deep copy of the CubeVoxel instance.

        Returns:
        - CubeVoxel: A new instance of CubeVoxel with identical properties.
        """
        # Use deepcopy to ensure all mutable objects are also copied.
        return deepcopy(self)

    # Modification methods

    def translate(self, translation_vector):
        """
        Translates the voxel by a given vector.

        Parameters:
        - translation_vector (list or np.array): The vector by which to translate the voxel.
        """
        self.nodesR += np.array(translation_vector)
        self.center = np.mean(self.nodesR, axis=0)

    def scale(self, scale_factors):
        """
        Scales the voxel by given scale factors along each axis.

        Parameters:
        - scale_factors (list or np.array): Scale factors for the x, y, and z axes.
        """
        self.nodesR *= np.array(scale_factors)
        self.center = np.mean(self.nodesR, axis=0)

    def rotate(self, angle, axis='z', center_of_rotation=None):
        """
        Rotates the voxel around a specified axis by a given angle. Optionally, a center of rotation can be specified.

        Parameters:
        - angle (float): Rotation angle in degrees.
        - axis (str): Axis of rotation ('x', 'y', or 'z').
        - center_of_rotation (list or np.array, optional): The center point around which to rotate. Defaults to the voxel's center.
        """
        angle_rad = np.radians(angle)
        rotation_center = self.center if center_of_rotation is None else np.array(center_of_rotation)
        rotation_matrix = CubeVoxel.get_rotation_matrix(angle_rad, axis)

        self.nodesR = np.dot(self.nodesR - rotation_center, rotation_matrix) + rotation_center
        self.center = np.mean(self.nodesR, axis=0)

    @staticmethod
    def get_rotation_matrix(angle_rad, axis):
        # Generate a unique key for the cache based on angle and axis
        cache_key = (round(angle_rad, 5), axis)  # Round to handle floating-point precision
        if cache_key in CubeVoxel.rotation_matrix_cache:
            # Move the used key to the end to show it's recently accessed
            CubeVoxel.rotation_matrix_cache.move_to_end(cache_key)
            return CubeVoxel.rotation_matrix_cache[cache_key]
        
        # Compute the rotation matrix if not in cache
        rotation_matrix = CubeVoxel._compute_rotation_matrix(angle_rad, axis)

        # Add the new rotation matrix to the cache
        CubeVoxel.rotation_matrix_cache[cache_key] = rotation_matrix
        # Ensure the cache does not exceed the specified limit
        if len(CubeVoxel.rotation_matrix_cache) > CubeVoxel.cache_limit:
            CubeVoxel.rotation_matrix_cache.popitem(last=False)  # Remove the oldest item
        return rotation_matrix

    @staticmethod
    def _compute_rotation_matrix(angle_rad, axis):
        # Your existing logic to compute the rotation matrix based on the axis
        if axis == 'x':
            return np.array([[1, 0, 0],
                             [0, np.cos(angle_rad), -np.sin(angle_rad)],
                             [0, np.sin(angle_rad), np.cos(angle_rad)]])
        elif axis == 'y':
            return np.array([[np.cos(angle_rad), 0, np.sin(angle_rad)],
                             [0, 1, 0],
                             [-np.sin(angle_rad), 0, np.cos(angle_rad)]])
        elif axis == 'z':
            return np.array([[np.cos(angle_rad), -np.sin(angle_rad), 0],
                             [np.sin(angle_rad), np.cos(angle_rad), 0],
                             [0, 0, 1]])
        else:
            raise ValueError("Axis must be 'x', 'y', or 'z'.")

    def transform_element(self, transformation_function, S):
        """
        Applies a transformation function to the entire voxel.

        Parameters:
        - transformation_function (callable): A function that takes and returns an np.array of node positions.
        """
        self.nodesR = transformation_function(self.nodesR, S)
        self.center = np.mean(self.nodesR, axis=0)

    def transform_node(self, transformation_function, S):
        """
        Applies a transformation function to each node of the voxel individually.

        Parameters:
        - transformation_function (callable): A function that takes and returns a single node position.
        """
        for i, node in enumerate(self.nodesR):
            self.nodesR[i] = transformation_function(node, S)
        self.center = np.mean(self.nodesR, axis=0)

    # Overload the arithmetic operators for easy linear matrix manipulation 

    def __mul__(self, other):
            """
            Supports multiplication of the CubeVoxel with a transformation matrix.

            Parameters:
            - other (np.array): A 3x3 or 4x4 transformation matrix for linear transformations.

            Returns:
            - CubeVoxel: A new instance of CubeVoxel with transformed nodes.
            """
            if not isinstance(other, np.ndarray) or other.shape not in [(3, 3), (4, 4)]:
                raise ValueError("Multiplication is supported only with a 3x3 or 4x4 transformation matrix.")

            # If a 3x3 matrix is provided, apply it directly.
            # If a 4x4 matrix is provided, assume it includes translation and apply it to homogeneous coordinates.
            if other.shape == (3, 3):
                transformed_nodes = np.dot(self.nodesR, other.T)
            else:  # other.shape == (4, 4)
                # Convert nodes to homogeneous coordinates, apply transformation, and convert back.
                nodes_homogeneous = np.hstack([self.nodesR, np.ones((self.nodesR.shape[0], 1))])
                transformed_nodes_homogeneous = np.dot(nodes_homogeneous, other.T)
                transformed_nodes = transformed_nodes_homogeneous[:, :3] / transformed_nodes_homogeneous[:, [3]]

            # Create a new CubeVoxel instance with the transformed nodes.
            return CubeVoxel(self.id, transformed_nodes, self.faceConnectNodes, self.colourF, self.colourE, self.faceAlpha)
    
    def _translate(self, vector, operation):
        """
        A helper method to support translation of the CubeVoxel by a vector.

        Parameters:
        - vector (np.array): A translation vector.
        - operation (callable): The numpy operation to apply (np.add or np.subtract).

        Returns:
        - CubeVoxel: A new instance of CubeVoxel with translated nodes.
        """
        if not isinstance(vector, (np.ndarray, list, tuple)) or len(vector) != 3:
            raise ValueError("Translation is supported only with a 3-element vector.")

        # Apply the translation operation to the nodes.
        translated_nodes = operation(self.nodesR, np.array(vector))

        # Create a new CubeVoxel instance with the translated nodes.
        return CubeVoxel(self.id, translated_nodes, self.faceConnectNodes, self.colourF, self.colourE, self.faceAlpha)

    def __add__(self, other):
        """
        Supports addition of the CubeVoxel with a translation vector.
        """
        return self._translate(other, np.add)

    def __sub__(self, other):
        """
        Supports subtraction of the CubeVoxel with a translation vector.
        """
        return self._translate(other, np.subtract)
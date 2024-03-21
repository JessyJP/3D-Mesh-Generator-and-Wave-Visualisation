import numpy as np
from copy import deepcopy

from cube_voxel import CubeVoxel

class PyramidVoxel(CubeVoxel):  # Inherit from CubeVoxel
    def __init__(self, id, colourF=[1.0, 1.0, 1.0], colourE=[0, 0, 0], faceAlpha=1.0):
        """
        Initializes a PyramidVoxel instance with specified properties or defaults.

        Parameters:
        - id (int): Identifier for the voxel.
        - colourF (list, optional): Face color of the voxel. Defaults to white.
        - colourE (list, optional): Edge color of the voxel. Defaults to black.
        - faceAlpha (float, optional): Transparency of the voxel faces. Defaults to 1.0 (opaque).
        """
        # Define nodes for a triangular pyramid
        nodesR = [[0, 0, 0], [1, 0, 0], [0.5, np.sqrt(3)/2, 0], [0.5, np.sqrt(3)/6, np.sqrt(6)/3]]
        # Define faces for a triangular pyramid (tetrahedron)
        faceConnectNodes = [[1, 2, 3], [1, 4, 2], [2, 4, 3], [3, 4, 1]]

        # Initialize the CubeVoxel (now PyramidVoxel) with pyramid-specific geometry
        super().__init__(id, nodesR, faceConnectNodes, colourF, colourE, faceAlpha)

    def copy(self):
        """
        Creates a deep copy of the PyramidVoxel instance.

        Returns:
        - PyramidVoxel: A new instance of PyramidVoxel with identical properties.
        """
        return deepcopy(self)

    # If there are any pyramid-specific methods or overrides needed, they can be added here.

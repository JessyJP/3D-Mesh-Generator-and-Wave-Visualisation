"""
MeshGenerator

Description:
    A Python class for generating and visualizing 2D/3D meshes. This class provides
    functionalities to set up mesh parameters, generate mesh elements based on these
    parameters, and visualize the generated mesh using Matplotlib.

Author:
    JessyJP (2022) - Original MATLAB version
    Conversion to Python - [Your Name or Identifier] (2024)

License:
    GPLv3 - This code is released under the GNU General Public License v3.0.
    See LICENSE.md for more details.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from cube_voxel import CubeVoxel

def handle_single_or_list(method):
    def wrapper(elements, *args, **kwargs):
        if isinstance(elements, list):
            return [method(element, *args, **kwargs) for element in elements]
        else:
            return method(elements, *args, **kwargs)
    return wrapper


class MeshGenerator:
    def __init__(self):
        """
        Constructor for the MeshGenerator class.
        Initializes mesh parameters, the cube template for mesh elements, plotting settings,
        and other internal properties.
        """
        # Private properties (indicated by underscore prefix)
        self._cubeTemplate = None  # Cube prototype
        self._elementAlignmentTypes = None  # Element alignment types
        self._elements = []  # Store the elements
        self._xyz = None  # The dimension and setup of the mesh

        # Public properties
        self.cube = None  # Working cube modified
        self.dimRemoveID = 0  # 0=>OFF, 1=>x, 2=>y, 3=>z
        self.Pplot = {'View': 3, 'gridON': True, 'scaleFactor': 1, 'axisEqual': True,
                      'limX': [0, 1], 'limY': [0, 1], 'limZ': [0, 1]} # Plot settings
        self.fig = None  # Figure for plotting

        # Initialize internal and public properties
        self._initialize_properties()

    def _initialize_properties(self):
        """
        Initializes both private and public properties of the class.
        """
        # XYZ properties initialization
        self._xyz = {'x': {'Label': 'x', 'ID': 1, 'Lim': [0, 1], 'N': 1, 'gapP': 0.1, 'ElementL': None, 'gap': None, 'Coords': None, 'mesh': None},
                    'y': {'Label': 'y', 'ID': 2, 'Lim': [0, 1], 'N': 1, 'gapP': 0.1, 'ElementL': None, 'gap': None, 'Coords': None, 'mesh': None},
                    'z': {'Label': 'z', 'ID': 3, 'Lim': [0, 1], 'N': 1, 'gapP': 0.1, 'ElementL': None, 'gap': None, 'Coords': None, 'mesh': None}}
        # % Constant IDs
        # % Lim:      Defautl Limits
        # % N:        Default Number of Elements
        # % gapP:     Element Gap factor between elements
        # % ElementL: Element side length
        # % gap:      Element side length gap
        # % Coords:   Mesh grid coordinates
        # % mesh:     

        for key, axis in self._xyz.items():
            axis['ElementL'] = self._calc_el_spacing(axis) * (1 - axis['gapP'])
            axis['gap'] = self._calc_el_spacing(axis) * axis['gapP']
            axis['Coords'] = []
            axis['mesh'] = None

        # Initialize cube template using CubeVoxel
        self._cubeTemplate = CubeVoxel(id=0)

        # Initialize element alignment types
        self._elementAlignmentTypes = ["Symmetric", "Center"]

        # Copy the cube template to the public cube property
        self.cube = self._cubeTemplate.copy()

    def _calc_el_spacing(self, axis):
        """
        Calculate the spacing between elements along a given axis.

        Parameters:
        - axis: A dictionary containing the mesh properties for a specific axis.

        Returns:
        - The calculated element spacing as a float.
        """
        return (axis['Lim'][1] - axis['Lim'][0]) / axis['N']

    def setup_mesh(self, xyz_=None):
        """
        Sets up or updates the mesh based on provided or existing xyz parameters.
        
        Parameters:
        - xyz_: Optional dictionary containing updated mesh parameters for x, y, and/or z axes.
        """
        # Update the mesh parameters if a new xyz_ dictionary is provided
        if xyz_ is not None:
            self._xyz = xyz_.copy()
            print('Table updated!')
        else:
            print('Table not updated!')

        # Update cube template colors and face alpha from the cube
        self._cubeTemplate.colourF = self.cube.colourF
        self._cubeTemplate.colourE = self.cube.colourE
        self._cubeTemplate.faceAlpha = self.cube.faceAlpha
        # Apply the updated template to the cube
        self.cube = self._cubeTemplate.copy()

        # Update viewing limits based on xyz settings and scaleFactor % Default viewing Limitts 
        for axis_key, axis_props in self._xyz.items():
            self.Pplot[f'lim{axis_key.upper()}'] = np.array(axis_props['Lim']) * self.Pplot['scaleFactor']

        # Generate coordinates for each axis
        for axis_key, axis_props in self._xyz.items():
            self._xyz[axis_key] = self.setup_coordinate(self._xyz[axis_key])
            # axis_props['Coords'] = np.linspace(axis_props['Lim'][0], axis_props['Lim'][1], axis_props['N'] + 1)[:-1]

        # Handle dimension removal and adjust cube nodes if necessary
        if self.dimRemoveID > 0:
            # Adjust the cube's nodes based on the removed dimension
            for axis_key, axis_props in self._xyz.items():
                if axis_props['ID'] == self.dimRemoveID:
                    # This effectively removes the dimension by setting all coordinates in that dimension to 0
                    self.cube['nodesR'][:, self.dimRemoveID - 1] = 0

            # Adjust the FaceConnectNodes based on the removed dimension
            if self.dimRemoveID == 1:
                self.cube['FaceConnectNodes'] = np.array([[1, 2, 4, 3]]) - 1
            elif self.dimRemoveID == 2:
                self.cube['FaceConnectNodes'] = np.array([[1, 2, 3, 4]]) - 1
            elif self.dimRemoveID == 3:
                self.cube['FaceConnectNodes'] = np.array([[1, 2, 4, 3]]) - 1

        # Scale the cube based on the element spacing and gap percentage
        scale_factors = np.array([self._calc_el_spacing(self._xyz['x']) * (1 - self._xyz['x']['gapP']),
                                  self._calc_el_spacing(self._xyz['y']) * (1 - self._xyz['y']['gapP']),
                                  self._calc_el_spacing(self._xyz['z']) * (1 - self._xyz['z']['gapP'])])
        self.cube.scale( scale_factors)

        # Prepare elements for plotting
        X, Y, Z = np.meshgrid(self._xyz['x']['Coords'], self._xyz['y']['Coords'], self._xyz['z']['Coords'], indexing='ij')
        # Assign the generated mesh grid back to the _xyz attribute
        self._xyz['x']['mesh'] = X
        self._xyz['y']['mesh'] = Y
        self._xyz['z']['mesh'] = Z
        for i, (x, y, z) in enumerate(zip(X.flatten(), Y.flatten(), Z.flatten())):
            el = self.cube.copy()
            el.id = i
            # Translate the element based on the generated coordinates
            el.translate(np.array([x, y, z]));
            el.faceConnectNodes_Fast = el.faceConnectNodes + (i)* (8*(self.dimRemoveID==0)+4*(self.dimRemoveID!=0))
            self._elements.append(el)

    def setup_coordinate(self, axis):
        """
        Set up the coordinate array for a given axis based on its mesh properties.

        Parameters:
        - axis: A dictionary containing the mesh properties for a specific axis.

        Returns:
        - The axis dictionary updated with calculated element length, gap, and coordinates.
        """
        # Calculate element spacing (length) and gap based on the axis properties
        element_spacing = self._calc_el_spacing(axis) * (1 - axis['gapP'])
        gap = self._calc_el_spacing(axis) * axis['gapP']

        # Update the axis dictionary with calculated values
        axis['ElementL'] = element_spacing
        axis['gap'] = gap

        # Adjust the limits if the dimension corresponds to the one being removed
        if self.dimRemoveID == axis['ID']:
            axis['Lim'] = [0, 0]  # This effectively removes the dimension

        # Generate the coordinates for the axis
        coords = np.linspace(axis['Lim'][0], axis['Lim'][1], axis['N'] + 1)
        coords = coords[:-1]  # Exclude the last point to match MATLAB behavior

        # Update the axis dictionary with the generated coordinates
        axis['Coords'] = coords

        return axis


    def _adjust_cube_for_removed_dimension(self):
        """
        Adjusts the cube's nodes and face connections if a dimension is removed.
        """
        if self.dimRemoveID == 1:
            self.cube['FaceConnectNodes'] = np.array([[1, 2, 4, 3]]) - 1
        elif self.dimRemoveID == 2:
            self.cube['FaceConnectNodes'] = np.array([[1, 2, 3, 4]]) - 1
        elif self.dimRemoveID == 3:
            self.cube['FaceConnectNodes'] = np.array([[1, 2, 4, 3]]) - 1
        # Note: Additional logic may be needed to adjust nodesR based on the removed dimension.

    def _generate_mesh_elements(self):
        """
        Generates mesh elements based on the current settings of the mesh.
        """
        elements = []
        X, Y, Z = np.meshgrid(self._xyz['x']['Coords'], self._xyz['y']['Coords'], self._xyz['z']['Coords'], indexing='ij')
        for i, (x, y, z) in enumerate(zip(X.flatten(), Y.flatten(), Z.flatten())):
            el = self.cube.copy()
            el['id'] = i
            # Translate the element based on the generated coordinates
            el['nodesR'] += np.array([x, y, z])
            elements.append(el)
        return elements

    def plot_elements(self, elements):
        """
        Visualize the generated mesh elements using Matplotlib in a 3D plot.

        Parameters:
        - elements: Optional list of elements to plot. If None, plots all generated elements.
        """

        # Initialize the figure and axes if they don't exist or if the figure was closed
        if not hasattr(self, 'fig') or self.fig is None or not plt.fignum_exists(self.fig.number):
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111, projection='3d')
            # Set the window style if needed; note that 'docked' style might not directly translate to matplotlib
        else:
            self.ax.clear()

        # Set plot limits and labels based on MeshGenerator properties
        self.ax.set_xlim(self.Pplot['limX'])
        self.ax.set_ylim(self.Pplot['limY'])
        self.ax.set_zlim(self.Pplot['limZ'])
        
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')

        # Enable or disable the grid based on MeshGenerator properties
        self.ax.grid(self.Pplot.get('gridON', True))

        # Ensure axis aspect ratio is equal if specified
        if self.Pplot.get('axisEqual', True):
            self._set_equal_aspect()

        if not isinstance(elements, list):
            elements = [elements]
        # Plot each element
        for el in elements:
            # Assuming 'FaceConnectNodes' are indices starting from 0 in Python
            verts = [el.nodesR[face] for face in el.faceConnectNodes]
            face_collection = Poly3DCollection(verts, facecolors=el.colourF, edgecolors=el.colourE, alpha=el.faceAlpha)
            self.ax.add_collection3d(face_collection)

        plt.draw()
        plt.show()

    def plot_mesh_fast(self, elements):
        """
        Visualizes the mesh elements efficiently using a single Poly3DCollection.
        """
        if not hasattr(self, 'fig') or not plt.fignum_exists(self.fig.number):
            self.fig = plt.figure()
            self.ax = self.fig.add_subplot(111, projection='3d')
        else:
            self.ax.clear()

        # Prepare a list for collecting vertices and faces across all elements
        all_faces = []
        all_vertices = np.vstack([el['nodesR'] for el in elements])

        # Offset for vertices indexing in faces due to concatenation
        offset = 0
        for el in elements:
            # Adjust face indices based on current offset
            faces = el['FaceConnectNodes'] + offset
            all_faces.extend(faces)
            offset += len(el['nodesR'])

        # Create a Poly3DCollection from all faces and vertices
        face_collection = Poly3DCollection(all_vertices[all_faces], 
                                           facecolors=self.cube['colourF'], 
                                           edgecolors=self.cube['colourE'], 
                                           alpha=self.cube['FaceAlpha'])

        self.ax.add_collection3d(face_collection)

        # Set plot limits and labels based on MeshGenerator properties
        self.ax.set_xlim(self.Pplot['limX'])
        self.ax.set_ylim(self.Pplot['limY'])
        self.ax.set_zlim(self.Pplot['limZ'])
        
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')

        # Enable or disable the grid based on MeshGenerator properties
        self.ax.grid(self.Pplot.get('gridON', True))

        # Ensure axis aspect ratio is equal if specified
        if self.Pplot.get('axisEqual', True):
            self._set_equal_aspect()

        plt.show()

    def _set_equal_aspect(self):
        """
        A workaround to set equal aspect ratio for the 3D plot, as Matplotlib does not support this directly.
        """
        extents = np.array([getattr(self.ax, f'get_{axis}lim')() for axis in 'xyz'])
        sz = extents[:,1] - extents[:,0]
        centers = np.mean(extents, axis=1)
        maxsize = max(abs(sz))
        r = maxsize/2
        for ctr, val in zip(centers, sz):
            self.ax.auto_scale_xyz(*[ctr - r, ctr + r], had_data=True)

    def update_colors_and_plot_properties(self, colorF, colorE, face_alpha, plot_scale_factor):
        """
        Updates the colors and plot properties of the mesh.

        Parameters:
        - colorF: Face color of the cube.
        - colorE: Edge color of the cube.
        - face_alpha: Transparency of the cube faces.
        - plot_scale_factor: Scale factor for the plot.
        """
        self.cube.colourF = colorF
        self.cube.colourE = colorE
        self.cube.faceAlpha = face_alpha
        self.Pplot['scaleFactor'] = plot_scale_factor

    def update_xyz_settings(self, xyz_settings):
        """
        Updates the xyz settings of the mesh based on the provided dictionary.

        Parameters:
        - xyz_settings: A dictionary containing the settings for x, y, and z axes.
        """
        for axis, settings in xyz_settings.items():
            if axis in self._xyz:
                for setting, value in settings.items():
                    if setting == 'Lim' and isinstance(value, list) and len(value) == 2:
                        self._xyz[axis]['Lim'] = value
                    elif setting in self._xyz[axis]:
                        self._xyz[axis][setting] = value
                # Recalculate element spacing, gap, and coordinates based on updated settings
                self._xyz[axis]['ElementL'] = self._calc_el_spacing(self._xyz[axis]) * (1 - self._xyz[axis]['gapP'])
                self._xyz[axis]['gap'] = self._calc_el_spacing(self._xyz[axis]) * self._xyz[axis]['gapP']
                self._xyz[axis]['Coords'] = self.setup_coordinate(self._xyz[axis])['Coords']
            else:
                print(f"Axis '{axis}' not found in mesh settings.")

    ## Modification methods
            
    @staticmethod
    @handle_single_or_list
    def translate(element, translation_vector):
        element.translate(translation_vector)
        return element

    @staticmethod
    @handle_single_or_list
    def scale(element, scale_factors):
        element.scale(scale_factors)
        return element

    @staticmethod
    @handle_single_or_list
    def rotate(element, angle, axis='z', center_of_rotation=None):
        element.rotate(angle, axis, center_of_rotation)
        return element

    @staticmethod
    @handle_single_or_list
    def transform_elements(element, transformation_function, S):
        element.transform_element(transformation_function, S)
        return element

    @staticmethod
    @handle_single_or_list
    def transform_nodes(element, transformation_function):
        element.transform_node(transformation_function)
        return element

    # Example usage of operator overloads
    def __mul__(self, matrix):
        for element in self.elements:
            element * matrix  # Assuming __mul__ is defined in CubeVoxel for matrix multiplication

    def __add__(self, vector):
        for element in self.elements:
            element + vector  # Assuming __add__ is defined in CubeVoxel for translation

    def __sub__(self, vector):
        for element in self.elements:
            element - vector  # Assuming __sub__ is defined in CubeVoxel for subtraction

    @property
    def xyz(self):
        """
        Property getter for the mesh dimensions.
        """
        return self._xyz

    def get_elements(self):
        return self._elements.copy()

    def get_matrix(self, elements, property):
        """
        Converts a specified property from a list of elements into a concatenated matrix.
        
        Parameters:
        - elements: List of elements to process.
        - property: The property of the elements to convert into a matrix.
        
        Returns:
        - A numpy array containing the concatenated values of the specified property from all elements.
        """
        return np.concatenate([el[property] for el in elements])

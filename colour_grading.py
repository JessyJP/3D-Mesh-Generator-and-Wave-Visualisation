import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def vals2colormap(vals, colormap='jet', crange=None):
    """
    Map values to a colormap.

    Parameters:
    - vals: A vector of values to map to a colormap.
    - colormap: A string name of a matplotlib colormap.
    - crange: The range of values to map to the minimum and maximum of the colormap.

    Returns:
    - rgb: An array of RGB colors corresponding to the vals mapped onto the colormap.
    """
    # Ensure vals is an array to handle both individual and sequences of values
    if np.isscalar(vals):
        vals = np.array([vals])
    elif isinstance(vals, tuple):
        vals = np.array(vals)

    if crange is None:
        crange = [np.min(vals), np.max(vals)]
    
    norm = mcolors.Normalize(vmin=crange[0], vmax=crange[1], clip=True)
    mapper = plt.cm.ScalarMappable(norm=norm, cmap=colormap)
    
    rgb = mapper.to_rgba(vals, bytes=True)[:, :3]  # Discard alpha channel
    return rgb

def displace_center_uxyz(E, S):
    """
    Displace the center of an element based on a transformation function.
    """
    rc = E.center
    T = S['Uxyz'](rc, S) - rc
    E.nodesR += T * np.ones((E.nodesR.shape[0], 1))

def colour_grade(E, Es_original, distM, cmapInd):
    """
    Apply colour grading to an element based on its displacement magnitude.
    """

    E0 = Es_original[E.id]  # Adjusted to use dot notation if E is a class instance
    R0 = E0.nodesR  # Assuming nodesR is the correct attribute name
    R = E.nodesR
    p = np.mean(np.sqrt(np.sum((R - R0)**2, axis=1))) / distM
    
    cmaps = ['parula', 'turbo', 'hsv', 'hot', 'cool', 'spring', 'summer',
             'autumn', 'winter', 'gray', 'bone', 'copper', 'pink', 'lines', 'jet',
             'colorcube', 'prism', 'flag']
    
    if 0 < cmapInd < len(cmaps):
        E.colourF = vals2colormap(p, colormap=cmaps[cmapInd], crange=[0, 1])
    else:
        E.colourF = np.array([238, 104, 104]) / 255 * p + np.array([204, 204, 0]) / 255 * (1 - p)

def colour_update(E, cube):
    """
    Update the colours of an element based on the cube's colours.
    """
    E.colourF = cube.colourF
    E.colourE = cube.colourE
    E.faceAlpha = cube.faceAlpha

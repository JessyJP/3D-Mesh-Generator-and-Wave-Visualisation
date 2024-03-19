import numpy as np

def love_wave_displacements(r_in, S):
    """
    Calculate displacements for Love Waves.
    """
    x, y, z = r_in[:, 0], r_in[:, 1], -r_in[:, 2]
    c, c1, c2, A1, A2, h, omega, t = S.values()

    Uy = np.zeros_like(y)
    mask_z_less_h = z < h
    mask_z_greater_equal_h = z >= h

    if np.any(mask_z_less_h):
        Uy[mask_z_less_h] = A1 * np.cos(omega * z[mask_z_less_h] * np.sqrt(1/c1**2 - 1/c**2)) * np.sin(omega * (t - x/c))
    if np.any(mask_z_greater_equal_h):
        Uy[mask_z_greater_equal_h] = A2 * np.exp(-omega * z[mask_z_greater_equal_h] * np.sqrt(1/c**2 - 1/c2**2)) * np.sin(omega * (t - x/c))

    R_out = np.vstack((x, y + Uy, -z)).T
    return R_out

def p_longitudinal_wave_displacements(r_in, S):
    """
    Calculate displacements for P Longitudinal Waves.
    """
    x, y, z = r_in.T
    t, A, k, omega = S.values()

    Ux = A * np.exp(1j * (omega * t - k * x))
    R_out = np.vstack((x + np.real(Ux), y, z)).T
    return R_out

def rayleigh_wave_displacements(r_in, S):
    """
    Calculate displacements for Rayleigh Waves.
    """
    x, y, z = r_in.T
    func, c, cl, ct, A, omega, t = S.values()

    s1 = -np.sqrt(c**2 - cl**2)
    s2 = -np.sqrt(c**2 - ct**2)
    Ux = A * c * (np.exp(-s1 * z) - (2 * s2 * s1 / ct**2) * np.exp(-s2 * z)) * np.exp(1j * (c * x - omega * t - np.pi/2))
    Uz = A * s1 * (np.exp(-s1 * z) - (2 * c**2 / (2 * c**2 - ct**2)) * np.exp(-s2 * z)) * np.exp(1j * (c * x - omega * t))

    R_out = np.vstack((np.real(Ux) + x, y, np.real(Uz) + z)).T
    return R_out

def s_transverse_wave_displacements(r_in, S):
    """
    Calculate displacements for S Transverse Waves.
    """
    x, y, z = r_in.T
    t, A, omega, k = S.values()

    Uy = A * np.cos(omega * t - k * z)
    Uz = A * np.sin(omega * t - k * x)

    R_out = np.vstack((x, y + Uy * 0, z + Uz)).T
    return R_out

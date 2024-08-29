import MDAnalysis as mda
import numpy as np
from scipy.spatial import ConvexHull

def load_xyz(file_path):
    """
    Load the XYZ file using MDAnalysis.
    """
    u = mda.Universe(file_path)
    return u

def calculate_surface_area_and_volume(u):
    """
    Calculate the surface area and volume of the molecule using ConvexHull.
    """
    positions = u.atoms.positions
    hull = ConvexHull(positions)

    # Surface area
    surface_area = hull.area

    # Volume
    volume = hull.volume

    return surface_area, volume

def calculate_pore_size(u, probe_radius=1.7):
    """
    Calculate the pore size distribution of the molecule.
    """
    positions = u.atoms.positions

    # Calculate pairwise distances
    dist_matrix = np.sqrt(np.sum((positions[:, np.newaxis, :] - positions[np.newaxis, :, :]) ** 2, axis=-1))

    # Determine pore sizes by finding distances larger than a given threshold (e.g., probe radius)
    pore_sizes = dist_matrix[dist_matrix > 2 * probe_radius]

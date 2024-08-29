# Molecular-Surface-and-Pore-Size-Analysis-using-MDAnalysis-and-SciPy
This repository contains a Python script designed to calculate the surface area, volume, and pore size distribution of molecules using atomic positions from XYZ files. The script leverages MDAnalysis for molecular data processing and SciPy's ConvexHull for geometric computations.
## Features
Surface Area and Volume Calculation: Computes the surface area and volume of a molecule by generating a convex hull around the atomic positions.
Pore Size Distribution: Calculates the distribution of pore sizes by analyzing pairwise distances between atoms and comparing them to a defined probe radius.
Customizable Probe Radius: Allows for adjusting the probe radius to match different molecular structures.
## Requirements
*Python 3.x
*MDAnalysis
*NumPy
*SciPy
## Example
python molecular_analysis.py
The script will print the calculated surface area, volume, and a list of pore sizes.
## Usage
1.Place your XYZ file in the working directory.
2.Update the file path in the script to point to your XYZ file.
3.Run the script to output the surface area, volume, and pore sizes.

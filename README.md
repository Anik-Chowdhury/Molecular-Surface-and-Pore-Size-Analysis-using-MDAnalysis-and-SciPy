# Molecular Surface and Pore Size Analysis using MDAnalysis and SciPy

This repository contains a Python script designed to calculate the surface area, volume, and pore size distribution of molecules using atomic positions from XYZ files. The script leverages `MDAnalysis` for molecular data processing and `SciPy`'s `ConvexHull` for geometric computations.

## Features
- **Surface Area and Volume Calculation**: Computes the surface area and volume of a molecule by generating a convex hull around the atomic positions.
- **Pore Size Distribution**: Calculates the distribution of pore sizes by analyzing pairwise distances between atoms and comparing them to a defined probe radius.
- **Customizable Probe Radius**: Allows for adjusting the probe radius to match different molecular structures.

## Requirements
- Python 3.x
- MDAnalysis
- NumPy
- SciPy
## Installation
To install the required dependencies, use the following command:

```bash
pip install MDAnalysis numpy scipy

```

## Usage
1. Place your XYZ file in the working directory.
2. Update the file path in the script to point to your XYZ file.
3. Run the script to output the surface area, volume, and pore sizes.

### Example
To run the script, use the following command:

```bash
python molecular_analysis.py

```
The script will output:

- Surface Area in Å²
- Volume in Å³
- A list of pore sizes larger than the specified probe radius

## Applications

This tool can be used for:
- Analyzing molecular structures in computational chemistry.
- Studying material properties such as porosity and surface area.
- Supporting research in nanotechnology, catalysis, and materials science.

## License

This project is licensed under the MIT License.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request if you have improvements or bug fixes.

## Acknowledgments

- [MDAnalysis](https://www.mdanalysis.org/) for molecular analysis.
- [SciPy](https://www.scipy.org/) for the computational tools.



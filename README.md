

# VBData_Gen: A method for generating single-bond valence bond property data of organic small molecules

This repository contains helper functions used in the creation and analysis of the valence bond property database


 

### Organization of the repository is as follows:
- `confsgen.py`: a python script that takes a SMILES string as input and generates multiple 3D conformations using RDKit.
- `mndogen.py`: a python script that generates MNDO input files from existing molecular coordinate XYZ files.
- `gasgen.sh`: a shell script that reads MNDO optimization results to generate Gaussian input files and the `checkgas.py` script used to check Gaussian optimization results.
- `xmvbgen.py`: a python script that reads Gaussian optimization results and generates input files for Xiamen valence bond.
- `checkxmo.py`: a python script that reads the valence properties from an xmvb output file and organizes them into a csv file.

`data/`:
- `res.csv`: a dataset of unique single-bond valence bond properties of small organic molecules.









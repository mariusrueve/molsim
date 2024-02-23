# molsim

Molecular similarity (molsim). A command line tool to find the most similar molecules in a database to molecules in an input file.
The similarity is calculated using the Tanimoto coefficient using Morgan fingerprints. The input file can be in SDF or SMILES format. The database can be in SDF or SMILES format.
The output is written to a CSV file with the following columns: `input_molecule, database_molecule, similarity_score`.

## Installation

```
conda create -p .conda
conda activate .conda/
conda install -c conda-forge cmake rdkit eigen boost spdlog gxx

mkdir build
cd build
cmake ..
make molsim
```

## Usage

### Help
```
./molsim --help
```

### Example
```
./molsim --input path/to/input.[sdf/smi] --database path/to/database.[sdf/smi]
```
The output will be written to best_matches.csv in the same directory where program is executed.

# PyProt
Automation scripts for protein model preparation and fetching data from UniProt.

## Project setup
First create the python package by running the command:
```
python -m build
```

Then move to `dist` subfolder and run:
```
pip install pyprot-0.0.1-py3-none-any.whl
```

Package requirements:
- bipython
- modeller
- pandas

## Usage
### PyProtein
Contains the tools for interfacing with UniProt. 
```
from Pyprot import PyProtein
```

`PyProtein.search(<text>)` performs a search on UniProt and returns the results in a specified format.

`PyProtein.seq_download(<name>)` takes a UniProt sequence name and returns the sequence as a biopython Seq object.

### Model
'Model' is the subpackage dedicated to construction of structure models. It is esentially the extrapolation of 'modeller' by Salilab with some functionality added in.

First, and instance of the `Input` class has to be created, as this instance acts as the holder of the settings. When invoking the class, we need to provide the path to the template structure, a list of tuples referencing pairs of template-model subunits (for the modelling reference), and optionally the names of the input and target structures and path to the pir file containing the alignment of at least the sequences involved in the modelling we wish to perform (template and model sequences!). 

Next, the `Input().chain_preparation()` is run to prepare the sequences required for the modelling.

After preparing the sequences `Input().ali_write()` is called for writing the sequence alignment, followed by `Input().pdb_clean()` for preparation of the PDB input file with the template structure.

The initial segment would thus look something like that:

```
settings = Model.Input(<input-structure-path>,
                       [("Alpha-1", "Alpha-6"), ("Beta-3", "Beta-3")],
                       <input-structure-name>,
                       <target-structure-name>,
                       <path-to-pir-file>)
settings.chain_preparation()
settings.ali_write()
settings.pdb_clean()
```

The settings get passed to the instance of class `Model` when it is created. After we create it we run the modelling with `Model().run_model()` and optionally clean up the created models with `Model().chain_cleanup()`. The last step opens individual created models, renames chains and renumbers residues to match the template's and finally saves the corrected structures.
In the end, the modelling code looks something like this:

```
model = Model.Model(settings)
model.run_model()
model.chain_cleanup()
```
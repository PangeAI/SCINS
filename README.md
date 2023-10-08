# SCINS

## Intended usage

Clone and install the package:

`git clone git@github.com:PangeAI/SCINS.git`

`cd SCINS/scins`

`pip install -e scins`

In your script:

```python
from rdkit import Chem
from rdkit.Chem.Scaffold.MurckoScaffold import GetScaffoldForMol, MakeScaffoldGeneric
from scins import scins

mol = Chem.MolFromSmiles('Cc1cc(C)nc(SCC(=O)Nc2ncc(Cc3ccccc3)s2)n1')
scaffold = GetScaffoldForMol(mol)
generic_scaffold = MakeScaffoldGeneric(scaffold)
scins_str = scins.mol_to_scins(mol)
```

## Important Note:

It is essential that you apply the function on
the generic scaffold. Otherwise, the result would not be meaningful.

## Second Important Note:

A slightly weird behaviour was observed with the code above even though
it is does implement what was described in the original paper. 
In the following example, notice how the carbonyl oxygen lead to a "separate chain". 
This does not correspond to my chemical intuition.

![alt text](assets/weird_scaffold_def.png "Figure 1")

The following code snippet shows how to avoid this issue:

```python
scaffold = scins.GetScaffoldForMol_edited(mol)
generic_scaffold = MakeScaffoldGeneric(scaffold)
scins_str = scins.mol_to_scins(mol)
```
For the same molecule, the function GetScaffoldForMol_edited in the 
package trims the carbonyl oxygen (and other bits that are sticking out).
Therefore, this definition of the scaffold harmonates with my notions of the world.

![alt text](assets/better_scaffold_def.png "Figure 2")

However, this function has not been tested thoroughly, so apply at your own risk.

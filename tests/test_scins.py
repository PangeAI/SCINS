from collections import Counter
import json
from pathlib import Path
import unittest

from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol, MakeScaffoldGeneric

from scins.scins import (
    mol_to_num_chain_assemblies,
    mol_to_chain_lengths,
    _mol_to_rings,
    get_ring_assemblies,
    GetScaffoldForMol_edited,
    generic_scaffold_mol_to_scins,
)

parent_dir = Path(__file__).parent
with open(parent_dir / "tests_data/test_data.json") as f:
    test_data = json.load(f)


class TestSCINS(unittest.TestCase):
    def test_mol_to_num_chain_assemblies(self):
        smi2num_chain_assemblies = test_data["num_chain_assemblies"]
        for smi, answer in smi2num_chain_assemblies.items():
            mol = Chem.MolFromSmiles(smi)
            gen_scaffold = MakeScaffoldGeneric(GetScaffoldForMol(mol))
            num_chain_assemblies = mol_to_num_chain_assemblies(gen_scaffold)
            self.assertEqual(num_chain_assemblies, answer)

    def test_mol_to_chain_lengths(self):
        smi2chain_assemblies_lengths = test_data["chain_assemblies_lengths"]
        for smi, answer in smi2chain_assemblies_lengths.items():
            mol = Chem.MolFromSmiles(smi)
            gen_scaffold = GetScaffoldForMol(MakeScaffoldGeneric(mol))
            chain_lengths = mol_to_chain_lengths(gen_scaffold)
            counter_chain_lengths = Counter(chain_lengths)
            counter_answer = Counter(answer)
            self.assertEqual(counter_chain_lengths, counter_answer)

    def test_get_num_ring_assemblies(self):
        smi2num_ring_assemblies = test_data["num_ring_assemblies"]
        for smi, answer in smi2num_ring_assemblies.items():
            mol = Chem.MolFromSmiles(smi)
            gen_scaffold = GetScaffoldForMol(MakeScaffoldGeneric(mol))
            rings = _mol_to_rings(gen_scaffold)
            ring_assemblies = get_ring_assemblies(rings)
            self.assertEqual(len(ring_assemblies), answer)

    def test_generic_scaffold_to_scins(self):
        smiles_to_scins = test_data["scins"]
        for smi, answer in smiles_to_scins.items():
            mol = GetScaffoldForMol(MakeScaffoldGeneric(Chem.MolFromSmiles(smi)))
            res = generic_scaffold_mol_to_scins(mol)
            self.assertEqual(res, answer)

    def test_GetScaffoldForMol_edited(self):
        from rdkit import RDLogger

        rdl = RDLogger.logger()
        # 🤫
        rdl.setLevel(RDLogger.CRITICAL)
        smiles_to_scins = test_data["scins"]
        for smi, _ in smiles_to_scins.items():
            gen_scaff1 = Chem.MolToSmiles(
                GetScaffoldForMol(MakeScaffoldGeneric(Chem.MolFromSmiles(smi))),
                canonical=True,
            )
            gen_scaff2 = Chem.MolToSmiles(
                MakeScaffoldGeneric(GetScaffoldForMol_edited(Chem.MolFromSmiles(smi))),
                canonical=True,
            )
            self.assertEqual(gen_scaff1, gen_scaff2)

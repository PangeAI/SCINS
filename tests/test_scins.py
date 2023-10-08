import unittest

from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol, MakeScaffoldGeneric

from scins.scins import (mol_to_num_chain_assemblies,
                         mol_to_chain_lengths,
                         get_rings_for_mol,
                         get_ring_assemblies,
                         GetScaffoldForMol_edited
                         )


class TestSCINS(unittest.TestCase):

    def test_mol_to_num_chain_assemblies(self):
        smi2num_chain_assemblies = {'Cc1cc(C)nc(SCC(=O)Nc2ncc(Cc3ccccc3)s2)n1': 2}
        for smi, answer in smi2num_chain_assemblies.items():
            mol = Chem.MolFromSmiles(smi)
            gen_scaffold = MakeScaffoldGeneric(GetScaffoldForMol(mol))
            num_chain_assemblies = mol_to_num_chain_assemblies(gen_scaffold)
            print(num_chain_assemblies)
            self.assertEqual(num_chain_assemblies, answer)

    def test_mol_to_chain_lengths(self):
        smi2num_chain_assemblies = {'Cc1cc(C)nc(SCC(=O)Nc2ncc(Cc3ccccc3)s2)n1': {3, 1, 2, 2}}
        for smi, answer in smi2num_chain_assemblies.items():
            mol = Chem.MolFromSmiles(smi)
            gen_scaffold = MakeScaffoldGeneric(GetScaffoldForMol(mol))
            chain_lengths = set(mol_to_chain_lengths(gen_scaffold))
            difference = chain_lengths ^ answer
            self.assertEqual(len(difference), 0)

    def test_mol_to_chain_lengths_with_trimmed_scaffold(self):
        smi2num_chain_assemblies = {'Cc1cc(C)nc(SCC(=O)Nc2ncc(Cc3ccccc3)s2)n1': {5, 2}}
        for smi, answer in smi2num_chain_assemblies.items():
            mol = Chem.MolFromSmiles(smi)
            gen_scaffold = MakeScaffoldGeneric(GetScaffoldForMol_edited(mol))
            chain_lengths = set(mol_to_chain_lengths(gen_scaffold))
            print(chain_lengths)
            difference = chain_lengths ^ answer
            self.assertEqual(len(difference), 0)

    def test_get_num_ring_assemblies(self):
        smi2num_ring_assemblies = {'Cc1cc(C)nc(SCC(=O)Nc2ncc(Cc3ccccc3)s2)n1': 3}
        for smi, answer in smi2num_ring_assemblies.items():
            mol = Chem.MolFromSmiles(smi)
            gen_scaffold = MakeScaffoldGeneric(GetScaffoldForMol(mol))
            rings = get_rings_for_mol(gen_scaffold)
            ring_assemblies = get_ring_assemblies(rings)
            self.assertEqual(len(ring_assemblies), answer)

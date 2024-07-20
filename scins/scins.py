#
#  Copyright (c) 2024 Kamen Petrov
#  All rights reserved.
#
#  This file is part of the SCINS project.
#  The contents are covered by the terms of the MIT license
#  which is included in the file LICENSE, found at the root
#  of the source tree.
#  The orginal paper that first introduced SCINS
#  by Schuffenhauer et al.: https://doi.org/10.1021/ci6004004

import logging
from copy import deepcopy
from collections import Counter

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

ERROR_SCINS = "ERROR_SCINS"


def _rdkit_mol_warning(mol):
    logging.warning(
        "Input should be an RDKit molecule object, but got an object of type %s. Therefore returning ERROR_SCINS"
        % type(mol)
    )


def _update_mol_graph_dict_inplace(
    atom_key: int, atom_bound: int, non_ring_mol_graph: dict[int, list[int]]
) -> None:
    if atom_key in non_ring_mol_graph:
        non_ring_mol_graph[atom_key].append(atom_bound)
    else:
        non_ring_mol_graph[atom_key] = [atom_bound]


def _mol_to_non_ring_mol_graph(mol: Chem.rdchem.Mol) -> dict[int, list[int]]:
    """
    Obtains the non-ring molecular graph (i.e. chain only atom graph) from
    an rdkit mol object (intended to be the generic scaffold of a molecule).
    :param mol: Chem.rdchem.Mol - The input generic scaffold molecule
    :return: Dict[int, List[int]] - The non-ring molecular graph
    """
    # Create a dictionary to track the atoms connected by bonds
    non_ring_mol_graph = {}

    # Iterate through all bonds in the molecule
    for bond in mol.GetBonds():
        if not bond.IsInRing():
            atom1_idx = bond.GetBeginAtomIdx()
            atom2_idx = bond.GetEndAtomIdx()

            # Check if both atoms are already in the connected_atoms dictionary's keys
            _update_mol_graph_dict_inplace(atom1_idx, atom2_idx, non_ring_mol_graph)
            _update_mol_graph_dict_inplace(atom2_idx, atom1_idx, non_ring_mol_graph)
    return non_ring_mol_graph


def _non_ring_mol_graph_to_num_chain_assemblies(
    non_ring_mol_graph: dict[int, list[int]]
) -> int:
    """
    Finds the number of chain assemblies.
    :param non_ring_mol_graph: Dict[int, List[int]] - The non-ring molecular graph for a generic scaffold.
    :return: num_chain_assemblies: int - The number of chain assemblies.
    """
    visited = set()
    num_chain_assemblies = 0

    # Depth-First Search (DFS) to find the number of chains
    for atom in non_ring_mol_graph:
        if atom not in visited:
            num_chain_assemblies += 1
            stack = [atom]
            while stack:
                node = stack.pop()
                visited.add(node)

                for neighbor in non_ring_mol_graph[node]:
                    if neighbor not in visited:
                        stack.append(neighbor)
    return num_chain_assemblies


def _non_ring_mol_graph_to_chain_lengths(
    non_ring_mol_graph: Dict[int, List[int]]
) -> List[int]:
    visited = set()
    # chains = 0
    chain_lengths = []

    # Depth-First Search (DFS) to find unbranched chains
    for atom in non_ring_mol_graph:
        # if atom not in ring_atoms: # want to start the chain with a ring atom
        #     continue
        if atom not in visited:
            # chains += 1
            chain_length = 0
            stack = [atom]
            while stack:
                node = stack.pop()
                visited.add(node)

                # Check if the current node has more than one neighbor
                neighbors = non_ring_mol_graph[node]
                if len(neighbors) > 2:
                    continue

                for neighbor in non_ring_mol_graph[node]:
                    # this looks inefficient, but is here to catch edge cases
                    if len([n for n in non_ring_mol_graph[neighbor]]) >= 3:
                        chain_length += 1
                        if neighbor not in visited:
                            stack.append(neighbor)
                    elif neighbor not in visited:
                        stack.append(neighbor)
                        chain_length += 1

            chain_lengths.append(chain_length)

    return chain_lengths


def _mol_to_rings(mol: Chem.rdchem.Mol) -> list[list[int]]:
    ring_info = mol.GetRingInfo()
    rings = [list(ring) for ring in ring_info.AtomRings()]
    return rings


def get_ring_assemblies(rings: list[list[int]]) -> list[list[int]]:
    """
    Merge rings with common atoms to produce ring assemblies. Return the list of merged rings.
    :param rings: a list of rings, where each ring is represented by the list of atoms that make up the ring.
    :return: ring assemblies: a list of ring assemblies, where each ring is represented by the list of atoms that make up the ring assembly.
    """
    rings_list_copy = deepcopy(rings)

    # Function to merge rings with a common atom
    def merge_rings(ring_list):
        merged = False
        for i in range(len(ring_list)):
            for j in range(i + 1, len(ring_list)):
                common_atoms = set(ring_list[i]) & set(ring_list[j])
                if common_atoms:
                    ring_list[i].extend(ring_list[j])
                    ring_list[i] = list(set(ring_list[i]))  # Remove duplicates
                    ring_list.pop(j)
                    merged = True
                    return merged
        return merged

    # Merge rings iteratively until no more merging is possible
    while merge_rings(rings_list_copy):
        pass
    return rings_list_copy


def get_num_bridgehead_atoms(mol: Chem.rdchem.Mol) -> int:
    return rdMolDescriptors.CalcNumBridgeheadAtoms(mol)


def _ring_list_to_num_macrocycles(ring_list: list[list[int]]) -> int:
    ring_sizes = [len(ring) for ring in ring_list]
    return 1 if any(size > 12 for size in ring_sizes) else 0


def _rings_list_to_atom2ring_idx(ring_list: list[list[int]]) -> dict[int, int]:
    atom_assemblies = {}
    for i, assembly in enumerate(ring_list):
        for atom in assembly:
            atom_assemblies[atom] = i
    return atom_assemblies


def _get_ring_mapping(
    atom2ring_idx_dict: dict[int, int], atom2ring_assembly_dict: dict[int, int]
) -> dict[int, int]:
    ring_mapping = {}  # Initialize an empty dictionary to store the mapping

    for atom_index, initial_ring_index in atom2ring_idx_dict.items():
        final_ring_index = atom2ring_assembly_dict.get(atom_index, None)
        if final_ring_index is not None:
            ring_mapping[initial_ring_index] = final_ring_index

    return ring_mapping


def _count_keys_corresponding_to_values(ring_mapping: dict[int, int]) -> Counter[int]:
    value_counts = Counter(Counter(ring_mapping.values()).values())
    return value_counts


def _bin(num: int) -> int:
    """Binning scheme."""
    d = {0: 1, 1: 1, 2: 2, 3: 3, 4: 3, 5: 4, 6: 4}
    if num in d.keys():
        return d[num]
    else:
        # if the chain is longer than 6, the bin is 7
        return 7


def _get_the_four_smallest_values_binned(
    chain_lengths: list[int],
) -> tuple[int, int, int, int]:
    chain_lengths = sorted(chain_lengths)
    if len(chain_lengths) >= 4:
        return (
            _bin(chain_lengths[0]),
            _bin(chain_lengths[1]),
            _bin(chain_lengths[2]),
            _bin(chain_lengths[3]),
        )
    elif len(chain_lengths) == 3:
        return _bin(chain_lengths[0]), _bin(chain_lengths[1]), _bin(chain_lengths[2]), 0
    elif len(chain_lengths) == 2:
        return _bin(chain_lengths[0]), _bin(chain_lengths[1]), 0, 0
    elif len(chain_lengths) == 1:
        return _bin(chain_lengths[0]), 0, 0, 0
    else:
        return 0, 0, 0, 0


def mol_to_num_chain_assemblies(mol: Chem.rdchem.Mol) -> int | str:
    """
    Calculate the number of chain assemblies in a molecule.
    :param mol: Chem.rdchem.Mol
    :return: num chain assemblies or ERROR_SCINS if you didn't provide a rdkit.Chem.rdchem.Mol.
    """
    if not isinstance(mol, Chem.rdchem.Mol):
        _rdkit_mol_warning(mol)
        return ERROR_SCINS
    non_ring_mol_graph = _mol_to_non_ring_mol_graph(mol)
    num_chain_assemblies = _non_ring_mol_graph_to_num_chain_assemblies(
        non_ring_mol_graph
    )
    return num_chain_assemblies


def mol_to_chain_lengths(mol: Chem.rdchem.Mol) -> list[int] | str:
    """
    Calculate the chain lengths in a molecule.
    :param mol: Chem.rdchem.Mol
    :return: chain lengths or ERROR_SCINS if you didn't provide a rdkit.Chem.rdchem.Mol.
    """
    if not isinstance(mol, Chem.rdchem.Mol):
        _rdkit_mol_warning(mol)
        return ERROR_SCINS
    non_ring_mol_graph = _mol_to_non_ring_mol_graph(mol)
    chain_lengths = _non_ring_mol_graph_to_chain_lengths(non_ring_mol_graph)
    return chain_lengths


def mol_to_num_ring_assemblies(mol: Chem.rdchem.Mol) -> Counter[int]:
    """
    Calculate the number of ring assemblies in a molecule.
    :param mol: Chem.rdchem.Mol
    :return: dict[size_of_ring_assembly: count], where size of the ring assembly is the number of rings in a ring assembly.
    """
    rings_list = _mol_to_rings(mol)
    ring_assemblies_list = get_ring_assemblies(rings_list)
    atom2ring_idx = _rings_list_to_atom2ring_idx(rings_list)
    atom2ring_assembly_idx = _rings_list_to_atom2ring_idx(ring_assemblies_list)
    ring_mapping = _get_ring_mapping(atom2ring_idx, atom2ring_assembly_idx)
    num_rings_in_assemblies = _count_keys_corresponding_to_values(ring_mapping)
    return num_rings_in_assemblies


def generic_scaffold_mol_to_scins(mol: Chem.rdchem.Mol) -> str:
    """
    Calculate the SCINS for the generic scaffold mol of a molecule.
    Currently you have to decide whether you want to use the generic scaffold of the molecule
    as originally defined by Bemis and Murcko (1: GetScaffoldForMol, 2: MakeScaffoldGeneric),
    or use the scaffold that is trimmed (1: MakeScaffoldGeneric, 2: GetScaffoldForMol).
    Hence, this function directly takes the Chem.rdchem.Mol scaffold object as input.
    :param mol: The generic Murcko scaffold of a molecule.
    :return: str: the SCINS string.
    """
    if not isinstance(mol, Chem.rdchem.Mol):
        _rdkit_mol_warning(mol)
        return ERROR_SCINS
    rings_list = _mol_to_rings(mol)

    non_ring_mol_graph = _mol_to_non_ring_mol_graph(mol)
    num_chain_assemblies = _non_ring_mol_graph_to_num_chain_assemblies(
        non_ring_mol_graph
    )
    chain_lengths = _non_ring_mol_graph_to_chain_lengths(non_ring_mol_graph)

    ring_assemblies_list = get_ring_assemblies(rings_list)
    num_bridgehead_atoms = get_num_bridgehead_atoms(mol)

    part1 = "_".join(
        [
            str(num_chain_assemblies),
            str(len(chain_lengths)),
            str(len(rings_list)),
            str(len(ring_assemblies_list)),
            str(num_bridgehead_atoms),
        ]
    )

    atom2ring_idx = _rings_list_to_atom2ring_idx(rings_list)
    atom2ring_assembly_idx = _rings_list_to_atom2ring_idx(ring_assemblies_list)
    ring_index2ring_assembly_index = _get_ring_mapping(
        atom2ring_idx, atom2ring_assembly_idx
    )
    num_rings_in_assemblies = _count_keys_corresponding_to_values(
        ring_index2ring_assembly_index
    )
    num_assemblies_with_one_ring = num_rings_in_assemblies.get(1, 0)
    num_assemblies_with_two_rings = num_rings_in_assemblies.get(2, 0)
    num_assemblies_with_three_or_more_rings = 0  # instantiate variable for summation
    for k in num_rings_in_assemblies.keys():
        if k > 2:
            num_assemblies_with_three_or_more_rings += num_rings_in_assemblies[k]

    num_macrocycles = _ring_list_to_num_macrocycles(rings_list)
    part2 = "_".join(
        [
            str(num_assemblies_with_one_ring),
            str(num_assemblies_with_two_rings),
            str(num_assemblies_with_three_or_more_rings),
            str(num_macrocycles),
        ]
    )

    four_shortest_chains = _get_the_four_smallest_values_binned(chain_lengths)
    part3 = "_".join([str(i) for i in four_shortest_chains])

    return part1 + "-" + part2 + "-" + part3


def generic_scaffold_smiles_to_scins(smiles: str) -> str:
    """
    Calculate the SCINS for the generic scaffold SMILES of a molecule.
    Currently you have to decide whether you want to use the generic scaffold of the molecule
    as originally defined by Bemis and Murcko (1: GetScaffoldForMol, 2: MakeScaffoldGeneric),
    or use the scaffold that is trimmed (1: MakeScaffoldGeneric, 2: GetScaffoldForMol).
    Hence, this function directly takes the Chem.rdchem.Mol scaffold object as input.
    :param smiles: The generic Murcko scaffold SMILES of a molecule.
    :return: str: the SCINS string.
    """
    if not isinstance(smiles, str):
        logging.warning(
            "Input should be a SMILES str, but got an object of type %s. Therefore returning ERROR_SCINS."
            % type(smiles)
        )
        return ERROR_SCINS
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        _rdkit_mol_warning(mol)
        return ERROR_SCINS
    return generic_scaffold_mol_to_scins(mol)


def GetScaffoldForMol_edited(mol):
    """
    A modified version of GetScaffoldForMol in rdkit that also removes carbonyl (and similar groups bound by a double bond or triple bond).
    """
    murckoPatts = [
        "[!#1;D3;$([D3]-[!#1])](=[AD1])=[AD1]",
        "[!#1;D2;$([D2]-[!#1])]=,#[AD1]",
        "[!#1;D1;$([D1]~[!#1;!n])]",  ## KP slight modification here compared to rdkit
        "[3H]",
        "[2H]",
        "[O;D1;$([D1]=[D3])]",
    ]
    murckoQ_sm = "[" + ",".join(["$(%s)" % x for x in murckoPatts]) + "]"
    murckoQ = Chem.MolFromSmarts(murckoQ_sm)
    murckoPatts = [Chem.MolFromSmarts(x) for x in murckoPatts]

    while mol.HasSubstructMatch(murckoQ):
        for patt in murckoPatts:
            patt.UpdatePropertyCache()
            Chem.GetSymmSSSR(patt)
            patt.GetRingInfo()
            mol = Chem.DeleteSubstructs(mol, patt)
            mol.UpdatePropertyCache()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            if atom.GetNoImplicit() and atom.GetExplicitValence() < 4:
                atom.SetNoImplicit(False)
    h = Chem.MolFromSmiles("[H]")
    mol = Chem.ReplaceSubstructs(mol, Chem.MolFromSmarts("[D1;$([D1]-n)]"), h, True)[0]
    mol = Chem.RemoveHs(mol, sanitize=False)
    if mol.GetNumAtoms() < 3:
        return Chem.MolFromSmiles("")
    return mol

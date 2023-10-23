import logging
from collections import Counter
import itertools

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

EMPTY_SCINS = '00000-0000-0000'


def rdkit_mol_warning(mol):
    logging.warning(
        "Input should be an RDKit molecule object, but got an object of type %s. Therefore returning empty scins" % type(
            mol))


def _update_mol_graph_dict_inplace(atom_key, atom_bound, non_ring_mol_graph):
    if atom_key in non_ring_mol_graph:
        non_ring_mol_graph[atom_key].append(atom_bound)
    else:
        non_ring_mol_graph[atom_key] = [atom_bound]


def mol_to_non_ring_mol_graph(mol):
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


def _non_ring_mol_graph_to_num_chain_assemblies(non_ring_mol_graph):
    visited = set()
    num_chains = 0

    # Depth-First Search (DFS) to find the number of chains
    for atom in non_ring_mol_graph:
        if atom not in visited:
            num_chains += 1
            stack = [atom]
            while stack:
                node = stack.pop()
                visited.add(node)

                for neighbor in non_ring_mol_graph[node]:
                    if neighbor not in visited:
                        stack.append(neighbor)
    return num_chains


def _non_ring_mol_graph_to_chain_lengths(non_ring_mol_graph, rings):
    ring_atoms = list(itertools.chain.from_iterable(rings))
    visited = set()
    # chains = 0
    chain_lengths = []

    tertiary_c = []
    for atom in non_ring_mol_graph:
        nei = [n for n in non_ring_mol_graph[atom]]
        if len(nei) > 2:
            tertiary_c.append(atom)

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
                # neighbor_count = len([n for n in non_ring_mol_graph[node]]) # if n not in visited and n not in ring_atoms
                neighbors = non_ring_mol_graph[node]
                if len(neighbors) > 2: ## or n not in ring_atoms if we want to split the chains there
                    # If more than one neighbor, it's a branching point; terminate this chain
                    # if len(visited) == 1:
                    # if this is the first atom that we start with in the graph
                    # we need to correct for the number of chains
                    # chains = 0
                    continue

                for neighbor in non_ring_mol_graph[node]:
                    ## this looks inefficient, but actually needs to be like that to capture all edge cases
                    if len([n for n in non_ring_mol_graph[neighbor]]) >= 3:
                        chain_length += 1
                        if neighbor not in visited:
                            stack.append(neighbor)
                    elif neighbor not in visited:
                        stack.append(neighbor)
                        chain_length += 1

            chain_lengths.append(chain_length)

    return chain_lengths


def get_rings_for_mol(mol):
    ring_info = mol.GetRingInfo()
    rings = [list(ring) for ring in ring_info.AtomRings()]
    return rings


def get_num_ring_assemblies(rings):
    from copy import deepcopy
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


def get_num_bridgehead_atoms(mol):
    return rdMolDescriptors.CalcNumBridgeheadAtoms(mol)


def ring_list_to_num_macrocycles(ring_list):
    ring_sizes = [len(ring) for ring in ring_list]
    num_macrocycles = sum(1 for size in ring_sizes if size > 12)
    return num_macrocycles


def _rings_list_to_atom2ring_idx(ring_list):
    atom_assemblies = {}
    for i, assembly in enumerate(ring_list):
        for atom in assembly:
            atom_assemblies[atom] = i
    return atom_assemblies


def _get_ring_mapping(atom2ring_idx_dict, atom2ring_assembly_dict):
    ring_mapping = {}  # Initialize an empty dictionary to store the mapping

    for atom_index, initial_ring_index in atom2ring_idx_dict.items():
        final_ring_index = atom2ring_assembly_dict.get(atom_index, None)
        if final_ring_index is not None:
            ring_mapping[initial_ring_index] = final_ring_index

    return ring_mapping


def _count_keys_corresponding_to_values(ring_mapping):
    value_counts = Counter(Counter(ring_mapping.values()).values())
    return value_counts


def _bin(num):
    d = {
        0: 1, 1: 1, 2: 2, 3: 3, 4: 3, 5: 4, 6: 4
    }
    if num in d.keys():
        return d[num]
    else:
        # if the chain is longer than 6, the bin is 7
        return 7


def _get_the_four_smallest_values_binned(chain_lengths):
    chain_lengths = sorted(chain_lengths)
    print(chain_lengths)
    if len(chain_lengths) >= 4:
        return _bin(chain_lengths[0]), _bin(chain_lengths[1]), _bin(chain_lengths[2]), _bin(chain_lengths[3])
    elif len(chain_lengths) == 3:
        return _bin(chain_lengths[0]), _bin(chain_lengths[1]), _bin(chain_lengths[2]), 0
    elif len(chain_lengths) == 2:
        return _bin(chain_lengths[0]), _bin(chain_lengths[1]), 0, 0
    elif len(chain_lengths) == 1:
        return _bin(chain_lengths[0]), 0, 0, 0
    else:
        return 0, 0, 0, 0


def mol_to_num_chain_assemblies(mol):
    if not isinstance(mol, Chem.rdchem.Mol):
        rdkit_mol_warning(mol)
        return EMPTY_SCINS
    non_ring_mol_graph = mol_to_non_ring_mol_graph(mol)
    num_chain_assemblies = _non_ring_mol_graph_to_num_chain_assemblies(non_ring_mol_graph)
    return num_chain_assemblies


def mol_to_chain_lengths(mol):
    if not isinstance(mol, Chem.rdchem.Mol):
        rdkit_mol_warning(mol)
        return EMPTY_SCINS
    non_ring_mol_graph = mol_to_non_ring_mol_graph(mol)
    chain_lengths = _non_ring_mol_graph_to_chain_lengths(non_ring_mol_graph)
    return chain_lengths


def mol_to_num_ring_assemblies(mol):
    rings_list = get_rings_for_mol(mol)
    ring_assemblies_list = get_num_ring_assemblies(rings_list)
    atom2ring_idx = _rings_list_to_atom2ring_idx(rings_list)
    atom2ring_assembly_idx = _rings_list_to_atom2ring_idx(ring_assemblies_list)
    ring_mapping = _get_ring_mapping(atom2ring_idx, atom2ring_assembly_idx)
    num_rings_in_assemblies = _count_keys_corresponding_to_values(ring_mapping)
    return num_rings_in_assemblies



def scaffold_mol_to_scins(mol):
    """
    Currently you have to decide whether you want to use the generic scaffold of the molecule
    as defined in rdkit, or Use Kamen's recommendation to use the scaffold that is trimmed further.
    Hence, this function directly takes the Chem.rdchem.Mol representation of the scaffold as input.
    :param mol: The Generic Murcko Scaffold of the molecule
    :return: str: the SCINS string
    """
    if not isinstance(mol, Chem.rdchem.Mol):
        rdkit_mol_warning(mol)
        return EMPTY_SCINS
    rings_list = get_rings_for_mol(mol)

    non_ring_mol_graph = mol_to_non_ring_mol_graph(mol)
    num_chain_assemblies = _non_ring_mol_graph_to_num_chain_assemblies(non_ring_mol_graph)
    chain_lengths = _non_ring_mol_graph_to_chain_lengths(non_ring_mol_graph, rings_list)

    ring_assemblies_list = get_num_ring_assemblies(rings_list)
    num_bridgehead_atoms = get_num_bridgehead_atoms(mol)

    part1 = str(num_chain_assemblies) + str(len(chain_lengths)) + str(len(rings_list)) + str(
        len(ring_assemblies_list)) + str(num_bridgehead_atoms)

    atom2ring_idx = _rings_list_to_atom2ring_idx(rings_list)
    atom2ring_assembly_idx = _rings_list_to_atom2ring_idx(ring_assemblies_list)
    ring_index2ring_assembly_index = _get_ring_mapping(atom2ring_idx, atom2ring_assembly_idx)
    num_rings_in_assemblies = _count_keys_corresponding_to_values(ring_index2ring_assembly_index)
    num_assemblies_with_one_ring = num_rings_in_assemblies.get(1, 0)
    num_assemblies_with_two_rings = num_rings_in_assemblies.get(2, 0)
    num_assemblies_with_three_rings = num_rings_in_assemblies.get(3, 0)
    num_macrocycles = ring_list_to_num_macrocycles(rings_list)
    part2 = (str(num_assemblies_with_one_ring) + str(num_assemblies_with_two_rings) +
             str(num_assemblies_with_three_rings) + str(num_macrocycles))

    part3 = ''.join([str(i) for i in _get_the_four_smallest_values_binned(chain_lengths)])

    return part1 + '-' + part2 + '-' + part3


def smiles_to_scins(smiles):
    if not isinstance(smiles, str):
        logging.warning(
            "Input should be a SMILES str, but got an object of type %s. Therefore returning empty scins" % type(
                smiles))
        return EMPTY_SCINS
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        rdkit_mol_warning(mol)
        return EMPTY_SCINS
    return scaffold_mol_to_scins(mol)


def GetScaffoldForMol_edited(mol):
    murckoPatts = [
        '[!#1;D3;$([D3]-[!#1])](=[AD1])=[AD1]', '[!#1;D2;$([D2]-[!#1])]=,#[AD1]',
        '[!#1;D1;$([D1]~[!#1;!n])]'  ## KP slight modification here compared to rdkit
    ]
    murckoQ_sm = '[' + ','.join(['$(%s)' % x for x in murckoPatts]) + ']'
    murckoQ = Chem.MolFromSmarts(murckoQ_sm)
    murckoPatts = [Chem.MolFromSmarts(x) for x in murckoPatts]

    while mol.HasSubstructMatch(murckoQ):
        for patt in murckoPatts:
            mol = Chem.DeleteSubstructs(mol, patt)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            if atom.GetNoImplicit() and atom.GetExplicitValence() < 4:
                atom.SetNoImplicit(False)
    h = Chem.MolFromSmiles('[H]')
    mol = Chem.ReplaceSubstructs(mol, Chem.MolFromSmarts('[D1;$([D1]-n)]'), h, True)[0]
    mol = Chem.RemoveHs(mol, sanitize=False)
    return mol


if __name__ == "__main__":
    from rdkit.Chem.Scaffolds.MurckoScaffold import MakeScaffoldGeneric
    gen_scaffold = MakeScaffoldGeneric(GetScaffoldForMol_edited(Chem.MolFromSmiles("Cc(cccc1)c1NC(CCC(Nc(c(OC)c1)cc(OC)c1NC(c1ccco1)=O)=O)=O")))
    print(scaffold_mol_to_scins(gen_scaffold))


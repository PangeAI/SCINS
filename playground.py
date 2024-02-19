def get_ring_assemblies(rings_list):
    merged_rings = _merge_rings(rings_list)
    return merged_rings

def _merge_rings(ring_list):
    # Function to merge rings with a common atom
    merged_rings = []
    for i in range(len(ring_list)):
        merged = False
        for j in range(i + 1, len(ring_list)):
            common_atoms = set(ring_list[i]) & set(ring_list[j])
            if common_atoms:
                merged_ring = list(set(ring_list[i] + ring_list[j]))  # Merge and remove duplicates
                merged_rings.append(merged_ring)
                merged = True
        if not  merged:
            merged_rings.append(ring_list[i])
    return merged_rings


def mol_to_num_bridge_bonds(mol):
    rings_list = get_rings_for_mol(mol)
    num_bridge_bonds = _get_num_bridge_bonds(rings_list)
    return num_bridge_bonds

def get_rings_for_mol(mol):
    ring_info = mol.GetRingInfo()
    rings = [list(ring) for ring in ring_info.AtomRings()]
    return rings

def _get_num_bridge_bonds(ring_list):
    # can combine with one of the other ring iterators
    num_bridge_bonds = 0
    for i in range(len(ring_list)):
        for j in range(i + 1, len(ring_list)):
            ring1_atoms = set(ring_list[i])
            ring2_atoms = set(ring_list[j])
            common_atoms = ring1_atoms.intersection(ring2_atoms)
            if len(common_atoms) > 2:
                num_bridge_bonds += len(common_atoms) - 1
    return num_bridge_bonds


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
    if mol.GetNumAtoms() < 3:
        return Chem.MolFromSmiles('')
    return mol

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
                if len(neighbors) > 2:  ## or n not in ring_atoms if we want to split the chains there
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


## the method for calculating modified Murcko scaffolds
def GetScaffoldForMol_edited(mol):
    murckoPatts = [
        '[!#1;D3;$([D3]-[!#1])](=[AD1])=[AD1]',
        '[!#1;D2;$([D2]-[!#1])]=,#[AD1]',
        '[!#1;D1;$([D1]~[!#1;!n])]',  ## KP slight modification here compared to rdkit
        '[3H]', '[2H]', '[O;D1;$([D1]=[D3])]'
    ]
    murckoQ_sm = '[' + ','.join(['$(%s)' % x for x in murckoPatts]) + ']'
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
    h = Chem.MolFromSmiles('[H]')
    mol = Chem.ReplaceSubstructs(mol, Chem.MolFromSmarts('[D1;$([D1]-n)]'), h, True)[0]
    mol = Chem.RemoveHs(mol, sanitize=False)
    if mol.GetNumAtoms() < 3:
        return Chem.MolFromSmiles('')
    return mol

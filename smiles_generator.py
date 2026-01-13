from rdkit import Chem
from rdkit.Chem import AllChem
import random
import argparse

# R-Group Fragments for Addition
R_GROUP_FRAGMENTS = [
    "[*:1]C",        # Methyl
    "[*:1]N",        # Amine
    "[*:1]O",        # Hydroxyl
    "[*:1]C(=O)O",   # Carboxylic Acid
    "[*:1]F",        # Fluoro
    "[*:1]c1ccccc1", # Phenyl
    "[*:1]CC",       # Ethyl
]

H_QUERY = Chem.MolFromSmarts('[#1]') 

def generate_valid_random_molecule(core_smiles, num_modifications=3):
    mol = Chem.MolFromSmiles(core_smiles)
    if mol is None:
        print(f"Error: Invalid core SMILES provided: {core_smiles}")
        return None
    
    mol = Chem.AddHs(mol)

    for _ in range(num_modifications):
        h_matches = mol.GetSubstructMatches(H_QUERY)
        if not h_matches:
            break 

        target_h_idx = random.choice(h_matches)[0]
        frag_smiles_with_dummy = random.choice(R_GROUP_FRAGMENTS)
        r_group_pure_smiles = frag_smiles_with_dummy.replace("[*:1]", "")
        fragment_mol_to_add = Chem.MolFromSmiles(r_group_pure_smiles)

        if fragment_mol_to_add is None or fragment_mol_to_add.GetNumAtoms() == 0:
            continue

        try:
            rwmol = Chem.RWMol(mol)
            h_atom = rwmol.GetAtomWithIdx(target_h_idx)
            parent_atom = h_atom.GetNeighbors()[0]
            parent_idx = parent_atom.GetIdx()

            rwmol.RemoveAtom(target_h_idx)

            fragment_atom_map = {}
            for atom in fragment_mol_to_add.GetAtoms():
                fragment_atom_map[atom.GetIdx()] = rwmol.AddAtom(atom)

            if 0 in fragment_atom_map:
                rwmol.AddBond(parent_idx, fragment_atom_map[0], Chem.BondType.SINGLE)
            
            for bond in fragment_mol_to_add.GetBonds():
                rwmol.AddBond(fragment_atom_map[bond.GetBeginAtomIdx()], 
                              fragment_atom_map[bond.GetEndAtomIdx()], 
                              bond.GetBondType())

            mol = rwmol.GetMol()
            Chem.SanitizeMol(mol) 
            
        except Exception:
            return None 
                
    try:
        mol = Chem.RemoveHs(mol)
        Chem.SanitizeMol(mol)
    except:
        return None
        
    return mol

def main():
    DEFAULT_CORE_SMILES = "C1=CC=C2C(=C1)NC=C2" 

    parser = argparse.ArgumentParser(description="Generate random molecules and save to a file.")
    
    # ADDED: Required output argument
    parser.add_argument(
        '--out', 
        type=str, 
        required=True, 
        help='The filename to save the generated SMILES (e.g., results.smi).'
    )
    
    parser.add_argument(
        '--core', 
        type=str, 
        default=DEFAULT_CORE_SMILES, 
        help=f'Starting core SMILES (default: {DEFAULT_CORE_SMILES}).'
    )
    parser.add_argument(
        '--count', 
        type=int, 
        default=1, 
        help='Number of molecules to generate.'
    )
    parser.add_argument(
        '--max_attempts', 
        type=int, 
        default=1000, 
        help='Maximum attempts before stopping.'
    )
    
    args = parser.parse_args()
  
    CORE_SMILES = args.core
    TARGET_COUNT = args.count
    MAX_ATTEMPTS = args.max_attempts
    OUTPUT_FILE = args.out
    
    found_mols = []
    attempts = 0

    print(f"Starting generator from core: {CORE_SMILES}")
    print(f"Targeting: {TARGET_COUNT} molecule(s) -> Saving to: {OUTPUT_FILE}")

    while len(found_mols) < TARGET_COUNT:
        if attempts >= MAX_ATTEMPTS:
            print(f"\nWARNING: Reached maximum attempts ({MAX_ATTEMPTS}).")
            break 

        mol = generate_valid_random_molecule(CORE_SMILES, num_modifications=random.randint(2, 5))
        if mol is not None:
            found_mols.append(mol)
        attempts += 1

    if found_mols:
        final_smiles_list = [Chem.MolToSmiles(m) for m in found_mols]
        
        try:
            with open(OUTPUT_FILE, "w") as f:
                f.write("\n".join(final_smiles_list) + "\n")
            
            print("\n-------------------------------------------------")
            print(f"Success! {len(found_mols)} molecules written to: {OUTPUT_FILE}")
            print(f"Total attempts: {attempts}.")
            print("-------------------------------------------------")
            
        except Exception as e:
            print(f"Error writing file {OUTPUT_FILE}: {e}")
    else:
        print(f"\nFAILURE. No molecules generated after {MAX_ATTEMPTS} attempts.")

if __name__ == "__main__":
    main()
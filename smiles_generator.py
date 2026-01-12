from rdkit import Chem
from rdkit.Chem import AllChem
import random
import argparse

# --- R-Group Fragments for Addition (Molecular Growth) ---
# Each fragment uses a dummy atom [*:1] as the attachment point to the core.
R_GROUP_FRAGMENTS = [
    "[*:1]C",        # Methyl
    "[*:1]N",        # Amine
    "[*:1]O",        # Hydroxyl
    "[*:1]C(=O)O",   # Carboxylic Acid
    "[*:1]F",        # Fluoro
    "[*:1]c1ccccc1", # Phenyl
    "[*:1]CC",       # Ethyl
]

# Query for finding the replaceable Hydrogen atom
H_QUERY = Chem.MolFromSmarts('[#1]') 

def generate_valid_random_molecule(core_smiles, num_modifications=3):
    """
    Generates a valid random molecule by performing robust R-group additions 
    on a core structure, which is the proper way to achieve molecular growth.
    
    This version uses RWMol for reliable fragment addition (R-group enumeration).
    """
    mol = Chem.MolFromSmiles(core_smiles)
    if mol is None:
        print(f"Error: Invalid core SMILES provided: {core_smiles}")
        return None
    
    # 1. Add Hydrogens to the core to identify available substitution sites
    mol = Chem.AddHs(mol)

    # 2. Iterate and apply random modifications (R-group addition)
    for _ in range(num_modifications):
        
        # 2a. Find all substitution sites (Hydrogen atoms)
        h_matches = mol.GetSubstructMatches(H_QUERY)
        
        if not h_matches:
            # No Hydrogens left to replace, stop modifications
            break 

        # Select a random Hydrogen atom index to replace
        target_h_idx = random.choice(h_matches)[0]
        
        # 2b. Choose a random R-group fragment
        frag_smiles_with_dummy = random.choice(R_GROUP_FRAGMENTS)
        
        # Get the SMILES without the dummy atom [*:1]
        r_group_pure_smiles = frag_smiles_with_dummy.replace("[*:1]", "")
        
        # Get the fragment's molecule structure
        fragment_mol_to_add = Chem.MolFromSmiles(r_group_pure_smiles)

        if fragment_mol_to_add is None or fragment_mol_to_add.GetNumAtoms() == 0:
            continue # Skip if the fragment is invalid

        try:
            # --- R-Group Addition (Molecular Growth) ---
            rwmol = Chem.RWMol(mol)
            
            # Find the H atom and its parent heavy atom (the attachment point)
            h_atom = rwmol.GetAtomWithIdx(target_h_idx)
            parent_atom = h_atom.GetNeighbors()[0]
            parent_idx = parent_atom.GetIdx()

            # 1. Delete the H atom and its bond
            rwmol.RemoveAtom(target_h_idx)

            # 2. Transfer all atoms and bonds from the fragment to the core (rwmol)
            fragment_atom_map = {}
            for atom in fragment_mol_to_add.GetAtoms():
                # Add the atom to the RWMol and record the new index
                fragment_atom_map[atom.GetIdx()] = rwmol.AddAtom(atom)

            # 3. Bond the fragment's primary attachment atom (index 0) to the core's parent atom
            if 0 in fragment_atom_map:
                rwmol.AddBond(parent_idx, fragment_atom_map[0], Chem.BondType.SINGLE)
            
            # 4. Transfer all internal bonds within the fragment
            for bond in fragment_mol_to_add.GetBonds():
                rwmol.AddBond(fragment_atom_map[bond.GetBeginAtomIdx()], 
                              fragment_atom_map[bond.GetEndAtomIdx()], 
                              bond.GetBondType())

            # Convert back to a standard Mol object and sanitize
            mol = rwmol.GetMol()
            Chem.SanitizeMol(mol) 
            
        except Exception:
            # If the complex graph manipulation fails, discard
            return None 
                
    # 3. Final cleanup and return
    try:
        # Remove explicit Hs and sanitize the final structure
        mol = Chem.RemoveHs(mol)
        Chem.SanitizeMol(mol)
    except:
        return None
        
    return mol


def main():
    # 1. Define Default Core SMILES
    # Using the aromatic indole core for drug-like molecule generation
    DEFAULT_CORE_SMILES = "C1=CC=C2C(=C1)NC=C2" # A standard Indole SMILES

    # 2. Setup argparse
    parser = argparse.ArgumentParser(description="Generate a random, valid molecule starting from a core SMILES scaffold.")
    parser.add_argument(
        '--core', 
        type=str, 
        default=DEFAULT_CORE_SMILES, 
        help=f'Starting core SMILES pattern (default: {DEFAULT_CORE_SMILES} - Indole).'
    )
    # ARGUMENT: Number of molecules to generate
    parser.add_argument(
        '--count', 
        type=int, 
        default=1, 
        help='Number of valid molecules to generate (default: 1).'
    )
    # NEW ARGUMENT: Maximum attempts before giving up
    parser.add_argument(
        '--max_attempts', 
        type=int, 
        default=1000, 
        help='Maximum attempts before stopping the generation process (default: 1000).'
    )
    args = parser.parse_args()
    
    # Use the argument values
    CORE_SMILES = args.core
    TARGET_COUNT = args.count
    MAX_ATTEMPTS = args.max_attempts
    
    # List to store successfully generated molecules
    found_mols = []
    attempts = 0

    print(f"Starting generator from core: {CORE_SMILES} (Targeting {TARGET_COUNT} molecule(s))")

    # --- Running the Generator Loop ---
    # Loop until the target count is reached
    while len(found_mols) < TARGET_COUNT:
        if attempts >= MAX_ATTEMPTS:
            print(f"\nWARNING: Reached maximum attempts ({MAX_ATTEMPTS}) before generating {TARGET_COUNT} molecules.")
            break # Exit the loop if max attempts are hit

        # Generate 2 to 5 random modifications
        mol = generate_valid_random_molecule(CORE_SMILES, num_modifications=random.randint(2, 5))
        
        if mol is not None:
            found_mols.append(mol) # Append the molecule to the list
            print(f"  --> Found valid molecule #{len(found_mols)}...")
            
        attempts += 1

    # --- Final Output ---
    if found_mols: # Check if the list is non-empty
        # Convert all generated molecules to SMILES strings
        final_smiles_list = [Chem.MolToSmiles(m) for m in found_mols]
        
        # Save the SMILES patterns to the desired file
        try:
            # Write all SMILES strings separated by newlines
            with open("random_mol.smi", "w") as f:
                f.write("\n".join(final_smiles_list) + "\n")
            
            print("\n-------------------------------------------------")
            print(f"✅ SUCCESS! Generated {len(found_mols)} valid molecule(s).")
            print(f"Total attempts: {attempts}.")
            print(f"SMILES patterns written to random_mol.smi")
            print("-------------------------------------------------")
            
        except Exception as e:
            print(f"Error writing file random_mol.smi: {e}")
    else:
        print(f"\nFAILURE. Could not generate any valid molecules after {MAX_ATTEMPTS} attempts.")

if __name__ == "__main__":
    main()

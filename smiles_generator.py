from rdkit import Chem
from rdkit.Chem import AllChem
import random
import argparse

# R-Group Fragments for Addition
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
    Generates a valid random molecule by performing R-group additions 
    on a core structure.
    
    This version uses RWMol for fragment addition (R-group enumeration).
    """
    mol = Chem.MolFromSmiles(core_smiles)
    if mol is None:
        print(f"Error: Invalid core SMILES provided: {core_smiles}")
        return None
    
    #Add Hydrogens to the core to identify available substitution sites
    mol = Chem.AddHs(mol)

    #Iterate and apply random modifications (R-group addition)
    for _ in range(num_modifications):
        
        # Find all substitution sites (Hydrogen atoms)
        h_matches = mol.GetSubstructMatches(H_QUERY)
        
        if not h_matches:
            # No Hydrogens left to replace, stop modifications
            break 

        # Select a random Hydrogen atom index to replace
        target_h_idx = random.choice(h_matches)[0]
        
        #Choose a random R-group fragment
        frag_smiles_with_dummy = random.choice(R_GROUP_FRAGMENTS)
        
        # Get the SMILES without the dummy atom [*:1]
        r_group_pure_smiles = frag_smiles_with_dummy.replace("[*:1]", "")
        
        # Get the fragment's molecule structure
        fragment_mol_to_add = Chem.MolFromSmiles(r_group_pure_smiles)

        if fragment_mol_to_add is None or fragment_mol_to_add.GetNumAtoms() == 0:
            continue # Skip if the fragment is invalid

        try:
            # R-Group Addition
            rwmol = Chem.RWMol(mol)
            
            # Find the H atom and its parent heavy atom (the attachment point)
            h_atom = rwmol.GetAtomWithIdx(target_h_idx)
            parent_atom = h_atom.GetNeighbors()[0]
            parent_idx = parent_atom.GetIdx()

            #Delete the H atom and its bond
            rwmol.RemoveAtom(target_h_idx)

            #Transfer all atoms and bonds from the fragment to the core (rwmol)
            fragment_atom_map = {}
            for atom in fragment_mol_to_add.GetAtoms():
                # Add the atom to the RWMol and record the new index
                fragment_atom_map[atom.GetIdx()] = rwmol.AddAtom(atom)

            #Bond the fragment's primary attachment atom (index 0) to the core's parent atom
            if 0 in fragment_atom_map:
                rwmol.AddBond(parent_idx, fragment_atom_map[0], Chem.BondType.SINGLE)
            
            #Transfer all internal bonds within the fragment
            for bond in fragment_mol_to_add.GetBonds():
                rwmol.AddBond(fragment_atom_map[bond.GetBeginAtomIdx()], 
                              fragment_atom_map[bond.GetEndAtomIdx()], 
                              bond.GetBondType())

            #Convert back to a standard Mol object and sanitize
            mol = rwmol.GetMol()
            Chem.SanitizeMol(mol) 
            
        except Exception:
            # If the graph manipulation fails, discard
            return None 
                
    #Final cleanup and return
    try:
        # Remove explicit Hs and sanitize the final structure
        mol = Chem.RemoveHs(mol)
        Chem.SanitizeMol(mol)
    except:
        return None
        
    return mol


def main():
    # Using the aromatic indole core for drug-like molecule generation.
    # user can pass different core in args
    DEFAULT_CORE_SMILES = "C1=CC=C2C(=C1)NC=C2" # A standard Indole SMILES

    parser = argparse.ArgumentParser(description="Generate a random, valid molecule starting from a core SMILES scaffold.")
    parser.add_argument(
        '--core', 
        type=str, 
        default=DEFAULT_CORE_SMILES, 
        help=f'Starting core SMILES pattern (default: {DEFAULT_CORE_SMILES} - Indole).'
    )
    parser.add_argument(
        '--count', 
        type=int, 
        default=1, 
        help='Number of valid molecules to generate (default: 1).'
    )
    #allow setting maximum attempts before giving up
    parser.add_argument(
        '--max_attempts', 
        type=int, 
        default=1000, 
        help='Maximum attempts before stopping the generation process (default: 1000).'
    )
    args = parser.parse_args()
  
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

    # Final output
    if found_mols: # Check if the list is non-empty
        # Convert all generated molecules to SMILES strings
        final_smiles_list = [Chem.MolToSmiles(m) for m in found_mols]
        
        # Save the SMILES patterns to the desired file
        try:
            # Write all SMILES strings separated by newlines
            with open("random_mol.smi", "w") as f:
                f.write("\n".join(final_smiles_list) + "\n")
            
            print("\n-------------------------------------------------")
            print(f"Generated {len(found_mols)} valid molecule(s).")
            print(f"Total attempts: {attempts}.")
            print(f"SMILES patterns written to random_mol.smi")
            print("-------------------------------------------------")
            
        except Exception as e:
            print(f"Error writing file random_mol.smi: {e}")
    else:
        print(f"\nFAILURE. Could not generate any valid molecules after {MAX_ATTEMPTS} attempts.")

if __name__ == "__main__":
    main()

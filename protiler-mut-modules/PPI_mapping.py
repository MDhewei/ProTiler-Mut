import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


from Bio.PDB.PDBParser import PDBParser
import pandas as pd

def get_interface_contacts(pdb_file: str, target_chain_id: str, residue_number: int, distance_threshold: float = 5.0) -> list:
    """
    For a given PDB file, target chain, and residue number, return a list of contacting residues
    (from other chains) that have at least one atom within the given distance threshold of any atom
    in the target residue.
    
    Each contact is recorded as a dictionary with keys:
      - 'chain': chain identifier of the contacting residue,
      - 'residue_number': the residue number,
      - 'residue_name': the three-letter residue name.
    
    Parameters:
        pdb_file (str): Path to the PDB file.
        target_chain_id (str): Chain identifier where the target mutation is located.
        residue_number (int): Residue number of the mutation.
        distance_threshold (float): Distance threshold (in Å) for contact (default 5.0 Å).
    
    Returns:
        list: A list of dictionaries for contacting residues.
    """
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = structure[0]

    try:
        target_chain = model[target_chain_id]
    except KeyError:
        print(f"Chain {target_chain_id} not found in {pdb_file}")
        return []

    target_residue = None
    for res in target_chain.get_residues():
        if res.id[1] == residue_number:
            target_residue = res
            break
    if target_residue is None:
        print(f"Residue {residue_number} not found in chain {target_chain_id} of {pdb_file}")
        return []

    contacts = []
    for other_chain in model:
        if other_chain.id == target_chain_id:
            continue
        for other_res in other_chain.get_residues():
            found_contact = False
            for atom1 in target_residue.get_atoms():
                for atom2 in other_res.get_atoms():
                    if atom1 - atom2 < distance_threshold:
                        found_contact = True
                        break
                if found_contact:
                    break
            if found_contact:
                contact_info = {
                    'chain': other_chain.id,
                    'residue_number': other_res.id[1],
                    'residue_name': other_res.get_resname()
                }
                if contact_info not in contacts:
                    contacts.append(contact_info)
    return contacts


def build_ppi_interface_table(df_gene: pd.DataFrame, pdb_files: list,chain_list: list, distance_threshold: float = 5.0, output_folder:str=None) -> pd.DataFrame:
    """
    Build a table that records, for each mutation (given as a residue number and corresponding chain)
    and each PDB file, whether the mutation is at the protein–protein interface, the number of contacting residues,
    and details of those contacts.
    
    Parameters:
        mut_aa_list (list): List of mutation residue numbers (integers).
        chain_list (list): List of chain identifiers (strings) corresponding to each mutation.
        pdb_files (list): List of PDB file paths representing protein–protein complex structures.
        distance_threshold (float): Distance threshold (in Å) for considering a contact.
    
    Returns:
        pd.DataFrame: A table with columns:
            ['PDB_File', 'Mutation_Chain', 'Mutation', 'At_Interface', 'Contact_Count', 'Contact_Details'].
            'At_Interface' is True if Contact_Count > 0.
    """
    
    mut_aa_list = list(df_gene['AA'])
    
    records = []
    for pdb in pdb_files:
        for mut, chain in zip(mut_aa_list, chain_list):
            contacts = get_interface_contacts(pdb, chain, mut, distance_threshold)
            contact_count = len(contacts)
            at_interface = contact_count > 0
            contact_details = "; ".join(
                [f"{c['chain']}:{c['residue_number']}-{c['residue_name']}" for c in contacts]
            )
            records.append({
                'PDB_File': pdb,
                'Mutation_Chain': chain,
                'Mutation': mut,
                'At_Interface': at_interface,
                'Contact_Count': contact_count,
                'Contact_Details': contact_details
            })
     
    df_interface = pd.DataFrame(records)
    
    if output_folder is not None:
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        
        csv_path = os.path.join(os.getcwd(), output_folder, "ppi_interface_mutations.csv")
        df_interface.to_csv(csv_path, index=False)
        print(f"Interface table saved to: {csv_path}")
    
    return df_interface


# Example usage:
if __name__ == "__main__":
    # List of mutation residue numbers (e.g., [45, 34, 66])
    mutation_numbers = [45, 34, 66]
    # Corresponding list of chain identifiers (e.g., ['A', 'B', 'A'])
    mutation_chains = ['A', 'B', 'A']

    # List of PDB file paths for protein–protein complexes
    pdb_list = ["complex1.pdb", "complex2.pdb", "complex3.pdb"]

    # Build the interface table.
    interface_table = build_ppi_interface_table(mutation_numbers, mutation_chains, pdb_list, distance_threshold=5.0)

    # Optionally, filter to include only entries where the mutation is at an interface.
    interface_hits = interface_table[interface_table['At_Interface'] == True]

    print("Full interface table:")
    print(interface_table)
    print("\nMutations at interface:")
    print(interface_hits)
    
    # Save the table to CSV if desired.
    interface_table.to_csv("ppi_interface_mutations.csv", index=False)

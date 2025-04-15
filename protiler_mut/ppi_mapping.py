import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


from Bio.PDB.PDBParser import PDBParser
import pandas as pd
import os
from pymol import cmd

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


def build_ppi_interface_table(gene,mut_aa_list:list, pdb_files: list,chain_list: list, distance_threshold: float = 5.0) -> pd.DataFrame:
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
    
   
    records = []
    for pdb,chain in zip(pdb_files,chain_list):
        for mut in mut_aa_list:
            contacts = get_interface_contacts(pdb, chain, int(mut), distance_threshold)
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
    df_interface['Gene'] = [gene]*df_interface.shape[0]
    
#     if output_folder is not None:
#         if not os.path.exists(output_folder):
#             os.makedirs(output_folder)
        
#         csv_path = os.path.join(os.getcwd(), output_folder, gene+"_ppi_interface_mutations.csv")
#         df_interface.to_csv(csv_path, index=False)
#         print(f"Interface table saved to: {csv_path}")
    
    return df_interface

def visualize_interfaces_for_pdb(group_df: pd.DataFrame, output_folder: str):
    """
    For a given PDB file (grouped by the interface table), load the structure in PyMOL
    and visualize all true interface items (mutations with contacts).
    
    For each row in group_df, the anchor (mutated) residue will be highlighted in magenta,
    and all contacting residues will be shown in red. The overall cartoon is shown with 0.8 transparency.
    
    Parameters:
        group_df (pd.DataFrame): Subset of the interface table for a single PDB file.
                                 Must include columns: 'Mutation_Chain', 'Mutation', and 'Contact_Details'.
        output_folder (str): Directory in which the PyMOL session will be saved.
    """
    # Reinitialize the PyMOL session.
    cmd.reinitialize()
    
    # Get the current PDB file (should be the same for all rows in the group).
    pdb_file = group_df['PDB_File'].iloc[0]
    
    # Load the structure.
    cmd.load(pdb_file, "protein")
    cmd.hide("everything", "protein")
    cmd.bg_color('white')
    cmd.show("cartoon", "protein")
    cmd.set("cartoon_transparency", 0.8)
    
    # Color chains: target chain from anchor mutations will be colored blue; others green.
    # (Assumes all anchors in this file belong to the same target chain.)
    target_chain = group_df['Mutation_Chain'].iloc[0]
    cmd.color("red", f"chain {target_chain}")
    cmd.color("blue", f"not chain {target_chain}")
    
    # Iterate over each interface record in the group.
    for idx, row in group_df.iterrows():
        mutation = row['Mutation']
        contact_details = row['Contact_Details']
        
        # Highlight the anchor residue for this record.
        anchor_sel = f"chain {target_chain} and resi {mutation}"
        anchor_name = f"anchor_{target_chain}_{mutation}_{idx}"
        cmd.select(anchor_name, anchor_sel)
        cmd.show("sticks", anchor_name)
        cmd.color("red", anchor_name)
        
        # Parse and highlight contacts.
        if contact_details:
            contacts = [x.strip() for x in contact_details.split(';') if x.strip()]
            for contact in contacts:
                try:
                    chain_part, res_part = contact.split(":")
                    res_num, _ = res_part.split("-")
                    sel_name = f"contact_{chain_part.strip()}_{res_num.strip()}_{idx}"
                    sel_str = f"chain {chain_part.strip()} and resi {res_num.strip()}"
                    cmd.select(sel_name, sel_str)
                    cmd.show("sticks", sel_name)
                    cmd.color("blue", sel_name)
                except Exception as e:
                    print(f"Error parsing contact detail '{contact}': {e}")
    
    # Save the session for this PDB file.
    session_name = os.path.basename(pdb_file) + "_interfaces.pse"
    output_name = os.path.join(os.getcwd(), output_folder, session_name)
    cmd.save(output_name)
    print(f"Saved PyMOL session for {pdb_file} to {output_name}")
    cmd.reinitialize()


def visualize_interfaces(interface_table: pd.DataFrame, output_folder: str):
    """
    Group the interface table by PDB file and visualize all interface mutations for each PDB.
    Only rows where At_Interface is True will be processed.
    
    Parameters:
        interface_table (pd.DataFrame): DataFrame with interface mapping results.
            Must include columns: 'PDB_File', 'Mutation_Chain', 'Mutation', 'At_Interface', 'Contact_Details'.
        output_folder (str): Directory to save the PyMOL sessions.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Filter for rows with At_Interface True.
    true_df = interface_table[interface_table['At_Interface'] == True]
    
    # Group by PDB file.
    grouped = true_df.groupby('PDB_File')
    for pdb_file, group_df in grouped:
        visualize_interfaces_for_pdb(group_df, output_folder)


# Example usage:
if __name__ == "__main__":
    
    mutation_numbers = set([32, 14, 37, 39, 47, 47, 49, 44, 43, 51, 43, 43, 68, 68, 36, 35])

    pdb_list = ["./Input_smaples/PNKP_HUMAN__XRCC1_HUMAN__1154aa_unrelaxed_af2mv3_model_1.pdb",
                "./Input_smaples/PNKP_HUMAN__XRCC1_HUMAN__1154aa_unrelaxed_af2mv3_model_2.pdb", 
                "./Input_smaples/PNKP_HUMAN__XRCC1_HUMAN__1154aa_unrelaxed_af2mv3_model_4.pdb"]

    mutation_chains = ['A', 'A', 'A']

    interface_table = build_ppi_interface_table(mutation_numbers,  pdb_list, mutation_chains,distance_threshold=5.0)
    visualize_interfaces(interface_table, 'PPI_3D_visualization')

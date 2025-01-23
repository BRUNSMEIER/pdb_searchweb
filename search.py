import os
from Bio.PDB import MMCIFParser
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.Data.IUPACData import protein_letters_3to1

# Define output folder
base_output_dir = os.path.abspath("pdb_files")
os.makedirs(base_output_dir, exist_ok=True)

# Representative PDB sequences
representatives = ["6KSHD", "4ctaA", "2x14A"]

for pdb_code in representatives:
    pdb_id = pdb_code[:4].lower()  # Extract PDB ID
    chain_id = pdb_code[4].upper()  # Extract chain ID

    # Determine CIF file path
    cif_file = os.path.join(base_output_dir, f"{pdb_id}.cif")
    if not os.path.exists(cif_file):
        print(f"File not found: {cif_file}")
        continue

    # Parse CIF file
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, cif_file)
    except Exception as e:
        print(f"Error parsing CIF file {cif_file}: {e}")
        continue

    # Extract sequence from the chain
    sequence = []
    found_chain = False
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                found_chain = True
                print(f"Processing chain {chain.id}")
                for residue in chain:
                    if residue.id[0] == " ":  # Exclude non-amino acids
                        residue_name = residue.resname.strip().capitalize()  # Convert to capitalize format
                        try:
                            if residue_name in protein_letters_3to1:
                                sequence.append(protein_letters_3to1[residue_name])
                            else:
                                print(f"Skipping unknown residue: {residue_name}")
                        except KeyError:
                            print(f"Skipping unknown residue: {residue_name}")
                break

    if not found_chain:
        print(f"Chain {chain_id} not found in PDB {pdb_id}")
        continue

    # Save sequence as FASTA file
    if sequence:
        fasta_sequence = "".join(sequence)
        fasta_file = os.path.join(base_output_dir, f"{pdb_id}_{chain_id}.fasta")
        record = SeqRecord(Seq(fasta_sequence), id=f"{pdb_id}_{chain_id}", description="")
        with open(fasta_file, "w") as f:
            FastaIO.FastaWriter(f).write_record(record)
        print(f"FASTA file saved: {fasta_file}")
    else:
        print(f"Sequence for {pdb_id}_{chain_id} is empty, no FASTA file saved")

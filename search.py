import os
import gzip
from Bio.PDB import PDBList
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO

# Use the HTTPS server to avoid old FTP mirror issues
pdb_list = PDBList(server="https://files.rcsb.org/")

base_output_dir = os.path.abspath("pdb_files")
os.makedirs(base_output_dir, exist_ok=True)

representatives = ["6KSHD", "4ctaA", "2x14A"]

for pdb_code in representatives:
    pdb_id = pdb_code[:4].lower()  # e.g. "6ksh", "4cta", "2x14"
    chain = pdb_code[4].upper()    # e.g. "D", "A", "A"

    try:
        pdb_file = pdb_list.retrieve_pdb_file(
            pdb_id,
            file_format="mmCif",
            pdir=base_output_dir,
            overwrite=True
        )
        print(f"Biopython returned path: {pdb_file}")
    except Exception as e:
        print(f"Error downloading {pdb_id}: {e}")
        continue

    # Check if the file or gz version exists
    if not os.path.exists(pdb_file):
        if os.path.exists(pdb_file + ".gz"):
            pdb_file += ".gz"
        else:
            print(f"File not found after download: {pdb_file}")
            continue

    # Parse mmCIF
    parser = MMCIFParser(QUIET=True)
    try:
        if pdb_file.endswith(".gz"):
            with gzip.open(pdb_file, "rt") as handle:
                structure = parser.get_structure(pdb_id, handle)
        else:
            structure = parser.get_structure(pdb_id, pdb_file)

        # Extract chain sequence
        for model in structure:
            for chain_obj in model:
                if chain_obj.id == chain:
                    sequence = "".join(
                        residue.resname
                        for residue in chain_obj.get_residues()
                        if residue.id[0] == " "
                    )
                    record = SeqRecord(Seq(sequence),
                                       id=f"{pdb_id}_{chain}",
                                       description="")
                    fasta_file = os.path.join(base_output_dir,
                                              f"{pdb_id}_{chain}.fasta")
                    with open(fasta_file, "w") as out_fasta:
                        FastaIO.FastaWriter(out_fasta).write_record(record)
                    print(f"Saved FASTA file: {fasta_file}")

    except Exception as e:
        print(f"Error processing {pdb_id}: {e}")

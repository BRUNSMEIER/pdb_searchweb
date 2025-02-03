from Bio.PDB import MMCIFParser
from Bio.Data.IUPACData import protein_letters_3to1
# .cif path
cif_file = "D:\\phas0077\\pdb_files\\6ksh.cif"  # 替换为目标文件路径

# analyze
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("6ksh", cif_file)

# print chain info
print("Chain ID INFO:")
for model in structure:
    for chain in model:
        print(f"Chain ID: {chain.id}")
print(protein_letters_3to1)
print(protein_letters_3to1["ALA"])  # output: A
print(protein_letters_3to1["GLY"])  # output: G
print(protein_letters_3to1["LYS"])  # output: K


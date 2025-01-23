from Bio.PDB import MMCIFParser
from Bio.Data.IUPACData import protein_letters_3to1
# 确定 .cif 文件路径
cif_file = "D:\\phas0077\\pdb_files\\6ksh.cif"  # 替换为目标文件路径

# 解析 .cif 文件
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure("6ksh", cif_file)

# 打印链信息
print("文件中的链 ID:")
for model in structure:
    for chain in model:
        print(f"链 ID: {chain.id}")
print(protein_letters_3to1)
print(protein_letters_3to1["ALA"])  # 输出: A
print(protein_letters_3to1["GLY"])  # 输出: G
print(protein_letters_3to1["LYS"])  # 输出: K


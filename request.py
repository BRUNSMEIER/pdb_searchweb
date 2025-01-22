import requests
import os

# 定义提交的URL和参数
url = "http://biomine.cs.vcu.edu/servers/NsitePred/"
fasta_files = [f for f in os.listdir() if f.endswith(".fasta")]

for fasta_file in fasta_files:
    with open(fasta_file, 'r') as f:
        fasta_content = f.read()

    # 模拟表单提交
    response = requests.post(url, files={"sequence": fasta_content})

    if response.status_code == 200:
        result_file = fasta_file.replace(".fasta", "_result.html")
        with open(result_file, "w") as f:
            f.write(response.text)
        print(f"Result saved for {fasta_file}")
    else:
        print(f"Failed to process {fasta_file}")

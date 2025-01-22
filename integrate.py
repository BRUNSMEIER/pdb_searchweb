import os
import requests

# Define directories
base_dir = os.path.abspath("pdb_files")
fasta_files = [os.path.join(base_dir, f) for f in os.listdir(base_dir) if f.endswith(".fasta")]

# URL for form submission
url = "http://biomine.cs.vcu.edu/servers/NsitePred/"

# Define your email address
your_email = "schonnbiukmit1@gmail.com"  # Replace this with your actual email address

# Function to submit the fasta files
def submit_fasta(file_path, email):
    with open(file_path, 'r') as f:
        fasta_content = f.read()

    # Submit the form with both the sequence and the email
    response = requests.post(url, data={"sequence": fasta_content, "email": email})

    if response.status_code == 200:
        print(f"Submission successful for {file_path}. Check your email for results.")
    else:
        print(f"Failed to process {file_path}. Status code: {response.status_code}")

# Process each fasta file
for fasta_file in fasta_files:
    print(f"Processing {fasta_file}...")
    submit_fasta(fasta_file, your_email)

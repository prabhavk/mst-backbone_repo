from Bio import Entrez, SeqIO
import time

# Set your email (required by NCBI)
Entrez.email = "prabhavk.gmail.com"

# Search for HA gene sequences of influenza A H3N2
search_term = '(("H3N2 subtype"[Organism] OR H3N2[All Fields]) AND HA[All Fields] AND gene[All Fields]) AND (ddbj_embl_genbank[filter] AND ("1701"[SLEN] : "1701"[SLEN]))'

# Search in NCBI
handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=30000)
record = Entrez.read(handle)
handle.close()
ids = record["IdList"]
print(f"Found {len(ids)} sequences.")

min_length = 1701
batch_size = 1000  # Adjust as needed
pause_duration = 5  # seconds

# Process in batches
for i in range(0, len(ids), batch_size):
    batch_ids = ids[i:i + batch_size]
    batch_sequences = []

    # Fetch sequences for the current batch
    handle = Entrez.efetch(db="nucleotide", id=batch_ids, rettype="gb", retmode="text")
    records = SeqIO.parse(handle, "genbank")
    for record in records:
        if len(record.seq) >= min_length:
            batch_sequences.append(record)
    handle.close()

    # Save current batch to a FASTA file
    batch_filename = f"HA_sequences_batch_{i//batch_size +1}.fasta"
    with open(batch_filename, "w") as out_handle:
        SeqIO.write(batch_sequences, out_handle, "fasta")
    print(f"Saved {len(batch_sequences)} sequences to {batch_filename}")

    # Pause between batches, except after the last one
    if i + batch_size < len(ids):
        time.sleep(pause_duration)

print("Processing complete.")

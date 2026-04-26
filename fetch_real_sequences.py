"""
fetch_real_sequences.py

Fetches real, publicly available reference sequences from NCBI GenBank.
We use influenza A (H1N1) segments; sequences that have been publicly
available in scientific literature for 20+ years and are standard
reference material in biosecurity research.

These are NOT dangerous sequences. They are published reference genomes
used by biosecurity tools including SecureDNA and IBBIS for validation.
"""

import time
from Bio import Entrez, SeqIO
from pathlib import Path

# NCBI requires an email for API access
Entrez.email = "hackathon.research@example.com"

# ─────────────────────────────────────────────────────────────
# NCBI ACCESSION NUMBERS
# These are publicly available, published H1N1 influenza segments
# from the NCBI GenBank database. Standard biosecurity reference data.
# Source: NCBI Influenza Virus Resource
# ─────────────────────────────────────────────────────────────

REFERENCE_ACCESSIONS = {
    "H1N1_HA_segment": "AF117241",   # Hemagglutinin — A/South Carolina/1/18
    "H1N1_NP_segment": "AF116575",   # Nucleoprotein — A/Brevig Mission/1/1918
    "H1N1_NA_segment": "AF250356",   # Neuraminidase — A/Brevig Mission/1/1918
}

def fetch_sequence(accession, retries=3):
    """Fetch a single sequence from NCBI GenBank by accession number."""
    for attempt in range(retries):
        try:
            print(f"  Fetching {accession}...")
            handle = Entrez.efetch(
                db="nucleotide",
                id=accession,
                rettype="fasta",
                retmode="text"
            )
            record = SeqIO.read(handle, "fasta")
            handle.close()
            time.sleep(0.4)  # NCBI rate limit: max 3 requests/second
            return str(record.seq), record.description
        except Exception as e:
            print(f"  Attempt {attempt+1} failed for {accession}: {e}")
            time.sleep(2)
    return None, None

def fetch_all_references():
    Path("data/references").mkdir(parents=True, exist_ok=True)
    sequences = {}

    print("\nFetching real reference sequences from NCBI GenBank...")
    print("(These are publicly available published sequences)\n")

    for name, accession in REFERENCE_ACCESSIONS.items():
        seq, desc = fetch_sequence(accession)
        if seq:
            sequences[name] = {
                "accession": accession,
                "description": desc,
                "sequence": seq,
                "length": len(seq)
            }
            # Save individual FASTA files
            fasta_path = f"data/references/{name}.fasta"
            with open(fasta_path, "w") as f:
                f.write(f">{accession} {desc}\n{seq}\n")
            print(f"  ✅ {name}: {len(seq)} bp")
        else:
            print(f"  ❌ Failed to fetch {accession}")

    return sequences

def summarize(sequences):
    print("\n" + "="*60)
    print("  REAL REFERENCE SEQUENCES FETCHED")
    print("="*60)
    for name, data in sequences.items():
        print(f"\n  {name}")
        print(f"    Accession  : {data['accession']}")
        print(f"    Length     : {data['length']} bp")
        print(f"    Description: {data['description'][:60]}...")
    print("\n  Saved to data/references/")
    print("="*60)

if __name__ == "__main__":
    sequences = fetch_all_references()
    if sequences:
        summarize(sequences)
        print("\n✅ Real sequences ready. Run real_dataset_generator.py next.")
    else:
        print("\n❌ No sequences fetched. Check your internet connection.")
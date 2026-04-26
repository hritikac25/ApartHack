"""
debug_scorer.py
Diagnoses why FARS scores are capping at 0.150 on real sequences.
"""

import json
import pandas as pd
from scorer import align_fragment, kmer_set, compute_coverage, compute_contiguity, REFERENCE_SEQUENCE, REF_LEN

# Load one real suspicious order
df = pd.read_csv("data/real_fragment_orders.csv")

# Pick the first high-coverage suspicious order
sus = df[df.label_name == "suspicious_high"].iloc[0]
fragments = json.loads(sus.fragments)

print(f"Order: {sus.order_id}")
print(f"Target segment: {sus.target_segment}")
print(f"Target accession: {sus.target_accession}")
print(f"Coverage declared: {sus.coverage}")
print(f"N fragments: {len(fragments)}")
print(f"Fragment lengths: {[len(f) for f in fragments]}")

print(f"\nFARS reference sequence length: {REF_LEN} bp")
print(f"FARS reference (first 60bp): {REFERENCE_SEQUENCE[:60]}")

print(f"\nFirst fragment (first 60bp): {fragments[0][:60]}")

# Test alignment of first fragment
start, end, score = align_fragment(fragments[0], REFERENCE_SEQUENCE)
print(f"\nAlignment of fragment[0] against FARS fictional reference:")
print(f"  Best position: {start}-{end}")
print(f"  Similarity score: {score:.4f}")

# Now load the REAL reference and test against it
from Bio import SeqIO
real_refs = {}
import glob
for f in glob.glob("data/references/*.fasta"):
    record = SeqIO.read(f, "fasta")
    real_refs[record.id] = str(record.seq).upper()

print(f"\nReal references loaded: {list(real_refs.keys())}")

# Test alignment against real reference
for acc, real_seq in real_refs.items():
    start, end, score = align_fragment(fragments[0], real_seq)
    print(f"\nAlignment of fragment[0] against {acc} ({len(real_seq)} bp):")
    print(f"  Best position: {start}-{end}")
    print(f"  Similarity score: {score:.4f}")
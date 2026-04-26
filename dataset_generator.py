"""
dataset_generator.py
Generates a synthetic dataset of DNA fragment orders; both benign and
suspicious (fragments that could assemble into a dangerous sequence).
This is the foundation of the Fragment Assembly Risk Scorer.
"""

import random
import json
import pandas as pd
from pathlib import Path

random.seed(42)

# ──────────────────────────────────────────────
# REFERENCE SEQUENCES
# These are simplified, fictional sequences used
# for research simulation. Not real pathogens.
# ──────────────────────────────────────────────

# A fictional "sequence of concern" — 600bp total
# In a real tool this would reference a secure DB
FICTIONAL_DANGEROUS_SEQUENCE = (
    "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTT"
    "AATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAA"
    "TTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTC"
    "AATGCTTTTCAAGATACCCAGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTA"
    "TGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAA"
    "GGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGAC"
    "ACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAA"
    "GTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATA"
    "CTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAA"
)

# Benign reference sequences
BENIGN_SEQUENCES = [
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGT",
    "GCTAGCATGACTGGTGGACAGCAAATGGGTCGGGATCTGTACGACGATGACGATAAGGATCCGGCTGCT",
    "AACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGAC",
    "TTCGAGCAAGAGATCGAGCACAGTGGCGGCCGCTCGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTC",
    "CCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGC",
    "GATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCC",
    "ACCCTCGTGACCACCTTCGGCTACGGCCTGCAGTGCTTCGCCCGCTACCCCGACCACATGAAGCAGCAC",
    "GACTTCTTCAAGAGCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGC",
]

def generate_fragment(sequence, start, length, add_noise=False):
    """Cut a fragment from a sequence, optionally adding point mutations."""
    end = min(start + length, len(sequence))
    fragment = sequence[start:end]
    if add_noise:
        # Simulate slight mutations (evading naive exact-match screening)
        fragment = list(fragment)
        n_mutations = max(1, int(len(fragment) * 0.03))  # 3% mutation rate
        for _ in range(n_mutations):
            pos = random.randint(0, len(fragment) - 1)
            fragment[pos] = random.choice("ATCG")
        fragment = "".join(fragment)
    return fragment

def generate_benign_order(n_fragments=None):
    """Generate a benign order: random fragments from housekeeping genes."""
    if n_fragments is None:
        n_fragments = random.randint(1, 4)
    fragments = []
    for _ in range(n_fragments):
        seq = random.choice(BENIGN_SEQUENCES)
        start = random.randint(0, max(0, len(seq) - 50))
        length = random.randint(30, 80)
        fragments.append({
            "sequence": generate_fragment(seq, start, length),
            "length": length,
            "claimed_purpose": random.choice([
                "GFP expression vector",
                "CRISPR guide RNA",
                "Vaccine antigen cloning",
                "Protein expression tag",
                "PCR primer design",
                "Reporter gene construct"
            ])
        })
    return fragments

def generate_suspicious_order(coverage_target=0.85, add_noise=False):
    """
    Generate a suspicious order: overlapping fragments that together
    cover a large portion of the dangerous reference sequence.
    This simulates a split-order attack.
    """
    fragments = []
    seq = FICTIONAL_DANGEROUS_SEQUENCE
    seq_len = len(seq)
    target_len = int(seq_len * coverage_target)

    pos = 0
    while pos < target_len:
        length = random.randint(80, 150)
        overlap = random.randint(10, 30)  # fragments overlap slightly
        fragment_seq = generate_fragment(seq, pos, length, add_noise=add_noise)
        fragments.append({
            "sequence": fragment_seq,
            "length": len(fragment_seq),
            "claimed_purpose": random.choice([
                "Synthetic gene study",
                "Academic research",
                "Protein engineering",
                "Gene therapy vector",
                "Metabolic pathway study"
            ])
        })
        pos += length - overlap

    return fragments

def build_dataset(n_benign=150, n_suspicious=100, n_noisy_suspicious=50):
    """
    Build the full labeled dataset.
    Labels:
        0 = benign order
        1 = suspicious (clean fragments)
        2 = suspicious (mutated fragments — evading naive screening)
    """
    records = []

    print(f"Generating {n_benign} benign orders...")
    for i in range(n_benign):
        frags = generate_benign_order()
        records.append({
            "order_id": f"BEN-{i:04d}",
            "label": 0,
            "label_name": "benign",
            "n_fragments": len(frags),
            "total_bp": sum(f["length"] for f in frags),
            "fragments": json.dumps([f["sequence"] for f in frags]),
            "purposes": json.dumps([f["claimed_purpose"] for f in frags]),
        })

    print(f"Generating {n_suspicious} suspicious orders (clean)...")
    for i in range(n_suspicious):
        coverage = random.uniform(0.75, 0.98)
        frags = generate_suspicious_order(coverage_target=coverage, add_noise=False)
        records.append({
            "order_id": f"SUS-{i:04d}",
            "label": 1,
            "label_name": "suspicious_clean",
            "n_fragments": len(frags),
            "total_bp": sum(f["length"] for f in frags),
            "fragments": json.dumps([f["sequence"] for f in frags]),
            "purposes": json.dumps([f["claimed_purpose"] for f in frags]),
        })

    print(f"Generating {n_noisy_suspicious} suspicious orders (mutated)...")
    for i in range(n_noisy_suspicious):
        coverage = random.uniform(0.75, 0.95)
        frags = generate_suspicious_order(coverage_target=coverage, add_noise=True)
        records.append({
            "order_id": f"MUT-{i:04d}",
            "label": 2,
            "label_name": "suspicious_mutated",
            "n_fragments": len(frags),
            "total_bp": sum(f["length"] for f in frags),
            "fragments": json.dumps([f["sequence"] for f in frags]),
            "purposes": json.dumps([f["claimed_purpose"] for f in frags]),
        })

    random.shuffle(records)
    df = pd.DataFrame(records)
    return df

if __name__ == "__main__":
    Path("data").mkdir(exist_ok=True)
    print("Building synthetic dataset...")
    df = build_dataset()
    output_path = "data/fragment_orders.csv"
    df.to_csv(output_path, index=False)
    print(f"\n✅ Dataset saved to {output_path}")
    print(f"   Total orders: {len(df)}")
    print(f"   Benign:              {len(df[df.label == 0])}")
    print(f"   Suspicious (clean):  {len(df[df.label == 1])}")
    print(f"   Suspicious (mutated):{len(df[df.label == 2])}")
    print(f"   Avg fragments/order: {df.n_fragments.mean():.1f}")
    print(f"   Avg total bp/order:  {df.total_bp.mean():.0f}")
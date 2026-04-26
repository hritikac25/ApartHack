"""
hard_dataset_generator.py
Generates a harder, more realistic synthetic dataset where:
  - Benign sequences share some motifs with the dangerous sequence
  - Suspicious mutated fragments have higher mutation rates (8-12%)
  - Partial orders only cover 40-70% of the reference (ambiguous cases)
  - Some benign orders have many fragments (stress-testing false positive rate)
"""

import random
import json
import pandas as pd
from pathlib import Path

random.seed(99)

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

# Benign sequences that SHARE some short motifs with dangerous sequence
# making them harder to distinguish purely by k-mer overlap
HARD_BENIGN_SEQUENCES = [
    # Shares some codons with dangerous seq but is a real housekeeping-style gene
    "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATCCTCGTCGACCTCGACGGCGACGTG"
    "AACGGCCACAAGTTCAGCGTCAGCGGCGAGGGCGAGGGCGACGCCACCTACGGCAAGCTGACCCTGAAG",
    "AATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGCAAGCTGCCCGTGCCC"
    "TGGCCCACCCTCGTGACCACCTTCGGCTACGGCCTGCAGTGCTTCGCCCGCTACCCCGACCACATGAAG",
    "TTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTC"
    "GATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCC",
    "GGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGAC"
    "GCTAGCATGACTGGTGGACAGCAAATGGGTCGGGATCTGTACGACGATGACGATAAGGATCCGGCTGCT",
    "ACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAA"
    "AACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGAC",
    # Pure random-ish benign sequences
    "TTCGAGCAAGAGATCGAGCACAGTGGCGGCCGCTCGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTC"
    "CCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGC",
    "GATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCC"
    "ACCCTCGTGACCACCTTCGGCTACGGCCTGCAGTGCTTCGCCCGCTACCCCGACCACATGAAGCAGCAC",
    "GACTTCTTCAAGAGCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGC"
    "AACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGC",
]

PURPOSES = [
    "GFP expression vector", "CRISPR guide RNA", "Vaccine antigen cloning",
    "Protein expression tag", "PCR primer design", "Reporter gene construct",
    "Synthetic gene study", "Academic research", "Protein engineering",
    "Gene therapy vector", "Metabolic pathway study", "Enzyme characterization",
    "Antibody humanization", "Cell line development", "Diagnostic probe",
]

def mutate(sequence, rate):
    """Apply random point mutations at a given rate."""
    seq = list(sequence)
    for i in range(len(seq)):
        if random.random() < rate:
            seq[i] = random.choice("ATCG")
    return "".join(seq)

def generate_fragment(sequence, start, length):
    end = min(start + length, len(sequence))
    return sequence[start:end]

def generate_hard_benign_order():
    """
    Benign order using sequences that share motifs with dangerous seq.
    Some orders have many fragments to stress-test false positive rate.
    """
    n_fragments = random.choices(
        [1, 2, 3, 4, 5, 6, 7, 8],
        weights=[10, 20, 25, 20, 10, 7, 5, 3]
    )[0]
    fragments = []
    for _ in range(n_fragments):
        seq = random.choice(HARD_BENIGN_SEQUENCES)
        start = random.randint(0, max(0, len(seq) - 50))
        length = random.randint(30, 100)
        frag = generate_fragment(seq, start, length)
        # Occasionally mutate benign fragments too (natural variation)
        if random.random() < 0.3:
            frag = mutate(frag, rate=0.05)
        fragments.append({
            "sequence": frag,
            "length": len(frag),
            "claimed_purpose": random.choice(PURPOSES)
        })
    return fragments

def generate_hard_suspicious_order(order_type="high"):
    """
    Three tiers of suspicious orders:
      high:    85-98% coverage, low mutation (0-2%) — should be caught
      medium:  55-84% coverage, medium mutation (5-8%) — harder
      partial: 30-54% coverage, high mutation (9-13%) — ambiguous
    """
    seq = FICTIONAL_DANGEROUS_SEQUENCE

    if order_type == "high":
        coverage = random.uniform(0.85, 0.98)
        mutation_rate = random.uniform(0.00, 0.02)
    elif order_type == "medium":
        coverage = random.uniform(0.55, 0.84)
        mutation_rate = random.uniform(0.05, 0.08)
    else:  # partial
        coverage = random.uniform(0.30, 0.54)
        mutation_rate = random.uniform(0.09, 0.13)

    target_len = int(len(seq) * coverage)
    fragments = []
    pos = 0
    while pos < target_len:
        length = random.randint(60, 130)
        overlap = random.randint(8, 25)
        raw = generate_fragment(seq, pos, length)
        frag = mutate(raw, mutation_rate)
        fragments.append({
            "sequence": frag,
            "length": len(frag),
            "claimed_purpose": random.choice(PURPOSES)
        })
        pos += length - overlap

    return fragments, order_type

def build_hard_dataset(
    n_benign=200,
    n_high=60,
    n_medium=60,
    n_partial=60
):
    records = []

    print(f"Generating {n_benign} hard benign orders...")
    for i in range(n_benign):
        frags = generate_hard_benign_order()
        records.append({
            "order_id": f"BEN-{i:04d}",
            "label": 0,
            "label_name": "benign",
            "n_fragments": len(frags),
            "total_bp": sum(f["length"] for f in frags),
            "fragments": json.dumps([f["sequence"] for f in frags]),
            "purposes": json.dumps([f["claimed_purpose"] for f in frags]),
        })

    print(f"Generating {n_high} high-coverage suspicious orders...")
    for i in range(n_high):
        frags, otype = generate_hard_suspicious_order("high")
        records.append({
            "order_id": f"HIG-{i:04d}",
            "label": 1,
            "label_name": "suspicious_high",
            "n_fragments": len(frags),
            "total_bp": sum(f["length"] for f in frags),
            "fragments": json.dumps([f["sequence"] for f in frags]),
            "purposes": json.dumps([f["claimed_purpose"] for f in frags]),
        })

    print(f"Generating {n_medium} medium-coverage suspicious orders...")
    for i in range(n_medium):
        frags, otype = generate_hard_suspicious_order("medium")
        records.append({
            "order_id": f"MED-{i:04d}",
            "label": 2,
            "label_name": "suspicious_medium",
            "n_fragments": len(frags),
            "total_bp": sum(f["length"] for f in frags),
            "fragments": json.dumps([f["sequence"] for f in frags]),
            "purposes": json.dumps([f["claimed_purpose"] for f in frags]),
        })

    print(f"Generating {n_partial} partial/ambiguous suspicious orders...")
    for i in range(n_partial):
        frags, otype = generate_hard_suspicious_order("partial")
        records.append({
            "order_id": f"PAR-{i:04d}",
            "label": 3,
            "label_name": "suspicious_partial",
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
    print("Building hard dataset...")
    df = build_hard_dataset()
    output_path = "data/hard_fragment_orders.csv"
    df.to_csv(output_path, index=False)
    print(f"\n✅ Dataset saved to {output_path}")
    print(f"   Total orders     : {len(df)}")
    print(f"   Benign           : {len(df[df.label == 0])}")
    print(f"   Suspicious high  : {len(df[df.label == 1])}")
    print(f"   Suspicious medium: {len(df[df.label == 2])}")
    print(f"   Suspicious partial: {len(df[df.label == 3])}")
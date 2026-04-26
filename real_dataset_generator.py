"""
real_dataset_generator.py

Generates a dataset using REAL publicly available H1N1 reference sequences
fetched from NCBI GenBank. Simulates split-order attacks against each
of the three segments independently, then cross-segment combinations.
"""

import random
import json
import pandas as pd
from pathlib import Path
from Bio import SeqIO

random.seed(42)

PURPOSES = [
    "GFP expression vector", "CRISPR guide RNA", "Vaccine antigen cloning",
    "Protein expression tag", "PCR primer design", "Reporter gene construct",
    "Synthetic gene study", "Academic research", "Protein engineering",
    "Gene therapy vector", "Metabolic pathway study", "Enzyme characterization",
]

def load_real_references():
    """Load the real sequences fetched from NCBI GenBank."""
    refs = {}
    ref_dir = Path("data/references")
    for fasta_file in ref_dir.glob("*.fasta"):
        name = fasta_file.stem
        record = SeqIO.read(fasta_file, "fasta")
        refs[name] = {
            "sequence": str(record.seq).upper(),
            "accession": record.id,
            "description": record.description,
            "length": len(record.seq)
        }
        print(f"  Loaded {name}: {len(record.seq)} bp ({record.id})")
    return refs

def mutate(sequence, rate):
    seq = list(sequence)
    for i in range(len(seq)):
        if random.random() < rate:
            seq[i] = random.choice("ATCG")
    return "".join(seq)

def generate_fragment(sequence, start, length):
    end = min(start + length, len(sequence))
    return sequence[start:end]

def generate_benign_order(refs):
    """
    Benign order: random short fragments from non-target regions,
    or scrambled sequences that don't assemble into anything meaningful.
    """
    n_fragments = random.randint(1, 5)
    fragments = []

    # Pick a random reference but take only tiny non-overlapping pieces
    # from random positions — not enough to assemble
    ref = random.choice(list(refs.values()))
    seq = ref["sequence"]

    used_positions = set()
    for _ in range(n_fragments):
        # Pick a random short fragment from scattered positions
        while True:
            start = random.randint(0, max(0, len(seq) - 40))
            if start not in used_positions:
                used_positions.add(start)
                break
        length = random.randint(20, 50)
        frag = generate_fragment(seq, start, length)
        # Apply high mutation to make it non-assembling
        frag = mutate(frag, rate=random.uniform(0.15, 0.30))
        fragments.append({
            "sequence": frag,
            "length": len(frag),
            "claimed_purpose": random.choice(PURPOSES),
            "source_ref": "benign_derived"
        })
    return fragments

def generate_split_order(ref_name, ref_data, order_type="high"):
    """
    Generate a split-order attack against a real reference segment.
    order_type controls coverage and mutation rate.
    """
    seq = ref_data["sequence"]

    if order_type == "high":
        coverage = random.uniform(0.80, 0.95)
        mutation_rate = random.uniform(0.00, 0.02)
    elif order_type == "medium":
        coverage = random.uniform(0.50, 0.79)
        mutation_rate = random.uniform(0.04, 0.08)
    else:  # partial
        coverage = random.uniform(0.25, 0.49)
        mutation_rate = random.uniform(0.08, 0.14)

    target_len = int(len(seq) * coverage)
    fragments = []
    pos = 0

    while pos < target_len:
        length = random.randint(60, 150)
        overlap = random.randint(10, 30)
        raw = generate_fragment(seq, pos, length)
        frag = mutate(raw, mutation_rate)
        fragments.append({
            "sequence": frag,
            "length": len(frag),
            "claimed_purpose": random.choice(PURPOSES),
            "source_ref": ref_name,
            "target_accession": ref_data["accession"]
        })
        pos += length - overlap

    return fragments, coverage, mutation_rate

def build_real_dataset(refs, n_benign=150, n_per_type_per_ref=30):
    """
    Build dataset across all three real reference segments.
    Each segment gets its own split-order attack scenarios.
    """
    records = []

    print(f"\nGenerating {n_benign} benign orders...")
    for i in range(n_benign):
        frags = generate_benign_order(refs)
        records.append({
            "order_id": f"BEN-{i:04d}",
            "label": 0,
            "label_name": "benign",
            "target_segment": "none",
            "target_accession": "none",
            "n_fragments": len(frags),
            "total_bp": sum(f["length"] for f in frags),
            "fragments": json.dumps([f["sequence"] for f in frags]),
            "purposes": json.dumps([f["claimed_purpose"] for f in frags]),
            "coverage": 0.0,
            "mutation_rate": 0.0
        })

    label_counter = 1
    for ref_name, ref_data in refs.items():
        for order_type, label_name in [
            ("high", "suspicious_high"),
            ("medium", "suspicious_medium"),
            ("partial", "suspicious_partial")
        ]:
            print(f"Generating {n_per_type_per_ref} {label_name} orders "
                  f"for {ref_name} ({ref_data['accession']})...")
            for i in range(n_per_type_per_ref):
                frags, cov, mut = generate_split_order(
                    ref_name, ref_data, order_type
                )
                records.append({
                    "order_id": f"{ref_name[:3].upper()}-{order_type[:3].upper()}-{i:04d}",
                    "label": label_counter,
                    "label_name": label_name,
                    "target_segment": ref_name,
                    "target_accession": ref_data["accession"],
                    "n_fragments": len(frags),
                    "total_bp": sum(f["length"] for f in frags),
                    "fragments": json.dumps([f["sequence"] for f in frags]),
                    "purposes": json.dumps([f["claimed_purpose"] for f in frags]),
                    "coverage": round(cov, 3),
                    "mutation_rate": round(mut, 3)
                })
            label_counter += 1

    random.shuffle(records)
    df = pd.DataFrame(records)
    return df

if __name__ == "__main__":
    Path("data").mkdir(exist_ok=True)
    print("Loading real NCBI GenBank reference sequences...")
    refs = load_real_references()

    print(f"\nLoaded {len(refs)} real reference segments:")
    for name, data in refs.items():
        print(f"  {data['accession']}: {data['description'][:55]}...")

    df = build_real_dataset(refs)

    output_path = "data/real_fragment_orders.csv"
    df.to_csv(output_path, index=False)

    print(f"\n✅ Real dataset saved to {output_path}")
    print(f"   Total orders     : {len(df)}")
    print(f"   Benign           : {len(df[df.label == 0])}")
    print(f"   Suspicious orders: {len(df[df.label > 0])}")
    print(f"   Reference segments tested: {len(refs)}")
    for name, data in refs.items():
        print(f"     {data['accession']} — {name} ({data['length']} bp)")
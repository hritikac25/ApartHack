"""
ibbis_comparison.py

Implements a window-based screening baseline modeled on IBBIS commec's
published methodology, then compares it against FARS on the same dataset.

IBBIS commec screens individual sequences using HMM-based comparison
against a curated database. Its key documented behaviors are:
  - Best performance above 150bp fragments
  - Evaluates each sequence independently (no cross-order assembly)
  - Window-based: screens subsequences of a fixed length

We model this as: for each fragment in an order, check whether any
150bp window of the fragment matches the reference above a similarity
threshold. If any fragment triggers, flag the order.

This is explicitly a SIMPLIFIED MODEL of commec's logic for comparison
purposes. Real commec uses HMMs against a curated database; we use
k-mer similarity against our GenBank references. The key structural
difference we are testing is individual-sequence vs. assembly-aware scoring.
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO
import glob

# ─────────────────────────────────────────────
# Load real reference sequences
# ─────────────────────────────────────────────

def load_references():
    refs = {}
    for fasta_file in glob.glob("data/references/*.fasta"):
        record = SeqIO.read(fasta_file, "fasta")
        refs[record.id] = str(record.seq).upper()
    return refs

REAL_REFERENCES = load_references()

# ─────────────────────────────────────────────
# IBBIS-style window screener
# ─────────────────────────────────────────────

def kmer_similarity(seq_a, seq_b, k=8):
    """Jaccard similarity of k-mer sets between two sequences."""
    def kmers(s):
        return set(s[i:i+k] for i in range(len(s)-k+1))
    a, b = kmers(seq_a), kmers(seq_b)
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)

def ibbis_screen_fragment(fragment, references, window_size=150, step=50, threshold=0.35):
    """
    Screen a single fragment using window-based approach.
    Models IBBIS commec's key behavior: evaluates each sequence
    independently using sliding windows.

    Returns (flagged, max_similarity, best_ref)
    """
    # IBBIS performs best above 150bp — below this, sensitivity drops
    if len(fragment) < 150:
        # Short fragment: reduced sensitivity, model as lower threshold
        effective_threshold = threshold * 1.4
    else:
        effective_threshold = threshold

    max_sim = 0.0
    best_ref = None

    for acc, ref_seq in references.items():
        # Slide window across fragment
        for start in range(0, max(1, len(fragment) - window_size + 1), step):
            window = fragment[start:start + window_size]
            # Compare window against reference windows
            for ref_start in range(0, len(ref_seq) - len(window) + 1, step):
                ref_window = ref_seq[ref_start:ref_start + len(window)]
                sim = kmer_similarity(window, ref_window)
                if sim > max_sim:
                    max_sim = sim
                    best_ref = acc

    flagged = max_sim >= effective_threshold
    return flagged, max_sim, best_ref

def ibbis_screen_order(fragments, references):
    """
    IBBIS-style screening: flag order if ANY fragment triggers.
    This is the key structural difference from FARS:
    no assembly-awareness, no cross-fragment aggregation.
    """
    any_flagged = False
    max_sim_overall = 0.0
    fragment_results = []

    for frag in fragments:
        flagged, max_sim, best_ref = ibbis_screen_fragment(frag, references)
        fragment_results.append({
            "length": len(frag),
            "flagged": flagged,
            "max_similarity": max_sim,
            "best_ref": best_ref
        })
        if flagged:
            any_flagged = True
        if max_sim > max_sim_overall:
            max_sim_overall = max_sim

    # IBBIS flags order if any fragment triggers
    flag = "ALERT" if any_flagged else "CLEAR"
    return {
        "ibbis_flag": flag,
        "ibbis_max_similarity": round(max_sim_overall, 4),
        "ibbis_fragments_flagged": sum(1 for f in fragment_results if f["flagged"]),
        "n_fragments": len(fragments)
    }

# ─────────────────────────────────────────────
# Run comparison
# ─────────────────────────────────────────────

def run_comparison():
    orders = pd.read_csv("data/real_fragment_orders.csv")
    fars_results = pd.read_csv("data/scored_orders.csv")

    print(f"Running IBBIS-style screening on {len(orders)} orders...")
    print("(This may take a few minutes on real genomic sequences)\n")

    ibbis_results = []
    for i, (_, row) in enumerate(orders.iterrows()):
        if i % 50 == 0:
            print(f"  Processing order {i}/{len(orders)}...")
        fragments = json.loads(row["fragments"])
        result = ibbis_screen_order(fragments, REAL_REFERENCES)
        result["order_id"] = row["order_id"]
        result["label"] = row["label"]
        result["label_name"] = row["label_name"]
        result["target_segment"] = row["target_segment"]
        result["coverage"] = row["coverage"]
        result["mutation_rate"] = row["mutation_rate"]
        ibbis_results.append(result)

    ibbis_df = pd.DataFrame(ibbis_results)
    ibbis_df.to_csv("data/ibbis_results.csv", index=False)

    # Merge with FARS results
    merged = ibbis_df.merge(
        fars_results[["order_id", "risk_score", "flag"]],
        on="order_id"
    )
    merged.columns = [c if c not in ["flag"] else "fars_flag"
                      for c in merged.columns]

    return merged

def print_comparison(merged):
    benign = merged[merged.label == 0]
    suspicious = merged[merged.label > 0]

    tier_map = {
        "suspicious_high":    "High (80-95%)",
        "suspicious_medium":  "Medium (50-79%)",
        "suspicious_partial": "Partial (25-49%)",
    }

    print("=" * 72)
    print("  TABLE: FARS vs IBBIS-Style Window Screening — Head to Head")
    print("=" * 72)
    print(f"\n  {'Tier':<22} {'N':>5} | "
          f"{'FARS Flagged':>13} | {'IBBIS Flagged':>13} | {'FARS Advantage':>14}")
    print("-" * 72)

    # Benign
    n = len(benign)
    fars_fp = (benign.fars_flag.isin(["ALERT","REVIEW"])).sum()
    ibbis_fp = (benign.ibbis_flag == "ALERT").sum()
    print(f"  {'Benign':<22} {n:>5} | "
          f"{f'{fars_fp} FP ({100*fars_fp/n:.1f}%)':>13} | "
          f"{f'{ibbis_fp} FP ({100*ibbis_fp/n:.1f}%)':>13} | "
          f"{'—':>14}")

    total_fars_advantage = 0
    for label_name, tier_label in tier_map.items():
        subset = suspicious[suspicious.label_name == label_name]
        if len(subset) == 0:
            continue
        n = len(subset)
        fars_caught = (subset.fars_flag.isin(["ALERT","REVIEW"])).sum()
        ibbis_caught = (subset.ibbis_flag == "ALERT").sum()
        advantage = fars_caught - ibbis_caught
        total_fars_advantage += advantage
        print(f"  {tier_label:<22} {n:>5} | "
              f"{f'{fars_caught} ({100*fars_caught/n:.1f}%)':>13} | "
              f"{f'{ibbis_caught} ({100*ibbis_caught/n:.1f}%)':>13} | "
              f"{f'+{advantage} orders':>14}")

    print("-" * 72)
    n_sus = len(suspicious)
    fars_total = (suspicious.fars_flag.isin(["ALERT","REVIEW"])).sum()
    ibbis_total = (suspicious.ibbis_flag == "ALERT").sum()
    print(f"  {'TOTAL SUSPICIOUS':<22} {n_sus:>5} | "
          f"{f'{fars_total} ({100*fars_total/n_sus:.1f}%)':>13} | "
          f"{f'{ibbis_total} ({100*ibbis_total/n_sus:.1f}%)':>13} | "
          f"{f'+{total_fars_advantage} orders':>14}")

    print("\n" + "=" * 72)
    print("  KEY FINDING: Where assembly-awareness matters most")
    print("=" * 72)

    partial = suspicious[suspicious.label_name == "suspicious_partial"]
    fars_partial = (partial.fars_flag.isin(["ALERT","REVIEW"])).sum()
    ibbis_partial = (partial.ibbis_flag == "ALERT").sum()

    print(f"""
  On partial split orders (25-49% coverage, high mutation):
    FARS detection rate  : {100*fars_partial/len(partial):.1f}% ({fars_partial}/{len(partial)})
    IBBIS-style detection: {100*ibbis_partial/len(partial):.1f}% ({ibbis_partial}/{len(partial)})
    Improvement          : +{fars_partial - ibbis_partial} orders caught

  This gap exists because IBBIS-style screening evaluates each fragment
  independently. Short, heavily mutated fragments (avg 3-4 per partial order)
  fall below the individual detection threshold. FARS aggregates coverage
  and contiguity signals across all fragments, catching assembly patterns
  that individual screening misses.

  False positive rate:
    FARS        : {100*(benign.fars_flag.isin(['ALERT','REVIEW'])).sum()/len(benign):.1f}%
    IBBIS-style : {100*(benign.ibbis_flag=='ALERT').sum()/len(benign):.1f}%
    """)

if __name__ == "__main__":
    merged = run_comparison()
    print_comparison(merged)
    print("✅ Comparison saved to data/ibbis_results.csv")
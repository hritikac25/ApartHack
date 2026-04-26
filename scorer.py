"""
scorer.py
Fragment Assembly Risk Scorer — core detection engine.

For each order, computes:
  - Coverage: what % of the reference sequence do the fragments cover?
  - Overlap graph: do fragments connect into a contiguous assembly?
  - Risk score: 0.0 (safe) to 1.0 (high risk)
  - Flag: CLEAR, REVIEW, or ALERT
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from colorama import init, Fore, Style

init(autoreset=True)

# ──────────────────────────────────────────────────────────
# REFERENCE DATABASE — Real NCBI GenBank sequences
# Loaded dynamically from data/references/
# ──────────────────────────────────────────────────────────

from Bio import SeqIO
import glob

def load_real_references():
    refs = {}
    for fasta_file in glob.glob("data/references/*.fasta"):
        record = SeqIO.read(fasta_file, "fasta")
        refs[record.id] = str(record.seq).upper()
    return refs

REAL_REFERENCES = load_real_references()

# Keep fictional reference for backward compatibility
REFERENCE_SEQUENCE = (
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
REF_LEN = len(REFERENCE_SEQUENCE)

# ──────────────────────────────────────────────────────────
# ALIGNMENT: find best match position of a fragment
# Uses a sliding window similarity score (k-mer based)
# ──────────────────────────────────────────────────────────

def kmer_set(seq, k=10):
    """Return the set of all k-mers in a sequence."""
    return set(seq[i:i+k] for i in range(len(seq) - k + 1))

def align_fragment(fragment, reference, k=8, step=3):
    """
    Slide the fragment across the reference in steps.
    At each position, compute Jaccard similarity of k-mer sets.
    Uses smaller k and step for real genomic sequences.
    Returns (best_start, best_end, best_score).
    """
    frag_kmers = kmer_set(fragment, k)
    if not frag_kmers:
        return None, None, 0.0

    best_score = 0.0
    best_start = 0
    frag_len = len(fragment)

    for start in range(0, len(reference) - frag_len + 1, step):
        window = reference[start:start + frag_len]
        window_kmers = kmer_set(window, k)
        union = frag_kmers | window_kmers
        if not union:
            continue
        score = len(frag_kmers & window_kmers) / len(union)
        if score > best_score:
            best_score = score
            best_start = start

    best_end = best_start + frag_len
    return best_start, best_end, best_score

# ──────────────────────────────────────────────────────────
# COVERAGE: given aligned fragments, what % of ref is covered?
# ──────────────────────────────────────────────────────────

def compute_coverage(alignments, ref_len):
    """
    alignments: list of (start, end, score) tuples
    Returns coverage as a float 0.0–1.0
    """
    covered = np.zeros(ref_len, dtype=bool)
    for start, end, score in alignments:
        if start is not None and score > 0.15:  # minimum alignment quality
            covered[start:end] = True
    return covered.sum() / ref_len

# ──────────────────────────────────────────────────────────
# CONTIGUITY: are the fragments connected into one assembly?
# ──────────────────────────────────────────────────────────

def compute_contiguity(alignments, ref_len, gap_tolerance=50):
    """
    Check whether aligned fragments form a contiguous block
    (with allowed gaps <= gap_tolerance bp).
    Returns contiguity score 0.0–1.0
    """
    valid = sorted(
        [(s, e) for s, e, sc in alignments if s is not None and sc > 0.25],
        key=lambda x: x[0]
    )
    if not valid:
        return 0.0

    # Merge overlapping/adjacent intervals
    merged = [valid[0]]
    for start, end in valid[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + gap_tolerance:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))

    # Largest contiguous block as fraction of ref
    largest = max(e - s for s, e in merged)
    return largest / ref_len

# ──────────────────────────────────────────────────────────
# RISK SCORE
# ──────────────────────────────────────────────────────────

def compute_risk_score(coverage, contiguity, n_fragments):
    """
    Combine signals into a single risk score 0.0–1.0.
    Weights tuned to reward both breadth (coverage) and
    assembly potential (contiguity).
    """
    # Penalize orders with suspiciously many fragments
    frag_penalty = min(1.0, n_fragments / 10.0)

    score = (
        0.50 * coverage +
        0.35 * contiguity +
        0.15 * frag_penalty
    )
    return round(min(score, 1.0), 4)

def flag_order(risk_score):
    if risk_score >= 0.60:
        return "ALERT"
    elif risk_score >= 0.35:
        return "REVIEW"
    else:
        return "CLEAR"

# ──────────────────────────────────────────────────────────
# SCORE A SINGLE ORDER
# ──────────────────────────────────────────────────────────

def score_order(fragments):
    """
    Score an order against ALL real reference sequences.
    Returns the highest risk score across all references.
    This reflects real screening: flag if dangerous against ANY known sequence.
    """
    best_result = None
    best_score = -1

    for accession, ref_seq in REAL_REFERENCES.items():
        ref_len = len(ref_seq)
        alignments = []

        for frag in fragments:
            start, end, sim = align_fragment(frag, ref_seq)
            alignments.append((start, end, sim))

        coverage = compute_coverage(alignments, ref_len)
        contiguity = compute_contiguity(alignments, ref_len)
        risk_score = compute_risk_score(coverage, contiguity, len(fragments))
        flag = flag_order(risk_score)

        if risk_score > best_score:
            best_score = risk_score
            best_result = {
                "coverage": round(coverage, 4),
                "contiguity": round(contiguity, 4),
                "risk_score": risk_score,
                "flag": flag,
                "n_fragments": len(fragments),
                "best_match_accession": accession
            }

    return best_result

# ──────────────────────────────────────────────────────────
# SCORE FULL DATASET
# ──────────────────────────────────────────────────────────

def score_dataset(csv_path="data/real_fragment_orders.csv"):
    df = pd.read_csv(csv_path)
    results = []

    print(f"Scoring {len(df)} orders...\n")

    for _, row in df.iterrows():
        fragments = json.loads(row["fragments"])
        result = score_order(fragments)
        result["order_id"] = row["order_id"]
        result["true_label"] = row["label"]
        result["true_label_name"] = row["label_name"]
        results.append(result)

    results_df = pd.DataFrame(results)
    results_df.to_csv("data/scored_orders.csv", index=False)
    return results_df

# ──────────────────────────────────────────────────────────
# PRETTY PRINT SUMMARY
# ──────────────────────────────────────────────────────────

def print_summary(results_df):
    df = pd.read_csv("data/real_fragment_orders.csv")
    print("=" * 60)
    print("  FRAGMENT ASSEMBLY RISK SCORER — RESULTS SUMMARY")

    label_names = df[['label','label_name']].drop_duplicates().sort_values('label')
    for _, row in label_names.iterrows():
        label_id = row['label']
        label_name = row['label_name']
        subset = results_df[results_df.true_label == label_id]
        alerts = (subset.flag == "ALERT").sum()
        reviews = (subset.flag == "REVIEW").sum()
        clears = (subset.flag == "CLEAR").sum()
        avg_score = subset.risk_score.mean()

        print(f"\n── {label_name.upper()} (n={len(subset)}) ──")
        print(f"   Avg risk score : {avg_score:.3f}")
        print(f"   🔴 ALERT       : {alerts}  ({100*alerts/len(subset):.1f}%)")
        print(f"   🟡 REVIEW      : {reviews} ({100*reviews/len(subset):.1f}%)")
        print(f"   🟢 CLEAR       : {clears}  ({100*clears/len(subset):.1f}%)")

    print("\n" + "=" * 60)
    print("  DETECTION PERFORMANCE")
    print("=" * 60)

    # Treat ALERT+REVIEW as "flagged"
    suspicious = results_df[results_df.true_label > 0]
    benign = results_df[results_df.true_label == 0]

    true_positives = (suspicious.flag.isin(["ALERT", "REVIEW"])).sum()
    false_negatives = (suspicious.flag == "CLEAR").sum()
    true_negatives = (benign.flag == "CLEAR").sum()
    false_positives = (benign.flag.isin(["ALERT", "REVIEW"])).sum()

    sensitivity = true_positives / len(suspicious)
    specificity = true_negatives / len(benign)
    precision = true_positives / max(1, true_positives + false_positives)

    print(f"\n   True Positives  (caught suspicious) : {true_positives}/{len(suspicious)}")
    print(f"   False Negatives (missed suspicious) : {false_negatives}/{len(suspicious)}")
    print(f"   True Negatives  (cleared benign)    : {true_negatives}/{len(benign)}")
    print(f"   False Positives (flagged benign)    : {false_positives}/{len(benign)}")
    print(f"\n   Sensitivity (recall) : {sensitivity:.1%}")
    print(f"   Specificity          : {specificity:.1%}")
    print(f"   Precision            : {precision:.1%}")
    print("=" * 60)

if __name__ == "__main__":
    results_df = score_dataset()
    print_summary(results_df)
    print(f"\n✅ Full results saved to data/scored_orders.csv")
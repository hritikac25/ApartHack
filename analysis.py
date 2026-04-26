"""
analysis.py
Deeper analysis of scorer results
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path

def load_results():
    df = pd.read_csv("data/scored_orders.csv")
    return df

def per_category_summary(df):
    print("\n" + "="*65)
    print("  TABLE 1: Per-Category Detection Summary")
    print("="*65)
    print(f"{'Category':<25} {'N':>5} {'Avg Score':>10} {'ALERT%':>8} {'REVIEW%':>8} {'CLEAR%':>8}")
    print("-"*65)

    for label_id, label_name in [
        (0, "Benign"),
        (1, "Suspicious-High"),
        (2, "Suspicious-Medium"),
        (3, "Suspicious-Partial")
    ]:
        subset = df[df.true_label == label_id]
        n = len(subset)
        avg = subset.risk_score.mean()
        alert_pct = 100 * (subset.flag == "ALERT").sum() / n
        review_pct = 100 * (subset.flag == "REVIEW").sum() / n
        clear_pct = 100 * (subset.flag == "CLEAR").sum() / n
        print(f"{label_name:<25} {n:>5} {avg:>10.3f} {alert_pct:>7.1f}% {review_pct:>7.1f}% {clear_pct:>7.1f}%")

def threshold_analysis(df):
    """
    Show how sensitivity and specificity change as we vary
    the ALERT threshold from 0.1 to 0.9.
    This finds the optimal operating point.
    """
    print("\n" + "="*65)
    print("  TABLE 2: Sensitivity vs Specificity at Different Thresholds")
    print("="*65)
    print(f"{'Threshold':>10} {'Sensitivity':>13} {'Specificity':>13} {'Precision':>11} {'F1':>8}")
    print("-"*65)

    benign = df[df.true_label == 0]
    suspicious_all = df[df.true_label.isin([1, 2, 3])]

    best_f1 = 0
    best_threshold = 0

    for threshold in np.arange(0.10, 0.95, 0.05):
        tp = (suspicious_all.risk_score >= threshold).sum()
        fn = (suspicious_all.risk_score < threshold).sum()
        tn = (benign.risk_score < threshold).sum()
        fp = (benign.risk_score >= threshold).sum()

        sensitivity = tp / len(suspicious_all)
        specificity = tn / len(benign)
        precision = tp / max(1, tp + fp)
        f1 = 2 * precision * sensitivity / max(0.001, precision + sensitivity)

        if f1 > best_f1:
            best_f1 = f1
            best_threshold = threshold

        marker = " <-- optimal" if abs(threshold - best_threshold) < 0.001 else ""
        print(f"{threshold:>10.2f} {sensitivity:>12.1%} {specificity:>12.1%} {precision:>10.1%} {f1:>8.3f}{marker}")

    return best_threshold

def coverage_vs_detection(df):
    """
    For suspicious orders only: show how coverage correlates
    with detection rate. Key finding for the paper.
    """
    print("\n" + "="*65)
    print("  TABLE 3: Coverage Level vs Detection Rate (Suspicious Only)")
    print("="*65)
    print(f"{'Coverage Tier':<25} {'Avg Risk Score':>15} {'Detection Rate':>15}")
    print("-"*65)

    sus = df[df.true_label.isin([1, 2, 3])].copy()

    tier_map = {1: "High (85-98%)", 2: "Medium (55-84%)", 3: "Partial (30-54%)"}
    for label_id, tier_name in tier_map.items():
        subset = sus[sus.true_label == label_id]
        avg_score = subset.risk_score.mean()
        detected = (subset.flag.isin(["ALERT", "REVIEW"])).sum()
        detection_rate = detected / len(subset)
        print(f"{tier_name:<25} {avg_score:>15.3f} {detection_rate:>14.1%}")

def false_negative_analysis(df):
    """
    Examine the orders that were missed (false negatives).
    What do they have in common?
    """
    print("\n" + "="*65)
    print("  TABLE 4: Missed Suspicious Orders (False Negatives)")
    print("="*65)

    suspicious = df[df.true_label.isin([1, 2, 3])]
    missed = suspicious[suspicious.flag == "CLEAR"]

    print(f"  Total suspicious orders : {len(suspicious)}")
    print(f"  Missed (CLEAR)          : {len(missed)}")
    print(f"  Overall miss rate       : {100*len(missed)/len(suspicious):.1f}%")

    if len(missed) > 0:
        print(f"\n  Missed orders by category:")
        for label_id, name in [(1,"High"),(2,"Medium"),(3,"Partial")]:
            n_missed = len(missed[missed.true_label == label_id])
            n_total = len(suspicious[suspicious.true_label == label_id])
            print(f"    {name:<12}: {n_missed}/{n_total} missed ({100*n_missed/n_total:.1f}%)")

        print(f"\n  Avg risk score of missed orders : {missed.risk_score.mean():.3f}")
        print(f"  Avg fragments in missed orders  : {missed.n_fragments.mean():.1f}")

def key_findings(df, best_threshold):
    print("\n" + "="*65)
    print("  KEY FINDINGS FOR REPORT")
    print("="*65)

    benign = df[df.true_label == 0]
    sus_high = df[df.true_label == 1]
    sus_med = df[df.true_label == 2]
    sus_partial = df[df.true_label == 3]
    sus_all = df[df.true_label.isin([1,2,3])]

    tp = (sus_all.risk_score >= best_threshold).sum()
    tn = (benign.risk_score < best_threshold).sum()
    fp = (benign.risk_score >= best_threshold).sum()
    fn = (sus_all.risk_score < best_threshold).sum()

    sensitivity = tp / len(sus_all)
    specificity = tn / len(benign)
    precision = tp / max(1, tp+fp)
    f1 = 2*precision*sensitivity / max(0.001, precision+sensitivity)

    print(f"""
  At optimal threshold ({best_threshold:.2f}):
    Sensitivity    : {sensitivity:.1%}
    Specificity    : {specificity:.1%}
    Precision      : {precision:.1%}
    F1 Score       : {f1:.3f}

  High-coverage split orders (85-98%):
    Detection rate : {100*(sus_high.flag.isin(['ALERT','REVIEW'])).sum()/len(sus_high):.1f}%
    Avg risk score : {sus_high.risk_score.mean():.3f}

  Medium-coverage + mutated (55-84%, 5-8% mutation):
    Detection rate : {100*(sus_med.flag.isin(['ALERT','REVIEW'])).sum()/len(sus_med):.1f}%
    Avg risk score : {sus_med.risk_score.mean():.3f}

  Partial orders (30-54%, 9-13% mutation) — key gap:
    Detection rate : {100*(sus_partial.flag.isin(['ALERT','REVIEW'])).sum()/len(sus_partial):.1f}%
    Avg risk score : {sus_partial.risk_score.mean():.3f}
    → These represent the detection boundary of the current tool.
    → A shared cross-device database would be needed to catch these.

  False positive rate on benign orders:
    {100*fp/len(benign):.1f}% — meaning {fp} benign orders were incorrectly flagged.
    """)

if __name__ == "__main__":
    df = load_results()
    per_category_summary(df)
    best_threshold = threshold_analysis(df)
    coverage_vs_detection(df)
    false_negative_analysis(df)
    key_findings(df, best_threshold)
    print("✅ Analysis complete. Use these tables directly in your report.")
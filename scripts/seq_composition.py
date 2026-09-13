"""Sequence-composition penalty for the RL reward.

Motivation: ProteinMPNN at low sampling temperature can collapse into
low-complexity sequence (long homopolymer runs, one residue dominating).
AlphaFold2 scores such sequences confidently -- alanine is the strongest
helix former, so an Ala-rich bundle folds cleanly in silico -- but those
designs are poor candidates for synthesis. Raw ipTM as an RL reward
therefore rewards the degeneracy instead of penalising it.

Observed case that motivated this: an IL-7Ra design scored ipTM 0.858
(pLDDT 91, interface PAE 5.8) while being 40% alanine with 7-residue
poly-Ala runs.

This module returns a penalty in ipTM units, so the RL reward becomes:

    reward = iptm - composition_penalty(seq)

Thresholds are deliberately set so that ordinary RFdiffusion/ProteinMPNN
helical binders (which are legitimately E/K/A-rich) incur only a small
penalty, while genuinely degenerate sequence is pushed well down.
"""
import math
from collections import Counter

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

# --- thresholds: no penalty at or below these, penalty grows beyond ---
MAX_RUN_FREE = 4        # natural proteins rarely exceed 4 identical residues in a row
MAX_FRAC_FREE = 0.25    # designed helical bundles legitimately reach ~25% for E or K
MIN_ENTROPY_FREE = 3.2  # bits; natural globular protein ~4.1, helical design ~3.3-3.6

# --- weights: chosen so a 40%-Ala sequence with long runs loses ~0.3 ipTM ---
W_RUN = 0.04            # per residue of run length beyond MAX_RUN_FREE
W_FRAC = 1.20           # per unit of single-AA fraction beyond MAX_FRAC_FREE
W_ENTROPY = 0.15        # per bit of entropy below MIN_ENTROPY_FREE

PENALTY_CAP = 0.50      # never subtract more than this, so reward stays informative


def shannon_entropy(seq: str) -> float:
    """Shannon entropy of the residue distribution, in bits."""
    if not seq:
        return 0.0
    counts = Counter(seq)
    n = len(seq)
    return -sum((c / n) * math.log2(c / n) for c in counts.values())


def max_homopolymer_run(seq: str) -> int:
    """Length of the longest run of a single repeated residue."""
    if not seq:
        return 0
    best = run = 1
    for prev, cur in zip(seq, seq[1:]):
        run = run + 1 if cur == prev else 1
        best = max(best, run)
    return best


def max_aa_fraction(seq: str) -> float:
    """Fraction of the sequence taken by its single most common residue."""
    if not seq:
        return 0.0
    return max(Counter(seq).values()) / len(seq)


#              Below this length the composition statistics are not meaningful
#              (a 5-residue peptide is trivially "100% one residue" at worst and
#              cannot reach the entropy threshold at best), so no penalty is
#              applied. Real binder designs are 50-120 aa.
MIN_ASSESSABLE_LENGTH = 20


def composition_penalty(seq: str, detail: bool = False):
    """Penalty in ipTM units for low-complexity / degenerate composition.

    Returns a float, or (float, dict) when detail=True. Sequences shorter than
    MIN_ASSESSABLE_LENGTH get no penalty -- they carry too little signal to
    judge, and silently penalising them would disguise upstream parse errors.
    """
    run = max_homopolymer_run(seq)
    frac = max_aa_fraction(seq)
    ent = shannon_entropy(seq)

    if len(seq) < MIN_ASSESSABLE_LENGTH:
        run_pen = frac_pen = ent_pen = 0.0
    else:
        run_pen = W_RUN * max(0, run - MAX_RUN_FREE)
        frac_pen = W_FRAC * max(0.0, frac - MAX_FRAC_FREE)
        # Entropy is bounded above by log2(len), so for shorter sequences the
        # fixed threshold is unreachable and would penalise them unfairly.
        ent_ceiling = min(MIN_ENTROPY_FREE, math.log2(len(seq)))
        ent_pen = W_ENTROPY * max(0.0, ent_ceiling - ent)

    total = min(PENALTY_CAP, run_pen + frac_pen + ent_pen)

    if not detail:
        return total
    return total, {
        "length": len(seq),
        "assessable": len(seq) >= MIN_ASSESSABLE_LENGTH,
        "max_run": run,
        "max_aa_fraction": frac,
        "dominant_aa": Counter(seq).most_common(1)[0][0] if seq else None,
        "entropy_bits": ent,
        "run_penalty": run_pen,
        "fraction_penalty": frac_pen,
        "entropy_penalty": ent_pen,
        "total_penalty": total,
    }


def shaped_reward(iptm: float, seq: str, detail: bool = False):
    """RL reward: ipTM minus the composition penalty, floored at 0."""
    pen, info = composition_penalty(seq, detail=True)
    reward = max(0.0, iptm - pen)
    if not detail:
        return reward
    info["iptm"] = iptm
    info["shaped_reward"] = reward
    return reward, info

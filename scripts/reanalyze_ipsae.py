#!/usr/bin/env python3
"""Re-run this project's key analyses using ipSAE instead of ipTM.

Motivation (Dunbrack 2025, "Res ipSAE loquuntur: What's wrong with AlphaFold's
ipTM score and how to fix it"): ipTM is computed over whole chains, so unequal
chain lengths shift it spuriously -- which is exactly the size-asymmetry
artifact recorded in this project's own notes, where the true ACE2/RBD complex
scored lower than an off-target pair. ipSAE restricts the calculation to residue
pairs below a pAE threshold, removing that dependence.

Overath et al. 2025 (meta-analysis, 3766 experimentally characterised binders)
found ipSAE_min > 0.61 outperformed ipAE by ~1.4x average precision.

ipSAE is directional (A-B and B-A); ipSAE_min takes the worse direction, which
is the conservative choice the meta-analysis recommends. Every number here comes
from JSON files already on disk -- no new folding.
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def metrics(path):
    d = json.load(open(path))

    def mn(x):
        if isinstance(x, dict):
            v = [float(i) for i in x.values()]
            return min(v) if v else np.nan
        return float(x) if x is not None else np.nan

    plddt = np.array(d.get("plddt", [0.0]))
    return {
        "ipTM": d.get("iptm", np.nan),
        "pTM": d.get("ptm", np.nan),
        "ipSAE_min": mn(d.get("ipsae")),
        "pDockQ2_min": mn(d.get("pdockq2")),
        "pLDDT_all": plddt.mean(),
        "n_res": len(plddt),
    }


def collect(fold_dir, key_from_name=lambda n: n):
    rows = {}
    for f in glob.glob(os.path.join(fold_dir, "*_scores_rank_001_*.json")):
        n = os.path.basename(f).split("_scores_rank_001_")[0]
        rows[key_from_name(n)] = metrics(f)
    return rows


# ---------------------------------------------------------------- analyses

def variance_decomposition(kit):
    print("=" * 74)
    print("1. VARIANCE DECOMPOSITION -- ipTM vs ipSAE_min")
    print("=" * 74)
    print("   (within-backbone share = fraction of variance from WHICH SEQUENCE was drawn)\n")
    out = []
    for t in ("il7ra", "rbd", "pdl1"):
        d = os.path.join(kit, "results", "variance_study", t, "folds")
        man_p = os.path.join(kit, "results", "variance_study", t, "manifest.csv")
        if not os.path.exists(man_p):
            continue
        man = pd.read_csv(man_p)
        keyed = {f"{r['backbone']}__{r['sample']}": r["backbone"] for _, r in man.iterrows()}
        rows = []
        for n, m in collect(d).items():
            if n in keyed:
                rows.append({"backbone": keyed[n], **m})
        if not rows:
            continue
        df = pd.DataFrame(rows)
        rec = {"target": t, "n": len(df)}
        for metric in ("ipTM", "ipSAE_min", "pDockQ2_min"):
            if df[metric].isna().all():
                continue
            g = df.groupby("backbone")[metric]
            between = g.mean().var(ddof=1)
            within = g.var(ddof=1).mean()
            tot = between + within
            rec[metric] = within / tot if tot > 0 else np.nan
        out.append(rec)
    r = pd.DataFrame(out)
    print(r.to_string(index=False, float_format=lambda x: f"{x:.1%}" if x < 1 else f"{x:.0f}"))
    print()


def reward_channel(kit):
    print("=" * 74)
    print("2. CHEAP-SCREEN SELECTION QUALITY -- does the metric pick good sequences?")
    print("=" * 74)
    man = pd.read_csv(os.path.join(kit, "results/variance_study/rbd/manifest.csv"))
    keyed = {f"{r['backbone']}__{r['sample']}": r["backbone"] for _, r in man.iterrows()}
    hi = collect(os.path.join(kit, "results/variance_study/rbd/folds"))

    for label, cheap_dir in [("fast screen", "results/queue/rewardchannel_rbd"),
                             ("hybrid MSA ", "results/queue/hybrid_rbd_folds")]:
        cheap = collect(os.path.join(kit, cheap_dir))
        common = [n for n in cheap if n in hi and n in keyed]
        if not common:
            continue
        print(f"\n  --- selecting by {label} ---")
        for metric in ("ipTM", "ipSAE_min", "pDockQ2_min"):
            picks, rnd, oracle, rhos = [], [], [], []
            for bb in sorted({keyed[n] for n in common}):
                grp = [n for n in common if keyed[n] == bb]
                cv = np.array([cheap[n][metric] for n in grp], dtype=float)
                hv = np.array([hi[n]["ipTM"] for n in grp], dtype=float)
                if np.isnan(cv).any():
                    continue
                picks.append(hv[int(np.nanargmax(cv))])
                rnd.append(hv.mean()); oracle.append(hv.max())
                if len(set(cv)) > 1:
                    rhos.append(spearmanr(cv, hv).correlation)
            if not picks:
                continue
            edge = np.mean(picks) - np.mean(rnd)
            print(f"    by {metric:12s}  picked={np.mean(picks):.3f}  random={np.mean(rnd):.3f}  "
                  f"oracle={np.mean(oracle):.3f}  edge={edge:+.3f}  mean_rho={np.mean(rhos):+.3f}")
    print()


def specificity(kit):
    print("=" * 74)
    print("3. SPECIFICITY RANKING -- ipTM margin vs ipSAE_min margin")
    print("=" * 74)
    on = collect(os.path.join(kit, "results/il7ra_seqrefine/folds"))
    off = {}
    for f in glob.glob(os.path.join(kit, "results/il7ra_decoys/*_scores_rank_001_*.json")):
        base = os.path.basename(f).split("_scores_rank_001_")[0]
        if "--vs--" not in base:
            continue
        cand, decoy = base.split("--vs--")
        off.setdefault(cand, {})[decoy] = metrics(f)

    rows = []
    for cand, dv in off.items():
        if cand not in on:
            continue
        r = {"name": cand}
        for metric in ("ipTM", "ipSAE_min"):
            worst = max(v[metric] for v in dv.values())
            r[f"on_{metric}"] = on[cand][metric]
            r[f"off_{metric}"] = worst
            r[f"margin_{metric}"] = on[cand][metric] - worst
        rows.append(r)
    if not rows:
        print("  (no decoy folds matched)\n"); return
    df = pd.DataFrame(rows)
    a = df.sort_values("margin_ipTM", ascending=False)["name"].tolist()
    b = df.sort_values("margin_ipSAE_min", ascending=False)["name"].tolist()
    print(df.sort_values("margin_ipSAE_min", ascending=False)
            .to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\n  top by ipTM margin      : {a[0]}")
    print(f"  top by ipSAE_min margin : {b[0]}")
    print(f"  ranking agrees: {a == b}")
    print()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kit", default="/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit")
    args = ap.parse_args()
    pd.set_option("display.width", 200)
    variance_decomposition(args.kit)
    reward_channel(args.kit)
    specificity(args.kit)


if __name__ == "__main__":
    main()

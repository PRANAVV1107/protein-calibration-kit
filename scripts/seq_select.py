#!/usr/bin/env python3
"""Select sequences from a multi-sample ProteinMPNN run by composition quality.

Folding is the expensive step (~5 min per high-accuracy fold), sequence design
is cheap (~1.5 s). So generate many sequences per backbone, then spend fold
budget only on the ones whose composition is already acceptable.

Two modes:

  filter  Read ProteinMPNN .fa output, rank each backbone's samples by
          composition penalty, keep the best K, write a FASTA for folding.

  select  Read the resulting fold scores and rank by shaped reward
          (ipTM minus composition penalty), i.e. the same objective the RL
          loop optimises.
"""
import argparse
import glob
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seq_composition import composition_penalty, shaped_reward  # noqa: E402


def read_mpnn_samples(seqs_dir):
    """Parse ProteinMPNN .fa files that contain multiple samples per backbone.

    Layout: line 0 = backbone header, line 1 = original sequence, then
    alternating (sample header, sample sequence) pairs.
    """
    records = []
    for fa in sorted(glob.glob(os.path.join(seqs_dir, "*.fa"))):
        backbone = os.path.splitext(os.path.basename(fa))[0]
        with open(fa) as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        # samples start at index 2, header/sequence alternating
        for i in range(2, len(lines) - 1, 2):
            header, seq = lines[i], lines[i + 1]
            sample_id = "s0"
            for field in header.lstrip(">").split(","):
                if "sample=" in field:
                    sample_id = "s" + field.split("=")[1].strip()
            records.append({"backbone": backbone, "sample": sample_id, "sequence": seq})
    return pd.DataFrame(records)


def cmd_filter(args):
    df = read_mpnn_samples(args.seqs_dir)
    if df.empty:
        sys.exit(f"[FAIL] no sequences parsed from {args.seqs_dir}")

    if args.backbones:
        wanted = {b.strip() for b in args.backbones.split(",") if b.strip()}
        missing = wanted - set(df["backbone"])
        if missing:
            sys.exit(f"[FAIL] requested backbones not found: {sorted(missing)}")
        df = df[df["backbone"].isin(wanted)]
        print(f"Restricted to {len(wanted)} requested backbone(s)")

    detail = df["sequence"].apply(lambda s: composition_penalty(s, detail=True))
    df["penalty"] = [d[0] for d in detail]
    df["max_run"] = [d[1]["max_run"] for d in detail]
    df["max_aa_frac"] = [d[1]["max_aa_fraction"] for d in detail]
    df["dominant_aa"] = [d[1]["dominant_aa"] for d in detail]
    df["length"] = df["sequence"].str.len()

    print(f"Parsed {len(df)} sequences across {df['backbone'].nunique()} backbones")
    print(f"Composition penalty: mean={df['penalty'].mean():.3f}  "
          f"median={df['penalty'].median():.3f}  max={df['penalty'].max():.3f}")
    print(f"Sequences exceeding --max-penalty {args.max_penalty}: "
          f"{(df['penalty'] > args.max_penalty).sum()}/{len(df)}")

    eligible = df[df["penalty"] <= args.max_penalty]
    dropped_backbones = set(df["backbone"]) - set(eligible["backbone"])
    if dropped_backbones:
        print(f"[WARN] {len(dropped_backbones)} backbone(s) had no sequence under the "
              f"penalty threshold and are excluded: {sorted(dropped_backbones)[:5]}")

    # keep the K lowest-penalty samples per backbone
    keep = (eligible.sort_values(["backbone", "penalty"])
                    .groupby("backbone", as_index=False)
                    .head(args.per_backbone))

    target_seq = open(args.target_seq_file).read().strip() if args.target_seq_file else args.target_seq
    if not target_seq:
        sys.exit("[FAIL] need --target-seq or --target-seq-file")

    with open(args.out_fasta, "w") as f:
        for _, r in keep.iterrows():
            name = f"{r['backbone']}__{r['sample']}"
            f.write(f">{name}\n{target_seq}:{r['sequence']}\n")

    keep.to_csv(args.out_manifest, index=False)
    print(f"\n[OK] kept {len(keep)} sequences "
          f"({args.per_backbone} per backbone x {keep['backbone'].nunique()} backbones)")
    print(f"[OK] fold input : {args.out_fasta}")
    print(f"[OK] manifest   : {args.out_manifest}")
    est = len(keep) * args.minutes_per_fold
    print(f"     estimated fold time: ~{est:.0f} min ({est/60:.1f} h) "
          f"at {args.minutes_per_fold} min/fold")


def cmd_select(args):
    manifest = pd.read_csv(args.manifest)
    seq_by_name = {f"{r.backbone}__{r.sample}": r.sequence for r in manifest.itertuples()}

    rows = []
    for f in glob.glob(os.path.join(args.fold_dir, "*_scores_rank_001_*.json")):
        name = os.path.basename(f).split("_scores_rank_001_")[0]
        if name not in seq_by_name:
            continue
        d = json.load(open(f))
        seq = seq_by_name[name]
        plddt = d.get("plddt", [0])
        reward, info = shaped_reward(d["iptm"], seq, detail=True)
        rows.append({
            "name": name,
            "backbone": name.split("__")[0],
            "ipTM": d["iptm"],
            "pTM": d.get("ptm"),
            "mean_pLDDT": sum(plddt) / len(plddt),
            "penalty": info["total_penalty"],
            "shaped_reward": reward,
            "max_run": info["max_run"],
            "max_aa_frac": round(info["max_aa_fraction"], 3),
            "length": len(seq),
        })

    if not rows:
        sys.exit(f"[FAIL] no fold results in {args.fold_dir} matched the manifest")

    df = pd.DataFrame(rows).sort_values("shaped_reward", ascending=False)
    pd.set_option("display.width", 220)
    print(f"\nRanked by shaped reward (ipTM - composition penalty), {len(df)} folds:\n")
    print(df.to_string(index=False))

    best_per_bb = df.sort_values("shaped_reward", ascending=False).groupby("backbone", as_index=False).head(1)
    print(f"\nBest sequence per backbone ({len(best_per_bb)} backbones):\n")
    print(best_per_bb.sort_values("shaped_reward", ascending=False).to_string(index=False))

    df.to_csv(args.out_csv, index=False)
    print(f"\n[OK] full ranking written to {args.out_csv}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("filter", help="rank MPNN samples by composition, emit fold FASTA")
    f.add_argument("--seqs-dir", required=True, help="ProteinMPNN out_folder/seqs directory")
    f.add_argument("--per-backbone", type=int, default=3, help="how many samples to keep per backbone")
    f.add_argument("--backbones", default="", help="comma-separated backbone names to restrict to (default: all)")
    f.add_argument("--max-penalty", type=float, default=0.05, help="reject samples above this composition penalty")
    f.add_argument("--target-seq", default="", help="target sequence inline")
    f.add_argument("--target-seq-file", default="", help="file containing the target sequence")
    f.add_argument("--out-fasta", required=True)
    f.add_argument("--out-manifest", required=True)
    f.add_argument("--minutes-per-fold", type=float, default=5.0)
    f.set_defaults(func=cmd_filter)

    s = sub.add_parser("select", help="rank fold results by shaped reward")
    s.add_argument("--fold-dir", required=True)
    s.add_argument("--manifest", required=True)
    s.add_argument("--out-csv", required=True)
    s.set_defaults(func=cmd_select)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

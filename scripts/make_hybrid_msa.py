#!/usr/bin/env python3
"""Build "hybrid" MSAs: real alignment for the target, single sequence for the binder.

WHY
---
The fast screen (--msa-mode single_sequence) was found to carry no ranking
signal. The mechanism is not subtle: with no MSA the TARGET does not fold.
Measured on 60 RBD complexes, same sequences both ways:

    fast screen   target pLDDT 27.2   binder pLDDT 70.1
    high accuracy target pLDDT 84.3   binder pLDDT 72.3

pLDDT 27 is random coil. The fast screen was scoring an interface against an
unfolded target. Meanwhile the binder folds fine either way (70 -> 72), because
a de novo sequence has no homologs for an MSA to find.

So the MSA matters for the target and is nearly worthless for the binder -- and
the target is IDENTICAL across every design in a campaign. Compute it once,
reuse it for all N designs, and keep single-sequence for the binder.

FORMAT
------
A ColabFold multimer a3m is:

    #<len1>,<len2>\t1,1
    >101\t102
    <chain1 query><chain2 query>
    >101
    <chain1 query><gaps>
    >homolog                       (chain-1 block, chain-2 region gapped)
    ...
    >102
    <gaps><chain2 query>           (chain-2 block)

a3m uses lowercase for insertions, which do NOT consume alignment columns, so
splitting a row by chain has to count only uppercase and gap characters.
"""
import argparse
import os
import sys


def split_a3m_row(row, n_cols_first):
    """Split an a3m row after n_cols_first alignment columns.

    Lowercase characters are insertions and do not count as columns. An
    insertion run belongs to the chain whose column it follows.
    """
    cols = 0
    for i, ch in enumerate(row):
        if ch.islower():
            continue            # insertion: no column consumed
        if cols == n_cols_first:
            return row[:i], row[i:]
        cols += 1
    return row, ""


def read_a3m(path):
    header = None
    entries = []
    name = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                header = line
            elif line.startswith(">"):
                name = line
            elif name is not None:
                entries.append((name, line))
                name = None
    return header, entries


def extract_target_msa(a3m_path):
    """Pull the target-only alignment out of an existing complex a3m."""
    header, entries = read_a3m(a3m_path)
    if not header or not header.startswith("#"):
        sys.exit(f"[FAIL] {a3m_path} has no '#len1,len2' header")
    lens = header[1:].split()[0].split(",")
    if len(lens) != 2:
        sys.exit(f"[FAIL] expected a 2-chain a3m, got header {header!r}")
    t_len, b_len = int(lens[0]), int(lens[1])

    target_rows = []
    seen = set()
    for name, seq in entries:
        first, _ = split_a3m_row(seq, t_len)
        # keep rows that actually align to the target; drop the binder-only block
        if first.replace("-", "").strip() == "":
            continue
        key = first
        if key in seen:
            continue
        seen.add(key)
        target_rows.append((name.split("\t")[0], first))
    return t_len, b_len, target_rows


def build_hybrid(target_rows, t_len, target_seq, binder_seq, name):
    """One multimer a3m: full target MSA, binder as a lone sequence."""
    b_len = len(binder_seq)
    out = [f"#{t_len},{b_len}\t1,1"]
    out.append(">101\t102")
    out.append(target_seq + binder_seq)
    # chain-1 block: target query + all its homologs, binder region gapped
    for rname, rseq in target_rows:
        out.append(rname)
        out.append(rseq + "-" * b_len)
    # chain-2 block: binder alone, no homologs (there are none)
    out.append(">102")
    out.append("-" * t_len + binder_seq)
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-a3m", required=True,
                    help="any existing complex a3m for this target (donates the target MSA)")
    ap.add_argument("--binders", required=True,
                    help="FASTA of binder sequences, or a CSV with a 'sequence' column")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args()

    t_len, b_len_src, target_rows = extract_target_msa(args.source_a3m)
    target_seq = target_rows[0][1].replace("-", "").upper() if target_rows else ""
    # the first row is the query itself; recover it cleanly
    _, entries = read_a3m(args.source_a3m)
    q_first, _ = split_a3m_row(entries[0][1], t_len)
    target_seq = "".join(c for c in q_first if not c.islower()).replace("-", "")

    print(f"Target MSA from {os.path.basename(args.source_a3m)}")
    print(f"  target length : {t_len}  ({len(target_seq)} recovered)")
    print(f"  MSA depth     : {len(target_rows)} unique rows (reused for every design)")

    # load binders
    binders = []
    if args.binders.endswith(".csv"):
        import pandas as pd
        df = pd.read_csv(args.binders)
        col = "sequence"
        for i, s in enumerate(df[col].tolist()):
            binders.append((f"design_{i:03d}", s))
    else:
        name = None
        for line in open(args.binders):
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
            elif line and name:
                binders.append((name, line.split(":")[-1]))
                name = None
    if args.limit:
        binders = binders[:args.limit]

    os.makedirs(args.out_dir, exist_ok=True)
    for name, seq in binders:
        a3m = build_hybrid(target_rows, t_len, target_seq, seq, name)
        with open(os.path.join(args.out_dir, f"{name}.a3m"), "w") as f:
            f.write(a3m)

    print(f"[OK] wrote {len(binders)} hybrid a3m files -> {args.out_dir}")
    print("     fold with: colabfold_batch <out-dir> <results>  (no --msa-mode; a3m inputs)")


if __name__ == "__main__":
    main()

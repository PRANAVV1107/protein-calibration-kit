#!/usr/bin/env python3
"""Fetch and VERIFY candidate decoy structures for specificity screening.

Every structure is checked by reading its own TITLE/COMPND records before use.
This project previously designed a full campaign against PDB 1ILR believing it
was IL-7Ra; it is interleukin-1 receptor antagonist. Nothing downstream catches
that, so verification happens here or not at all.

Decoy selection rationale
-------------------------
The informative decoys are structural HOMOLOGS, not random proteins. IL-7Ra is a
type-I cytokine receptor with fibronectin-III ectodomains; a binder that also
hits another family member is a real specificity failure, whereas one that hits a
viral protein with an unrelated fold is a curiosity. So the panel is tiered:

  tier 1  cytokine-receptor family -- the hard, meaningful cases
  tier 2  abundant serum proteins  -- what a therapeutic actually meets in vivo
  tier 3  unrelated folds          -- easy controls (PD-L1, RBD; already held)
"""
import os
import re
import sys
import urllib.request

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
OUT = os.path.join(KIT, "benchmark", "decoys")

# pdb_id -> (what we expect to find, tier, note)
CANDIDATES = {
    "2B5I": ("interleukin-2", "tier1", "IL-2 quaternary complex: IL-2Ra/IL-2Rb/gamma-c"),
    "1IAR": ("interleukin-4", "tier1", "IL-4 with IL-4R alpha ectodomain"),
    "3TGX": ("interleukin-21", "tier1", "IL-21 receptor"),
    "1AO6": ("albumin", "tier2", "human serum albumin, most abundant serum protein"),
    "1FC1": ("immunoglobulin", "tier2", "IgG Fc fragment"),
}

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
       "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
       "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def fetch(pdb_id, path):
    url = "https://files.rcsb.org/download/%s.pdb" % pdb_id
    try:
        with urllib.request.urlopen(url, timeout=60) as r, open(path, "wb") as f:
            f.write(r.read())
        return True
    except Exception as e:
        print("    [FAIL] download: %s" % e)
        return False


def header(path):
    title, compnds = [], []
    for line in open(path, errors="ignore"):
        if line.startswith("TITLE"):
            title.append(line[10:].strip())
        elif line.startswith("COMPND") and "MOLECULE:" in line:
            compnds.append(line.split("MOLECULE:")[1].strip().rstrip(";"))
        elif line.startswith(("ATOM", "HETATM")):
            break
    return " ".join(title), compnds


def chains(path):
    """chain id -> (sequence, n_residues), ATOM records only."""
    out = {}
    seen = {}
    for line in open(path, errors="ignore"):
        if not line.startswith("ATOM"):
            continue
        ch = line[21]
        key = line[22:27]
        rn = line[17:20].strip()
        if (ch, key) in seen or rn not in AA3:
            continue
        seen[(ch, key)] = 1
        out.setdefault(ch, []).append(AA3[rn])
    return {c: ("".join(s), len(s)) for c, s in out.items()}


def main():
    os.makedirs(OUT, exist_ok=True)
    ok, bad = [], []
    for pdb_id, (expect, tier, note) in CANDIDATES.items():
        path = os.path.join(OUT, "%s.pdb" % pdb_id)
        print("\n=== %s  [%s] ===" % (pdb_id, tier))
        print("    expected: %s   (%s)" % (expect, note))
        if not os.path.exists(path) and not fetch(pdb_id, path):
            bad.append((pdb_id, "download failed")); continue
        if os.path.getsize(path) < 2000:
            bad.append((pdb_id, "file too small")); continue

        title, mols = header(path)
        print("    TITLE   : %s" % title[:110])
        for m in mols[:5]:
            print("      MOL   : %s" % m[:95])

        blob = (title + " " + " ".join(mols)).lower()
        verified = expect.lower() in blob
        print("    VERIFIED: %s" % ("yes" if verified else "NO -- expected text not found"))
        if not verified:
            bad.append((pdb_id, "header does not mention %r" % expect)); continue

        ch = chains(path)
        for c, (seq, n) in sorted(ch.items(), key=lambda kv: -kv[1][1])[:4]:
            print("      chain %s: %3d res  %s..." % (c, n, seq[:42]))
        ok.append(pdb_id)

    print("\n" + "=" * 66)
    print("verified: %s" % (", ".join(ok) if ok else "none"))
    for pid, why in bad:
        print("rejected: %s -- %s" % (pid, why))
    print("\nNext: choose the specific chain per structure, then extract it as a")
    print("decoy target. Nothing is used until its chain is picked deliberately.")


if __name__ == "__main__":
    main()

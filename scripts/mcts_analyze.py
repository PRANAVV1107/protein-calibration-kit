#!/usr/bin/env python3
"""Pre-registered outcomes for the search run.

A subtlety that produced a bogus result on the first attempt: MCTS backprop
propagates `best` UP the tree, so a node's `best` field in state.json is the
maximum over its entire subtree, not its own score. The root reads 0.810
while its own evaluation was 0.150. Comparing a child's `best` against its
parent's is therefore structurally guaranteed to fail -- the parent's value
already contains the child's -- and returns 0% every time.

Own scores exist only in the run log ("child N -> 0.650"), so they are parsed
from there and joined to the tree structure in state.json.
"""
import json, re, sys
import numpy as np

D = sys.argv[1] if len(sys.argv) > 1 else "results/rl/mcts_run1"
LOG = sys.argv[2] if len(sys.argv) > 2 else "results/mcts_run1.log"
BAR = 0.87          # pre-registered: blind sampling, 60 folds
ROOT_BASE = 0.52    # pre-registered: binder_47 best-of-6

own = {}
for line in open(LOG, errors="ignore"):
    m = re.search(r"child (\d+) -> (?:value )?([0-9.]+)", line)
    if m:
        own[int(m.group(1))] = float(m.group(2))
    m = re.search(r"root score ([0-9.]+)", line)
    if m:
        own[0] = float(m.group(1))

st = json.load(open(D + "/state.json"))
N = {n["id"]: n for n in st["nodes"]}

print("=" * 70)
print("SEARCH -- pre-registered outcomes")
print("=" * 70)
print("  folds used %d   nodes %d   with own score %d"
      % (st["folds_used"], len(N), len(own)))

wins = tot = 0
for nid, sc in own.items():
    p = N.get(nid, {}).get("parent")
    if p is None or p not in own:
        continue
    tot += 1
    wins += int(sc > own[p])
print("\nSECONDARY -- children beating their own parent")
if tot:
    f = wins / tot
    se = (f * (1 - f) / tot) ** 0.5
    lo, hi = 100 * (f - 1.96 * se), 100 * (f + 1.96 * se)
    verdict = ("ABOVE chance" if lo > 50 else
               "BELOW chance" if hi < 50 else "not distinguishable from chance")
    print("  %d of %d = %.1f%%   95%% CI [%.1f%%, %.1f%%]   %s"
          % (wins, tot, 100 * f, lo, hi, verdict))

print("\nOWN score by depth")
for d in sorted({N[i]["depth"] for i in own if i in N}):
    s = [own[i] for i in own if i in N and N[i]["depth"] == d]
    print("  depth %d: n=%2d  mean %.3f  max %.3f" % (d, len(s), np.mean(s), max(s)))

best_id = max(own, key=own.get)
best = own[best_id]
print("\nPRIMARY")
print("  root own score        %.3f" % own.get(0, float("nan")))
print("  best found            %.3f  (node %d, depth %d)"
      % (best, best_id, N[best_id]["depth"]))
print("  bar (blind, 60 folds) %.3f  -> %s"
      % (BAR, "CLEARED" if best > BAR else "NOT cleared"))
print("  root best-of-6        %.3f  -> %s"
      % (ROOT_BASE, "improved on root" if best > ROOT_BASE else "did not improve on root"))

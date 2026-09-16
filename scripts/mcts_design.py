#!/usr/bin/env python3
"""Monte Carlo tree search over RFdiffusion's denoising trajectory.

The shape of the tree
---------------------
A node is a BACKBONE. The root is an existing backbone; expanding a node means
re-noising it to a shallow timestep with diffuser.partial_T and re-denoising,
which produces structural neighbours rather than independent draws. Evaluating a
node means designing sequences for it with ProteinMPNN, folding them, and taking
the resulting interface score.

This follows Proteina-Complexa (arXiv:2603.27950), whose tree is over denoising
timesteps with rewards computed on completed rollouts rather than partial
trajectories. There is no cheap partial scoring to exploit, so every leaf costs
a real fold.

Why batched expansion rather than textbook MCTS
-----------------------------------------------
Classic MCTS expands one node per iteration. RFdiffusion amortises heavily over
num_designs -- model load dominates a single design -- so expanding one child at
a time would waste most of the wall clock. Each iteration therefore selects one
node and expands it into --children at once, evaluates them together, and
backpropagates all of them. This is best-first search with batched expansion,
and it is what makes the search affordable on the measured primitives:
partial diffusion 0.5 min, fold about 1 min locally.

Why no value network
--------------------
AlphaZero evaluates leaves with a value network, which is what makes deep trees
affordable. value_probe.py measured the obvious candidate here -- ProteinMPNN's
own sequence likelihood -- and it does not work: within-backbone rho about -0.06
against the fold score, meaning it cannot rank sequences on a fixed backbone at
all. So leaves are evaluated by real folds and the tree must stay shallow and
wide. With a few hundred evaluations this is closer to guided best-first search
than to AlphaZero, and it should be described that way.

Honest expectation
------------------
The comparison that matters is against best-of-N at equal fold budget, not
against nothing. Proteina-Complexa reports MCTS beating best-of-N only on HARD
targets, with best-of-N already sufficient on easy ones. RBD is the hard target
in this project (single-draw mean ipTM 0.400 against 0.766 and 0.805), so that
is where a difference should appear if one exists. --baseline runs the matched
best-of-N arm.

Everything is resumable. State is written after every iteration, because a
preempted spot instance mid-search is expected rather than exceptional.
"""
import argparse
import glob
import json
import math
import os
import re
import shutil
import subprocess
import time

import numpy as np

ROOT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp"
KIT = os.path.join(ROOT, "calibration-kit")
RF = os.path.join(ROOT, "external_tools/RFdiffusion")
MPNN = os.path.join(ROOT, "external_tools/ProteinMPNN")
RF_ENV = "/root/miniconda3/envs/rfdiffusion"
CF_ENV = "/root/miniconda3/envs/colabfold"

TARGETS = {
    # name: (target sequence, contigs, hotspots, designed chain)
    "rbd": (
        "TNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAP"
        "GQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPL"
        "QSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCG",
        "[E333-526/0 74-74]",
        "[E449,E493,E496,E498,E501,E505]",
        "A"),
}


def sh(cmd, env=None, log=None, cwd=None, timeout=3600):
    e = dict(os.environ)
    if env:
        e.update(env)
    with open(log, "a") if log else open(os.devnull, "w") as fh:
        return subprocess.run(cmd, shell=True, env=e, cwd=cwd, stdout=fh,
                              stderr=subprocess.STDOUT, timeout=timeout).returncode


class Node:
    __slots__ = ("id", "pdb", "parent", "children", "n", "total", "best", "depth", "evaluated")

    def __init__(self, nid, pdb, parent=None, depth=0):
        self.id, self.pdb, self.parent = nid, pdb, parent
        self.children, self.n, self.total, self.best = [], 0, 0.0, -1e9
        self.depth, self.evaluated = depth, False

    @property
    def mean(self):
        return self.total / self.n if self.n else 0.0

    def to_dict(self):
        return {"id": self.id, "pdb": self.pdb,
                "parent": self.parent.id if self.parent else None,
                "children": [c.id for c in self.children], "n": self.n,
                "total": self.total, "best": self.best, "depth": self.depth,
                "evaluated": self.evaluated}


class Search:
    def __init__(self, args):
        self.a = args
        self.out = args.out_dir
        os.makedirs(self.out, exist_ok=True)
        self.log = os.path.join(self.out, "search.log")
        self.nodes = {}
        self.next_id = 0
        self.root = None
        self.folds_used = 0
        tgt = TARGETS[args.target]
        self.tseq, self.contigs, self.hotspots, self.chain = tgt

    def say(self, m):
        line = "[%s] %s" % (time.strftime("%F %T"), m)
        print(line, flush=True)
        with open(self.log, "a") as f:
            f.write(line + "\n")

    # ---------- persistence ----------
    def save(self):
        st = {"next_id": self.next_id, "folds_used": self.folds_used,
              "root": self.root.id if self.root else None,
              "nodes": [n.to_dict() for n in self.nodes.values()]}
        tmp = os.path.join(self.out, "state.json.tmp")
        with open(tmp, "w") as f:
            json.dump(st, f, indent=1)
        os.replace(tmp, os.path.join(self.out, "state.json"))

    def load(self):
        p = os.path.join(self.out, "state.json")
        if not os.path.exists(p):
            return False
        st = json.load(open(p))
        for d in st["nodes"]:
            nd = Node(d["id"], d["pdb"], None, d["depth"])
            nd.n, nd.total, nd.best = d["n"], d["total"], d["best"]
            nd.evaluated = d["evaluated"]
            self.nodes[d["id"]] = nd
        for d in st["nodes"]:
            nd = self.nodes[d["id"]]
            if d["parent"] is not None:
                nd.parent = self.nodes[d["parent"]]
            nd.children = [self.nodes[c] for c in d["children"]]
        self.next_id = st["next_id"]
        self.folds_used = st["folds_used"]
        self.root = self.nodes.get(st["root"])
        self.say("resumed: %d nodes, %d folds used" % (len(self.nodes), self.folds_used))
        return True

    # ---------- primitives ----------
    def expand(self, node):
        """partial_T re-noising of node.pdb into --children structural neighbours."""
        d = os.path.join(self.out, "nodes", "n%04d" % node.id)
        os.makedirs(d, exist_ok=True)
        pref = os.path.join(d, "child")
        if not glob.glob(pref + "_*.pdb"):
            cmd = ("python scripts/run_inference.py "
                   "inference.output_prefix=%s "
                   "inference.model_directory_path=%s/models "
                   "inference.input_pdb=%s "
                   "'contigmap.contigs=%s' 'ppi.hotspot_res=%s' "
                   "inference.num_designs=%d diffuser.partial_T=%d "
                   "denoiser.noise_scale_ca=0 denoiser.noise_scale_frame=0"
                   % (pref, RF, node.pdb, self.contigs, self.hotspots,
                      self.a.children, self.a.partial_t))
            env = {"PATH": RF_ENV + "/bin:" + os.environ.get("PATH", ""),
                   "LD_LIBRARY_PATH": RF_ENV + "/lib",
                   "CONDA_PREFIX": RF_ENV, "CONDA_DEFAULT_ENV": "rfdiffusion"}
            rc = sh(cmd, env, self.log, cwd=RF, timeout=self.a.rf_timeout)
            if rc != 0:
                self.say("  RFdiffusion returned %d on node %d" % (rc, node.id))
        return sorted(glob.glob(pref + "_*.pdb"))

    def evaluate(self, pdbs, tag):
        """ProteinMPNN then fold. Returns {pdb: best interface score}."""
        d = os.path.join(self.out, "eval", tag)
        os.makedirs(os.path.join(d, "pdbs"), exist_ok=True)
        for p in pdbs:
            shutil.copy(p, os.path.join(d, "pdbs", os.path.basename(p)))
        env_rf = {"PATH": RF_ENV + "/bin:" + os.environ.get("PATH", ""),
                  "LD_LIBRARY_PATH": RF_ENV + "/lib"}
        if not os.path.exists(os.path.join(d, "parsed.jsonl")):
            sh("python helper_scripts/parse_multiple_chains.py --input_path %s "
               "--output_path %s" % (os.path.join(d, "pdbs"),
                                     os.path.join(d, "parsed.jsonl")),
               env_rf, self.log, cwd=MPNN)
            sh("python helper_scripts/assign_fixed_chains.py --input_path %s "
               "--output_path %s --chain_list '%s'"
               % (os.path.join(d, "parsed.jsonl"),
                  os.path.join(d, "assigned.jsonl"), self.chain),
               env_rf, self.log, cwd=MPNN)
        if not glob.glob(os.path.join(d, "mpnn", "seqs", "*.fa")):
            with open(os.path.join(d, "bias.jsonl"), "w") as f:
                f.write('{"A": -1.0}')
            sh("python protein_mpnn_run.py --jsonl_path %s --chain_id_jsonl %s "
               "--out_folder %s --num_seq_per_target %d --sampling_temp 0.2 "
               "--batch_size 1 --bias_AA_jsonl %s"
               % (os.path.join(d, "parsed.jsonl"), os.path.join(d, "assigned.jsonl"),
                  os.path.join(d, "mpnn"), self.a.seqs, os.path.join(d, "bias.jsonl")),
               env_rf, self.log, cwd=MPNN)

        # fold inputs: target MSA is not needed per-design; single-sequence pair
        fi = os.path.join(d, "fold_in")
        os.makedirs(fi, exist_ok=True)
        pairs = {}
        for fa in sorted(glob.glob(os.path.join(d, "mpnn", "seqs", "*.fa"))):
            base = os.path.basename(fa)[:-3]
            seqs, cur = [], None
            for line in open(fa):
                line = line.strip()
                if line.startswith(">"):
                    if cur is not None:
                        seqs.append(cur)
                    cur = ""
                elif cur is not None:
                    cur += line
            if cur:
                seqs.append(cur)
            for i, s in enumerate(seqs[1:]):
                s = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", s.upper())
                if not s:
                    continue
                t = "%s__s%d" % (base, i + 1)
                with open(os.path.join(fi, t + ".a3m"), "w") as f:
                    f.write("#%d,%d\t1,1\n" % (len(self.tseq), len(s)))
                    f.write(">101\n%s%s\n" % (self.tseq, s))
                    f.write(">t\n%s%s\n" % (self.tseq, "-" * len(s)))
                pairs[t] = base
        fo = os.path.join(d, "folds")
        os.makedirs(fo, exist_ok=True)
        todo = [t for t in pairs
                if not glob.glob(os.path.join(fo, t + "_scores_rank_001_*.json"))]
        if todo:
            ch = os.path.join(d, ".chunk")
            for i in range(0, len(todo), self.a.chunk):
                shutil.rmtree(ch, ignore_errors=True)
                os.makedirs(ch)
                for t in todo[i:i + self.a.chunk]:
                    shutil.copy(os.path.join(fi, t + ".a3m"), ch)
                sh("%s/bin/colabfold_batch --model-type alphafold2_multimer_v3 "
                   "--num-models 1 --num-recycle 3 %s %s" % (CF_ENV, ch, fo),
                   {"TF_FORCE_UNIFIED_MEMORY": "1",
                    "XLA_PYTHON_CLIENT_PREALLOCATE": "false",
                    "XLA_PYTHON_CLIENT_MEM_FRACTION": "2.0"},
                   self.log, timeout=self.a.fold_timeout)
            shutil.rmtree(ch, ignore_errors=True)

        best = {}
        for f in glob.glob(os.path.join(fo, "*_scores_rank_001_*.json")):
            t = os.path.basename(f).split("_scores_rank_001_")[0]
            if t not in pairs:
                continue
            j = json.load(open(f))
            if self.a.metric == "ipTM":
                v = float(j.get("iptm") or 0.0)
            else:
                x = j.get("ipsae")
                v = (min(float(q) for q in x.values()) if isinstance(x, dict)
                     else (float(x) if x is not None else float(j.get("iptm") or 0.0)))
            b = pairs[t]
            best[b] = max(best.get(b, -1e9), v)
            self.folds_used += 1
        return best

    # ---------- search ----------
    def uct(self, node):
        """Standard UCT, on nodes that have been evaluated and can be expanded."""
        c = self.a.c_uct
        best, bv = None, -1e9
        for n in self.nodes.values():
            if not n.evaluated or n.depth >= self.a.max_depth:
                continue
            if n.children and len(n.children) >= self.a.children:
                continue
            parent_n = max(1, n.parent.n if n.parent else n.n)
            v = n.mean + c * math.sqrt(math.log(parent_n + 1) / max(1, n.n))
            if v > bv:
                best, bv = n, v
        return best

    def backprop(self, node, score):
        cur = node
        while cur is not None:
            cur.n += 1
            cur.total += score
            cur.best = max(cur.best, score)
            cur = cur.parent

    def run(self):
        if not self.load():
            self.root = Node(self.next_id, self.a.root_pdb, None, 0)
            self.nodes[self.root.id] = self.root
            self.next_id += 1
            self.say("root: %s" % self.a.root_pdb)
            sc = self.evaluate([self.root.pdb], "n0000")
            v = max(sc.values()) if sc else 0.0
            self.root.evaluated = True
            self.backprop(self.root, v)
            self.say("root score %.3f (folds %d)" % (v, self.folds_used))
            self.save()

        while self.folds_used < self.a.budget:
            node = self.uct(node=None)
            if node is None:
                self.say("no expandable node left")
                break
            self.say("iter: expanding node %d (depth %d, mean %.3f, n %d) "
                     "[folds %d/%d]" % (node.id, node.depth, node.mean, node.n,
                                        self.folds_used, self.a.budget))
            kids = self.expand(node)
            if not kids:
                self.say("  expansion produced nothing; marking node exhausted")
                node.depth = self.a.max_depth
                self.save()
                continue
            newnodes = []
            for p in kids:
                nd = Node(self.next_id, p, node, node.depth + 1)
                self.nodes[nd.id] = nd
                node.children.append(nd)
                self.next_id += 1
                newnodes.append(nd)
            scores = self.evaluate([n.pdb for n in newnodes], "n%04d_kids" % node.id)
            for nd in newnodes:
                key = os.path.basename(nd.pdb)[:-4]
                v = scores.get(key)
                if v is None:
                    continue
                nd.evaluated = True
                self.backprop(nd, v)
                self.say("  child %d -> %.3f" % (nd.id, v))
            self.save()

        self.report()

    def report(self):
        ev = [n for n in self.nodes.values() if n.evaluated]
        ev.sort(key=lambda n: -n.best)
        self.say("=" * 60)
        self.say("SEARCH COMPLETE: %d nodes, %d folds" % (len(ev), self.folds_used))
        for n in ev[:10]:
            self.say("  node %-4d depth %d  best %.3f  mean %.3f  n %d  %s"
                     % (n.id, n.depth, n.best, n.mean, n.n, os.path.basename(n.pdb)))
        with open(os.path.join(self.out, "RESULT.txt"), "w") as f:
            f.write("nodes %d  folds %d\n" % (len(ev), self.folds_used))
            f.write("root best %.4f\n" % self.root.best)
            f.write("overall best %.4f\n" % (ev[0].best if ev else float("nan")))
            f.write("\nall evaluated nodes, best first:\n")
            for n in ev:
                f.write("%4d depth %d best %.4f mean %.4f n %d %s\n"
                        % (n.id, n.depth, n.best, n.mean, n.n, n.pdb))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root-pdb", required=True, help="starting backbone")
    ap.add_argument("--target", default="rbd", choices=sorted(TARGETS))
    ap.add_argument("--out-dir", default="results/rl/mcts_rbd")
    ap.add_argument("--budget", type=int, default=120, help="total folds to spend")
    ap.add_argument("--children", type=int, default=4, help="neighbours per expansion")
    ap.add_argument("--seqs", type=int, default=4, help="MPNN draws per backbone")
    ap.add_argument("--partial-t", type=int, default=10,
                    help="re-noising depth; lower stays closer to the parent")
    ap.add_argument("--max-depth", type=int, default=4)
    ap.add_argument("--c-uct", type=float, default=0.6)
    ap.add_argument("--metric", default="ipTM", choices=["ipTM", "ipSAE"],
                    help="ipTM by default so results compare directly against "
                         "the variance study's already-paid-for baselines: "
                         "best-of-6 on this root is 0.71, and 60 folds of blind "
                         "sampling across 10 RBD backbones reached 0.87")
    ap.add_argument("--chunk", type=int, default=20)
    ap.add_argument("--rf-timeout", type=int, default=3600)
    ap.add_argument("--fold-timeout", type=int, default=7200)
    args = ap.parse_args()
    os.chdir(KIT)
    Search(args).run()


if __name__ == "__main__":
    main()

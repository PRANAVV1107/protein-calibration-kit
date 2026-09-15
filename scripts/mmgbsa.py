#!/usr/bin/env python3
"""MM-GBSA binding free energy for predicted binder complexes.

STATUS: PARKED, NOT VALIDATED -- do not use these numbers
---------------------------------------------------------
The code runs end to end and splits chains correctly, but the one test done so
far returned dG = +264,742 kJ/mol on adaptyv_145 (target 201 res, binder 27).
A binding free energy must be negative and of order -50 to -500 kJ/mol, so that
figure is meaningless. The cause is understood:

We feed this the *unrelaxed* ColabFold output. ColabFold's Amber relaxation
step is what removes AF2's small clashes, and it never ran here --
--use-gpu-relax fails because the colabfold env's OpenMM has no CUDA platform,
and no CPU-relax fallback was run. So every PDB carries clash energy that
swamps the interface term: E_before_min was 397,873 kJ/mol and 20 minimisation
steps only reached 359,953. Even the isolated target came out at +85,630.

Before trusting anything from this script:
  1. Relax the structures first (ColabFold Amber relax in the openmm-gpu env),
     or minimise to convergence here rather than for a fixed step count. The
     relaxed PDBs are reusable by md_stability.py, which argues for doing it
     upstream once rather than inside this script 165 times.
  2. Re-run and confirm dG lands in a physical range on a known complex before
     computing it across the set.
  3. Time it per design on CUDA. Cost is half the decision: a reward term that
     needs an energy minimisation per candidate is expensive inside an RL loop
     specifically, where it is paid on every candidate at every step. "Predictive
     but too costly to use as a reward" is a legitimate outcome.

CPU is not viable regardless: 44 s for 20 steps on the smallest complex in the
set (1724 atoms, 4 threads), and it competes with ColabFold's CPU-bound MSA
stage.

Why this is worth running
-------------------------
ipTM, ipSAE, pDockQ2 and pLDDT are all the same network's opinion re-expressed,
so they cannot disagree with one another in an informative way. A force-field
energy is a physically different question, which is the only reason it could
add predictive power over metrics we already have.

That "could" is the point. This project is about to put dG into a reinforcement
learning reward alongside ipTM, ipSAE and a specificity margin. Every one of
those terms except dG has now been checked against the 165 Adaptyv designs with
MEASURED binding outcomes. Weighting an unvalidated term is how a reward
function learns to optimise something that does not predict binding, so this
script exists to give dG the same test the others got.

Method -- single-trajectory MM-GBSA, no entropy term
----------------------------------------------------
    dG_bind = E_complex - E_target - E_binder

All three energies come from the SAME minimised complex coordinates, the
standard single-trajectory approximation. It cancels the intramolecular terms
that dominate and would otherwise swamp the interface signal, at the cost of
ignoring conformational relaxation on unbinding -- acceptable here because the
comparison is between designs against one shared target, not between systems.

No configurational entropy (-T*dS). Normal-mode or quasi-harmonic entropy
roughly doubles the cost and typically adds more noise than signal for ranking
within a congeneric set, which is exactly this case. So the number below is an
enthalpy-like score for RANKING, not an absolute affinity. It must never be
reported in kcal/mol as though comparable to an experimental dG.

Minimisation first, always. AF2 structures carry small clashes that produce
enormous positive Lennard-Jones energies; a single-point evaluation on the raw
prediction is dominated by them and is meaningless. A few hundred steps of
L-BFGS removes them without materially moving the interface.

Implicit solvent (GBn2) supplies the desolvation penalty that makes this a
binding free energy rather than a vacuum interaction energy.
"""
import argparse
import glob
import json
import os
import sys
import time

import numpy as np

try:
    import openmm
    from openmm import app, unit
    from pdbfixer import PDBFixer
except Exception as e:  # pragma: no cover
    sys.exit("[FAIL] OpenMM/pdbfixer unavailable: %s" % e)

FF = ("amber14-all.xml", "implicit/gbn2.xml")
KJ = unit.kilojoule_per_mole


def prepare(pdb_path):
    """Add hydrogens and any missing heavy atoms. Never build in residues."""
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.missingResidues = {}          # modelling unresolved residues is invention
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    return fixer.topology, fixer.positions


def platform_for(name, threads):
    if name == "CUDA":
        try:
            p = openmm.Platform.getPlatformByName("CUDA")
            return p, {"Precision": "mixed"}
        except Exception:
            pass
    p = openmm.Platform.getPlatformByName("CPU")
    return p, {"Threads": str(threads)}


def energy_of(topology, positions, ff, platform, props):
    """Potential energy of exactly this topology at these coordinates."""
    system = ff.createSystem(topology,
                             nonbondedMethod=app.NoCutoff,
                             constraints=None,
                             rigidWater=False)
    integ = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    ctx = openmm.Context(system, integ, platform, props)
    ctx.setPositions(positions)
    e = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KJ)
    del ctx, integ, system
    return float(e)


def subset(topology, positions, keep_chain_ids):
    """A new topology/positions holding only the named chains."""
    modeller = app.Modeller(topology, positions)
    drop = [c for c in modeller.topology.chains() if c.id not in keep_chain_ids]
    modeller.delete(drop)
    return modeller.topology, modeller.positions


def run_one(pdb_path, platform_name="CPU", threads=2, min_steps=500):
    name = os.path.basename(pdb_path).split("_unrelaxed_")[0]
    t0 = time.time()
    topology, positions = prepare(pdb_path)

    chains = list(topology.chains())
    if len(chains) < 2:
        return {"name": name, "error": "single chain"}
    # target is the larger chain; ids are what Modeller.delete matches on
    sizes = [(len(list(c.atoms())), c.id) for c in chains]
    sizes.sort(reverse=True)
    target_id, binder_id = sizes[0][1], sizes[1][1]

    ff = app.ForceField(*FF)
    platform, props = platform_for(platform_name, threads)

    # ---- minimise the complex; everything downstream uses these coordinates ----
    system = ff.createSystem(topology, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds, rigidWater=False)
    integ = openmm.LangevinMiddleIntegrator(300 * unit.kelvin,
                                            1.0 / unit.picosecond,
                                            0.002 * unit.picoseconds)
    sim = app.Simulation(topology, system, integ, platform, props)
    sim.context.setPositions(positions)
    e_before = sim.context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KJ)
    sim.minimizeEnergy(maxIterations=min_steps)
    st = sim.context.getState(getPositions=True, getEnergy=True)
    minpos = st.getPositions()
    e_complex = st.getPotentialEnergy().value_in_unit(KJ)
    del sim, integ, system

    # ---- the two partners, at those same minimised coordinates ----
    t_top, t_pos = subset(topology, minpos, {target_id})
    b_top, b_pos = subset(topology, minpos, {binder_id})
    e_target = energy_of(t_top, t_pos, ff, platform, props)
    e_binder = energy_of(b_top, b_pos, ff, platform, props)

    dg = e_complex - e_target - e_binder
    return {
        "name": name,
        "dG_kj": dg,
        "dG_kcal": dg / 4.184,
        "E_complex": e_complex,
        "E_target": e_target,
        "E_binder": e_binder,
        "E_before_min": e_before,
        "min_steps": min_steps,
        "n_res_target": t_top.getNumResidues(),
        "n_res_binder": b_top.getNumResidues(),
        "platform": platform.getName(),
        "seconds": round(time.time() - t0, 1),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+", help="complex PDB files or a glob")
    ap.add_argument("--out-dir", default="results/mmgbsa")
    ap.add_argument("--platform", default="CPU", choices=["CPU", "CUDA"])
    ap.add_argument("--threads", type=int, default=2,
                    help="CPU threads; keep low while ColabFold is running")
    ap.add_argument("--min-steps", type=int, default=500)
    args = ap.parse_args()

    files = []
    for pat in args.inputs:
        files.extend(sorted(glob.glob(pat)) if any(c in pat for c in "*?[") else [pat])
    os.makedirs(args.out_dir, exist_ok=True)

    done = fail = 0
    for f in files:
        name = os.path.basename(f).split("_unrelaxed_")[0]
        out = os.path.join(args.out_dir, name + "_mmgbsa.json")
        if os.path.exists(out):
            continue
        try:
            r = run_one(f, args.platform, args.threads, args.min_steps)
        except Exception as e:
            r = {"name": name, "error": "%s: %s" % (type(e).__name__, e)}
        with open(out, "w") as fh:
            json.dump(r, fh, indent=2)
        if "error" in r:
            fail += 1
            print("[FAIL] %-16s %s" % (name, r["error"]), flush=True)
        else:
            done += 1
            print("[OK]   %-16s dG %9.1f kJ/mol (%7.1f kcal)  %5.1fs  %s"
                  % (name, r["dG_kj"], r["dG_kcal"], r["seconds"], r["platform"]),
                  flush=True)
    print("\n%d computed, %d failed -> %s" % (done, fail, args.out_dir))


if __name__ == "__main__":
    main()

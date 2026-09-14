#!/usr/bin/env python3
"""Molecular-dynamics stability screen for predicted binder complexes.

Why this is worth running
-------------------------
ipTM, ipSAE and pDockQ2 are all the same network's opinion re-expressed. They
cannot disagree with each other in an informative way. Molecular dynamics asks a
physically different question -- does the predicted interface survive when the
system is propagated under a force field? -- so it can reject a pose that every
confidence metric liked.

This is the only orthogonal evidence available without a wet lab.

Method
------
Implicit solvent (GBn2) rather than explicit water: roughly an order of
magnitude cheaper, and sufficient for the question being asked, which is gross
interface stability rather than precise energetics. Amber14 protein force field,
Langevin thermostat at 300 K, 2 fs timestep with hydrogen bonds constrained.

Platforms. The colabfold env's OpenMM exposes only Reference and CPU -- which is
also the real reason --use-gpu-relax fails there, rather than a CUDA version
mismatch. A separate openmm-gpu env carries the CUDA platform, which OpenMM
rates at speed 100 against the CPU platform's 1.

That difference decides how this is used. On CPU, 1 ns on one complex took over
90 minutes and was still unfinished, which buys a screen that can only reject
the grossest failures. On CUDA the same work is minutes, so 50-100 ns runs
become affordable and the result can actually support a stability claim rather
than merely a rejection. Use --platform CUDA once the GPU is free of folding.

Measured per frame
------------------
  interface_contacts  heavy-atom residue pairs across the interface within 5 A
  binder_rmsd         binder RMSD after superposing on the TARGET, so it
                      reports movement of the binder RELATIVE to its target
                      rather than overall tumbling
  com_distance        centre-of-mass separation of the two chains

A complex that holds has roughly flat contacts and a bounded binder RMSD. One
that dissociates shows contacts decaying toward zero and com_distance climbing.
"""
import argparse
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


def prepare(pdb_path, out_prefix):
    """Add missing atoms and hydrogens; return a fixed topology/positions."""
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    # do not try to build in unresolved terminal residues -- for predicted
    # structures there are none, and modelling them would be invention
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    with open(out_prefix + "_fixed.pdb", "w") as f:
        app.PDBFile.writeFile(fixer.topology, fixer.positions, f)
    return fixer.topology, fixer.positions


def chain_atom_indices(topology):
    chains = list(topology.chains())
    if len(chains) < 2:
        return None, None
    idx = [[a.index for a in c.atoms()] for c in chains]
    # target is the larger chain
    order = sorted(range(len(idx)), key=lambda i: -len(idx[i]))
    return idx[order[0]], idx[order[1]]


def heavy(topology, indices):
    s = set(indices)
    return [a.index for a in topology.atoms()
            if a.index in s and a.element is not None and a.element.symbol != "H"]


def metrics(pos, t_heavy, b_heavy, t_ref, b_ref, cutoff=0.5):
    """cutoff in nm (0.5 nm = 5 A)."""
    T = pos[t_heavy]
    B = pos[b_heavy]
    d = np.linalg.norm(T[:, None, :] - B[None, :, :], axis=-1)
    contacts = int((d < cutoff).sum())
    com = float(np.linalg.norm(T.mean(axis=0) - B.mean(axis=0)))

    # superpose on target, then measure binder displacement (Kabsch)
    P, Q = T - T.mean(0), t_ref - t_ref.mean(0)
    V, S, W = np.linalg.svd(P.T @ Q)
    dsign = np.sign(np.linalg.det(V @ W))
    D = np.diag([1.0, 1.0, dsign])
    R = V @ D @ W
    b_moved = (B - T.mean(0)) @ R
    b_start = b_ref - t_ref.mean(0)
    rmsd = float(np.sqrt(((b_moved - b_start) ** 2).sum(axis=1).mean()))
    return contacts, com, rmsd


def make_simulation(topology, system, integrator, platform_name, threads):
    """CUDA where available (about 100x the CPU platform), else CPU.

    The colabfold env's OpenMM build exposes only Reference and CPU -- which is
    the real reason --use-gpu-relax fails there, rather than a CUDA version
    mismatch. A separate openmm-gpu env carries the CUDA platform.
    """
    if platform_name == "CUDA":
        plat = openmm.Platform.getPlatformByName("CUDA")
        return app.Simulation(topology, system, integrator, plat,
                              {"Precision": "mixed"})
    plat = openmm.Platform.getPlatformByName("CPU")
    return app.Simulation(topology, system, integrator, plat, {"Threads": str(threads)})


def run(pdb_path, out_dir, ns=2.0, threads=4, equil_ps=100.0, report_ps=50.0,
        platform_name="CPU"):
    name = os.path.splitext(os.path.basename(pdb_path))[0]
    os.makedirs(out_dir, exist_ok=True)
    prefix = os.path.join(out_dir, name)

    topology, positions = prepare(pdb_path, prefix)
    t_idx, b_idx = chain_atom_indices(topology)
    if t_idx is None:
        return {"name": name, "error": "single chain, nothing to measure"}
    t_heavy, b_heavy = heavy(topology, t_idx), heavy(topology, b_idx)

    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = ff.createSystem(topology, nonbondedMethod=app.NoCutoff,
                             constraints=app.HBonds, hydrogenMass=1.5 * unit.amu)
    integrator = openmm.LangevinMiddleIntegrator(300 * unit.kelvin,
                                                 1.0 / unit.picosecond,
                                                 0.002 * unit.picoseconds)
    sim = make_simulation(topology, system, integrator, platform_name, threads)
    sim.context.setPositions(positions)

    sim.minimizeEnergy(maxIterations=2000)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin)
    sim.step(int(equil_ps / 0.002))                      # equilibration

    start = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    t_ref, b_ref = start[t_heavy], start[b_heavy]
    c0, com0, _ = metrics(start, t_heavy, b_heavy, t_ref, b_ref)

    n_steps = int(ns * 1000 / 0.002)
    stride = int(report_ps / 0.002)
    frames = []
    t0 = time.time()
    for _ in range(0, n_steps, stride):
        sim.step(stride)
        p = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        c, com, rmsd = metrics(p, t_heavy, b_heavy, t_ref, b_ref)
        frames.append({"contacts": c, "com_nm": round(com, 3), "binder_rmsd_nm": round(rmsd, 3)})

    cs = [f["contacts"] for f in frames]
    res = {
        "name": name,
        "ns": ns,
        "platform": platform_name,
        "wall_min": round((time.time() - t0) / 60, 1),
        "contacts_start": c0,
        "contacts_mean": float(np.mean(cs)),
        "contacts_final": cs[-1],
        "contacts_retained": round(cs[-1] / c0, 3) if c0 else None,
        "com_start_nm": round(com0, 3),
        "com_final_nm": frames[-1]["com_nm"],
        "binder_rmsd_final_nm": frames[-1]["binder_rmsd_nm"],
        "binder_rmsd_max_nm": round(max(f["binder_rmsd_nm"] for f in frames), 3),
        "dissociated": bool(cs[-1] < 0.25 * c0) if c0 else None,
        "frames": frames,
    }
    with open(prefix + "_md.json", "w") as f:
        json.dump(res, f, indent=2)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdbs", nargs="+")
    ap.add_argument("--out-dir", default="results/md")
    ap.add_argument("--ns", type=float, default=2.0)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--platform", default="CPU", choices=["CPU", "CUDA"],
                    help="CUDA needs the openmm-gpu env; it is ~100x the CPU platform")
    ap.add_argument("--equil-ps", type=float, default=100.0,
                    help="equilibration before production; 20 is enough for a "
                         "gross-stability screen on the CPU platform")
    args = ap.parse_args()

    for p in args.pdbs:
        print("[%s] starting %s" % (time.strftime("%H:%M:%S"), os.path.basename(p)), flush=True)
        try:
            r = run(p, args.out_dir, ns=args.ns, threads=args.threads,
                    equil_ps=args.equil_ps, platform_name=args.platform)
        except Exception as e:
            print("   FAILED: %s" % e, flush=True)
            continue
        if "error" in r:
            print("   %s" % r["error"], flush=True)
            continue
        print("   contacts %d -> %d (%.0f%% retained)   binder RMSD %.2f nm   "
              "COM %.2f -> %.2f nm   %s   [%.1f min]"
              % (r["contacts_start"], r["contacts_final"], 100 * r["contacts_retained"],
                 r["binder_rmsd_final_nm"], r["com_start_nm"], r["com_final_nm"],
                 "DISSOCIATED" if r["dissociated"] else "stable", r["wall_min"]), flush=True)


if __name__ == "__main__":
    main()

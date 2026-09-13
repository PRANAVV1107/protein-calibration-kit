#!/usr/bin/env python3
"""Generate the IL-7Ra campaign report as a self-contained HTML page.

Embeds backbone-only PDB coordinates for the three control folds so the
structures are viewable without any external fetch (the artifact CSP blocks
those anyway).
"""
import json
import os

KIT = "/mnt/c/Users/prana/Downloads/proteinfoldingexp/calibration-kit"
OUT = os.path.join(KIT, "il7ra_report.html")

STRUCTURES = {
    "ontarget": "results/win_final/win_5model_relaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb",
    "scrambled": "results/win_validate/win_scrambled_unrelaxed_rank_001_alphafold2_multimer_v3_model_2_seed_000.pdb",
    "offtarget": "results/win_validate/win_offtarget_PDL1_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb",
}

KEEP_ATOMS = {"N", "CA", "C", "O", "CB"}


def slim_pdb(path):
    """Keep backbone + CB only: enough for cartoon rendering, ~60% smaller."""
    out = []
    for line in open(path):
        if line.startswith(("ATOM", "TER", "END")):
            if line.startswith("ATOM"):
                if line[12:16].strip() not in KEEP_ATOMS:
                    continue
            out.append(line.rstrip())
    return "\n".join(out)


def main():
    pdbs = {k: slim_pdb(os.path.join(KIT, v)) for k, v in STRUCTURES.items()}
    for k, v in pdbs.items():
        print(f"  {k}: {len(v)/1024:.0f} KB")

    html = TEMPLATE.replace("__PDB_DATA__", json.dumps(pdbs))
    with open(OUT, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"[OK] {OUT}  ({os.path.getsize(OUT)/1024:.0f} KB)")


TEMPLATE = r"""<title>IL-7Rα Binder Campaign</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Serif:wght@500;600&display=swap">
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
<style>
:root{
  --bg:#F6F7F9; --surface:#FFFFFF; --surface-2:#EFF2F6;
  --ink:#10141C; --ink-2:#4C5769; --ink-3:#7A8598;
  --line:#DCE1E9; --line-2:#C3CBD8;
  --accent:#0B5FD0; --accent-soft:#E6EFFC;
  --good:#157F3D; --good-soft:#E4F2E9;
  --warn:#A8690B; --warn-soft:#FBF0DC;
  --crit:#B3261E; --crit-soft:#FBE7E5;
  --plddt-vhigh:#0053D6; --plddt-high:#65CBF3; --plddt-low:#FFDB13; --plddt-vlow:#FF7D45;
  --measure:68ch;
}
@media (prefers-color-scheme:dark){
  :root:not([data-theme="light"]){
    --bg:#0C1016; --surface:#141A23; --surface-2:#1B222D;
    --ink:#E9EDF4; --ink-2:#99A5B8; --ink-3:#6E7B8E;
    --line:#242C39; --line-2:#323C4C;
    --accent:#5AA0FF; --accent-soft:#152538;
    --good:#4FBE7B; --good-soft:#13251A;
    --warn:#E0A33C; --warn-soft:#2A2013;
    --crit:#F0736A; --crit-soft:#2C1715;
  }
}
:root[data-theme="dark"]{
  --bg:#0C1016; --surface:#141A23; --surface-2:#1B222D;
  --ink:#E9EDF4; --ink-2:#99A5B8; --ink-3:#6E7B8E;
  --line:#242C39; --line-2:#323C4C;
  --accent:#5AA0FF; --accent-soft:#152538;
  --good:#4FBE7B; --good-soft:#13251A;
  --warn:#E0A33C; --warn-soft:#2A2013;
  --crit:#F0736A; --crit-soft:#2C1715;
}
*{box-sizing:border-box}
body{
  margin:0; background:var(--bg); color:var(--ink);
  font-family:"IBM Plex Sans",system-ui,-apple-system,sans-serif;
  font-size:16px; line-height:1.6; -webkit-font-smoothing:antialiased;
}
.wrap{max-width:1180px; margin:0 auto; padding:0 28px 96px}
.col{max-width:var(--measure)}
h1,h2,h3{font-family:"IBM Plex Serif",Georgia,serif; text-wrap:balance; margin:0}
h1{font-size:clamp(30px,4.4vw,46px); font-weight:600; line-height:1.12; letter-spacing:-.015em}
h2{font-size:25px; font-weight:600; line-height:1.25; margin-bottom:6px}
h3{font-size:18px; font-weight:600; margin-bottom:4px}
p{margin:0}
a{color:var(--accent)}
.mono{font-family:"IBM Plex Mono",ui-monospace,monospace; font-variant-numeric:tabular-nums}
.eyebrow{
  font-family:"IBM Plex Mono",monospace; font-size:11.5px; font-weight:500;
  letter-spacing:.13em; text-transform:uppercase; color:var(--ink-3);
}
.lede{font-size:18px; color:var(--ink-2); line-height:1.62}
.muted{color:var(--ink-2)}
.small{font-size:14px}

/* ---- header ---- */
header{padding:64px 0 40px; border-bottom:1px solid var(--line)}
.titlerow{display:flex; flex-wrap:wrap; gap:20px 32px; align-items:flex-end; justify-content:space-between}
.meta{display:flex; flex-direction:column; gap:3px; font-size:13px; color:var(--ink-3)}
.meta b{color:var(--ink-2); font-weight:500}

/* ---- sections ---- */
section{padding:52px 0; border-bottom:1px solid var(--line)}
section:last-of-type{border-bottom:0}
.sec-head{display:flex; flex-direction:column; gap:6px; margin-bottom:26px}
.stack{display:flex; flex-direction:column; gap:18px}

/* ---- result banner ---- */
.headline{
  display:grid; grid-template-columns:repeat(auto-fit,minmax(132px,1fr)); gap:1px;
  background:var(--line); border:1px solid var(--line); border-radius:3px; overflow:hidden;
}
.stat{background:var(--surface); padding:16px 18px; display:flex; flex-direction:column; gap:3px}
.stat .k{font-family:"IBM Plex Mono",monospace; font-size:10.5px; letter-spacing:.11em;
  text-transform:uppercase; color:var(--ink-3)}
.stat .v{font-family:"IBM Plex Mono",monospace; font-size:26px; font-weight:500;
  letter-spacing:-.02em; font-variant-numeric:tabular-nums}
.stat .v.good{color:var(--good)} .stat .v.warn{color:var(--warn)}

/* ---- structure grid ---- */
.structs{display:grid; grid-template-columns:repeat(auto-fit,minmax(288px,1fr)); gap:16px; margin-top:8px}
.sv{border:1px solid var(--line); border-radius:3px; background:var(--surface); overflow:hidden;
  display:flex; flex-direction:column}
.sv-top{padding:13px 15px 11px; border-bottom:1px solid var(--line); display:flex;
  flex-direction:column; gap:2px}
.sv-name{font-size:14.5px; font-weight:600}
.sv-sub{font-size:12.5px; color:var(--ink-3); line-height:1.45}
.viewer{position:relative; height:270px; background:var(--surface-2)}
.sv-bot{padding:11px 15px; border-top:1px solid var(--line); display:flex; gap:18px;
  align-items:baseline; flex-wrap:wrap}
.sv-bot .n{font-family:"IBM Plex Mono",monospace; font-size:19px; font-weight:500;
  font-variant-numeric:tabular-nums}
.sv-bot .lbl{font-family:"IBM Plex Mono",monospace; font-size:10.5px; letter-spacing:.1em;
  text-transform:uppercase; color:var(--ink-3)}
.verdict{margin-left:auto; font-size:12px; font-weight:500; padding:3px 9px; border-radius:2px}
.v-pass{background:var(--good-soft); color:var(--good)}
.v-fail{background:var(--crit-soft); color:var(--crit)}
.v-warn{background:var(--warn-soft); color:var(--warn)}
.legend{display:flex; gap:16px; flex-wrap:wrap; margin-top:14px; font-size:12.5px; color:var(--ink-2)}
.legend span{display:inline-flex; align-items:center; gap:6px}
.sw{width:22px; height:9px; border-radius:1px; display:inline-block}
.hint{font-size:12.5px; color:var(--ink-3); margin-top:10px}

/* ---- tables ---- */
.tw{overflow-x:auto; border:1px solid var(--line); border-radius:3px; background:var(--surface)}
table{border-collapse:collapse; width:100%; font-size:14.5px}
th,td{padding:10px 15px; text-align:left; border-bottom:1px solid var(--line); white-space:nowrap}
thead th{font-family:"IBM Plex Mono",monospace; font-size:10.5px; letter-spacing:.1em;
  text-transform:uppercase; color:var(--ink-3); font-weight:500; background:var(--surface-2)}
tbody tr:last-child td{border-bottom:0}
td.num{font-family:"IBM Plex Mono",monospace; font-variant-numeric:tabular-nums; text-align:right}
tr.hl td{background:var(--accent-soft)}
td .bar{display:inline-block; height:7px; border-radius:1px; background:var(--accent);
  vertical-align:middle; margin-right:8px}

/* ---- callouts ---- */
.note{border-left:2px solid var(--line-2); padding:2px 0 2px 18px; color:var(--ink-2); font-size:15px}
.note.crit{border-color:var(--crit)}
.note.warn{border-color:var(--warn)}
.note.good{border-color:var(--good)}
.note b{color:var(--ink)}

/* ---- sequence ---- */
.seq{font-family:"IBM Plex Mono",monospace; font-size:13px; line-height:1.85;
  word-break:break-all; background:var(--surface); border:1px solid var(--line);
  border-radius:3px; padding:14px 16px; color:var(--ink-2)}

/* ---- misdirection list ---- */
.errs{display:flex; flex-direction:column; gap:0; border:1px solid var(--line);
  border-radius:3px; background:var(--surface); overflow:hidden}
.err{display:grid; grid-template-columns:auto 1fr; gap:0 16px; padding:17px 19px;
  border-bottom:1px solid var(--line)}
.err:last-child{border-bottom:0}
.err .tag{font-family:"IBM Plex Mono",monospace; font-size:11px; letter-spacing:.08em;
  text-transform:uppercase; color:var(--crit); padding-top:3px; white-space:nowrap}
.err h3{margin-bottom:3px}
.err p{font-size:14.5px; color:var(--ink-2)}
.err code{font-family:"IBM Plex Mono",monospace; font-size:13px; background:var(--surface-2);
  padding:1px 5px; border-radius:2px; color:var(--ink)}

/* ---- chart ---- */
.chart{border:1px solid var(--line); border-radius:3px; background:var(--surface); padding:18px}
svg{display:block; width:100%; height:auto; overflow:visible}
.ax{stroke:var(--line-2); stroke-width:1}
.gr{stroke:var(--line); stroke-width:1}
.tk{font-family:"IBM Plex Mono",monospace; font-size:10.5px; fill:var(--ink-3)}
.sl{font-family:"IBM Plex Sans",sans-serif; font-size:12px; font-weight:500}

footer{padding:40px 0 0; color:var(--ink-3); font-size:13px}
@media (prefers-reduced-motion:reduce){*{animation:none!important; transition:none!important}}
</style>

<div class="wrap">

<header>
  <div class="titlerow">
    <div class="col" style="flex:1 1 420px">
      <div class="eyebrow" style="margin-bottom:12px">De novo binder campaign · Report</div>
      <h1>IL-7Rα Binder Campaign</h1>
      <p class="lede" style="margin-top:16px">
        Fifty RFdiffusion backbones against the interleukin-7 receptor α ectodomain.
        One validated candidate at ipTM&nbsp;0.88 — and two null results that removed
        a pipeline stage and reframed what the pipeline was optimising.
      </p>
    </div>
    <div class="meta">
      <div><b>Target</b> PDB 3DI3 chain B</div>
      <div><b>Residues</b> 17–209 (193 aa)</div>
      <div><b>Hotspots</b> S31 · K77 · K138 · Y192</div>
      <div><b>Dates</b> 10–13 Sep 2026</div>
      <div><b>Hardware</b> RTX 3060 laptop, 6 GB</div>
    </div>
  </div>
</header>

<section>
  <div class="sec-head col">
    <div class="eyebrow">Result</div>
    <h2>il7ra_binder_6__s2</h2>
    <p class="muted">Validated under a 5-model ensemble with CPU Amber relaxation. Four of five
    models agree within 0.03 ipTM; one outlier at 0.31.</p>
  </div>

  <div class="headline">
    <div class="stat"><span class="k">ipTM</span><span class="v good">0.88</span></div>
    <div class="stat"><span class="k">pTM</span><span class="v">0.89</span></div>
    <div class="stat"><span class="k">mean pLDDT</span><span class="v">94.8</span></div>
    <div class="stat"><span class="k">interface PAE</span><span class="v">5.8 Å</span></div>
    <div class="stat"><span class="k">length</span><span class="v">83 aa</span></div>
    <div class="stat"><span class="k">off-target</span><span class="v warn">0.61</span></div>
  </div>

  <div class="stack col" style="margin-top:24px">
    <div class="seq">MEKKEEIKKLLEETEKRMEKVCKEAGEKGNKELIDKCLEARQEVLGCFENASYYIKKGDLEKAMEEAKKAQEIVDKLEEIVKK</div>
    <p class="note warn"><b>The off-target score is the open flaw.</b> This binder scores 0.61
    against PD-L1, an unrelated protein — only 1.4× discrimination. The selection criterion
    (ipTM minus a composition penalty) contains no specificity term, so it optimised for
    “binds” and got exactly that.</p>
    <p class="note"><b>This is a scoring result, not a binding result.</b> ipTM is a structure-prediction
    confidence metric. No binding assay has been run on any sequence in this campaign.</p>
  </div>
</section>

<section>
  <div class="sec-head col">
    <div class="eyebrow">Controls</div>
    <h2>Is the score real?</h2>
    <p class="muted">The same 83-residue binder under three conditions. If the score came from
    composition or from AlphaFold's familiarity with the target, the scramble would hold up.
    It does not.</p>
  </div>

  <div class="structs">
    <div class="sv">
      <div class="sv-top">
        <div class="sv-name">On-target</div>
        <div class="sv-sub">Binder + IL-7Rα · relaxed, rank 1 of 5 models</div>
      </div>
      <div class="viewer" id="v-ontarget"></div>
      <div class="sv-bot">
        <span><span class="n" style="color:var(--good)">0.88</span> <span class="lbl">ipTM</span></span>
        <span class="verdict v-pass">designed</span>
      </div>
    </div>
    <div class="sv">
      <div class="sv-top">
        <div class="sv-name">Scrambled</div>
        <div class="sv-sub">Identical residue composition, order destroyed</div>
      </div>
      <div class="viewer" id="v-scrambled"></div>
      <div class="sv-bot">
        <span><span class="n" style="color:var(--crit)">0.26</span> <span class="lbl">ipTM</span></span>
        <span class="verdict v-pass">control passes</span>
      </div>
    </div>
    <div class="sv">
      <div class="sv-top">
        <div class="sv-name">Off-target</div>
        <div class="sv-sub">Same binder against PD-L1</div>
      </div>
      <div class="viewer" id="v-offtarget"></div>
      <div class="sv-bot">
        <span><span class="n" style="color:var(--warn)">0.61</span> <span class="lbl">ipTM</span></span>
        <span class="verdict v-warn">weak specificity</span>
      </div>
    </div>
  </div>

  <div class="legend">
    <span><i class="sw" style="background:var(--plddt-vhigh)"></i> pLDDT &gt; 90</span>
    <span><i class="sw" style="background:var(--plddt-high)"></i> 70–90</span>
    <span><i class="sw" style="background:var(--plddt-low)"></i> 50–70</span>
    <span><i class="sw" style="background:var(--plddt-vlow)"></i> &lt; 50</span>
  </div>
  <p class="hint">Coloured by per-residue pLDDT, the AlphaFold confidence convention. Drag to rotate, scroll to zoom.</p>

  <div class="stack col" style="margin-top:26px">
    <p class="note good"><b>Scramble result settles it.</b> Same 83 residues, same alanine
    content, shuffled order → 0.26. A generic alanine-rich helix scored 0.12 against the same
    target. Sequence <em>order</em> carries the signal, so the score is not a composition artifact.</p>
  </div>
</section>

<section>
  <div class="sec-head col">
    <div class="eyebrow">Null result</div>
    <h2>The fast screen does nothing</h2>
    <p class="muted">Standard practice screens designs with single-sequence ColabFold, then
    rescores the top ranks at high accuracy. We folded the top ten <em>and</em> the bottom ten.</p>
  </div>

  <div class="chart col" style="max-width:none">
    <svg viewBox="0 0 720 200" role="img" aria-label="High-accuracy ipTM for top-ten versus bottom-ten fast-screen designs. The distributions overlap almost exactly.">
      <line class="ax" x1="120" y1="170" x2="700" y2="170"></line>
      <g>
        <line class="gr" x1="120" y1="170" x2="120" y2="30"></line>
        <line class="gr" x1="265" y1="170" x2="265" y2="30"></line>
        <line class="gr" x1="410" y1="170" x2="410" y2="30"></line>
        <line class="gr" x1="555" y1="170" x2="555" y2="30"></line>
        <line class="gr" x1="700" y1="170" x2="700" y2="30"></line>
      </g>
      <text class="tk" x="120" y="188" text-anchor="middle">0.0</text>
      <text class="tk" x="265" y="188" text-anchor="middle">0.25</text>
      <text class="tk" x="410" y="188" text-anchor="middle">0.50</text>
      <text class="tk" x="555" y="188" text-anchor="middle">0.75</text>
      <text class="tk" x="700" y="188" text-anchor="middle">1.0</text>
      <text class="tk" x="410" y="16" text-anchor="middle">high-accuracy ipTM</text>

      <text class="sl" x="112" y="64" text-anchor="end" fill="var(--ink)">top 10</text>
      <text class="tk" x="112" y="79" text-anchor="end">fast 0.10–0.12</text>
      <g fill="var(--accent)" opacity="0.75">
        <circle cx="618" cy="58" r="6"></circle><circle cx="613" cy="58" r="6"></circle>
        <circle cx="607" cy="58" r="6"></circle><circle cx="601" cy="58" r="6"></circle>
        <circle cx="596" cy="58" r="6"></circle><circle cx="555" cy="58" r="6"></circle>
        <circle cx="532" cy="58" r="6"></circle><circle cx="404" cy="58" r="6"></circle>
        <circle cx="294" cy="58" r="6"></circle><circle cx="265" cy="58" r="6"></circle>
      </g>
      <line x1="508" y1="44" x2="508" y2="72" stroke="var(--ink)" stroke-width="2"></line>
      <text class="tk" x="508" y="38" text-anchor="middle" fill="var(--ink)">mean 0.670</text>

      <text class="sl" x="112" y="134" text-anchor="end" fill="var(--ink)">bottom 10</text>
      <text class="tk" x="112" y="149" text-anchor="end">fast 0.06–0.08</text>
      <g fill="var(--warn)" opacity="0.75">
        <circle cx="625" cy="128" r="6"></circle><circle cx="625" cy="128" r="6"></circle>
        <circle cx="601" cy="128" r="6"></circle><circle cx="596" cy="128" r="6"></circle>
        <circle cx="590" cy="128" r="6"></circle><circle cx="561" cy="128" r="6"></circle>
        <circle cx="485" cy="128" r="6"></circle><circle cx="416" cy="128" r="6"></circle>
        <circle cx="311" cy="128" r="6"></circle><circle cx="282" cy="128" r="6"></circle>
      </g>
      <line x1="509" y1="114" x2="509" y2="142" stroke="var(--ink)" stroke-width="2"></line>
      <text class="tk" x="509" y="108" text-anchor="middle" fill="var(--ink)">mean 0.671</text>
    </svg>
  </div>

  <div class="stack col" style="margin-top:22px">
    <p>Difference in means: <b class="mono">0.001</b>. Mann–Whitney <b class="mono">p = 0.970</b>.
    The two groups are indistinguishable.</p>
    <p class="note crit"><b>The screen actively discards good designs.</b>
    <span class="mono">pair_033</span> ranked <b>last of fifty</b> on the fast screen and scored
    <b class="mono">0.870</b> at high accuracy — tied for best in the campaign. The stage costs
    roughly two hours per campaign and returns noise.</p>
  </div>
</section>

<section>
  <div class="sec-head col">
    <div class="eyebrow">Variance study · 168 folds · 28 backbones · 3 targets</div>
    <h2>The sequence is in control, not the backbone</h2>
    <p class="muted">Holding geometry fixed and redrawing the sequence six times, then decomposing
    where the ipTM variance actually lives.</p>
  </div>

  <div class="tw col" style="max-width:none">
    <table>
      <thead><tr>
        <th>Target</th><th>Backbones</th><th style="text-align:right">Between-backbone</th>
        <th style="text-align:right">Within-backbone (sequence)</th><th style="text-align:right">Mean spread</th>
      </tr></thead>
      <tbody>
        <tr><td>IL-7Rα</td><td class="num">10</td><td class="num">41.0%</td>
          <td class="num"><span class="bar" style="width:59px"></span>59.0%</td><td class="num">0.228</td></tr>
        <tr><td>PD-L1</td><td class="num">8</td><td class="num">28.2%</td>
          <td class="num"><span class="bar" style="width:72px"></span>71.8%</td><td class="num">0.302</td></tr>
        <tr class="hl"><td>SARS-CoV-2 RBD</td><td class="num">10</td><td class="num">21.7%</td>
          <td class="num"><span class="bar" style="width:78px"></span>78.3%</td><td class="num">0.540</td></tr>
      </tbody>
    </table>
  </div>
  <p class="hint">Three independent targets, consistent direction. A protocol that generates many
  backbones and one sequence each is optimising the minor term.</p>

  <div class="stack col" style="margin-top:30px">
    <h3>What you gain by drawing more sequences</h3>
    <p class="muted small">Bootstrapped expected ipTM from the best of N draws, same backbones throughout.</p>
  </div>

  <div class="chart col" style="max-width:none; margin-top:14px">
    <svg viewBox="0 0 720 260" role="img" aria-label="Expected ipTM from best-of-N sequence draws, for three targets. RBD rises steeply from 0.386 at one draw to 0.644 at six; IL-7Ra and PD-L1 rise modestly from about 0.78 to about 0.86.">
      <g>
        <line class="gr" x1="90" y1="215" x2="690" y2="215"></line>
        <line class="gr" x1="90" y1="170" x2="690" y2="170"></line>
        <line class="gr" x1="90" y1="125" x2="690" y2="125"></line>
        <line class="gr" x1="90" y1="80"  x2="690" y2="80"></line>
        <line class="gr" x1="90" y1="35"  x2="690" y2="35"></line>
      </g>
      <line class="ax" x1="90" y1="215" x2="690" y2="215"></line>
      <text class="tk" x="82" y="219" text-anchor="end">0.3</text>
      <text class="tk" x="82" y="174" text-anchor="end">0.45</text>
      <text class="tk" x="82" y="129" text-anchor="end">0.6</text>
      <text class="tk" x="82" y="84"  text-anchor="end">0.75</text>
      <text class="tk" x="82" y="39"  text-anchor="end">0.9</text>
      <text class="tk" x="20" y="125" text-anchor="middle" transform="rotate(-90 20 125)">expected ipTM</text>
      <g class="tk" text-anchor="middle">
        <text x="90" y="236">1</text><text x="210" y="236">2</text><text x="330" y="236">3</text>
        <text x="450" y="236">4</text><text x="570" y="236">5</text><text x="690" y="236">6</text>
      </g>
      <text class="tk" x="390" y="254" text-anchor="middle">sequence draws per backbone (N)</text>

      <polyline fill="none" stroke="var(--good)" stroke-width="2.5"
        points="90,75 210,63 330,58 450,56 570,55 690,54"></polyline>
      <polyline fill="none" stroke="var(--accent)" stroke-width="2.5"
        points="90,69 210,52 330,45 450,43 570,41 690,41"></polyline>
      <polyline fill="none" stroke="var(--crit)" stroke-width="2.5"
        points="90,189 210,153 330,135 450,124 570,117 690,112"></polyline>

      <circle cx="90" cy="189" r="4" fill="var(--crit)"></circle>
      <circle cx="690" cy="112" r="4" fill="var(--crit)"></circle>
      <text class="sl" x="700" y="116" fill="var(--crit)">RBD</text>
      <text class="sl" x="700" y="45" fill="var(--accent)">PD-L1</text>
      <text class="sl" x="700" y="58" fill="var(--good)">IL-7Rα</text>
      <text class="tk" x="98" y="204" fill="var(--crit)">0.386</text>
      <text class="tk" x="672" y="104" text-anchor="end" fill="var(--crit)">0.644</text>
    </svg>
  </div>

  <div class="stack col" style="margin-top:22px">
    <p>On the one hard target, a single draw averages <b class="mono">0.386</b> while the same
    backbones reach <b class="mono">0.644</b> at six draws — <b>+0.258 with zero new geometry</b>.
    Most of the gain lands by N&nbsp;=&nbsp;3.</p>
    <p class="note"><b>Where this is thin.</b> RBD is the only hard target here, so the direction
    is consistent across three targets but the magnitude rests on n&nbsp;=&nbsp;1. PD-L1 contributes
    8 backbones rather than 10 — the public MSA server throttled the run. The bootstrap resamples
    with replacement from six observed draws, so the N&nbsp;=&nbsp;6 endpoint is biased low.</p>
  </div>
</section>

<section>
  <div class="sec-head col">
    <div class="eyebrow">Corrections</div>
    <h2>Three false starts</h2>
    <p class="muted">Each survived multiple downstream stages without tripping anything. They are
    the most transferable part of this campaign.</p>
  </div>

  <div class="errs col" style="max-width:none">
    <div class="err">
      <div class="tag">wrong target</div>
      <div>
        <h3>PDB 1ILR is not IL-7Rα</h3>
        <p>It is interleukin-1 receptor antagonist — an unrelated protein. The error passed through
        structure download, design, folding and scoring untouched, because every stage accepts any
        valid PDB. Caught by reading the file header. The correct entry is <code>3DI3</code>.</p>
      </div>
    </div>
    <div class="err">
      <div class="tag">fabricated input</div>
      <div>
        <h3>The first hundred “designs” were random strings</h3>
        <p>A generator stub printed <code>[SKIP] Skipping actual RFdiffusion run</code> and
        substituted <code>random.choice(aa)</code> sequences. One reached ipTM 0.615 and advanced
        through three validation stages before the stub was noticed.</p>
      </div>
    </div>
    <div class="err">
      <div class="tag">wrong input file</div>
      <div>
        <h3>Two hours folding the previous campaign</h3>
        <p>A hard-coded benchmark path sent the IL-7Rα run at 158 pairs of the earlier RBD/PD-L1
        data. Fixed by adding explicit <code>--benchmark</code> and <code>--output-dir</code> arguments.</p>
      </div>
    </div>
  </div>

  <p class="note crit col" style="margin-top:20px"><b>The common thread:</b> every one of these
  produced plausible-looking numbers. None was caught by a score looking wrong — they were caught by
  reading headers, reading source, and running controls.</p>
</section>

<section>
  <div class="sec-head col">
    <div class="eyebrow">Standing</div>
    <h2>What is and isn't established</h2>
  </div>
  <div class="col stack">
    <p class="note good"><b>Established.</b> The fast screen carries no selective signal
    (p = 0.970, n = 20). Within-backbone sequence variance exceeds between-backbone variance on all
    three targets. The 0.88 candidate reproduces across a 5-model ensemble and passes a scramble control.</p>
    <p class="note warn"><b>Not established.</b> That the candidate binds anything — no assay has been
    run. That it is specific — off-target is 0.61. That the sampling effect size generalises — one hard
    target. That the effect holds at default sampling settings — measured only with an alanine bias at
    temperature 0.2.</p>
  </div>
</section>

<footer>
  <p>Pipeline: RFdiffusion → ProteinMPNN → ColabFold (AF2-multimer v3) → composition filter →
  selection on shaped reward. Structures rendered from relaxed and unrelaxed rank-1 predictions,
  backbone atoms only.</p>
</footer>

</div>

<script>
const PDB = __PDB_DATA__;
function plddtColor(a){const b=a.b;
  return b>90?"#0053D6":b>70?"#65CBF3":b>50?"#FFDB13":"#FF7D45";}
function mount(id,key){
  const el=document.getElementById(id); if(!el||!window.$3Dmol) return;
  const css=getComputedStyle(document.body).getPropertyValue("--surface-2").trim();
  const v=$3Dmol.createViewer(el,{backgroundColor:css||"#eff2f6"});
  v.addModel(PDB[key],"pdb");
  v.setStyle({},{cartoon:{colorfunc:plddtColor,thickness:0.3,arrows:true}});
  v.zoomTo(); v.render();
  return v;
}
const viewers=[];
window.addEventListener("load",()=>{
  viewers.push(mount("v-ontarget","ontarget"));
  viewers.push(mount("v-scrambled","scrambled"));
  viewers.push(mount("v-offtarget","offtarget"));
});
// keep viewer grounds in step with the host theme
const mq=window.matchMedia("(prefers-color-scheme: dark)");
mq.addEventListener&&mq.addEventListener("change",()=>{
  const css=getComputedStyle(document.body).getPropertyValue("--surface-2").trim();
  viewers.forEach(v=>{if(v){v.setBackgroundColor(css);v.render();}});
});
</script>
"""

if __name__ == "__main__":
    main()

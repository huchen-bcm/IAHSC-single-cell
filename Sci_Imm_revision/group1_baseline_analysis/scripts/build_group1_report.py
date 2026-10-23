import base64

G1_ASSET_DIR = "/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group1_baseline_analysis/assets"

def b64(path):
    with open(path, "rb") as f:
        return base64.b64encode(f.read()).decode("ascii")

g1_heatmap = b64(f"{G1_ASSET_DIR}/four_condition_heatmap.png")
g1_phase = b64(f"{G1_ASSET_DIR}/cellcycle_phase_composition.png")
g1_rankrank = b64(f"{G1_ASSET_DIR}/rank_rank_infection_response.png")

TEMPLATE = r"""<title>Group 1 — Steady-State Baseline (IA-HSC Revision)</title>
<style>
  :root {
    color-scheme: light;
    --bg: #f9f9f7;
    --surface: #ffffff;
    --surface-2: #f3f2ee;
    --ink: #0b0b0b;
    --ink-2: #52514e;
    --ink-3: #898781;
    --border: rgba(11,11,11,0.11);
    --accent: #4a3aa7;
    --accent-soft: #efecfb;
    --accent-ink: #2f2470;
    --hl-bg: rgba(237,161,0,0.15);
    --hl-border: #b8790a;
    --hl-tag: #8a5a00;
    --res-fg: #075c1f; --res-bg: #e3f5df;
    --figure-paper: #fcfcfb;
    --figure-border: #e1e0d9;
    --shadow: 0 1px 2px rgba(11,11,11,0.04), 0 6px 20px rgba(11,11,11,0.05);
    --font-display: "Iowan Old Style","Palatino Linotype",Palatino,"URW Palladio L",Georgia,serif;
    --font-body: -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
    --font-mono: ui-monospace,SFMono-Regular,Menlo,Consolas,"Liberation Mono",monospace;
  }
  @media (prefers-color-scheme: dark) {
    :root:where(:not([data-theme="light"])) {
      color-scheme: dark;
      --bg: #0c0c0d; --surface: #17171b; --surface-2: #1d1d21;
      --ink: #f4f4f2; --ink-2: #c3c2b7; --ink-3: #918f88;
      --border: rgba(255,255,255,0.10);
      --accent: #9c8fef; --accent-soft: #221d3d; --accent-ink: #cabfff;
      --hl-bg: rgba(201,133,0,0.18); --hl-border: #e0a83a; --hl-tag: #f0c374;
      --res-fg: #6fdb8f; --res-bg: rgba(12,163,12,0.14);
      --shadow: none;
    }
  }
  :root[data-theme="dark"] {
    color-scheme: dark;
    --bg: #0c0c0d; --surface: #17171b; --surface-2: #1d1d21;
    --ink: #f4f4f2; --ink-2: #c3c2b7; --ink-3: #918f88;
    --border: rgba(255,255,255,0.10);
    --accent: #9c8fef; --accent-soft: #221d3d; --accent-ink: #cabfff;
    --hl-bg: rgba(201,133,0,0.18); --hl-border: #e0a83a; --hl-tag: #f0c374;
    --res-fg: #6fdb8f; --res-bg: rgba(12,163,12,0.14);
    --shadow: none;
  }
  :root[data-theme="light"] {
    color-scheme: light;
    --bg: #f9f9f7; --surface: #ffffff; --surface-2: #f3f2ee;
    --ink: #0b0b0b; --ink-2: #52514e; --ink-3: #898781;
    --border: rgba(11,11,11,0.11);
    --accent: #4a3aa7; --accent-soft: #efecfb; --accent-ink: #2f2470;
    --hl-bg: rgba(237,161,0,0.15); --hl-border: #b8790a; --hl-tag: #8a5a00;
    --res-fg: #075c1f; --res-bg: #e3f5df;
    --shadow: 0 1px 2px rgba(11,11,11,0.04), 0 6px 20px rgba(11,11,11,0.05);
  }

  * { box-sizing: border-box; }
  html { background: var(--bg); }
  body {
    margin: 0; background: var(--bg); color: var(--ink);
    font-family: var(--font-body); font-size: 16px; line-height: 1.55;
    -webkit-font-smoothing: antialiased;
  }
  ::selection { background: var(--accent-soft); color: var(--accent-ink); }
  a { color: var(--accent-ink); text-decoration-color: var(--border); text-underline-offset: 2px; }

  h1, h2, h3 { font-family: var(--font-display); font-weight: 600; text-wrap: balance; color: var(--ink); }
  code, .mono { font-family: var(--font-mono); }

  .wrap { max-width: 760px; margin: 0 auto; padding: 56px 20px 96px; }

  header.hero { padding-bottom: 28px; border-bottom: 1px solid var(--border); margin-bottom: 36px; }
  .eyebrow {
    font-family: var(--font-mono); font-size: 12px; letter-spacing: 0.09em; text-transform: uppercase;
    color: var(--accent-ink); margin: 0 0 14px; display: flex; gap: 10px; align-items: center;
  }
  .eyebrow .dot { width: 5px; height: 5px; border-radius: 50%; background: var(--accent); display: inline-block; }
  h1.title { font-size: 30px; line-height: 1.2; margin: 0 0 14px; letter-spacing: -0.01em; }
  .dek { font-size: 17px; color: var(--ink-2); max-width: 62ch; margin: 0 0 18px; }
  .meta-row { display: flex; flex-wrap: wrap; gap: 8px 22px; font-size: 13.5px; color: var(--ink-3); }
  .meta-row strong { color: var(--ink-2); font-weight: 600; }

  .pill {
    display: inline-flex; align-items: center; gap: 6px; font-size: 12px; font-weight: 600;
    padding: 3px 10px 3px 8px; border-radius: 999px; white-space: nowrap;
  }
  .pill::before { content: ""; width: 6px; height: 6px; border-radius: 50%; }
  .pill.resolved { color: var(--res-fg); background: var(--res-bg); }
  .pill.resolved::before { background: var(--res-fg); }

  .group-card {
    border: 1px solid var(--border); border-radius: 14px; background: var(--surface);
    box-shadow: var(--shadow); padding: 30px 32px 32px; margin-bottom: 28px;
  }
  .subhead {
    font-size: 11.5px; text-transform: uppercase; letter-spacing: 0.08em; color: var(--ink-3);
    font-weight: 600; margin: 22px 0 10px;
  }
  .subhead:first-of-type { margin-top: 4px; }

  blockquote.quote {
    margin: 0 0 12px; padding: 12px 16px; background: var(--hl-bg); border-left: 3px solid var(--hl-border);
    border-radius: 0 8px 8px 0; font-size: 14.5px; line-height: 1.55; color: var(--ink);
  }
  blockquote.quote .who {
    display: block; font-family: var(--font-mono); font-size: 11px; letter-spacing: 0.05em;
    text-transform: uppercase; color: var(--hl-tag); margin-bottom: 5px; font-weight: 600;
  }
  blockquote.quote p { margin: 0; }
  blockquote.quote b { background: rgba(0,0,0,0.06); padding: 0 3px; border-radius: 3px; }

  .prose p { margin: 0 0 12px; color: var(--ink); }
  .prose p:last-child { margin-bottom: 0; }
  .muted { color: var(--ink-2); }
  .prose code { background: var(--surface-2); padding: 0.5px 5px; border-radius: 4px; font-size: 0.89em; border: 1px solid var(--border); }

  .callout {
    margin: 14px 0; padding: 14px 16px; border-radius: 10px;
    background: var(--accent-soft); border: 1px solid color-mix(in srgb, var(--accent) 28%, transparent);
    font-size: 14.5px;
  }
  .callout .lbl {
    display: block; font-size: 11px; text-transform: uppercase; letter-spacing: 0.07em;
    color: var(--accent-ink); font-weight: 700; margin-bottom: 6px;
  }
  .callout p { margin: 0 0 8px; color: var(--ink); }
  .callout p:last-child { margin-bottom: 0; }
  .callout.warn { border-color: color-mix(in srgb, var(--hl-border) 55%, transparent); background: var(--hl-bg); }
  .callout.warn .lbl { color: var(--hl-tag); }

  figure.fig {
    margin: 18px 0; border-radius: 10px; overflow: hidden;
    border: 1px solid var(--figure-border); background: var(--figure-paper); box-shadow: var(--shadow);
  }
  figure.fig img { display: block; width: 100%; height: auto; }
  figure.fig figcaption {
    padding: 10px 16px 13px; font-size: 13px; color: #52514e; background: var(--figure-paper);
    border-top: 1px solid var(--figure-border); line-height: 1.5;
  }

  .tbl-wrap { overflow-x: auto; margin: 14px 0; border: 1px solid var(--border); border-radius: 10px; }
  table.data { width: 100%; border-collapse: collapse; font-size: 13.5px; min-width: 480px; }
  table.data th {
    text-align: left; padding: 9px 14px; background: var(--surface-2); color: var(--ink-2);
    font-weight: 600; font-size: 12px; text-transform: uppercase; letter-spacing: 0.04em;
    border-bottom: 1px solid var(--border); white-space: nowrap;
  }
  table.data td { padding: 8px 14px; border-bottom: 1px solid var(--border); font-variant-numeric: tabular-nums; white-space: nowrap; }
  table.data tr:last-child td { border-bottom: none; }
  table.data td.cl { font-family: var(--font-mono); font-weight: 600; color: var(--accent-ink); }
  table.data td.ns { color: var(--ink-3); }

  footer {
    margin-top: 40px; padding-top: 20px; border-top: 1px solid var(--border);
    font-size: 12.5px; color: var(--ink-3); display: flex; flex-wrap: wrap; gap: 6px 20px;
  }

  @media (max-width: 600px) {
    h1.title { font-size: 25px; }
    .group-card { padding: 22px 18px 24px; }
    .wrap { padding: 40px 14px 80px; }
  }
</style>

<div class="wrap">

  <header class="hero">
    <p class="eyebrow"><span class="dot"></span>Science Immunology &middot; Revision Round 1 &middot; Group 1 of the bioinformatics review</p>
    <h1 class="title">Is IA-HSC "immune-primed" constitutively, or infection-activated?</h1>
    <p class="dek">Cao <i>et al.</i>, "Infection-activated hematopoietic stem cells harbor central trained immunity." A focused writeup answering the reviewers' steady-state-baseline questions &mdash; the four-condition pathway comparison, cell cycle, and the LT-HSC infection-response cross-check.</p>
    <div class="meta-row">
      <span><strong>Status</strong> <span class="pill resolved">Resolved</span></span>
      <span><strong>Manuscript</strong> aei2812_CombinedPDF_v2.pdf (submitted version)</span>
      <span><strong>Analysis folder</strong> Sci_Imm_revision/group1_baseline_analysis</span>
    </div>
  </header>

  <div class="group-card">
    <p class="subhead">Reviewer quotes</p>
    <blockquote class="quote"><span class="who">Reviewer 1 · point 2</span><p>"RNA-seq in Figure 4 shows increased 'myeloid cell activation' of IA-HSCs compared to LT-HSCs one month post <i>M. avium</i>. (2.1) Is there a difference in cell cycle score between these populations? (2.2) Can these pathways also be compared in IA-HSC and LT-HSC at steady-state (without bacterial exposure)? Specifically, showing relevant pathways across four conditions (LT-HSC at steady-state, IA-HSC at steady-state, LT-HSC post <i>M. avium</i>, IA-HSC post <i>M. avium</i>) would be helpful."</p></blockquote>
    <blockquote class="quote"><span class="who">Reviewer 1 · point 3</span><p>"...the conclusion... 'that IA-HSC represent an immune-primed subset poised for rapid innate immune gene activation', is poorly supported without showing the steady-state level of innate immune gene expression: if IA-HSCs already have higher innate immune gene expression at steady-state, the relative activation may not differ from LT-HSCs."</p></blockquote>
    <blockquote class="quote"><span class="who">Reviewer 2 · major 3, 2nd bullet</span><p>"A comparison of differentially expressed genes (DEGs) between LT-HSCs from uninfected vs. infected mice should be included and discussed in relation to the DEGs presented in Figure 4A. This would help clarify whether the transcriptional changes in IA-HSCs are unique or shared with other HSC populations under infection."</p></blockquote>

    <div class="prose"><p class="muted">Maps to Results section "IA-HSCs undergo transcriptional reprogramming..." (Fig. 4A&ndash;F, fig. S4A&ndash;B in the actual submitted manuscript); text "IA-HSC represent an immune-primed subset poised for rapid innate immune gene activation" (page 12).</p></div>

    <p class="subhead">Confirmed before analysis</p>
    <div class="prose">
      <p>Bulk RNA-seq (Fig. 4A/4B) was only ever run in infected mice &mdash; <i>"we conducted bulk RNA sequencing in IA-HSCs and LT-HSCs from mice that had been infected with M. avium one month prior"</i> &mdash; there's no steady-state bulk arm to fall back on, which is exactly why a scRNA-based approach is needed here.</p>
      <p>Across every notebook in the repo, the only true single-cell DE (<code>singlecell_de</code>, a per-cell t-test &mdash; not pseudobulk/DESeq2) ever run was <b>IA-HSC Infected vs. Naive</b>. The one file that looked like an "IA-HSC vs. other-HSC" scRNA comparison (<code>resultsDEG_HSC.IAHSC_vs_otherHSC.tsv</code>) was actually generated via pseudobulk/DESeq2, confirmed by its <code>baseMean</code>/<code>lfcSE</code>/<code>stat</code> columns, which <code>singlecell_de</code> output never has.</p>
    </div>

    <p class="subhead">What we ran</p>
    <div class="prose"><p>Ran the missing <strong>LT-HSC infected-vs-naive</strong> scRNA contrast (same pipeline, thresholds, and h5ad already used for the published IA-HSC contrast), cross-tabulated it against the existing IA-HSC infected-vs-naive and Fig. 4A bulk IA-HSC-vs-LT-HSC gene lists, and built the exact 4-condition pathway panel R1.2.2 requested &mdash; the primary evidence below, aggregated to n=4 mice/group rather than tested cell-by-cell, to avoid pseudoreplication. Cell-cycle scores were already computed in the h5ad. We also attempted a more direct gene-level scRNA "IA-HSC vs. LT-HSC" DE comparison (Finding 4) &mdash; found unreliable, not used as evidence.</p></div>

    <div class="callout">
      <span class="lbl">Finding 1 &middot; R1.3 confirmed, with a consequence for the text</span>
      <p>For all three named Fig. 4B pathways (type&nbsp;II IFN, chemotaxis, myeloid activation), IA-HSC already scores significantly higher than LT-HSC <b>at steady-state</b> (naive, p&nbsp;=&nbsp;0.003&ndash;0.03), and barely changes with infection within IA-HSC itself (all p&nbsp;&gt;&nbsp;0.65) &mdash; while LT-HSC's own score for the same pathways shifts significantly with infection (p&nbsp;&le;&nbsp;0.0014). These read as a constitutive IA-HSC identity feature, not an infection-triggered activation. The one exception &mdash; a small "innate defense mediator" panel (Irf1/Irf7/Stat1/Gbp6/Gbp7) &mdash; is genuinely infection-induced in both cell types and isn't IA-HSC-specific (LT-HSC's induced level is actually higher).</p>
    </div>

    <figure class="fig">
      <img src="data:image/png;base64,__G1_HEATMAP__" alt="Heatmap of five pathway and cell-cycle scores across four conditions: LT-HSC naive, IA-HSC naive, LT-HSC infected, and IA-HSC infected, showing IA-HSC elevated for interferon, chemotaxis, and myeloid activation scores at both timepoints with little change upon infection, while LT-HSC shifts more with infection." />
      <figcaption>Row-wise z-scored pathway/cell-cycle scores across the four conditions (numbers = mean UCell / cell-cycle score, n=4 mice/group). Columns 1&amp;2 vs. 3&amp;4 compare naive to infected within each cell type.</figcaption>
    </figure>

    <div class="tbl-wrap">
      <table class="data">
        <thead><tr><th>Signature</th><th>IA vs LT, naive</th><th>IA vs LT, infected</th><th>LT: inf. vs naive</th><th>IA: inf. vs naive</th></tr></thead>
        <tbody>
          <tr><td class="cl">Type II IFN</td><td>0.0043</td><td>0.0051</td><td>0.0014</td><td class="ns">0.655</td></tr>
          <tr><td class="cl">Chemotaxis</td><td>&lt;0.0001</td><td>0.0003</td><td>0.0002</td><td class="ns">0.831</td></tr>
          <tr><td class="cl">Myeloid activation</td><td>0.0265</td><td>0.0001</td><td>0.0003</td><td class="ns">0.712</td></tr>
          <tr><td class="cl">Innate defense mediators</td><td class="ns">0.221</td><td>0.0011</td><td>&lt;0.0001</td><td>&lt;0.0001</td></tr>
        </tbody>
      </table>
    </div>

    <div class="callout">
      <span class="lbl">Finding 2 &middot; cell cycle (R1.2.1)</span>
      <p>IA-HSC is ~47&ndash;57% G2M at both timepoints vs. LT-HSC's ~19&ndash;22% &mdash; present already at baseline and essentially flat with infection in IA-HSC, while LT-HSC's G2M fraction <i>drops</i> with infection (becoming more quiescent). Complements (doesn't replace) the manuscript's existing BrdU proliferation data.</p>
    </div>

    <figure class="fig">
      <img src="data:image/png;base64,__G1_PHASE__" alt="Stacked bar chart of cell cycle phase composition (G1, S, G2M) across the four conditions, showing IA-HSC has roughly half its cells in G2M phase at both naive and infected timepoints, compared to about one fifth for LT-HSC." />
      <figcaption>Cell-cycle phase composition (from <code>scanpy.tl.score_genes_cell_cycle</code>, already computed in the object) across the four conditions.</figcaption>
    </figure>

    <div class="callout">
      <span class="lbl">Finding 3 &middot; R2 Major-3: is the Fig. 4A signature just generic infection response?</span>
      <p>LT-HSC has far more infection-responsive genes than IA-HSC (842 vs. 115). Only <b>9.0%</b> of the Fig. 4A "IA-HSC identity" signature (1296 genes) overlaps with LT-HSC's own infection response &mdash; significantly above the 3.1% background rate (2.93&times; enrichment, p&nbsp;=&nbsp;7.2&times;10<sup>&minus;26</sup>), but the other 91% is <i>not</i> explained this way. Reassuring for the paper's core claim, with a small caveat worth stating explicitly.</p>
    </div>

    <figure class="fig">
      <img src="data:image/png;base64,__G1_RANKRANK__" alt="Scatter plot comparing log2 fold change upon infection in LT-HSC versus IA-HSC, with Fig 4A identity genes highlighted in purple and genes significant in both infection contrasts circled in red, showing only weak correlation (Spearman rho 0.15) between the two cell types' infection responses." />
      <figcaption>LT-HSC vs. IA-HSC infection response, gene-by-gene. The 34-gene shared core (red) sits along the diagonal; most Fig. 4A identity genes (purple) don't track either axis strongly.</figcaption>
    </figure>

    <div class="callout warn">
      <span class="lbl">Finding 4 &middot; a more direct approach was tried, and found unreliable &mdash; not used as evidence</span>
      <p>Also ran two new gene-level <code>singlecell_de</code> contrasts, IA-HSC vs. LT-HSC within Naive and within Infected, to build a real steady-state DEG list instead of leaning on pathway scores. Before trusting either, checked the Infected version against the Fig. 4A bulk gold-standard (same cell types, same treatment, just scRNA vs. bulk) &mdash; <b>it doesn't reproduce it</b>: only 17.2% of the 1230 Fig. 4A bulk-significant genes are also significant and concordant in scRNA, 10.4% are significant but flip direction, and 64&ndash;67% of everything scanpy calls "significant" in these two contrasts has an implausible |log2FC|&gt;15 (a near-zero-count artifact expected when testing a 90&ndash;107-cell population cell-by-cell against a 4,000+-cell one). GSEA on these contrasts also failed to independently recover the named pathways. <b>Conclusion: this gene-level comparison isn't reliable enough to use as evidence</b> &mdash; Findings 1&ndash;3 above (aggregated to n=4 mice/group) remain the basis for every conclusion in this writeup.</p>
    </div>

    <p class="subhead">Recommended text fix</p>
    <div class="prose"><p>Replace <i>"IA-HSC represent an immune-primed subset poised for rapid innate immune gene activation"</i> with language describing a constitutive, steady-state-present identity difference rather than an infection-triggered one, e.g.: <i>"IA-HSCs constitutively express elevated type&nbsp;II interferon, leukocyte chemotaxis, and myeloid-activation gene signatures relative to LT-HSCs &mdash; a difference already present at steady-state and largely maintained rather than amplified by infection &mdash; consistent with a stable, lineage-intrinsic identity rather than an infection-triggered activation state."</i></p></div>

    <p class="subhead">Caveat</p>
    <div class="prose"><p>The three pathway gene sets used for Findings 1&ndash;3 are our own best-matching GO terms from re-running GSEA on the bulk contrast (leading-edge genes) &mdash; one reached FDR&lt;0.05, the other two were directionally correct but FDR 0.07&ndash;0.14 in this re-run. This is a gene-set-construction exercise, not a re-confirmation of Fig. 4B's own reported significance; swapping in the original Fig. 4B leading-edge lists, if on file, would directly improve on this.</p></div>
  </div>

  <footer>
    <span>Prepared for internal co-author review &mdash; not for distribution outside the author group.</span>
    <span>Generated 2026-07-21</span>
  </footer>

</div>
"""

html = (TEMPLATE.replace("__G1_HEATMAP__", g1_heatmap)
                 .replace("__G1_PHASE__", g1_phase)
                 .replace("__G1_RANKRANK__", g1_rankrank))

out_path = "/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group1_baseline_analysis/group1_steady_state_report.html"
with open(out_path, "w") as f:
    f.write(html)
print("wrote", out_path, len(html), "bytes")

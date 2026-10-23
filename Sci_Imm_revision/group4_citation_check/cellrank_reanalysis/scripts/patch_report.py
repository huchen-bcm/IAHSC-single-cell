import re

SRC = '/home/huc/.claude/projects/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/tool-results/artifact-717798e7-1788733159-fdce.html'
with open(SRC, 'r', encoding='utf-8') as f:
    html = f.read()

with open('/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/img_v1_b64.txt') as f:
    img_v1_b64 = f.read()

# 1. Update header status pill
old_pill = ('<span><strong>Status</strong> <span class="pill mixed">'
            "All 15 donors run — median/magnitude match, significance doesn't</span></span>")
assert old_pill in html, "pill not found"
new_pill = ('<span><strong>Status</strong> <span class="pill mixed">'
            'Does not reproduce — tested on both CellRank 2.x and period-matched CellRank 1.x software</span></span>')
html = html.replace(old_pill, new_pill)

# 2. Insert PART 6 after PART 5's group-card, before SUMMARY section-head
marker = """  <div class="section-head">
    <span class="num">SUMMARY</span>
    <h2>What's done, what's next</h2>
  </div>"""
assert marker in html, "summary marker not found"

part6 = """  <div class="section-head">
    <span class="num">PART 6</span>
    <h2>Testing software-version drift as the explanation, empirically</h2>
    <p class="sub">Sun et al.'s code was written against CellRank 1.x. The Part 5 run used CellRank 2.0.6 (ported API). Re-ran all 30 samples on a reconstructed CellRank 1.5.1 / scVelo 0.2.5 stack, matching the original code line-for-line, to test whether the version gap explains the non-reproduction.</p>
  </div>

  <div class="group-card">
    <div class="prose">
      <p>Rebuilt the original environment from scratch: CellRank 1.5.1, scVelo 0.2.5, and the exact API calls from Sun et al.'s own scripts (<code>cellrank.tl.kernels</code>, <code>GPCCA.compute_absorption_probabilities</code>, <code>adata.obsm['to_terminal_states']</code>) &mdash; no porting, verbatim. This required resolving several environment conflicts along the way (matplotlib/pandas/numpy pins for scanpy 1.9.1 compatibility) and, notably, a Schur-decomposition solver problem: without <code>petsc4py</code>/<code>slepc4py</code> installed, CellRank silently falls back to a dense (&ldquo;<code>brandts</code>&rdquo;) solver that produced wildly unstable per-donor results (0% to 74% swings). Installing <code>petsc4py</code>/<code>slepc4py</code> via conda-forge initially pulled in a too-new PETSc (3.23) whose Python API had dropped an attribute (<code>ksp.converged</code>) that CellRank 1.5.1's code depends on &mdash; pinning to petsc4py/slepc4py 3.17.x (contemporaneous with CellRank 1.5.1, with openmpi rather than mpich) fixed this and gave numerically stable, properly-normalized fate probabilities (verified: probabilities sum to 1.0 per cell, std &lt; 1.2e-7).</p>
    </div>

    <p class="subhead">Full 30-sample batch, CellRank 1.5.1 + properly-configured PETSc/SLEPc</p>
    <div class="prose">
      <p>All 30 samples (15 donors &times; 2 timepoints) completed with 0 failures.</p>
    </div>
    <figure class="fig">
      <img src="data:image/png;base64,__IMG_V1_PANELC__" alt="Panel C reproduction, CellRank 1.5.1, CMP_b bias at D90, full donor set">
      <figcaption>% of each donor's HSCs with CMP_b (= CMP c2) as top predicted fate, at D90/Tm3, all 15 donors, using the period-matched CellRank 1.5.1 stack.</figcaption>
    </figure>
    <div class="callout resolve">
      <span class="lbl">Result: does not reproduce, and version drift is not the explanation</span>
      <p><b>BCG median = 2.7% (n=11), CTL median = 2.6% (n=4), Mann-Whitney p = 0.84.</b> This is not significant &mdash; and it is a <i>worse</i> match to the paper's reported values (BCG &asymp;9&ndash;10%) than the CellRank 2.x run in Part 5 was. Using the original, period-matched software stack did not recover the paper's result; if anything the BCG estimate moved further from the paper's reported value. <b>Bottom line, stated plainly: with the same underlying data and the authors' own published code, run first on a modern CellRank/scVelo stack and then on a reconstructed period-matched CellRank 1.5.1/scVelo 0.2.5/PETSc 3.17 stack, this reproduction does not recover the paper's reported significant BCG vs. CTL difference in Day-90 CMP_b fate bias, under either software configuration.</b></p>
    </div>
  </div>

"""
html = html.replace(marker, part6 + marker)
html = html.replace('__IMG_V1_PANELC__', img_v1_b64)

# 3. Update SUMMARY table's Figure 3C row and add a CellRank 1.x row
old_row = ('<tr><td class="wrap-cell">Figure 3C reproduction (% CMP_b/CMP-c2-biased HSCs, D90, BCG vs. CTL)</td>'
           '<td class="wrap-cell">All 15 donors run (11 BCG / 4 CTL). <b>Medians match the paper closely</b> '
           "(BCG 9.4% vs. paper's &asymp;9&ndash;10%; CTL 2.6% vs. &asymp;3%). <b>Wilcoxon significance does not reproduce</b> "
           '(p=0.66; p=0.17 excluding one CTL outlier donor, S13) &mdash; see Part 5</td></tr>')
assert old_row in html, "old summary row not found"
new_rows = old_row + ('\n          <tr><td class="wrap-cell">Same, re-run on period-matched CellRank 1.5.1/scVelo 0.2.5/PETSc 3.17 stack</td>'
                       '<td class="wrap-cell">All 15 donors run. <b>BCG median 2.7%, CTL median 2.6%, p=0.84</b> '
                       '&mdash; does not reproduce, and is a worse match to the paper than the CellRank 2.x run &mdash; see Part 6</td></tr>')
html = html.replace(old_row, new_rows)

# 4. Replace the stale "Next step" callout with a plain closing statement
old_next = """    <div class="callout">
      <span class="lbl">Next step</span>
      <p>Run the remaining candidate donor pairs (2 more BCG, 2 CTL) through the same pipeline, then reproduce Sun's own aggregation step (<code>02_terminal_fate_analysis.R</code>: per-donor Td0&rarr;Tm3 shift, Wilcoxon BCG vs. CTL) to see whether c2/c3/c4/c7's fate bias shifts with vaccination &mdash; and whether that shift differs between BCG and CTL, which is the actual claim in Sun's Figure 3.</p>
    </div>"""
assert old_next in html, "next-step callout not found"
new_next = """    <div class="callout">
      <span class="lbl">Bottom line</span>
      <p>The pipeline, data, sample mapping, and terminal-fate logic were all validated against Sun et al.'s own code and confirmed correct. Two independent full-cohort runs &mdash; one on a modern CellRank 2.x stack, one on a reconstructed period-matched CellRank 1.5.1 stack &mdash; both land on the same outcome: the reported BCG vs. CTL significance in Figure 3C does not reproduce from this data and code.</p>
    </div>"""
html = html.replace(old_next, new_next)

OUT = '/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/group4_cellrank_report_v2.html'
with open(OUT, 'w', encoding='utf-8') as f:
    f.write(html)
print("wrote", OUT, len(html))

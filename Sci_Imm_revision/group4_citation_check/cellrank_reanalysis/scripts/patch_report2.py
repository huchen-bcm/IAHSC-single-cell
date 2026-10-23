SRC = '/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/group4_cellrank_report_v2.html'
with open(SRC, 'r', encoding='utf-8') as f:
    html = f.read()

with open('/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/img_fig3d_v2_b64.txt') as f:
    img_b64 = f.read()

# 1. Part 4 sub-header
old_sub = ('<p class="sub">Fig. 3D is a PHATE map, not a boxplot: HSCs colored by predicted terminal fate, '
           'placebo (top) vs. BCG (bottom). Ran one more donor (S18, CTL/placebo, both timepoints) '
           'specifically to fill out this structure.</p>')
assert old_sub in html
new_sub = ('<p class="sub">Fig. 3D is a PHATE map, not a boxplot: HSCs colored by predicted terminal fate, '
           'placebo (top) vs. BCG (bottom). <b>Correction: the two columns are not Td0 vs. Tm3</b> &mdash; '
           'the panel headers in the manuscript read &ldquo;D90&rdquo; for both columns in each row, so the '
           'two columns are two different example donors at the same timepoint, not two timepoints of one donor. '
           'An earlier version of this reproduction misread the columns as Td0/Tm3; corrected below.</p>')
html = html.replace(old_sub, new_sub)

# 2. figure caption
old_cap = ('<figcaption>Reproduction of Fig. 3D\'s structure: rows = placebo/CTL (donor S18) vs. BCG (donor S20), '
           'columns = Td0 vs. Tm3, colored by predicted terminal fate.</figcaption>')
assert old_cap in html
new_cap = ('<figcaption>Corrected reproduction of Fig. 3D\'s structure: rows = CTL vs. BCG, both at D90 &mdash; '
           'columns = two different example donors per group (CTL: S11, S18; BCG: S1, S20 &mdash; picked as '
           'plausible visual matches to the manuscript panel, not confirmed to be the same donors), colored '
           'by predicted terminal fate.</figcaption>')
html = html.replace(old_cap, new_cap)

# 3. swap the embedded image (find the <figure class="fig"> block right before the old caption)
# locate the img tag preceding the caption we just replaced
marker_before_img = '<figure class="fig">\n      <img src="data:image/png;base64,'
idx = html.find(marker_before_img)
assert idx != -1
start = idx + len(marker_before_img)
end = html.find('"', start)
html = html[:start] + img_b64 + html[end:]

# 4. callout text
old_callout = ('<span class="lbl">Reading &mdash; structurally consistent, but not yet a real comparison</span>\n'
               '      <p>The PHATE shape itself looks like a genuine developmental trajectory &mdash; a curved arc '
               'with CMPs/GMPs (green) concentrated at one end and MLPs/PreBNK (blue) + MEPs (orange) mixed at the '
               'other, consistent across all four panels. This is a meaningful sanity check that velocity&rarr;'
               'CellRank&rarr;fate-probability&rarr;PHATE-projection is behaving correctly end to end. <b>But this is '
               'one CTL donor vs. one BCG donor, not yet a real placebo-vs-BCG comparison</b> &mdash; any visual '
               'difference between rows here (e.g., the denser orange/MEP patch in BCG&ndash;Td0) is just as likely to '
               'be donor-to-donor variation as a vaccination effect, and shouldn\'t be read as a finding yet.</p>')
assert old_callout in html
new_callout = ('<span class="lbl">Reading &mdash; structurally consistent, but still not a real comparison, and '
               'donor identity is a guess</span>\n'
               '      <p>The PHATE shape itself looks like a genuine developmental trajectory &mdash; a curved arc '
               'with CMPs/GMPs (green) concentrated at one end and MLPs/PreBNK (blue) + MEPs (orange) mixed at the '
               'other, consistent across all four panels. This is a meaningful sanity check that velocity&rarr;'
               'CellRank&rarr;fate-probability&rarr;PHATE-projection is behaving correctly end to end. The two BCG '
               'panels (and the two CTL panels) do resemble each other more than they resemble the opposite group, '
               'consistent with the manuscript\'s point. <b>But which specific donors the manuscript actually used is '
               'not stated anywhere in the paper</b> &mdash; S11/S18/S1/S20 were picked from all 15 available donors '
               'by eye, as the closest match in color composition to the manuscript screenshot; this is a plausible '
               'illustration, not a confirmed identification, and still isn\'t a substitute for the full 15-donor '
               'statistical comparison in Part 5&ndash;6.</p>')
html = html.replace(old_callout, new_callout)

OUT = '/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/group4_cellrank_report_v3.html'
with open(OUT, 'w', encoding='utf-8') as f:
    f.write(html)
print('wrote', OUT, len(html))

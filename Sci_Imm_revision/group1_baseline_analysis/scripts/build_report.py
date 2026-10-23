import json
with open('/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/img_b64.json') as f:
    IMG = json.load(f)

html = open('/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/report_template.html').read()

html = html.replace('__IMG_FOUR_COND__', IMG['four_condition_heatmap.png'])
html = html.replace('__IMG_CELLCYCLE__', IMG['cellcycle_phase_composition.png'])
html = html.replace('__IMG_RANKRANK__', IMG['rank_rank_infection_response.png'])
html = html.replace('__IMG_ROBUST__', IMG['robustness_alt_genesets_barplot.png'])
html = html.replace('__IMG_PROLIF__', IMG['proliferation_multiangle_summary.png'])

path = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group1_baseline_analysis/group1_steady_state_report.html'
with open(path, 'w') as f:
    f.write(html)
print("written", len(html), "bytes")

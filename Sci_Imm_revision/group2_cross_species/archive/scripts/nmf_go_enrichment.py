import json
import gseapy as gp
import pandas as pd

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species/nmf_metaprogram'

mt = json.load(open(f'{OUT}/mouse_mp_top_genes.json'))
ht = json.load(open(f'{OUT}/human_mp_top_genes.json'))
all_programs = {f'mouse_{k}': v for k, v in mt.items()}
all_programs.update({f'human_{k}': v for k, v in ht.items()})

results = {}
top_terms = {}
for name, genes in all_programs.items():
    try:
        enr = gp.enrichr(gene_list=list(genes), gene_sets=['GO_Biological_Process_2025'],
                          organism='human', outdir=None)
        df = enr.results
        df = df[df['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value')
        results[name] = df
        top_terms[name] = df['Term'].head(5).tolist()
        print(name, '-> top terms:', top_terms[name][:3])
    except Exception as e:
        print(name, 'FAILED', e)
        top_terms[name] = []

with open(f'{OUT}/mp_go_top_terms.json', 'w') as f:
    json.dump(top_terms, f, indent=1)

full = pd.concat([df.assign(program=name) for name, df in results.items()], ignore_index=True)
full.to_csv(f'{OUT}/mp_go_enrichment_full.csv', index=False)
print('done')

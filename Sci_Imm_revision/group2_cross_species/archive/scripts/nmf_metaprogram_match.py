"""
Cross-species NMF meta-program matching pipeline.
Run AFTER cNMF consensus() has been called for both mouse_cnmf and human_cnmf at chosen K's.
Steps:
  1. Load gene_spectra_score matrices (genes x programs) for mouse and human.
  2. For each program, take top-N genes by score -> gene set.
  3. Jaccard similarity matrix between all mouse programs x all human programs.
  4. Match programs with Jaccard > threshold (or best mutual match) -> cross-species meta-programs.
  5. GO/Reactome enrichment (gseapy Enrichr, human gene symbols) per matched meta-program's gene union.
  6. Score per-cell usage (already in consensus usage matrix) -> compare mean usage by
     mouse HSC_subclass and human hsc_clusters_lab for each matched MP.
"""
import numpy as np
import pandas as pd
import json

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species/nmf_metaprogram'
TOP_N = 50


def load_spectra(name, k, dt):
    f = f'{OUT}/{name}/{name}.gene_spectra_score.k_{k}.dt_{dt}.txt'
    df = pd.read_csv(f, sep='\t', index_col=0)
    return df  # programs x genes


def load_usage(name, k, dt):
    f = f'{OUT}/{name}/{name}.usages.k_{k}.dt_{dt}.consensus.txt'
    df = pd.read_csv(f, sep='\t', index_col=0)
    return df  # cells x programs


DT = '0_2'


def top_genes(spectra_row, n=TOP_N):
    return set(spectra_row.sort_values(ascending=False).head(n).index)


def jaccard(a, b):
    return len(a & b) / len(a | b)


def run(mouse_k, mouse_dt, human_k, human_dt):
    mouse_spectra = load_spectra('mouse_cnmf', mouse_k, mouse_dt)
    human_spectra = load_spectra('human_cnmf', human_k, human_dt)
    mouse_usage = load_usage('mouse_cnmf', mouse_k, mouse_dt)
    human_usage = load_usage('human_cnmf', human_k, human_dt)

    mouse_top = {f'mouseMP{i}': top_genes(mouse_spectra.loc[i]) for i in mouse_spectra.index}
    human_top = {f'humanMP{i}': top_genes(human_spectra.loc[i]) for i in human_spectra.index}

    jmat = pd.DataFrame(index=mouse_top.keys(), columns=human_top.keys(), dtype=float)
    for mk, mg in mouse_top.items():
        for hk, hg in human_top.items():
            jmat.loc[mk, hk] = jaccard(mg, hg)
    jmat.to_csv(f'{OUT}/jaccard_matrix.csv')
    print(jmat.round(3))

    # best mutual match
    matches = []
    for mk in jmat.index:
        best_h = jmat.loc[mk].idxmax()
        best_h_val = jmat.loc[mk, best_h]
        best_m_for_h = jmat[best_h].idxmax()
        mutual = (best_m_for_h == mk)
        matches.append(dict(mouse_mp=mk, human_mp=best_h, jaccard=best_h_val, mutual_best=mutual))
    matches = pd.DataFrame(matches).sort_values('jaccard', ascending=False)
    matches.to_csv(f'{OUT}/mp_matches.csv', index=False)
    print(matches)

    # gene sets for matched pairs (union) for GO enrichment
    match_genesets = {}
    for _, row in matches.iterrows():
        pair_name = f"{row['mouse_mp']}__{row['human_mp']}"
        genes = sorted(mouse_top[row['mouse_mp']] | human_top[row['human_mp']])
        match_genesets[pair_name] = genes
    with open(f'{OUT}/matched_mp_genesets.json', 'w') as f:
        json.dump(match_genesets, f, indent=1)

    mouse_top_json = {k: sorted(v) for k, v in mouse_top.items()}
    human_top_json = {k: sorted(v) for k, v in human_top.items()}
    with open(f'{OUT}/mouse_mp_top_genes.json', 'w') as f:
        json.dump(mouse_top_json, f, indent=1)
    with open(f'{OUT}/human_mp_top_genes.json', 'w') as f:
        json.dump(human_top_json, f, indent=1)

    mouse_usage.columns = [f'mouseMP{c}' for c in mouse_usage.columns]
    human_usage.columns = [f'humanMP{c}' for c in human_usage.columns]
    mouse_usage.to_csv(f'{OUT}/mouse_usage.csv')
    human_usage.to_csv(f'{OUT}/human_usage.csv')

    return jmat, matches


if __name__ == '__main__':
    import sys
    mouse_k, mouse_dt, human_k, human_dt = sys.argv[1:5]
    run(int(mouse_k), mouse_dt, int(human_k), human_dt)

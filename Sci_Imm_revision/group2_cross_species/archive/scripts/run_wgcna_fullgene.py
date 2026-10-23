import sys
import pandas as pd
import numpy as np
import anndata as ad
import PyWGCNA

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'
species = sys.argv[1]

expr = pd.read_csv(f'{OUT}/metacells_{species}_fullgene_expr.csv', index_col=0)

a = ad.AnnData(X=expr.values.astype(np.float32))
a.var_names = expr.columns
a.obs_names = expr.index

wgcna = PyWGCNA.WGCNA(name=f'{species}_full', species='mouse' if species == 'mouse' else 'human',
                       anndata=a, TPMcutoff=0, save=False, outputPath=f'{OUT}/wgcna_{species}_full',
                       minModuleSize=30, MEDissThres=0.15)
wgcna.preprocess()
wgcna.findModules(kwargs_function={'cutreeHybrid': {'deepSplit': 3, 'pamRespectsDendro': False}})

mod_df = wgcna.datExpr.var[['moduleColors']].copy()
mod_df.to_csv(f'{OUT}/wgcna_{species}_fullgene_gene_modules.csv')

me = wgcna.MEs
me.index = expr.index
me.to_csv(f'{OUT}/wgcna_{species}_fullgene_module_eigengenes.csv')

print(species, 'done. modules:')
print(mod_df['moduleColors'].value_counts())

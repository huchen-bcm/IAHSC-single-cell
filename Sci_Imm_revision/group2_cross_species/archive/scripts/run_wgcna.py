import sys
import pandas as pd
import PyWGCNA

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'
species = sys.argv[1]  # 'mouse' or 'human'

expr = pd.read_csv(f'{OUT}/metacells_{species}_expr.csv', index_col=0)
# PyWGCNA expects a geneExp df with genes as columns... but also needs sample ids as a column typically for its internal GeneExp wrapper.
# Use anndata interface instead, which is simplest / most robust.
import anndata as ad
import numpy as np

a = ad.AnnData(X=expr.values.astype(np.float32))
a.var_names = expr.columns
a.obs_names = expr.index

wgcna = PyWGCNA.WGCNA(name=f'{species}_mc', species='mouse' if species == 'mouse' else 'human',
                       anndata=a, TPMcutoff=0, save=False, outputPath=f'{OUT}/wgcna_{species}',
                       minModuleSize=20, MEDissThres=0.15)
wgcna.preprocess()
wgcna.findModules(kwargs_function={'cutreeHybrid': {'deepSplit': 3, 'pamRespectsDendro': False}})
wgcna.updateSampleInfo(sampleInfo=pd.DataFrame(index=expr.index))
wgcna.saveWGCNA()

# export modules
mod_df = wgcna.datExpr.var[['moduleColors']].copy()
mod_df.to_csv(f'{OUT}/wgcna_{species}_gene_modules.csv')

me = wgcna.MEs
me.index = expr.index
me.to_csv(f'{OUT}/wgcna_{species}_module_eigengenes.csv')

print(species, 'done. modules:', mod_df['moduleColors'].value_counts())

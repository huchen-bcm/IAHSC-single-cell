from cnmf import cNMF
import numpy as np

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species/nmf_metaprogram'

cnmf_obj = cNMF(output_dir=OUT, name='mouse_cnmf')
cnmf_obj.prepare(counts_fn=f'{OUT}/mouse_ortholog_raw.h5ad',
                  components=np.arange(5, 12),
                  n_iter=40, seed=14, num_highvar_genes=2000)
cnmf_obj.factorize_multi_process(total_workers=12)
cnmf_obj.combine()
cnmf_obj.k_selection_plot(close_fig=True)
print('mouse cNMF factorize+combine done')

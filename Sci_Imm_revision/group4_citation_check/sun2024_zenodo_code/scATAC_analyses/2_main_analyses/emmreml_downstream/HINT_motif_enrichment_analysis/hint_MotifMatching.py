inputs = [
"HSC_homer_DE_peaks.bed",
"CMP_homer_DE_peaks.bed",
"GMP_homer_DE_peaks.bed",
"MEP_homer_DE_peaks.bed",
"CLP_homer_DE_peaks.bed",
"PreBNK_homer_DE_peaks.bed",
"homer_background_peaks.bed"
]


motifmatching = "/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HINT_motif_enrichment_0.05/hint_MotifMatching.sbatch"

import os
for i in range(0,7):
    os.system("sbatch " + motifmatching + " " + inputs[i] )

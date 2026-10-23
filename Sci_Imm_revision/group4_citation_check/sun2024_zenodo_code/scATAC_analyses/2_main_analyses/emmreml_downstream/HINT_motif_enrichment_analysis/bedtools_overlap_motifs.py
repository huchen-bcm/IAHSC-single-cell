
input1 = [
"HSC_homer_DE_peaks.bed",
"CMP_homer_DE_peaks.bed",
"GMP_homer_DE_peaks.bed",
"MEP_homer_DE_peaks.bed",
"CLP_homer_DE_peaks.bed",
"PreBNK_homer_DE_peaks.bed"
]

input2 = [
"HSC_HINT_mpbs.bed",
"CMP_HINT_mpbs.bed",
"GMP_HINT_mpbs.bed",
"MEP_HINT_mpbs.bed",
"CLP_HINT_mpbs.bed",
"PreBNK_HINT_mpbs.bed"
]

output = [
"HSC_HINT_out.bed",
"CMP_HINT_out.bed",
"GMP_HINT_out.bed",
"MEP_HINT_out.bed",
"CLP_HINT_out.bed",
"PreBNK_HINT_out.bed"
]


bedtools = "/project2/lbarreiro/users/Sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HINT_motif_enrichment_analysis/bedtools_overlap_motifs.sbatch"

import os
for i in range(0,6):
    os.system("sbatch " + bedtools + " " + input1[i] + " " + input2[i] + " " + output[i])
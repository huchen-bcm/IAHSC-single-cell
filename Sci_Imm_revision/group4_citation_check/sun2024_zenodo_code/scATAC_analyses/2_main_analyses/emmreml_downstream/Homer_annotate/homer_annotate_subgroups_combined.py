
inputs = [
"CMP_homer_DE_peaks.bed",
"HSC_homer_DE_peaks.bed",
"MEP_homer_DE_peaks.bed",
"PreBNK_homer_DE_peaks.bed",
"GMP_homer_DE_peaks.bed",
"CLP_homer_DE_peaks.bed"
]


outputs_go = [
"CMP_homer_annotate_GO",
"HSC_homer_annotate_GO",
"MEP_homer_annotate_GO",
"PreBNK_homer_annotate_GO",
"GMP_homer_annotate_GO",
"CLP_homer_annotate_GO"
  ]

output = [
"CMP_homer_annotate_out.txt",
"HSC_homer_annotate_out.txt",
"MEP_homer_annotate_out.txt",
"PreBNK_homer_annotate_out.txt",
"GMP_homer_annotate_out.txt",
"CLP_homer_annotate_out.txt"
]


homer_annotate = "/scATAC_analyses/2_main_analyses/emmreml_downstream/HOMER_annotate/homer_annotate_subgroups_combined.sbatch"

import os
for i in range(0,6):
  os.system("sbatch " + homer_annotate + " " + inputs[i] + " " + outputs_go[i] + " " + output[i])


inputs = [
"CMP1_homer_DE_peaks.bed",
"CMP2_homer_DE_peaks.bed",
"HSC1_homer_DE_peaks.bed",
"HSC2_homer_DE_peaks.bed",
"MEP1_homer_DE_peaks.bed",
"Pre-BNK_homer_DE_peaks.bed",
"GMP1_homer_DE_peaks.bed",
"CLP_homer_DE_peaks.bed",
"MEP2_homer_DE_peaks.bed",
"MEP3_homer_DE_peaks.bed",
"CMP3_homer_DE_peaks.bed",
"MEP4_homer_DE_peaks.bed",
"CMP4_homer_DE_peaks.bed",
"CMP5_homer_DE_peaks.bed",
"GMP2_homer_DE_peaks.bed",
"GMP3_homer_DE_peaks.bed"
  ]


outputs_go = [
"CMP1_homer_annotate_GO",
"CMP2_homer_annotate_GO",
"HSC1_homer_annotate_GO",
"HSC2_homer_annotate_GO",
"MEP1_homer_annotate_GO",
"Pre-BNK_homer_annotate_GO",
"GMP1_homer_annotate_GO",
"CLP_homer_annotate_GO",
"MEP2_homer_annotate_GO",
"MEP3_homer_annotate_GO",
"CMP3_homer_annotate_GO",
"MEP4_homer_annotate_GO",
"CMP4_homer_annotate_GO",
"CMP5_homer_annotate_GO",
"GMP2_homer_annotate_GO",
"GMP3_homer_annotate_GO"
  ]

output = [
"CMP1_homer_annotate_out.txt",
"CMP2_homer_annotate_out.txt",
"HSC1_homer_annotate_out.txt",
"HSC2_homer_annotate_out.txt",
"MEP1_homer_annotate_out.txt",
"Pre-BNK_homer_annotate_out.txt",
"GMP1_homer_annotate_out.txt",
"CLP_homer_annotate_out.txt",
"MEP2_homer_annotate_out.txt",
"MEP3_homer_annotate_out.txt",
"CMP3_homer_annotate_out.txt",
"MEP4_homer_annotate_out.txt",
"CMP4_homer_annotate_out.txt",
"CMP5_homer_annotate_out.txt",
"GMP2_homer_annotate_out.txt",
"GMP3_homer_annotate_out.txt"
  ]


homer_annotate = "/scATAC_analyses/2_main_analyses/emmreml_downstream/HOMER_annotate/homer_annotate.sbatch"

import os
for i in range(0,16):
  os.system("sbatch " + homer_annotate + " " + inputs[i] + " " + outputs_go[i] + " " + output[i])

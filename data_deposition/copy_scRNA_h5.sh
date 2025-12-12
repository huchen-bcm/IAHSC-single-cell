for i in `grep WT /home/huc/projects/Ruoqiong_Cao.IAHSC/meta/samples.scRNA.meta.tsv | cut -f1`;
do
	sample=`echo $i | cut -d"_" -f2-4 | awk 'BEGIN{OFS=""}{print "scRNA_"$1,$2,$3}'`
	echo $sample
	cp ~/projects/Ruoqiong_Cao.IAHSC/scRNA/1.mapping/cellranger_output/${i}/outs/filtered_feature_bc_matrix.h5 ./GEO_upload/scRNA/${sample}.filtered_feature_bc_matrix.h5
done


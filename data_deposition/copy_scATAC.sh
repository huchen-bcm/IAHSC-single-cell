for i in `grep WT /home/huc/projects/Ruoqiong_Cao.IAHSC/meta/samples.scATAC.csv | cut -f1`;
do
	sample=`echo $i | awk 'BEGIN{OFS=""}{print "scATAC_"$1}'`
	echo $sample

	# copy fragments
	cp /home/huc/projects/Ruoqiong_Cao.IAHSC/scATAC/1.mapping/cellranger_output/${i}/outs/fragments.tsv.gz ./GEO_upload/scATAC/${sample}.fragments.tsv.gz

	# copy peaks
	cp /home/huc/projects/Ruoqiong_Cao.IAHSC/scATAC/1.mapping/cellranger_output/${i}/outs/peaks.bed ./GEO_upload/scATAC/${sample}.peaks.bed
done


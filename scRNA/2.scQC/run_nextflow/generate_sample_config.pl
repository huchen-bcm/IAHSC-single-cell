my $input = "/home/hu/projects/Ruoqiong_Cao.IAHSC/meta/samples.scRNA.csv";
open(INPUT, "<$input");
open(OUTPUT, ">samples.tsv");
print OUTPUT "sample_id\tinput_dir\n";
my $header = <INPUT>;
while(<INPUT>){
	chomp;
	my ($id, $fq) = split(",", $_);
	print OUTPUT "${id}\t/home/hu/projects/Ruoqiong_Cao.IAHSC/scRNA/1.mapping/cellranger_output/${id}/outs/raw_feature_bc_matrix.h5\n";
}
close INPUT;
close OUTPUT;

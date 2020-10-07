#!/usr/bin/env perl

@seed_genes = ('CHP1','GPAT4','ACACA','FASN','GPI','PGP','LSS','ERO1A','SLC2A1');
#@seed_genes = ('CHP1');

#
# first pass: get all neighbors of seed genes at PVAL< 0.001
#

foreach $s (@seed_genes) {
	$seed{$s}=1;
}


open(IN, "cc-table-aml15cells-exBlood530cells-dPCC-Pval-17407genes-filtered-absZgt3.txt") || die "fail 1\n";
$skip = <IN>;
while(<IN>) {
	chomp;
	($g1, $g2, $aml, $other, $dPCC, $pval) = split(/\t/);
	if ( ($seed{$g1} || $seed{$g2}) && ($pval==0) ) {
		$keep{$g1}=1;
		$keep{$g2}=1;
#		print "$_\n";
	}
}
close(IN);

#
# repeat to get edges between neighbors
# - this allows determination of clustering and clust coeff of seeds 
# - loosen pval restriction to 0.005
#
#

open(IN, "cc-table-aml15cells-exBlood530cells-dPCC-Pval-17407genes-filtered-absZgt3.txt") || die "fail 2\n";
open(OUT, ">cc-table-aml15cells-exBlood530cells-dPCC-Pval-17407genes-filtered-absZgt3-seed_neighbors_pval.txt");

$header = <IN>;
print OUT $header;
while(<IN>) {
	chomp;
	($g1, $g2, $aml, $other, $dPCC, $pval) = split(/\t/);
	if ( $keep{$g1} && $keep{$g2} && $pval<=0.005 ) {
		print OUT "$_\n";
	}
}
close(IN);

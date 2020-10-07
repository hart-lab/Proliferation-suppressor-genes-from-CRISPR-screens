#!/usr/bin/env perl

open(Z, "zscore_aml15screens_summary_countsZgt3.txt") || die "fail 1\n";
$skip = <Z>;
while(<Z>) {
	chomp;
	($g, $ess, $tsg) = split(/\t/);
	$keep{$g} = 1 if ( ($ess>=2) || ($tsg>=2));
}
close(Z);

#
# now use this list to filter the cc-table-with-pvals
#

open(IN, "cc-table-aml15cells-exBlood530cells-dPCC-Pval-17407genes.txt") || die "fail 2\n";
$header = <IN>;
print $header;
while(<IN>) {
	chomp;
	@data = split(/\t/);
	if ( $keep{$data[0]} && $keep{ $data[1] } ) {
		print "$_\n";
	}
}
close(IN);
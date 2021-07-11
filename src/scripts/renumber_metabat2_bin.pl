#!/usr/bin/perl

if (@ARGV != 1) {
	print STDERR ("usage: $0 metabat2_out.tsv(--saveCls)\n");
	exit(0);
}

open($in, $ARGV[0]);
while (chomp($l = <$in>)) {
	($name, $id) = split(/\t/, $l);
	$id_flag{$id} = 1;
}
close $in;

$count = 0;
for $id (sort{$1 <=> $b} keys %id_flag) {
	if ($id == 0) {
		$new_id{$id} = 0;
	}
	else {
		++$count;
		$new_id{$id} = $count;
	}
}

open($in, $ARGV[0]);
while (chomp($l = <$in>)) {
	($name, $id) = split(/\t/, $l);
	print(join("\t", ($name, $new_id{$id})), "\n");
}
close $in;

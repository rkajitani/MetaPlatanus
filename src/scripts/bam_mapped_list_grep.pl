#!/usr/bin/perl

(@ARGV != 2) and die "usage: $0 seq_name_list in.bam\n";
$list_name = $ARGV[0];
$bam_name = $ARGV[1];

open IN, $list_name;
while (<IN>) {
	chomp;
	$target{$_} = 1;
}
close IN;

open IN, "samtools view -h $bam_name |";
open OUT, "| samtools view -Sb -";

while ($l1 = <IN>) {
	if ($l1 =~ /^@/) {
		if (not($l1 =~ /^\@SQ/)) {
			print OUT $l1;
		}
		elsif ($l1 =~ /SN:(\S+)\s/ and $target{$1}) {
			print OUT $l1;
		}
		next;
	}

	$l2 = <IN>;

	if ($target{(split(/\t/, $l1))[2]} and $target{(split(/\t/, $l2))[2]}) {
		print OUT $l1;
		print OUT $l2;
	}
}
close IN;
close OUT;

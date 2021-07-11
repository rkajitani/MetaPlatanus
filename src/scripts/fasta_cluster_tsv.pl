#!/usr/bin/perl

if (@ARGV == 0) {
	print STDERR ("usage: $0 in1.fa [in2.fa ...]\n");
	exit(0);
}

for $file (@ARGV) {
	open($in, $file);
	while ($l = <$in>) {
		if ($l =~ /^>(\S+)/) {
			$name = $1;
			if ($name =~ /_cluster(\d+)/ and $1 ne "0") {
				print("$name\t$1\n");
			}
		}
	}
	close $in;
}

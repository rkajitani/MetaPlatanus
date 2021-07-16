#!/usr/bin/env perl

if (@ARGV == 0) {
	print STDERR ("usage: $0 in1.fa [in2.fa ...]\n");
	exit(0);
}

for $file (@ARGV) {
	open($in, $file);
	while (($name, $seq) = fasta_get($in)) {
		$name =~ s/_read\d+_maxK\d+//;
		$name =~ s/_np1212$//;
		if ($name =~ /(^\S+_len)\d+(\S+)/) {
			$name = $1 . length($seq) . $2;
		}
		$seq = uc($seq);
		print(">$name\n");
		print_seq($seq, 80);
	}
	close $in;
}


sub fasta_get
{
	my $in = shift;
	my($l, $name, $seq);

	while ($l = <$in>) {
		if ($l =~ /^>(\S*)/) {
			$name = $1;
			last;
		}
	}
	if ($l eq '') {
		return ();
	}

	while ($l = <$in>) {
		if (substr($l, 0, 1) eq '>') {
			seek($in, -length($l), 1);
			return ($name, $seq);
		}
		chomp $l;
		$seq .= $l 
	}
	return ($name, $seq);
}

sub print_seq
{
	my $seq = shift;
	my $line_len = shift;
	my $seq_len = length $seq;

	unless ($line_len) {
		$line_len = 80;
	}

	for (my $i = 0; $i < $seq_len; $i += $line_len) {
		print(substr($seq, $i, $line_len) . "\n");
	}
}


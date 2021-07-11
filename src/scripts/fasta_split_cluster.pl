#!/usr/bin/perl

if (@ARGV != 3) {
	print STDERR ("usage: $0 scf.fa clusters.tsv out_dir\n");
	exit(0);
}

open($in, $ARGV[1]);
while (chomp($l = <$in>)) {
	@f = split(/\t/, $l);
	$name2cl{$f[0]} = $f[1];
}
close $in;

open($in, $ARGV[0]);
while (($name, $seq) = fasta_get($in)) {
	push(@{$names[$name2cl{$name}]}, $name);
	push(@{$seqs[$name2cl{$name}]}, $seq);
}
close $in;

$out_dir = $ARGV[2];
mkdir($out_dir) or die("error: cannot make directory $out_dir\n");
for $i (0..$#names) {
	if ($i == 0) {
		open($out, ">$out_dir/unclassified.fa");
	}
	else {
		open($out, ">$out_dir/cluster$i.fa");
	}
	
	for $j (0..(@{$names[$i]} - 1)) {
		printf($out ">%s_cluster%d\n", $names[$i][$j], $i);
		print_scf($out, $seqs[$i][$j], 80);
	}

	close $out;
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

sub print_scf
{
	my $out = shift;
	my $seq = shift;
	my $line_len = shift;
	my $seq_len = length $seq;

	unless ($line_len) {
		$line_len = 80;
	}

	for (my $i = 0; $i < $seq_len; $i += $line_len) {
		print($out substr($seq, $i, $line_len) . "\n");
	}
}

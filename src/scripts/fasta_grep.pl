#! /usr/bin/perl

use Getopt::Long;
use FindBin;
use lib "$FindBin::Bin";
use my_bio;

%opt = ();
GetOptions(\%opt, qw/e f v x m=i help/);

if (@ARGV < 2 or exists $opt{help}) {
	print_usage();
	exit;
}

$ptn = shift @ARGV;;
if (@ARGV == 0) {
	push(@files, *STDIN);
}
else {
	@files = @ARGV;
}

$template =  <<'EOS';

%s
for $file (@files) {
	open($fh, $file) or warn "WARNING: cannot open $file\n";
	while (($name, $seq) = fasta_nonwhite_get($fh)) {
		if (%s) {
			print ">$name\n";
			print_seq($seq, 80);
			%s
		}
	}
	close $fh;
}

EOS

if (exists $opt{m}) {
	$first_exp = "\$max_count = $opt{m};";
	$last_exp = '++$count; last if ($count >= $max_count);';
}

if (exists $opt{v}) {
	$match_cond = 'not ';
}
if (exists $opt{e}) {
	$match_cond .= "\$name =~ /$ptn/";
}
elsif (exists $opt{f}) {
	open IN, $ptn;
	while (<IN>) {
		chomp;
		$is_match{$_} = 1;
	}
	close IN;
	$match_cond .= "\$is_match\{\$name\}";
}
elsif (exists $opt{x}) {
	$match_cond .= "\$name eq \'$ptn\'";
}
else {
	$match_cond .= '$name =~ /$ptn/';
}

eval sprintf($template, $first_exp, $match_cond, $last_exp);



sub print_usage
{
	
	print
<<'USAGE';
Usage: faste_grep.pl [OPTION]... PATTERN FASTA ...
Search for PATTERN in each FASTA.
Example: grep gi000000 in.fasta

Regexp selection and interpretation:
  -e         PATTERN is a Perl regular expression
  -f         OR exact match to PATTERNS in the file
  -x         force PATTERN to match only whole title

Miscellaneous:
  -v         select non-matching lines

Output control:
  -m         stop after NUM matches
  --help     display this help and exit
USAGE
}

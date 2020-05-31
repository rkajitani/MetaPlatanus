#!/usr/bin/perl

use strict;
use File::Spec;
use FindBin;
use Getopt::Long qw/:config posix_default no_ignore_case/;

my $version = 'v1.2.2';

printf("meta_platanus.pl version %s\n\n", $version, "\n");

my $cmd = shift(@ARGV);

if ($cmd eq 'cons_asm') {
	cons_asm_main(@ARGV);
}
elsif ($cmd eq 'phase_asm') {
	phase_asm_main(@ARGV);
}
else {
	print_usage();
}



sub cons_asm_main
{
	my @argv = @_;

	my $help_flag = 0;
	my $overwrite_flag = 0;

	my $sub_bin = ($FindBin::Bin . '/sub_bin');
	my $mem = get_mem_total() * 0.75;
	my $out = 'out';
	my $n_thread = 1;
	my $tmp_dir = '.';

	my $megahit_min_len = 500;
	my @lib_file;
	my @lib_direction;
	my @binning_lib_file;
	my @pac_file = ();
	my @ont_file = ();
	my @tenx_file_interleaved = ();
	my @tenx_file_paired = ();

	my %opt_fmt = (
		'h' => \$help_flag,
		'help' => \$help_flag,
		'overwrite' => \$overwrite_flag,

		'm=i' => \$mem,
		'o=s' => \$out,
		't=i' => \$n_thread,
		'tmp=s' => \$tmp_dir,
		'sub_bin=s' => \$sub_bin,

		'megahit_min_len=i' => \$megahit_min_len,
		'p=s{1,}' => \@pac_file,
		'ont=s{1,}' => \@ont_file,
		'x=s{1,}' => \@tenx_file_interleaved,
		'X=s{1,}' => \@tenx_file_paired,
	);

	for my $arg (@argv) {
		if ($arg =~ /-([IO])P(\d+)/) {
			my $lib_i = $2 - 1;
			$lib_file[$lib_i] = [];
			if (defined($lib_direction[$lib_i]) and $lib_direction[$lib_i] ne $1) {
				die(sprintf("Error: Confliction of directions of short-reads libraries (%s).\n", $arg));
			}
			$lib_direction[$lib_i] = $1;
			$opt_fmt{sprintf('%sP%s=s{1,}', $lib_direction[$lib_i], $lib_i + 1)} = $lib_file[$lib_i];
		}
		elsif ($arg =~ /-binning_IP(\d+)/) {
			my $lib_i = $1 - 1;
			$binning_lib_file[$lib_i] = [];
			$opt_fmt{sprintf('binning_IP%s=s{1,}', $lib_i + 1)} = $binning_lib_file[$lib_i];
		}
	}

	GetOptions(%opt_fmt);

	if (@argv == 0 or  $help_flag) {
		print_cons_asm_usage();
		return;
	}

	$sub_bin = File::Spec->rel2abs($sub_bin);
	$ENV{"PATH"} = ($sub_bin . ':' . $ENV{"PATH"});

	for (my $i = 0; $i < @lib_file; ++$i) {
		rel2abs_multi($lib_file[$i]);
	}
	for (my $i = 0; $i < @binning_lib_file; ++$i) {
		rel2abs_multi($binning_lib_file[$i]);
	}
	rel2abs_multi(\@pac_file);
	rel2abs_multi(\@ont_file);
	rel2abs_multi(\@tenx_file_interleaved);
	rel2abs_multi(\@tenx_file_paired);

	my $contig_asm_dir = sprintf("%s_contig_asm", $out);
	if ($overwrite_flag or not(-e "$contig_asm_dir/done")) {
		contig_asm(
			$contig_asm_dir,
			$sub_bin,
			$n_thread,
			$mem,
			$tmp_dir,
			$out,
			\@lib_file,
			\@lib_direction,
		);
	}

	my $contig_asm_megahit_dir = sprintf("%s_contig_asm_megahit", $out);
	if ($overwrite_flag or not(-e "$contig_asm_megahit_dir/done")) {
		contig_asm_megahit(
			$contig_asm_megahit_dir,
			$sub_bin,
			$n_thread,
			$tmp_dir,
			\@lib_file,
			\@lib_direction,
		);
	}

	my $iterative_asm_dir = sprintf("%s_iterative_asm", $out);
	if ($overwrite_flag or not(-e "$iterative_asm_dir/done")) {
		iterative_asm(
			$iterative_asm_dir,
			$contig_asm_dir,
			$contig_asm_megahit_dir,
			$sub_bin,
			$n_thread,
			$mem,
			$tmp_dir,
			$out,
			$megahit_min_len,
			\@lib_file,
			\@lib_direction,
			\@pac_file,
			\@ont_file,
			\@tenx_file_paired,
			\@tenx_file_interleaved,
		);
	}

	my $binning_metabat_dir = sprintf("%s_binning_metabat", $out);
	if ($overwrite_flag or not(-e "$binning_metabat_dir/done")) {
		binning_metabat(
			$binning_metabat_dir,
			$iterative_asm_dir,
			$sub_bin,
			$n_thread,
			$mem,
			$tmp_dir,
			$out,
			\@lib_file,
			\@lib_direction,
			\@binning_lib_file,
		);
	}

	my $re_scaffold_dir = sprintf("%s_re_scaffold", $out);
	if ($overwrite_flag or not(-e "$re_scaffold_dir/done")) {
		re_scaffold(
			$re_scaffold_dir,
			$binning_metabat_dir,
			$sub_bin,
			$n_thread,
			$tmp_dir,
			$out,
			\@lib_file,
			\@lib_direction,
			\@pac_file,
			\@ont_file,
			\@tenx_file_paired,
			\@tenx_file_interleaved,
		);
	}
	
	system("mv $re_scaffold_dir/${out}_finalClusters_all.fa .");
	system("mv $re_scaffold_dir/${out}_finalClusters .");

	printf("Final results: %s, %s\n", "${out}_finalClusters_all.fa", "${out}_finalClusters");
	print("meta_platanus.pl cons_asm finished.\n");
	exit(0);
}


sub contig_asm
{
	my $work_dir = shift;
	my $sub_bin = shift;
	my $n_thread = shift;
	my $mem = shift;
	my $tmp_dir = shift;
	my $out = shift;
	my $lib_file = shift;
	my $lib_direction = shift;

	mkdir($work_dir);
	chdir($work_dir) or die(sprintf("Error: can not make or enter %s (directory).\n$!", $work_dir));

	my $file_input_str = '';
	for (my $i = 0; $i < @{$lib_file}; ++$i) {
		if (not(defined($lib_file->[$i]))) {
			next;
		}
		for (my $j = 0; $j < @{$lib_file->[$i]}; ++$j) {
			if ($lib_direction->[$i] eq 'I') {
				$file_input_str .= ($lib_file->[$i][$j] . " ");
			}
		}
	}

	my $script_fh;
	open($script_fh, ">cmd.bash");

	print $script_fh <<_EOS;
#!/usr/bin/bash

t=$n_thread
m=$mem
tmp=$tmp_dir
o=$out
sub_bin=$sub_bin

/usr/bin/time \$sub_bin/meta_platanus assemble -tmp \$tmp -t \$t -m \$m -f $file_input_str -o \${o} 2>\${o}.assembleLog
_EOS

	system("bash cmd.bash >cmd.stdout 2>cmd.stderr") == 0 or die "Error: contig_assembly error.\n$!\n";

	system('touch done');
	chdir('..');
}


sub contig_asm_megahit
{
	my $work_dir = shift;
	my $sub_bin = shift;
	my $n_thread = shift;
	my $tmp_dir = shift;
	my $lib_file = shift;
	my $lib_direction = shift;

	mkdir($work_dir);
	chdir($work_dir) or die(sprintf("Error: can not make or enter %s (directory).\n$!", $work_dir));

	my @file_input = ();
	for (my $i = 0; $i < @{$lib_file}; ++$i) {
		if (not(defined($lib_file->[$i]))) {
			next;
		}
		for (my $j = 0; $j < @{$lib_file->[$i]}; ++$j) {
			if ($lib_direction->[$i] eq 'I') {
				push(@file_input, $lib_file->[$i][$j]);
			}
		}
	}

	my @file_input1 = ();
	my @file_input2 = ();
	for (my $i = 0; $i + 1 < @file_input; $i += 2) {
		push(@file_input1, $file_input[$i]);
		push(@file_input2, $file_input[$i + 1]);
	}
	my $file_input_str1 = join(',', @file_input1);
	my $file_input_str2 = join(',', @file_input2);

	my $script_fh;
	open($script_fh, ">cmd.bash");

	print $script_fh <<_EOS;
#!/usr/bin/bash

t=$n_thread
tmp=$tmp_dir
sub_bin=$sub_bin
PE1=$file_input_str1
PE2=$file_input_str2

/usr/bin/time \$sub_bin/megahit -t \$t --tmp-dir \$tmp -1 \$PE1 -2 \$PE2 >megahit.stdout 2>megahit.stderr
rm -r megahit_out/intermediate_contigs
_EOS

	system("bash cmd.bash >cmd.stdout 2>cmd.stderr") == 0 or die "Error: megahit_contig_assembly error.\n$!\n";

	system('touch done');
	chdir('..');
}


sub iterative_asm
{
	my $work_dir = shift;
	my $metaplatanus_contig_dir = shift;
	my $megahit_contig_dir = shift;
	my $sub_bin = shift;
	my $n_thread = shift;
	my $mem = shift;
	my $tmp_dir = shift;
	my $out = shift;
	my $megahit_min_len = shift;
	my $lib_file = shift;
	my $lib_direction = shift;
	my $pac_file = shift;
	my $ont_file = shift;
	my $tenx_file_paired = shift;
	my $tenx_file_interleaved = shift;

	$metaplatanus_contig_dir = File::Spec->rel2abs($metaplatanus_contig_dir);
	$megahit_contig_dir = File::Spec->rel2abs($megahit_contig_dir);

	mkdir($work_dir);
	chdir($work_dir) or die(sprintf("Error: can not make or enter %s (directory).\n$!", $work_dir));

	my $file_input_str = '';
	for (my $i = 0; $i < @{$lib_file}; ++$i) {
		if (not(defined($lib_file->[$i]))) {
			next;
		}
		$file_input_str .= sprintf("-%sP%d ", $lib_direction->[$i], $i + 1);
		for (my $j = 0; $j < @{$lib_file->[$i]}; ++$j) {
			$file_input_str .= ($lib_file->[$i][$j] . " ");
		}
	}

	if (ref($pac_file) eq 'ARRAY' and @{$pac_file} > 0) {
		$file_input_str .= " -p ";
		for (my $i = 0; $i < @{$pac_file}; ++$i) {
			$file_input_str .= ($pac_file->[$i] . " ");
		}
	}

	if (ref($ont_file) eq 'ARRAY' and @{$ont_file} > 0) {
		$file_input_str .= " -ont ";
		for (my $i = 0; $i < @{$ont_file}; ++$i) {
			$file_input_str .= ($ont_file->[$i] . " ");
		}
	}

	if (ref($tenx_file_interleaved) eq 'ARRAY' and @{$tenx_file_interleaved} > 0) {
		$file_input_str .= " -x ";
		for (my $i = 0; $i < @{$tenx_file_interleaved}; ++$i) {
			$file_input_str .= ($tenx_file_interleaved->[$i] . " ");
		}
	}

	if (ref($tenx_file_paired) eq 'ARRAY' and @{$tenx_file_paired} > 0) {
		$file_input_str .= " -X ";
		for (my $i = 0; $i < @{$tenx_file_paired}; ++$i) {
			$file_input_str .= ($tenx_file_paired->[$i] . " ");
		}
	}

	my $script_fh;
	open($script_fh, ">cmd.bash");

	print $script_fh <<_EOS;
#!/usr/bin/bash

t=$n_thread
m=$mem
tmp=$tmp_dir
o=$out
megahit_min_len=$megahit_min_len
metaplatanus_contig_dir=$metaplatanus_contig_dir
megahit_contig_dir=$megahit_contig_dir
sub_bin=$sub_bin

ln -fs \$metaplatanus_contig_dir/\${o}_contig.fa >/dev/null 2>&1
ln -fs \$metaplatanus_contig_dir/\${o}_junctionKmer.fa >/dev/null 2>&1
ln -fs \$metaplatanus_contig_dir/\${o}_kmer_occ.bin >/dev/null 2>&1

$sub_bin/fasta_length_filter.pl \$megahit_contig_dir/megahit_out/final.contigs.fa \$megahit_min_len >megahit.fa

/usr/bin/time \$sub_bin/meta_platanus divide -recalc_cov -tmp \$tmp -k \${o}_kmer_occ.bin -c \${o}_contig.fa \${o}_junctionKmer.fa megahit.fa -o initial 2>initial_recalc.log
/usr/bin/time \$sub_bin/meta_platanus merge -fixed_k_len -tmp \$tmp -m \$m -f initial_recalc.fa -o initial 2>initial_merge.log
/usr/bin/time \$sub_bin/meta_platanus iterate -tmp \$tmp -c initial_merged.fa initial_mergedJunctionKmer.fa -k \${o}_kmer_occ.bin -t \$t -m \$m $file_input_str -o \${o} 2>\${o}.iteLog
_EOS

	system("bash cmd.bash >cmd.stdout 2>cmd.stderr") == 0 or die "Error: iterative_assembly error.\n$!\n";

	system('touch done');
	chdir('..');
}


sub binning_metabat 
{
	my $work_dir = shift;
	my $iterative_asm_dir = shift;
	my $sub_bin = shift;
	my $n_thread = shift;
	my $mem = shift;
	my $tmp_dir = shift;
	my $out = shift;
	my $lib_file = shift;
	my $lib_direction = shift;
	my $binning_lib_file = shift;

	$iterative_asm_dir = File::Spec->rel2abs($iterative_asm_dir);

	mkdir($work_dir);
	chdir($work_dir) or die(sprintf("Error: can not make or enter %s (directory).\n$!", $work_dir));

	my @file_input = ();
	for (my $i = 0; $i < @{$binning_lib_file}; ++$i) {
		if (not(defined($binning_lib_file->[$i]))) {
			next;
		}
		for (my $j = 0; $j < @{$binning_lib_file->[$i]}; ++$j) {
			push(@{$file_input[$i]}, $binning_lib_file->[$i][$j]);
		}
	}

	if (@file_input == 0) {
		for (my $i = 0; $i < @{$lib_file}; ++$i) {
			if (not(defined($lib_file->[$i]))) {
				next;
			}
			for (my $j = 0; $j < @{$lib_file->[$i]}; ++$j) {
				if ($lib_direction->[$i] eq 'I') {
					push(@{$file_input[0]}, $lib_file->[$i][$j]);
				}
			}
		}
	}
		
	my $mem_per_thread = (int($mem / $n_thread) . 'g');
	my $map_cmd_str = '';
	my $bam_str = '';
	for (my $i = 0; $i < @file_input; ++$i) {
		if (not(defined($file_input[$i]))) {
			next;
		}
		my @file_input1 = ();
		my @file_input2 = ();
		for (my $j = 0; $j + 1 < @{$file_input[$i]}; $j += 2) {
			push(@file_input1, $file_input[$i][$j]);
			push(@file_input2, $file_input[$i][$j + 1]);
		}

		my $lib_id = ("lib_" . ($i + 1));
		$map_cmd_str .= sprintf(
			"\$sub_bin/bwa mem -t \$t base.fa <(cat %s) <(cat %s) 2>>%s_bwa_mem.log | \$sub_bin/samtools view -b - >%s_raw.bam 2>%s_bwa_mem.log\n",
			join(' ', @file_input1),
			join(' ', @file_input2),
			$lib_id,
			$lib_id,
			$lib_id,
		);
			
		$map_cmd_str .= sprintf(
			"\$sub_bin/samtools sort -m %s -\@ \$t -o %s_sorted.bam %s_raw.bam >%s_samtools_sort.log 2>&1\n",
			$mem_per_thread,
			$lib_id,
			$lib_id,
			$lib_id,
		);

		$bam_str .= sprintf("%s_sorted.bam ", $lib_id);
	}
	
	my $script_fh;
	open($script_fh, ">cmd.bash");

	print $script_fh <<_EOS;
#!/usr/bin/bash

t=$n_thread
tmp=$tmp_dir
o=$out
iterative_asm_dir=$iterative_asm_dir
sub_bin=$sub_bin

ln -fs \$iterative_asm_dir/\${o}_iterativeAssembly.fa base.fa
\$sub_bin/bwa index base.fa >bwa_index.log 2>&1

$map_cmd_str

BAMS="$bam_str"
/usr/bin/time \$sub_bin/jgi_summarize_bam_contig_depths --outputDepth depth.txt \$BAMS >jgi_summarize_bam_contig_depths.log 2>&1
/usr/bin/time \$sub_bin/metabat2 -t \$t -i base.fa -a depth.txt -o metabat2_out.tsv -l -v --saveCls --noBinOut >metabat2.stdout 2>metabat2.stderr
\$sub_bin/renumber_metabat2_bin.pl metabat2_out.tsv >metabat2_out_renumbered.tsv

rm *.bam
_EOS

	system("bash cmd.bash >cmd.stdout 2>cmd.stderr") == 0 or die "Error: binning_metabat2 error.\n$!\n";

	system('touch done');
	chdir('..');
}


sub re_scaffold
{
	my $work_dir = shift;
	my $binning_metabat_dir = shift;
	my $sub_bin = shift;
	my $n_thread = shift;
	my $tmp_dir = shift;
	my $out = shift;
	my $lib_file = shift;
	my $lib_direction = shift;
	my $pac_file = shift;
	my $ont_file = shift;
	my $tenx_file_paired = shift;
	my $tenx_file_interleaved = shift;

	$binning_metabat_dir = File::Spec->rel2abs($binning_metabat_dir);

	mkdir($work_dir);
	chdir($work_dir) or die(sprintf("Error: can not make or enter %s (directory).\n$!", $work_dir));

	my $scaffold_input_str = '';
	for (my $i = 0; $i < @{$lib_file}; ++$i) {
		if (not(defined($lib_file->[$i]))) {
			next;
		}
		$scaffold_input_str .= sprintf("-%sP%d ", $lib_direction->[$i], $i + 1);
		for (my $j = 0; $j < @{$lib_file->[$i]}; ++$j) {
			$scaffold_input_str .= ($lib_file->[$i][$j] . " ");
		}
	}
	my $fill_input_str = $scaffold_input_str;

	if (ref($pac_file) eq 'ARRAY' and @{$pac_file} > 0) {
		$scaffold_input_str .= " -p ";
		for (my $i = 0; $i < @{$pac_file}; ++$i) {
			$scaffold_input_str .= ($pac_file->[$i] . " ");
		}
	}

	if (ref($ont_file) eq 'ARRAY' and @{$ont_file} > 0) {
		$scaffold_input_str .= " -ont ";
		for (my $i = 0; $i < @{$ont_file}; ++$i) {
			$scaffold_input_str .= ($ont_file->[$i] . " ");
		}
	}

	if (ref($tenx_file_interleaved) eq 'ARRAY' and @{$tenx_file_interleaved} > 0) {
		$scaffold_input_str .= " -x ";
		for (my $i = 0; $i < @{$tenx_file_interleaved}; ++$i) {
			$scaffold_input_str .= ($tenx_file_interleaved->[$i] . " ");
		}
	}

	if (ref($tenx_file_paired) eq 'ARRAY' and @{$tenx_file_paired} > 0) {
		$scaffold_input_str .= " -X ";
		for (my $i = 0; $i < @{$tenx_file_paired}; ++$i) {
			$scaffold_input_str .= ($tenx_file_paired->[$i] . " ");
		}
	}

	my $script_fh;
	open($script_fh, ">cmd.bash");

	print $script_fh <<_EOS;
#!/usr/bin/bash

t=$n_thread
tmp=$tmp_dir
o=$out
binning_metabat_dir=$binning_metabat_dir
sub_bin=$sub_bin

ln -s \$binning_metabat_dir/metabat2_out_renumbered.tsv >/dev/null 2>&1
ln -s \$binning_metabat_dir/base.fa >/dev/null 2>&1

/usr/bin/time \$sub_bin/meta_platanus solve_DBG -unphase -tmp \$tmp -g metabat2_out_renumbered.tsv -c base.fa $scaffold_input_str -t \$t -o \$o 2>\${o}.rescaffoldLog
cat \${o}_extendedClusters/*.fa >\${o}_extendedClusters_all.fa
/usr/bin/time \$sub_bin/meta_platanus cluster_fill -tmp \$tmp -c \${o}_extendedClusters_all.fa $fill_input_str -t \$t -o \$o 2>\${o}.fillLog
mv \${o}_filledClusters \${o}_finalClusters 
mv \${o}_filledClusters_all.fa \${o}_finalClusters_all.fa
_EOS

	system("bash cmd.bash >cmd.stdout 2>cmd.stderr") == 0 or die "Error: re_scaffold error.\n$!\n";

	system('touch done');
	chdir('..');
}


sub print_usage
{
	print

<<'USAGE';
Usage:
	meta_platanus.pl [Command] ...

Commands:
    cons_asm   : construct consensus sequences using contig-assembly, scaffolding and binning processes.
    phase_asm  : construct phased haplotype sequences using contig-assembly and scaffolding processes.

USAGE

}


sub print_cons_asm_usage
{
	print

<<'USAGE';
Usage:
    meta_platanus.pl cons_asm -IP1 short_R1.fastq(a) short_R2.fastq(a) [Options] ...

Options:
    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq; at least one library required)
    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq; aka mate-pairs or jumping-library)
    -binning_IP{INT} FWD1 REV1 ...     : lib_id inward_pair_files for binning process. (reads in 2 files, fasta or fastq; the data are usually from another sample)
    -p FILE1 [FILE2 ...]               : PacBio long-read file (fasta or fastq)
    -ont FILE1 [FILE2 ...]             : Oxford Nanopore long-read file (fasta or fastq)
    -x PAIR1 [PAIR2 ...]               : barcoded_pair_files (10x Genomics) (reads in 1 file, interleaved, fasta or fastq)
    -X FWD1 REV1 [FWD2 REV2 ...]       : barcoded_pair_files (10x Genomics) (reads in 2 files, fasta or fastq)
    -t INT                             : number of threads (<= 1; default, 1)
    -m INT                             : memory limit for making kmer distribution (unit, GB; default, 0.75 * available_memory))
    -o STR                             : prefix of output files (default "out")
    -tmp DIR                           : directory for temporary files (default, ".")
    -sub_bin DIR                       : directory for sub-executables, such as mata_plantaus and minimap2 (default, directory-of-this-script/sub_bin)
    -megahit_min_len                   : minimum length of contigs of MEGAHIT (default, 500)
    -overwrite                         : overwrite the previous results, not re-start (default, off)
    -h, -help                          : display usage

USAGE
}


sub phase_asm_main
{
	my @argv = @_;

	my $help_flag = 0;
	my $overwrite_flag = 0;

	my $sub_bin = ($FindBin::Bin . '/sub_bin');
	my $mem = get_mem_total() * 0.75;
	my $out = 'out';
	my $n_thread = 1;
	my $tmp_dir = '.';

	my @lib_file;
	my @lib_direction;
	my @pac_file = ();
	my @ont_file = ();
	my @tenx_file_interleaved = ();
	my @tenx_file_paired = ();

	my %opt_fmt = (
		'h' => \$help_flag,
		'help' => \$help_flag,
		'overwrite' => \$overwrite_flag,

		'm=i' => \$mem,
		'o=s' => \$out,
		't=i' => \$n_thread,
		'tmp=s' => \$tmp_dir,
		'sub_bin=s' => \$sub_bin,

		'p=s{1,}' => \@pac_file,
		'ont=s{1,}' => \@ont_file,
		'x=s{1,}' => \@tenx_file_interleaved,
		'X=s{1,}' => \@tenx_file_paired,
	);

	for my $arg (@argv) {
		if ($arg =~ /-([IO])P(\d+)/) {
			my $lib_i = $2 - 1;
			$lib_file[$lib_i] = [];
			if (defined($lib_direction[$lib_i]) and $lib_direction[$lib_i] ne $1) {
				die(sprintf("Error: Confliction of directions of short-reads libraries (%s).\n", $arg));
			}
			$lib_direction[$lib_i] = $1;
			$opt_fmt{sprintf('%sP%s=s{1,}', $lib_direction[$lib_i], $lib_i + 1)} = $lib_file[$lib_i];
		}
	}

	GetOptions(%opt_fmt);

	if (@argv == 0 or $help_flag) {
		print_phase_asm_usage();
		return;
	}

	$sub_bin = File::Spec->rel2abs($sub_bin);
	$ENV{"PATH"} = ($sub_bin . ':' . $ENV{"PATH"});

	for (my $i = 0; $i < @lib_file; ++$i) {
		rel2abs_multi($lib_file[$i]);
	}
	rel2abs_multi(\@pac_file);
	rel2abs_multi(\@ont_file);
	rel2abs_multi(\@tenx_file_interleaved);
	rel2abs_multi(\@tenx_file_paired);

	my $contig_asm_dir = sprintf("%s_contig_phase_asm", $out);
	if ($overwrite_flag or not(-e "$contig_asm_dir/done")) {
		contig_phase_asm(
			$contig_asm_dir,
			$sub_bin,
			$n_thread,
			$mem,
			$tmp_dir,
			$out,
			\@lib_file,
			\@lib_direction,
		);
	}

	my $phase_asm_dir = sprintf("%s_phase_asm", $out);
	if ($overwrite_flag or not(-e "$phase_asm_dir/done")) {
		phase_asm(
			$phase_asm_dir,
			$contig_asm_dir,
			$sub_bin,
			$n_thread,
			$mem,
			$tmp_dir,
			$out,
			\@lib_file,
			\@lib_direction,
			\@pac_file,
			\@ont_file,
			\@tenx_file_paired,
			\@tenx_file_interleaved,
		);
	}
	
	system("mv $phase_asm_dir/${out}_primaryBubble.fa .");
	system("mv $phase_asm_dir/${out}_secondaryBubble.fa .");
	system("mv $phase_asm_dir/${out}_nonBubbleOther.fa .");

	printf("Main results: %s, %s\n", "${out}_primaryBubble.fa", "${out}_secondaryBubble.fa", "${out}_nonBubbleOther.fa");
	print("meta_platanus.pl phase_asm finished.\n");
	exit(0);
}


sub contig_phase_asm
{
	my $work_dir = shift;
	my $sub_bin = shift;
	my $n_thread = shift;
	my $mem = shift;
	my $tmp_dir = shift;
	my $out = shift;
	my $lib_file = shift;
	my $lib_direction = shift;

	mkdir($work_dir);
	chdir($work_dir) or die(sprintf("Error: can not make or enter %s (directory).\n$!", $work_dir));

	my $file_input_str = '';
	for (my $i = 0; $i < @{$lib_file}; ++$i) {
		if (not(defined($lib_file->[$i]))) {
			next;
		}
		for (my $j = 0; $j < @{$lib_file->[$i]}; ++$j) {
			if ($lib_direction->[$i] eq 'I') {
				$file_input_str .= ($lib_file->[$i][$j] . " ");
			}
		}
	}

	my $script_fh;
	open($script_fh, ">cmd.bash");

	print $script_fh <<_EOS;
#!/usr/bin/bash

t=$n_thread
m=$mem
tmp=$tmp_dir
o=$out
sub_bin=$sub_bin

/usr/bin/time \$sub_bin/meta_platanus_phase assemble -u 0 -l 1.0 -tmp \$tmp -t \$t -m \$m -f $file_input_str -o \${o} 2>\${o}.assembleLog
_EOS

	system("bash cmd.bash >cmd.stdout 2>cmd.stderr") == 0 or die "Error: contig_assembly error.\n$!\n";

	system('touch done');
	chdir('..');
}


sub phase_asm
{
	my $work_dir = shift;
	my $contig_dir = shift;
	my $sub_bin = shift;
	my $n_thread = shift;
	my $mem = shift;
	my $tmp_dir = shift;
	my $out = shift;
	my $lib_file = shift;
	my $lib_direction = shift;
	my $pac_file = shift;
	my $ont_file = shift;
	my $tenx_file_paired = shift;
	my $tenx_file_interleaved = shift;

	$contig_dir = File::Spec->rel2abs($contig_dir);

	mkdir($work_dir);
	chdir($work_dir) or die(sprintf("Error: can not make or enter %s (directory).\n$!", $work_dir));

	my $file_input_str = '';
	for (my $i = 0; $i < @{$lib_file}; ++$i) {
		if (not(defined($lib_file->[$i]))) {
			next;
		}
		$file_input_str .= sprintf("-%sP%d ", $lib_direction->[$i], $i + 1);
		for (my $j = 0; $j < @{$lib_file->[$i]}; ++$j) {
			$file_input_str .= ($lib_file->[$i][$j] . " ");
		}
	}

	if (ref($pac_file) eq 'ARRAY' and @{$pac_file} > 0) {
		$file_input_str .= " -p ";
		for (my $i = 0; $i < @{$pac_file}; ++$i) {
			$file_input_str .= ($pac_file->[$i] . " ");
		}
	}

	if (ref($ont_file) eq 'ARRAY' and @{$ont_file} > 0) {
		$file_input_str .= " -p ";
		for (my $i = 0; $i < @{$ont_file}; ++$i) {
			$file_input_str .= ($ont_file->[$i] . " ");
		}
	}

	if (ref($tenx_file_interleaved) eq 'ARRAY' and @{$tenx_file_interleaved} > 0) {
		$file_input_str .= " -x ";
		for (my $i = 0; $i < @{$tenx_file_interleaved}; ++$i) {
			$file_input_str .= ($tenx_file_interleaved->[$i] . " ");
		}
	}

	if (ref($tenx_file_paired) eq 'ARRAY' and @{$tenx_file_paired} > 0) {
		$file_input_str .= " -X ";
		for (my $i = 0; $i < @{$tenx_file_paired}; ++$i) {
			$file_input_str .= ($tenx_file_paired->[$i] . " ");
		}
	}

	my $script_fh;
	open($script_fh, ">cmd.bash");

	print $script_fh <<_EOS;
#!/usr/bin/bash

t=$n_thread
m=$mem
tmp=$tmp_dir
o=$out
sub_bin=$sub_bin
contig_dir=$contig_dir

ln -fs \$contig_dir/\${o}_contig.fa >/dev/null 2>&1
ln -fs \$contig_dir/\${o}_junctionKmer.fa >/dev/null 2>&1
ln -fs \$contig_dir/\${o}_kmer_occ.bin >/dev/null 2>&1

/usr/bin/time \$sub_bin/meta_platanus_phase phase -tmp \$tmp -k \${o}_kmer_occ.bin -c \${o}_contig.fa \${o}_junctionKmer.fa $file_input_str -t \$t -o \$o 2>\${o}.phaseLog >\${o}.stdout
/usr/bin/time \$sub_bin/meta_platanus_phase consensus -tmp \$tmp -t \$t -c \${o}_primaryBubble.fa \${o}_nonBubbleOther.fa -k \${o}_kmer_occ.bin $file_input_str -o \$o 2>\${o}.consensusLog >>\${o}.stdout
_EOS

	system("bash cmd.bash >cmd.stdout 2>cmd.stderr") == 0 or die "Error: phase_assembly error.\n$!\n";

	system('touch done');
	chdir('..');
}


sub print_phase_asm_usage
{
	print

<<'USAGE';
Usage:
    meta_platanus.pl phase_asm -IP1 short_R1.fastq(a) short_R2.fastq(a) [Options] ...

Options:
    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id inward_pair_files (reads in 2 files, fasta or fastq; at least one library required)
    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq; aka mate-pairs or jumping-library)
    -p FILE1 [FILE2 ...]               : PacBio long-read file (fasta or fastq)
    -ont FILE1 [FILE2 ...]             : Oxford Nanopore long-read file (fasta or fastq)
    -x PAIR1 [PAIR2 ...]               : barcoded_pair_files (10x Genomics) (reads in 1 file, interleaved, fasta or fastq)
    -X FWD1 REV1 [FWD2 REV2 ...]       : barcoded_pair_files (10x Genomics) (reads in 2 files, fasta or fastq)
    -t INT                             : number of threads (<= 1; default, 1)
    -m INT                             : memory limit for making kmer distribution (unit, GB; default, 0.75 * available_memory))
    -o STR                             : prefix of output files (default "out")
    -tmp DIR                           : directory for temporary files (default, ".")
    -sub_bin DIR                       : directory for sub-executables, such as mata_plantaus and minimap2 (default, directory-of-this-script/sub_bin)
    -overwrite                         : overwrite the previous results, not re-start (default, off)
    -h, -help                          : display usage

USAGE
}




sub get_n_cpu
{
	my $n_cpu = 1;

	$n_cpu = int(`getconf _NPROCESSORS_ONLN`);

	return($n_cpu);
}


sub get_mem_total
{
	my $mem = 1;
	my $in;

	open($in, "/proc/meminfo");
	while (<$in>) {
		if (/MemTotal:\s+(\d+)/) {
			$mem = $1 / 1024 / 1024;
		}
	}
	close($in);

	return(int($mem + 0.5));
}


sub rel2abs_multi
{
	my $files_ref = shift;

	if (ref($files_ref) ne 'ARRAY') {
		return;
	}
	
	for (@{$files_ref}) {
		if (defined($_)) {
			$_ = File::Spec->rel2abs($_);
		}
	}
}

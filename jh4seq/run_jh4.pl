#/usr/bin/perl
use warnings;
use strict;

use 5.26.1;
use Data::Dumper;
use XML::Hash::XS;
use List::Util;
use Bio::SeqIO;

### edit the following ###
my $design_root = '/camp/home/hungm/working/Matthew/project/oscar/OA_DN19023/input/20240606_OA_DN19023.csv';
my $script_root = '/camp/home/hungm/working/Matthew/project/oscar/OA_DN19023/output/20240618/';
my $strictness = 'default'; #xstrict, ustrict, vstrict, strict, default, loose, vloose, uloose, xloose
my $reference_fasta = '/camp/home/hungm/scratch/hungm/reference/alignment/JH4_ref/GRCm39_C57BL6J.fasta';
my $reference_id = '12:113428106-113428538';
#my $min_overlap = '25';
#my $min_overlap0 = '20';

# dont use kmer because this is amplicon sequencing with low diversity
# reduce strictness when data is low quality (often happens in low-diversity amplicon sequencing using long reads)

#########################

system('mkdir '.$script_root.'/scripts') unless (-d $script_root.'/scripts'); # makes a folder that we'll save all the scripts that we'll make unless it exists
system('mkdir '.$script_root.'/slurm_logs') unless (-d $script_root.'/slurm_logs'); # makes a folder that will save all the slurm logs in case we need to check them unless it exists 
system('mkdir '.$script_root.'/merged_reads') unless (-d $script_root.'/merged_reads'); # makes a folder that will save all the slurm logs in case we need to check them unless it exists 
system('mkdir '.$script_root.'/merged_read_fastq') unless (-d $script_root.'/merged_read_fastq'); # makes a folder that will save all the slurm logs in case we need to check them unless it exists 
system('mkdir '.$script_root.'/collapsed') unless (-d $script_root.'/collapsed'); # makes a folder that will save all the slurm logs in case we need to check them unless it exists 
system('mkdir '.$script_root.'/mafft') unless (-d $script_root.'/mafft'); # makes a folder that will save all the slurm logs in case we need to check them unless it exists 


## collect software versions

my $software_versions_file = '/camp/home/hungm/template_pips/jh4/software_versions.csv';  ## this file lists the different camp modules that we'll need.  Here we just check it exists (that I created it properly)
if (!-e $software_versions_file) {
	say  "Software mapping file must exist!"; 
	exit; 
};

open (my $sfh,'<',$software_versions_file);  # opening that file
chomp (my @software_array = <$sfh>); # read the content into an array, which is like a vector in R (or a one dimensional list)
close $sfh; # close that file
shift @software_array;  # shift removes the first element of the array, it gets rid of the header
my $software_map = {}; # this creates a hash, which is a datastructure like a key->value pair
for my $lines (@software_array) { # loops over each of the elements of the array
	if ($lines =~ /\w/) {  # checks that there is a character in that line (we don't want blank lines)
		my ($param, $value) = split(",",$lines); # splits the line into two values, the param (which is the module name) and the value (which is the module definition on CAMP)
		$software_map->{$param} = $value; # builds our key value pair hash that we'll use later
	};
};

 
open (my $dfh, '<', $design_root); # reads in the design file we made
chomp (my @design_array = <$dfh>); # slurps it into an array, 
close $dfh; # close the file we opened

shift @design_array; # gets rid of the header

## here we are assuming that each sample is only present once in the list, as in, we don't need to merge any bam files.  This is fine for now, but if these were re-sequenced / sequenced on something else than a MiSeq this wouldn't work

## but now we know that this isn't always going to be the case, so we need to fix this, maybe we can just zcat the files together

## let's check for duplicate sample_labels

my $counter = {};

for my $lines (@design_array) {
	my ($sample_name, $fastq_1, $fastq_2) = split(',', $lines); # each line has 3 elements, separated by a ,.  here we split them out so we can access them independently
	$counter->{$sample_name}->{'count'}++;
	if (exists $counter->{$sample_name}->{'fq_1'}) {
		my $arrayref = $counter->{$sample_name}->{'fq_1'};
		my @array = @$arrayref;
		push(@array, $fastq_1);
		$counter->{$sample_name}->{'fq_1'}= \@array;
	} else {
		$counter->{$sample_name}->{'fq_1'} = [$fastq_1];
	};
	if (exists $counter->{$sample_name}->{'fq_2'}) {
		my $arrayref = $counter->{$sample_name}->{'fq_2'};
		my @array = @$arrayref;
		push(@array, $fastq_1);
		$counter->{$sample_name}->{'fq_2'}= \@array;
	} else {
		$counter->{$sample_name}->{'fq_2'} = [$fastq_2];
	};
};

say Dumper $counter;

say 'if there are duplicate sample names we need to deal with these';

<STDIN>;

my $map_queue = {};  ## here we're just setting up another key value map that will help us keep track of the jobs we submitted to CAMP to see if they're done or not
my $merge_queue = {};  ## here we're just setting up another key value map that will help us keep track of the jobs we submitted to CAMP to see if they're done or not

print Dumper(\@design_array);

for my $lines (@design_array) {
	if ($lines =~ /\w/) {
		my ($sample_name, $fastq_1, $fastq_2) = split(',', $lines); # each line has 3 elements, separated by a ,.  here we split them out so we can access them independently
		merge_reads($sample_name, $fastq_1, $fastq_2);
	};
};

sleep 10 if (scalar( keys %{$merge_queue} ) > 0); # if there are any keys in the queue hash (as in, if we're doing anything in this step), sleep for 10 seconds
check_status($merge_queue);

## next, run fastx collapser

my $fastx_queue = {};

for my $lines (@design_array) {
	my ($sample_name, $fastq_1, $fastq_2) = split(',', $lines); # each line has 3 elements, separated by a ,.  here we split them out so we can access them independently
	system('zcat '.$script_root.'/merged_reads/'.$sample_name.'.fastq.gz > '.$script_root.'/merged_read_fastq/'.$sample_name.'.fastq');
	fastx($sample_name);
};
sleep 10 if (scalar( keys %{$fastx_queue} ) > 0); # if there are any keys in the queue hash (as in, if we're doing anything in this step), sleep for 10 seconds
check_status($fastx_queue);

# next, run MAFFT

my $mafft_queue = {};

for my $lines (@design_array) {
	my ($sample_name, $fastq_1, $fastq_2) = split(',', $lines); # each line has 3 elements, separated by a ,.  here we split them out so we can access them independently
	mafft($sample_name);
};

sleep 10 if (scalar( keys %{$mafft_queue} ) > 0); # if there are any keys in the queue hash (as in, if we're doing anything in this step), sleep for 10 seconds
check_status($mafft_queue);


my $results = {};
my @res = ('sample	rank	count	snps	deletions	insertions	total	length	sequence	conversions');																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																											

for my $lines (@design_array) {
	my $sequences = {};
	my ($sample_name, $fastq_1, $fastq_2) = split(',', $lines); # each line has 3 elements, separated by a ,.  here we split them out so we can access them independently
	my $inseq = Bio::SeqIO->new( -file => $script_root.'/mafft/'.$sample_name.'.aln', -format => 'fasta' );
	while (my $seq = $inseq->next_seq) {
   		$sequences->{$seq->primary_id}->{'sequence'} = $seq->seq;
	};
	my $reference = $sequences->{$reference_id}->{'sequence'};
	delete $sequences->{$reference_id};
	my @array = split('', $reference);
	for my $keys (sort keys %{$sequences} ) {
		$sequences->{$keys}->{'deletion'} = 0;
		$sequences->{$keys}->{'insertion'} = 0;
		$sequences->{$keys}->{'snp'} = 0;
		$sequences->{$keys}->{'variants'}->{'A:T'} = 0;
		$sequences->{$keys}->{'variants'}->{'A:C'} = 0;
		$sequences->{$keys}->{'variants'}->{'A:G'} = 0;
		$sequences->{$keys}->{'variants'}->{'T:A'} = 0;
		$sequences->{$keys}->{'variants'}->{'T:C'} = 0;
		$sequences->{$keys}->{'variants'}->{'T:G'} = 0;
		$sequences->{$keys}->{'variants'}->{'C:A'} = 0;
		$sequences->{$keys}->{'variants'}->{'C:G'} = 0;
		$sequences->{$keys}->{'variants'}->{'C:T'} = 0;
		$sequences->{$keys}->{'variants'}->{'G:A'} = 0;
		$sequences->{$keys}->{'variants'}->{'G:C'} = 0;
		$sequences->{$keys}->{'variants'}->{'G:T'} = 0;
		my @test = split('', $sequences->{$keys}->{'sequence'});
		my $i = 0;
		my $len = 0;
		for my $position (@array) {
			$len++ if ($test[$i] =~ /\w/);
			if ($test[$i] eq $position) {
		#		say 'positions match! nothing to do!';
			} else {
				$sequences->{$keys}->{'deletion'}++ if ($position =~ /\w/ && $test[$i] eq '-');
				$sequences->{$keys}->{'insertion'}++ if ($position eq '-' && $test[$i] =~ /\w/);
				$sequences->{$keys}->{'snp'}++ if ($position =~ /\w/ && $test[$i] =~ /\w/);
				$sequences->{$keys}->{'variants'}->{'A:T'}++ if( $position eq 'a' && $test[$i] eq 't');
				$sequences->{$keys}->{'variants'}->{'A:C'}++ if( $position eq 'a' && $test[$i] eq 'c');
				$sequences->{$keys}->{'variants'}->{'A:G'}++ if( $position eq 'a' && $test[$i] eq 'g');
				$sequences->{$keys}->{'variants'}->{'T:A'}++ if( $position eq 't' && $test[$i] eq 'a');
				$sequences->{$keys}->{'variants'}->{'T:C'}++ if( $position eq 't' && $test[$i] eq 'c');
				$sequences->{$keys}->{'variants'}->{'T:G'}++ if( $position eq 't' && $test[$i] eq 'g');
				$sequences->{$keys}->{'variants'}->{'C:A'}++ if( $position eq 'c' && $test[$i] eq 'a');
				$sequences->{$keys}->{'variants'}->{'C:G'}++ if( $position eq 'c' && $test[$i] eq 'g');
				$sequences->{$keys}->{'variants'}->{'C:T'}++ if( $position eq 'c' && $test[$i] eq 't');
				$sequences->{$keys}->{'variants'}->{'G:A'}++ if( $position eq 'g' && $test[$i] eq 'a');
				$sequences->{$keys}->{'variants'}->{'G:C'}++ if( $position eq 'g' && $test[$i] eq 'c');
				$sequences->{$keys}->{'variants'}->{'G:T'}++ if( $position eq 'g' && $test[$i] eq 't');
			};
			$i++;
		};
		$sequences->{$keys}->{'length'} = $len;
	};	
	for my $keys ( sort { ($a =~ /(\d+)/)[0] <=> ($b =~ /(\d+)/)[0] } keys %{$sequences}) { 
		my ($rank, $count) = split('-', $keys);	
		my @mut_array = ();
		for my $mut (sort keys %{$sequences->{$keys}->{'variants'}}) {
			push(@mut_array, $mut.'='.$sequences->{$keys}->{'variants'}->{$mut});
		};
		my $total = $sequences->{$keys}->{'snp'} + $sequences->{$keys}->{'deletion'} + $sequences->{$keys}->{'insertion'};
		push(@res,$sample_name.'	'.$rank.'	'.$count.'	'.$sequences->{$keys}->{'snp'}.'	'.$sequences->{$keys}->{'deletion'}.'	'.$sequences->{$keys}->{'insertion'}.'	'.$total.'	'.$sequences->{$keys}->{'length'}.'	'.$sequences->{$keys}->{'sequence'}.'	'.join(';', @mut_array));
	};
};

print join("\n", @res);

open(OF, '>', $script_root.'/mutation_counts.tsv');
print OF join("\n", @res);
close OF;


## subs 


sub mafft {
	my $sample = shift;
	unless (-e $script_root.'/mafft/'.$sample.'.aln') {
		my $align = 'mafft --auto --add '.$script_root.'/collapsed/'.$sample.'.fasta '.$reference_fasta.' > '.$script_root.'/mafft/'.$sample.'.aln';
		my $script = '#!/usr/bin/bash'."\n".'module purge'."\n"."\n".'module load '.$software_map->{'mafft'}."\n"."\n".$align;
		open(OF, '>', $script_root.'/scripts/'.$sample.'.mafft.sh');
		say OF $script;
		close OF;
		my $sbatch = 'sbatch -c 8 --mem=36G -t 1000 -o '.$script_root.'/slurm_logs/'.$sample.'.mafft.out '.$script_root.'/scripts/'.$sample.'.mafft.sh';
		my $cmd = `$sbatch`;
		chomp $cmd;
		$cmd =~ s/Submitted batch job //g;
		$mafft_queue->{$sample.'-mafft'} = $cmd;
	};
}

sub fastx {
	my $sample = shift;
	unless (-e $script_root.'/collapsed/'.$sample.'.fasta') {
		say $sample.' to collapse';
		my $collapse = 'fastx_collapser -v -Q 33 -i'.$script_root.'/merged_read_fastq/'.$sample.'.fastq -o '.$script_root.'/collapsed/'.$sample.'.fasta';
		my $script = '#!/usr/bin/bash'."\n".'module purge'."\n"."\n".'module load '.$software_map->{'fastx'}."\n"."\n".$collapse;
		open(OF, '>', $script_root.'/scripts/'.$sample.'.fastx.sh');
		say OF $script;
		close OF;
		my $sbatch = 'sbatch -c 4 --mem=4G -t 1000 -o '.$script_root.'/slurm_logs/'.$sample.'.fastx.out '.$script_root.'/scripts/'.$sample.'.fastx.sh';
		my $cmd = `$sbatch`;
		chomp $cmd;
		$cmd =~ s/Submitted batch job //g;
		$fastx_queue->{$sample.'-fastx'} = $cmd;
	};
};


sub merge_reads {
	my ($sample, $f1, $f2) = @_;
	unless (-e $script_root.'/merged_reads/'.$sample.'.fastq.gz') {
		say 'running merging things for '.$sample; # just prints this out to the screen to let us know things are moving
		my $merge = 'bbmerge.sh in='.$f1.' '.'in2='.$f2.' '.'out='.$script_root.'/merged_reads/'.$sample.'.fastq.gz '.$strictness;
		#my $merge = 'bbmerge.sh in='.$f1.' '.'in2='.$f2.' '.'minoverlap='.$min_overlap.' '.'minoverlap0='.$min_overlap0.' '.'out='.$script_root.'/merged_reads/'.$sample.'.fastq.gz '.$strictness;
		my $script = '#!/usr/bin/bash'."\n".'module purge'."\n"."\n".'module load '.$software_map->{'conda'}."\n"."\n".'eval "$(conda shell.bash hook)"'."\n".'conda activate bbtools'."\n".$merge;
		open(OF, '>', $script_root.'/scripts/'.$sample.'.merge.sh');
		say OF $script;
		close OF;
		my $sbatch = 'sbatch -c 2 --mem=8G -t 240 -o '.$script_root.'/slurm_logs/'.$sample.'.merge.out '.$script_root.'/scripts/'.$sample.'.merge.sh';
		my $cmd = `$sbatch`;
		chomp $cmd;
		$cmd =~ s/Submitted batch job //g;
		$merge_queue->{$sample.'-merge'} = $cmd;
	};

};




### general subroutines that we use just for checking CAMP.  I won't comment all this as it's not really important, it just helps us monitor the job status (queued, running, failed etc) on CAMP


sub check_status {
	my $queue = shift;
	my $completed = {};
	my $failed = {};
	do {
		for my $keys (keys %{$queue}) {
			my $check = 'sacct -j '.$queue->{$keys}.' |grep -v "^[0-9]*\."';
			my $check_cmd = `$check`;
			my @check_array = split("\n", $check_cmd);
			if (scalar(@check_array) > 1) {  ## if there is a row for this job id
				my $row = $check_array[2];
				$row =~ s/\s+/ /g;
				my @row_array = split(' ', $row);
				say $keys.' : '.$row_array[5];
				if ($row_array[5] ne 'PENDING') {
					if ($row_array[5] ne 'RUNNING' ) {
						if ($row_array[5] eq 'COMPLETED') {
							$completed->{$keys}++;
							delete $queue->{$keys};
						} else {
							$failed->{$keys}++;
						delete $queue->{$keys};
						};
					};
				};
			};
		};
		say Dumper $queue;
		sleep 60 if (scalar( keys %{$queue} ) > 0);;
	} until (scalar( keys %{$queue} ) == 0);
	say 'complete:'. Dumper $completed;
	say 'failed:'. Dumper $failed;
	if (scalar (keys %{$failed}) > 0 ) {
		say 'pipeline encountered failed processes';
		exit;
	};
};



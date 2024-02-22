use strict;

my $matrix_file = $ARGV[0]; #The matrix file

#For this script to work, the matrix file must look like this:
#	genome1	genome2	genome3
#genefamilyA	1	2	2
#genefamilyB	1	0	1
#genefamilyC	2	0	1
#genefamilyD	0	1	1


#The first task, counting the genomes

my $header= `head -n1 $matrix_file`;
my @name_genomes = split(/\t/,$header);
my $number_genomes = scalar(@name_genomes)-1;
#print("$colnum\n");
my %score_ubiquity;
my %mode_modifier;

print("Score for every gene family/gene cluster:\n");
open (FILE,$matrix_file) or die $!;
	foreach my $line (<FILE>){
		my @genefam_idNvalues;
		#print($line);
		#The next line discards the first line, which is the header with the genome names
		#We just need the values for the gene families/clusters not the names of the genomes
		if ($line!~/^\t.*/){
			@genefam_idNvalues=split(/\t/,$line);
			#print("$genefam_idNvalues[0]\n");
			my $weigth=0;
			my $counter=1;
			for my $i (1..$number_genomes){
				#print("$i\t$genefam_idNvalues[$i]\n");
				#These line are for assign a weight for each row according to the value in each row. If a row has a value of more than 0 in each row, its weigth is 1. The value decreases with each 0.
				if ($genefam_idNvalues[$i]>0){
					$weigth=$weigth+$counter;
				}	
			}
			my $percentweigth=$weigth/$number_genomes;
			$score_ubiquity{$genefam_idNvalues[0]}=$percentweigth;
			print("$genefam_idNvalues[0]\t$score_ubiquity{$genefam_idNvalues[0]}\n");
			#Now let's get the mode and frequencies of each value and the penalties for each one.
			my %freq;
			my @numberset=splice(@genefam_idNvalues,1);
			foreach (@numberset) {$freq{$_}++};
			my @sorted_array = sort { $freq{$a} <=> $freq{$b} } keys %freq;
			my $mode = pop(@sorted_array);
			foreach (sort keys %freq){
				my $modevalue=1-(abs($freq{$_}-$freq{$mode})/(($freq{$_}+$freq{$mode})/2));
				$mode_modifier{$genefam_idNvalues[0]}{$_}=$modevalue;
			}			
		}
	}
close FILE;

my $jobID=`head -n1 $matrix_file | tr '\t' '\n' | sed '/^\$/d'`;
my @rastID=split(/\n/,$jobID);
#print("@rastID\n");
foreach my $ID(@rastID){
	my $value=0;
	my @filas=(1);
	my $number =`head -n1 $matrix_file| tr '\t' '\n' | cat -n | grep -w $ID | cut -f1 | sed -e 's/^[ \t]*//'`;
	chomp $number;
	push (@filas,$number);
	my $cutFeed=join(",",@filas);
	my $minifile=`cut -f $cutFeed $matrix_file | tail -n +2`;
	my @arrayminifile=split(/\n/,$minifile);
	#print("$arrayminifile[0]");
	foreach my $llave (sort keys %score_ubiquity){	
		my @extract = (grep { $_ =~ /$llave\t/ } @arrayminifile);
		#print("$llave\t@extract\n");
		#print("$extract[0]\n");		
		my @breaks=split("\t",$extract[0]);
		#print("$breaks[1]\n")
		#print("$score_ubiquity{$llave}\t$breaks[1]\n");
		$value=$value+($score_ubiquity{$llave}*$mode_modifier{$breaks[0]}{$breaks[1]});		
	}
	print("$ID\t$value\n");	
}




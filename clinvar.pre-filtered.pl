#!/usr/bin/env perl
use warnings;
use List::MoreUtils qw/ uniq /;
use Sort::Fields;

################# opening ##################
open CLINVAR, "<./variant_summary.txt" or die "$!";
# Poner fecha
open OUT, ">./clinvar.feb19" or die "$!";
open ERROR, ">./variant_summary.error" or die "$!";

my @clinvar = <CLINVAR>;
@file = grep /\tGRCh37\t/, @clinvar;

################# clean Clinvar ##################

foreach my $h (0..$#file) {
	chomp($h);
	my ($id, $type, $mut, $gene, $chr, $ini, $fin, $ref, $alt, $dis, $origin, $sig, $rev) = (split /\t/, $file[$h])[0, 1, 2, 4, 18, 19, 20, 21, 22, 13, 14, 6, 24]; # mutacion, gen
	if ($h == 0) {
		print OUT "GeneSymbol\tType\tClinicalSignificance\tOriginSimple\tPhenotypeList\tName\tcvid\tkey\ttranscrip\tconsequence\tref_aa\talt_aa\tpos_aa\treview\n";
		next;
	}
	# defitions of extras
	my $consequence = "NA";
	my $aa_ref = "NA";
	my $pos_aa = "NA";
	my $aa_alt = "NA";
	my $transcript = "NA";
	my $var = "NA";
	my $var2 = "NA";
	# Estructura de la entrada Mut.
	if ($mut =~ /:/) {
		($transcript, $var) = split /:/, $mut, 2;
		$transcript =~ s/\(.*\)//g;
		if ($var =~ /\(*p\./) {
			$var =~ s/.*\(p\.//g;
			$var =~ s/\)//g;
			$var =~ s/\*/Ter/g;
			$aa_ref = $var;
			$aa_ref =~ s/[0-9]+.*//g;
			$var2 = $var;
			$var2 =~ s/[0-9]+/\t/;
			$aa_alt = (split /\t/, $var2)[1];
			$pos_aa = $var;
			$pos_aa =~ s/[a-zA-Z=]+//g;
			if (! defined $aa_ref ) {
				$aa_ref = "NA";
			} 
			if (! defined $aa_alt) {
				$aa_alt = "NA";
			} 
			if ($aa_alt eq "=") {
				$aa_alt = $aa_ref;
				$consequence = "Synonyomous";
			} elsif ($aa_alt =~ "fs") {
				$consequence = "Frameshift";
			} elsif ($aa_alt =~ "del|ins|dup") {
				$consequence = "In frame indel"; 
			} elsif ($aa_alt ne $aa_ref) {
				$consequence = "Missense";
			} 
			if ($aa_alt =~ "Ter") {
				$consequence = "Stop gain";
			}
		} elsif ($var =~ /\Ac\.-/) {
			if ($var =~ /\Ac\.-[0-9]+\+[123][ACTG]/) {
				$consequence = "Start loss";
			} else {
				$consequence = "5-UTR";				
			}
		} elsif ($var =~ /\Ac\.\*/) {
			if ($var =~ /\Ac\.\*[0-9]+-[123][ACTG]/) {
				$consequence = "Stop loss";
			} else {
				$consequence = "3-UTR";				
			}
		} elsif ($var =~ /\Ac\.[0-9]/) {
			if ($var =~ /\Ac\.[0-9]+\+[12][ACTG]/) {
				$consequence = "Splice-D/A";
			} elsif ($var =~ /\Ac\.[0-9]+-[12][ACTG]/) {
				$consequence = "Splice-D/A";
			} else {
				$consequence = "Intronic";				
			}
		} elsif ($var =~ /\Ag\..*del/) {
			$consequence = "NA";
		} elsif ($var =~ /\Ag\..*dup/) {
			$consequence = "NA";
		} elsif ($type =~ /\Acopy number/) {
			$consequence = "NA";
			$pos_aa = $var;
		} elsif ($var =~ /\Ag\./) {
			$consequence = "Genomic";
			$pos_aa = $var;
		} elsif ($var =~ /\An\./) {
			$consequence = "Non-coding";
		} elsif ($var =~ /\Am\./) {
			$consequence = "Mitocondrial";
		} elsif ($type =~ /\Adeletion/) {
			$consequence = "NA";
			$pos_aa = $var;
		} elsif ($type =~ /\Aduplication/) {
			$consequence = "NA";
			$pos_aa = $var;
		} elsif ($type =~ /\Acomplex/) {
			$consequence = "NA";
			$pos_aa = $var;
		} else {
				print ERROR $type."\t".$var."\n";
		}

	} else {
			print ERROR "===>".$type."\t".$var."\n";
	}
	my $key = $chr.";".$ini.";".$fin.";".$ref.";".$alt;
	print OUT "$gene\t$type\t$sig\t$origin\t$dis\t$mut\tCV:$id\t$key\t$transcript\t$consequence\t$aa_ref\t$aa_alt\t$pos_aa\t$rev\n";
#                0        1        2         3      4          5          6        7        8       9      10    11               
}
close CLINVAR;
close OUT;
close ERROR;


# Agregar serie de gsubs que marie hizo con el otro codigo de R. 
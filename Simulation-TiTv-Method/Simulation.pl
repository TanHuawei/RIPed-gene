use strict;

# Author:		TAN, Huawei
# Affiliation:	GATTACA Lab, Nanjing University
# Email:		tanhuawei1991@163.com

use Data::Dump qw(dump);
use List::Util 'shuffle';
use IO::Handle;
use POSIX qw(strftime);


# 2025-04-02
# 2025-04-05 Update
# 2025-04-12 Update
# 2025-04-21 Update
# 2025-05-01 Update to SimuW4
# 2025-05-04 Update  GeneScan DupReg RIPed deRIP
# 2025-05-12 Update and Integrate statistics
# 2025-05-13 Update to SimuW6$rea
# 2025-05-20 Update to SimuW7
# 2025-05-28 Update to SimuW8
# 2025-06-03 Simulation-V1.pl

# If you feel hard to read this Perl script, you can just follow the Workflow.sh to run

my ($RunPart, $SpeciesSet, $prefix, $pct, $repiSet, $folderout) = @ARGV;
# $RunPart	IntroSNP  Stat
# $repi is the counter of replicates. Also set 'Each2' means only set $repi = 2; Set '5' means set $repi = 2, 3, 4, 5 in Batch mode

# Here we set the parameters of RIP-type SNPs ratio.
# For example, 20 means 20% RIP-type SNPs, with 80% other type SNPs.
my @round = qw/0 20 40 60 80 100/;

# when $pct == 0, we only calculate the RIP index, but not introduce SNPs
if ( $pct == 0 ) { @round = qw/0/; } 

# When we run a large number of replicated, set $BatchNum = 5 or other numbers could accelerate the speed.
my $BatchNum = 5;
my @repi;
if ( $repiSet=~/^Each(\d+)$/ ) {
	@repi = ($1);
} elsif ( $repiSet=~/^(\d+)$/ ) {
	if    ( $repiSet  < $BatchNum ) { @repi = ($repiSet); }
	elsif ( $repiSet == $BatchNum ) { @repi = 2 .. $BatchNum; }
	elsif ( $repiSet >= $BatchNum*2 and $repiSet % $BatchNum == 0 ) { @repi = $repiSet - $BatchNum + 1 .. $repiSet; }
	else { @repi = ($repiSet); }
}


# Set the parameters
my $RIPRateSet      = $pct / 10 ** 2;
my $clusterlen		= 1000;
#my $RIPmodelfile 	= 'SimuW1.RIPmutNondupModel.ini';
my $RIPmodelfile 	= 'Simu.RIPMod.Tetrad4Spores.nondup.ini';
my $othermodelfile 	= 'SimuW2.MutDirt.g6722.ini';

# Nature 2003 paper parameters
my %BinucProb = qw/Ca 0.3 tG 0.3 Ct 0.05 aG 0.05 Cg 0.01 cG 0.01 Cc 0.009 gG 0.009/;
if ( $pct >= 15 ) { %BinucProb = qw/Ca 0.3 tG 0.3 Ct 0.1 aG 0.1 Cg 0.1 cG 0.1 Cc 0.05 gG 0.05/; }
elsif ( $pct >= 20 ) { %BinucProb = qw/Ca 0.3 tG 0.3 Ct 0.15 aG 0.15 Cg 0.15 cG 0.15 Cc 0.1 gG 0.1/; }

my ($InitFasFile, $InitGeneFasFile, $cdscoofile, $geneusefile, $repeatcoofile, $repeatusefile, $ChrHead, $UseTotalLen);


# This simulation is work on Sordaria macrospora, a RIP-lack fungi, and we also provide other choices
if     ( $SpeciesSet eq 'SorMacSN1693' ) { 
	$InitFasFile		= 'SorMacSN1693.genome.fa.gz';
	$cdscoofile 		= 'SorMacSN1693.cds.coo';
	$repeatcoofile		= 'SorMacSN1693.dupblock2.coo';
	$geneusefile		= 'SorMacSN1693-ortho.gene.coo';
	$repeatusefile		= 'SorMacSN1693.dupblock2.coo';
	$ChrHead			= 'CP083';
} elsif   ( $SpeciesSet eq 'FGSC2489' ) { 
	$InitFasFile		= 'thwNC12v42-contig20Mt.fas';
	$cdscoofile			= 'thwNC12v42-cds-coo.txt';
	$geneusefile 		= 'thwNC12v42-cds+intron.nuclear.coo';
	$repeatcoofile		= 'FGSC2489.dupblock2.coo';
	$repeatusefile		= 'FGSC2489.dupblock2.coo';
	$ChrHead			= 'Supercontig';
} elsif ( $SpeciesSet eq 'FGSC2225' ) { 
	$InitFasFile		= 'ref2225v03.nuclear.fa';
	$cdscoofile 		= 'ref2225v03-cds-coo.txt';
	$geneusefile		= 'ref2225v03.gene.nuclear.coo';
	$repeatcoofile		= 'FGSC2225.dupblock2.coo';
	$repeatusefile		= 'FGSC2225.dupblock2.coo';
	$ChrHead			= 'Chr';
} else { die "That species not found\n"; }


if ( $RunPart eq 'GeneScan' ) {

	# Initiate the Genome from the fasta file
	my %GeneFas = &InitFasV0($InitGeneFasFile, '', 'Hash');
	open STDOUT, ">$folderout/0-$SpeciesSet.gene.txt";
	open STDERR, ">$folderout/0-$SpeciesSet.region.txt";

	print STDOUT join ("\t" => qw/Project	Species	GeneID	GeneLength	GeneNumCA	GeneNumTG	GeneNumCATG	GeneCATGdivergence	RIPindex	RegLen	RegCATGcontain	RegCATGNeedMinGene	RegCATGdivergence	RegCATGNeedMinGeneRatio/), "\n";
	my $repi = $pct;
#	print STDERR join ("\t" => $repi, ''), "\n";
	my (%RIPindexThisHash, %MutCountThisHash);
	for my $gene ( sort keys %GeneFas ) {
		my $genenewseq = $GeneFas{$gene};
#		print STDERR join ("\t" => 'genenewseq', $gene, $genenewseq), "\n";
		my ($typereg, $RegLen, $RegCATGcontain, $RegCATGNeedMinGene, $RIPinfoHashRF) = &RIPindexV8GeneScan($genenewseq, 1000, $gene);
		$RIPindexThisHash{$gene} = $typereg;

		my ($numCA, $numTG);
		my $numN = $genenewseq =~tr/N/N/;
		if ($genenewseq=~/CA/) { $numCA = () = $genenewseq=~/CA/g; } else { $numCA = 0; }
		if ($genenewseq=~/TG/) { $numTG = () = $genenewseq=~/TG/g; } else { $numTG = 0; }
		my $numCATG = $numCA + $numTG;
		my $genenewlen = length($genenewseq);
		my $CATGdivergence = $numCATG/$genenewlen;

		my $RegCATGdivergence = 'NA';			if ( $RegLen > 0 )         { $RegCATGdivergence       = $RegCATGNeedMinGene/$RegLen; }
		my $RegCATGNeedMinGeneRatio = 'NA';		if ( $RegCATGcontain > 0 ) { $RegCATGNeedMinGeneRatio = $RegCATGNeedMinGene/$RegCATGcontain; }
		print STDOUT join ("\t" => $repi, $SpeciesSet, $gene, $genenewlen, $numCA, $numTG, $numCATG, $CATGdivergence, $typereg, $RegLen, $RegCATGcontain, $RegCATGNeedMinGene, $RegCATGdivergence, $RegCATGNeedMinGeneRatio ), "\n";
		STDOUT->flush(); STDERR->flush();
	}
	close STDOUT; close STDERR;
	exit;
}


elsif ( $RunPart =~/IntroSNP/ ) {


################################################################################ Eliminate SNPs that cause pre-stop of peptides
	my %NewStopHash;
	open FILE, "gzip -dc SorMacSN1693-prestop.txt.gz | " or die;
	while (<FILE>) { s/\r//g; chomp;
		my ( $gene, $iref, $chr, $site, $refbase, $altbase, $dnacodonref, $dnacodonalt, $genetic_code_ref, $genetic_code_alt ) = split /\t/;
		if ( $genetic_code_ref ne '*' and $genetic_code_alt eq '*' ) {
			$NewStopHash{$chr}{$site}{$altbase} = 'prestop';
		}
	}
	close FILE;

################################################################################ Initiate the RIP SNP Model

	my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime; print STDERR "START $current_time\n";

	my ($UseGeneHashRF, $UseRepeatHashRF, $UseAllHashRF) = &GeneUseInitV8($geneusefile, $repeatusefile);
	my ($HashClassRF, $UseDupLen)   = &HashClassdupV8($repeatcoofile);
	my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime; print STDERR "HashClassRF Init done at $current_time\n";

	# Initiate the Genome from the fasta file
	my %GenomeFas = &InitFasV0($InitFasFile, '', 'Hash');

	my ($GenomeEachRF, $HashClassRF, $BinucArrHashRF, $BinucCountHashRF, $chrsiteclass3HashRF, $chrsiteArrhashRF, $UseTotalLen ) = &GenomeEachV8(\%GenomeFas, $HashClassRF, 'IntroSNP');
	my $chrsiteclass3ArrRF;
	my ($genecoohashRF, $cdscoohashRF, $cdssitehashRF, $site2genehashRF ) = &InitGeneCOOV8($UseAllHashRF, $cdscoofile, $repeatcoofile);
	my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime; print STDERR "InitGeneCOOV8 Init done at $current_time\n";

	my @UseGeneArr   = sort keys %$UseGeneHashRF;
	my @UseRepeatArr = sort keys %$UseRepeatHashRF;

	my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	print STDERR "Total ", scalar @UseGeneArr, " genes. ", scalar @UseRepeatArr, " Repeats. Init done $current_time\n";
################################################################################ Initiate END

#### Set min SNPs number in each cluster, only useful in cls InSNP method
	my $MinSNPInOneCluster = 3;
	my %modetetradhash;
	my @LenSNPArrMode;
	open MODE, "<$RIPmodelfile" or die;
	while (<MODE>) { s/\r//g; chomp;
		my ( $gid, $countall, $DUP, $GFullLen, $sitenewARR, $dirt, $countdirt3ARR, $combnewARR ) = split /\t/;
		if ( $countall >= $MinSNPInOneCluster ) { } else { next; }
		my $tetrad = (split /\t/)[0];
	#	push @{ $modetetradhash{$tetrad} }, $_ ;
		push @LenSNPArrMode, "$countall $GFullLen $sitenewARR" ;
	}
	close MODE;


	### Introduce Other SNP information
	my @OtherDirtArr;
	my %DirtSet;
	my $DirtSetTotalNum;
	open MODE, "<$othermodelfile" or die;
	while (<MODE>) { s/\r//g; chomp;
		my ($dirt, $num) = split /\t/;
		$DirtSet{$dirt} = $num;
		$DirtSetTotalNum += $num;
		# Add 10 folds of these mutations for sufficient shuffle
	#	for my $n ( 1 .. $num * 10 ) { push @OtherDirtArr, $dirt; }
		for my $n ( 1 .. $num      ) { push @OtherDirtArr, $dirt; }
	}
	close MODE;


	for my $repi ( @repi ) {

	# Output files
	open SITELIST1S, " | gzip > $folderout/$prefix/${pct}pct.$repi.sitelist.txt.gz";	# short information
	print SITELIST1S join ("\t" => qw/repi round Process ClusterID SNPNumFin chr site refbase altbase dirt selection/), "\n";

		if ( $repi == 1 ) {
			open SITELISTIN, " | gzip > $folderout/$prefix/${pct}pct.$repi.InClusterInfo.txt.gz" if $prefix=~/cluster/ or $prefix=~/^cls/;	# introduced clustered SNPs information
			open SITEINSWAY, " | gzip > $folderout/$prefix/${pct}pct.$repi.InClusterSway.txt.gz" if $prefix=~/^cls[23]/;	# introduced clustered SNPs information

			print SITELISTIN join ("\t" => qw/repi round Process ClusterID chrsiteNum chrsite/), "\n";
		}

	for my $round ( @round ) {

	# Define the number of Simulated RIP type SNPs and other type SNPs
		my $RIPRatioSet     = $round/100;
		my $RIPNumSet       = int( $RIPRateSet * $UseTotalLen *    $RIPRatioSet  + 0.5 );
		my $OthNumSet       = int( $RIPRateSet * $UseTotalLen * (1-$RIPRatioSet) + 0.5 );

	############################################################################################################## First Introduce RIP Model SNPs
	# Prepare the Mode for RIP type SNPs
		my $RIPfold = int($RIPNumSet / scalar @LenSNPArrMode) + 1;
		my @LenSNPArr;
		for my $itmp ( 1 .. $RIPfold ) { push @LenSNPArr, @LenSNPArrMode; }
		@LenSNPArr = shuffle @LenSNPArr;

	############################################################################################################## First Initiate RIP type SNPs
		# Summary %hash 

		my (%UsedChrSiteHash);
		# shuffle and obtain the initial site
		my @chrsiteclass3shuf       = shuffle (@$chrsiteclass3ArrRF); # Actually not run

		my $RIPNumFin  = 0;
		my $RIPshufi   = 0;
		my $Process    = 'InRIP';
		my $ClusterID  = 0;

############################################## introduce clustered model
		if ( $prefix=~/^cls4/ ) {

			my $RIPDupNumSet    = int( $RIPRateSet * $UseDupLen   *    $RIPRatioSet  + 0.5 );
			my $RIPDupNumFin    = 0;

			my $swaySet = 100;
			my @sway = (0);
			for my $i (1 .. $swaySet) { push @sway, -$i, $i; }

			my $RF1 = $BinucArrHashRF->{"Ca dup"};	my $RF2 = $BinucArrHashRF->{"tG dup"};
			my @chrsiteclass3shuf       = shuffle (@$RF1, @$RF2);
			while ( $RIPDupNumFin < $RIPDupNumSet ) {
				$RIPshufi ++;
				my $chrsite = $chrsiteclass3shuf[$RIPshufi-1];
				my ($chr, $site) = split /\s+/, $chrsite;
				my ( $clusterRIPnum, $RIPNewFlankLen, $sitenewARR ) = split /\s+/, $LenSNPArr[$RIPshufi-1];

				my $dinuc = (shuffle qw/CA TG/)[0];
				my $base3;
				my %chrsitenew;

				# Simulate SNPs accroding to the site Model
				for my $sitecoo ( split /\;/, $sitenewARR ) {
					my $flag = 'NotFound';
					my $sitenew;
					my $swaynew;
					while ( $flag eq 'NotFound' ) {
						for my $sway ( @sway ) {
						#	print STDOUT join ("\t" => $mode, $sitecoo, 'flag', $flag, $sway), "\n";
						
							$swaynew = $sway;
							$sitenew = $site + $sway + $sitecoo;
							$base3 = $GenomeEachRF->{$chr}[$sitenew-1].$GenomeEachRF->{$chr}[$sitenew].$GenomeEachRF->{$chr}[$sitenew+1];

							if ( $base3=~/$dinuc/ and length($base3) == 3 and $chrsiteclass3HashRF->{$chr}{$sitenew} =~/$ChrHead/ and not $UsedChrSiteHash{$chr}{$sitenew} ) {
								$chrsitenew{"$chr $sitenew"} ++;
								$flag = 'meet';
							}
							if ( $flag eq 'meet' or abs($swaynew) >= $swaySet ) { last; }
						}
						if ( $flag eq 'meet' or abs($swaynew) >= $swaySet ) { last; }
					}
					print SITEINSWAY join ("\t" => $repi, $round, $dinuc, $base3, $sitenewARR, $sitecoo, $swaynew, $sitenew, $flag), "\n" if $repi == 1;
				}

				my @chrsitenew = sort keys %chrsitenew;
				if ( scalar @chrsitenew >= $MinSNPInOneCluster ) { $ClusterID ++; } else { next; }
				print SITELISTIN join ("\t" => $repi, $round, $Process, $ClusterID, scalar @chrsitenew, @chrsitenew), "\n" if $repi == 1;

				for my $chrsitenew ( @chrsitenew ) {
					$RIPNumFin ++;
					my ($chrnew, $sitenew) = split /\s+/, $chrsitenew;
		#			delete $chrsiteclass3HashRF->{$chr}{$sitenew};
					# store the mutated sites
					my $refbase = $GenomeEachRF->{$chr}[$site];
					my $altbase;	if ( $refbase eq 'C' ) { $altbase = 'T'; } elsif ( $refbase eq 'G' ) { $altbase = 'A'; } else { die join ("\t" => $chrsite, $refbase, $altbase), "\n"; }
					$UsedChrSiteHash{$chr}{$sitenew} = $altbase;

					### simple information
					my $selection = $NewStopHash{$chr}{$sitenew}{$altbase} || 'normal';
					my $dirt = "$refbase>$altbase";
					my $dup = $HashClassRF->{$chr}{$sitenew};
					if ( $dup eq 'dup' ) { $RIPDupNumFin ++; }

					my $print1 = join ("\t" => $repi, $round, $Process, $ClusterID, $RIPNumFin, $chr, $sitenew, $refbase, $altbase, $dirt, $dup, $selection);
					print SITELIST1S "$print1\n";
				}
			} # End Dup CA/TG shuf in Dup


#		my $RIPNonDupNumSet = $RIPNumFin - $RIPDupNumSet + 0;
#			$RIPshufi = 0;  # Initiate again, not used

			my $RF1 = $BinucArrHashRF->{"Ca nondup"};	my $RF2 = $BinucArrHashRF->{"tG nondup"};
			my @chrsiteclass3shuf       = shuffle (@$RF1, @$RF2);

			while ( $RIPNumFin < $RIPNumSet ) {  # Exclude dup mut,then we supply only nondup mut
				$RIPshufi ++;
				my $chrsite = $chrsiteclass3shuf[$RIPshufi-1];
				my ($chr, $site) = split /\s+/, $chrsite;
				my ( $clusterRIPnum, $RIPNewFlankLen, $sitenewARR ) = split /\s+/, $LenSNPArr[$RIPshufi-1];

				my $dinuc = (shuffle qw/CA TG/)[0];
				my $base3;
				my %chrsitenew;

				# Simulate SNPs accroding to the site Model
				for my $sitecoo ( split /\;/, $sitenewARR ) {
					my $flag = 'NotFound';
					my $sitenew;
					my $swaynew;
					while ( $flag eq 'NotFound' ) {
						for my $sway ( @sway ) {
						#	print STDOUT join ("\t" => $mode, $sitecoo, 'flag', $flag, $sway), "\n";
						
							$swaynew = $sway;
							$sitenew = $site + $sway + $sitecoo;
							$base3 = $GenomeEachRF->{$chr}[$sitenew-1].$GenomeEachRF->{$chr}[$sitenew].$GenomeEachRF->{$chr}[$sitenew+1];

							if ( $base3=~/$dinuc/ and length($base3) == 3 and $chrsiteclass3HashRF->{$chr}{$sitenew} =~/$ChrHead/ and not $UsedChrSiteHash{$chr}{$sitenew} and $HashClassRF->{$chr}{$sitenew} eq 'nondup' ) {    # Only supply with nondup sites
								$chrsitenew{"$chr $sitenew"} ++;
								$flag = 'meet';
							}
							if ( $flag eq 'meet' or abs($swaynew) >= $swaySet ) { last; }
						}
						if ( $flag eq 'meet' or abs($swaynew) >= $swaySet ) { last; }
					}
					print SITEINSWAY join ("\t" => $repi, $round, $dinuc, $base3, $sitenewARR, $sitecoo, $swaynew, $sitenew, $flag), "\n" if $repi == 1;
				}

				my @chrsitenew = sort keys %chrsitenew;
				if ( scalar @chrsitenew >= $MinSNPInOneCluster ) { $ClusterID ++; } else { next; }
				print SITELISTIN join ("\t" => $repi, $round, $Process, $ClusterID, scalar @chrsitenew, @chrsitenew), "\n" if $repi == 1;

				for my $chrsitenew ( @chrsitenew ) {
					$RIPNumFin ++;
					my ($chrnew, $sitenew) = split /\s+/, $chrsitenew;
		#			delete $chrsiteclass3HashRF->{$chr}{$sitenew};
					# store the mutated sites
					my $refbase = $GenomeEachRF->{$chr}[$site];
					my $altbase;	if ( $refbase eq 'C' ) { $altbase = 'T'; } elsif ( $refbase eq 'G' ) { $altbase = 'A'; } else { die join ("\t" => $chrsite, $refbase, $altbase), "\n"; }
					$UsedChrSiteHash{$chr}{$sitenew} = $altbase;

					### simple information
					my $selection = $NewStopHash{$chr}{$sitenew}{$altbase} || 'normal';
					my $dirt = "$refbase>$altbase";
					my $dup = $HashClassRF->{$chr}{$sitenew};

					my $print1 = join ("\t" => $repi, $round, $Process, $ClusterID, $RIPNumFin, $chr, $sitenew, $refbase, $altbase, $dirt, $dup, $selection);
					print SITELIST1S "$print1\n";
				}
			} # End Dup CA/TG shuf in Nondup

		}
	############################################################################################################## First Introduce RIP SNPs randomly on genomic CpA


		######################## Introduce novel Mutations only due to substitution rate

		elsif ( $prefix=~/single/ or $prefix=~/^sgl1/ ) {
			while ( $RIPNumFin < $RIPNumSet ) {
					$RIPshufi ++;
					my $chrsite = $chrsiteclass3shuf[$RIPshufi-1];
					my ($chr, $site) = split /\s+/, $chrsite;


					$RIPNumFin ++;
					my ($chrnew, $sitenew) = ($chr, $site);
		#			delete $chrsiteclass3HashRF->{$chr}{$sitenew};
					# store the mutated sites
					my $refbase = $GenomeEachRF->{$chr}[$site];
					my $altbase;	if ( $refbase eq 'C' ) { $altbase = 'T'; } elsif ( $refbase eq 'G' ) { $altbase = 'A'; } else { die join ("\t" => $chrsite, $refbase, $altbase), "\n"; }
					$UsedChrSiteHash{$chr}{$sitenew} = $altbase;
					$ClusterID ++;
		#			print SITELISTIN join ("\t" => $repi, $round, $Process, $ClusterID, 1, $chrsite), "\n";

					### simple information
					my $selection = $NewStopHash{$chr}{$sitenew}{$altbase} || 'normal';
					my $dirt = "$refbase>$altbase";
					my $dup = $HashClassRF->{$chr}{$sitenew};

					my $print1 = join ("\t" => $repi, $round, $Process, $ClusterID, $RIPNumFin, $chr, $sitenew, $refbase, $altbase, $dirt, $dup, $selection);
					print SITELIST1S "$print1\n";
			}
		}

		######################## Introduce novel Mutations due to substitution rate on specific di-nucleotide context. This is used in the paper
                               # The genome sequence of the filamentous fungus Neurospora crassa, Nature, 2003
                               # According to the probabilities in the previous report: (CpA = 0.3, CpT = 0.05, CpG = 0.01, CpC = 0.009) (Galagan et al., 2003), we incorporated these probabilities during each round of simulation involving 1% divergence for RIP (CpA = 0.059843, CpT = 0.009969, CpG = 0.001989, CpC = 0.00179), every cytosine-containing dinucleotide on one strand (selected with equal probabilities) was mutated.

		elsif ( $prefix=~/^sgl3/ ) {
			my $BinucMuttedSiteAllNum = 0;
			for my $dup ( qw/dup nondup/) {
				for my $Binuc ( sort keys %BinucProb ) {
					my $BinucSiteNum = $BinucCountHashRF->{"$Binuc $dup"};
					my $BinucMuttedSiteNum = $BinucSiteNum * $BinucProb{$Binuc};
					$BinucMuttedSiteAllNum += $BinucMuttedSiteNum;
					print STDERR join ("\t" => $repi, $round, $Binuc, $dup, $BinucSiteNum, $BinucMuttedSiteNum ), "\n";
				}
			}
			my $FoldSet = $RIPRateSet * $UseTotalLen * $RIPRatioSet / $BinucMuttedSiteAllNum;
			print STDERR join ("\t" => $repi, $round, 'Total', $BinucMuttedSiteAllNum, $UseTotalLen, $BinucMuttedSiteAllNum / $UseTotalLen, $FoldSet ), "\n";
			for my $dup ( qw/dup nondup/) {
				for my $Binuc ( sort keys %BinucProb ) {
					my ($refbase, $altbase);
					if    ( $Binuc=~/^C[atgc]$/ ) { ($refbase, $altbase) = qw/C T/; }
					elsif ( $Binuc=~/^[atgc]G$/ ) { ($refbase, $altbase) = qw/G A/; }
					my $dirt = "$refbase>$altbase";
					$RIPshufi = 0;
					my $BinucSiteNum = $BinucCountHashRF->{"$Binuc $dup"};
					my $BinucMuttedSiteNumSet = int($BinucSiteNum * $BinucProb{$Binuc} * $FoldSet + 0.5);
					print STDERR join ("\t" => $repi, $round, $Binuc, $dup, $BinucSiteNum, $BinucMuttedSiteNumSet ), "\n";

					$ClusterID = 0;
					my $RF = $BinucArrHashRF->{"$Binuc $dup"};
					for my $chrsite ( (shuffle @$RF)[0 .. $BinucMuttedSiteNumSet-1] ) {
						my ($chr, $site) = split /\s+/, $chrsite;
						my $sitenew = $site;
						$RIPNumFin ++;
				#		delete $chrsiteclass3HashRF->{$chr}{$sitenew};
						# store the mutated sites
						$UsedChrSiteHash{$chr}{$sitenew} = $altbase;
						$ClusterID ++;
				#		print SITELISTIN join ("\t" => $repi, $round, $Process, $ClusterID, 1, $chrsite), "\n";

							### simple information
						my $selection = $NewStopHash{$chr}{$sitenew}{$altbase} || 'normal';
						my $dup = $HashClassRF->{$chr}{$sitenew};

						my $print1 = join ("\t" => $repi, $round, $Process, "$Binuc-$ClusterID", $RIPNumFin, $chr, $sitenew, $refbase, $altbase, $dirt, $dup, $selection);
						print SITELIST1S "$print1\n";

					}

				}
			}
		}
	############################################################################################################## Secondly Introduce other SNPs

			$Process   = 'InOther';
	#		my $OthNumSet       = int( $RIPRateSet * $UseTotalLen * (1-$RIPRatioSet) + 0.5 );
			my $FoldSet = $OthNumSet / $DirtSetTotalNum;
			my %chrsiteArr2hash;
			for my $chr ( sort keys %GenomeFas ) {
				my $chrlen = length($GenomeFas{$chr});
				for my $site ( 1 .. $chrlen ) {
					my $thisbase = uc( $GenomeEachRF->{$chr}[$site] );
					push @{ $chrsiteArr2hash{$thisbase} }, "$chr $site" if not exists $UsedChrSiteHash{$chr}{$site};
				}
			}

			for my $refbase (keys %chrsiteArr2hash) {
			    $chrsiteArr2hash{$refbase} = [shuffle @{$chrsiteArr2hash{$refbase}}];
			}


			my $OthNumFin  = 0;
			my $Othshufi   = 0;
			my %RefBasei   = qw/A 0 C 0 G 0 T 0/;

			for my $dirt ( sort keys %DirtSet ) {
				my $DirtSetNum = int($DirtSet{$dirt} * $FoldSet + 0.5);
				my ($refbase, $altbase) = split /\>/, $dirt;

				my ($istart, $istop) = ($RefBasei{$refbase}, $RefBasei{$refbase} + $DirtSetNum - 1);
				$RefBasei{$refbase} += $DirtSetNum;

				for my $ii ( $istart .. $istop ) {
					my $chrsite = $chrsiteArr2hash{$refbase}[$ii];
					my ($chr, $sitenew) = split /\s+/, $chrsite;
					
					# store the mutated sites
					$UsedChrSiteHash{$chr}{$sitenew} = $altbase;

					$OthNumFin ++;
					$ClusterID = $OthNumFin;
			#		print SITELISTIN join ("\t" => $repi, $round, $Process, $ClusterID, 1, $chrsite), "\n";

					### simple information
					my $selection = $NewStopHash{$chr}{$sitenew}{$altbase} || 'normal';
#					my $dirt = "$refbase>$altbase";
					my $dup = $HashClassRF->{$chr}{$sitenew};

					my $print1 = join ("\t" => $repi, $round, $Process, $ClusterID, $OthNumFin, $chr, $sitenew, $refbase, $altbase, $dirt, $dup, $selection);
					print SITELIST1S "$print1\n";

				}
			}
###########################################################################
	} # END round
	close SITELIST1S; close SITELISTIN;	close SITEINSWAY;
	}
}

if ( $RunPart =~/Stat/ ) {

################################################################################ Initiate START
	my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime; print STDERR "START $current_time\n";

	my ($UseGeneHashRF, $UseRepeatHashRF, $UseAllHashRF) = &GeneUseInitV8($geneusefile, $repeatusefile);
	my $HashClassRF;

	# Initiate the Genome from the fasta file
	my %GenomeFas = &InitFasV0($InitFasFile, '', 'Hash');

	my ($GenomeEachRF, $HashClassRF, $chrsiteclass2RF, $chrsiteclass3ArrRF, $chrsiteclass3HashRF, $chrsiteArrhashRF ) = &GenomeEachV8(\%GenomeFas, $HashClassRF, 'Stat');
	my ($genecoohashRF, $cdscoohashRF, $cdssitehashRF, $site2genehashRF ) = &InitGeneCOOV8($UseAllHashRF, $cdscoofile, $repeatcoofile);
	my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime; print STDERR "InitGeneCOOV8 Init done at $current_time\n";

	# $chrsiteclass3ArrRF is the array of gene chr sites for shuffle
	# $chrsiteclass3HashRF is the hash of gene chr sites for get seq


	my @UseGeneArr   = sort keys %$UseGeneHashRF;
	my @UseRepeatArr = sort keys %$UseRepeatHashRF;

	my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime;
	print STDERR "Total ", scalar @UseGeneArr, " genes. ", scalar @UseRepeatArr, " Repeats. Init done $current_time\n";
################################################################################ Initiate END

#	@round = qw/0 20 40 60 80 100/;
	for my $repi ( @repi ) {
		my (%MuttedGeneHash, %RIPedGeneHash, %RIPIntroSNPHash, %CTtypehash, %UsedChrSiteRdHash, %StatSNPnumHash);
		my $line  = 0;

		open SITELIST1S, "gzip -dc $folderout/$prefix/${pct}pct.$repi.sitelist.txt.gz | ";	# short information
		while (<SITELIST1S>) { s/\r//g; chomp;
				my ($repi, $round, $Process, $ClusterID, $RIPNumFin, $chr, $sitenew, $refbase, $altbase, $dirt, $dup, $selection) = split /\t/;
#				my ($repiSet, $round, $gid, $countall, $DUP, $fulllen, $sitenewARR, $dirt, $countdirt3ARR, $combnewARR, $flag, $sitecoo, $swaynew, $chr, $sitenew, $base3, $altbase, $dup, $dirt3, $ref5base, $alt5base, $geneARR ) = split /\t/;
				$line ++;	if ($line % 10000000 == 0 ) { my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime; print STDERR "$line lines round $round done $current_time\n"; }

	#			my $selection = $NewStopHash{$chr}{$sitenew}{$altbase} || 'normal';
				if ( $selection eq 'prestop' ) { next; }

				$UsedChrSiteRdHash{$round}{$chr}{$sitenew} = $altbase;
				if ( $Process eq 'InRIP' ) { $StatSNPnumHash{'RIPNumFin'}{$round} ++; } else { $StatSNPnumHash{'OthNumFin'}{$round} ++; }

				# obtain the mutated sites related genes
				my $geneARR = $site2genehashRF->{$chr}{$sitenew};
				# record which gene is altered in this round
				for my $genetmp ( @$geneARR ) {
					if ( $Process eq 'InRIP' ) { $RIPedGeneHash{$round}{$genetmp}{$ClusterID} ++; $RIPIntroSNPHash{$round}{$genetmp} ++; }
					$MuttedGeneHash{$round}{$genetmp} ++ ;
				}

				# Distinct the Ti Tv
				my $CTtype = &CTtypeDetect($dirt);
				$CTtypehash{$round}{$chr}{$sitenew} = $CTtype;
		}
		close SITELIST1S;

		my $current_time = strftime "%Y-%m-%d %H:%M:%S", localtime; print STDERR "All lines done $current_time\n";

################################################################################ Summary of Genes
		for my $SeqClass ( qw/Gene Repeat/ ) {
			open OUTROUNDSTAT, " | gzip > $folderout/$prefix-stat/${pct}pct-$SeqClass-NoprestopSite.$repi.roundstat.txt.gz" or die;
			open OUTROUNDSEQ,  " | gzip > $folderout/$prefix-stat/${pct}pct-$SeqClass-NoprestopSite.$repi.roundseq.txt.gz" or die;
			open OUTCOUNT,     "        > $folderout/$prefix-stat/${pct}pct-$SeqClass-NoprestopSite.$repi.count.txt" or die;

			my $UseThisHashRF; if ( $SeqClass eq 'Gene' ) { $UseThisHashRF = $UseGeneHashRF; } elsif ( $SeqClass eq 'Repeat' ) { $UseThisHashRF = $UseRepeatHashRF; } else { die; }
				my @TiTvSet = qw/10 3 5 7 15/;
				my @Matrix;
				for my $IndexMethod ( qw/RIPindex TiTv10 TiTv3 TiTv5 TiTv7 TiTv15/ ) {
					push @Matrix, "$IndexMethod-TP", "$IndexMethod-FP", "$IndexMethod-TN", "$IndexMethod-FN", "$IndexMethod-TPR", "$IndexMethod-FPR", "$IndexMethod-Recall", "$IndexMethod-Accuracy", "$IndexMethod-Precision";
				}
				print OUTCOUNT join ("\t" => qw/repi round RIPRateSet RIPRatioSet MuttedSeqNum RIPedSeqNum SNPNumOnRIPedSeq/, @Matrix, qw/SNP10RIPedSeqNum SNP1_9RIPedSeqNum SNPNum RIPNumFin OthNumFin SNP10SeqNum SNP1_9SeqNum SNP0SeqNum/ ), "\n";

##################### Detect the Ti and Tv on mutated genes in this round

			for my $round ( @round ) {
				my $RIPRatioSet = $round/100;
				my (%MutCountThisHash, %RIPindexThisHash);

				for my $gene ( sort keys %$UseThisHashRF ) {
					my ($MutCount, $TiTvCount) = qw/0 0/;
					my %MutCountGene = qw/ti 0 tv 0/;
					my ($gchr, $gstrand, $gstart, $gstop);
					my ($geneoriseq, $genenewseq);
					my ($gene2, @coo) = split /\t/, $cdscoohashRF->{$gene};
					for my $gcoo (@coo) {
						($gchr, $gstrand, $gstart, $gstop) = split /\s+/, $gcoo;
				#		$geneoriseq .= substr($GenomeFas{$gchr}, $gstart - 1, $gstop - $gstart + 1);
						for my $gsite ( $gstart .. $gstop ) {
							# summary mutated sites on each gene
							# obtain new seq
							if ( exists $UsedChrSiteRdHash{$round}{$gchr}{$gsite} )  { $MutCount ++; $genenewseq .= $UsedChrSiteRdHash{$round}{$gchr}{$gsite}; }
							else { $genenewseq .= $GenomeEachRF->{$gchr}[$gsite]; }

							# Stat of Ti, Tv on each gene in this round
							if ( exists $CTtypehash{$round}{$gchr}{$gsite} ) { my $CTtype = $CTtypehash{$round}{$gchr}{$gsite}; $MutCountGene{$CTtype} ++; }
						}
					}

					my ($titv, $SNPnum, $RIPtype) = &MutSum($MutCountGene{'ti'}, $MutCountGene{'tv'} );		
			#		print STDERR join ("\t" => $repi, $round, 'titv', $gene, $gene2, $titv, $SNPnum, $RIPtype), "\n";

					# double check
					if ( $SNPnum != $MutCount ) { print STDERR "Mut information not match SNP: $SNPnum   MutCount $MutCount\n"; }


					my ($typereg, $RegLen, $RegCATGcontain, $RegCATGNeedMinGene, $RIPinfoHashRF) = &RIPindexV8SeqInit($genenewseq, 1000, $gene, 'GB2020', $SeqClass );

					$RIPindexThisHash{$gene} = $typereg;
					$MutCountThisHash{'SNPnum'}{$gene} = $MutCount;
				#	print STDERR join ("\t" => $repi, $round, 'geneseq', $gene, $gene2, $typereg, $MutCount, @coo, $genenewseq), "\n";

					$MutCountThisHash{'titv'}{$gene} = $titv;

################ RIP index part START
			#		if ( $round % 10 == 0 ) {
						print OUTROUNDSEQ join ("\t" => $round, $SeqClass, $gene, $typereg, $RIPtype, $titv, $MutCount, $MutCountGene{'ti'}, $MutCountGene{'tv'}, $RIPIntroSNPHash{$round}{$gene} + 0 ), "\n";

						for my $infoi ( sort { $a<=>$b} keys %$RIPinfoHashRF ) {
							my $RF2 = $RIPinfoHashRF->{$infoi};
							print OUTROUNDSTAT join ("\t" => $round, $SeqClass, $gene, $typereg, $infoi, @$RF2 ), "\n";
						}
			#		}
################ RIP index part END
				}
		#		dump %RIPindexThisHash;
#################### Summary on all genes
				my ($RIPedGeneNum, $UseGeneNum, $TotalSNPonGenenum, $TotalSNPonRIPedGenenum) = qw/0 0 0 0/;
				my %countgene = qw/10 0 1-9 0 0 0 RIPed10 0 RIPed1-9 0/;
				my %StatHash;
				for my $TypeSet ( qw/TP FP TN FN/) {
						my $IndexMethod = 'RIPindex';
						$StatHash{$IndexMethod}{$TypeSet} = 0;
					for my $TiTvSet ( @TiTvSet ) {
						my $IndexMethod = "TiTv$TiTvSet";
						for my $TypeSet ( qw/TP FP TN FN/) { $StatHash{$IndexMethod}{$TypeSet} = 0; }
					}
				}

				for my $gene ( sort keys %$UseThisHashRF ) {
					$UseGeneNum ++;
					my $RIPintro;
					if ( exists $RIPedGeneHash{$round}{$gene} ) {
						$RIPintro = 'Y';
						$RIPedGeneNum ++;
						$TotalSNPonRIPedGenenum += $MutCountThisHash{'SNPnum'}{$gene};

						if    ( $MutCountThisHash{'SNPnum'}{$gene} >= 10 ) { $countgene{'RIPed10'} ++; }
						elsif ( $MutCountThisHash{'SNPnum'}{$gene} >=  1 ) { $countgene{'RIPed1-9'}  ++; }
				#		for my $titv ( @titvSet ) { $TitvTPNumHash{$titv} ++ if $MutCountThisHash{'titv'}{$gene} >= $titv; }
					} else {
						$RIPintro = 'N';
					}

					my $Type;
						my $IndexMethod = 'RIPindex';
						if    ( $RIPintro eq 'Y' and $RIPindexThisHash{$gene} eq 'RIPed' ) { $Type = 'TP'; }
						elsif ( $RIPintro eq 'N' and $RIPindexThisHash{$gene} eq 'RIPed' ) { $Type = 'FP'; }
						elsif ( $RIPintro eq 'N' and $RIPindexThisHash{$gene} ne 'RIPed' ) { $Type = 'TN'; }
						elsif ( $RIPintro eq 'Y' and $RIPindexThisHash{$gene} ne 'RIPed' ) { $Type = 'FN'; }
						else { die "$IndexMethod $Type\n"; }
						$StatHash{$IndexMethod}{$Type} ++;

					for my $TiTvSet ( @TiTvSet ) {
						my $IndexMethod = "TiTv$TiTvSet";
						if    ( $RIPintro eq 'Y' and $MutCountThisHash{'titv'}{$gene} >= $TiTvSet ) { $Type = 'TP'; }
						elsif ( $RIPintro eq 'N' and $MutCountThisHash{'titv'}{$gene} >= $TiTvSet ) { $Type = 'FP'; }
						elsif ( $RIPintro eq 'N' and $MutCountThisHash{'titv'}{$gene}  < $TiTvSet ) { $Type = 'TN'; }
						elsif ( $RIPintro eq 'Y' and $MutCountThisHash{'titv'}{$gene}  < $TiTvSet ) { $Type = 'FN'; }
						else { die "$IndexMethod $Type\n"; }
						$StatHash{$IndexMethod}{$Type} ++;
					}

					if    ( $MutCountThisHash{'SNPnum'}{$gene} >= 10 ) { $countgene{'10'} ++; }
					elsif ( $MutCountThisHash{'SNPnum'}{$gene} >=  1 ) { $countgene{'1-9'}  ++; }
					else                                               { $countgene{'0'} ++; }
					$TotalSNPonGenenum += $MutCountThisHash{'SNPnum'}{$gene};

				}
		#		dump %StatHash;

				my $MuttedGeneNum = $countgene{'10'} + $countgene{'1-9'} + 0;

				my @Matrix;
				for my $IndexMethod ( qw/RIPindex TiTv10 TiTv3 TiTv5 TiTv7 TiTv15/ ) { 
					$StatHash{$IndexMethod}{'TPR'}       =   $StatHash{$IndexMethod}{'TP'} / $UseGeneNum;
					$StatHash{$IndexMethod}{'FPR'}       =   $StatHash{$IndexMethod}{'FP'} / $UseGeneNum;
					$StatHash{$IndexMethod}{'Precision'} =   ($StatHash{$IndexMethod}{'TP'} + $StatHash{$IndexMethod}{'FP'}) > 0  ? $StatHash{$IndexMethod}{'TP'} / ( $StatHash{$IndexMethod}{'TP'} + $StatHash{$IndexMethod}{'FP'} ) : $StatHash{$IndexMethod}{'TP'};
					$StatHash{$IndexMethod}{'Recall'}    =   ($StatHash{$IndexMethod}{'TP'} + $StatHash{$IndexMethod}{'FN'}) > 0  ? $StatHash{$IndexMethod}{'TP'} / ( $StatHash{$IndexMethod}{'TP'} + $StatHash{$IndexMethod}{'FN'} ) : $StatHash{$IndexMethod}{'TP'};
					$StatHash{$IndexMethod}{'Accuracy'}  = ( $StatHash{$IndexMethod}{'TP'} +   $StatHash{$IndexMethod}{'TN'} ) / $UseGeneNum ;
					push @Matrix, $StatHash{$IndexMethod}{'TP'}, $StatHash{$IndexMethod}{'FP'}, $StatHash{$IndexMethod}{'TN'}, $StatHash{$IndexMethod}{'FN'}, $StatHash{$IndexMethod}{'TPR'}, $StatHash{$IndexMethod}{'FPR'}, $StatHash{$IndexMethod}{'Recall'}, $StatHash{$IndexMethod}{'Accuracy'}, $StatHash{$IndexMethod}{'Precision'};
				}

				print OUTCOUNT join ("\t" => $repi, $round, $RIPRateSet, $RIPRatioSet, $MuttedGeneNum, $RIPedGeneNum, $TotalSNPonRIPedGenenum, @Matrix, $countgene{'RIPed10'}, $countgene{'RIPed1-9'}, $TotalSNPonGenenum, $StatSNPnumHash{'RIPNumFin'}{$round} + 0, $StatSNPnumHash{'OthNumFin'}{$round} + 0, $countgene{'10'}, $countgene{'1-9'}, $countgene{'0'}), "\n";

			#	if ($round % 1000 == 0 ) { print STDERR "$repi $round done\n"; }
				OUTCOUNT->flush();	OUTROUNDSTAT->flush();	OUTROUNDSEQ->flush();
			} # END each round

			close OUTCOUNT;
			close OUTROUNDSTAT;
			close OUTROUNDSEQ;
			sleep (5);

		}	# END $SeqClass

	} # END each $repi
}


elsif ( $RunPart eq 'Gather' ) {
#	perl Simulation-V1.pl Gather SorMacSN1693 sgl3 '0.1;0.2;0.5;1;2;5;10' 1-100 SimuOutDir

	my @RateSet = split /\;/, $pct;
	my @tmp     = split /\-/, $repiSet;
	my @RepiSet = $tmp[0] .. $tmp[1];

	open OUT1, ">$folderout/$prefix.count.txt";
	for my $RateSet ( @RateSet ) {
		for my $SeqClass ( qw/Gene Repeat/ ) {
			for my $repi ( @RepiSet ) {
				open FILE, "<$folderout/$prefix-stat/${RateSet}pct-$SeqClass-NoprestopSite.$repi.count.txt" or die;
				while (<FILE>) { s/\r//g; chomp;
					if ( /round/ ) {
						print OUT1 join ("\t" => 'SeqClass', $_), "\n" if $repi == 1 and $SeqClass eq 'Gene';
					} else {
						print OUT1 join ("\t" => $SeqClass, $_), "\n";
					}
				}
				close FILE;
			}
		}
	}
	close OUT1;

	open OUT3, ">$folderout/$prefix.TiTvParameters.txt";
	print OUT3 join ("\t" => qw/InModel SeqClass RIPpct TiTvSet Classification GeneNum/), "\n";

	my %Hash;
	my $InModel = $prefix;
	my $pct = 1;
	for my $RepiSet ( @RepiSet ) {
		for my $SeqClass ( qw/Gene Repeat/ ) {
			open FILE, "gzip -dc $folderout/$InModel-stat/${pct}pct-$SeqClass-NoprestopSite.$RepiSet.roundseq.txt.gz | " or die; 
			while (<FILE>) { s/\r//g; chomp;
	#			my ( $repi, $round, $gene, $typereg, $RIPtype, $titv, $MutCount, $ti, $tv, $RIPIntroSNP ) = split /\t/;
				my ( $round, $SeqClass2, $gene, $typereg, $RIPtype, $titv, $MutCount, $ti, $tv, $RIPIntroSNP ) = split /\t/;
				$RIPIntroSNP = $MutCount + 0;
				if    ( $RIPIntroSNP >= 20 ) { $RIPIntroSNP = 20; }
				elsif ( $RIPIntroSNP >= 10 ) { $RIPIntroSNP = 10; }
				elsif ( $RIPIntroSNP >=  1 ) { $RIPIntroSNP =  1; }
				else                         { $RIPIntroSNP =  0; }

				for my $TiTvSet ( qw/3 5 7 10 15/) {
					my $comb;
					if ($round == 0) {
						if    ( $titv >= $TiTvSet ) { $comb = "FP"; }
						elsif ( $titv <  $TiTvSet ) { $comb = "TN"; }
					} elsif ($round == 100 or $round == 80) {
						if    ( $RIPIntroSNP >  0 and $titv >= $TiTvSet ) { $comb = "TP$RIPIntroSNP"; }
						elsif ( $RIPIntroSNP == 0 and $titv >= $TiTvSet ) { $comb = "FP"; }
						elsif ( $RIPIntroSNP == 0 and $titv <  $TiTvSet ) { $comb = "TN"; }
						elsif ( $RIPIntroSNP >  0 and $titv <  $TiTvSet ) { $comb = "FN"; }
					}
					$Hash{$SeqClass}{$round}{$TiTvSet}{$comb} ++;
				}
			}
			close FILE;

		}
	}

		for my $SeqClass ( qw/Gene Repeat/ ) {
			for my $round ( qw/0 80 100/ ) {
				for my $comb ( qw/TP1 TP10 TP20 TN FP FN/) {
					for my $TiTvSet ( qw/3 5 7 10 15/ ) {
						my $comb2 = $comb;
						$comb2=~s/^TP1$/TP SNP:1-9/;
						$comb2=~s/^TP10$/TP SNP:10-19/;
						$comb2=~s/^TP20$/TP SNP:20+/;
						
						my $class = sprintf "RIP%03dpct", ${round};
						$class=~s/^RIP000pct$/RIP 0% All/;
						$class=~s/^RIP080pct$/RIP 80% All/;
						$class=~s/^RIP100pct$/RIP100% All/;
						print OUT3 join ("\t" => $InModel, $SeqClass, $class, $TiTvSet, $comb2, $Hash{$SeqClass}{$round}{$TiTvSet}{$comb} / scalar @RepiSet ), "\n";
					}
				}
			}
		}
}


# Classification of SNPs into transition / transversion
sub MutSum {
# my ($titv, $SNPnum, $RIPtype) = &MutSum($MutCountGene{'ti'}, $MutCountGene{'tv'} );
	my ($TiNum, $TvNum) = @_;
	my $titv;
	if ( $TvNum == 0 ) { $titv = $TiNum; } else { $titv = $TiNum/$TvNum; }
	my $RIPtype = 'deRIP';
	if ( $titv >= 10 ) { $RIPtype = 'RIPed'; }
	my $SNPnum = ($TiNum + $TvNum) || 0;
	return ($titv, $SNPnum, $RIPtype);
}

sub CTtypeDetect {
	my ( $dirt ) = @_;

	my $CTtype;
	if    ($dirt=~/^(C>T)|(G>A)|(T>C)|(A>G)$/) { $CTtype = 'ti'; }
	elsif ($dirt=~/^(A|C|G|T)\>(A|C|G|T)$/i) { $CTtype = 'tv'; }
	else  { $CTtype = 'oth'; }
	return ($CTtype);
}


sub RIPindexV8SeqInit {
# Used only in Stat RunPart
# my ($typereg, $RegLen, $CATGcontain, $CATGNeedMinGene, $RIPinfoHashRF) = &RIPindexV8seq($genenewseq, 1000, $gene);
	my ($SEQ, $lenreg, $tag, $Parameters, $SeqClass) = @_;
	my @RIPSumArr;
	my $step = int($lenreg/2);
	my $minN = 0.05 * $lenreg;
	my $minTA = 0; # 6;
	my $minACGT = 0; # 20;

	my $typereg = 'deRIP';
	my @bestinfo;
	my (%RIPinfoHash, %RIPcooHash);
	my ($CATGNeedMinGene, $CATGcontain, $RegLen );

	for (my $iset = 0; $iset <= length($SEQ); $iset += $step) {
		my $i = $iset;
		my $wini = "win$i";
		my $seq = substr($SEQ, $i, $lenreg);
		my $length = length($seq);
		my $stop = $i + $length;

		if ( $i == 500 and length($SEQ) < $lenreg ) { next; }

		# if the full seq shorter than $lenreg
		if ( length($SEQ) < $lenreg ) {
			# change nothing
		} elsif ( $i == 0 ) {
			# change nothing
		# if the last seq is shorter than $lenreg, re-set the $seq
		} elsif ( $length < $lenreg ) {
			$i = length($SEQ) - $lenreg;
			$seq = substr($SEQ, $i, $lenreg);
			$length = length($seq);
		}

		if ( exists $RIPcooHash{$i}{$stop} ) { next; }
		$RIPcooHash{$i}{$stop} ++;

		my @info = &RIPindexV8RIPper( $wini, $seq, $minACGT, $minN, $i, $stop, $length, $Parameters, $SeqClass );
		$RIPinfoHash{$i} = \@info;
		if ( $i == 0 ) { 
			@bestinfo = @info;
			$CATGNeedMinGene = $info[16];
			$CATGcontain = $info[10] + $info[11];
			$RegLen = $length;
		} else {
			if ( $CATGNeedMinGene > $info[16] ) {
				$CATGNeedMinGene = $info[16];
				$CATGcontain = $info[10] + $info[11];
				$RegLen = $length;
			}
		}

		if ( $info[1] eq 'RIPed' ) { $typereg = 'RIPed'; @bestinfo = @info; }
#		print STDERR join ("\t" => 'indexEach', $tag, $i, $length, $typereg, @info ), "\n";
	}

	return ($typereg, $RegLen, $CATGcontain, $CATGNeedMinGene, \%RIPinfoHash);
}


sub RIPindexV8RIPper {
	my ( $wini, $seq, $minACGT, $minN, $start0, $stop, $length, $Parameters, $SeqClass ) = @_;
	my ($numAT, $numTA, $numCA, $numTG, $numAC, $numGT, $numN, $CATGfold, $Composite);
	my $minN = 0;
	my $minACGT = 0;
	my $typereg = 'deRIP';
	$numN = $seq =~tr/N/N/;
	my ($Product, $Substrate) = qw/na na/;

#	$numAT = ($seq =~ /AT/g) ? () = $seq =~ /AT/g : 0;
	if ($seq=~/AT/) { $numAT = () = $seq=~/AT/g; } else { $numAT = 0;  }
	if ($seq=~/TA/) { $numTA = () = $seq=~/TA/g; } else { $numTA = 0;  }
	if ($seq=~/CA/) { $numCA = () = $seq=~/CA/g; } else { $numCA = 0; }
	if ($seq=~/TG/) { $numTG = () = $seq=~/TG/g; } else { $numTG = 0; }
	if ($seq=~/AC/) { $numAC = () = $seq=~/AC/g; } else { $numAC = 0; }
	if ($seq=~/GT/) { $numGT = () = $seq=~/GT/g; } else { $numGT = 0; }
	$numN = $seq =~tr/N/N/;
	my ($Product, $Substrate) = qw/na na/;

	if ($numN <= $minN) { } else { next; }

	if ( $numAT > 0 )        { $Product    = sprintf "%0.3f", $numTA/$numAT; }                   else { $Product   = sprintf "%0.3f", $numTA; }
	if ( $numAC+$numGT > 0 ) { $Substrate  = sprintf "%0.3f", ($numCA+$numTG)/($numAC+$numGT); } else { $Substrate = sprintf "%0.3f", ($numCA+$numTG); }
	if ( $numTA > 0 )        { $CATGfold   = sprintf "%0.3f", ($numCA+$numTG)/$numTA; }          else { $CATGfold  = sprintf "%0.3f", ($numCA+$numTG); }
	my $Composite = sprintf "%0.3f", ($Product - $Substrate);

	my $TANeed = 1.15 * $numAT - $numTA;
	my $CANeed = $numCA + $numTG - 0.75 * ( $numAC + $numGT );
	my $CATANeedMin;

	if ( $Parameters eq 'Nature2003' ) {
		if    ( $SeqClass eq 'Gene'   ) {  if ( $Product > 2 or $Substrate < 0.7 ) { $typereg = 'RIPed'; } $CATANeedMin = ($TANeed < $CANeed) ? $TANeed : $CANeed; }
		elsif ( $SeqClass eq 'Repeat' ) {  if ( $Product > 2                     ) { $typereg = 'RIPed'; } $CATANeedMin = $TANeed; }
		else { die; }
	} elsif ( $Parameters eq 'GB2020' ) {
		if ( $Product > 2 or $Substrate < 0.7 ) { $typereg = 'RIPed'; }
	} elsif ( $Parameters eq 'PeerjOr2And0' ) {
		if ( ($Product >  1.15  or $Substrate <= 0.75 ) and $Composite >  0 ) { $typereg = 'RIPed'; }
		if ( $TANeed < $CANeed ) { $CATANeedMin = $TANeed; } else { $CATANeedMin = $CANeed; }
	}

#                                  0       1       2        3        4        5       6        7       8       9      10       11     12      13      14         15         16
	my $head2 = join ("\t" => qw/wini   typereg  start0    stop   length    Product   Substrate  Composite  numAT   numTA   numCA   numTG   numAC   numGT   numN   CATGfold    CATANeedMin/);
	return (                    $wini, $typereg, $start0, $stop, $length, $Product, $Substrate, $Composite, $numAT, $numTA, $numCA, $numTG, $numAC, $numGT, $numN, $CATGfold, $CATANeedMin );
}

###############################################################################################


sub InitFasV0 {
	my ($item, $type, $ReturnInfo) = @_;
	$type = 'SingleFile' unless $type;
	$ReturnInfo = 'HashRf' unless $ReturnInfo; 
	my %fas; my $id; 
	my @file; if ($type eq 'ArrayRf') { @file = @$item; } else { push @file, $item; }
	for my $file (@file) {
		if ($file=~/\.gz$/) { open FAS, "gzip -dc $file | " or die "no gziped $file\n"; }
		else                { open FAS, "<$file" or die "no $file\n"; }
		while (<FAS>) { s/\r//g; chomp;
			s/^\s+$//;
			if   (/^>(\S+)/) { $id = $1; }
			else             { $fas{$id} .= $_; }
		}
		close FAS;
	}
	if ($ReturnInfo eq 'Hash') { return %fas; } else { return \%fas; }
}

sub GeneUseInitV8 {
#my ($UseGeneHashRF, $UseRepeatHashRF, $UseAllHashRF) = &GeneUseInitV8($geneusefile, $repeatusefile);
	my ($geneusefile, $repeatusefile) = @_;
	my (%UseGeneHash, %UseRepeatHash, %UseAllHash);
	open MODE, "<$geneusefile" or die;
	while (<MODE>) { s/\r//g; chomp;
		next if /^geneID/;
		my ($gene, @t) = split /\t/;
		$UseGeneHash{$gene} ++;
		$UseAllHash{$gene} ++;
	}
	close MODE;
	open MODE, "<$repeatusefile" or die;
	while (<MODE>) { s/\r//g; chomp;
		next if /^geneID/;
		my ($gene, @t) = split /\t/;
		$UseRepeatHash{$gene} ++;
		$UseAllHash{$gene} ++;
	}
	close MODE;
	return (\%UseGeneHash, \%UseRepeatHash, \%UseAllHash);
}

sub HashClassdupV8 {
	my ($repeatcoofile) = @_;
	my (%HashClass);
	my $UseDupLen = 0;
	open FILE, "<$repeatcoofile" or die;
	while (<FILE>) { s/\r//g; chomp;
		my ($id, $chr, $strand, $start, $stop) = split /\s+/;
		my $len = $stop - $start + 1;
		$UseDupLen += $len;
		for my $site ( $start .. $stop ) {
#			push @chrsiteclass1, "$chr $site";
			$HashClass{$chr}{$site} = 'dup';
		}
	}
	close FILE;
	return (\%HashClass, $UseDupLen);
}

# Initiate each site, including CpA and consider cluster
sub GenomeEachV8 {
	my ($GenomeFasRF, $HashClassRF, $InitPart) = @_;
	my (%GenomeEach, %chrsiteclass3Hash, %BinucArrHash, %BinucCountHash);
	my %chrsiteArrhash;
# @chrsiteclass2 stores nondup chr sites which is not used
#	my (@chrsiteclass2);
	my $UseTotalLen = 0;

	if ( $InitPart eq 'Stat' ) { ### only initiate part function
		for my $chr ( sort keys %$GenomeFasRF ) {
			$UseTotalLen += length( $GenomeFasRF->{$chr});
			my @Each = ( '', split //, $GenomeFasRF->{$chr});
			$GenomeEach{$chr} = \@Each;
		}
	} elsif ( $InitPart eq 'IntroSNP' ) { ### initiate full function
		for my $chr ( sort keys %$GenomeFasRF ) {
			$UseTotalLen += length( $GenomeFasRF->{$chr});
			my @Each = ( '', split //, $GenomeFasRF->{$chr});
			$GenomeEach{$chr} = \@Each;
			for my $site ( 1 .. $#Each + 1 ) {
				if ( not exists $HashClassRF->{$chr}{$site} ) { $HashClassRF->{$chr}{$site} = 'nondup'; }

				# CA / TG
				my ($prebase, $thisbase, $nextbase) = ( lc($Each[$site-1]), uc($Each[$site]), lc($Each[$site+1]) );
				push @{ $chrsiteArrhash{$thisbase} }, "$chr $site";

				my $dup = $HashClassRF->{$chr}{$site};
				if    ( $thisbase eq 'C' ) { push @{ $BinucArrHash{"$thisbase$nextbase $dup"} }, "$chr $site"; $BinucCountHash{"$thisbase$nextbase $dup"} ++; $chrsiteclass3Hash{$chr}{$site} = "$chr $site"; }
				elsif ( $thisbase eq 'G' ) { push @{ $BinucArrHash{"$prebase$thisbase $dup"} }, "$chr $site"; $BinucCountHash{"$prebase$thisbase $dup"} ++; $chrsiteclass3Hash{$chr}{$site} = "$chr $site"; }
			}
		}
	} else { die "sub GenomeEachV8 Fault\n"; }

	return (\%GenomeEach, $HashClassRF, \%BinucArrHash, \%BinucCountHash, \%chrsiteclass3Hash, \%chrsiteArrhash, $UseTotalLen);
}

# Initiate the gene coordinates
sub InitGeneCOOV8 {
# my ($genecoohashRF, $cdscoohashRF, $cdssitehashRF, $site2genehashRF ) = &InitGeneCOOV8($UseAllHashRF, $cdscoofile, $repeatcoofile);
	my ($UseAllHashRF, $cdscoofile, $repeatcoofile) = @_;
	my (%genecoohash, %cdscoohash, %cdssitehash, %site2genehash);

	# Initiate the CDS coordinates
	for my $coofile ( $cdscoofile, $repeatcoofile ) {
		open COO, "<$coofile" or die;
		while (<COO>) { s/\r//g; chomp;
			my ($gene, @gcoo) = split /\t/;
			my ($gchr1, $gstrand1, $gstart1, $gstop1) = split /\s/, $gcoo[0]; 
			my ($gchrx, $gstrandx, $gstartx, $gstopx) = split /\s/, $gcoo[-1]; 

			# only use selected genes in that list
			if ( exists $UseAllHashRF->{$gene} ) { } else { next; }

			# Mark the Intron
	#		for my $gsite ( $gstart1 .. $gstopx ) { $cdssitehash{$gchr}{$gsite} = 'Intron'; }

			# Mark start codon and stop codon
	#		for my $gsite ( $gstart1 .. $gstart1 + 2, $gstopx - 2 .. $gstopx ) { $cdssitehash{$gchr}{$gsite} = 'StartStopCodon'; }

			# Initiate the gene coordinates
			$genecoohash{$gene} = join ("\t" => $gene, $gchr1, $gstrand1, $gstart1, $gstopx);

			for my $gcoo (@gcoo) {
				my ($gchr, $gstrand, $gstart, $gstop) = split /\s/, $gcoo; 
				for my $gsite ( $gstart .. $gstop ) {
					push @{ $site2genehash{$gchr}{$gsite} }, $gene;
	#				$cdssitehash{$gchr}{$gsite} = 'CDS';
	#				$cdssitehash{$gchr}{$gsite} = 'CDS' if $cdssitehash{$gchr}{$gsite} ne 'StartStopCodon';
				}
			}
			$cdscoohash{$gene} = $_;
		}
		close COO;
	}
	return (\%genecoohash, \%cdscoohash, \%cdssitehash, \%site2genehash);
}




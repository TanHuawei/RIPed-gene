#################################################################################################### 04. Gene Annotation

# conda create -n busco_env
# conda activate busco_env
# mamba install busco orthofinder
# BUSCO 5.8.1
# BUSCO v5 database  https://busco-data.ezlab.org/v5/data/lineages/      11-Jan-2024

cd /home/thw/Neurospora-2225/50.GenePredict/
cp -p -a /mnt/san1/usr/thw/info/Neurospora-2225/ref2225v03/ref2225v03.pep.fa .

mkdir /home/thw/Neurospora-2225/50.GenePredict/41.busco
cd /home/thw/Neurospora-2225/50.GenePredict/41.busco/

# genome BUSCO
for database in ascomycota_odb10 sordariomycetes_odb10
do
  for sample in CCN51h1 CCN51h2
  do
    busco -l /home/thw/tansoft7/busco/database-202401/$database -i $sample.genomic.fa -o $sample--geno--$database -m geno --cpu 16 --offline &> $sample--geno--$database.log
  done
done


#
ln -s ../*.pep.fa .

# protein BUSCO
for database in ascomycota_odb10
do
  for sample in ref2225v03 2489 4200 4830 Gtetra Ndiscreta Nhispaniola Ntetra2508A Ntetra2509a Podan ref2225v02 Sbrevicollis Shumana Sordaria trp3
  do
    busco -l /home/thw/tansoft7/busco/database-202401/$database -i $sample.pep.fa -o $sample--prot--$database -m prot --cpu 16 --offline &> $sample--prot--$database.log
  done
done


####### test the protein completeness
mkdir BUSCO_summaries-batch1/
# protein BUSCO
cp -p -a ./*/short_summary.*.txt BUSCO_summaries-batch1/

	#database='eudicots_odb10'
	#for sample in xxx
	#	cp $sample--prot--eudicots_odb10/short_summary.*.txt 
	#done

# draw the plots
generate_plot.py -wd BUSCO_summaries-batch1/

Rscript BUSCO_summaries-batch1/busco_figure.R


####### test the protein completeness
mkdir BUSCO_summaries-batch2/
# protein BUSCO
# cp -p -a ./*/short_summary.*.txt BUSCO_summaries-batch2/

database='ascomycota_odb10'
for sample in ref2225v03 2489 trp3 4830 Gtetra Ndiscreta Nhispaniola Ntetra2508A Sordaria Podan
do
	cp -p -a $sample--prot--$database/short_summary.*.txt  BUSCO_summaries-batch2/
done

# draw the plots
generate_plot.py -wd BUSCO_summaries-batch2/

Rscript BUSCO_summaries-batch2/busco_figure2.R

#################################################################################################### orthofinder



#################################################################################################### Annotation ####################################################################################################


#################################################################################################### prepareDB  coffee Project
# 2024-12-09
## downlaod nr database
# nr.gz 2024-02-07 10:05  186G	0c481846c1c6c1cc807f4a8e2dafb622  nr.gz
curl -C - -L -O  https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
curl -C - -L -O  https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz.md5
md5sum nr.gz

## downlaod uniprot database
curl -C - -L -O  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
curl -C - -L -O  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

## prepare diamond database
cd /home/thw/tandb2/
singularity exec /home/thw/singularity/GENETools202309.sif diamond makedb  --threads 24 --tmpdir ./  --in nr-202402.gz            --db nr-202402        >nr-202402.makedb.log 2>&1
singularity exec /home/thw/singularity/GENETools202309.sif diamond makedb  --threads 24 --tmpdir ./  --in uniprot_sprot.fasta.gz  --db  uniprot_sprot   >uniprot_sprot.makedb.log 2>&1
singularity exec /home/thw/singularity/GENETools202309.sif diamond makedb  --threads 24 --tmpdir ./  --in uniprot_trembl.fasta.gz --db  uniprot_trembl  >uniprot_trembl.makedb.log 2>&1
# diamond v2.1.8.162

## install interrproscan
curl -C - -L -O  https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.72-103.0/interproscan-5.72-103.0-64-bit.tar.gz
curl -C - -L -O  https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.72-103.0/interproscan-5.72-103.0-64-bit.tar.gz.md5

cd /home/thw/tandb2
tar -zxvf interproscan-5.72-103.0-64-bit.tar.gz
cd interproscan-5.72-103.0/
python3.10 setup.py -f interproscan.properties

## prepare for eggnog database
# on centos7v1
mkdir /home/thw/tandb2/eggnog-mapper-data/
singularity exec /home/thw/singularity/emapper-2.1.12.sif download_eggnog_data.py -P -M -y --data_dir /home/thw/tandb2/eggnog-mapper-data
# http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz

#################################################################################################### Annotation 01.check_fasta

cd /home/thw/Neurospora-2225/50.GenePredict/

## checking the genome assembly
awk  '$1 !~ /^>/ &&  $1 !~ /^[A-Z]+$/ {print NR": "$0}' ref2225v03.pep.fa > ref2225v03.out

protfile='/home/thw/Neurospora-2225/50.GenePredict/ref2225v03.pep.fa'

#################################################################################################### Annotation 02.run_emapper

mkdir /home/thw/Neurospora-2225/50.GenePredict/51.emapper/
cd /home/thw/Neurospora-2225/50.GenePredict/51.emapper/
singularity exec /home/thw/singularity/emapper-2.1.12.sif emapper.py -m mmseqs -i $protfile \
     --evalue 1e-4 --output emapper_out --data_dir /home/thw/tandb2/eggnog-mapper-data --cpu 24 >emapper.log 2>&1

#################################################################################################### Annotation 02.run_interpro
# conda activate anno_env
# mamba install -c conda-forge openjdk=11

# remove  * character inside ../ref2225v03.pep.fa or you got failure

mkdir /home/thw/Neurospora-2225/50.GenePredict/52.interpro/
cd /home/thw/Neurospora-2225/50.GenePredict/52.interpro/
/home/thw/tandb2/interproscan-5.72-103.0/interproscan.sh -cpu 24 -i $protfile -dp -goterms -b interproscan_out > interproscan.log 2>&1

#################################################################################################### Annotation 03.run_nr

mkdir /home/thw/Neurospora-2225/50.GenePredict/53.nr/
cd /home/thw/Neurospora-2225/50.GenePredict/53.nr/
singularity exec /home/thw/singularity/GENETools202309.sif diamond blastp --outfmt 6 --evalue 1e-5 --threads 24 --max-target-seqs 1 \
	--db /home/thw/tandb2/nr-202402 --query $protfile --out nr_diamond.tsv 1>nr_diamond.log 2>&1


#################################################################################################### Annotation 04.run_uniprot

mkdir /home/thw/Neurospora-2225/50.GenePredict/54.trembl/
cd /home/thw/Neurospora-2225/50.GenePredict/54.trembl/

singularity exec /home/thw/singularity/GENETools202309.sif diamond blastp --outfmt 6 --evalue 1e-5 --threads 24 --max-target-seqs 1 --db /home/thw/tandb2/uniprot_sprot --query $protfile -o  sprot_diamond.tsv 1>swissprot_diamond.log 2>&1

singularity exec /home/thw/singularity/GENETools202309.sif diamond blastp --outfmt 6 --evalue 1e-5 --threads 24 --max-target-seqs 1 --db /home/thw/tandb2/uniprot_trembl --query $protfile -o  trembl_diamond.tsv 1>trembl_diamond.log 2>&1


#################################################################################################### Annotation 05.statistics
conda activate anno_env
# mamba install r-ggplot2 r-dplyr r-tibble r-sf r-tidyverse r-ggVennDiagram r-argparser

mkdir /home/thw/Neurospora-2225/50.GenePredict/62.stat-gene/
cd /home/thw/Neurospora-2225/50.GenePredict/62.stat-gene/
cp -p -a ../5*/*.tsv .
cp -p -a ../51.emapper/emapper_out.emapper.annotations .

cut -f 1  trembl_diamond.tsv | sort -u > ann.trembl.id
cut -f 1  nr_diamond.tsv     | sort -u > ann.nr.id
cut -f 1  sprot_diamond.tsv  | sort -u > ann.sprot.id
cut -f 1  interproscan_out.tsv   | sort -u > ann.interpro.id

#awk 'NF>1'  user_ko.txt > kegg_out.tsv
#awk 'NF>1 {print $1}'  kegg_out.tsv | sort -u > ann.kegg.id

grep -v "^#" emapper_out.emapper.annotations | cut -f 1| sort -u > ann.emapper.id

ls ann.*.id |awk -F "."  '{print $2"\t"$0}'  > ann_venn.list

ln -s ../ref2225v03.pep.fa protein.sf.fasta
num=`grep -c ">" protein.sf.fasta`

# below is wrong
# Rscript ../script/draw_annVenn.R  ann_venn.list  $num  stat

# emapper kegg & GO
conda install conda-forge::mamba
mamba install bioconda::bioconductor-annotationforge
mamba install bioconda::bioconductor-clusterprofiler
Rscript  ../script/emcp/emapperx.R   emapper_out.emapper.annotations  protein.sf.fasta

# interproscan GO
awk  -F "\t" '$NF ~ /GO:/{print $1"\t"$NF-1}' interproscan_out.tsv > interproscan.GO.tsv
# awk -F "\t" '{for(i=1;i<=NF;i++) if($i ~ /^GO:/) print $1"\t"$i}' interproscan_out.tsv > interproscan.GO.tsv
# awk -F "\t" '{for(i=1;i<=NF;i++) if($i ~ /^GO:/) print $1"\t"$i}' interproscan_out.tsv | sed 's/(InterPro)//g' | sed 's/(PANTHER)//g' | sed 's/|/ /g' > interproscan.GO.tsv
Rscript ../script/GOmapperx.R interproscan.GO.tsv 
 

wc -l *.id
  8841 ann.emapper.id
  8765 ann.interpro.id
  9616 ann.nr.id
  5135 ann.sprot.id
  9428 ann.trembl.id

cat ann.*.id | sort | uniq | wc -l
9664

cat protein.sf.fasta | grep ">" | wc -l
9726

#################################################################################################### Annotation 06.summury
# 2025-03-04

mkdir /home/thw/Neurospora-2225/50.GenePredict/62.stat-gene/
cd /home/thw/Neurospora-2225/50.GenePredict/62.stat-gene/

cut -f1 nr_diamond.tsv | sort | uniq | wc -l
	9616
cut -f1 sprot_diamond.tsv | sort | uniq | wc -l
	5135
cut -f1 trembl_diamond.tsv | sort | uniq | wc -l
	9428
cat *.tsv | grep -vP "^\#" | cut -f1 | sort | uniq | wc -l
	9664



#################################################################################################### ncRNA
###########
# 2025-03-02

########### 00.Rfam_prepare
cd /home/thw/tandb2/

# Rfam v2024.09.16
# first download into centos7v1, upload to iyun32 and unzipped
wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz
wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin
wget http://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/database_files/family.txt.gz

gzip -d family.txt.gz 
gzip -d Rfam.cm.gz

## construct index
singularity exec /home/thw/singularity/ncRNATools202309.sif cmpress Rfam.cm
		Pressed and indexed 4178 CMs and p7 HMM filters (4178 names and 4178 accessions).
		Covariance models and p7 filters pressed into binary file:  Rfam.cm.i1m
		SSI index for binary covariance model file:                 Rfam.cm.i1i
		Optimized p7 filter profiles (MSV part)  pressed into:      Rfam.cm.i1f
		Optimized p7 filter profiles (remainder) pressed into:      Rfam.cm.i1p

mkdir /home/thw/tandb2/Rfam-202409/
mv Rfam.* /home/thw/tandb2/Rfam-202409/

########### 01.ncRNA_shell
# split the genome
mkdir /home/thw/Neurospora-2225/50.GenePredict/55.ncRNA/
cd /home/thw/Neurospora-2225/50.GenePredict/55.ncRNA/

genomefile='../ref2225v03.genomic.fa'
singularity exec /home/thw/singularity/ncRNATools202309.sif tRNAscan-SE -E -o tRNA.results --gff tRNA.gff3 -f tRNA.structures -m tRNA.summary $genomefile 1>tRNAscan-SE.log 2>&1


############################ all barrnap rRNA analysis failed
#conda activate anno_env
#mamba install barrnap
# barrnap v0.9      hmmer v3.4     bedtools v2.31.1
#barrnap --kingdom euk --threads 24 --outseq rRNA.fasta $genomefile 1>rRNA.gff 2>barrnap.log

#I got this issue for a genome that started with a telomere repeat and did not have all four bases in the first few hundred characters. I got around it by replacing the first four characters of each sequence with GATC and then running on the temporary file:
sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' ../CCN51hm.genomic.fa > CCN51hm-changeheadtail.genomic.fa

genomefile=ref2225v03-changeheadtail.genomic.fa
singularity exec /home/thw/singularity/ncRNATools202309.sif barrnap --kingdom euk --threads 24 --outseq rRNA.fasta < $genomefile 1>rRNA.gff 2>barrnap.log

############################ siRNA & miRNA
conda activate anno_env

genomefile='../ref2225v03.genomic.fa'
seqkit split -p 200 -O split $genomefile

ls split/*.fa  | while read fa
do
	size=`seqtk size $fa |awk '{printf $2/1000000}' `
	echo "singularity exec  /home/thw/singularity/ncRNATools202309.sif cmscan --cut_ga --rfam --nohmmonly \
	 --cpu  1 \
	 -Z $size \
	 --tblout $fa.tblout \
	 --fmt 2 \
	 --clanin  /opt/rfam14.9/Rfam.clanin \
	 /opt/rfam14.9/Rfam.cm \
	 $fa > $fa.cmscan 2>$fa.cmscan.log " 
done  > ncRNA_cmd.list

#run_parallel
parallel  -j 20  < ncRNA_cmd.list

########### 03.filter
cat split/*.tblout  > genome_rfam.tblout

grep -v " = "  genome_rfam.tblout > genome_rfam.deoverlapped.tblout

# change to gff3, add family information
perl ../script/tblout2gff3.pl ../script/family.txt genome_rfam.deoverlapped.tblout > genome_rfam.gff3
  
awk -F "\t" '$NF!~ /tRNA/ && $NF!~ /rRNA/ ' genome_rfam.gff3 >  genome_rfam.other.gff3

# obtain miRNA snRNA snoRNA
grep "miRNA" genome_rfam.gff3 > genome_rfam.miRNA.gff3
grep "snoRNA" genome_rfam.gff3 > genome_rfam.snoRNA.gff3
grep "snRNA" genome_rfam.gff3 | grep -v 'snoRNA' > genome_rfam.snRNA.gff3


############################ ncRNA statistics
cat rRNA.gff | grep -vP "^\#" | wc -l
	108

wc -l genome_rfam.*.gff3
   0 genome_rfam.miRNA.gff3
  39 genome_rfam.other.gff3
  17 genome_rfam.snoRNA.gff3
  10 genome_rfam.snRNA.gff3


cat tRNA.summary | grep 'Total tRNAs:'
Total tRNAs:                                446


#################################################################################################### Genome stat ####################################################################################################
# 2025-03-03
mkdir /home/thw/Neurospora-2225/50.GenePredict/61.stat-genome/
cd /home/thw/Neurospora-2225/50.GenePredict/61.stat-genome/

genomefile='../ref2225v03.genomic.fa'

############################# basic
## chop genome into contig
# chmod u+x ../script/breakGCscaf
../script/breakGCscaf -m 1 $genomefile shuffle-contig.fasta

## N50 
singularity exec /home/thw/singularity/Assembly202306.sif assembly-stats $genomefile shuffle-contig.fasta > genome.N50_stats

## CG content
singularity exec /home/thw/singularity/Assembly202306.sif seqkit fx2tab -g -H -l -n -i $genomefile > genome.all_stats

############################# tgs_mapping
cp -p -a /home/thw/Neurospora-2225/21.DNABAM--ref2225v02/09.fq-3Gseq/FGSC2225-Pacbio.0.fq.gz /home/thw/Neurospora-2225/50.GenePredict/FGSC2225-Pacbio.0.fq.gz

genomekeyword=ref2225v03
tgsfile='/home/thw/Neurospora-2225/50.GenePredict/FGSC2225-Pacbio.0.fq.gz'
#tgsfile=CCN51ExcludeCpMtOth.fq.gz
tgskeyword='FGSC2225-Pacbio'

conda activate hifi_env
minimap2 -t 16 -ax map-pb --secondary=no $genomefile $tgsfile | samtools sort - -o genome.aln.bam
samtools faidx     genome.aln.bam > genome.aln.bai
samtools flagstats genome.aln.bam > genome.aln.flagstas
samtools coverage  genome.aln.bam > genome.aln.coverage
awk 'NR > 1 {sum+=$3; rd +=$4; cov += $5; dp+=$7*$5 }END{print "Total: "sum"\nCovbase: "cov"\ncoverage: "cov/sum"\nmeandepth: "dp/sum}' genome.aln.coverage > genome.aln.coverage_all


################################################# Quast evaluate genome assembly

cd /home/thw/Neurospora-2225/50.GenePredict/61.stat-genome/
conda activate hifi_env
mamba install bioconda::quast
# QUAST v5.3.0
mamba install bioconda::gridss
# including Dfam-RepeatMasker.lib.gz

# 2025-03-03

genomekeyword=ref2225v03
quast.py -o ${genomekeyword}.quast ../${genomekeyword}.genomic.fa -t 24 --eukaryote --k-mer-stats --circos --conserved-genes-finding --pacbio ${tgsfile}

############################# 02.Merqury
singularity exec /home/thw/singularity/Assembly202306.sif best_k.sh 41000000
	#	genome: 41000000
	#	tolerable collision rate: 0.001
	#	17.6267


mkdir /home/thw/Neurospora-2225/50.GenePredict/61.stat-genome/merqury-ref2225v03/
cd /home/thw/Neurospora-2225/50.GenePredict/61.stat-genome/merqury-ref2225v03/

conda activate hifi_env
mamba install bioconda::meryl bioconda::merqury
# meryl v1.4.1   merqury  v1.1

## kmer statistics
meryl k=17 memory=256G threads=24 count $tgsfile output $tgskeyword.meryl &> $tgskeyword.meryl.log

## run merquery
genomekeyword=ref2225v03
merqury.sh $tgskeyword.meryl ../../$genomekeyword.genomic.fa $genomekeyword.$tgskeyword.merqury_out
mv completeness.stats $genomekeyword.$tgskeyword.merqury_out.completeness.stats


############################# 04.centromere
# quartet v1.2.0


conda activate quartet_env
mamba install quartet      # this brings the quartet1.0, not use

# install newest v1.2.5
cd /home/thw/tansoft22.i32/
unzip quarTeT-v1.2.5.zip
mv quarTeT-main/ quarTeT-v1.2.5/

# on iyun32

conda activate quartet_env
# mkdir /home/thw/Neurospora-2225/50.GenePredict/61.stat-genome/
cd /home/thw/Neurospora-2225/50.GenePredict/61.stat-genome/
genomekeyword=ref2225v03

python /home/thw/tansoft22.i32/quarTeT-v1.2.5/quartet.py TeloExplorer -i ../$genomekeyword.genomic.fa -c animal -m 10 -p $genomekeyword.Telo
mv $genomekeyword.Telo.telo.* tmp/
mv tmp/ $genomekeyword.Telo/


mkdir $genomekeyword.Centro/
nohup python /home/thw/tansoft22.i32/quarTeT-v1.2.5/quartet.py CentroMiner -i ../$genomekeyword.genomic.fa -p $genomekeyword.Centro -t 8 --TE ../../16.TE/FGSC2225/genome.fa.mod.EDTA.anno/genome.fa.mod.EDTA.TEanno.gff3 --gene ../$genomekeyword.gene.gff > $genomekeyword.Centro.log
mv Candidates/ TandemRepeat/ $genomekeyword.Centro.log $genomekeyword.Centro/

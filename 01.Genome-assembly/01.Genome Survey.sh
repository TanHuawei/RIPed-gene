#################################################################################################### 01.Genome Survey


###################################################################### kmer
cd /home/thw/Neurospora-2225/50.GenePredict/


# FGSC2225-Pacbio.0.fq.gz is the Pacbio sequencing data
ls FGSC2225-Pacbio.0.fq.gz | sed 's/.0.fq.gz//' > 1-sample-se.list

mkdir /home/thw/Neurospora-2225/50.GenePredict/30.kmer/
cd /home/thw/Neurospora-2225/50.GenePredict/30.kmer/

conda activate denovo_env
# mamba install r r-devtools        # r-4.4.2
# install.packages("pracma")
# install.packages("fGarch")
# devtools::install_github("schneebergerlab/findGSE")
# mamba install genomescope2 smudgeplot

# jellyfish v2.3.0


mkdir findGSEhet-k19/  findGSEhet-k21/ genomescope2-k19/ genomescope2-k21/

# kmer=19       # test results indicate kmer=19 is a bad parameter
kmer=21

# For single-ends data, including Pacbio data
for samp in `cat ../1-sample-se.list`
do
if [[ -e ../$samp.0.fq.gz && ! -e $samp.k$kmer.histo && ! -e genomescope2-k$kmer/$samp/summary.txt ]]; then
	echo $samp
	jellyfish count -t 8 -C -m $kmer -s 1G -o $samp.k$kmer.jf -G 2 <(zcat ../$samp.0.fq.gz)
	jellyfish stats -o $samp.k$kmer.stat $samp.k$kmer.jf
	jellyfish histo -v  -t 4 -h 10000000 -o $samp.k$kmer.histo $samp.k$kmer.jf
fi
done


# For paired-ends data
for samp in `cat xx/01.fq/1-sample-pe.list`
do
if [[ -e xx/01.fq/$samp.1.fq.gz && ! -e $samp.k$kmer.histo && ! -e genomescope2-k$kmer/$samp/summary.txt ]]; then
	echo $samp
	jellyfish count -t 8 -C -m $kmer -s 1G -o $samp.k$kmer.jf -G 2 <(zcat xx/01.fq/$samp.1.fq.gz) <(zcat xx/01.fq/$samp.2.fq.gz)
	jellyfish stats -o $samp.k$kmer.stat $samp.k$kmer.jf
	jellyfish histo -v  -t 4 -h 10000000 -o $samp.k$kmer.histo $samp.k$kmer.jf
fi
done


### findGSE
for samp in `ls *.k19.histo | sed 's/.k19.histo//'`
do
if [ ! -e findGSEhet-k$kmer/$samp.R ]; then
	echo "setwd('/home/thw/Neurospora-2225/50.GenePredict/30.kmer')" > findGSEhet-k$kmer/$samp.R
	echo "library(findGSE)" >> findGSEhet-k$kmer/$samp.R
	echo "findGSE(histo=\"$samp.k$kmer.histo\", sizek=$kmer, exp_hom = 140, outdir=\"findGSEhet-k$kmer/$samp\")" >> findGSEhet-k$kmer/$samp.R
	Rscript findGSEhet-k$kmer/$samp.R &> findGSEhet-k$kmer/$samp.log
fi
done


### GenomeScope2  设置 -p 1 为单倍体
for samp in `ls *.k19.histo | sed 's/.k19.histo//'`
do
if [ ! -e genomescope2-k$kmer/$samp/summary.txt ]; then
	genomescope2 -i $samp.k$kmer.histo -o genomescope2-k$kmer/$samp -k $kmer -p 1
fi
done

### GenomeScope2  设置 -p 1 为单倍体
for samp in `ls *.k21.histo | sed 's/.k21.histo//'`
do
if [ ! -e genomescope2-k$kmer/$samp/summary.txt ]; then
	genomescope2 -i $samp.k$kmer.histo -o genomescope2-k$kmer/$samp -k $kmer -p 1
fi
done


### Smudgeplot
# https://github.com/KamilSJaron/smudgeplot

mkdir smudgeplot smudgeplot-k27 smudgeplot-k31

# ref2225v03 is the ID of FGSC2225 genome assembly
species=ref2225v03
kmer=27
kmer=31
for samp in `cat ../1-sample-se.list`
do if [ ! -e smudgeplot-k$kmer/$samp/run.log ]; then
	mkdir smudgeplot-k$kmer/$samp/
	nohup FastK -v -t4 -k$kmer -M64 -T16 ../$samp.0.fq.gz -Nsmudgeplot-k$kmer/$samp/FastK_Table > smudgeplot-k$kmer/$samp/run.log
	rm -f /tmp/$samp.[12].fq
fi
done


ylim=100
for samp in `cat ../1-sample-se.list`
do
  for L in 8 12 16 20
  do
  if [ ! -e smudgeplot-k$kmer/$samp/$samp.L$L.${ylim}_smudgeplot_log10.pdf ]; then
    echo $kmer $samp L$L ${ylim}
    nohup smudgeplot.py hetmers -L $L -t 8 -o smudgeplot-k$kmer/$samp/kmerpairs.L$L --verbose smudgeplot-k$kmer/$samp/FastK_Table >> smudgeplot-k$kmer/$samp/run.log
    nohup smudgeplot.py all -t $samp -o smudgeplot-k$kmer/$samp/$samp.L$L.$ylim smudgeplot-k$kmer/$samp/kmerpairs.L${L}_text.smu -ylim $ylim >> smudgeplot-k$kmer/$samp/run.log
  fi
  done
done


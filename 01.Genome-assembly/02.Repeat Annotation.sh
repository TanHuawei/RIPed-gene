#################################################################################################### 02. Repeat Annotation


######################################## Preparation

mkdir /home/thw/Neurospora-2225/16.TE/
cp -p -a /home/thw/cacao/16.TE/sine_line_base.fa /home/thw/Neurospora-2225/16.TE/


### thwNC12v42-contig20Mt.fas is the genome sequence of FGSC2489, NC12
mkdir /home/thw/Neurospora-2225/16.TE/FGSC2489/
cd    /home/thw/Neurospora-2225/16.TE/FGSC2489/
sed 's/>Supercontig_12./>tig/' /home/thw/Neurospora-2225/16.TE/thwNC12v42-contig20Mt.fas > genome.fa

### ref2225v03.genomic.fa is the genome sequence of FGSC2225
mkdir /home/thw/Neurospora-2225/16.TE/FGSC2225/
cd    /home/thw/Neurospora-2225/16.TE/FGSC2225/
cp -p -a /home/thw/Neurospora-2225/16.TE/ref2225v03.genomic.fa genome.fa


######################################## Let's run EDTA

conda activate EDTA2_env
BuildDatabase -name genome genome.fa

# sine_line_base.fa is the additional database that downloaded from SINEBASE

SINEbase='../sine_line_base.fa'
# EDTA v2.2
	perl /home/thw/tansoft22.i32/EDTA/EDTA.pl --genome genome.fa --species others --sensitive 1 --anno 1 --threads 8 --curatedlib $SINEbase | tee pipe3.1-EDTA.log &

##########
# DeepTE	Based on the Deep learning of known TEs ( Also useful for partial TEs)
# TEsorter	Based on the structure of TEs ( Good for intact TEs )

# conda activate DeepTE_env
# mamba install bioconda::seqtk
# DeepTE version 10/22/2022

#################### DeepTE2 style
conda activate DeepTE_env

mv genome.fa.mod.EDTA.TEanno.sum genome.fa.mod.EDTA.TEanno.sum.bak
mv genome.fa.mod.EDTA.TEanno.gff3 genome.fa.mod.EDTA.TEanno.gff3.bak
mv genome.fa.mod.EDTA.final/genome.fa.mod.EDTA.TElib.fa genome.fa.mod.EDTA.final/genome.fa.mod.EDTA.TElib.fa.bak
mv genome.fa.mod.EDTA.anno/ genome.fa.mod.EDTA.anno-firstrun/

mkdir DeepTE2/
cd DeepTE2/
mkdir DeepTE-LTR/ DeepTE-Other/

TElibfile=../genome.fa.mod.EDTA.TElib.fa
grep "#LTR/unknown" $TElibfile | sed 's/>//' | seqtk subseq $TElibfile - > LTR_unknown.fa
grep "#unknown" $TElibfile | sed 's/>//' | seqtk subseq $TElibfile - > Other_unknown.fa
grep -v "#LTR/unknown" $TElibfile | grep -v "#unknown" | sed 's/>//' | seqtk subseq $TElibfile - > All_known.fa

python /home/thw/tansoft22.i32/DeepTE/DeepTE.py -i LTR_unknown.fa -sp F -m_dir /home/thw/tansoft22.i32/DeepTE_Model/Fungi_model -fam LTR 1>DeepTE-LTR.log 2>&1
mv opt_* DeepTE-LTR.log store_temp_opt_dir DeepTE-LTR/

python /home/thw/tansoft22.i32/DeepTE/DeepTE.py -i Other_unknown.fa -sp F -m_dir /home/thw/tansoft22.i32/DeepTE_Model/Fungi_model 1>DeepTE-Other.log 2>&1
mv opt_* DeepTE-Other.log store_temp_opt_dir DeepTE-Other/

cat DeepTE-LTR/opt_DeepTE.fasta DeepTE-Other/opt_DeepTE.fasta > DeepTE-merge.fasta
cat $TElibfile | grep ">" | sed 's/>//' | awk '{print $1}' > fam-rename.EDTAuse.list
cat DeepTE-merge.fasta | grep ">" | sed 's/>//' | awk '{print $1}' > fam-rename.DeepTEuse.list


perl ../../fam-rename.pl DeepTE-merge.fasta ../../fam-rename.DeepTEuse.new.list > DeepTE-merge.rename.fa


cat All_known.fa DeepTE-merge.rename.fa > All_merged.TElib.fa
cp -p -a All_merged.TElib.fa ../genome.fa.mod.EDTA.final/genome.fa.mod.EDTA.TElib.fa
cp -p -a All_merged.TElib.fa ../genome.fa.mod.EDTA.TElib.fa
cd ../

conda activate EDTA2_env
perl /home/thw/tansoft22.i32/EDTA/EDTA.pl --genome genome.fa --step anno --overwrite 1 --anno 1


################################################################################################## Mclintock pipeline for TEs in re-sequencing data

######## Prepare the database
# Use FGSC2489 NC12 as the reference

cd /home/thw/Data/Neurospora-2225/53.TE/db/

perl -ne 'if (/>/) {s/(\#|\-|\:|\?|\/)/\_/g;} print;' \
    FGSC2489.EDTA2.TElib2.fa \
    > FGSC2489.EDTA2.TElib2.rename.fa

python3 /home/wl/Data/biosoft/mcclintock/mcclintock.py \
    -r FGSC2489.genomic.fa \
    -c FGSC2489.EDTA2.TElib2.rename.fa \
    --make_annotations \
    -o TE_annotation



######## Run the analysis
# ../23.DNAbam7--ref2225v03/01.fq-NC134/0-samp130pe.list is a list contain sample name of resequencing data, one sample in one line
# ../23.DNAbam7--ref2225v03/01.fq-NC134/ dir contain resequencing data 

cat ../23.DNAbam7--ref2225v03/01.fq-NC134/0-samp130pe.list | parallel -j 4 '
  echo "running {}"
  date

  read1=../23.DNAbam7--ref2225v03/01.fq-NC134/{}.1.fq.gz
  read2=../23.DNAbam7--ref2225v03/01.fq-NC134/{}.2.fq.gz

  python3 /home/wl/Data/biosoft/mcclintock/mcclintock.py \
            -r db/FGSC2489.genomic.fa \
            -c db/FGSC2489.EDTA2.TElib2.rename.fa \
            -1 ${read1} -2 ${read2} \
            -p 8 -k general -v sample -n {} --resume \
            -m ngs_te_mapper2,temp2,tebreak \
            -o mcclintock-out/ 2> {}.log



### Chech which samples are processed done
for sample in $(cat ../23.DNAbam7--ref2225v03/01.fq-NC134/0-samp130pe.list); do
    if [ ! -d "mcclintock-out/$sample/tmp" ] && [ -f "mcclintock-out/$sample/results/summary/data/run/summary_report.txt" ]; then
        echo "$sample"
    fi
done




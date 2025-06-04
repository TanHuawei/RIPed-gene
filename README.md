
Example for Simulation
############################
# Simulating introduce of SNPs, Run several times

# 0.1% sequence nucleotide divergence, 3 replicates
seq 1 3 | parallel -j 10 perl Simulation-V1.pl IntroSNP SorMacSN1693 sgl3 0.1 Each{} SimuOutDir

# 0.5% sequence nucleotide divergence, 3 replicates
seq 1 3 | parallel -j 10 perl Simulation-V1.pl IntroSNP SorMacSN1693 sgl3 0.5 Each{} SimuOutDir

# 1% sequence nucleotide divergence, 3 replicates
seq 1 3 | parallel -j 10 perl Simulation-V1.pl IntroSNP SorMacSN1693 sgl3 1 Each{} SimuOutDir

# 2% sequence nucleotide divergence, 3 replicates
seq 1 3 | parallel -j 10 perl Simulation-V1.pl IntroSNP SorMacSN1693 sgl3 2 Each{} SimuOutDir

# 5% sequence nucleotide divergence, 3 replicates
seq 1 3 | parallel -j 10 perl Simulation-V1.pl IntroSNP SorMacSN1693 sgl3 5 Each{} SimuOutDir

# 10% sequence nucleotide divergence, 3 replicates
seq 1 3 | parallel -j 10 perl Simulation-V1.pl IntroSNP SorMacSN1693 sgl3 10 Each{} SimuOutDir

# 15% sequence nucleotide divergence, 3 replicates
seq 1 3 | parallel -j 10 perl Simulation-V1.pl IntroSNP SorMacSN1693 sgl3 15 Each{} SimuOutDir

# 20% sequence nucleotide divergence, 3 replicates
seq 1 3 | parallel -j 10 perl Simulation-V1.pl IntroSNP SorMacSN1693 sgl3 20 Each{} SimuOutDir

############################
# Stats of simulated SNPs
perl Simulation-V1.pl Gather SorMacSN1693 sgl3 '0.1;0.5;1;2;5;10;15' 1-3 SimuOutDir








Please cite "Pervasive influence of repeat-induced point mutations on single-copy functional genes confirm its importance in the evolution of Neurospora crassa".

If you have any questions, feel free to contact TAN, Huawei ( tanhuawei1991@163.com ).

# Thalassiosirales chloroplast phylogenomics

--

Working directory: `/storage/ruck/Thal_orgs/Thal_plastid_genomes_gbf`

### 1. Create fasta files for each protein-coding gene

##### Run the script `extract-gbk-features.py`
```
module load python/3.6.0-anaconda

mkdir fasta-files

python extract-gbk-features.py gbf .

mv *.fna fasta-files/
mv *.faa fasta-files/
```

##### Remove tRNAs, rRNAs, ORFs, and difficult sequences
```
cd fasta-files

mkdir tRNAs
mv trn* tRNAs/
mv tRNA-* tRNAs/

mkdir rRNAs
mv rnl.* rRNAs/
mv rns.* rRNAs/
mv rrl.* rRNAs/
mv rrs.* rRNAs/
mv rrf.* rRNAs/
mv rrn5* rRNAs/
mv 5S* rRNAs/

mkdir orfs
mv ORF* orfs/
mv orf* orfs/

mkdir others
mv acpP1.* others/ # only in outgroups
mv acpP2.* others/ # only in outgroups
mv clpC.* others/ # difficult to align
mv rpoB.* others/ # difficult to align
mv tyrC.* others/ # only in one species
```

##### Remove pseudogenes
| gene | strain |
| --- | --- |
| petF | L549\_Minidiscus_variabilis |
| ycf66 | WR78\_Roundia_cardiophora |

### 2. Create sequence alignments for each protein-coding gene
##### Align using MACSE
```
module load java

for i in *.fna
do
java -jar /home/wader/src/macse_v2.06.jar \
-prog alignSequences \
-gc_def 11 \
-out_AA ${i%.fna}.align.faa \
-out_NT ${i%.fna}.align.fna \
-seq $i
done
```
###### `-gc_def 11` will use the Bacterial, Archael, and Plant Plastid genetic code

##### Remove the final stop codons 
```
for i in *.fna
do
java -jar /home/wader/src/macse_v2.06.jar \
-prog exportAlignment \
-align $i \
-codonForFinalStop --- \
-codonForInternalFS NNN
done

for i in *.align_NT.fna
do
mv $i ${i%.align_NT.fna}.align.fna
mv ${i%.align_NT.fna}.align_AA.fna ${i%_NT.fna}.align.faa
done
```
##### Copy the alignments to different directory
```
cd /storage/wader/organelles-project/alignments-dna/

for i in /storage/ruck/Thal_orgs/Thal_plastid_genomes_gbf/fasta-files/*.align.fna
do
ln -s $i
done

cd ../storage/wader/organelles-project/alignments-amino/

for i in /storage/ruck/Thal_orgs/Thal_plastid_genomes_gbf/fasta-files/*.align.faa
do
ln -s $i
done
```
##### Trim the alignments using ClipKit
```
conda activate clipkit

for i in *.faa
do
clipkit $i -o $i.clipkit \
-m kpic-gappy -g 0.5 \
--log
done

cd ../storage/wader/organelles-project/alignments-dna/

for i in *.fna
do clipkit $i -o $i.clipkit \
-m gappy -g 0.5 \
--log
done

conda deactivate
``` 
###### `-m kpic-gappy` keeps parsimony informative and constant sites
###### `-g 0.5` specifies gap threshold of 50%

##### Remove the gene names from thee fasta headers
```
for i in *.clipkit 
do
cat $i | rev | cut -d "_" -f2- | rev > $i.2
mv $i.2 $i
done

cd ../storage/wader/organelles-project/alignments-amino/

for i in *.clipkit 
do
cat $i | rev | cut -d "_" -f2- | rev > $i.2
mv $i.2 $i
done
```

##### Concatenate the alignments
```
module load python

python /home/wader/src/AMAS/amas/AMAS.py concat \
-i *.clipkit \
-f fasta \
-d aa \
-p plastid.amino.concat.model \
-t plastid.amino.concat.fasta \
-u fasta \
-y raxml

cd ../storage/wader/organelles-project/alignments-dna/

python /home/wader/src/AMAS/amas/AMAS.py concat \
-i *.clipkit \
-f fasta \
-d dna \
-p plastid.codon.concat.model \
-t plastid.dna.concat.fasta \
-u fasta \
-y raxml \
-n 123

python /home/wader/src/AMAS/amas/AMAS.py concat \
-i *.clipkit \
-f fasta \
-d dna \
-p plastid.gene.concat.model \
-t plastid.dna.concat.fasta \
-u fasta \
-y raxml
```

### 3. Phylogenetic analyses
##### AA sequences partitioned by gene
```
cd ../storage/wader/organelles-project/plastid-phylogeny-amino/

conda activate iqtree2

iqtree2 \
-s plastid.amino.concat.fasta \
-p plastid.amino.concat.model \
-m TEST --mset cpREV \
-T AUTO -B 10000 --runs 5 \
--prefix plastid.amino.gene_partitioned
```
##### AA sequences with PartitionFinder to merge gene partitions
```
iqtree2 \
-s plastid.amino.concat.fasta \
-p plastid.amino.concat.model \
--merge rclusterf --rclusterf 10 \
-m TEST --mset cpREV \
-T AUTO -B 10000 --runs 5 \
--prefix plastid.amino.gene_merged
```
##### Model test the AA sequences for the GHOST heterotachy model
```
iqtree2 \
-s plastid.amino.concat.fasta \
-m TESTONLY --mset cpREV \
--mrate H,*H --cmin 2 --cmax 8 \
-T AUTO --prefix plastid.amino.ghost_modeltest
```
| criterion | model selected |
| --- | --- |
| AIC | cpREV+F*H5 |
| AICc | cpREV+F*H5 |
| BIC | cpREV+F+H4|

##### AA sequences under GHOST model with 5 rate classes, unlinked cpREV parameters
```
iqtree2 \
-s plastid.amino.concat.fasta \
-m cpREV+F*H6 \
-T AUTO -B 10000 --runs 5 \
--prefix plastid.amino.ghost_heterotachy
```

##### AA sequences using maximum parsimony
```
mpboot \
-s plastid.amino.concat.fasta \
-st AA -v -seed 12345 -bb 10000 \
-pre plastid.amino.mpboot
```

##### AA sequences using Bayesian inference
```
/home/wader/src/mrbayes-3.2.7/src/mb \
plastid.amino.concat.nexus
```
###### Partition by merged genes from PartitionFinder
###### Model is cpREV+I+G
###### Run for 5,000,000 generations


##### DNA sequences partitioned by gene
```
cd ../storage/wader/organelles-project/plastid-phylogeny-dna/

iqtree2 \
-s plastid.dna.concat.fasta \
-p plastid.gene.concat.model \
-m TEST -T AUTO -B 10000 --runs 5 \
--prefix plastid.dna.gene_partitioned
```
##### DNA sequences with PartitionFinder to merge gene partitions
```
iqtree2 \
-s plastid.dna.concat.fasta \
-p plastid.gene.concat.model \
--merge rclusterf --rclusterf 10 \
-m TEST -T AUTO -B 10000 --runs 5 \
--prefix plastid.dna.gene_merged
```
##### DNA sequences partitioned by gene and codon position, merged with PartitionFinder
```
iqtree2 \
-s plastid.dna.concat.fasta \
-p plastid.codon.concat.model \
--merge rclusterf --rclusterf 10 \
-m TEST -T AUTO -B 10000 --runs 5 \
--prefix plastid.codon.gene_codon_merged
```
##### DNA sequences partitioned by gene and codon position
```
iqtree2 \
-s plastid.dna.concat.fasta \
-p plastid.codon.concat.fasta \
-m TEST -T AUTO -B 10000 --runs 5 \
--prefix plastid.codon.codon_partitioned
```

##### Model test the DNA sequences for the GHOST heterotachy model
```
iqtree2 \
-s plastid.dna.concat.fasta \
-m TESTONLY --mset GTR \
--mrate H,*H --cmin 2 --cmax 8 \
-T AUTO --prefix plastid.dna.ghost_modeltest
```
| criterion | model selected |
| --- | --- |
| AIC | GTR*H9 |
| AICc | GTR*H9 |
| BIC | GTR*H9 |

##### DNA sequences under the GHOST model with 4 rate classes
```
iqtree2 \
-s plastid.dna.concat.fasta \
-m GTR*H9 \
-T AUTO -B 10000 --runs 5 \
--prefix plastid.dna.ghost_heterotachy
```
##### DNA sequences using maximum parsimony
```
mpboot \
-s plastid.codon.concat.fasta \
-st DNA -v -seed 12345 -bb 10000 \
-pre plastid.dna.mpboot
```
##### DNA sequences using Bayesian inference
```
/home/wader/src/mrbayes-3.2.7/src/mb \
plastid.dna.concat.nexus
```
###### Partition by merged genes from PartitionFinder
###### Model is GTR+F+I+G
###### Run for 5,000,000 generations

### 4. Estimate synonymous and nonsynonymous rates
```
cd ..
mkdir dNdS-analysis
cd dNdS-analysis/

for i in ../alignments-dna/*.clipkit
do
ln -s $i
done

ln -s ../plastid-phylogeny-dna/plastid.dna.gene_merged.treefile
```
##### Estimate dN/dS using HyPhy 
###### Standard MG94
###### Global and local models
```
conda activate hyphy

for i in *.clipkit 
do
/home/wader/src/hyphy-analyses/FitMG94/FitMG94.bf \
--alignment $i \
--tree plastid.dna.gene_merged.treefile \
--type global \
--output ${i%.clipkit}.hyphy.global.json
done

for i in *.clipkit
do
/home/wader/src/hyphy-analyses/FitMG94/FitMG94.bf \
--alignment $i \
--tree plastid.dna.gene_merged.treefile \
--type local \
--output ${i%.clipkit}.hyphy.local.json
done
```
##### Summarize the estimated dN/dS for the global model
```

```

### 5. Estimate divergence times
##### Prepare the alignment files, tree files, and control files for MCMCtree.
###### Create a new working directory.
```
cd ..
mkdir divergence-times
cd divergence-times/
```

###### Split the `plastid.dna.concat.fasta` alignment according to the partitions output by PartitionFinder.
```
mkdir alignments
cd alignments

ln -s ../../plastid-phylogeny-dna/plastid.dna.gene_merged.best_scheme
ln -s ../../alignments-dna/plastid.dna.concat.fasta

module load python

python /home/wader/src/AMAS/amas/AMAS.py \
split -i plastid.dna.concat.fasta -f fasta -d dna \
-l plastid.dna.gene_merged.best_scheme -j -u phylip

cd ..
```

###### Open the newick file `plastid.dna.gene_merged.treefile` in a text editor to root the tree, remove branch lengths, and add fossil calibrations to the appropriate branches. Save the edited newick tree as a new file and view in FigTree to make sure the calibrations are placed correctly.
| Node | Calibration |
| --- | --- |
| Root | B(1.31,1.435) |
| *Lithodesmium* | L(0.2996,0.1,1) |
| Thalassiosirales | B(0.75,1.18) |
| *Porosira glacialis* | L(0.09,0.1,1) |
| *Bacterosira* | L(0.0835,0.1,1) |
| *Shionodiscus* | L(0.0615,0.1,2) |
| *Cyclostephanos* | L(0.05,0.1,1) |

###### Combine all the partition phylip alignments into a single file.
```
mkdir prep
cd prep/

for i in ../alignments/*.phy
do
ln -s $i
done

for i in *.phy
do
cat $i >> plastid.dna.gene_merged.concat.phy
done

rm *-out.phy
```

###### Create a control file (`mcmctree_run-prep.ctl`) for MCMCtree with the following lines.
```
seed = -1
seqfile = plastid.dna.gene_merged.concat.phy
treefile = plastid.dna.gene_merged.rooted.tree
mcmcfile = mcmc-prep.txt
outfile = out-prep.txt

ndata = 23
seqtype = 0
usedata = 3
clock = 3
RootAge = <1.0

model = 7
alpha = 0.5
ncatG = 5

cleandata = 0

BDparas = 1 1 0.1
kappa_gamma = 6 2
alpha_gamma = 1 1

rgene_gamma = 2 20 1
sigma2_gamma = 1 10 1

finetune = 1: 0.1 0.1 0.1 0.1 0.1 0.1

print = 1
burnin = 200000
sampfreq = 50
nsample = 20000
```

###### Submit the control file.
```
module load PAML/4.9e

mcmctree mcmctree_run-prep.ctl
```

###### Rename `out.BV` to `in.BV`.
```
mv out.BV in.BV
```

##### Run0 - no data.
###### Copy the alignment, tree file, in.BV, and the control file to a new directory.
```
cd ..
mkdir run0-priors
cd run0-priors

cp ../prep/plastid.dna.gene_merged.concat.phy .
cp ../prep/plastid.dna.gene_merged.rooted.tree .
cp ../prep/in.BV .
cp ../prep/mcmctree_run-prep.ctl mcmctree_run0.ctl
```

###### Edit the following lines in `mcmctree_run0.ctl`. We are running MCMCtree without the data to only estimate the prior distributions.
```
mcmcfile = mcmc-run0.txt
outfile = out-run0.txt
usedata = 0
```

###### Submit `mcmctree_run0.ctl`.
```
mcmctree mcmctree_run0.ctl
```

##### Run1 - autocorrelated rates.
###### Copy the alignment, tree file, in.BV, and the control file to a new directory.
```
cd ..
mkdir run1-autocorr
cd run1-autocorr/

cp ../prep/plastid.dna.gene_merged.concat.phy .
cp ../prep/plastid.dna.gene_merged.rooted.tree .
cp ../prep/in.BV .
cp ../prep/mcmctree_run-prep.ctl mcmctree_run1.ctl
```
###### Edit the following lines in `mcmctree_run1.ctl`. We are running MCMCtree with the molecular data and autocorrelated rates clock model.
```
mcmcfile = mcmc-run1.txt
outfile = out-run1.txt
usedata = 2
clock = 3
```

###### Submit `mcmctree_run1.ctl`.
```
mcmctree mcmctree_run1.ctl
```

##### Run2 - autocorrelated rates.
###### Repeat the steps for Run1 in a new directory called `run2-autocorr`.

##### Run3 - independent rates.
###### Copy the alignment, tree file, in.BV, and the control file to a new directory.
```
cd ..
mkdir run3-indep
cd run3-indep/

cp prep/plastid.dna.gene_merged.concat.phy run3-indep/
cp prep/plastid.dna.gene_merged.rooted.tree run3-indep/
cp prep/in.BV run3-indep/
cp prep/mcmctree_run-prep.ctl run3-indep/mcmctree_run3.ctl

cd run3-indep/
```
###### Edit the following lines in `mcmctree_run3.ctl`. We are running MCMCtree with the molecular data and independent rates clock model.
```
mcmcfile = mcmc-run3.txt
outfile = out-run3.txt
usedata = 2
clock = 2
```

###### Submit `mcmctree_run3.ctl`.
```
mcmctree mcmctree_run3.ctl
```

##### Run4 - independent rates.
###### Repeat the steps for Run3 in a new directory called `run4-indep`.

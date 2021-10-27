# Scripts for inferring the Mitogenome 

The code below consists of the steps I used for the mtDNA phylogenetic analyses. I used Beast v1.10.4 to construct the phylogenetics. 

The fasta file of the two mitogenome datasets I used for final analyses: 

All samples with all positions except Dloop: /home/hennelly/projects/Mitogenome_Aug52020/FINAL_MITOGENOME_DATASET_JAN152021/canids_AWD_August52020_aligned_adjusted_withoutDloop_aligned_edited_FINALUSEDINMANUSCRIPT_reducedsamples.fasta
3rd codon positions to estimate divergence times: /home/hennelly/projects/Mitogenome_Aug52020/codingregions_Jan132021/3rdposition/canids_Jan142021_withoutDloop_3rdposition_combined_final.fasta


The steps are below: 

Step One: Aligning the Mitogenomes

Step Two: Partitioning the mitogenomes into coding and noncoding regions 

Step Three: Extracting the 1st, 2nd, and 3rd codon 
- Step 3.1: Move sequences to Farm
- Step 3.2: Linearize the multifasta file
- Step 3.3: Extract codons and names of sequences 
- Step 3.4: Extract the 1st, 2nd, and 3rd positions 
- Step 3.5: Past the names and sequence files together
- Step 3.6: Cat all files together to merge each of the protein sequence into one multifasta file
- Step 3.7: Paste the sequences of each protein file into one line, determined by identical names of sequences

Step Four: Using Beast for estimating divergence time with the 3rd codon 
- Step 4.1: Make Beauti file 
- Step 4.2: Run BEAST with Beauti file on the farm cluster
- Step 4.3: Move output to desktop to check with Tracer and Figtree


## Step One: Aligning the Mitogenomes 

After creating the fasta file with adding all the individual mitogenomes, these mitogenomes now need to be aligned. I used MUSCLE v3.8.31 muscle3.8.31_i86linux64 to align. 
```
#### 1st step: Move to Farm ####
rsync -avz -e "ssh -p 2022" --progress ~/Desktop/canids_AWD_August52020_withoutDloop_noancientwolves.fasta hennelly@farm.cse.ucdavis.edu:/home/hennelly/projects/Mitogenome_Aug52020/alignedfasta

#### 2nd step: Align with muscle ####
/home/hennelly/bin/muscle3.8.31_i86linux64 -in canids_AWD_August52020_withoutDloop_noancientwolves.fasta -out canids_AWD_August52020_withoutDloop_noancientwolves_aligned.fasta

#### 3rd step: Move back to desktop ####

scp -P 2022 -r hennelly@farm.cse.ucdavis.edu:/home/hennelly/projects/Mitogenome_Aug52020/script/canids_AWD_August52020_aligned_adjusted_withoutDloop_aligned.fasta ~/Desktop
```
After obtaining the aligned fasta file from muscle, I then loaded the fasta file into Sequencher. In Sequencer, I manually checked the alignment and deleted any gaps or "N's that 
were present in the reference fasta (the domestic dog and Ben's African wild dog). 

## Step Two: Partitioning the mitogenomes into coding and noncoding regions 

I used ClustalX 2.1 to create new multi-sample fasta with a specific range of sequence length from the whole mitogenome files, which excluded the Dloop region. For partitioning
the mitogenome, I used the non-coding (tRNA, rRNA) and coding regions of the domestic dog mitogenome from Kim et al. with NCBI reference NC_002008.4.

Kim KS, Lee SE, Jeong HW, Ha JH. 1998. The complete nucleotide sequence of the domestic dog (Canis familiaris) mitochondrial genome. Mol. Phylogenet. Evol. 10(2):210-220

Here are the positions of the non coding and coding regions: 
```
tRNA 1-69 (tRNA1)
rRNA 70-1023 (rRNA1)
tRNA 1024-1090 (tRNA2)
rRNA 1091-2670 (rRNA2)
tRNA 2671-2744 (tRNA3)
CDS (ND1) 2747-3702 (CDS) - ND1 gene
tRNA 3703-3771 (tRNA4)
tRNA 3768-3842 (tRNA5)
tRNA 3844-3913 (tRNA6)
ND2 3914-4955
tRNA 4956-5023 (tRNA7)
tRNA 5037-5105 (tRNA8)
tRNA 5107-5178 (tRNA9)
tRNA 5212-5279 (tRNA10)
tRNA 5280-5347 (tRNA11)
COX1 5349-6893 (COX1)
tRNA 6891-6961 (tRNA12)
tRNA 6966-7033 (tRNA13)
COX2 7034-7717 (COX2)
tRNA 7735-7801 (tRNA14)
ATP8 7803-8006 (ATP8)
ATP6 7964-8644 (ATP6)
COX3 8644-9427 (COX3)
tRNA 9428-9495 (tRNA15)
ND3 9496-9841 (ND3)
tRNA 9842-9910 (tRNA16)
ND4L 9911-10207 (ND4L)
ND4 10201-11578 (ND4)
tRNA 11579-11647 (tRNA17)
tRNA 11648-11707 (tRNA18)
tRNA 11708-11777 (tRNA19)
ND5 11778-13598 (ND5)
ND6 13582-14109 (ND6)
tRNA 14110-14178 (tRNA20)
CYTB 14183-15322 (CYTB)
tRNA 15323-15437 (tRNA21)
```

## Step Three: Extracting the 1st, 2nd, and 3rd codon

For this step three, I extracted the 1st, 2nd, and 3rd codons, as well as the non coding regions within Unix. Its a bit tedious, but the script works well. 

### Step 3.1 

First, I need to move my sequences into farm.
```
rsync -avz -e "ssh -p 2022" --progress ~/Desktop/canids_AWD_August52020_withoutDloop_ND5.fasta hennelly@farm.cse.ucdavis.edu:/home/hennelly/projects/Mitogenome_Aug52020/codingregions/ND5
```

### Step 3.2 

Next, I need to linearlize my multifasta, which I can do the code below in unix.
```
1.) Linearize multiline fasta

cat canids_AWD_August52020_withoutDloop_ND5.fasta | awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0);next;} {printf("%s",$0);}END{printf("\n");}' > canids_AWD_August52020_withoutDloop_ND5_2.fasta

awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' canids_AWD_August52020_withoutDloop_ND5_2.fasta > canids_AWD_August52020_withoutDloop_ND5_linear.fasta

2.) Create each tab as a new line: 
  
sed -e 'y/\t/\n/' canids_AWD_August52020_withoutDloop_ND5_linear.fasta > canids_AWD_August52020_withoutDloop_ND5_linear2.fasta

3.) Remove empty lines: 
  
sed '/^$/d' canids_AWD_August52020_withoutDloop_ND5_linear2.fasta > canids_AWD_August52020_withoutDloop_ND5_linear_Final.fasta
```

### Step 3.3 

In order to just remove specific codons from the sequences, I need to remove and save the names of the sequences into another file. After I extract the codon positions, I can 
then paste back the names of each sequence. The NR % 2 command in unix copy and pastes every other line into a new file. To grab the names of the sequences, I need to add a 
extra line in the original file.
```
#Remove name lines from file and save as a seperate file: 

awk 'NR % 2 == 0' canids_AWD_August52020_withoutDloop_ND5_linear_Final.fasta > canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences.fasta #save sequences in another file

#add a line in the original file to grab the names of each sequence

awk 'NR % 2 == 0' canids_AWD_August52020_withoutDloop_ND5_linear_Final.fasta > canids_AWD_August52020_withoutDloop_ND5_linear_Final_names.fasta #save sequences in another file
```

### Step 3.4 

The script below is used to extract every nth position in the sequences. Because the codons are triplets and coded as every 3rd position, overall I will need to extract one position and then delete the next two positions across all lines. For example, if the sequence was 
GRAYWOLF, I will need to remove RA, WO, F to keep only the 1st letter (G) and then skip two letters and keep the third letter (GYL). For extracting the 2nd and 3rd codon positions, 
I ended up cutting the first character of each line for grabbing 2nd codon positions, and then cut the first two characters to grab the 3rd codon positions. 
```
cut -c 2- < canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences.fasta > canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_nofirstline.fasta # remove first character for extracting 2nd codon positions

cut -c 2- < canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_nofirstline.fasta > canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_nofirstandsecondline.fasta # remove first character again on the output above to have the file start on the 3rd position

### First position ###
sed 's/\(.\{1\}\)\(.\{2\}\)/\1/g' canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences.fasta > canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_1stposition_FINISHED.fasta #keeps the first position of each line and deletes the next two position for each line of the whole file. 
### Second position ###
sed 's/\(.\{1\}\)\(.\{2\}\)/\1/g' canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_nofirstline.fasta > canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_2ndposition_FINISHED.fasta #keeps the first position of each line and deletes the next two position for each line of the whole file. 
### third position ###
sed 's/\(.\{1\}\)\(.\{2\}\)/\1/g' canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_nofirstandsecondline.fasta > canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_3rdposition_FINISHED.fasta #keeps the first position of each line and deletes the next two position for each line of the whole file. 
```

### Step 3.5 

Finally, I need to paste together the names file and sequence files in every other row: 
```
paste -d \\n canids_AWD_August52020_withoutDloop_ND5_linear_Final_names.fasta canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_1stposition_FINISHED.fasta > canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_1stposition_FINISHED_namesandsequences.fasta
paste -d \\n canids_AWD_August52020_withoutDloop_ND5_linear_Final_names.fasta canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_2ndposition_FINISHED.fasta > canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_2ndposition_FINISHED_namesandsequences.fasta
paste -d \\n canids_AWD_August52020_withoutDloop_ND5_linear_Final_names.fasta canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_3rdposition_FINISHED.fasta > canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_3rdposition_FINISHED_namesandsequences.fasta
```

### Step 3.6 

Next, I need to cat all files together to prepare to merge each of the protein sequences into one multifasta file.
```
cat canids_AWD_August52020_withoutDloop_ATP6_linear_Final_sequences_3rdposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ATP8_linear_Final_sequences_3rdposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_COX1_linear_Final_sequences_3rdposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_COX2_linear_Final_sequences_3rdposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_COX3_linear_Final_sequences_3rdposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_CYTB_linear_Final_sequences_3rdposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND1_linear_Final_sequences_3rdposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND2_linear_Final_sequences_3rdposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND3_linear_Final_sequences_3rdposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND4_linear_Final_sequences_3rdposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND4L_linear_Final_sequences_3rdposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_3rdposition_FINISHED_namesandsequences.fasta > canids_AWD_August52020_withoutDloop_3rdposition.fasta 

cat canids_AWD_August52020_withoutDloop_ATP6_linear_Final_sequences_2ndposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ATP8_linear_Final_sequences_2ndposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_COX1_linear_Final_sequences_2ndposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_COX2_linear_Final_sequences_2ndposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_COX3_linear_Final_sequences_2ndposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_CYTB_linear_Final_sequences_2ndposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND1_linear_Final_sequences_2ndposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND2_linear_Final_sequences_2ndposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND3_linear_Final_sequences_2ndposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND4_linear_Final_sequences_2ndposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND4L_linear_Final_sequences_2ndposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_2ndposition_FINISHED_namesandsequences.fasta > canids_AWD_August52020_withoutDloop_2ndposition.fasta

cat canids_AWD_August52020_withoutDloop_ATP6_linear_Final_sequences_1stposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ATP8_linear_Final_sequences_1stposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_COX1_linear_Final_sequences_1stposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_COX2_linear_Final_sequences_1stposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_COX3_linear_Final_sequences_1stposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_CYTB_linear_Final_sequences_1stposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND1_linear_Final_sequences_1stposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND2_linear_Final_sequences_1stposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND3_linear_Final_sequences_1stposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND4_linear_Final_sequences_1stposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND4L_linear_Final_sequences_1stposition_FINISHED_namesandsequences.fasta canids_AWD_August52020_withoutDloop_ND5_linear_Final_sequences_1stposition_FINISHED_namesandsequences.fasta > canids_AWD_August52020_withoutDloop_1stposition.fasta

cat canids_AWD_August52020_withoutDloop_ND6complement.fasta canids_AWD_August52020_withoutDloop_rRNA2.fasta canids_AWD_August52020_withoutDloop_rRNA3.fasta canids_AWD_August52020_withoutDloop_tRNA10.fasta canids_AWD_August52020_withoutDloop_tRNA11.fasta canids_AWD_August52020_withoutDloop_tRNA12.fasta canids_AWD_August52020_withoutDloop_tRNA13.fasta canids_AWD_August52020_withoutDloop_tRNA14.fasta canids_AWD_August52020_withoutDloop_tRNA16.fasta canids_AWD_August52020_withoutDloop_tRNA17.fasta canids_AWD_August52020_withoutDloop_tRNA18.fasta canids_AWD_August52020_withoutDloop_tRNA1.fasta canids_AWD_August52020_withoutDloop_tRNA20.fasta canids_AWD_August52020_withoutDloop_tRNA21.fasta canids_AWD_August52020_withoutDloop_tRNA22.fasta canids_AWD_August52020_withoutDloop_tRNA2.fasta canids_AWD_August52020_withoutDloop_tRNA4.fasta canids_AWD_August52020_withoutDloop_tRNA5.fasta canids_AWD_August52020_withoutDloop_tRNA6.fasta canids_AWD_August52020_withoutDloop_tRNA7.fasta canids_AWD_August52020_withoutDloop_tRNA8.fasta canids_AWD_August52020_withoutDloop_tRNA9.fasta > canids_AWD_August52020_noncoding.fasta
```

### Step 3.7 

Finally, I can paste the sequences of each protein file into one line, determined by having identical names of sequences:
```{r}
perl -nle '/^>/?($n=$_):($s{$n}.=$_);}{print"$_\n$s{$_}"for keys%s' canids_AWD_August52020_withoutDloop_3rdposition.fasta > canids_AWD_August52020_withoutDloop_3rdposition_combined.fasta 

perl -nle '/^>/?($n=$_):($s{$n}.=$_);}{print"$_\n$s{$_}"for keys%s' canids_AWD_August52020_withoutDloop_2ndposition.fasta > canids_AWD_August52020_withoutDloop_2ndposition_combined.fasta 

perl -nle '/^>/?($n=$_):($s{$n}.=$_);}{print"$_\n$s{$_}"for keys%s' canids_AWD_August52020_withoutDloop_1stposition.fasta > canids_AWD_August52020_withoutDloop_1stposition_combined.fasta 

perl -nle '/^>/?($n=$_):($s{$n}.=$_);}{print"$_\n$s{$_}"for keys%s' canids_AWD_August52020_noncoding.fasta > canids_AWD_August52020_noncoding_combined.fasta 

```

## Using Beast for estimating divergence time with the 3rd codon position 

### Step One: make Beauti file 

Move 3rd position to desktop
```
scp -P 2022 -r hennelly@farm.cse.ucdavis.edu:/home/hennelly/projects/Mitogenome_Aug52020/codingregions/3rdposition/canids_AWD_August52020_withoutDloop_3rdposition_combined.fasta ~/Desktop
```

I used the parameters below to generate the BEAUTI file on the Beast version 1 gui within my local computer: 
```
-- "Speciation:Birth-Death Process" for the tree prior, 
-- relaxed lognormal clock, 
-- We used the same partition scheme and substitution models previously determined for gray wolves: HKY+1 for 1st codon and non-coding, TN93+1 for 2nd codon position, and TN93 + G for 3rd codon position (Loog et al. 2020)
-- For the 3rd codon tree, we used a normal prior for the TMRCA of the African wild dog at 3.9 Ma (SD=0.3Ma) based on previous fossil and genetic analyses (Chavez et al. 2019). 
-- Each analylsis was run for 50 million MCMC cycles, sampling every 5,00 and discarding the first 5 million states as burnin
-- Parameter kappa1: logNormal [1,1.25], initial=2
-- Parameter kappa2: logNormal [1,1.25], initial=2
-- Parameter frequencies: Uniform [0,1], initial=0.25
-- Parameter alpha: exponential [0.5], initial=0.5
-- Parameter ucld.stdev: Exponential [0.333333], initial=0.3333333
-- Parameter ucld.mean: LogNormal [-18.42068, 1.5], initial=1
-- Parameter treeModel.rootHeight: Normal [3.7E6,3E5] in [1.9E6, 8E6], initial =4E6
-- Parameter birthDeath.meanGrowthRate: Uniform [0,1E5], initial = 820
-- Parameter birthDeath.relativeDeathRate: Uniform [0,1], initial 0.5
-- Parameter meanRate: Indirectly specified through other Parameters 
-- Parameter covariance: Indirectly specified through other Parameters 
-- Parameter coefficientOfVariation: Indirectly specified through other Parameters
```
After generating the Beauti file, I then transfered the Beauti file to Farm: 
```
rsync -avz -e "ssh -p 2022" --progress ~/Desktop/canids_AWD_3rdposition_10millionMCMC_lognormalprior_TN93_noancientwolves.xml hennelly@farm.cse.ucdavis.edu:/home/hennelly/projects/Mitogenome_Aug52020/BEAST
```

### Step Two: Run BEAST with BEAUTI file on the cluster 

First, I need to create a conda environment and load the program beastie2. Below, I use `conda activate` to load the conda environment containing Beat. 
```
conda activate beastie1:
```
Script: 
```
#!/usr/bin/env bash
#SBATCH --job-name=beast
#SBATCH --time 24:00:00
#SBATCH --mem=15GB
#SBATCH -p high
#SBATCH -o beast.out
#SBATCH -o beast.err
#SBATCH --exclude=c10-96,c10-69,c11-76,c11-93,c11-86

conda init beastie1

beast canids_AWD_3rdposition_10millionMCMC_lognormalprior_TN93_noancientwolves.xml
```

Finally, I can move the results to my desktop 
```
scp -P 2022 -r hennelly@farm.cse.ucdavis.edu:/home/hennelly/projects/Mitogenome_Aug52020/BEAST/10millionMCMC_normalprior_noancientwolves/canids_AWD_August52020_withoutDloop_3rdposition_noancientwolves_combined.trees ~/Desktop
```

### Step Three: Check output on Tracer and Figtree 

With the output on my local computer, I then checked the Beast run with Tracer and Figtree. 






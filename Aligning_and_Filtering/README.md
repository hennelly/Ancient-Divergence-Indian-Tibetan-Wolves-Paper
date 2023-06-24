# Attaching the Y chromosome to the CanFam 3.1 Dog Assembly

Code for attaching the Y chromosome on the Canfam 3.1 Dog Assembly

                                 

## Downloading the Wolf Genome that has the Y chromosome 

I used the de novo assembled gray wolf reference genome that was a male wolf. This was from the reference: 

> S. Gopalakrishnan, J.A.S. Castuita, M.H.S. Sinding, L.F.K. Kuderna, J. Rakkonen, B. Petersen, T. Sicheritz-Ponten, G. Larson, L. Orlando, T. Marques-Bonet, A.J. Hansen, L. Dalen, 
M.T.P. Gilbert. The wolf reference genome sequence (Canis lupus lupus) and its implications for Canis spp. population genomics. BMC Genomics 18: 294 (2017). 

I downloaded the de novo wolf genome from this website:
```
wget https://sid.erda.dk/share_redirect/f1ppDgUPQG/L.Dalen_14_wolf.scf.fasta
wget https://sid.erda.dk/share_redirect/f1ppDgUPQG/L.Dalen_14_wolf.scf.fasta.fai
```                                    

## Attaching the  Y chromosome to Dog Reference Genome                        

For creating a reference genome with the Y chromosome, I need to extract Y chromosome from the wolf genome, add the wolf Y chromosome onto the CanFam3.1 genome, and then index the new canFam3.1 genome.

The canFam3_withY.fa within the directory consists of the CanFam3.1 dog genome and the Y chromosome of the gray wolf genome (Gopalakrishnan et al. 2017)

directory: /home/hennelly/fastqfiles/DogRefwithY/genomes/canFam3_withY.fa

## Step One: Splitting the wolf genome into scaffolds

First, I used this paper to identify the scaffolds that are on the Y chromosome:

https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fmec.15054&file=mec15054-sup-0001-Supinfo.pdf

The wolf genome is full of scaffolds instead of chromosomes. The scaffold of the Y chromosome were assembled based from table S4 in the supplemental materials of Smed et al. 2019 :

```
awk '/^>/{split($1,a,"[|.]")}{print >> a[2]".fa"}' L.Dalen_14_wolf.scf_copy.fasta 
```

## Step Two: Removing PAR from scaffold 242

scaffold242.fa, which turned out to be the largest putatively Y-linked scaffold, consisted of ~425-kb-long region that met the threshold for Y-linkage, followed by a 1.74-Mb-long segment with equal male and female coverage. The 1.74-Mb-long sequence aligned to the dog pseudo-autosomal region (PAR). In the Y chromosome study, they only kept 425-kb MSY region of the scaffold. 
 
 In the supplementary materials, it shows that the PAR region starts at 425kb. First, I'll convert it to a single line file (only letters) and then keep a specific amount of columns in the file (exactly 425,000 columns):
 
```
###Convert scaffold_242 into only letters, which makes the fasta file into a single line
 tr -cd '[:alpha:]' < scaffold_242.fa > scaffold_242_noletters.fa
 ## Keep only characters between 0 to 425000 (corresponding to the Y chromosome region)
cut -c-425000 scaffold_242_noletters.fa > scaffold_242_noletters_fixed.fa
```

## Step Three: Combining all the scaffolds into the Y chromosome

Combining all scaffolds by subsequently "cat"ed together the scaffolds into one file. 

```
cat scaffold_1620.fa scaffold_1839.fa scaffold_2073.fa scaffold_2091.fa scaffold_2320.fa scaffold_2336.fa scaffold_2411.fa scaffold_2415.fa scaffold_2416.fa scaffold_242_noletters_fixed.fa scaffold_2444.fa scaffold_2510.fa scaffold_2542.fa scaffold_2571.fa scaffold_2578.fa scaffold_2580.fa scaffold_2587.fa scaffold_2711.fa scaffold_2714.fa scaffold_2775.fa scaffold_2802.fa scaffold_2915.fa scaffold_2930.fa scaffold_2937.fa scaffold_2955.fa scaffold_2992.fa scaffold_3003.fa scaffold_3040.fa scaffold_3047.fa scaffold_3133.fa scaffold_3144.fa scaffold_3150.fa scaffold_3205.fa scaffold_3210.fa scaffold_3254.fa scaffold_3272.fa scaffold_3276.fa > Wolf_chrY_Feb26.fa

cat Wolf_chrY_Feb26.fa scaffold_3285.fa scaffold_3306.fa scaffold_3311.fa scaffold_3333.fa scaffold_3353.fa scaffold_3462.fa scaffold_3512.fa scaffold_3525.fa scaffold_3536.fa scaffold_3549.fa scaffold_3551.fa scaffold_3558.fa scaffold_3563.fa scaffold_3588.fa scaffold_3616.fa scaffold_3623.fa scaffold_3663.fa scaffold_3693.fa scaffold_3738.fa scaffold_3747.fa scaffold_3834.fa scaffold_3889.fa scaffold_3892.fa scaffold_3911.fa scaffold_3930.fa scaffold_4031.fa scaffold_4051.fa scaffold_4057.fa scaffold_4073.fa scaffold_4146.fa scaffold_4159.fa scaffold_4197.fa scaffold_4251.fa scaffold_4268.fa scaffold_4292.fa scaffold_4310.fa scaffold_4313.fa scaffold_4489.fa scaffold_4675.fa scaffold_4710.fa scaffold_4814.fa scaffold_4936.fa scaffold_4990.fa > Wolf_chrY_2_Feb26.fa

cat Wolf_chrY_2_Feb26.fa scaffold_5031.fa scaffold_5052.fa scaffold_5098.fa scaffold_5267.fa scaffold_5293.fa scaffold_5347.fa scaffold_5497.fa scaffold_5554.fa scaffold_5579.fa scaffold_5584.fa scaffold_5774.fa scaffold_5911.fa scaffold_5929.fa scaffold_5970.fa scaffold_6061.fa scaffold_6087.fa scaffold_6207.fa scaffold_6307.fa scaffold_6405.fa scaffold_6473.fa scaffold_6529.fa scaffold_6535.fa scaffold_6875.fa scaffold_7108.fa scaffold_7290.fa scaffold_7505.fa scaffold_7525.fa scaffold_7553.fa scaffold_7632.fa scaffold_7714.fa scaffold_7725.fa scaffold_7838.fa scaffold_7934.fa scaffold_8081.fa scaffold_8085.fa scaffold_8244.fa scaffold_8316.fa scaffold_8404.fa scaffold_8549.fa scaffold_8598.fa > Wolf_chrY_3_Feb26.fa
```

## Step Four: Adjust the format of the Y chromosome contig

Next, I need to remove the lines with ">" that indicate the different scaffolds, since it should be a continuous file. After removing the ">", I then added a first line written as ">ChrY" using cat. 

```
awk '!/>/' Wolf_chrY_3_Feb26.fa > Wolf_chrY_3_Feb26_fixed.fa
```

Finally, I need to remove all blank spaces from the fasta file. The script below removes everyone other than the letters. 

```
tr -cd '[:alpha:]' < Wolf_chrY_3_Feb26_fixed.fa > Wolf_chrY_3_Feb26_fixed_2.fa
### Add ">ChrY" to the first line
```      

## Step Five: Add the Y chromosome at the end of the dog genome

Add Y chromosome to the end of the dog genome
```
cat /home/hennelly/fastqfiles/DogRef/canFam3.fa /home/hennelly/fastqfiles/wolfrefgenome/Wolf_chrY_3_Feb26_fixed_2.fa > /home/hennelly/fastqfiles/DogRefwithY/canFam3_withY.fa


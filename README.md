# Merging BAM files with Reference Panels

*By Pavel Salazar-Fernandez (epsalazarf@gmail.com), in collaboration with Maria Avila (maricugh@gmail.com) and Andres Moreno-Estrada (amoreno@langebio.cinvestav.mx)*

*Human Population and Evolutionary Genomics Lab | LANGEBIO*

## About
This document explains the procedure to merge data from a sequencing data file (BAM) with the genotyping data from populations of interest. The merged data can be used in a principal components analysis (smartPCA) and individual ancestry estimation (admixture).

>This pipeline was optimized by Maria Avila to filter out possible triallelic sites (most likely cause by sequencing errors) and exclude those with possibility of damage. BAM files processed by this pipeline should have the base qualities recalibrated by mapDamage2, which reduces the quality of mismatches likely caused by damages.

## Preparing the Reference Panel
Before the merging, the Reference Panel(s) must comply to all these requirements:
* **All panels have the same coordinates system:** The coordinates of every SNP are updated with each Human Reference Genome build released. Please check the version used for the manufacture of each genotyping platform. If the build versions utilized differ, the older versions must be updated to match the latest version.
* **All samples have the same ID for each polymorphism:** Especially important when merging data originated from different genotyping platforms. Depending on the manufacturer, some SNP IDs may be private to their platform and thus create problems. If this is condition is not true, the best option is to replace the ID with a compound ID made from the chromosome number and the position of the SNP (Chr_Position, eg. 1_12345). Check the "Renaming SNP IDs" section for the code.
* **No duplicated positions:** Sometimes, more than one rsID or SNP ID are listed at the same position on the same chromosome. When this happens, it is preferable to eliminate these SNPs from further analysis. See the "Finding and Removing Duplicated SNPs" section explaining how to fix this.
* **All SNPs must be located on the Forward Strand:** Since aligned sequencing data is usually provided in a forward strand fashion, the reference panels must also be in the same format, as SNPs that are in the reverse strand will originated false results. If you are not sure if the SNPs provided are in the forward strand format, pick some random rsIDs (skip the ones with alleles G/C or A/T) check their strand info allocated at a SNP database such as [dbSNP](http://www.ncbi.nlm.nih.gov/SNP/).
* **No ambiguous allele orientation:** Those SNPs where the alleles are complimentary to each other (A/T and G/C) are preferably omitted since a correct determination of the orientation of an allele is difficult to achieve. The pipeline itself will ignore sites that follow this behaviour.

If the references panel already comply to all these conditions, you may skip the following sub-sections.

### Finding and removing duplicated SNPs

*NOTE: When modifying a dataset, please make sure to create renamed copies, and leave the original files untouched.*

Use the following terminal instruction to check if more than one SNP ID falls in the same position:

`cut -f 1,4 [FILE.bim] | sort | uniq -d | wc -l`

The script pulls the chromosome number (1) and the position (4) columns from the file, sorts the list and identifies those positions present more than once. If the result is greater than 0, it means that some SNP IDs are repeated two or more times. To generate a list of these duplicated IDs, use the following instructions:

`cut -f 1,4 [FILE.bim] | sort | uniq -d | cut -f 2 > repeatedPos.txt`

`while read line; do grep -wF "$line"  [FILE.bim]; done < repeatedPos.txt | cut -f 2 > repeatedIDs.txt`

The last file created, `repeatedIDs.txt`, contains a list of  IDs that is then feeded to the following code:

`plink --bfile [FILE.bim] --exclude [repeatedIDs.txt] --make-bed --out [FILE2] `

This generates a .bim file with no redundant positions.

### Renaming SNP IDs
To replace the default SNP IDs with the compound name that includes chromosome number and SNP position, use this code:

`awk '{print$1,$1"_"$4,$3,$4,$5,$6}' [FILE.bim] > [FILE.ChrPos.bim]`

Since the modified BIM file needs BED and FAM files named exactly the same, we can create symbolic links renamed to fit the new file with these instructions:

`ln -s [FILE.bed] [FILE.ChrPos.bed]`

`ln -s [FILE.fam] [FILE.ChrPos.fam]`

### Merging reference panels
Now that our reference panels are clean, we can use `plink` to merge them.

`plink --bfile [FILE-1.ChrPos] --bmerge [FILE-2.ChrPos.bin] [FILE-1.ChrPos.bed] [FILE-1.ChrPos.fam] --make-bed --out [OUTPUTNAME] --allow-extra-chr --allow-no-sex`

If the previous steps were performed correctly, the program should not show any error warnings. A succesful merging indicates that the reference panels comply to the previously stated conditions and are ready to use.

### Random allele selection
If the BAM file has a low coverage and genotypes cannot be called without a high degree of uncertainty, it is recommended that the reference panels are forced to homozygosity using a random sampling.

First, we set a seed for the process of randomization:

`RANDOM=0.1`

Then, a random allele for each SNP is selected randomly:

`awk -v seed=$RANDOM 'BEGIN{srand(seed);}  {printf $1" "$2" "$3" "$4" "$5" "$6; for (i=7; i<=NF; i=i+2){x=int(.5 + rand()); printf (" "$(i+x)" "$(i+x));} print "" }' [ReferencePanel.ped] > [ReferencePanel.rand.ped]`

We can copy or soft link the associated MAP file:

`ln -s [ReferencePanel.map] [ReferencePanel.rand.map]`

Finally, we save all the positions:

`awk '{print $1"\t"$4}' [ReferencePanel.bim] > [ReferencePanel.pos]`

## Merging BAM files with the Reference Panel
### Preparations
To do the merging, it is necessary to have installed the following programs, and all the involved files and scripts located or linked in the same directory.

Programs:
* samtools
* plink
* smartpca

Scripts (Provided):
* make_evect_per_sample_mD2_filter.Q20.sh
* overlap_mpileup_map.pl
* overlap_mpileup_tped.pl
* get_random_allele_30.pl

Data (files or symlinks):
* SAMPLE.bam
* REFERENCE_PANEL (.bed,.bim,.fam,.map,.ped,.rand.map,rand.ped,.pos)

### Running
Once all the requirements are met, simply run the main script to perform the analysis:

`make_evect_per_sample_hg19_pavel_mD2_filters.Q20.sh [SAMPLE.bam] [REFERENCE_PANEL]`

Many output files will be created at the running directory, and are identifiable by the tag `MERGED` at the end of the file names.

### Re-Merging merged files
From two merged files, it is possible to produce a third set of merged files that contains all samples from both original files without redundancy:

`plink --file [FILE1] --merge [FILE2.ped] [FILE1.ped] --out [FILE3 name] --allow-no-sex --geno 0.1`

*NOTE*: This command will produce a binary-type file, with no complementary PED and MAP files. Please refer to the following section to learn how to create them.

This merging will produce a **union**, increasing the number of variant sites. If a stricter **intersection** is desired, reduce the `--geno` parameter down to 0; but this could greatly reduce the total sites.

If the merging of more than two datasets is required, repeat the mentioned procedure with the resulting merged datasets and a new one:
> ( fileA + fileB = fileAB; fileAB + fileC = fileABC...)

### Generate PED and MAP files
While some programs can read BED files, others may require PED and MAP files (such as *smartPCA*, used later). To create these files from the binary fileset, run:

`plink --bfile [FILE] --recode --tab --out [output file name]`

### About popinfo
A `popinfo` is a simple text file that contains additional information about the samples, usually provided with the reference panel.

The usual columns, defined by spaces or tabs, are:
* FAMID - *(Required)* Family identification number, utilized for grouping members of a family in association studies. Can be the same as ID.
* ID - *(Required)* Unique sample name, used for matching sample names between different files.
* PID - Parental ID, indicating that this sample is offspring of another male sample in the list. Defaults to 0.
* MID - Maternal ID, same as PID but for a female sample. Defaults to 0.
* SEX - Indicates sex of the sample (0= unknown, 1 = male, 2 = female).
* PHENO - Phenotype, used when doing case-control studies. Defaults to 1.
* POP - *(Required)* Tag for the population of origin, used for grouping together several members from the same population.
* REGION - Geographical region designation that can be used to group many populations together.

Other columns may be added to provide more information (population extended name, linguistic family, geographical coordinates, etc.), but are optional for this procedure.

Merging datasets will cause some samples to be not enlisted in the original popinfo files. Thus, a new popinfo file must be created either adding required rows manually or via copy and paste. All rows must have values in all columns; if no data is available, default the value for that column to either *0* or *NA*. While a popinfo may contain rows for samples not present in the fileset, it is preferrable that only the samples used are enlisted.

##Principal Components Analysis (PCA)
*Documentation: [EIGENSOFT (GitHub)](http://github.com/chrchang/eigensoft/blob/master/POPGEN/README) *

During the merging process, the required files for PCA plotting (.eval and .evec) are automatically produced. But if a composite dataset was created as described before (see **Re-Merging merged files**), it is necessary to identify the principal components again for this dataset.

### Performing a PCA
*NOTE:* If you already have the EVAC and EVAL files for your dataset, you may skip this section.

Program:
* smartPCA

Required files:
* SAMPLE.ped
* SAMPLE.map

First, the user must create a PAR file, that is the input for `smartPCA`. It is a simple text file containing the names for the input and output files. You can use the following text as a template:

```
genotypename:	file.ped  
snpname:         file.map  
indivname:       file.ped  
evecoutname:	 file.evec  
evaloutname:     file.eval  
familynames:     NO  
numoutlieriter:  0
```
Save the file, preferably with a `.par` termination.

To run the analysis, simply type the command:

`smartPCA -p [PAR file]`

The program will generate an EVAL and an EVEC file, both required for plotting.

### Plotting

Program:
* R

Required files:
* Bam+RefPanel_PCA_Plot.R (Provided)
* SAMPLE.eval (Output of smartPCA)
* SAMPLE.evec (Output of smartPCA)
* popinfo.txt


## ADMIXTURE
*Documentation: [ADMIXTURE (UCLA)](https://www.genetics.ucla.edu/software/admixture/)*

### Preparing the input file
Program:
* plink

Required files:
* SAMPLE.bed
* SAMPLE.bim
* SAMPLE.fam

### Linkage Disequilibrium
Due to linkage disequilibrium, not all SNPs are informative. To thin the marker set for linkage disequilibrium, run the following command:

`plink --bfile SAMPLE --indep-pairwise 50 10 0.1`

This command will produce two files:
* plink.prune.in (SNPs targeted for inclusion)
* plink.prune.out (SNPs targeted for exclusion)

To aplly the filter to the dataset and generate a new cleaned dataset, run now the following command:

`plink --noweb --bfile SAMPLE --mind 0.1 --geno 0.05 --extract plink.prune.in --alleleACGT --make-bed --out SAMPLE.clean`

The options used are:
* `--mind 0.1`: filter samples with >10% missing SNPs. Ancient DNA samples may be excluded due to this filter, so this option may be omitted to prevent that.
* `--geno 0.05`: filters SNPs with >5% missing samples.
* `--extract plink.prune.in`: filters SNPs with >0.1 LD in 50-SNP windows.
* `--alleleACGT`: recodes genotypes from numbers to letters. 
* `--out SAMPLE.clean`: The output file name. A name different from the original file is recommended to prevent overwriting (a tag can be added as in the example).

### Running ADMIXTURE

Program:
* admixture

Required files:
* SAMPLE.clean.bed
* SAMPLE.clean.bim
* SAMPLE.clean.fam

To run ADMIXTURE, modify the following command:

`admixture --cv SAMPLE.clean.bed [Number of Ks] > SAMPLE.clean.[Number of Ks].log &`

ADMIXTURE will produce a set of files for each analysis:
* .Q file: ancestry fractions for each individual.
* .P file: allele frequencies in the inferred populations.
* log file: Details on parameters used and other values.

Each analysis for a given number of Ks must be run independently. Recommended values range all numbers between 3 and 10, but depends on the dataset. Keep in mind that higher Ks will require greater computing time.

The `--cv` flag estimates the cross-validation error for each run at the bottom of each log file. The lower the CV error, the more accurate the model involving that *K* number of clusters. To quickly view all CV errors, run:

`grep -h CV *.K*.log`

### Plotting Results

Program:
* R

Required files:
* plot-admixture.R (Provided)
* All `.Q`-type files
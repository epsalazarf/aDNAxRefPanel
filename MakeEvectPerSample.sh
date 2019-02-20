###
#
# MAKE PCA plot based on bam file and reference panel
# November 2012
# maricugh@gmail.com
# Update:
# This pipeline was optimized July 2013 by Maria to 
# only filter out possible triallelic sites (most likely caused by sequence error)
# and ignore the possibility of damage 
# ideally the bam files processed by this pipeline should have the base qualities
# recalibrated by mapDamage2 which reduces the quality of mismatches likely caused by damages
#  
###

#awk '{printf $1" "$2" "$3" "$4" "$5" "$6; for (i=7; i<=NF; i=i+2){x=int(.5 + rand()); printf (" "$(i+x)" "$(i+x));} print "" }' $2.ped > $2.rand.ped
#cp $2.map $2.rand.map

file="$1" 

samtools index $file
base=`basename $file .bam` 

samtools mpileup -Q 0 -f  /data/reference_panels/human_g1k_v37/human_g1k_v37.fasta  -l $2.pos $file |get_random_allele.pl > $base.Q20.rand.mpileup


##Make tfam, tped and map files for sample

overlap_mpileup_tped.pl  $2.rand.map $base.Q20.rand.mpileup > $base.$2.tped 
overlap_mpileup_map.pl $2.rand.map $base.Q20.rand.mpileup  > $base.$2.map 

echo "$base $base 0 0 3 1" > $base.$2.tfam

cut -f2 $base.$2.tped > $base.$2.tped.ids

plink2 --bfile $2 --extract $base.$2.tped.ids --make-bed --out $base.preFilt.sites

 
paste $base.preFilt.sites.bim $base.$2.tped | awk '{if ((($5 ~ /[A]/ || $6 ~ /[A]/) && $11 ~ /[A]/) || ( ($5 ~ /[C]/ || $6 ~ /[C]/) && $11 ~ /C/) || (($5 ~ /[G]/ || $6 ~ /[G]/) && $11 ~ /G/ ) || (($5 ~ /[T]/ || $6 ~ /[T]/) && $11 ~ /[T]/)) print $1,$2,$3,$4}' > $base.filt_sites.map 

##for id in `cut -f2  $base.$2.tped `; do  grep "$id	" $2.bim ; done | paste - $base.$2.tped |awk '{if ((($5 ~ /[A]/ || $6 ~ /[A]/) && $11 ~ /[A]/) || ( ($5 ~ /[C]/ || $6 ~ /[C]/) && $11 ~ /C/) || (($5 ~ /[G]/ || $6 ~ /[G]/) && $11 ~ /G/ ) || (($5 ~ /[T]/ || $6 ~ /[T]/) && $11 ~ /[T]/)) print $1,$2,$3,$4}' > $base.filt_sites.map

## Extract the overlapping positions in the reference panel and in the sample

cut -f2 -d ' ' $base.filt_sites.map > $base.filt_sites.map.ids

plink2 --tfile $base.$2 --extract $base.filt_sites.map.ids --recode --out $base.$2.filt_sites 
plink2 --file $2.rand --extract $base.filt_sites.map.ids --recode --out $base\_extract 

plink2 --file $base.$2.filt_sites --merge $base\_extract.ped $base\_extract.map --recode --out $base.filtSites\_MERGED --allow-no-sex

### BUILD PARFILES and run pca file 


mapfile=$base.filtSites\_MERGED.map
mapbase=`basename $mapfile .map`
sed 's/ -9 / 1 /g' $mapbase.ped > tmp
mv tmp $mapbase.ped 
samplename=`echo $base |cut -f1 -d '.'|cut -f2 -d'_'`
sed 's/'$base'/'$samplename'/g' $mapbase.ped > tmp
mv tmp $mapbase.ped 

echo -e  genotypename: $mapbase.ped"\n"snpname: $mapfile"\n"indivname: $mapbase.ped"\n"evecoutname: $mapbase.evec"\n"evaloutname: $mapbase.eval"\n"numoutlieriter: 0 > $mapbase.par 
#echo -e  genotypename: $mapbase.ped"\n"snpname: $mapfile"\n"indivname: $mapbase.ped"\n"evecoutname: $mapbase.noOut.evec"\n"evaloutname: $mapbase.noOut.eval > $mapbase.noOut.par
echo "smartpca -p $mapbase.par > $mapbase.out" 
smartpca -p $mapbase.par > $mapbase.out
#smartpca -p $mapbase.noOut.par > $mapbase.noOut.out

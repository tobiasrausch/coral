#!/bin/bash

# GRCh38
wget 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz'
bcftools view --min-ac 3 -v snps -m2 -M2 ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz | bcftools annotate -x ^FORMAT/GT -x ^INFO/AC - | grep -v "^##bcftools_" | bcftools view -O b -o GRCh38.snps.bcf -
bcftools index GRCh38.snps.bcf
rm ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz

# GRCh37
wget 'http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz'
bcftools view --min-ac 3 -v snps -m2 -M2 ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz | bcftools annotate -x ^FORMAT/GT -x ^INFO/AC - | grep -v "^##bcftools_" | grep -v "^##ALT=" | bcftools view -O b -o GRCh37.snps.noChr.bcf -
bcftools index GRCh37.snps.noChr.bcf
rm ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz

# Rename chromosomes
rm -f rename.fwd.chrs rename.rev.chrs
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
do
    echo chr${CHR} ${CHR} >> rename.fwd.chrs
    echo ${CHR} chr${CHR} >> rename.rev.chrs
done
bcftools annotate -O b -o GRCh38.snps.noChr.bcf --rename-chrs rename.fwd.chrs GRCh38.snps.bcf
bcftools index GRCh38.snps.noChr.bcf
bcftools annotate -O b -o GRCh37.snps.bcf --rename-chrs rename.rev.chrs GRCh37.snps.noChr.bcf
bcftools index GRCh37.snps.bcf
rm rename.fwd.chrs rename.rev.chrs

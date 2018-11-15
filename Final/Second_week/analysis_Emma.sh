#!/bin/bash

# Pour me connecter sur la machine de Thibault, ou est aussi connectée justine
emma$ ssh -A -X -p 22 ubuntu@bedtool
wget https://tinyurl.com/2018-tp-polymorphism

#download du chr et unzipage
gunzip ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz


#pour pouvoir utiliser les fichiers pythons depuis n'importe où
for PYTHON_SCRIPT in ./src/*.py
do
chmod a+x ${PYTHON_SCRIPT}
done
echo 'export PATH='$(pwd)':${PATH}' >> ~/.bashrc
source ~/.bashrc

#la on a trier le fichier gtf, pour ensuite utiliser bedtools : la mnt on a un fichier propre
gtf_to_bed.py -g Homo_sapiens.GRCh38.94.chr.gtf

#on a modifié la fonction python pour enlever les sequences qui ne sont pas multiples de 3 

#bedtools 
bedtools sort -i Homo_sapiens.GRCh38.94.chr.bed > Homo_sorted.bed
bedtools merge -c 4 -o distinct -i Homo_sorted.bed > Homo_merge.bed
bedtools intersect -header -wb -a ALL.chr20_GRCh38.genotypes.20170504.vcf -b Homo_merge.bed > Chr20_intersect.vcf

 
#Maintenant on veut comparer les deux types de variances
#Comment marche vcf_analysis.py? on a un tableau, pour chaque personnes, pour chaque SNIP, on a une valeur : soit 0 (pas de mutations), 1 (hétérozygotes), 2 (homozygotes)
##ATTENTION dans l'article ils ont décidé de ne prendre en compte que les snp  qui sont présents sur 1 seule personne dans la pop, pas plus !! Mais si ça se trouve c'était une mutation ponctuelle qui n'a aucun impact sur la protéine 
#génotypage =on est capable de dire, l'individu a cet allele la . On s'intéresse à certaines positions seulement
#sequencage = quand on sequence on genotype, mais pas inversement. Quand on sequence on génotype toutes les positions

## A partir d'ici j'ai repris le fichier analysis.sh de Justine car je n'étais pas présente mercredi
# Classification des mutations synonymes, non-synonymes ou stop
# Téléchargement des transcrits
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
gunzip Homo_sapiens.GRCh38.cds.all.fa.gz
vcf_coding_polymorphism.py -f Homo_sapiens.GRCh38.cds.all.fa -g Homo_sapiens.GRCh38.94.chr.gtf -v Chr20_intersect.vcf

# Filtrage par population + analyse de la variance
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
# Fichier qui contient la population pour chaque individu
for CLASS in "Stop" "Syn" "NonSyn"
do
    for POP in "EUR" "AFR" "EAS" "AMR" "SAS" "GBR" "FIN" "CHS" "PUR" "CDX" "JPT" "CLM" "IBS" "PEL" "PJL" "KHV" "LWK" "ACB" "GWD" "ESN" "BEB" "MSL" "MXL" "STU" "ITU" "CEU" "YRI" "CHB" "ASW" "TSI" "GIH"
    do
    extract_pop.py -p integrated_call_samples_v3.20130502.ALL.panel -k ${POP}
    vcftools --vcf Chr20_intersect.${CLASS}.vcf --keep ${POP}.txt --recode --out ${CLASS}_${POP}
    done
done

for POP in "EUR" "AFR" "EAS" "AMR" "SAS" "GBR" "FIN" "CHS" "PUR" "CDX" "JPT" "CLM" "IBS" "PEL" "PJL" "KHV" "LWK" "ACB" "GWD" "ESN" "BEB" "MSL" "MXL" "STU" "ITU" "CEU" "YRI" "CHB" "ASW" "TSI" "GIH"
    do
    vcf_meta_analysis.py -o Stop_${POP}.recode.vcf -s Syn_${POP}.recode.vcf -n NonSyn_${POP}.recode.vcf -c 200
done

# Au final, on obtient des graphes indiquant la distribution de la variance pour les mutations synonymes et non-synonymes, ainsi que le ratio des variances pour les mutations LoF. On peut ainsi comparer ces ratios (LoF vs synonymes) afin de savoir à quel type d'interaction sont soumises les LoF.




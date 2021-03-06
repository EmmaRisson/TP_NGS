# MA VERSION
#hey c'est un test pour voir si ca passe via git
#2eme motif
#!/bin/bash
mkdir -p /mnt/data/variant_calling
cd /mnt/data/variant_calling

########################################################################################################################
# Requirements:
#	Java (version 8)
#	FastQC (version 0.11.7)
#	BWA-MEM (version 0.7.17-r1194-dirty)
#	SAMtools (version 1.9)
#	IGV (version 2.4.14)
#	GATK (version 3.3)
########################################################################################################################

java -version
fastqc -version
bwa
samtools
java -jar ${GATK} --help   # pour la détection de variants
java -jar ${PICARD}

##########################################################
## Download, extract and index the reference chromosome ##
##########################################################

# Download the reference Human chromosome (chromosome 20) from Ensembl
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: compressed reference sequence (.fa.gz)
wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz -O Homo_sapiens.Chr20.fa.gz

# Extract the reference chromosome
# Command: gunzip
# Input: compressed reference sequence (.fa.gz)
# Ouput: reference sequence (remove .gz file) ## Donc ici on décompresse
gunzip Homo_sapiens.Chr20.fa.gz

# Index the reference chromosome # c’est comme un sommaire
# Command: bwa index
# Input: reference (.fa)
# Ouput: indexed reference (.fa.amb, .fa.ann, .fa.bwt, fa.pac, .fa.sa)
bwa index Homo_sapiens.Chr20.fa

######################################################
## Mapping of a family trio to the reference genome ##
######################################################

# The sequences are from an East Asian (Kinh Vietnamese) family forming a trio : daughter/mother/father
# Data available at http://www.internationalgenome.org/data-portal/sample/HG02024
# Daughter:
#       StudyId: SRP004063
#       SampleName: HG02024
#       Library: Pond-206419
#       ExperimentID: SRX001595
#       RunId: SRR822251
#       PlatformUnit: C1E0PACXX121221.6.tagged_373
#       InstrumentModel: Illumina HiSeq 2000
#       InsertSize: 160
# Mother:
#       StudyId: SRP004063
#       SampleName: HG02025
# Father:
#       SampleName: HG02026

#############################
## Mapping of the daughter ##
#############################

# Download paired sequencing reads for the daughter (SampleName: HG02024, RunId: SRR822251)
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: compressed sequencing reads (.fastq.gz)
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02024/sequence_read/SRR822251_1.filt.fastq.gz
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


# Map the paired sequencing reads against the reference Human chromosome 20
# Command: bwa mem
# Options: -M (Mark shorter split hits as secondary for GATK compatibility)
#          -t [number of CPU] (multi-threading)
# Input: indexed reference (.fa), and compressed sequencing reads (.fastq.gz)
# Ouput: alignment (.sam)
bwa mem -t 4 -M Homo_sapiens.Chr20.fa SRR822251_1.filt.fastq.gz SRR822251_2.filt.fastq.gz > SRR822251.sam

# (Optional)
# Compute summary statistics of the alignment
# Command: samtools flagstats
# Input: alignment (.sam)
# Ouput: text file (human and computer readable)
samtools flagstat SRR822251.sam > SRR822251.sam.flagstats

# Compress the alignment and filter unaligned reads
# Command: samtools view
# Options: -@ [number of CPU] (multi-threading)
#	   -S (input format is auto-detected)
# 	   -b (output BAM)
#	   -h (include header in output)
#          -f [flag] (include reads with all  of the FLAGs in INT present) # pour filtrer pour ne regarder que les pairs 
# 	        flag=3 for read paired & read mapped in proper pair, see
#	        https://broadinstitute.github.io/picard/explain-flags.html
# Input: alignment (.sam)
# Ouput: compressed alignment (.bam)
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx > SRR822251.bam
samtools view -@ 4 -Sbh  -f 3  SRR822251.sam> SRR822251.bam # Samtools c’est un programme qui fait des trucs pour nous



# Sort the alignment
# Command: samtools sort
# Input: compressed alignment (.bam)
# Ouput: sorted and compressed alignment (.bam)
samtools sort -@ 4 SRR822251.bam > SRR822251.sorted.bam

# Add Read group (cf https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups) # ça sert à stocker les infos, pour qu’ensuite le programme puisse comprendre d’ou viennent les biais
# Command: gatk AddOrReplaceReadGroups
# Input: alignment (.bam) and read group (Read group identifier, DNA preparation library identifier, Platform, Platform Unit, Sample)
# Ouput: annotated alignment (.bam)
java -jar ${PICARD} AddOrReplaceReadGroups I=SRR822251.sorted.bam \
                                         O=daughter.bam \
                                         RGID=SRR822251 RGLB=Pond-206419 RGPL=illumina \
                                         RGPU=C1E0PACXX121221.6.tagged_373 RGSM=HG02024 RGPI=160

# (Optional)
# Compute statistics of the alignment
# Command: samtools-stats
# Input: alignment (.bam)
# Ouput: text file (human and computer readable)
samtools stats daughter.bam > daughter.bam.stats

# (Optional)
# Plot statistics of the alignment
# Command: plot-bamstats
# Input: statistics text file (output of samtools-stats)
# Ouput: plots (.png)
plot-bamstats -p ~/TP-mardi/plots/ daughter.bam.stats

# Index the alignment
# Command: samtools index
# Input: alignment (.bam)
# Ouput: indexed alignment (.bam.bai)
samtools index daughter.bam


###########################
## Mapping of the mother ##
# Mother:
#       StudyId: SRP004063
#       SampleName: HG02025
# Father:
#       SampleName: HG02026
###########################
# on a toutes les infos 
# On a RunID = SRR361100

head -n 1 20130502.phase3.analysis.sequence.index > mere_info_importante # la on prend que la premiere ligne
grep "SRR361100" 20130502.phase3.analysis.sequence.index >> mere_info_importante  # la je prends que les lignes du gros dossier qui m'intéresse


# Variables definition
FTP_SEQ_FOLDER=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3 # Ftp folder from 1000Genomes project
RUN_ID=SRR361100 # 
SAMPLE_NAME=HG02025 # Sample
INSTRUMENT_PLATFORM=ILLUMINA # Platform/technology used to produce the read
LIBRARY_NAME=Catch-88584 # DNA preparation library identifier
RUN_NAME=BI.PE.110902_SL-HBC_0182_AFCD046MACXX.2.tagged_851.srf # Platform Unit
INSERT_SIZE=96 # Insert size

# Download paired sequencing reads for the mother
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: compressed sequencing reads (.fastq.gz)
wget ${FTP_SEQ_FOLDER}/data/${SAMPLE_NAME}/sequence_read/${RUN_ID}_1.filt.fastq.gz
wget ${FTP_SEQ_FOLDER}/data/${SAMPLE_NAME}/sequence_read/${RUN_ID}_2.filt.fastq.gz


## JE VAIS ESSAYER D'ADAPTER POUR LA MERE MAIS EN AUTOMATISANT
bwa mem -t 4 -M Homo_sapiens.Chr20.fa ${RUN_ID}_1.filt.fastq.gz ${RUN_ID}_2.filt.fastq.gz | samtools view -@ 4 -Sbh  -f 3 | samtools sort -@ 4 > ${RUN_ID}.sorted.bam
# on a enlevé les chevrons pour éviter de créer les dossiers intermédiaires à chaque fois. On crée que celui à la fin
#samtools view ca compresse



# Add Read group
# Command: gatk AddOrReplaceReadGroups
# Input: alignment (.bam) and read group
# Ouput: alignment (.bam)
java -jar ${PICARD} AddOrReplaceReadGroups I=${RUN_ID}.sorted.bam O=mother.bam \
                                         RGID=${RUN_ID} RGLB=${LIBRARY_NAME} RGPL=${INSTRUMENT_PLATFORM} \
                                         RGPU=${RUN_NAME} RGSM=${SAMPLE_NAME} RGPI=${INSERT_SIZE}

# Index the alignment
# Command: samtools index
# Input: alignment (.bam)
# Ouput: indexed alignment (.bam.bai)
samtools index mother.bam

###########################
## Mapping of the father ##
###########################

# Variables definition
SAMPLE_NAME=HG02026 # Sample

# Download index file containing sequencing runs information
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: text file (.index)
wget ${FTP_SEQ_FOLDER}/20130502.phase3.analysis.sequence.index -O 20130502.phase3.index

# Filter paired exome sequencing runs related to father (HG02026)
# Command: grep && grep -v
# Input: tab-separated values file (.index)
# Ouput: filtered comma-separated values file (.index)
grep ${SAMPLE_NAME} 20130502.phase3.analysis.sequence.index | grep "exome" | grep 'PAIRED' | grep -v 'Catch-88526' | grep -v 'Solexa' | grep -v 'from blood' | grep -v '_1.filt.fastq.gz' | grep -v '_2.filt.fastq.gz' | sed 's/\t/,/g' > father.index

# File containing the list of alignments (each line is a .bam file)
# This file is necessary to merge multiple alignments into a single alignment.
# Command: touch
# Input: file name
# Ouput: empty file (.bamlist)
touch father.bamlist

# for each sequencing run (the first 10), align to the reference, sort, add read group and index
head -6 father.index | while IFS="," read FASTQ_FILE MD5 RUN_ID STUDY_ID STUDY_NAME CENTER_NAME SUBMISSION_ID SUBMISSION_DATE SAMPLE_ID SAMPLE_NAME POPULATION EXPERIMENT_ID INSTRUMENT_PLATFORM INSTRUMENT_MODEL LIBRARY_NAME RUN_NAME RUN_BLOCK_NAME INSERT_SIZE LIBRARY_LAYOUT PAIRED_FASTQ WITHDRAWN WITHDRAWN_DATE COMMENT READ_COUNT BASE_COUNT ANALYSIS_GROUP
do

    # Variables definition
    FASTQ_FILE_1=${FASTQ_FILE/.filt.fastq.gz/_1.filt.fastq.gz} # Path of the fasta file in the FTP folder
    FASTQ_FILE_2=${FASTQ_FILE/.filt.fastq.gz/_2.filt.fastq.gz} # Path of the fasta file in the FTP folder (pairing file)

    # Download paired sequencing reads for the father
    # Command: wget
    # Input: url (http:// or ftp://)
    # Ouput: compressed sequencing reads (.fastq.gz)
    wget ${FTP_SEQ_FOLDER}/${FASTQ_FILE_1} -O ${SAMPLE_NAME}_${RUN_ID}_1.filt.fastq.gz
    wget ${FTP_SEQ_FOLDER}/${FASTQ_FILE_2} -O ${SAMPLE_NAME}_${RUN_ID}_2.filt.fastq.gz

    # Map, filter, and sort the paired reads of the sequencing run against the reference genome
    # Command: bwa mem && samtools view && samtools sort
    # Input: indexed reference (.fa), and compressed sequencing reads (.fastq.gz)
    # Ouput: sorted alignment (.bam)
    bwa mem -t 4 -M ../Homo_sapiens.Chr20.fa ${SAMPLE_NAME}_${RUN_ID}_1.filt.fastq.gz ${SAMPLE_NAME}_${RUN_ID}_2.filt.fastq.gz | samtools view -@ 4 -Sbh  -f 3 | samtools sort -@ 4 > ${SAMPLE_NAME}_${RUN_ID}.sorted.bam
    
# Add Read group
    # Command: gatk AddOrReplaceReadGroups
    # Input: alignment (.bam) and read group
    # Ouput: alignment (.bam)
    java -jar ${PICARD} AddOrReplaceReadGroups I=${SAMPLE_NAME}_${RUN_ID}.sorted.bam O=${SAMPLE_NAME}_${RUN_ID}.sorted.RG.bam \
                                         RGID=${RUN_ID} RGLB=${LIBRARY_NAME} RGPL=${INSTRUMENT_PLATFORM} \
                                         RGPU=${RUN_NAME} RGSM=${SAMPLE_NAME} RGPI=${INSERT_SIZE}

    # Index the alignment
    # Command: samtools index
    # Input: alignment (.bam)
    # Ouput: indexed alignment (.bam.bai)
    samtools index ${SAMPLE_NAME}_${RUN_ID}.sorted.bam


    # Append the file name (.bam) to the list of alignments that will be merged
    echo ${SAMPLE_NAME}_${RUN_ID}.sorted.RG.bam >> father.bamlist
done

# Merge the list of alignments into a single file
# Command: samtools merge
# Input: file containing the list of alignments (each line is a .bam file)
# Ouput: alignment (.bam)
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
samtools merge father.bamlist


# Index the alignment
# Command: samtools index
# Input: alignment (.sam or .bam)
# Ouput: indexed alignment (.sam.bai or .bam.bai)
samtools index father.bam

## PIPE = Alt maj L


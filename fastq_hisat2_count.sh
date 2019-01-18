#RNAseq Analysis

#Prepare Files
#awk -f foo.awk dict GCF_000001895.5_Rnor_6.0_genomic.gff  > genes.gff
#gffread Annotation/genes.gff -T -o Annotation/genes.gtf
#/mnt/b/said/NGI-RNAseq/bin/gtf2bed genes.gtf > genes.bed
hisat2_extract_splice_sites.py Annotation/genes.gtf > Annotation/splicesites.txt
hisat2_extract_exons.py Annotation/genes.gtf > Annotation/exons.txt
hisat2_extract_snps_haplotypes_VCF.py Sequence/genome.fa Annotation/Rattus_norvegicus.vcf Annotation/snps

hisat2-build --ss Annotation/splicesites.txt --exon Annotation/exons.txt --snp Annotation/snps Sequence/genome.fa Sequence/rn6_snp


#Concatening forward reads
cat 20171003_D12R1C1/Files/*R1_001.fastq.gz > D12R1C1.R1.fastq.gz
cat 20171003_D12R2/Files/*R1_001.fastq.gz > D12R2.R1.fastq.gz
cat 20171003_D12R3/Files/*R1_001.fastq.gz > D12R3.R1.fastq.gz
cat 20171003_D12RK2/Files/*R1_001.fastq.gz > D12RK2.R1.fastq.gz

cat 20171003_D18R1/Files/*R1_001.fastq.gz > D18R1.R1.fastq.gz
cat 20171003_D18R2CD2/Files/*R1_001.fastq.gz > D18R2CD2.R1.fastq.gz
cat 20171003_D18R3C1/Files/*R1_001.fastq.gz > D18R3C1.R1.fastq.gz
cat 20171003_D18R2C1/Files/*R1_001.fastq.gz > D18R2C1.R1.fastq.gz
cat 20171003_D18R3C2/Files/*R1_001.fastq.gz > D18R3C2.R1.fastq.gz


cat 20171003_NLR4/Files/*R1_001.fastq.gz > NLR4.R1.fastq.gz
cat 20171003_NLR6/Files/*R1_001.fastq.gz > NLR6.R1.fastq.gz
cat 20171003_NLRS/Files/*R1_001.fastq.gz > NLRS.R1.fastq.gz


#Concatening reverse reads
cat 20171003_D12R1C1/Files/*R2_001.fastq.gz > D12R1C1.R2.fastq.gz
cat 20171003_D12R2/Files/*R2_001.fastq.gz > D12R2.R2.fastq.gz
cat 20171003_D12R3/Files/*R2_001.fastq.gz > D12R3.R2.fastq.gz
cat 20171003_D12RK2/Files/*R2_001.fastq.gz > D12RK2.R2.fastq.gz

cat 20171003_D18R1/Files/*R2_001.fastq.gz > D18R1.R2.fastq.gz
cat 20171003_D18R2CD2/Files/*R2_001.fastq.gz > D18R2CD2.R2.fastq.gz
cat 20171003_D18R3C1/Files/*R2_001.fastq.gz > D18R3C1.R2.fastq.gz
cat 20171003_D18R2C1/Files/*R1_001.fastq.gz > D18R2C1.R1.fastq.gz
cat 20171003_D18R3C2/Files/*R2_001.fastq.gz > D18R3C2.R2.fastq.gz

cat 20171003_NLR4/Files/*R2_001.fastq.gz > NLR4.R2.fastq.gz
cat 20171003_NLR6/Files/*R2_001.fastq.gz > NLR6.R2.fastq.gz
cat 20171003_NLRS/Files/*R2_001.fastq.gz > NLRS.R2.fastq.gz


#fbef64b73570359e75eaaf68e3aed5f5  D12R1C1.R1.fastq.gz
#b32b8220c895eb39cf6537adb202b299  D12R1C1.R2.fastq.gz
#99a7fd6cacef99534cbf9feb35f3fc11  D12R2.R1.fastq.gz
#430771b5694c4d2458ad5d87a3def4d4  D12R2.R2.fastq.gz
#e15dc42055b76d737e071aa9ef438a85  D12R3.R1.fastq.gz
#9e00264c1612e7722e1c5f886c19b34c  D12R3.R2.fastq.gz
#e2da9e2791b843e522c50409ade4c5fe  D12RK2.R1.fastq.gz
#c2865a0f1e76601377d419f9df9f16c0  D12RK2.R2.fastq.gz
#c25e36393b279d3f7135bee4cca7d787  D18R1.R1.fastq.gz
#5cf6348072360310d2cd69bfe6fae98e  D18R1.R2.fastq.gz
#b11be5749c1d0318cab3c0c2a4d44af0  D18R2C1.R1.fastq.gz
#38388c51b8efad985a45ed2246e69863  D18R2C1.R2.fastq.gz
#c321b3fa22cb22532888d86d1f578071  D18R2CD2.R1.fastq.gz
#c54c9fc9ccf7409bb4e03cb0e3865094  D18R2CD2.R2.fastq.gz
#70719c08a4aee4a7f773efd9fba7a7f9  D18R3C1.R1.fastq.gz
#1b429e9b5f7133d5b595453fdb348338  D18R3C1.R2.fastq.gz
#d4fce6d6cdffd55f23fdf2563570e743  D18R3C2.R1.fastq.gz
#cb0a01a7a927e4538e6de617c14ed8e5  D18R3C2.R2.fastq.gz
#f956725a01ff1491e059e2bbe648173a  NLR4.R1.fastq.gz
#c940cc1bcf3de37b7f3ab7093d2d7435  NLR4.R2.fastq.gz
#698c10050f4a89d90ccd87267dc62be0  NLR6.R1.fastq.gz
#f49185dc363ee0e7647dcf51f00401b9  NLR6.R2.fastq.gz
#d0b033f9edb5245d6523f4690e020d29  NLRS.R1.fastq.gz
#4c273251589ff38395d162bf920d5523  NLRS.R2.fastq.gz


# STEP 1 - FastQC

fastqc *.fastq.gz 

multiqc .

fastqc --version
#FastQC v0.11.5

multiqc --version
#multiqc, version 1.0.dev0


#STEP 2 - Trim Galore!

trim_galore --paired --gzip D12R1C1.R?.fastq.gz &> D12R1C1.trimgalore.txt
trim_galore --paired --gzip D12R3.R?.fastq.gz &> D12R3.trimgalore.txt
trim_galore --paired --gzip D12RK2.R?.fastq.gz &> D12RK2.trimgalore.txt
trim_galore --paired --gzip D18R1.R?.fastq.gz &> D18R1.trimgalore.txt
trim_galore --paired --gzip D18R2C1.R?.fastq.gz &> D18R2C1.trimgalore.txt
trim_galore --paired --gzip D18R2CD2.R?.fastq.gz &> D18R2CD2.trimgalore.txt
trim_galore --paired --gzip D12R2.R?.fastq.gz &>D12R2.trimgalore.txt
trim_galore --paired --gzip D18R3C1.R?.fastq.gz &> D18R3C1.trimgalore.txt
trim_galore --paired --gzip D18R3C2.R?.fastq.gz &> D18R3C2.trimgalore.txt
trim_galore --paired --gzip NLR4.R?.fastq.gz &> NLR4.trimgalore.txt
trim_galore --paired --gzip NLR6.R?.fastq.gz &> NLR6.trimgalore.txt
trim_galore --paired --gzip NLRS.R?.fastq.gz &> NLRS.trimgalore.txt


# STEP 3 - align with HISAT2


#wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/rn6.tar.gz
hisat2index="/mnt/b/said/Ref/rn6/NCBI/Sequence/rn6_snp"


hisat2 --rna-strandness FR -q -p 4 -x $hisat2index -1 data/trimmed/D12R1C1.R1_val_1.fq.gz -2 data/trimmed/D12R1C1.R2_val_2.fq.gz | samtools view -bS -F 4 -F 256 - > bam/D12R1C1.bam
samtools sort -@ 3 bam/D12R1C1.bam  -o bam/D12R1C1.sorted.bam

hisat2 --rna-strandness FR -q -p 4 -x $hisat2index -1 data/trimmed/D12R2.R1_val_1.fq.gz -2 data/trimmed/D12R2.R2_val_2.fq.gz | samtools view -bS -F 4 -F 256 - > bam/D12R2.bam
samtools sort -@ 3 bam/D12R2.bam  -o bam/D12R2.sorted.bam

hisat2 --rna-strandness FR -q -p 4 -x $hisat2index -1 data/trimmed/D12R3.R1_val_1.fq.gz -2 data/trimmed/D12R3.R2_val_2.fq.gz | samtools view -bS -F 4 -F 256 - > bam/D12R3.bam
samtools sort -@ 3 bam/D12R3.bam  -o bam/D12R3.sorted.bam

hisat2 --rna-strandness FR -q -p 4 -x $hisat2index -1 data/trimmed/D12RK2.R1_val_1.fq.gz -2 data/trimmed/D12RK2.R2_val_2.fq.gz | samtools view -bS -F 4 -F 256 - > bam/D12RK2.bam
samtools sort -@ 3 bam/D12RK2.bam  -o bam/D12RK2.sorted.bam



hisat2 --rna-strandness FR -q -p 4 -x $hisat2index -1 data/trimmed/D18R1.R1_val_1.fq.gz -2 data/trimmed/D18R1.R2_val_2.fq.gz | samtools view -bS -F 4 -F 256 - > bam/D18R1.bam
samtools sort -@ 3 bam/D18R1.bam  -o bam/D18R1.sorted.bam

hisat2 --rna-strandness FR -q -p 4 -x $hisat2index -1 data/trimmed/D18R2C1.R1_val_1.fq.gz -2 data/trimmed/D18R2C1.R2_val_2.fq.gz | samtools view -bS -F 4 -F 256 - > bam/D18R2C1.bam
samtools sort -@ 3 bam/D18R2C1.bam  -o bam/D18R2C1.sorted.bam

hisat2 --rna-strandness FR -q -p 4 -x $hisat2index -1 data/trimmed/D18R3C2.R1_val_1.fq.gz -2 data/trimmed/D18R3C2.R2_val_2.fq.gz | samtools view -bS -F 4 -F 256 - > bam/D18R3C2.bam
samtools sort -@ 3 bam/D18R3C2.bam  -o bam/D18R3C2.sorted.bam

hisat2 --rna-strandness FR -q -p 4 -x $hisat2index -1 data/trimmed/D18R2CD2.R1_val_1.fq.gz -2 data/trimmed/D18R2CD2.R2_val_2.fq.gz | samtools view -bS -F 4 -F 256 - > bam/D18R2CD2.bam
samtools sort -@ 3 bam/D18R2CD2.bam  -o bam/D18R2CD2.sorted.bam




hisat2 --rna-strandness FR -q -p 4 -x $hisat2index -1 data/trimmed/D18R3C1.R1_val_1.fq.gz -2 data/trimmed/D18R3C1.R2_val_2.fq.gz | samtools view -bS -F 4 -F 256 - > bam/D18R3C1.bam
samtools sort -@ 3 bam/D18R3C1.bam  -o bam/D18R3C1.sorted.bam


hisat2 --rna-strandness FR -q -p 8 -x $hisat2index -1 data/trimmed/NLR4.R1_val_1.fq.gz -2 data/trimmed/NLR4.R2_val_2.fq.gz | samtools view -bS -F 4 -F 256 - > bam/NLR4.bam
samtools sort -@ 3 bam/NLR4.bam  -o bam/NLR4.sorted.bam

hisat2 --rna-strandness FR -q -p 8 -x $hisat2index -1 data/trimmed/NLR6.R1_val_1.fq.gz -2 data/trimmed/NLR6.R2_val_2.fq.gz | samtools view -bS -F 4 -F 256 - > bam/NLR6.bam
samtools sort -@ 3 bam/NLR6.bam  -o bam/NLR6.sorted.bam

hisat2 --rna-strandness FR -q -p 8 -x $hisat2index -1 data/trimmed/NLRS.R1_val_1.fq.gz -2 data/trimmed/NLRS.R2_val_2.fq.gz | samtools view -bS -F 4 -F 256 - > bam/NLRS.bam
samtools sort -@ 3 bam/NLRS.bam  -o bam/NLRS.sorted.bam



hisat2 --version
#/home/said/miniconda2/bin/hisat2-align-s version 2.1.0
#64-bit
#Built on login-node03
#Wed Jun  7 15:53:42 EDT 2017
#Compiler: gcc version 4.8.2 (GCC) 
#Options: -O3 -m64 -msse2 -funroll-loops -g3 -DPOPCNT_CAPABILITY
#Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}



# STEP 4 - RSeQC analysis
#wget ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Rattus_norvegicus/UCSC/rn6/Rattus_norvegicus_UCSC_rn6.tar.gz
baseName="D12R1C1"
baseName="D12RK2"
baseName="D12R2"
baseName="D18R2C1"
baseName="D12R3"
baseName="D18R1"
baseName="D18R2CD2"
baseName="D18R3C1"
baseName="D18R3C2"
baseName="NLR4"
baseName="NLR6"
baseName="NLRS"


bed12="/mnt/b/said/Ref/rn6/NCBI/Annotation/genes.bed"

bam_rseqc="bam/${baseName}.sorted.bam"

samtools index $bam_rseqc
infer_experiment.py -i $bam_rseqc -r $bed12 > bam/stats/${baseName}.infer_experiment.txt
junction_annotation.py -i $bam_rseqc -o algnm/stats/${baseName}.rseqc -r $bed12
bam_stat.py -i $bam_rseqc 2> bam/stats/${baseName}.bam_stat.txt
junction_saturation.py -i $bam_rseqc -o bam/stats/${baseName}.rseqc -r $bed12 2> bam/stats/${baseName}.junction_annotation_log.txt
inner_distance.py -i $bam_rseqc -o bam/stats/${baseName}.rseqc -r $bed12
read_distribution.py -i $bam_rseqc -r $bed12 > bam/stats/${baseName}.read_distribution.txt
read_duplication.py -i $bam_rseqc -o bam/stats/${baseName}.read_duplication

echo "Filename $bam_rseqc RseQC version: " $(read_duplication.py --version)


#Step 4.1 Rseqc genebody_coverage

cat <(samtools view -H ${bam_rseqc}) \
    <(samtools view ${bam_rseqc} | shuf -n 1000000) \
    | samtools sort - -o bam/${baseName}_subsamp_sorted.bam 

samtools index bam/${baseName}_subsamp_sorted.bam
geneBody_coverage.py -i bam/${baseName}_subsamp_sorted.bam -o bam/${baseName}.rseqc -r $bed12

# STEP 5 - preseq analysis

bam_preseq="bam/${baseName}.sorted.bam"

preseq lc_extrap -v -B $bam_preseq -o bam/${baseName}.ccurve.txt
echo "File name: $bam_preseq  preseq version: " $(preseq)

# STEP 6 Mark duplicates
#baseName="D12R3"
#baseName="D18R1"
#baseName="D18R2CD2"
#baseName="D18R3C1"
#baseName="D18R3C2"
#baseName="NLR4"
#baseName="NLR6"
#baseName="NLRS"

bam_markduplicates="bam/${baseName}.sorted.bam"


picard AddOrReplaceReadGroups I="$bam_markduplicates" \
   O="bam/${baseName}.gp.bam" \
   SO=coordinate \
   RGID=$baseName RGLB=$baseName RGPL=illumina RGPU=nextseq RGSM=$baseName 

#picard MarkDuplicates I=.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics 


java -Xms512m -Xmx30g -jar /home/said/miniconda2/share/picard-1.126-4/picard.jar MarkDuplicates \
    INPUT="bam/${baseName}.gp.bam"  \
    OUTPUT=bam/${baseName}.markDups.bam \
    METRICS_FILE=bam/${baseName}.markDups_metrics.txt \
    REMOVE_DUPLICATES=false \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT

echo "File name: $bam_markduplicates Picard version " $(picard MarkDuplicates --version 2>&1)
samtools index bam/${baseName}.markDups.bam

ref=/mnt/b/said/Ref/rn6/NCBI/Sequence/genome.fa

gatk \
                SplitNCigarReads \
                -R ${ref} \
                -I "${sample}.markDups.bam" \
                -O "${sample}.gatk.bam" 

vcf="/mnt/b/said/Ref/rn6/NCBI/Annotation/rn6.vcf"

gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R ${ref} \
            -I "${sample}.gatk.bam" \
            --use-original-qualities \
            -O ${sample}.recal_output_file \
            -known-sites ${vcf} 

gatk \
            --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R ${ref} \
            -I ${sample}.gatk.bam \
            --use-original-qualities \
            -O ${sample}.hc.bam \
            --bqsr-recal-file ${sample}.recalibration_report

bed12="/mnt/b/said/Ref/rn6/NCBI/Annotation/genes.bed"


gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    HaplotypeCaller \
    -R ${ref} \
    -I ${sample}.hc.bam \
    -L $bed12 \
    -O ${sample}.vcf.gz \
    -dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20 \
    --dbsnp ${vcf}

gatk \
        VariantFiltration \
      --R ${ref} \
      --V ${sample}.vcf.gz \
      --window 35 \
      --cluster 3 \
      --filter-name "FS" \
      --filter "FS > 30.0" \
      --filter-name "QD" \
      --filter "QD < 2.0" \
      -O ${sample}.filtro.vcf.gz

gatk -Xmx30G -T SplitNCigarReads -R $ref -I "bam/${baseName}.markDups.bam" -o "bam/${baseName}.gatk.bam" -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

#5. Base Recalibration
#picard SortVcf \
#    I="/mnt/b/said/Ref/rn6/NCBI/Annotation/Rattus_norvegicus.vcf" \
#    O="/mnt/b/said/Ref/rn6/NCBI/Annotation/rn6.vcf" \
#    SEQUENCE_DICTIONARY=/mnt/b/said/Ref/rn6/NCBI/Sequence/genome.dict
#awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' no_chr.vcf > with_chr.vcf

vcf="/mnt/b/said/Ref/rn6/NCBI/Annotation/rn6.vcf"


gatk -T BaseRecalibrator -nct 3 -R $ref -I "bam/${baseName}.gatk.bam" -knownSites "$vcf" -o "bam/${baseName}_recal.table"
gatk -T BaseRecalibrator -nct 2 -R $ref -I "bam/${baseName}.gatk.bam" -knownSites "$vcf" -BQSR "bam/${baseName}_recal.table" -o "bam/${baseName}_post_recal.table"	
gatk -T AnalyzeCovariates -R $ref -before "bam/${baseName}_recal.table" -after "bam/${baseName}_post_recal.table" -plots "bam/${baseName}_recalibration_plots.pdf"	
gatk -Xmx5G -T PrintReads -nct 2 -R $ref -I "bam/${baseName}.gatk.bam" -BQSR "bam/${baseName}_recal.table" -o "bam/${baseName}.hc.bam"

gatk -Xmx5G -T HaplotypeCaller -nct 2 -R $ref -I "bam/${baseName}.hc.bam" -dontUseSoftClippedBases -stand_call_conf 20.0 -o "vcf/${baseName}.raw.vcf"
gatk -T VariantFiltration -R $ref -V "vcf/${baseName}.raw.vcf" -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o "vcf/${baseName}.vcf"





#  STEP 7 - dupRadar

gtf="/mnt/b/said/Ref/rn6/NCBI/Annotation/genes.gtf"
bam_md="bam/${baseName}.hc.bam"

/mnt/b/said/NGI-RNAseq/bin/dupRadar.r $bam_md $gtf TRUE


# STEP 8 Feature counts

#==> bam/stats/D12R1C1.infer_experiment.txt <==

#This is PairEnd Data
#Fraction of reads failed to determine: 0.0000
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.0629
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.9371

baseName="D12R1C1"
bam_featurecounts="bam/${baseName}.sorted.bam"
featureCounts_direction="1"
featureCounts -a $gtf -g transcript_id -T 10 -o analysis/${baseName}_exon.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  

featureCounts -a $gtf -g gene_id -T 3 -o analysis/${baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  
featureCounts -a $gtf -g gene_biotype -o analysis/${baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 analysis/${baseName}_biotype.featureCounts.txt > analysis/${baseName}_biotype_counts.txt


#==> bam/stats/D12R2.infer_experiment.txt <==


#This is PairEnd Data
#Fraction of reads failed to determine: 0.0000
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.0814
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.9186

baseName="D12R2"
bam_featurecounts="bam/${baseName}.sorted.bam"
featureCounts_direction="1"
featureCounts -a $gtf -g transcript_id -T 10 -o analysis/${baseName}_exon.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  

featureCounts -a $gtf -g gene_id -o analysis/${baseName}_gene.featureCounts.txt -p -T 10 -s $featureCounts_direction $bam_featurecounts  
featureCounts -a $gtf -g gene_biotype -o analysis/${baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 analysis/${baseName}_biotype.featureCounts.txt > analysis/${baseName}_biotype_counts.txt


#==> bam/stats/D12R3.infer_experiment.txt <==


#This is PairEnd Data
#Fraction of reads failed to determine: 0.0000
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.8559
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.1441
gtf="/mnt/b/said/Ref/rn6/Annotation/Genes/genes.gtf"

baseName="D12R3"
bam_featurecounts="bam/${baseName}.sorted.bam"
featureCounts_direction="1"
featureCounts -a $gtf -g transcript_id -T 10 -o analysis/${baseName}_exon.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  

featureCounts -a $gtf -g gene_id -o analysis/${baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  
featureCounts -a $gtf -g gene_biotype -o analysis/${baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 analysis/${baseName}_biotype.featureCounts.txt > analysis/${baseName}_biotype_counts.txt



#==> bam/stats/D12RK2.infer_experiment.txt <==


#This is PairEnd Data
#Fraction of reads failed to determine: 0.0000
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.3448
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.6552

baseName="D12RK2"
bam_featurecounts="bam/${baseName}.sorted.bam"
featureCounts_direction="1"
featureCounts -a $gtf -g transcript_id -T 10 -o analysis/${baseName}_exon.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  

featureCounts -a $gtf -g gene_id -o analysis/${baseName}_gene.featureCounts.txt -p -T 10 -s $featureCounts_direction $bam_featurecounts  
featureCounts -a $gtf -g gene_biotype -o analysis/${baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 analysis/${baseName}_biotype.featureCounts.txt > analysis/${baseName}_biotype_counts.txt


#==> bam/stats/D18R1.infer_experiment.txt <==


#This is PairEnd Data
#Fraction of reads failed to determine: 0.0000
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.9434
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.0566

gtf="/mnt/b/said/Ref/rn6/Annotation/Genes/genes.gtf"

baseName="D18R1"
bam_featurecounts="bam/${baseName}.sorted.bam"
featureCounts_direction="1"
featureCounts -a $gtf -g transcript_id -T 10 -o analysis/${baseName}_exon.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  

featureCounts -a $gtf -g gene_id -o analysis/${baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  
featureCounts -a $gtf -g gene_biotype -o analysis/${baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 analysis/${baseName}_biotype.featureCounts.txt > analysis/${baseName}_biotype_counts.txt


#==> bam/stats/D18R2C1.infer_experiment.txt <==


#This is PairEnd Data
#Fraction of reads failed to determine: 0.0000
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.2715
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.7285

baseName="D18R2C1"
bam_featurecounts="bam/${baseName}.sorted.bam"
featureCounts_direction="1"
featureCounts -a $gtf -g transcript_id -T 10 -o analysis/${baseName}_exon.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  


featureCounts -a $gtf -g gene_id -o analysis/${baseName}_gene.featureCounts.txt -p -T 10 -s $featureCounts_direction $bam_featurecounts  
featureCounts -a $gtf -g gene_biotype -o analysis/${baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 analysis/${baseName}_biotype.featureCounts.txt > analysis/${baseName}_biotype_counts.txt


#==> bam/stats/D18R2CD2.infer_experiment.txt <==


#This is PairEnd Data
#Fraction of reads failed to determine: 0.0000
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.8847
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.1153

gtf="/mnt/b/said/Ref/rn6/Annotation/Genes/genes.gtf"

baseName="D18R2CD2"
bam_featurecounts="bam/${baseName}.sorted.bam"
featureCounts_direction="1"
featureCounts -a $gtf -g transcript_id -T 10 -o analysis/${baseName}_exon.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  


featureCounts -a $gtf -g gene_id -o analysis/${baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  
featureCounts -a $gtf -g gene_biotype -o analysis/${baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 analysis/${baseName}_biotype.featureCounts.txt > analysis/${baseName}_biotype_counts.txt


#==> bam/stats/D18R3C1.infer_experiment.txt <==


#This is PairEnd Data
#Fraction of reads failed to determine: 0.0000
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.8094
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.1906

gtf="/mnt/b/said/Ref/rn6/Annotation/Genes/genes.gtf"

baseName="D18R3C1"
bam_featurecounts="bam/${baseName}.sorted.bam"
featureCounts_direction="1"
featureCounts -a $gtf -g transcript_id -T 10 -o analysis/${baseName}_exon.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  


featureCounts -a $gtf -g gene_id -o analysis/${baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  
featureCounts -a $gtf -g gene_biotype -o analysis/${baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 analysis/${baseName}_biotype.featureCounts.txt > analysis/${baseName}_biotype_counts.txt


#==> bam/stats/D18R3C2.infer_experiment.txt <==


#This is PairEnd Data
#Fraction of reads failed to determine: 0.0000
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.8198
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.1802
gtf="/mnt/b/said/Ref/rn6/Annotation/Genes/genes.gtf"

baseName="D18R3C2"
bam_featurecounts="bam/${baseName}.sorted.bam"
featureCounts_direction="1"
featureCounts -a $gtf -g transcript_id -T 10 -o analysis/${baseName}_exon.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  


featureCounts -a $gtf -g gene_id -o analysis/${baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  
featureCounts -a $gtf -g gene_biotype -o analysis/${baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 analysis/${baseName}_biotype.featureCounts.txt > analysis/${baseName}_biotype_counts.txt


#==> bam/stats/NLR4.infer_experiment.txt <==


#This is PairEnd Data
#Fraction of reads failed to determine: 0.0000
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.7365
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.2635
gtf="/mnt/b/said/Ref/rn6/Annotation/Genes/genes.gtf"

baseName="NLR4"
bam_featurecounts="bam/${baseName}.sorted.bam"
featureCounts_direction="1"
featureCounts -a $gtf -g transcript_id -T 10 -o analysis/${baseName}_exon.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  


featureCounts -a $gtf -g gene_id -o analysis/${baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  
featureCounts -a $gtf -g gene_biotype -o analysis/${baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 analysis/${baseName}_biotype.featureCounts.txt > analysis/${baseName}_biotype_counts.txt


#==> bam/stats/NLR6.infer_experiment.txt <==


#This is PairEnd Data
#Fraction of reads failed to determine: 0.0000
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.9358
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.0642


baseName="NLR6"
bam_featurecounts="bam/${baseName}.sorted.bam"
featureCounts_direction="1"


gtf="/mnt/b/said/Ref/rn6/NCBI/Annotation/genes.gtf"
featureCounts -a $gtf -g gene_name -T 10 -o analysis/${baseName}_gene_name.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  


featureCounts -a $gtf -g gene_id -o analysis/${baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  
featureCounts -a $gtf -g gene_biotype -o analysis/${baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 analysis/${baseName}_biotype.featureCounts.txt > analysis/${baseName}_biotype_counts.txt


#==> bam/stats/NLRS.infer_experiment.txt <==


#This is PairEnd Data
#Fraction of reads failed to determine: 0.0000
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.8608
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.1392
gtf="/mnt/b/said/Ref/rn6/Annotation/Genes/genes.gtf"

baseName="NLRS"
bam_featurecounts="bam/${baseName}.sorted.bam"
featureCounts_direction="1"
featureCounts -a $gtf -g exon_id -o analysis/${baseName}_exon.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  


featureCounts -a $gtf -g gene_id -o analysis/${baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  
featureCounts -a $gtf -g gene_biotype -o analysis/${baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 analysis/${baseName}_biotype.featureCounts.txt > analysis/${baseName}_biotype_counts.txt





bam_featurecounts="bam/${baseName}.sorted.bam"
featureCounts_direction="1"

featureCounts -a $gtf -g gene_id -o analysis/${baseName}_gene.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  

featureCounts -a $gtf -g exon_id -o analysis/${baseName}_exon.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts  



featureCounts -a $gtf -g gene_biotype -o analysis/${baseName}_biotype.featureCounts.txt -p -s $featureCounts_direction $bam_featurecounts
cut -f 1,7 analysis/${baseName}_biotype.featureCounts.txt > analysis/${baseName}_biotype_counts.txt


bash mapQC.sh "D12R1C1" &> analysis/D12R1C1.log
bash mapQC.sh "D12RK2" &> analysis/D12RK2.log
bash mapQC.sh "D18R2C1" &> analysis/D18R2C1.log
bash mapQC.sh "D12R2" &> analysis/D12R2.log


bash mapQC.sh "D12R3" &> analysis/D12R3.log
bash mapQC.sh "D18R1" &> analysis/D18R1.log
bash mapQC.sh "D18R2CD2" &> analysis/D18R2CD2.log
bash mapQC.sh "D18R3C1" &> analysis/D18R3C1.log
bash mapQC.sh "D18R3C2" &> analysis/D18R3C2.log
bash mapQC.sh "NLR4" &> analysis/NLR4.log
bash mapQC.sh "NLR6" &> analysis/NLR6.log
bash mapQC.sh "NLRS" &> analysis/NLRS.log



# STEP 9 - Merge featurecounts
/mnt/b/said/NGI-RNAseq/bin/merge_featurecounts.py -o analysis/merged_gene_name_counts.txt -i analysis/*_gene_name.featureCounts.txt
/mnt/b/said/NGI-RNAseq/bin/merge_featurecounts.py -o analysis/merged_transcripts_counts.txt -i analysis/*exon.featureCounts.txt


# STEP 10 - stringtie FPKM

baseName="D12RK2"
#baseName="D12R1C1"
#baseName="D12R2"
#baseName="D18R2C1"
#baseName="D12R3"
#baseName="D18R1"
#baseName="D18R2CD2"
#baseName="D18R3C1"
#baseName="D18R3C2"
#baseName="NLR4"
#baseName="NLR6"
#baseName="NLRS"


bam_stringtieFPKM="bam/${baseName}.sorted.bam"

gtf="/mnt/b/said/Ref/rn6/Annotation/Genes/genes.gtf"


stringtie $bam_stringtieFPKM \
        -o analysis/${baseName}_transcripts.gtf \
        -G $gtf \
        -A analysis/${baseName}.gene_abund.txt \
        -C analysis/${baseName}.cov_refs.gtf \
        -e \
        -b analysis/${baseName}_ballgown



stringtie --merge -p 32 -G $gtf \
     -o stringtie_merged.gtf  mergelist.txt

gffcompare –r $gtf –G –o merged stringtie_merged.gtf


baseName="D12RK2"
#baseName="D12R1C1"
#baseName="D12R2"
#baseName="D18R2C1"
#baseName="D12R3"
#baseName="D18R1"
#baseName="D18R2CD2"
#baseName="D18R3C1"
#baseName="D18R3C2"
#baseName="NLR4"
#baseName="NLR6"
#baseName="NLRS"
bam_stringtieFPKM="bam/${baseName}.sorted.bam"

stringtie $bam_stringtieFPKM \
        -o analysis/${baseName}_stringtie_transcripts.gtf \
        -G analysis/stringtie_merged.gtf \
        -p 3 \
        -A analysis/${baseName}.stringtie.gene_abund.txt \
        -C analysis/${baseName}.stringtie.cov_refs.gtf \
        -e \
        -b analysis/ballgown/${baseName}_stringtie_ballgown





# STEP 11 - edgeR MDS and heatmap
/mnt/b/said/NGI-RNAseq/bin/edgeR_heatmap_MDS.r  bam/*sorted.bam




# STEP 12 MultiQC
multiqc -f $rtitle $rfilename --config $multiqc_config .








# STEP 13 GATK

baseName="D12RK2"
#baseName="D12R1C1"
#baseName="D12R2"
#baseName="D18R2C1"
#baseName="D12R3"
#baseName="D18R1"
#baseName="D18R2CD2"
#baseName="D18R3C1"
#baseName="D18R3C2"
#baseName="NLR4"
#baseName="NLR6"
#baseName="NLRS"

bam_markduplicates="bam/${baseName}.markDups.bam"

ref=/mnt/b/said/Ref/rn6/rn6.fa




picard AddOrReplaceReadGroups I=$bam_markduplicates \
   O="bam/${baseName}.gatk.bam" \
   SO=coordinate \
   RGID=$baseName RGLB=$baseName RGPL=illumina RGPU=nextseq RGSM=$baseName 

#picard MarkDuplicates I=.bam O=dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics 
gatk \
  SplitNCigarReads \
  -R ${ref_fasta} \
  -I ${input_bam} \
  -O ${base_name}.bam 

gatk -Xmx30G -T SplitNCigarReads -R $ref -I "bam/${baseName}.gatk.bam" -o "bam/${baseName}.gatk2.bam" -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

#5. Base Recalibration
vcf="/mnt/b/said/Ref/rn6/Rattus_norvegicus.vcf.wo"
#awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' no_chr.vcf > with_chr.vcf

gatk -T BaseRecalibrator -nct 12 -R $ref -I "bam/${baseName}.gatk2.bam" -knownSites "$vcf" -o "bam/${baseName}_recal.table"
gatk -T BaseRecalibrator -nct 12 -R $ref -I "bam/${baseName}.gatk2.bam" -knownSites "$vcf" -BQSR "bam/${baseName}_recal.table" -o "bam/${baseName}_post_recal.table"	
gatk -T AnalyzeCovariates -R $ref -before "bam/${baseName}_recal.table" -after "bam/${baseName}_post_recal.table" -plots "bam/${baseName}_recalibration_plots.pdf"	
gatk -Xmx30G -T PrintReads -nct 12 -R $ref -I "bam/${baseName}.gatk2.bam" -BQSR "bam/${baseName}_recal.table" -o "bam/${baseName}.hc.bam"

gatk -Xmx30G -T HaplotypeCaller -nct 12 -R $ref -I "bam/${baseName}.hc.bam" -dontUseSoftClippedBases -stand_call_conf 20.0 -o "vcf/${baseName}.raw.vcf"
gatk -T VariantFiltration -R $ref -V "vcf/${baseName}.raw.vcf" -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o "vcf/${baseName}.vcf"





sample="D12RK2"
sample="D12R1C1"
sample="D12R2"
sample="D18R2C1"
sample="D12R3"
sample="D18R1"
sample="D18R2CD2"
sample="D18R3C1"
sample="D18R3C2"
sample="NLR4"
sample="NLR6"
sample="NLRS"

cd bam/
source activate gatk4

ref=/mnt/b/said/Ref/rn6/NCBI/Sequence/genome.fa

gatk \
                SplitNCigarReads \
                -R ${ref} \
                -I "${sample}.markDups.bam" \
                -O "${sample}.gatk.bam" 

vcf="/mnt/b/said/Ref/rn6/NCBI/Annotation/rn6.vcf"

gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R ${ref} \
            -I "${sample}.gatk.bam" \
            --use-original-qualities \
            -O ${sample}.recal_output_file \
            -known-sites ${vcf} 

gatk \
            --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R ${ref} \
            -I ${sample}.gatk.bam \
            --use-original-qualities \
            -O ${sample}.hc.bam \
            --bqsr-recal-file ${sample}.recalibration_report

bed12="/mnt/b/said/Ref/rn6/NCBI/Annotation/genes.bed"


gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    HaplotypeCaller \
    -R ${ref} \
    -I ${sample}.hc.bam \
    -L $bed12 \
    -O ${sample}.vcf.gz \
    -dont-use-soft-clipped-bases \
    --standard-min-confidence-threshold-for-calling 20 \
    --dbsnp ${vcf}

gatk \
        VariantFiltration \
      --R ${ref} \
      --V ${sample}.vcf.gz \
      --window 35 \
      --cluster 3 \
      --filter-name "FS" \
      --filter "FS > 30.0" \
      --filter-name "QD" \
      --filter "QD < 2.0" \
      -O ${sample}.filtro.vcf.gz




#!/usr/bin/env bash
# Merge all reads (3 samples) from the same strain
## Reads 1
raw="./Cauris/raw"
input=("${raw}/F1_185_1" "${raw}/F4_188_1" "${raw}/F7_189_1" "${raw}/F10_191_1" "${raw}/S4_48_1")
for strain in "${input[@]}"; do
    strain=$(basename "${strain}")
    echo "${strain}"
    ### Strain 185
    if [[ ${strain} =~ "185" ]]; then
        cat \
            "${raw}/F1_185_1/F1_185_1_1.fq.gz" \
            "${raw}/F2_185_2/F2_185_2_1.fq.gz" \
            "${raw}/F3_185_3/F3_185_3_1.fq.gz" \
            > "${raw}/merged_reads/strain_185_R1.fq.gz"
    ### Strain 188
    elif [[ ${strain} =~ "188" ]]; then
        cat \
            "${raw}/F4_188_1/F4_188_1_1.fq.gz" \
            "${raw}/F5_188_2/F5_188_2_1.fq.gz" \
            "${raw}/F6_188_3/F6_188_3_1.fq.gz" \
            > "${raw}/merged_reads/strain_188_R1.fq.gz"
    ### Strain 189
    elif [[ ${strain} =~ "189" ]]; then
        cat \
            "${raw}/F7_189_1/F7_189_1_1.fq.gz" \
            "${raw}/F8_189_2/F8_189_2_1.fq.gz" \
            "${raw}/F9_189_3/F9_189_3_1.fq.gz" \
            > "${raw}/merged_reads/strain_189_R1.fq.gz"
    ### Strain 191
    elif [[ ${strain} =~ "191" ]]; then
        cat \
            "${raw}/F10_191_1/F10_191_1_1.fq.gz" \
            "${raw}/F11_191_2/F11_191_2_1.fq.gz" \
            "${raw}/F12_191_3/F12_191_3_1.fq.gz" \
            > "${raw}/merged_reads/strain_191_R1.fq.gz"
    ### Strain 48
    elif [[ ${strain} =~ "48" ]]; then
        cat \
            "${raw}/S4_48_1/S4_48_1_1.fq.gz" \
            "${raw}/S5_48_2/S5_48_2_1.fq.gz" \
            "${raw}/S6_48_3/S6_48_3_1.fq.gz" \
            > "${raw}/merged_reads/strain_48_R1.fq.gz"
    fi
done

## Reads 2
for strain in "${input[@]}"; do
    strain=$(basename "${strain}")
    echo "${strain}"
    ### Strain 185
    if [[ ${strain} =~ "185" ]]; then
        cat \
            "${raw}/F1_185_1/F1_185_1_2.fq.gz" \
            "${raw}/F2_185_2/F2_185_2_2.fq.gz" \
            "${raw}/F3_185_3/F3_185_3_2.fq.gz" \
            > "${raw}/merged_reads/strain_185_R2.fq.gz"
    ### Strain 188
    elif [[ ${strain} =~ "188" ]]; then
        cat \
            "${raw}/F4_188_1/F4_188_1_2.fq.gz" \
            "${raw}/F5_188_2/F5_188_2_2.fq.gz" \
            "${raw}/F6_188_3/F6_188_3_2.fq.gz" \
            > "${raw}/merged_reads/strain_188_R2.fq.gz"
    ### Strain 189
    elif [[ ${strain} =~ "189" ]]; then
        cat \
            "${raw}/F7_189_1/F7_189_1_2.fq.gz" \
            "${raw}/F8_189_2/F8_189_2_2.fq.gz" \
            "${raw}/F9_189_3/F9_189_3_2.fq.gz" \
            > "${raw}/merged_reads/strain_189_R2.fq.gz"
    ### Strain 191
    elif [[ ${strain} =~ "191" ]]; then
        cat \
            "${raw}/F10_191_1/F10_191_1_2.fq.gz" \
            "${raw}/F11_191_2/F11_191_2_2.fq.gz" \
            "${raw}/F12_191_3/F12_191_3_2.fq.gz" \
            > "${raw}/merged_reads/strain_191_R2.fq.gz"
    ### Strain 48
    elif [[ ${strain} =~ "48" ]]; then
        cat \
            "${raw}/S4_48_1/S4_48_1_2.fq.gz" \
            "${raw}/S5_48_2/S5_48_2_2.fq.gz" \
            "${raw}/S6_48_3/S6_48_3_2.fq.gz" \
            > "${raw}/merged_reads/strain_48_R2.fq.gz"
    fi
done

# Paths set-up
REF_DIR="refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Sequence/WholeGenomeFasta"
RG="Cauris_merged_output/RG_SplitReads/RG"
MD="Cauris_merged_output/MarkDuplicates"
SPLIT="Cauris_merged_output/RG_SplitReads/SplitReads"
BQSR="/mnt/rdisk/dminh/Cauris/Cauris_merged_output/BaseRecalibrator/recal_table"
BQSR2="Cauris_merged_output/BaseRecalibrator/recal_table2"
RECAL_BAM="Cauris_merged_output/BaseRecalibrator/recalibrated"
COV="Cauris_merged_output/BaseRecalibrator/covariate"
GVCF="Cauris_merged_output/HaplotypeCaller/gvcf"
COMBINED="Cauris_merged_output/HaplotypeCaller/combinegvcfs"
GENO="Cauris_merged_output/HaplotypeCaller/genotypegvcfs"
SNP="Cauris_merged_output/HaplotypeCaller/selectedsnps"
SNP_FILTERED="Cauris_merged_output/HaplotypeCaller/selectedsnpsfiltered"
RE_BAM="/mnt/rdisk/dminh/Cauris/Cauris_merged_output/HaplotypeCaller/realigned_bam"
BAM_COV="Cauris_merged_output/BaseRecalibrator/recalibrated/mosdepth"

# FastQC
fastqc \
    --quiet \
    --threads 10 \
    -o "Cauris_merged_output/FastQC/"
    "raw/merged_reads/*.fq.gz"

# Trim Galore
for read1 in raw/merged_reads/*_R1.fq.gz; do
    read1=$(basename "${read1}")
    read2=${read1/R1/R2}
    trim_galore \
        -j 8 \
        --paired \
        --fastqc \
        --gzip \
        --stringency 1 \
        --length 20 \
        "${read1} ${read2}"
done

# STAR align
for trimmed1 in Cauris_merged_output/Trim_Galore/*R1_val_1*; do
    trimmed1=$(basename "${trimmed1}")
    trimmed2=${trimmed1/R1_val_1/R2_val_2}
    sample_name=${trimmed1/_R1_val_1.fq.gz/}
    STAR --genomeDir STARIndex \
        --readFilesIn "${trimmed1} ${trimmed2}" \
        --runThreadN 32 \
        --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 64324509440 \
        --readFilesCommand zcat \
        --outSAMmultNmax 1 \
        --twopassMode Basic \
        --outFileNamePrefix "${sample_name}" \
        --peOverlapNbasesMin 10 \
        --peOverlapMMp 0.01
    
    samtools index "${sample_name}Aligned.sortedByCoord.out.bam"

done

# Adding read group (RG) information
for file in $(ls ${MD}/*.bam); do
    sample=$(basename "${file}")
    echo "${sample}"
    gatk --java-options "-Xmx128G" AddOrReplaceReadGroups \
        -I "${MD}/${sample}" \
        -O "${RG}/${sample/.markDups.bam/.md.RG.bam}" \
        --RGLB "lib1" \
        --RGPL "Illumina" \
        --RGPU "unit1" \
        --RGSM "${sample}"
done

# Handling Splicing Events in RNASeq Data using SplitNCigarReads: Split reads into exon and intron segments
for bam in $(ls ${RG}); do 
    sample=$(basename "${bam}")
    gatk SplitNCigarReads \
        --reference "${REF_DIR}/genome.fa" \
        --input "${RG}/${sample}" \
        --output "${SPLIT}/${sample/.md.RG.bam/.md.RG.SplitReads.bam}"
done

# BQSR using adjusted SNPs vcf files
## 1. BaseRecalibrator (First pass)
for file in $(ls ${SPLIT}/*.bam); do
    BAM=$(basename "${file}")
    gatk --java-options "-Xmx256G" BaseRecalibrator \
        --reference "${REF_DIR}/genome.fa" \
        --input "${SPLIT}/${BAM}" \
        --output "${BQSR}/${BAM/md.RG.SplitReads.bam/recal.table}" \
        --known-sites "./refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/VCF/B11221_CA05_contigs_rename.vcf.gz" \
        --known-sites "./refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/VCF/B11245_CA06_contigs_rename.vcf.gz" \
        --known-sites "./refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/VCF/cauris.consensus_contigs_rename_split_edited.vcf.gz"
done
## 2. ApplyBQSR
for file in $(ls ${SPLIT}/*.bam); do
    BAM=$(basename "${file}")
    gatk --java-options "-Xmx256G" ApplyBQSR \
        --reference "${REF_DIR}/genome.fa" \
        --input "${SPLIT}/${BAM}" \
        --bqsr-recal-file "${BQSR}/${BAM/md.RG.SplitReads.bam/recal.table}" \
        --output "${RECAL_BAM}/${BAM/md.RG.SplitReads/recal}"
done
## 3. Second pass BQSR on the recalibrated bam files
for file in $(ls ${RECAL_BAM}/*.bam); do
    BAM=$(basename "${file}")
    gatk --java-options "-Xmx256G" BaseRecalibrator \
        --reference "${REF_DIR}/genome.fa" \
        --input "${RECAL_BAM}/${BAM}" \
        --output "${BQSR2}/${BAM/recal.bam/recal.table2}" \
        --known-sites "./refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/VCF/B11221_CA05_contigs_rename.vcf.gz" \
        --known-sites "./refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/VCF/B11245_CA06_contigs_rename.vcf.gz" \
        --known-sites "./refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/VCF/cauris.consensus_contigs_rename_split_edited.vcf.gz"
done
## 4. Evaluate and compare base quality score recalibration tables
for TABLE in $(ls ${BQSR}); do
    gatk --java-options "-Xmx256G" AnalyzeCovariates \
        -before "${BQSR}/${TABLE}" \
        -after "${BQSR2}/${TABLE/table/table2}" \
        -plots "${COV}/${TABLE/recal.table/covariates.pdf}" \
        -csv "${COV}/${TABLE/recal.table/covariates.csv}"
done

#----------------------------------------------------------------------------------------------------------------#

# Germline variant calling with HaplotypeCaller in gVCF (Phred score confident = 30)
for file in $(ls ${RECAL_BAM}/*.bam); do
    BAM=$(basename "${file}")
    gatk --java-options "-Xmx256G" HaplotypeCaller \
        -R "${REF_DIR}/genome.fa" \
        -I "${RECAL_BAM}/${BAM}" \
        -O "${GVCF}/${BAM/recal.bam/haplotypecaller.g.vcf.gz}" \
        --dont-use-soft-clipped-bases \
        -stand-call-conf 30 \
        -ERC GVCF --sample-ploidy "1"
done

# GATK4_LOCALCOMBINEGVCFS: Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file
gatk --java-options "-Xmx256g" CombineGVCFs \
      -R "${REF_DIR}/genome.fa" \
      -O "${COMBINED}/cauris.combined.g.vcf.gz" \
      -V "${GVCF}/strain_185.haplotypecaller.g.vcf.gz" -V "${GVCF}/strain_188.haplotypecaller.g.vcf.gz" -V "${GVCF}/strain_189.haplotypecaller.g.vcf.gz" -V "${GVCF}/strain_191.haplotypecaller.g.vcf.gz" -V "${GVCF}/strain_48.haplotypecaller.g.vcf.gz"

# GATK4_GENOTYPEGVCFS: Perform joint genotyping on one or more samples pre-called with HaplotypeCaller
gatk --java-options "-Xmx256g" GenotypeGVCFs \
    -R "${REF_DIR}/genome.fa" \
    -V "${COMBINED}/cauris.combined.g.vcf.gz" \
    -O "${GENO}/cauris.combined.genotype.vcf.gz"

# Extract only SNP variants
gatk --java-options "-Xmx64G" SelectVariants \
    -R "${REF_DIR}/genome.fa" \
    -V "${GENO}/cauris.combined.genotype.vcf.gz" \
    -O "${SNP}/cauris.combined.genotype.snps.selectvariants.vcf.gz" \
    --select-type-to-include "SNP"

# SNPs hard filtering
gatk --java-options "-Xmx64G" VariantFiltration \
        -R "${REF_DIR}/genome.fa" \
        -V "${SNP}/cauris.combined.genotype.snps.selectvariants.vcf.gz" \
        -O "${SNP_FILTERED}/cauris.combined.genotype.snps.filtered.vcf.gz" \
        --filter-name "QD_filter" -filter "QD < 2.0" \
        --filter-name "FS_filter" -filter "FS > 60.0" \
        --filter-name "MQ_filter" -filter "MQ < 40.0" \
        --filter-name "SOR_filter" -filter "SOR > 3.0" \
        --filter-name "MQRankSum_filter" -filter "MQRankSum <-12.5" \
        --filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum <-8.0"

# Emit the assembled haplotypes and locally realigned reads
for file in $(ls ${RECAL_BAM}/*.bam); do
    BAM=$(basename "${file}")
    gatk --java-options "-Xmx256G" HaplotypeCaller \
        -R "${REF_DIR}/genome.fa" \
        -I "${RECAL_BAM}/${BAM}" \
        -O "${RE_BAM}/${BAM/recal.bam/haplotypecaller.g.vcf.gz}" \
        -bamout "${RE_BAM}/${BAM/recal/bamout}" \
        --dont-use-soft-clipped-bases \
        -stand-call-conf 30 \
        -ERC GVCF --sample-ploidy "1"
done

# BAM coverage QC
for bam in $(ls ${RECAL_BAM} | grep -vE "qualimap|bai|mosdepth"); do
    BAM=${bam/.recal.bam/}
    echo "${BAM}"
    MOSDEPTH_PRECISION=5 mosdepth \
        --threads 40 \
        --fasta "refs/Candida_auris_B8441/Genbank/Cand_auris_B8441_V2/Sequence/WholeGenomeFasta/genome.fa" \
        -n --fast-mode --by 414 \
        "${BAM_COV}/${BAM}" \
        "${RECAL_BAM}/${bam}"
done
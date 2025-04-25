#!/bin/bash
# RNA-seq preprocessing pipeline

# CREACIÓN DE ÁRBOL DE DIRECTORIOS
mkdir -p input/ out/fastqc/ out/fastp/ out/fastq_screen/ \
	out/multiqc/ \
	out/hisat2/ \
	out/htseq/ \
	log/fastp/ \
	log/hisat2/ \
	res/ \
	scripts/

# REVISAR ESTRUCTURA
tree

# CONTROL DE CALIDAD INICIAL
fastqc -o out/fastqc \
	input/*.fastq

# ELIMINACIÓN DE ADAPTADORES Y EXTREMOS CON BAJA CALIDAD
fastp -i input/SRR479052.chr21_1.fastq \
	-I input/SRR479052.chr21_2.fastq \
	-o out/fastp/SRR479052.chr21_1.fastp.fastq \
	-O out/fastp/SRR479052.chr21_2.fastp.fastq \
	--cut_tail --cut_mean_quality 28 \
	-j log/fastp/SRR479052.json \
	-h log/fastp/SRR479052.html

fastp -i input/SRR479054.chr21_1.fastq \
	-I input/SRR479054.chr21_2.fastq \
	-o out/fastp/SRR479054.chr21_1.fastp.fastq \
	-O out/fastp/SRR479054.chr21_2.fastp.fastq \
	--cut_tail --cut_mean_quality 28 \
	-j log/fastp/SRR479054.json \
	-h log/fastp/SRR479054.html

# CONTROL DE CALIDAD POSTPROCESAMIENTO
fastqc -o out/fastqc \
	out/fastp/*fastp.fastq

# ANÁLISIS DE CONTAMINACIÓN
fastq_screen --get_genomes \
	--outdir res/
# ⚠️ Recuerda modificar la ruta al alineador en fastq_screen.conf

fastq_screen --conf res/FastQ_Screen_Genomes/fastq_screen.conf \
	--outdir out/fastq_screen/ \
	out/fastp/*

# INFORME INTEGRADO DE CALIDAD
multiqc . -o out/multiqc_report/

# ALINEAMIENTO CON HISAT2
hisat2-build --seed 123 \
	-p 2 \
	res/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
	res/Homo_sapiens.GRCh38.dna.chromosome.21

hisat2 --new-summary \
	--summary-file log/hisat2/SRR479052.hisat2.summary  \
	--rna-strandness R \
	--seed 123 \
	--phred33 \
	-p 2 \
	-k 1 \
	-x res/Homo_sapiens.GRCh38.dna.chromosome.21 \
	-1 out/fastp/SRR479052.chr21_1.fastp.fastq \
	-2 out/fastp/SRR479052.chr21_2.fastp.fastq \
	-S out/hisat2/SRR479052.sam

hisat2 --new-summary \
	--summary-file log/hisat2/SRR479054.hisat2.summary \
	--rna-strandness R \
	--seed 123 \
	--phred33 \
	-p 2 \
	-k 1 \
	-x res/Homo_sapiens.GRCh38.dna.chromosome.21 \
	-1 out/fastp/SRR479054.chr21_1.fastp.fastq \
	-2 out/fastp/SRR479054.chr21_2.fastp.fastq \
	-S out/hisat2/SRR479054.sam

# CONVERSIÓN A BAM
samtools view -bS out/hisat2/SRR479052.sam \
	> out/hisat2/SRR479052.bam

samtools view -bS out/hisat2/SRR479054.sam \
	> out/hisat2/SRR479054.bam

# ORDENAMIENTO POR POSICIÓN GENÓMICA
samtools sort out/hisat2/SRR479052.bam \
	-o out/hisat2/SRR479052.sorted.bam

samtools sort out/hisat2/SRR479054.bam \
	-o out/hisat2/SRR479054.sorted.bam

# INDEXACIÓN DE BAM
samtools index out/hisat2/SRR479052.sorted.bam

samtools index out/hisat2/SRR479054.sorted.bam

# GENERACIÓN DE MATRIZ DE CUENTAS CRUDAS
htseq-count --format=bam \
	--order=pos \
	--stranded=reverse \
	--mode=intersection-nonempty \
	--minaqual=10 \
	--type=exon \
	--idattr=gene_id \
	--additional-attr=gene_name \
	out/hisat2/SRR479052.sorted.bam \
	out/hisat2/SRR479054.sorted.bam \
	res/Homo_sapiens.GRCh38.109.chr21.gtf \
	> out/htseq/raw_counts.htseq

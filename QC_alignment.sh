# Quality control of the Illumina raw reads using fastqc and multiqc to aggregate the files
fastqc /raw_reads/*.fastq.gz -o qc_reports

multiqc qc_reports -o qc_reports/multiqc

# Alignment using STAR
STAR --runThreadN 8 --genomeDir star_index \
--readFilesIn raw_reads/GSM8287690_1.fastq.gz raw_reads/GSM8287690_2.fastq.gz \
--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix star_out/GSM8287690_

STAR --runThreadN 8 --genomeDir star_index \
--readFilesIn raw_reads/GSM8287691_1.fastq.gz raw_reads/GSM8287691_2.fastq.gz \
--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix star_out/GSM8287691_

STAR --runThreadN 8 --genomeDir star_index \
--readFilesIn raw_reads/GSM8287692_1.fastq.gz raw_reads/GSM8287692_2.fastq.gz \
--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix star_out/GSM8287692_

STAR --runThreadN 8 --genomeDir star_index \
--readFilesIn raw_reads/GSM8287693_1.fastq.gz raw_reads/GSM8287693_2.fastq.gz \
--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix star_out/GSM8287693_

STAR --runThreadN 8 --genomeDir star_index \
--readFilesIn raw_reads/GSM8287694_1.fastq.gz raw_reads/GSM8287694_2.fastq.gz \
--readFilesCommand zcat --outFileNamePrefix star_out/GSM8287694_

STAR --runThreadN 8 --genomeDir star_index \
--readFilesIn raw_reads/GSM8287695_1.fastq.gz raw_reads/GSM8287695_2.fastq.gz \
--readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --outFileNamePrefix star_out/GSM8287695_

# Counting of genes to create mytable_features file using featurecount
featureCounts -T 8 -a annotation.gtf -o mytable_features.txt star_out/*.bam
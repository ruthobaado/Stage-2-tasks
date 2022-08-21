#Reproduce variant calling workflow tutorial using the E. coli:
#download the reference genome:
mkdir -p data/ref_genome
curl -L -o data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
gunzip data/ref_genome/ecoli_rel606.fasta.gz
head data/ref_genome/ecoli_rel606.fasta

#to download a set of trimmed files 
curl -L -o sub.tar.gz https://ndownloader.figshare.com/files/14418248
tar xvf sub.tar.gz
mv sub/ ~/data/trimmed_fastq_small
ls data/trimmed_fastq_small
mkdir -p results/sam results/bam results/bcf results/vcf

#to index the reference genome:
bwa index data/ref_genome/ecoli_rel606.fasta
bwa mem data/ref_genome/ecoli_rel606.fasta data/trimmed_fastq_small/SRR2584866_1.trim.sub.fastq data/trimmed_fastq_small/SRR2584866_2.trim.sub.fastq > results/sam/SRR2584866.aligned.sam

#to convert the SAM file to BAM format:
samtools view -S -b results/sam/SRR2584866.aligned.sam > results/bam/SRR2584866.aligned.bam
 
#to sort the BAM file:
samtools sort -o results/bam/SRR2584866.aligned.sorted.bam results/bam/SRR2584866.aligned.bam 
samtools flagstat results/bam/SRR2584866.aligned.sorted.bam

#variant calling:
conda install -c bioconda bcftools
bcftools mpileup -O b -o results/bcf/SRR2584866_raw.bcf --no-reference -f data/ref_genome/ecoli_rel606.fasta results/bam/SRR2584866.aligned.sorted.bam 
bcftools call --ploidy 1 -m -v -o results/vcf/SRR2584866_variants.vcf results/bcf/SRR2584866_raw.bcf 
vcfutils.pl varFilter results/vcf/SRR2584866_variants.vcf  > results/vcf/SRR2584866_final_variants.vcf
less -S results/vcf/SRR2584866_final_variants.vcf
grep -v "#" results/vcf/SRR2584866_final_variants.vcf | wc -l

#to assess or visualize the alignment:
samtools index results/bam/SRR2584866.aligned.sorted.bam
samtools tview results/bam/SRR2584866.aligned.sorted.bam data/ref_genome/ecoli_rel606.fasta
samtools tview ~/results/bam/SRR2584866.aligned.sorted.bam ~/data/ref_genome/ecoli_rel606.fasta


#Reproduce variant calling workflow tutorial using the datasets from normal tissue and tumor tissue:
#to download datasets:
mkdir datasets
cd datasets
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz

#download the reference genome:
mkdir ref_genome
cd ref_genomels 
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz
gunzip datasets/ref_genome/hg19.chr5_12_17.fa.gz
mkdir -p results/sam results/bam results/bcf results/vcf

#to check for quality of the datasets using fastqc:
mkdir Output
fastqc datasets/*.fastq.gz -o Output

#to aggregate all the html files into a single file using multiqc:
multiqc Output/ 

#to trim adapters and poor quality reads using fastp:
touch trim.sh
nano trim.sh
#copy the code template from the Github, edit the codes and paste in trim.sh 
cp trim.sh datasets/
cd datasets/
bash trim.sh
ls
mv qc_reads trimmed_reads
ls

#to aligns relatively short sequences to a sequence base using BWA:
bwa index datasets/ref_genome/hg19.chr5_12_17.fa
repair.sh in1=trimmed_reads/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz in2=trimmed_reads/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz out1=SLGFSK-N_231335_r1_rep.fastq.gz out2=SLGFSK-N_231335_r2_rep.fastq.gz outsingle=single.fq
repair.sh in1=trimmed_reads/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz in2=trimmed_reads/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz out1=SLGFSK-T_231336_r1_rep.fastq.gz out2=SLGFSK-T_231336_r2_rep.fastq.gz outsingle=single.fq

bwa mem datasets/ref_genomes/hg19.chr5_12_17.fa SLGFSK-N_231335_r1_rep.fastq.gz SLGFSK-N_231335_r2_rep.fastq.gz > results/sam/SLGFSK-N_231335.aligned.sam
bwa mem datasets/ref_genomes/hg19.chr5_12_17.fa SLGFSK-T_231336_r1_rep.fastq.gz SLGFSK-T_231336_r2_rep.fastq.gz > results/sam/SLGFSK-T_231336.aligned.sam

#to convert the SAM file to BAM format:
samtools view -S -b results/sam/SLGFSK-N_231335.aligned.sam > results/bam/SLGFSK-N_231335.aligned.bam
samtools view -S -b results/sam/SLGFSK-T_231336.aligned.sam > results/bam/SLGFSK-T_231336.aligned.bam
 
#to sort the BAM file:
samtools sort -o results/bam/SLGFSK-N_231335.aligned.sorted.bam results/bam/SLGFSK-N_231335.aligned.bam 
samtools sort -o results/bam/SLGFSK-T_231336.aligned.sorted.bam results/bam/SLGFSK-T_231336.aligned.bam 

samtools flagstat results/bam/SLGFSK-N_231335.aligned.sorted.bam
samtools flagstat results/bam/SLGFSK-T_231336.aligned.sorted.bam

#variant calling using bcftools:
bcftools mpileup -O b -o results/bcf/SLGFSK-N_231335_raw.bcf --no-reference -f datasets/ref_genome/hg19.chr5_12_17.fa results/bam/SLGFSK-N_231335.aligned.sorted.bam 
bcftools mpileup -O b -o results/bcf/SLGFSK-T_231336_raw.bcf --no-reference -f datasets/ref_genome/hg19.chr5_12_17.fa results/bam/SLGFSK-T_231336.aligned.sorted.bam 

bcftools call --ploidy 1 -m -v -o results/vcf/SLGFSK-N_231335_variants.vcf results/bcf/SLGFSK-N_231335_raw.bcf 
bcftools call --ploidy 1 -m -v -o results/vcf/SLGFSK-T_231336_variants.vcf results/bcf/SLGFSK-T_231336_raw.bcf 

vcfutils.pl varFilter results/vcf/SLGFSK-N_231335_variants.vcf  > results/vcf/SLGFSK-N_231335_final_variants.vcf
vcfutils.pl varFilter results/vcf/SLGFSK-T_231336_variants.vcf  > results/vcf/SLGFSK-T_231336_final_variants.vcf

grep -v "#" results/vcf/SLGFSK-N_231335_final_variants.vcf | wc -l
grep -v "#" results/vcf/SLGFSK-T_231336_final_variants.vcf | wc -l

#to assess or visualize the alignment:
samtools index results/bam/SLGFSK-N_231335.aligned.sorted.bam
samtools tview results/bam/SLGFSK-N_231335.aligned.sorted.bam datasets/ref_genome/hg19.chr5_12_17.fa


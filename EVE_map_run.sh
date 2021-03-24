!# /bin/bash

cd $1
genome=$2
GT_path=""
fastq_path=""
sam_path=""
bam_path=""
sorted_bam_path="" 
bed_path=""
tab_path=""


for i in *.fna; do bowtie2-build  --large-index -f${genome} ${genome}.bwt; done #smalls
for i in *.fna; do hisat2-build --large-index -p 25  ${genome} ${genome}.hst; done  #longs

for cond1 in F_G F_S M_G M_S; do 
cond2=$(echo $cond1| sed 's/_//')
bowtie2 --threads 25 -x ${GT_path}/${genome}_GT.fna.bwt -U ${fastq_path}/${genome}_${cond1}_sRNA_f_t.fq -S ${sam_path}/${genome}_${cond2}_small.SAM --no-unal 
done 

for cond1 in F_G F_S M_G M_S; do 
cond2=$(echo $cond1| sed 's/_//')
/mnt/65To/bin/hisat2-2.1.0/hisat2 -p 25 -x ${GT_path}/${genome}_GT.fna.hst -q -1 ${fastq_path}/${genome}_${cond1}_RNAseq_f_t.1.fq -2 ${fastq_path}/${genome}_${cond1}_RNAseq_f_t.2.fq --rna-strandness FR -S ${sam_path}/${genome}_${cond2}_long.SAM --no-unal 
done 

for l in small long; do  
  for cond1 in F_G F_S M_G M_S; do 
    /mnt/65To/bin/samtools/bin/samtools view -h -bS ${sam_path}/${genome}_${cond2}_${l}.SAM > ${bam_path}/${genome}_${cond2}_${l}.bam
    /mnt/65To/bin/samtools/bin/samtools sort ${bam_path}/${genome}_${cond2}_${l}.bam ${sorted_bam_path}/${genome}_${cond2}_${l}.bam
    /mnt/65To/bin/samtools/bin/samtools index ${sorted_bam_path}/${genome}_${cond2}_${l}.bam 
  done
done

for cond1 in FG FS MG MS; do 
python3 eve_mapping.py -f 80000 -g  -pb ${bed_path}/${genome}.bed -sB ${sorted_bam_path}/${genome}_${cond2}_small.bam -lB ${sorted_bam_path}/${genome}_${cond2}_long.bam -o ${tab_path}/${genomes}_${cond1}.tab
done


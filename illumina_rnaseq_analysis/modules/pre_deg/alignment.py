

#Reference: STAR: https://groups.google.com/g/rna-star/c/Cpsf-_rLK9I ;
#                 https://www.reneshbedre.com/blog/star-aligner-twopass-mode.html

#Good tutorials ; https://nbisweden.github.io/workshop-RNAseq/2011/lab_smallrna.html 
#                 https://faqs.lexogen.com/faq/how-can-i-analyze-my-small-rna-sequencing-data

import os

#STAR Alignment
def staraligner_run(sample_name, reads, threads, genomedir):
    bam_dir="Results/bamfiles/"
    os.system(f"mkdir -p {bam_dir}")
    os.system(f"mkdir -p Results/log_files/")
    if len(reads) > 1:
        os.system(f'''STAR --genomeDir {genomedir} \
                   --readFilesIn {reads[0]} {reads[1]} \
                   --readFilesCommand zcat \
                   --outSAMunmapped Within \
                   --outFileNamePrefix {bam_dir}/{sample_name}.\
                   --runThreadN {threads} > Results/log_files/alignment_star_{sample_name}.log 2>&1''')
    else:
        os.system(f'''STAR --genomeDir {genomedir} \
                   --readFilesIn {reads[0]} \
                   --readFilesCommand zcat \
                   --outSAMunmapped Within \
                   --outFileNamePrefix {bam_dir}/{sample_name}.\
                   --runThreadN {threads} > Results/log_files/alignment_star_{sample_name}.log 2>&1''')        


##Alignment using hisat2 ; https://daehwankimlab.github.io/hisat2/ 
def hisat2aligner_run(sample_name, reads, threads, fasta,strand="FR"):
    bam_dir="Results/bamfiles/"
    os.system(f"mkdir -p {bam_dir}")
    os.system(f"mkdir -p Results/log_files/")
    ID = f"HTYYCBBXX.1.{sample_name}"
    SM = sample_name
    LB = "Lib-1"
    PU = "HTYYCBBXX.1"
    PL = "ILLUMINA"
    if len(reads) > 1:
        os.system(f"""hisat2 \
                        -p {threads} \
                        --rg-id={ID} --rg SM:{SM} --rg LB:{LB} --rg PL:{PL} --rg PU:{PU} \
                        -x {fasta} \
                        --dta \
                        --rna-strandness {strand} \
                        -1 {reads[0]} \
                        -2 {reads[1]} \
                         2> Results/log_files/alignment_hisat2_{sample_name}.log \
                        | samtools sort -O BAM \
                        | tee Results/bamfiles/{sample_name}.bam \
                        | samtools index - Results/bamfiles/{sample_name}.bai """)
    else:
        os.system(f"""hisat2 \
                        -p {threads} \
                        --rg-id={ID} --rg SM:{SM} --rg LB:{LB} --rg PL:{PL} --rg PU:{PU} \
                        -x {fasta} \
                        --dta \
                        --rna-strandness {strand} \
                        -1 {reads[0]} \
                         2> Results/log_files/alignment_hisat2_{sample_name}.log \
                        | samtools sort -O BAM \
                        | tee Results/bamfiles/{sample_name}.bam \
                        | samtools index - Results/bamfiles/{sample_name}.bai""")


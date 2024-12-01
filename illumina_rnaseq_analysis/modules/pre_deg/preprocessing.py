import os
from typing import List
import pandas as pd
import glob 

## Fastqc command https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
def fastqc_run(sample_name:str, reads:List[str], threads:int,output_folder:str):
    fastqc_dir=f"Results/QC/fastq_qc/{output_folder}/"
    fastqc_dir = f"{fastqc_dir}/{sample_name}"
    os.system(f"mkdir -p {fastqc_dir}")
    os.system(f"mkdir -p Results/log_files/")
    if len(reads) > 1:
        os.system(f"fastqc -t {threads} -o {fastqc_dir} {reads[0]} {reads[1]} > Results/log_files/fastq-qc_fastqc_{sample_name}.log 2>&1 ")
    else:
        os.system(f"fastqc -t {threads} -o {fastqc_dir} {reads[0]}  > Results/log_files/fastq-qc_fastqc_{sample_name}.log 2>&1 ")


## Falco command https://github.com/smithlabcode/falco
def falco_run(sample_name:str, reads,output_folder):
    print(f"now sample {sample_name} running in falcon")
    falco_dir=f"Results/QC/fastq_qc/{output_folder}"
    falco_dir = f"{falco_dir}/{sample_name}"
    os.system(f"mkdir -p {falco_dir}")
    os.system(f"mkdir -p Results/log_files/")
    if len(reads) > 1:
        os.system(f"falco --outdir {falco_dir} \
                          -report-filename {falco_dir}/{sample_name}_R1.html \
                          -summary-filename {falco_dir}/{sample_name}_R1.txt \
                          -data-filename {falco_dir}/{sample_name}_R1_fastqc.txt {reads[0]} > Results/log_files/fastq-qc_falco_{sample_name}.log 2>&1")
        os.system(f"falco --outdir {falco_dir} \
                          -report-filename {falco_dir}/{sample_name}_R2.html \
                          -summary-filename {falco_dir}/{sample_name}_R2.txt \
                          -data-filename {falco_dir}/{sample_name}_R2_fastqc.txt {reads[1]} > Results/log_files/fastq-qc_falco_{sample_name}.log 2>&1")
    else:
        os.system(f"falco --outdir {falco_dir} \
                          -report-filename {falco_dir}/{sample_name}_R1.html \
                          -summary-filename {falco_dir}/{sample_name}_R1.txt \
                          -data-filename {falco_dir}/{sample_name}_R1_fastqc.txt {reads[0]} > Results/log_files/fastq-qc_falco_{sample_name}.log 2>&1")


##Fastp command ; https://github.com/OpenGene/fastp
def fastp_run(sample_name, reads, threads,qc_cutoff):
    print(f"now sample {sample_name} running in fastp")
    output_dir = f"Results/cleaned_fastQ_files/"
    output_dir2 = f"Results/QC/fastq_trimming_stat/"
    os.system(f"mkdir -p {output_dir}")
    os.system(f"mkdir -p {output_dir2} Results/log_files/")

    if len(reads) > 1:
        os.system(f''' fastp \
                    --in1 {reads[0]} --in2 {reads[1]} \
                    --out1 {output_dir}/{sample_name}_fastp_R1.gz --out2 {output_dir}/{sample_name}_fastp_R2.gz \
                    --unpaired1 {output_dir}/{sample_name}_unpaired_R1.gz --unpaired2 {output_dir}/{sample_name}_unpaired_R2.gz \
                    --cut_front cut_front_window_size 3 --cut_front_mean_quality {qc_cutoff} \
                    --cut_tail cut_tail_window_size 3 --cut_tail_mean_quality {qc_cutoff} \
                    --detect_adapter_for_pe --compression 4 \
                    --qualified_quality_phred {qc_cutoff} --length_required 30  \
                    --trim_poly_x  --correction --thread {threads} \
                    --json {output_dir2}{sample_name}_fastp.json \
                    --html {output_dir2}{sample_name}_fastp.html > Results/log_files/fastq-qc_fastp_{sample_name}.log 2>&1''')
    else:
        os.system(f'''fastp \
                    --in1 {reads[0]} \fastp
                    --out1 {output_dir}/{sample_name}_trimmed.fq.gz \
                    --cut_front cut_front_window_size 3 --cut_front_mean_quality {qc_cutoff} \
                    --cut_tail cut_tail_window_size 3 --cut_tail_mean_quality {qc_cutoff} \
                    --qualified_quality_phred {qc_cutoff} --length_required 30 \
                    --compression 4 \
                    --trim_poly_x  --thread {threads} \
                    --json {output_dir2}{sample_name}_fastp.json \
                    --html {output_dir2}{sample_name}_fastp.html > Results/log_files/fastq-qc_fastp_{sample_name}.log 2>&1''')


##trim_galore command ; https://github.com/FelixKrueger/TrimGalore/tree/master
def trim_galore_run(sample_name, reads, threads,qc_cutoff):
    trimmed_fq_dir = f"Results/cleaned_fastq_files/"
    trimming_stat_dir = f"Results/QC/fastq_trimming_stats/"
    os.system(f"mkdir -p {trimmed_fq_dir}")
    os.system(f"mkdir -p {trimming_stat_dir}")
    os.system(f"mkdir -p Results/log_files/")
    
    if len(reads) > 1:
        os.system(f'''trim_galore \
                        --quality {qc_cutoff} \
                        --gzip \
                        --length 16  \
                        --basename {sample_name} \
                        --cores {threads} \
                        --fastqc_args "--outdir {trimming_stat_dir}" 
                        --output_dir {trimmed_fq_dir} \
                        --polyA \
                        --paired \
                        {reads[0]} {reads[1]} > Results/log_files/fastq-qc_trim_galore_{sample_name}.log 2>&1''')
    else:
        os.system(f'''trim_galore \
                        --quality {qc_cutoff} \
                        --gzip \
                        --length 16  \
                        --basename {sample_name} \
                        --cores {threads} \
                        --fastqc_args "--outdir {trimming_stat_dir}" \
                        --output_dir {trimmed_fq_dir} \
                        {reads[0]} > Results/log_files/fastq-qc_trim_galore_{sample_name}.log 2>&1''')


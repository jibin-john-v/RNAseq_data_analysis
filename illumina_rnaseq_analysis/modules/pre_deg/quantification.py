
import os


def featurecount_run(gtf,feature,attribute,strand,threads,out_name):
    os.system(' mkdir -p Results/quantification')
    os.system(f'''featureCounts \
                        -t {feature} \
                        -g {attribute} \
                        -O \
                        -s {strand} \
                        -M \
                        -a {gtf} \
                        -T {threads} \
                        -o Results/quantification/{out_name} \
                        Results/Bam_Files/*.bam > \
                        Results/log_files/quantification_featurecout.log 2>&1''')

def htseqcount_run(gtf,feature,attribute,strand,threads,bam_file,sample_name):
    os.system(' mkdir -p Results/Quantification')
    os.system(f'''htseq-count \
                        --format bam \
                        --order pos \
                        --mode intersection-strict \
                        --stranded {strand} \
                        --minaqual 1 \
                        --type {feature} \
                        --idattr {attribute} \
                        --nprocesses {threads} \
                        {bam_file} \
                        {gtf}> Results/quantification/{sample_name}_gene.tsv > \
                        Results/log_files/quantification_htseq_{sample_name}.log 2>&1''')

#stringtie quantification ; https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
def quant_stringtie_run(sample_name:str,gtf,threads,bam_file):
    stringie_dir="Results/quantification/"
    os.system( f"mkdir -p {stringie_dir}")
    os.system(f'stringtie -p {threads} -G {gtf} -o {stringie_dir}/{sample_name}.gtf -l {sample_name} {bam_file}')



#stringtie quantification ; https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
def merge_stringtie_run(gtf,threads,merge_list):
    stringie_dir=f"Results/quantification/"
    os.system(f"mkdir -p  {stringie_dir}")
    os.system(f"stringtie --merge -p {threads} -G {gtf} -o {stringie_dir}/stringtie_merged.gtf {merge_list}")


#gffcompare ;https://ccb.jhu.edu/software/stringtie/gffcompare.shtml
def gffcompare_run(ref_gtf,merge_gtf):
    gffcompare_dir="Results/gffcompare_out/"
    os.system(f"mkdir -p  {gffcompare_dir}")
    os.system(f"gffcompare -r {ref_gtf} -G -o {gffcompare_dir}/merged {merge_gtf}")



#stringtie quantification ; https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual
def quant_stringtie_run2(sample_name:str,merged_gtf,threads,bam_file):
    stringie_dir=f"Results/ballgown/{sample_name}"
    os.system(f"mkdir -p {stringie_dir}")
    os.system(f'stringtie -e -B -p {threads} -G {merged_gtf} -o {stringie_dir}/{sample_name}_2nd.gtf -l {sample_name} {bam_file}')



def quant_salmon_run(salmon_index,sample_name, reads):
    salmon_dir=f"mkdir -p Results/quantification/"
    os.system( f"{salmon_dir}")

    if len(reads) > 1:
        os.system( f"salmon quant \
         -i {salmon_index} -l A \
         -1 {reads[0]} \
         -2 {reads[1]} \
         -p 8 --validateMappings \
         -o Results/Quantification/{sample_name}_quant")
    else:
        os.system( f"salmon quant \
         -i {salmon_index} -l A \
         -1 {reads[0]} \
         -p 8 --validateMappings \
         -o Results/Quantification/{sample_name}_quant")
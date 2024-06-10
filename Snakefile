# Defining files

work_dir = "/mnt/c/Users/sunil/Downloads/NGS_example_dataset/paeruginosa-reads/shell/script/"
R1 = "/mnt/c/Users/sunil/Downloads/NGS_example_dataset/paeruginosa-reads/shell/script/SRR396636.sra_1.fastq"
R2 = "/mnt/c/Users/sunil/Downloads/NGS_example_dataset/paeruginosa-reads/shell/script/SRR396636.sra_2.fastq"
Analysis_summary = "/mnt/c/Users/sunil/Downloads/NGS_example_dataset/paeruginosa-reads/shell/script/Analysis_summary.txt"

# Learning about basename
'''
basename can be used to hide complete path
'''
rule basename:
    shell: """
    echo "$(basename ${Analysis_summary})"
"""

# Learning about dirname
rule dirname:
    shell: """
    echo "$(dirname ${Analysis_summary})"
"""

# 0. Determining the number of reads in fastq file
rule fastq_reads_count:
    input:
        {R1}, {R2}

    shell: """
    echo "Number of reads in {R1} are: $(cat {R1} | grep @SRR | wc -l)" | tee -a {Analysis_summary}
    echo "Number of reads in {R2} are: $(cat {R2} | grep @SRR | wc -l)" | tee -a {Analysis_summary}
"""

# 1. Running fastqc
rule fastqc:
    input:
        {R1}, {R2}

    shell: """
    mkdir {work_dir}fastqc_results
    fastqc {input} -o {work_dir}fastqc_results
    echo "Fastqc is completed!"
"""

# 1.1
rule fastqc2:
    input:
        "{fastqc}_1.fastq"
        "{fastqc}_2.fastq"
    output:
        "{fastqc}_1.fastqc.html"
        "{fastqc}_2.fastqc.html"
    
    shell: '''
    mkdir fastqc_reports
    fastqc {input} -o fastqc_reports/ {output}
'''

    



# 2. Running trim_galore
rule trim_galore:
    input:
        {R1}, {R2}

    shell: """
    mkdir {work_dir}trim_galore {work_dir}trim_galore_fastqc
    trim_galore --paired {input} --output_dir {work_dir}trim_galore \
    --fastqc_args "--outdir {work_dir}trim_galore_fastqc"
    echo "Trimmomatic is completed!"
"""

import os
import glob
R1_trim = glob.glob(os.path.join(work_dir, "trim_galore", "*1.fq"))
R2_trim = glob.glob(os.path.join(work_dir, "trim_galore", "*2.fq"))

# 2.1. Determining the number of reads in fastq file after trimming
rule fastq_reads_count2:
    shell: """
    echo "{R1_trim}"
    echo "Number of reads in {R1_trim} are: $(cat {R1_trim} | grep @SRR | wc -l)" | tee -a {Analysis_summary}
    echo "Number of reads in {R2_trim} are: $(cat {R2_trim} | grep @SRR | wc -l)" | tee -a {Analysis_summary}
"""

# 3. Alignment

# 3.1 Generating index for the reference genome
rule reference_genome_index:
    input:
        {work_dir} reference_genome/GCA_000006765.1_ASM676v1_genomic.fna

    shell: """
    bwa index -p {input}
"""

# ## 3.2 Alignment of reads to the reference genome
# rule alignment:
# input:
#     "dir/fastq_files/{sample}.fastq",
#     "dir/reference_genome/{sample}.fna"
# output:
#     "dir/alignment/{sample}.sam"
# shell:
#     """
#     bwa mem {input} | samtools view -Sb - > {output}
#     """

# ## 3.3 Generating alignment index file
# rule bam.bai
# input:
#     "dir/alignement_files/{sample}.bam"
# output:
#     "dir/alignment_files/{sample}.bam.bai"
# shell:
#     '''
#     samtools index {input}
#     '''
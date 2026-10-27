################################################################################
########################## Ref download and generation #########################
################################################################################

from utils import *

rule curl_fa:
    threads: 1
    resources:
        mem_mb = 4096,
        partition = config['queuing']['standard']['partition'],
        account = config['queuing']['standard']['account'],
        gres = config['queuing']['standard']['gres'],
        time = '4:00:00'
    params:
        link = lambda wc: get_fa_link(wc, config),
    output:
        zip = temporary("{genotype}.zip")
    shell:
        """
        if [ "{wildcards.genotype}" == "B6J" ]; then
            wget -O {wildcards.genotype}.fa.gz {params.link}
            mkdir -p {wildcards.genotype}/ncbi_dataset/data/temp/
            gunzip {wildcards.genotype}.fa.gz
            mv {wildcards.genotype}.fa {wildcards.genotype}/ncbi_dataset/data/temp/temp.fna
            zip -r {wildcards.genotype}.zip {wildcards.genotype}
        else
            curl -OJX GET "{params.link}{output.zip}"
        fi
        """

rule fa_ref_fmt:
    input:
        zip = "{genotype}.zip"
    threads: 2
    resources:
        mem_mb = 4096,
        partition = config['queuing']['standard']['partition'],
        account = config['queuing']['standard']['account'],
        gres = config['queuing']['standard']['gres'],
        time = '1:00:00'
    output:
        fa = "ref/genomes/{genotype}.fa.gz"
    shell:
        """
        unzip {input.zip} -d {wildcards.genotype}
        gzip -cvf {wildcards.genotype}/ncbi_dataset/data/*/*fna > {output.fa}
        rm -r {wildcards.genotype}
        """

rule dl:
    threads: 1
    resources:
        mem_mb = 4096,
        partition = config['queuing']['standard']['partition'],
        account = config['queuing']['standard']['account'],
        gres = config['queuing']['standard']['gres'],
        time = '4:00:00'
    shell:
        "wget -O {output.out} {params.link}"

use rule dl as dl_annot with:
    params:
        link = config['ref']['annot_link'],
    output:
        out = 'ref/gencode.v{}.gtf.gz'.format(config['ref']['annot_version']),

use rule dl as dl_fa with:
    params:
        link = config['ref']['fa_link'],
    output:
        out = 'ref/{}.fa.gz'.format(config['ref']['fa_version']),

################################################################################
################################ Symlink fastqs ################################
################################################################################

rule symlink_fastq_r1:
    params:
        fastq = lambda wc:get_df_info(wc, df, 'fastq')
    threads: 1
    resources:
        mem_mb = 4096,
        partition = config['queuing']['standard']['partition'],
        account = config['queuing']['standard']['account'],
        gres = config['queuing']['standard']['gres'],
        time = '1:00:00'
    output:
        fastq = "fastq/{plate}_{subpool}_{run}_{lane}_R1_001.fastq.gz",
    shell:
        """
        ln -s {params.fastq} {output.fastq}
        """

rule symlink_fastq_r2:
    params:
        fastq = lambda wc:get_df_info(wc, df, 'fastq_r2')
    threads: 1
    resources:
        mem_mb = 4096,
        partition = config['queuing']['standard']['partition'],
        account = config['queuing']['standard']['account'],
        gres = config['queuing']['standard']['gres'],
        time = '1:00:00'
    output:
        fastq = "fastq/{plate}_{subpool}_{run}_{lane}_R2_001.fastq.gz",
    shell:
        """
        ln -s {params.fastq} {output.fastq}
        """




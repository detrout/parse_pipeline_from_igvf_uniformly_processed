import pandas as pd

from utils import *

BYTES_TO_MB = 1048576   # divide by this to convert from bytes to megabytes

######## Everything should be configured via the configfile ########
configfile: 'configs/config.yml'

config_tsv = config['analysis']['analysis_config']
sample_csv = config['analysis']['sample_metadata']

kit = config['analysis']['kit']
chemistry = config['analysis']['chemistry']

# read in config / analysis spec
df = parse_config(config_tsv)
bc_df = get_bc1_matches(kit, chemistry)
sample_df = parse_sample_df(sample_csv)

r1_RT_filename = f"ref/r1_RT_replace_{kit}_{chemistry}.txt"
create_r1_RT_replace(kit, chemistry, r1_RT_filename)

r2r2r3_filename = f"ref/r1r2r3_{kit}_{chemistry}.txt"
create_r1r2r3_file(kit, chemistry, r2r2r3_filename)

mult_genotype_1s = sample_df.loc[sample_df.mult_genotype_1.notnull(), 'mult_genotype_1'].unique().tolist()
mult_genotype_2s = sample_df.loc[sample_df.mult_genotype_2.notnull(), 'mult_genotype_2'].unique().tolist()

wildcard_constraints:
    plate='|'.join([re.escape(x) for x in df.plate.tolist()]),
    subpool='|'.join([re.escape(x) for x in df.subpool.tolist()]),
    lane='|'.join([re.escape(x) for x in df.lane.tolist()]),
    sample='|'.join([re.escape(x) for x in sample_df.Mouse_Tissue_ID.tolist()]),
    tissue='|'.join([re.escape(x) for x in sample_df.Tissue.tolist()]),
    genotype='|'.join([re.escape(x) for x in sample_df.Genotype.unique().tolist()]),
    mult_genotype_1='|'.join([re.escape(x) for x in mult_genotype_1s]),
    mult_genotype_2='|'.join([re.escape(x) for x in mult_genotype_2s])


# these patterns needed to be defined before the all rule so I can't
# use the rules.${rule}.output.${variable} convention because they're
# not defined yet.  but I still want the patterns in the snakemake and
# to not repeat myself.
KLUE_IND_OUTPUT = "ref/klue/{mult_genotype_1}_{mult_genotype_2}.idx"
KLUE_GET_GENOTYPE_COUNTS_OUTPUT = "pipeline_output_klue/{plate}/{subpool}/genotype_counts.tsv"
MAKE_TISSUE_ADATA_OUTPUT = "{plate}/{tissue}/adata.h5ad"
COPY_CELLBENDER_METRICS = "pipeline_output_cellbender/{plate}/{subpool}/adata_denoised_metrics.csv"

rule all:
    input:
        expand(KLUE_IND_OUTPUT,
               zip,
               mult_genotype_1=[g for g in mult_genotype_1s if g in get_founder_genotypes()],
               mult_genotype_2=[g for g in mult_genotype_2s if g in get_founder_genotypes()]),
        expand(KLUE_GET_GENOTYPE_COUNTS_OUTPUT,
                zip,
                plate=df.plate.tolist(),
                subpool=df.subpool.tolist()),
        expand(MAKE_TISSUE_ADATA_OUTPUT,
                plate=df.plate.tolist(),
                tissue=get_subset_tissues(df, sample_df)),
#        expand(COPY_CELLBENDER_METRICS, #rules.copy_cellbender_metrics.output.metrics_copy,
#               plate=df.plate.tolist(),
#               subpool=df.subpool.tolist())

################################################################################
########################## Ref download and generation #########################
################################################################################

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

rule kallisto_ind:
    input:
        annot = rules.dl_annot.output.out,
        fa = rules.dl_fa.output.out
    threads: 12
    resources:
        mem_mb = 32768,
        partition = config['queuing']['highmem']['partition'],
        account = config['queuing']['highmem']['account'],
        gres = config['queuing']['highmem']['gres'],
        time = '4:00:00'
    conda:
        'envs/kb_env.yaml' 
    output:
        t2g = "ref/kallisto/c57bl6j.t2g",
        ind = "ref/kallisto/c57bl6j.idx",
        fa = "ref/kallisto/c57bl6j.cdna.fasta",
        na = "ref/kallisto/c57bl6j.nascent.fasta",
        c1 = "ref/kallisto/c57bl6j.c1",
        c2 = "ref/kallisto/c57bl6j.c2",
    shell:
        """        
        kb ref \
            -t {threads} \
            --workflow=nac \
            -i {output.ind} \
            -g {output.t2g} \
            -c1 {output.c1} \
            -c2 {output.c2} \
            -f1 {output.fa} \
            -f2 {output.na} \
            {input.fa} \
            {input.annot}
        """

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



################################################################################
##################################### klue #####################################
################################################################################

rule klue_fa:
    input:
        fa_g1 = lambda wc: expand(rules.fa_ref_fmt.output.fa,
                    genotype=wc.mult_genotype_1)[0],
        fa_g2 = lambda wc: expand(rules.fa_ref_fmt.output.fa,
                    genotype=wc.mult_genotype_2)[0]
    params:
        temp = "ref/klue/{mult_genotype_1}_{mult_genotype_2}"
    conda:
        'envs/kb_env.yaml'
    threads: 32
    resources:
        mem_mb = 32768,
        partition = config['queuing']['highmem']['partition'],
        account = config['queuing']['highmem']['account'],
        gres = config['queuing']['highmem']['gres'],
        time = '4:00:00'
    output:
        fa = "ref/klue/{mult_genotype_1}_{mult_genotype_2}.fa",
        t2g = "ref/klue/{mult_genotype_1}_{mult_genotype_2}.t2g",
    shell:
        """
        klue distinguish \
            -o {output.fa} \
            -g {output.t2g} \
            -t {threads} \
            -r 61 \
            --all-but-one \
            --tmp {params.temp}/ \
            {input.fa_g1} {input.fa_g2}
        """

rule klue_ind:
    input:
        fa = rules.klue_fa.output.fa
    params:
        temp = rules.klue_fa.params.temp
    conda:
        'envs/kb_env.yaml'  
    threads: 32
    resources:
        mem_mb = 32768,
        partition = config['queuing']['highmem']['partition'],
        account = config['queuing']['highmem']['account'],
        gres = config['queuing']['highmem']['gres'],
        time = '4:00:00'
    output:
        ind = KLUE_IND_OUTPUT
    shell:
        """
        kb ref \
            --workflow=custom \
            --distinguish \
            -i {output.ind} \
            -t {threads} \
            --tmp {params.temp}/ \
            {input.fa}
        """

rule klue:
    input:
        r1_fastq = lambda wc:get_subpool_fastqs(wc, df, config, how='list', read='R1'),
        r2_fastq = lambda wc:get_subpool_fastqs(wc, df, config, how='list', read='R2'),
        t2g = rules.klue_fa.output.t2g,
        ind = rules.klue_ind.output.ind
    params:
        fastq_str = lambda wc:get_subpool_fastqs(wc, df, config, how='str'),
        odir = lambda output: output.cgb.split('counts_unfiltered_modified')[0],
    threads: 24
    resources:
        mem_mb = 32768,
        partition = config['queuing']['highmem']['partition'],
        account = config['queuing']['highmem']['account'],
        gres = config['queuing']['highmem']['gres'],
        time = '4:00:00'
    conda:
        'envs/kb_env.yaml'  
    output:
        cgb="klue/{plate}/{subpool}/{mult_genotype_1}_{mult_genotype_2}/counts_unfiltered_modified/cells_x_genes.barcodes.txt",
        cggn="klue/{plate}/{subpool}/{mult_genotype_1}_{mult_genotype_2}/counts_unfiltered_modified/cells_x_genes.genes.names.txt",
        cgg="klue/{plate}/{subpool}/{mult_genotype_1}_{mult_genotype_2}/counts_unfiltered_modified/cells_x_genes.genes.txt",
        adata="klue/{plate}/{subpool}/{mult_genotype_1}_{mult_genotype_2}/counts_unfiltered_modified/adata.h5ad",
        mtx="klue/{plate}/{subpool}/{mult_genotype_1}_{mult_genotype_2}/counts_unfiltered_modified/cells_x_genes.mtx",
        bus=temporary("klue/{plate}/{subpool}/{mult_genotype_1}_{mult_genotype_2}/output.bus"),
        bus_modified_unfilt=temporary("klue/{plate}/{subpool}/{mult_genotype_1}_{mult_genotype_2}/output_modified.unfiltered.bus"),
        bus_unfilt=temporary("klue/{plate}/{subpool}/{mult_genotype_1}_{mult_genotype_2}/output.unfiltered.bus"),
    shell:
        """
        kb count \
            --h5ad \
            --gene-names \
            --strand=forward \
            --mm \
            -r {r1_RT_filename} \
            -w {r2r2r3_filename} \
            --workflow=standard \
            -g {input.t2g} \
            -x SPLIT-SEQ \
            -i {input.ind} \
            -t {threads} \
            -o {params.odir} \
            --tmp {params.odir}/temp/ \
            --verbose \
            {params.fastq_str}
        """

# add metadata
rule make_adata_klue:
    threads: 4
    resources:
        mem_mb = lambda wildcards, input, attempt: (input.size // BYTES_TO_MB) * 2,
        partition = config['queuing']['highmem']['partition'],
        account = config['queuing']['highmem']['account'],
        gres = config['queuing']['highmem']['gres'],
        time = '1:00:00'
    input:
        mtx = rules.klue.output.mtx,
        cgb = rules.klue.output.cgb,
        cggn = rules.klue.output.cggn,
        cgg = rules.klue.output.cgg,
    output:
        adata = "pipeline_output_klue/{plate}/{subpool}/{mult_genotype_1}_{mult_genotype_2}/unfiltered_adata.h5ad",
    run:
        add_meta_klue(input.mtx,
                        input.cgb,
                        input.cggn,
                        input.cgg,
                        wildcards,
                        bc_df,
                        kit,
                        chemistry,
                        sample_df,
                        output.adata)

# get counts for each cell across all genotypes
rule klue_get_genotype_counts:
    input:
        adatas = lambda wc:get_subpool_adatas(df, sample_df, wc, rules.make_adata_klue.output.adata)
    threads: 2
    resources:
        mem_mb = lambda wildcards, input, attempt: (input.size // BYTES_TO_MB) * 3,
        partition = config['queuing']['standard']['partition'],
        account = config['queuing']['standard']['account'],
        gres = config['queuing']['standard']['gres'],
        time = '1:00:00'
    output:
        ofile = KLUE_GET_GENOTYPE_COUNTS_OUTPUT
    run:
        get_genotype_counts(input.adatas, output.ofile)


################################################################################
####################### Use precomputed unfilt_adata ###########################
################################################################################

rule link_existing_unfilt_adata:
    """
    Symlink the pre-generated unfilt_adata (.h5ad) whose path lives
    in the `igvf_adata` column of your igvf_b01_config.tsv.
    """
    input:
        # df was created at top via df = parse_config(config_tsv)
        h5ad = lambda wc: df.loc[
            (df.plate == wc.plate) &
            (df.subpool == wc.subpool),
            'igvf_adata'
        ].iloc[0]
    threads: 1
    resources:
        mem_mb = 4096,
        partition = config['queuing']['standard']['partition'],
        account = config['queuing']['standard']['account'],
        gres = config['queuing']['standard']['gres'],
        time = '1:00:00'

    output:
        unfilt_adata = "pipeline_output_kallisto/{plate}/{subpool}/unfiltered_adata.h5ad",
    shell:
        """
        mkdir -p $(dirname {output.unfilt_adata})
        ln -s {input.h5ad} {output.unfilt_adata}
        """



################################################################################
################################## cellbender ##################################
################################################################################


rule cellbender:
    input:
        raw_adata = rules.link_existing_unfilt_adata.output.unfilt_adata
    params:
        total_drops    = lambda wc: int(df.loc[df.subpool.eq(wc.subpool), 'droplets_included'].iloc[0]),
        learning_rate  = lambda wc: float(df.loc[df.subpool.eq(wc.subpool), 'learning_rate'].iloc[0]),
        expected_cells = lambda wc: int(df.loc[df.subpool.eq(wc.subpool), 'expected_cells'].iloc[0]),
        use_cuda = '--cuda' if config['queuing']['gpu']['has_cuda'] else '',
    threads: 1
    resources:
        mem_mb  = lambda wildcards, input, attempt: (input.size // BYTES_TO_MB) * 6,
        partition = config["queuing"]["gpu"]["partition"],
        account   = config["queuing"]["gpu"]["account"],
        gres      = config["queuing"]["gpu"]["gres"],
        time      = "8:00:00"
    output:
        filt_h5   = "cellbender/{plate}/{subpool}/adata_denoised_filtered.h5",
        unfilt_h5 = "cellbender/{plate}/{subpool}/adata_denoised.h5",
        metrics   = "cellbender/{plate}/{subpool}/adata_denoised_metrics.csv",
    log:
        "logs/cellbender/{plate}/{subpool}.log"
    conda:
        'envs/kb_env.yaml'
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.unfilt_h5})" "$(dirname {log})"

        set +u
        set -u

        INP="$(python - <<'PY'
from pathlib import Path
print(Path("{input.raw_adata}").resolve())
PY
)"

        RUN_DIR="$(dirname {output.unfilt_h5})"
        OUT_BASE="$(basename {output.unfilt_h5})"
        OUT_BASE="${{OUT_BASE%.h5}}"   # <-- escape braces so Snakemake doesn't parse them

        (
          cd "$RUN_DIR"
          cellbender remove-background \
            --input "$INP" \
            --output "$OUT_BASE" \
            --total-droplets-included {params.total_drops} \
            --learning-rate {params.learning_rate} \
            --expected-cells {params.expected_cells} \
            {params.use_cuda}
        ) 2>&1 | tee "{log}"

        python - <<'PY'
import csv, os
p = "{output.metrics}"
os.makedirs(os.path.dirname(p), exist_ok=True)
with open(p, "w", newline="") as f:
    w = csv.writer(f); w.writerow(["plate","subpool","total_drops","learning_rate","expected_cells"])
    w.writerow([{wildcards.plate!r}, {wildcards.subpool!r}, {params.total_drops}, {params.learning_rate}, {params.expected_cells}])
PY
        """



#        """
#        mkdir -p $(dirname {output.filt_h5})
#        cd $(dirname {output.filt_h5})

#        source ~/miniconda3/bin/activate cellbender
#        cellbender remove-background \
#          --input {input.raw_adata} \
#          --output {output.unfilt_h5} \
#          --total-droplets-included {params.total_drops} \
#          --learning-rate {params.learning_rate} \
#          --expected-cells {params.expected_cells} \
#          --cuda
#        """


#rule copy_cellbender_metrics:
#    input:
#        metrics = rules.cellbender.output.metrics
#    threads: 1
#    resources:
#        mem_mb = 4096,
#        partition = config['queuing']['standard']['partition'],
#        account = config['queuing']['standard']['account'],
#        gres = config['queuing']['standard']['gres'],
#        time = '1:00:00'
#    output:
#        metrics_copy = COPY_CELLBENDER_METRICS
#    shell:
#        """
#        mkdir -p $(dirname {output.metrics_copy})

#        cp {input.metrics} {output.metrics_copy}
#        """

#####################################################################################################
##################### Merge klue results with cellbender adata and run scrublet #####################
#####################################################################################################

rule make_filt_adata:
    threads: 4
    resources:
        mem_mb = lambda wildcards, input, attempt: (input.size // BYTES_TO_MB) * 4,
        partition = config['queuing']['highmem']['partition'],
        account = config['queuing']['highmem']['account'],
        gres = config['queuing']['highmem']['gres'],
        time = '1:00:00'
    input:
        filt_h5 = rules.cellbender.output.filt_h5,
        unfilt_adata = rules.link_existing_unfilt_adata.output.unfilt_adata,
        genotype_counts = rules.klue_get_genotype_counts.output.ofile,
    output:
        filt_adata = "pipeline_output_cellbender/{plate}/{subpool}/filtered_adata_denoised.h5ad",
    run:
        add_meta_filter(input.filt_h5,
                        input.unfilt_adata,
                        input.genotype_counts,
                        wildcards,
                        bc_df,
                        kit,
                        chemistry,
                        sample_df,
                        output.filt_adata)

################################################################################
################################ Combine adatas ################################
################################################################################

rule make_tissue_adata:
    input:
        adatas = lambda wc:get_tissue_adatas(df, sample_df, wc, rules.make_filt_adata.output.filt_adata)
    threads: 2
    resources:
        mem_mb = lambda wildcards, input, attempt: (input.size // BYTES_TO_MB) * 2,
        partition = config['queuing']['highmem']['partition'],
        account = config['queuing']['highmem']['account'],
        gres = config['queuing']['highmem']['gres'],
        time = '2:00:00'
    output:
        adata = MAKE_TISSUE_ADATA_OUTPUT,
    run:
        concat_large_adatas(input.adatas, output.adata)

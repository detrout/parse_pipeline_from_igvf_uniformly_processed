import pandas as pd

from utils import *

include:
    rules/config.smk

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

include:
    rules/download.smk
    


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


include:
    rules/cellbender.yml

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

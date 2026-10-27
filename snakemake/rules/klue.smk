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
    threads: 24
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
    threads: 24
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

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

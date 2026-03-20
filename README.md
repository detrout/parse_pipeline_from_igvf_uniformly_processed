# parse_pipeline_from_igvf_uniformly_processed

This repository is a fork of the original parse pipeline and is adapted for IGVF-uniformly processed single-cell data.

Original repository:
https://github.com/mortazavilab/parse_pipeline/tree/main

## What changed in this fork

1. The workflow in [snakemake/Snakefile_updated_from_igvf.smk](snakemake/Snakefile_updated_from_igvf.smk) can start from precomputed IGVF h5ad files for downstream processing.
2. FASTQ files are still required because klue-based genetic demultiplexing is run from FASTQs.
3. Inputs are driven by experiment config TSV files in [configs](configs), plus [configs/sample_metadata.csv](configs/sample_metadata.csv).

## Inputs required

1. A config TSV (see examples in [configs](configs)) with at least these columns used by this Snakefile:
    fastq, fastq_r2, subpool, plate, lane, run, droplets_included, learning_rate, expected_cells, igvf_adata.
2. A sample metadata CSV at [configs/sample_metadata.csv](configs/sample_metadata.csv).
3. Reference settings in [configs/config.yml](configs/config.yml).

## Notes

1. The current Snakefile header has hardcoded values for `config_tsv`, `sample_csv`, `kit`, and `chemistry`. Update those as needed per run.
2. For full setup details (environments, cluster execution patterns, and broader background), use the original repository README.

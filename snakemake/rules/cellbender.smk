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
        'envs/cellbender.yaml'
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


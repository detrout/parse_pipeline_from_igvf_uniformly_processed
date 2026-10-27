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


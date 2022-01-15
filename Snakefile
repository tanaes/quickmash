import pandas as pd
from os.path import join

configfile: 'config.yaml'

samples = pd.read_table(config['samples']).set_index('sample')

rule all:
    input:
        "output/merged/genomes.sig",
        "output/compare/genomes_cmp.matrix.pdf",
        "output/classify/genomes.tax.csv"

rule sketch:
    input:
        r1 = lambda wildcards: samples.loc[wildcards.sample,
                                            'r1'],
        r2 = lambda wildcards: samples.loc[wildcards.sample,
                                            'r2']
    output:
        cat = temp("output/sketches/{sample}.fq"),
        sig = "output/sketches/{sample}.sig"
    params: 
    log:
        "output/logs/sketch/sketch-{sample}.log"
    conda:
        "envs/sourmash.yaml"
    shell:
        """
        cat {input} > {output.cat}

        sourmash sketch dna -p scaled=1000,k=31 -o {output.sig} {input} 2> {log}
        """

rule compare:
    input:
        expand("output/sketches/{sample}.sig",
               sample=samples.index)
    output:
        cmp = "output/compare/genomes_cmp",
        tree = "output/compare/genomes_cmp.dendro.pdf",
        matrix = "output/compare/genomes_cmp.matrix.pdf",
        csv = "output/compare/genomes_cmp.matrix.csv"
    conda:
        "envs/sourmash.yaml"
    params:
        outdir = "output/compare"
    log:
        "output/logs/compare/compare.log"
    shell:
        """
        sourmash compare -p 8 {input} -o {output.cmp}

        sourmash plot --pdf --labels {output.cmp} \
            --output-dir {params.outdir} \
            --csv {output.csv}
        """

rule merge:
    input:
        expand("output/sketches/{sample}.sig",
               sample=samples.index) 
    output:
        sig = "output/merged/genomes.sig"
    conda:
        "envs/sourmash.yaml"
    params:
    log:
        "output/logs/merged/merge.log"
    shell:
        """
        sourmash sig merge -o {output.sig} {input}
        """

rule classify:
    input:
        expand("output/sketches/{sample}.sig",
               sample=samples.index) 
    output:
        tax = "output/classify/genomes.tax.csv"
    conda:
        "envs/sourmash.yaml"
    params:
        db = config['params']['classify']['db']
    log:
        "output/logs/classify/classify.log"
    shell:
        """
        sourmash lca classify \
        --db {params.db} \
        --query {input} -o {output}
        """

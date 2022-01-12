import pandas as pd
from os.path import join

configfile: 'config.yaml'

samples = pd.read_table(config['samples']).set_index('sample')

rule all:
    input:
        expand(rules.sketch.output.sig, sample=samples.index)

rule trim_low_abund:
    input:
        fq1 = lambda wildcards: samples.loc[wildcards.sample,
                                            'fq1'],
        fq2 = lambda wildcards: samples.loc[wildcards.sample,
                                            'fq2']
    output:
        cat = temp("output/trimmed/{sample}.fq"),
        trimmed = temp("output/trimmed/{sample}.fq.abundtrimmed")
    params: 
        C = config['params']['trim_low_abund']['C'],
        Z = config['params']['trim_low_abund']['Z'],
        M = config['params']['trim_low_abund']['M']
    conda:
        "envs/sourmash.yaml"
    shell:
        """
        cat {input} > {output.cat}

        trim-low-abund -C {params.C} -Z {params.Z} -M {params.M} {output.cat}
        """

rule sketch:
    input:
        rules.trim_low_abund.output.trimmed
    output:
        sig = "output/sketches/{sample}.sig"
    params:
    conda:
        "envs/sourmash.yaml"
    shell:
        """
        sourmash sketch dna -p scaled=1000,k=31 {input} -o {output}
        """

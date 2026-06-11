"""
QC rules: fastp adapter trimming + read QC for WES FASTQs.
Wildcard: {sample} (key in config['wes_samples']).
"""

import os

OUTDIR = config["output_dir"]


def _fastp_input(wildcards):
    s = config["wes_samples"][wildcards.sample]
    return {"r1": s["R1"], "r2": s["R2"]}


rule fastp:
    """
    Adapter auto-detection + quality trimming.
    """
    input:
        unpack(_fastp_input),
    output:
        r1=os.path.join(OUTDIR, "wes/fastp/{sample}/{sample}_R1.trimmed.fq.gz"),
        r2=os.path.join(OUTDIR, "wes/fastp/{sample}/{sample}_R2.trimmed.fq.gz"),
        html=os.path.join(OUTDIR, "wes/fastp/{sample}/{sample}.fastp.html"),
        json=os.path.join(OUTDIR, "wes/fastp/{sample}/{sample}.fastp.json"),
    log:
        os.path.join(OUTDIR, "logs/fastp/{sample}.log"),
    threads: 4
    params:
        extra=config["params"]["fastp"]["extra"],
    conda:
        "../envs/fastp.yaml"
    shell:
        r"""
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --detect_adapter_for_pe \
            -h {output.html} -j {output.json} \
            -w {threads} {params.extra} \
            > {log} 2>&1
        """

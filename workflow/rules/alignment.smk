"""
Alignment: BWA-MEM2 -> samtools sort -> samtools index (single shell pipe).
Wildcard: {sample}
"""

import os

OUTDIR = config["output_dir"]


rule bwa_mem2:
    """
    Align trimmed paired-end FASTQ to GRCh38 with BWA-MEM2.
    Read group is required by GATK downstream:
        @RG\tID:{sample}\tSM:{sample}\tLB:{sample}\tPL:ILLUMINA
    BWA-MEM2 -> samtools sort is connected by a shell pipe (no intermediate SAM).
    """
    input:
        r1=os.path.join(OUTDIR, "wes/fastp/{sample}/{sample}_R1.trimmed.fq.gz"),
        r2=os.path.join(OUTDIR, "wes/fastp/{sample}/{sample}_R2.trimmed.fq.gz"),
        ref=config["reference"]["genome"],
    output:
        bam=temp(os.path.join(OUTDIR, "wes/aligned/{sample}/{sample}.sorted.bam")),
        bai=temp(os.path.join(OUTDIR, "wes/aligned/{sample}/{sample}.sorted.bam.bai")),
    log:
        os.path.join(OUTDIR, "logs/bwa_mem2/{sample}.log"),
    threads: config["params"]["bwa"]["threads"]
    params:
        rg=lambda wc: r"@RG\tID:{s}\tSM:{s}\tLB:{s}\tPL:ILLUMINA".format(s=wc.sample),
        extra=config["params"]["bwa"]["extra"],
    conda:
        "../envs/alignment.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {output.bam})"
        ( bwa-mem2 mem \
              -t {threads} \
              -R "{params.rg}" \
              {params.extra} \
              {input.ref} {input.r1} {input.r2} \
          | samtools sort -@ {threads} -o {output.bam} - \
        ) 2> {log}
        samtools index -@ {threads} {output.bam}
        """

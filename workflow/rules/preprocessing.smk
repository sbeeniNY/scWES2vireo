"""
Preprocessing: MarkDuplicates (GATK4) + BQSR (BaseRecalibrator -> ApplyBQSR).
Wildcard: {sample}

`config['reference']['exome_targets_bed']` is optional. If empty, BQSR and
ApplyBQSR run genome-wide (no `-L`).
"""

import os

OUTDIR = config["output_dir"]

# Optional capture-kit intervals.
_BED = (config["reference"].get("exome_targets_bed") or "").strip()
_BED_INPUT = [_BED] if _BED else []
_L_FLAG = f"-L {_BED}" if _BED else ""


rule mark_duplicates:
    """
    Mark PCR/optical duplicates with GATK4 MarkDuplicatesSpark
    (uses Picard's MarkDuplicates engine; outputs sorted+indexed BAM).
    """
    input:
        bam=os.path.join(OUTDIR, "wes/aligned/{sample}/{sample}.sorted.bam"),
        bai=os.path.join(OUTDIR, "wes/aligned/{sample}/{sample}.sorted.bam.bai"),
    output:
        bam=temp(os.path.join(OUTDIR, "wes/dedup/{sample}/{sample}.dedup.bam")),
        bai=temp(os.path.join(OUTDIR, "wes/dedup/{sample}/{sample}.dedup.bam.bai")),
        metrics=os.path.join(OUTDIR, "wes/dedup/{sample}/{sample}.dedup.metrics.txt"),
    log:
        os.path.join(OUTDIR, "logs/mark_duplicates/{sample}.log"),
    threads: 4
    params:
        java_opts=config["params"]["gatk"]["java_options"],
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        # MarkDuplicatesSpark stages the BAM in a "{output}.parts" dir then commits.
        # (1) -M metrics is written by a plain writer that does NOT mkdir, so the
        #     output dir must exist. (2) MarkDuplicatesSpark REFUSES to run if a
        #     leftover ".parts" dir from a killed attempt exists ("Output directory
        #     ... .parts already exists"), so every retry fails until it is removed.
        #     Snakemake does not clean it because ".parts" is not a declared output.
        mkdir -p "$(dirname {output.bam})"
        rm -rf {output.bam}.parts {output.bam} {output.bai}
        gatk --java-options "{params.java_opts}" MarkDuplicatesSpark \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --conf 'spark.executor.cores={threads}' \
            > {log} 2>&1
        """


rule base_recalibrator:
    """
    Build BQSR recalibration table.
    Uses exome target intervals if provided; otherwise whole genome.
    """
    input:
        bam=os.path.join(OUTDIR, "wes/dedup/{sample}/{sample}.dedup.bam"),
        bai=os.path.join(OUTDIR, "wes/dedup/{sample}/{sample}.dedup.bam.bai"),
        ref=config["reference"]["genome"],
        dbsnp=config["reference"]["known_sites"]["dbsnp"],
        mills=config["reference"]["known_sites"]["mills"],
        g1k=config["reference"]["known_sites"]["g1000_indels"],
        bed=_BED_INPUT,
    output:
        recal_table=os.path.join(OUTDIR, "wes/recal/{sample}/{sample}.recal.table"),
    log:
        os.path.join(OUTDIR, "logs/base_recalibrator/{sample}.log"),
    threads: 4
    params:
        java_opts=config["params"]["gatk"]["java_options"],
        intervals=_L_FLAG,
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {output.recal_table})"
        gatk --java-options "{params.java_opts}" BaseRecalibrator \
            -R {input.ref} \
            -I {input.bam} \
            {params.intervals} \
            --known-sites {input.dbsnp} \
            --known-sites {input.mills} \
            --known-sites {input.g1k} \
            -O {output.recal_table} \
            > {log} 2>&1
        """


rule apply_bqsr:
    """
    Apply recalibration table -> final analysis-ready BAM.
    """
    input:
        bam=os.path.join(OUTDIR, "wes/dedup/{sample}/{sample}.dedup.bam"),
        bai=os.path.join(OUTDIR, "wes/dedup/{sample}/{sample}.dedup.bam.bai"),
        recal_table=os.path.join(OUTDIR, "wes/recal/{sample}/{sample}.recal.table"),
        ref=config["reference"]["genome"],
        bed=_BED_INPUT,
    output:
        bam=os.path.join(OUTDIR, "wes/recal/{sample}/{sample}.recal.bam"),
        bai=os.path.join(OUTDIR, "wes/recal/{sample}/{sample}.recal.bam.bai"),
    log:
        os.path.join(OUTDIR, "logs/apply_bqsr/{sample}.log"),
    threads: 4
    params:
        java_opts=config["params"]["gatk"]["java_options"],
        intervals=_L_FLAG,
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {output.bam})"
        gatk --java-options "{params.java_opts}" ApplyBQSR \
            -R {input.ref} \
            -I {input.bam} \
            {params.intervals} \
            --bqsr-recal-file {input.recal_table} \
            -O {output.bam} \
            > {log} 2>&1
        # GATK writes .bai (Picard style); also create .bam.bai for tools that expect it
        if [ -f {output.bam}.bai ]; then
            cp -f {output.bam}.bai {output.bai} || true
        else
            samtools index -@ {threads} {output.bam}
        fi
        """

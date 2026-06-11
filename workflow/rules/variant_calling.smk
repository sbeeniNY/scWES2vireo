"""
Variant calling:
  HaplotypeCaller (GVCF) -> CombineGVCFs -> GenotypeGVCFs -> VariantFiltration
  -> PASS-only cohort VCF.
Wildcards: {sample} (per-sample); cohort steps have no wildcard.

`config['reference']['exome_targets_bed']` is optional. If empty,
HaplotypeCaller runs genome-wide (no `-L`). Note this is significantly
slower than restricting to capture regions for WES data.
"""

import os

OUTDIR = config["output_dir"]
WES_SAMPLES = list(config["wes_samples"].keys())

# Optional capture-kit intervals.
_BED = (config["reference"].get("exome_targets_bed") or "").strip()
_BED_INPUT = [_BED] if _BED else []
_L_FLAG = f"-L {_BED}" if _BED else ""


rule haplotype_caller:
    """
    Per-sample GVCF. Uses exome target intervals if provided; otherwise
    runs genome-wide (slower but correct — HaplotypeCaller only emits
    variants where reads exist, so empty regions are skipped quickly).
    """
    input:
        bam=os.path.join(OUTDIR, "wes/recal/{sample}/{sample}.recal.bam"),
        bai=os.path.join(OUTDIR, "wes/recal/{sample}/{sample}.recal.bam.bai"),
        ref=config["reference"]["genome"],
        bed=_BED_INPUT,
    output:
        gvcf=os.path.join(OUTDIR, "wes/gvcf/{sample}/{sample}.g.vcf.gz"),
        tbi=os.path.join(OUTDIR, "wes/gvcf/{sample}/{sample}.g.vcf.gz.tbi"),
    log:
        os.path.join(OUTDIR, "logs/haplotype_caller/{sample}.log"),
    threads: 4
    params:
        java_opts=config["params"]["gatk"]["java_options"],
        intervals=_L_FLAG,
        extra=config["params"]["gatk"]["haplotypecaller_extra"],
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        mkdir -p "$(dirname {output.gvcf})"
        gatk --java-options "{params.java_opts}" HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            {params.intervals} \
            -ERC GVCF \
            -O {output.gvcf} \
            {params.extra} \
            > {log} 2>&1
        # HaplotypeCaller usually writes .tbi itself, but guarantee it: this GATK
        # build / GPFS sometimes leaves the gvcf unindexed -> MissingOutputException.
        if [ ! -f {output.tbi} ]; then
            gatk --java-options "{params.java_opts}" IndexFeatureFile -I {output.gvcf} >> {log} 2>&1
        fi
        """


rule combine_gvcfs:
    """
    Merge per-sample GVCFs into a multi-sample GVCF.
    """
    input:
        gvcfs=expand(
            os.path.join(OUTDIR, "wes/gvcf/{sample}/{sample}.g.vcf.gz"),
            sample=WES_SAMPLES,
        ),
        tbis=expand(
            os.path.join(OUTDIR, "wes/gvcf/{sample}/{sample}.g.vcf.gz.tbi"),
            sample=WES_SAMPLES,
        ),
        ref=config["reference"]["genome"],
    output:
        gvcf=os.path.join(OUTDIR, "wes/joint/combined.g.vcf.gz"),
        tbi=os.path.join(OUTDIR, "wes/joint/combined.g.vcf.gz.tbi"),
    log:
        os.path.join(OUTDIR, "logs/combine_gvcfs/cohort.log"),
    threads: 4
    params:
        java_opts=config["params"]["gatk"]["java_options"],
        variant_args=lambda wc, input: " ".join(f"-V {p}" for p in input.gvcfs),
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.java_opts}" CombineGVCFs \
            -R {input.ref} \
            {params.variant_args} \
            -O {output.gvcf} \
            > {log} 2>&1
        if [ ! -f {output.tbi} ]; then
            gatk --java-options "{params.java_opts}" IndexFeatureFile -I {output.gvcf} >> {log} 2>&1
        fi
        """


rule genotype_gvcfs:
    """
    Joint genotyping over the cohort.
    """
    input:
        gvcf=os.path.join(OUTDIR, "wes/joint/combined.g.vcf.gz"),
        tbi=os.path.join(OUTDIR, "wes/joint/combined.g.vcf.gz.tbi"),
        ref=config["reference"]["genome"],
    output:
        vcf=os.path.join(OUTDIR, "wes/joint/genotyped.vcf.gz"),
        tbi=os.path.join(OUTDIR, "wes/joint/genotyped.vcf.gz.tbi"),
    log:
        os.path.join(OUTDIR, "logs/genotype_gvcfs/cohort.log"),
    threads: 4
    params:
        java_opts=config["params"]["gatk"]["java_options"],
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.java_opts}" GenotypeGVCFs \
            -R {input.ref} \
            -V {input.gvcf} \
            -O {output.vcf} \
            > {log} 2>&1
        if [ ! -f {output.tbi} ]; then
            gatk --java-options "{params.java_opts}" IndexFeatureFile -I {output.vcf} >> {log} 2>&1
        fi
        """


rule select_snps:
    input:
        vcf=os.path.join(OUTDIR, "wes/joint/genotyped.vcf.gz"),
        tbi=os.path.join(OUTDIR, "wes/joint/genotyped.vcf.gz.tbi"),
        ref=config["reference"]["genome"],
    output:
        vcf=temp(os.path.join(OUTDIR, "wes/joint/cohort.snps.vcf.gz")),
        tbi=temp(os.path.join(OUTDIR, "wes/joint/cohort.snps.vcf.gz.tbi")),
    log:
        os.path.join(OUTDIR, "logs/select_variants/snps.log"),
    threads: 2
    params:
        java_opts=config["params"]["gatk"]["java_options"],
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.java_opts}" SelectVariants \
            -R {input.ref} -V {input.vcf} \
            --select-type-to-include SNP \
            -O {output.vcf} \
            > {log} 2>&1
        if [ ! -f {output.tbi} ]; then
            gatk --java-options "{params.java_opts}" IndexFeatureFile -I {output.vcf} >> {log} 2>&1
        fi
        """


rule select_indels:
    input:
        vcf=os.path.join(OUTDIR, "wes/joint/genotyped.vcf.gz"),
        tbi=os.path.join(OUTDIR, "wes/joint/genotyped.vcf.gz.tbi"),
        ref=config["reference"]["genome"],
    output:
        vcf=temp(os.path.join(OUTDIR, "wes/joint/cohort.indels.vcf.gz")),
        tbi=temp(os.path.join(OUTDIR, "wes/joint/cohort.indels.vcf.gz.tbi")),
    log:
        os.path.join(OUTDIR, "logs/select_variants/indels.log"),
    threads: 2
    params:
        java_opts=config["params"]["gatk"]["java_options"],
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.java_opts}" SelectVariants \
            -R {input.ref} -V {input.vcf} \
            --select-type-to-include INDEL \
            -O {output.vcf} \
            > {log} 2>&1
        if [ ! -f {output.tbi} ]; then
            gatk --java-options "{params.java_opts}" IndexFeatureFile -I {output.vcf} >> {log} 2>&1
        fi
        """


rule filter_snps:
    """
    GATK hard-filter for SNPs.
    """
    input:
        vcf=os.path.join(OUTDIR, "wes/joint/cohort.snps.vcf.gz"),
        tbi=os.path.join(OUTDIR, "wes/joint/cohort.snps.vcf.gz.tbi"),
        ref=config["reference"]["genome"],
    output:
        vcf=temp(os.path.join(OUTDIR, "wes/joint/cohort.snps.filtered.vcf.gz")),
        tbi=temp(os.path.join(OUTDIR, "wes/joint/cohort.snps.filtered.vcf.gz.tbi")),
    log:
        os.path.join(OUTDIR, "logs/variant_filtration/snps.log"),
    threads: 2
    params:
        java_opts=config["params"]["gatk"]["java_options"],
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.java_opts}" VariantFiltration \
            -R {input.ref} -V {input.vcf} \
            --filter-expression "QD < 2.0"               --filter-name "QD2" \
            --filter-expression "FS > 60.0"              --filter-name "FS60" \
            --filter-expression "MQ < 40.0"              --filter-name "MQ40" \
            --filter-expression "MQRankSum < -12.5"      --filter-name "MQRankSum-12.5" \
            --filter-expression "ReadPosRankSum < -8.0"  --filter-name "ReadPosRankSum-8" \
            -O {output.vcf} \
            > {log} 2>&1
        if [ ! -f {output.tbi} ]; then
            gatk --java-options "{params.java_opts}" IndexFeatureFile -I {output.vcf} >> {log} 2>&1
        fi
        """


rule filter_indels:
    """
    GATK hard-filter for INDELs.
    """
    input:
        vcf=os.path.join(OUTDIR, "wes/joint/cohort.indels.vcf.gz"),
        tbi=os.path.join(OUTDIR, "wes/joint/cohort.indels.vcf.gz.tbi"),
        ref=config["reference"]["genome"],
    output:
        vcf=temp(os.path.join(OUTDIR, "wes/joint/cohort.indels.filtered.vcf.gz")),
        tbi=temp(os.path.join(OUTDIR, "wes/joint/cohort.indels.filtered.vcf.gz.tbi")),
    log:
        os.path.join(OUTDIR, "logs/variant_filtration/indels.log"),
    threads: 2
    params:
        java_opts=config["params"]["gatk"]["java_options"],
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.java_opts}" VariantFiltration \
            -R {input.ref} -V {input.vcf} \
            --filter-expression "QD < 2.0"                --filter-name "QD2" \
            --filter-expression "FS > 200.0"              --filter-name "FS200" \
            --filter-expression "ReadPosRankSum < -20.0"  --filter-name "ReadPosRankSum-20" \
            -O {output.vcf} \
            > {log} 2>&1
        if [ ! -f {output.tbi} ]; then
            gatk --java-options "{params.java_opts}" IndexFeatureFile -I {output.vcf} >> {log} 2>&1
        fi
        """


rule merge_filtered_vcfs:
    """
    Merge filtered SNP + INDEL VCFs back together (FILTER column preserved).
    """
    input:
        snps=os.path.join(OUTDIR, "wes/joint/cohort.snps.filtered.vcf.gz"),
        snps_tbi=os.path.join(OUTDIR, "wes/joint/cohort.snps.filtered.vcf.gz.tbi"),
        indels=os.path.join(OUTDIR, "wes/joint/cohort.indels.filtered.vcf.gz"),
        indels_tbi=os.path.join(OUTDIR, "wes/joint/cohort.indels.filtered.vcf.gz.tbi"),
    output:
        vcf=os.path.join(OUTDIR, "wes/joint/cohort.filtered.vcf.gz"),
        tbi=os.path.join(OUTDIR, "wes/joint/cohort.filtered.vcf.gz.tbi"),
    log:
        os.path.join(OUTDIR, "logs/merge_vcfs/cohort.log"),
    threads: 2
    params:
        java_opts=config["params"]["gatk"]["java_options"],
    conda:
        "../envs/gatk.yaml"
    shell:
        r"""
        gatk --java-options "{params.java_opts}" MergeVcfs \
            -I {input.snps} -I {input.indels} \
            -O {output.vcf} \
            > {log} 2>&1
        # MergeVcfs creates .tbi already; guarantee it via gatk (tabix may be
        # absent in the gatk env).
        if [ ! -f {output.tbi} ]; then
            gatk --java-options "{params.java_opts}" IndexFeatureFile -I {output.vcf} >> {log} 2>&1
        fi
        """


rule select_pass:
    """
    Keep only PASS records to feed Vireo.
    """
    input:
        vcf=os.path.join(OUTDIR, "wes/joint/cohort.filtered.vcf.gz"),
        tbi=os.path.join(OUTDIR, "wes/joint/cohort.filtered.vcf.gz.tbi"),
    output:
        vcf=os.path.join(OUTDIR, "wes/joint/cohort.pass.vcf.gz"),
        tbi=os.path.join(OUTDIR, "wes/joint/cohort.pass.vcf.gz.tbi"),
    log:
        os.path.join(OUTDIR, "logs/select_pass/cohort.log"),
    threads: 2
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
        bcftools view -f PASS -Oz -o {output.vcf} {input.vcf} 2> {log}
        tabix -p vcf {output.vcf}
        """

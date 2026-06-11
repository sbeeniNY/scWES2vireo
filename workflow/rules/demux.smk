"""
Genetic demultiplexing per pool:
  cellSNP-lite (pool BAM + common SNPs) -> per-pool donor genotype subset -> Vireo.
Wildcard: {pool} (key in config['pools']).
"""

import os

OUTDIR = config["output_dir"]
POOLS = list(config["pools"].keys())


def _pool_bam(wc):
    cr_dir = config["pools"][wc.pool]["cellranger_dir"]
    cr_id = os.path.basename(os.path.dirname(cr_dir))
    return os.path.join(cr_dir, "per_sample_outs", cr_id, "sample_alignments.bam")


def _pool_barcodes(wc):
    return os.path.join(
        config["pools"][wc.pool]["cellranger_dir"],
        "filtered_feature_bc_matrix",
        "barcodes.tsv.gz",
    )


def _pool_donors_csv(pool):
    """Return the comma-joined donor list for a given pool, e.g. 'P029,P387'."""
    return ",".join(config["pools"][pool]["donors"])


rule cellsnp_lite:
    """
    Pileup pooled scRNA BAM at common SNP positions, per cell barcode.
    """
    input:
        bam=lambda wc: _pool_bam(wc),
        barcodes=lambda wc: _pool_barcodes(wc),
        common_vcf=config["reference"]["common_snps_vcf"],
    output:
        outdir=directory(os.path.join(OUTDIR, "demux/{pool}/cellsnp")),
        cells_vcf=os.path.join(OUTDIR, "demux/{pool}/cellsnp/cellSNP.cells.vcf.gz"),
        base_vcf=os.path.join(OUTDIR, "demux/{pool}/cellsnp/cellSNP.base.vcf.gz"),
    log:
        os.path.join(OUTDIR, "logs/cellsnp_lite/{pool}.log"),
    threads: config["params"]["cellsnp"]["threads"]
    params:
        min_maf=config["params"]["cellsnp"]["min_MAF"],
        min_count=config["params"]["cellsnp"]["min_count"],
        extra=config["params"]["cellsnp"]["extra"],
    conda:
        "../envs/cellsnp.yaml"
    shell:
        r"""
        mkdir -p {output.outdir}
        cellsnp-lite \
            -s {input.bam} \
            -b {input.barcodes} \
            -O {output.outdir} \
            -R {input.common_vcf} \
            --minMAF {params.min_maf} \
            --minCOUNT {params.min_count} \
            -p {threads} \
            --gzip \
            --genotype \
            {params.extra} \
            > {log} 2>&1
        """


rule subset_pool_donor_vcf:
    """
    Extract only the donors present in this pool from the cohort PASS VCF.
    Vireo MUST receive a VCF subset to the pool's donors.
    """
    input:
        vcf=os.path.join(OUTDIR, "wes/joint/cohort.pass.vcf.gz"),
        tbi=os.path.join(OUTDIR, "wes/joint/cohort.pass.vcf.gz.tbi"),
    output:
        vcf=os.path.join(OUTDIR, "demux/{pool}/donor_genotype.vcf.gz"),
        tbi=os.path.join(OUTDIR, "demux/{pool}/donor_genotype.vcf.gz.tbi"),
    log:
        os.path.join(OUTDIR, "logs/subset_pool_donor_vcf/{pool}.log"),
    threads: 2
    params:
        donors=lambda wc: _pool_donors_csv(wc.pool),
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
        # -s : sample subset, --force-samples avoids failure if any name mismatch
        # -a : trim alt alleles not seen in subset
        # -e 'GT[*]="mis"' : drop sites missing in any selected sample (Vireo-friendly)
        bcftools view -s {params.donors} --force-samples -a {input.vcf} 2> {log} \
          | bcftools view -e 'GT[*]="mis"' -Oz -o {output.vcf} 2>> {log}
        tabix -p vcf {output.vcf}
        """


rule vireo:
    """
    Assign each cell barcode to a donor (or doublet/unassigned) using
    the cellSNP cell pileups and the pool-specific donor genotype VCF.
    """
    input:
        cellsnp_dir=os.path.join(OUTDIR, "demux/{pool}/cellsnp"),
        donor_vcf=os.path.join(OUTDIR, "demux/{pool}/donor_genotype.vcf.gz"),
        donor_tbi=os.path.join(OUTDIR, "demux/{pool}/donor_genotype.vcf.gz.tbi"),
    output:
        outdir=directory(os.path.join(OUTDIR, "demux/{pool}/vireo")),
        donor_ids=os.path.join(OUTDIR, "demux/{pool}/vireo/donor_ids.tsv"),
    log:
        os.path.join(OUTDIR, "logs/vireo/{pool}.log"),
    threads: 4
    params:
        n_donor=lambda wc: config["pools"][wc.pool]["n_donor"],
        extra=config["params"]["vireo"]["extra"],
    conda:
        "../envs/vireo.yaml"
    shell:
        r"""
        # Ignore ~/.local site-packages: a user-site scipy/numpy there shadows
        # the conda env's pinned versions and breaks vireoSNP import
        # (ValueError: All ufuncs must have type numpy.ufunc).
        export PYTHONNOUSERSITE=1
        mkdir -p {output.outdir}
        vireo \
            -c {input.cellsnp_dir} \
            -d {input.donor_vcf} \
            -N {params.n_donor} \
            -o {output.outdir} \
            {params.extra} \
            > {log} 2>&1
        """


rule demux_summary:
    """
    Aggregate per-pool donor_ids.tsv into a single summary table:
      pool, donor_id, mapped_label, n_cells, doublet_rate, unassigned_rate
    Uses config['sample_labels'] to map donor_id -> human-readable label.
    """
    input:
        donor_ids=expand(
            os.path.join(OUTDIR, "demux/{pool}/vireo/donor_ids.tsv"),
            pool=POOLS,
        ),
        config_yaml="config/config.yaml",
    output:
        tsv=os.path.join(OUTDIR, "demux/demux_summary.tsv"),
    log:
        os.path.join(OUTDIR, "logs/demux_summary/all.log"),
    threads: 1
    params:
        script=os.path.abspath(os.path.join("workflow", "scripts", "demux_summary.py")),
    conda:
        "../envs/cellsnp.yaml"
    shell:
        r"""
        # Run as a plain CLI (not a Snakemake `script:` directive): the script:
        # machinery injects the host snakemake env's path, pulling a py3.11
        # pandas/numpy into this py3.10 env and crashing on import. unset
        # PYTHONPATH + PYTHONNOUSERSITE keep imports inside this conda env.
        unset PYTHONPATH
        export PYTHONNOUSERSITE=1
        python3 {params.script} \
            --donor-ids {input.donor_ids} \
            --config {input.config_yaml} \
            --out {output.tsv} \
            > {log} 2>&1
        """


rule combine_cell_labels:
    """
    Flatten every pool's Vireo per-cell assignment into one cohort table
    (pool, cell, donor, sample) for direct loading into Seurat. doublet /
    unassigned cells carry those literals through in both donor and sample.
    Depends on demux_summary so it runs last and verifies the pool set.
    """
    input:
        donor_ids=expand(
            os.path.join(OUTDIR, "demux/{pool}/vireo/donor_ids.tsv"),
            pool=POOLS,
        ),
        config_yaml="config/config.yaml",
        summary=os.path.join(OUTDIR, "demux/demux_summary.tsv"),
    output:
        tsv=os.path.join(OUTDIR, "demux/cell_donor_sample.tsv"),
    log:
        os.path.join(OUTDIR, "logs/combine_cell_labels/all.log"),
    threads: 1
    params:
        script=os.path.abspath(os.path.join("workflow", "scripts", "combine_cell_labels.py")),
    conda:
        "../envs/cellsnp.yaml"
    shell:
        r"""
        # Plain CLI under shell: (not script:) to avoid host-env path injection;
        # unset PYTHONPATH + PYTHONNOUSERSITE keep imports inside this conda env.
        unset PYTHONPATH
        export PYTHONNOUSERSITE=1
        python3 {params.script} \
            --donor-ids {input.donor_ids} \
            --config {input.config_yaml} \
            --summary {input.summary} \
            --out {output.tsv} \
            > {log} 2>&1
        """


rule unassigned_decompose:
    """
    Classify Vireo "unassigned" cells into data_limited / recoverable /
    true_ambiguous using per-cell n_vars and best-singlet posterior, so we can
    tell how much of the unassigned rate is recoverable vs genuinely ambiguous.
    Non-destructive: reads donor_ids.tsv, never edits Vireo output.
    """
    input:
        donor_ids=expand(
            os.path.join(OUTDIR, "demux/{pool}/vireo/donor_ids.tsv"),
            pool=POOLS,
        ),
    output:
        tsv=os.path.join(OUTDIR, "demux/qc_unassigned/unassigned_decomposition_per_pool.tsv"),
    log:
        os.path.join(OUTDIR, "logs/unassigned_decompose/all.log"),
    threads: 1
    params:
        demux_dir=os.path.join(OUTDIR, "demux"),
        out_dir=os.path.join(OUTDIR, "demux/qc_unassigned"),
        # Anchor to the working dir (repo root, where snakemake is launched and
        # config/cluster paths resolve). `workflow.snakefile` is unreliable here:
        # inside an included .smk it resolves to the .smk's dir, yielding the
        # wrong workflow/scripts/ path instead of the top-level scripts/.
        script=os.path.abspath(os.path.join("scripts", "16_unassigned_decompose.py")),
        min_vars=config["params"]["qc"]["min_vars"],
        recover_prob=config["params"]["qc"]["recover_prob"],
        assign_thr=config["params"]["qc"]["assign_thr"],
        pools=" ".join(POOLS),
    conda:
        "../envs/vireo.yaml"
    shell:
        r"""
        # See vireo rule: ignore ~/.local packages that shadow the conda env.
        export PYTHONNOUSERSITE=1
        python3 {params.script} \
            --demux-dir {params.demux_dir} \
            --out-dir {params.out_dir} \
            --pools {params.pools} \
            --min-vars {params.min_vars} \
            --recover-prob {params.recover_prob} \
            --assign-thr {params.assign_thr} \
            > {log} 2>&1
        """

# scWES2vireo
## WES → GATK joint genotyping → cellSNP-lite + Vireo

**WES-guided demultiplexing for pooled single-cell RNA-seq using known donor genotypes.**

`scWES2vireo` is an end-to-end Snakemake workflow for assigning pooled scRNA-seq
barcodes to donors when matched whole-exome sequencing (WES) data are available.
It calls donor germline variants from WES with the GATK best-practice germline
pipeline, performs joint genotyping across all donors, and then runs
`cellSNP-lite` + `Vireo` per pool. Because the WES cohort VCF provides **known
donor genotypes**, Vireo runs in supervised mode: each cell is assigned to a
named donor (or `doublet` / `unassigned`) rather than to an anonymous cluster.

## Why scWES2vireo?

When matched WES exists, supervised demultiplexing with known genotypes is
direct and interpretable — clusters never need a post-hoc mapping to donors. The
workflow is built around two practical realities of 10x 3′ data:

- **Self-contained WES variant calling.** No external variant-calling pipeline
  is required: FASTQ → cohort VCF is handled by GATK (BWA-MEM2, MarkDuplicates,
  BQSR, HaplotypeCaller, joint GenotypeGVCFs, hard filtering).
- **Coverage-aware cellSNP defaults.** Exome variants sit in coding regions
  while 10x 3′ reads are 3′-UTR biased, so only a fraction of WES sites are
  covered. The cellSNP defaults are tuned to retain sparse donor-discriminating
  sites instead of discarding them (see *Unassigned recovery*).
- **Honest unassigned accounting.** A companion step decomposes Vireo's
  `unassigned` cells into *data-limited* vs *recoverable* vs *truly ambiguous*,
  so a high unassigned rate can be interpreted rather than blindly force-called.
- **Seurat-ready output.** A single per-cell `pool / cell / donor / sample`
  table is emitted for direct loading into R.

### Compared with a Souporcell-based pipeline

`scWES2vireo` is a **supervised** demultiplexer: Vireo assigns each cell using
the WES donor genotypes directly. This is the key difference from a de novo
approach such as [`scWESdemux`](https://github.com/sbeeniNY/scWESdemux)
(nf-core/sarek multi-caller consensus → Souporcell), where cells are first
clustered without identity and clusters are then mapped to donors *post hoc* by
genotype concordance.

| Aspect | **scWES2vireo** (cellSNP-lite + Vireo) | scWESdemux (Sarek + Souporcell) |
| --- | --- | --- |
| Demux strategy | **Supervised** — Vireo uses WES known donor genotypes | **De novo** — Souporcell clusters, then clusters mapped to donors by concordance |
| Donor identity | **Direct**: each cell → named donor with a posterior; no cluster-mapping step | **Indirect**: numeric cluster → donor mapping, done after the run |
| Failure mode | Conservative `unassigned` on low-coverage cells (made explicit, recoverable) | Cluster↔donor mapping can swap/merge when donors are genetically close, a pool is small, or a cluster is split |
| WES variant calling | Single GATK best-practice chain (self-contained) | nf-core/sarek 3-caller consensus (DeepVariant + FreeBayes + HaplotypeCaller) |
| Software stack | Snakemake + conda only | Snakemake + **Nextflow** + Souporcell **container** |
| Extra modeling | Unassigned-cell decomposition (recoverable vs ambiguous) | Built-in ambient-RNA estimate + doublet model |

**When scWES2vireo is the better choice**

- **Donors are known and you want their identity tied in from the start.**
  Supervised assignment removes the cluster→donor mapping step entirely, so there
  is no risk of a whole cluster being assigned to the wrong donor — a real
  failure mode for de novo methods when donors are genetically similar, when a
  pool has few cells, or when Souporcell splits/merges a cluster.
- **Interpretability.** Every cell carries a per-donor posterior
  (`prob_max`, `prob_doublet`) and, for borderline cells, an explicit
  `data_limited` / `recoverable` / `true_ambiguous` label — so a high unassigned
  rate can be diagnosed instead of force-called.
- **Lighter, simpler stack.** No Nextflow and no Souporcell container: the whole
  pipeline is Snakemake + conda, which is easier to install, port, and debug.
- **Fewer moving parts in variant calling.** One GATK chain rather than a
  three-caller consensus, when you do not need multi-caller hardening.

**When the Souporcell pipeline is preferable**

- WES sites are very poorly covered by the scRNA reads, so de novo clustering on
  a genome-wide panel extracts signal that supervised assignment cannot.
- You want Souporcell's built-in ambient-RNA estimation and doublet model, or
  multi-caller consensus genotypes to guard against single-caller artifacts.

The two are complementary. When matched WES exists and its sites are reasonably
covered, supervised Vireo assignment is the most direct and interpretable
option; running both and checking concordance is a strong validation strategy.

## Workflow overview

```text
Matched WES FASTQs (per donor)
        │
        ▼
Part 1: GATK germline (Snakemake)
  - fastp trim
  - BWA-MEM2 alignment → samtools sort
  - MarkDuplicates
  - BQSR (BaseRecalibrator → ApplyBQSR)
  - HaplotypeCaller (per-donor GVCF)
  - CombineGVCFs → GenotypeGVCFs (joint)
  - hard-filter SNPs + INDELs → merge → PASS-only cohort VCF
        │
        ▼
cohort.pass.vcf.gz  (known donor genotypes)
        │
        ▼
Part 2: demultiplexing (Snakemake, per pool)
  - subset cohort VCF to the pool's donors
  - cellSNP-lite pileup on the pooled scRNA-seq BAM at common SNPs
  - Vireo with known donor genotypes (-N n_donor)
  - cohort summary + per-cell donor/sample table
  - unassigned-cell decomposition
```

Part 2 is triggered automatically once `cohort.pass.vcf.gz` exists — there is no
manual hand-off between the two parts.

## Pipeline structure

```text
scWES2vireo/
├── config/
│   └── config.yaml                 # central configuration (paths + parameters)
├── workflow/
│   ├── Snakefile                   # entry point, rule all
│   ├── rules/
│   │   ├── qc.smk                  # fastp
│   │   ├── alignment.smk           # BWA-MEM2 → samtools sort
│   │   ├── preprocessing.smk       # MarkDuplicates, BQSR
│   │   ├── variant_calling.smk     # HaplotypeCaller → joint genotyping → filter
│   │   └── demux.smk               # cellSNP-lite, donor subset, Vireo, summaries
│   ├── envs/                       # per-rule conda environments
│   └── scripts/
│       ├── demux_summary.py        # cohort donor cell-count summary
│       └── combine_cell_labels.py  # per-cell pool/cell/donor/sample table
├── scripts/
│   ├── 16_unassigned_decompose.py  # unassigned-cell decomposition (used by workflow)
│   └── demux_summary_standalone.py # standalone summary utility (no Snakemake)
├── cluster/
│   ├── config.yaml                 # Snakemake 8 cluster profile (LSF / bsub)
│   ├── lsf_submit.sh               # bsub wrapper — prints numeric job ID
│   └── lsf_status.sh               # bjobs wrapper — running/success/failed
└── README.md
```

---

## Requirements

### Software

- Snakemake (≥ 8) with `snakemake-executor-plugin-cluster-generic`
- Conda / Mamba (per-rule environments are created automatically with
  `--software-deployment-method conda`)
- The per-rule tools are pinned in `workflow/envs/`: `fastp`, `bwa-mem2`,
  `samtools`, `gatk4`, `bcftools`/`htslib`, `cellsnp-lite`, `vireosnp`.

The workflow was developed for an HPC environment using LSF, but the Snakemake
profile in `cluster/` can be adapted to other schedulers (or run locally).

### Reference files

All reference inputs are set in `config/config.yaml`.

1. **Reference genome FASTA.** Use the **same** FASTA and chromosome naming for
   WES and scRNA-seq. For Cell Ranger BAMs this is normally the FASTA inside the
   matching Cell Ranger reference directory. It must be BWA-MEM2 indexed, with a
   `.fai` and a `.dict`:

   ```bash
   bwa-mem2 index   /path/to/reference/genome.fa
   samtools faidx   /path/to/reference/genome.fa
   gatk CreateSequenceDictionary -R /path/to/reference/genome.fa
   ```

2. **GATK BQSR known-sites (hg38).** dbSNP, Mills indels, and 1000G known
   indels from the Broad public bucket; index each with `gatk IndexFeatureFile`:

   ```bash
   # https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0
   wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
   wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
   wget https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
   ```

3. **cellSNP-lite common-SNP panel.** A genome-wide common SNP list (1000G,
   MAF > 5%, hg38). It must use the **same** chromosome naming as the genome
   FASTA — rename to chr-prefixed if needed (see comments in `config.yaml`):

   ```bash
   wget -O genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz \
     "https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz/download"
   ```

4. **Exome targets BED — optional.** If `reference.exome_targets_bed` is left
   empty (`""`), BQSR and HaplotypeCaller run without interval restriction (see
   *Running without an exome targets BED*). Provide a capture-kit BED path to
   restrict with `-L`.

---

## Input data

Edit `config/config.yaml` to define:

- reference genome FASTA + GATK known-sites + cellSNP common-SNP panel
- WES R1/R2 FASTQ for each donor (`wes_samples`)
- scRNA-seq Cell Ranger `outs/` directory + donor IDs + `n_donor` per pool (`pools`)
- optional donor → human-readable label mapping (`sample_labels`)
- tool parameters and the output directory

Donor IDs must be **identical** across `wes_samples`, `pools.*.donors`, and
`sample_labels`. Example pool configuration:

```yaml
pools:
  Pool1:
    cellranger_dir: "/path/to/cellranger/Pool1/outs"
    donors: ["DonorA", "DonorB"]
    n_donor: 2
  Pool2:
    cellranger_dir: "/path/to/cellranger/Pool2/outs"
    donors: ["DonorB", "DonorC"]
    n_donor: 2
```

---

## Quick start

### 1. Install the cluster plugin (one-time)

```bash
conda activate <your-snakemake-env>
pip install snakemake-executor-plugin-cluster-generic
chmod +x cluster/lsf_submit.sh cluster/lsf_status.sh
```

Edit `cluster/config.yaml` to set your LSF project (`-P acc_YourProject`) and
queue, or adapt the profile to your scheduler.

### 2. Configure

Edit `config/config.yaml` (reference files, WES FASTQs, pools, labels, output
directory). Prepare the reference indices as in *Reference files* above.

### 3. Dry run / DAG

```bash
snakemake -n --profile cluster/                 # dry run
snakemake --dag | dot -Tpng > dag.png           # DAG visualization
```

### 4. Run

```bash
snakemake --profile cluster/ --jobs 20
```

Or locally (single node, no scheduler):

```bash
snakemake -j 16 --software-deployment-method conda
```

Resume after a partial run:

```bash
snakemake --profile cluster/ --jobs 20 --rerun-incomplete
```

The workflow automatically:

1. Trims WES reads (`fastp`) and aligns with BWA-MEM2 (piped to `samtools sort`)
2. Marks duplicates and applies BQSR
3. Calls per-donor GVCFs (HaplotypeCaller) and joint-genotypes the cohort
4. Hard-filters SNPs + INDELs, merges, and keeps PASS-only `cohort.pass.vcf.gz`
5. Subsets the cohort VCF to each pool's donors
6. Runs cellSNP-lite on the pooled scRNA-seq BAM at the common-SNP panel
7. Runs Vireo with known donor genotypes (`-N n_donor`) per pool
8. Aggregates `demux_summary.tsv`, the per-cell `cell_donor_sample.tsv`, and the
   unassigned-cell decomposition

---

## Outputs

```text
{output_dir}/
├── wes/joint/
│   ├── cohort.filtered.vcf.gz       # all variants, FILTER column annotated
│   └── cohort.pass.vcf.gz           # PASS-only, fed to Vireo
├── demux/
│   ├── {pool}/
│   │   ├── cellsnp/                 # cellSNP-lite per-cell pileups
│   │   ├── donor_genotype.vcf.gz    # cohort VCF subset to this pool's donors (+ .tbi)
│   │   └── vireo/donor_ids.tsv      # barcode -> donor_id, prob_max, prob_doublet, n_vars
│   ├── demux_summary.tsv            # cohort-wide donor cell counts + labels
│   ├── cell_donor_sample.tsv        # per-cell pool/cell/donor/sample table for Seurat
│   └── qc_unassigned/               # unassigned-cell decomposition
│       ├── unassigned_decomposition_per_pool.tsv
│       ├── unassigned_decomposition_cells.tsv
│       └── unassigned_decomposition_scatter.png
└── logs/                            # per-rule tool stdout/stderr
```

`demux_summary.tsv` columns:

```text
pool   donor_id   mapped_label   n_cells   doublet_rate   unassigned_rate   total_cells
```

`cell_donor_sample.tsv` — one row per cell across all pools, for attaching
labels in R/Seurat. `doublet` / `unassigned` cells carry those literals in both
`donor` and `sample`:

```text
pool   cell   donor   sample
```

Barcodes (`cell`) repeat across pools, so join on **(pool, cell)**:

```r
lab <- read.delim("demux/cell_donor_sample.tsv")
# build a per-pool key matching your Seurat cell names, then merge into meta.data
lab$key <- paste(lab$pool, lab$cell, sep = "_")
rownames(lab) <- lab$key
seu <- AddMetaData(seu, lab[colnames(seu), c("donor", "sample")])
```

---

## Running without an exome targets BED

The pipeline runs end-to-end **without a capture-kit BED** — leave
`reference.exome_targets_bed: ""` empty (the default).

| Rule                | With BED (`-L bed`)               | Without BED                                                         |
| ------------------- | --------------------------------- | ------------------------------------------------------------------ |
| `base_recalibrator` | recalibrates over capture intervals | recalibrates over all aligned reads                                |
| `apply_bqsr`        | applies on capture intervals      | applies on all aligned reads                                       |
| `haplotype_caller`  | calls only inside capture         | iterates the genome, skipping zero-coverage regions via active-region detection |

Other rules are unaffected. Off-target reads produce a few low-coverage calls
outside the kit, but joint hard-filtering removes most, and Vireo only uses
high-confidence common SNPs from cellSNP — so demultiplexing accuracy is not
degraded. Walltimes in `cluster/config.yaml` are set conservatively for the
no-BED case (`haplotype_caller: 36h`, `base_recalibrator/apply_bqsr: 12h`). Set
a BED path later with no other change required.

---

## Unassigned recovery (cellSNP tuning + decomposition)

A high Vireo `unassigned` rate is usually **not** a donor-calling error —
Vireo's singlet-vs-singlet accuracy is essentially perfect (donor-vs-donor
confusion ≈ 0). It reflects Vireo being conservative on low-coverage cells.
Force-assigning those cells lowers the unassigned rate but risks systematic
donor errors, so "fewer unassigned" is not automatically better. Two levers
address this honestly.

### 1. cellSNP coverage relaxation (default)

`config.params.cellsnp` defaults are recovery-oriented for sparse 10x 3′ data:

| Param       | Conservative | Default (this repo) | Why |
|-------------|--------------|---------------------|-----|
| `min_count` | 20           | **2**  | `min_count` is the *aggregate* read/UMI count per SNP across all cells. At 20, low-coverage donor-discriminating panel SNPs are dropped — exactly the SNPs that give a sparse cell its only signal. |
| `min_MAF`   | 0.05         | **0.0** | Vireo receives known donor genotypes, so allele-frequency-in-data filtering is unnecessary and can remove discriminating sites carried by a single donor. |

To A/B test, set the conservative values with a different `output_dir` and
re-run Part 2.

### 2. Unassigned decomposition

`scripts/16_unassigned_decompose.py` (run automatically as part of `rule all`)
splits every Vireo `unassigned` cell into three classes using per-cell `n_vars`
(informative SNPs) and best-singlet posterior `prob_max`:

| Class            | Rule | Meaning |
|------------------|------|---------|
| `data_limited`   | `n_vars < min_vars` | too few SNPs — genuinely unassignable |
| `recoverable`    | `n_vars ≥ min_vars` and `prob_max ≥ recover_prob` | confident best singlet just under Vireo's 0.9 cutoff — relax threshold / coverage to recover |
| `true_ambiguous` | `n_vars ≥ min_vars` and `prob_max < recover_prob` | enough data but signal split — likely doublet / ambient, correctly unassigned |

Thresholds live in `config.params.qc` (`min_vars`, `recover_prob`, `assign_thr`).
It is **non-destructive** — it never edits Vireo output and reports a
`proposed_donor` column for recoverable cells so you decide whether to adopt the
relaxed assignment. Standalone:

```bash
python3 scripts/16_unassigned_decompose.py \
    --demux-dir <output_dir>/demux \
    --pools Pool1 Pool2 \
    --min-vars 10 --recover-prob 0.7
```

If most unassigned cells are `data_limited`, the rate is honest. If many are
`recoverable`, relax cellSNP coverage or the assignment threshold. A large
`true_ambiguous` fraction indicates real doublets/ambient RNA.

---

## Key constraints (enforced)

1. All paths/parameters live in `config.yaml`; rule files hardcode nothing.
2. WES rules use the `{sample}` wildcard, demux rules use `{pool}` — never mixed.
3. The same `genome.fa` is used for WES and scRNA → consistent chromosome names.
4. BWA-MEM2 → samtools sort are connected by a shell pipe (no intermediate SAM).
5. Intermediate BAMs and split SNP/INDEL VCFs are `temp()`; recal BAMs, GVCFs,
   and cohort VCFs are kept. Index (`.tbi`/`.bai`) files are declared as rule
   inputs so they are never deleted before a downstream tool needs them.
6. Vireo always receives a per-pool donor-subset VCF, never the full cohort VCF.
7. Every VCF is tabix-indexed; per-rule conda envs pin all tool versions.

---

## Citation

If you use `scWES2vireo`, please cite this repository and the underlying tools.
At minimum, cite the workflow, `Vireo`, `cellSNP-lite`, `GATK`, the aligner, and
`Snakemake`.

### scWES2vireo

```bibtex
@software{scwes2vireo,
  title        = {scWES2vireo: WES-guided supervised demultiplexing for pooled single-cell RNA-seq with cellSNP-lite and Vireo},
  author       = {Cho, Subin},
  year         = {2026},
  url          = {https://github.com/sbeeniNY/scWES2vireo},
  version      = {v.1}
}
```

### Demultiplexing

- Huang Y, McCarthy DJ, Stegle O. **Vireo: Bayesian demultiplexing of pooled
  single-cell RNA-seq data without genotype reference.** *Genome Biology*.
  2019;20:273. doi: [10.1186/s13059-019-1865-2](https://doi.org/10.1186/s13059-019-1865-2).
  Used for: genotype-based scRNA-seq demultiplexing and doublet detection.

- Huang X, Huang Y. **Cellsnp-lite: an efficient tool for genotyping single
  cells.** *Bioinformatics*. 2021;37(23):4569–4571. doi:
  [10.1093/bioinformatics/btab358](https://doi.org/10.1093/bioinformatics/btab358).
  Used for: per-cell pileup at common SNPs feeding Vireo.

### WES variant calling

- Poplin R, Ruano-Rubio V, DePristo MA, Fennell TJ, Carneiro MO, Van der Auwera
  GA, et al. **Scaling accurate genetic variant discovery to tens of thousands
  of samples.** *bioRxiv*. 2018. doi: [10.1101/201178](https://doi.org/10.1101/201178).
  Used for: GATK `HaplotypeCaller` and joint germline genotyping.

- Van der Auwera GA, O'Connor BD. **Genomics in the Cloud: Using Docker, GATK,
  and WDL in Terra.** 1st ed. O'Reilly Media; 2020.
  Used for: GATK best-practice germline short-variant discovery.

- Li H. **Aligning sequence reads, clone sequences and assembly contigs with
  BWA-MEM.** *arXiv*. 2013. [arXiv:1303.3997](https://arxiv.org/abs/1303.3997).
  Used for: BWA-MEM2 short-read alignment.

- Chen S, Zhou Y, Chen Y, Gu J. **fastp: an ultra-fast all-in-one FASTQ
  preprocessor.** *Bioinformatics*. 2018;34(17):i884–i890. doi:
  [10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560).
  Used for: WES read trimming and adapter removal.

- Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, et al.
  **Twelve years of SAMtools and BCFtools.** *GigaScience*. 2021;10(2):giab008.
  doi: [10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008).
  Used for: `samtools`, `bcftools`, `tabix`, and HTSlib-based VCF/BAM processing.

### Workflow engine and resources

- Köster J, Rahmann S. **Snakemake—a scalable bioinformatics workflow engine.**
  *Bioinformatics*. 2012;28(19):2520–2522. doi:
  [10.1093/bioinformatics/bts480](https://doi.org/10.1093/bioinformatics/bts480).
  Used for: workflow orchestration.

- The 1000 Genomes Project Consortium. **A global reference for human genetic
  variation.** *Nature*. 2015;526:68–74. doi:
  [10.1038/nature15393](https://doi.org/10.1038/nature15393).
  Used for: the common variant resource used by cellSNP-lite.

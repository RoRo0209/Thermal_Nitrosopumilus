# Exploratory population differentiation and reciprocal validation

Set `THERMAL_PROJECT_ROOT` to a separate workspace containing `tajima/`. All data, intermediate files and outputs belong there. This code release includes no read pairs, sample identities/manifests, reference sequences, reference annotations, variant calls, selected-interval tables or observed allele frequencies. The exact original study inputs must be obtained separately.

The available SPOT01 workflow is `01_mapping_and_calling.sh` → `02_merge_populations.sh` → `06_fst_mapping.sh` → `07_fst_analysis.py` → `10_fst_outlier_gene_analysis.py`. `03_popgen_analysis.py` and `07c_pairwise_snv_distance.py` provide additional tables consumed by `figure_integrated_main.py`. `parse_kofam_output.py` parses the externally generated KOfam output. Run each script by its absolute code path; it switches to the workspace module. Linux/WSL is needed for Bash tools. For the legacy shell sources use paths without whitespace.

Required input families:

- `sample_manifest.tsv`: sample identifier, R1 and R2 paths; also `clade` for figure grouping. Respect positional fields used by individual scripts. The targeted caller specifically uses column 1 for sample ID and column 6 for ecotype.
- `refs/HighT_representative.fna`, `refs/LowT_representative.fna`, and matched `refs/LowT_representative.gff`; FASTA, GFF and VCF sequence identifiers must refer to the same coordinate system. Bowtie2 indexes/FAI may be generated in the workspace.
- QC-selected `HighT_samples_passing.txt` and `LowT_samples_passing.txt` from the merge step, and `vcf_shared/HighT_LowT_joint.vcf.gz` from shared-coordinate mapping.
- `LowT_kofam_clean.tsv` with gene IDs and annotation fields as emitted by `parse_kofam_output.py`; optional `candidate_genes.tsv` for overlap reporting. Annotation databases and their exact versions are not bundled.

The V4 SPOT01 scan is exploratory: it counts hard GT calls from pooled metagenomes, keeps sites biallelic within both groups, and therefore excludes sites fixed within either group. It estimates Hudson FST as a ratio of summed components using 10 kb windows and 5 kb steps. Empirical upper-tail windows are rankings, not formal evidence of positive selection. `07c_pairwise_snv_distance.py` measures consensus differences among retained variant sites; it is not genome-wide ANI.

One compatibility repair is explicit: the original `07_fst_analysis.py` called `allel.windowed_hudson_fst`, absent from the official scikit-allel API. This port uses `allel.hudson_fst` followed by `allel.windowed_statistic` sums of numerator and denominator, then division. Site filters and window boundaries are retained. The pairwise panel now handles non-estimable/constant groups rather than crashing a violin plot. Generic species-boundary interpretations and a literal historical FST annotation were removed. The original numerical outputs still require reconciliation with the repaired implementation. See the [official F-statistics API](https://scikit-allel.readthedocs.io/en/stable/stats/fst.html) and [window statistics API](https://scikit-allel.readthedocs.io/en/stable/stats/window.html).

`figure_integrated_main.py` is the available component source for the FST/selected-gene-context panels, not a complete final V4 assembly. It also reads `results/popgen_windowed.tsv`, `results/pairwise_snv_distance.tsv`, `results/fst_outlier_genes_top5.tsv`, the GFF and sample manifest. Its study-selected annotation constants have been externalized rather than copied into the release:

| Required variable | TSV schema |
|---|---|
| `THERMAL_FST_CLUSTERS` | `id`, `start`, `end`, `label`; six ordered intervals, IDs 1 through 6, on the matched reference |
| `THERMAL_FST_GENE_LABELS` | unique `gene_id`, `label` |

The selected intervals/labels need provenance and author-supplied values. They are input data, not newly inferred by the figure script. Final crop/layout and selection provenance remain gaps.

Targeted amoABC code is separate from the genome-wide GT scan. `workstation_amoabc_call_and_export.sh <BAM-directory-list> <reference-FASTA> <manifest> [outdir]` requires `THERMAL_AMO_REGION=contig:start-end` on the matched reference. `amoabc_site_differentiation.py --vcf ... --groups ... --ref ... --gff ... --outdir ...` uses FORMAT/AD when available. It requires `THERMAL_AMO_LOCI`, a TSV with `gene`, `gene_id`, `start`, `end`, `strand` defining amoA, amoB and amoC once each. The reference loci were externalized; do not reuse coordinates on another assembly. The available script is not by itself the complete variant-QC, topology and structural interpretation workflow.

Reciprocal SCM1 validation is available in `SCM1_reciprocal_FST/scripts/01_run_SCM1_mapping_calling.sh` and `02_analyze_SCM1_reciprocal.py`. Required module inputs are `sample_groups_26.tsv` (`sample`, `group`), `reference_manifest.tsv`, `SCM1_amoABC_loci.tsv`, plus the parent sample manifest/reference and the generated `results/SCM1_joint_26samples.vcf.gz`, `SCM1_mapping_QC.tsv`, `SCM1_amoABC_depth_by_sample.tsv`. Inspect each parser for its exact table columns. This analysis uses per-sample read-based allele frequencies and sample counts, and is a distinct estimator from the original GT-based scan. Exact reference and matched-sample membership remain required provenance.

Requirements: Bash, Bowtie2, samtools, FreeBayes, bcftools, tabix, GNU timeout/awk; Python NumPy, pandas, scikit-allel, matplotlib, seaborn, SciPy. No dependency lock from the original execution environment was found. Syntax was checked; reads/VCFs were not rerun, study values were not verified, and the full manuscript cannot be regenerated from this code-only release.

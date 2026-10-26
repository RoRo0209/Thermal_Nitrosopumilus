# Thermal Nitrosopumilus — V4 analysis code

Available analysis sources and reproduction documentation for *Thermal niche differentiation shapes climate response of ammonia-oxidizing archaea in the global ocean*, aligned to the V4 manuscript, figures and supplementary information. This release replaces the earlier documentation-only repository.

**Status: a curated collection of available source code, not a verified end-to-end reproduction package.** Missing original scripts, unavailable same-run outputs and unresolved source matching are documented. The Qingdao 60-week and Kuroshio sample-level data remain restricted and are not distributed. No manuscript files, real study datasets, fitted model objects or database snapshots are included.

## Start here

1. Read the [figure map](docs/FIGURE_MAP.md), which covers seven main and seventeen supplementary figures.
2. Review [data access and gaps](docs/DATA_ACCESS.md) and [known limitations](docs/KNOWN_GAPS.md).
3. Prepare dependencies using the [environment guide](docs/ENVIRONMENT.md). The package lists are installation templates, not recovered historical lockfiles.
4. Create a separate data/output workspace. Most scripts use `THERMAL_PROJECT_ROOT`; the selection and ALE event workflows use explicit path arguments. Follow each module README before running a stage.

```bash
python workflow/preflight.py --code-only
python workflow/preflight.py --workspace /path/to/analysis-workspace --module ALE_2
# The second command inventories requirements and exits 2: it is not an end-to-end readiness certificate.
```

PowerShell workspace configuration:

```powershell
$env:THERMAL_PROJECT_ROOT = 'D:/thermal-analysis-workspace'
```

Linux/WSL configuration:

```bash
export THERMAL_PROJECT_ROOT=/path/to/thermal-analysis-workspace
```

Keep inputs and generated results outside this code checkout. Do not place restricted datasets in a public fork. Repository ignore rules are an additional safeguard, not permission to redistribute inputs.

## Analysis modules

| Module | Scope |
|---|---|
| [abundance_wca_np](analysis/abundance_wca_np/README.md) | NP/WCA composition and global mapping |
| [1.phylogenetic](analysis/1.phylogenetic/README.md) | Predicted traits, phylogenetic signal and available ASR/display sources |
| [genome_character](analysis/genome_character/README.md) | Genome and functional traits, amino-acid composition |
| [60weeks](analysis/60weeks/README.md) | Restricted coastal time series and actual-week lag analyses |
| [kuroshio](analysis/kuroshio/README.md) | Restricted regional gradient analyses |
| [TEMP_experiment](analysis/TEMP_experiment/README.md) | Temperature-incubation comparisons |
| [instrain](analysis/instrain/README.md) | Population metrics and their provenance limitations |
| [tajima](analysis/tajima/README.md) | Differentiation scans, selected regions and reciprocal-reference checks |
| [ALE_2](analysis/ALE_2/README.md) | Same-run ALE node mapping/events and genome repertoires |
| [selection_analysis](analysis/selection_analysis/README.md) | Codon models, convergence checks and V4 Fig. 6c / S12 |
| [SSP585](analysis/SSP585/README.md) | Available climate-modeling, projection and niche sources |
| [ammonia_oxidation_rate_ESSD](analysis/ammonia_oxidation_rate_ESSD/README.md) | Published observation plots, conditional scaling and source-pair audits |

Historical figure numbers in script filenames are retained for traceability; use the V4 figure map to identify the intended panel. Empty files, redundant figure variants, private helper dependencies and unrelated archival scripts were omitted. Missing steps have not been reconstructed and passed off as original code.

The [source manifest](workflow/source_manifest.tsv) records original relative paths and source/release SHA-256 hashes. The [stage inventory](workflow/stages.json) lists inputs, outputs and qualifications. See [release changes](docs/RELEASE_NOTES.md) and [validation scope](docs/VALIDATION.md). The [methods overview](METHODS_OVERVIEW.md) provides conceptual context, not executable instructions.

No license has been selected for this repository. Public visibility does not itself grant a general reuse license; contact the repository owner for permission and citation guidance until a license and publication record are added.

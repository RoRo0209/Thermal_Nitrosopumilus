# ALE events, node identity and accessory repertoires

Only the designated ALE_2 analysis is authoritative for this release. Original same-run `.uml_rec` files were not found there. A derived numbered tree or mapping made from another run cannot independently establish these crown IDs. The supplied species tree has support labels, which must not be interpreted as ALE IDs. Consequently, **V4 Figure 6b and S11 cannot yet be independently verified from the available originals**. No HighT/LowT ALE branch numbers are hardcoded in the released event workflow.

## Same-run event workflow

Supply the full rooted input species tree, every original family `.ale.uml_rec` from that run, a verified `genome_id`/`group` TSV (HighT and LowT), the name/type metadata workbook, family functional annotations and the appropriate KO hierarchy JSON. Use original genome IDs, not display names. Each group must match a complete rooted clade in that full tree. The 108-genome comparison cohort is not interchangeable with the full 177-tip ALE tree.

```bash
Rscript analysis/ALE_2/verify_ale_mapping.R species_tree.nwk one_original.ale.uml_rec full_tree_groups.tsv verified_events
python analysis/ALE_2/extract_ale_events.py --uml-dir original_family_outputs --verified-dir verified_events --expected-families 5403
Rscript analysis/ALE_2/plot_branch_events.R species_tree.nwk one_original.ale.uml_rec metadata.xlsx verified_events figures
Rscript analysis/ALE_2/plot_functional_events.R verified_events family_annotations.xlsx ko_hierarchy.json figures
```

The verifier reads actual internal labels and terminal branch records from the same ALE output. It joins nodes by their complete descendant-tip sets, independently of R's integer node indices and tip ordering. It rejects differing rooted clades, duplicate IDs, incomplete branch records and groups that are not complete crowns. The extractor checks the numbered `S:` tree and complete branch coverage for every family, rejects duplicates/invalid event values, preserves all expected event counts and records input hashes. The drawing scripts consume this verified mapping. For V4 the plot checks expect 177 tips and 5,403 families; different datasets require a separately justified analysis configuration.

`node_mapping.tsv` distinguishes `node` (R plotting index) from `ALE_ID` (original ALE branch identifier). `crowns.tsv` records group membership and the recovered crown ID. Never join an R index directly to the DTLO `Branch` column. Branch counts and labels can change with pruning/rerooting; regenerate the mapping from the unmodified source tree.

## Figure 4e–f: genome repertoires

These calculations do not depend on ALE event numbering. Set `THERMAL_PROJECT_ROOT` to the separate workspace containing `ALE_2/`. Run `full_ko_metabolic_analysis.py`, then `plot_accessory_repertoire.py` by their code paths. Required inputs are:

- `ALE_2/genome_x_family_copy_number_kegg_annotated.tsv`, `gene_family_membership.tsv` and `tree-gen-name-gff.xlsx` (Sheet1).
- `THERMAL_GENOME_GROUPS_JSON`: an external JSON object with `HighT`, `LowT` and `excluded_from_108` arrays of original genome IDs. This replaces the embedded historical cohort lists. Checks require 11 HighT, 23 LowT and exactly 108 retained representative genomes.
- The drawing step also reads `ALE_2/genome_x_family_copy_number.tsv`, the generated `full_ko_analysis/genome_by_KO_copy_number.tsv`, and `genome_character/functional_gene_counts_all_np_processed.csv`.

The KO calculation drops families with ambiguous multiple-KO assignments, aggregates family copies at KO level and retains the source tests/FDR calculations. Family presence and KO copy number are different units. Input schemas are defined by the parsers; no study tables or KEGG database files are redistributed. Full data execution is pending. Synthetic tests of the new node-mapping/extraction guards are separate from biological validation.

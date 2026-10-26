# Environment and external tools

The release includes dependency **templates**, not a recovered historical lockfile. Exact original package/tool/database versions were not available for every stage. Installing current packages does not establish identical historical results.

Python source was syntax-checked with Python 3.12. R source was parsed with R 4.5.2. `environment/python-inspected.json` and `environment/r-inspected.tsv` record packages visible in the local packaging environment; absent packages are not represented as successfully tested. Most scientific Python dependencies were not installed in that runtime. No full analysis environment or dependency-resolution lock has been validated.

Create a separate Python environment, then use `python -m pip install -r environment/requirements.txt` as a starting template. Review package compatibility and record the resolved versions before scientific execution. pandas >=2.2 is required by one population-metric sensitivity script. Legacy scripts may need compatibility review with newly released packages. Plotly static exports also require a compatible Kaleido/Chrome installation and map resources.

For R, `environment/r-packages.txt` lists direct packages found in source. Install CRAN packages in your own library and use a Bioconductor release compatible with your R version for `ggtree`, `ggtreeExtra` and `treeio` when listed. Spatial packages such as terra/sf/ncdf4 require appropriate system libraries. Package managers resolve transitive dependencies. The code does not automatically install software.

Linux/WSL Bash workflows additionally require the relevant stages' external tools: BLAST+, DIAMOND, Bowtie2, samtools, bcftools, FreeBayes, tabix, inStrain, MAFFT, PAL2NAL, trimAl, IQ-TREE and PAML/codeml, plus GNU command-line utilities. Gene reconciliation requires the matching ALEml_undated output provenance. Some upstream methods (marker inference, TOME, MCMCTree and annotation databases) have incomplete execution recipes; installing them alone does not close the code/input gaps.

Use a fresh output workspace when changing references, models or parameters. Some legacy shell scripts reuse existing outputs and assume paths without whitespace. Match reference FASTA/GFF/VCF coordinates and gene-family identities before running. Arial/Helvetica substitution can change figure layout; no byte-identical PDF claim is made.

The synthetic ALE tests need only Python's standard library plus R with ape and jsonlite. Set `RSCRIPT` to the Rscript executable if it is not on PATH, and configure `R_LIBS_USER` if needed:

```bash
python workflow/test_ale_guards.py
```

After acquiring the authorized study inputs, capture `sessionInfo()`, Python environment versions, command-line tool versions, database releases and input hashes for the actual reproduction run. Keep private dataset identifiers outside public version control.

# Validation of this release

Inspection date: 2026-10-26. These checks concern the packaged source and new mapping guards, not a full scientific rerun.

| Check | Result and scope |
|---|---|
| Source inventory | 100 substantive analysis scripts in 12 modules; 13 empty original scripts recorded as unavailable rather than published as executable stages. |
| R syntax | All 53 R source files parsed successfully with R 4.5.2. Packages/data are not loaded by a parse check. |
| Python syntax | All 34 analysis and 2 workflow Python files parsed successfully with Python 3.12. |
| ALE synthetic tests | Eight tests passed: reordered tip labels preserve true ALE IDs and exact event sums; different rooted clades, non-clade group definitions, changed run numbering, missing branches across every family, duplicate records, changed extracted tables and mismatched plotting sources are rejected. Test data are generated in a temporary directory and are explicitly synthetic. |
| Bash | 13 wrappers retained and inspected as source. Bash/WSL was not accessible for a shell syntax/execution check in the packaging environment. Their external bioinformatics tools and real inputs were not run. |
| File integrity | CHECKSUMS.sha256 covers the released text files except itself; source_manifest.tsv separately records original/released analysis hashes. |
| Disclosure review | Staged paths/content were inspected for personal machine paths, credentials and unintended data files. No restricted observations, real sequence/read files, workbooks, fitted models or manuscripts are included. Generic workspace placeholders and public method identifiers remain. |
| Manuscript reproduction | Not performed: restricted inputs, missing original scripts/output provenance and unresolved stage matching remain. Numerical identity, full dependency resolution, model convergence and final figure rendering are not certified. |

To repeat file-integrity and Python syntax checks, run `python workflow/preflight.py --code-only`. To repeat synthetic ALE checks, use `python workflow/test_ale_guards.py` with R, ape and jsonlite available. A checksum mismatch after an intentional edit means the release has changed; it is not automatically a scientific error.

The first local synthetic-test attempt could not access the sandbox's default temporary directory. Re-running with a writable temporary directory completed all eight tests successfully. This was an environment issue; the final test results do not rely on waived assertions. For restricted environments, set TMP/TEMP to an existing writable directory.

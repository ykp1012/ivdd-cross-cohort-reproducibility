# Pre-specified Program Lock

The four program definitions in `config/program_modules.json` were fixed
before any expression result, module score, or differential analysis from the
discovery archives was examined. `data/derived/program_module_ledger.csv`
records the sorted gene lists and SHA-256 hashes.

The score direction is intentionally narrow: a higher score means higher
expression of the listed genes. It does not mean that the program is
beneficial, harmful, causal, or a treatment target. Genes absent from a given
platform are reported as mapping loss; no gene is replaced after looking at a
result. A comparable score requires at least 80% of the locked genes to be
measured, as specified in the analysis protocol.

The disc-matrix module is explicitly a fixed study-defined NP matrix core,
not a validated clinical signature. ECM and inflammatory modules likewise
represent broad transcriptional programs and are not interpreted as single
mechanisms.

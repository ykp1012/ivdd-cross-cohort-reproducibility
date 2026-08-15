# Analysis Environment

## Hardware snapshot

Resource detection was run on 2026-08-13 and saved as `tools/resource_snapshot.json`.

- CPU: 10 physical / 16 logical cores.
- Memory: 15.7 GB total, 7.91 GB available at detection.
- GPU: NVIDIA GeForce RTX 4060 Laptop GPU, 8 GB VRAM.
- Project disk: 328.9 GB free.

The analysis will process cohorts independently and favor sparse, disk-backed intermediate files because the raw archives exceed the available memory budget when expanded.

## Software

- R: 4.4.1 is available, with `MASS`, `Matrix`, `data.table`, `dplyr`, and `ggplot2`. `edgeR`, `limma`, and `DESeq2` are not currently installed.
- The exploratory S7-S9 random-effects synthesis was run with `metafor` 4.8.0,
  `digest` 0.6.37, and `jsonlite` 2.0.0 under R 4.4.1.
- Python: project-local virtual environment at `tools/python/venv`.
- Python package versions are frozen in `tools/python/requirements-lock.txt`.
- Primary single-cell object: AnnData H5AD with sparse matrices.
- The default analysis uses donor/library-level descriptive score summaries and
  does not require a count-model package. S7-S9 use `metafor::escalc` with
  `measure = "SMDH"`, random-effects REML, and Knapp-Hartung intervals. Their
  sensitivity analyses use pooled-SD Hedges *g* and Paule-Mandel tau-squared.
   S7-S9 set REML `maxiter = 10000`; all other `control` parameters retain
  `metafor` defaults. No package installation will be attempted during a
  rerun; a project-local R library is preferred over modifying the system
  library.

## Reproducibility rules

- Analysis scripts set an explicit random seed.
- Every generated result gets a script name, input manifest, and session information.
- Raw data are never overwritten.
- Software package installation outside this project is avoided.

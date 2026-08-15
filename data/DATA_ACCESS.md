# Public data access

The public release contains code, configuration, audit records, derived score
tables, figures, and manuscript files. It deliberately does not redistribute
the original GEO assets or Ensembl reference files. Those files are large,
remain subject to their upstream terms of use, and should be retrieved from the
source records rather than copied into a code release.

`public_data_manifest.tsv` is the source-of-truth acquisition manifest. It
records the expected local path, source URL, byte count, and SHA-256 checksum
for every external file used by the released analysis code.

## Download and verify

Run these commands from the repository root with Python 3.10 or newer:

```powershell
# Frozen default analysis inputs.
python .\scripts\fetch_public_data.py --group default

# Add the isolated exploratory and S8/S9 sensitivity inputs.
python .\scripts\fetch_public_data.py --group exploratory --group s8 --group s9

# Fetch every source file, including candidate-only and Ensembl audit inputs.
python .\scripts\fetch_public_data.py --all

# Verify files already downloaded without using the network.
python .\scripts\fetch_public_data.py --all --verify-only
```

The downloader writes each file beneath `data/raw/`, verifies its exact SHA-256
digest, and never overwrites an existing file. A partial download is retained
and resumed on a later invocation when the source server supports HTTP ranges.

## Release boundary

The GitHub repository and Zenodo snapshot should include the following:

- `config/`, `scripts/`, `docs/`, and `tools/python/requirements-lock.txt`;
- `data/derived/` and `results/`, which contain the small, generated records
  needed to inspect the released findings;
- the formal manuscript and release metadata.

They should exclude `data/raw/`, `tools/python/venv/`, transient logs, Python
bytecode, and regenerated manuscript working files. The raw-data manifest and
downloader preserve a reproducible route to every excluded public input.

## Scope labels

| Group | Contents |
| --- | --- |
| `default` | Frozen 20-effect default analysis inputs. |
| `exploratory` | Isolated GSE251686 sensitivity input. |
| `s8` | Post hoc external-expansion inputs. |
| `s9` | Source-family replacement sensitivity inputs. |
| `candidate` | Candidate-only datasets excluded from the primary synthesis. |
| `reference` | Ensembl release 113 files used only by the candidate probe-specificity audit. |

The scientific interpretation boundaries remain in `README.md`,
`docs/01_data_provenance.md`, and `docs/15_research_quality_audit.md`. The
manifest supports computational reproducibility; it does not change those
boundaries or confer permission to reinterpret the data as confirmatory.

# GSE56081 Probe Annotation Audit

This is a candidate-only extension audit. It does not modify the default IVDD result.

The raw tar contains ten Agilent Feature Extraction files and ten scan images. GPL15314 exposes 60,756 platform rows and probe sequences, but most coding-probe ENSEMBL_ID/ACCESSION_STRING fields are blank. The mapping table therefore records exact sequence matches to Ensembl GRCh38 canonical cDNA, in either orientation. This is a reproducible candidate mapping, not a manufacturer-certified annotation: genome-wide uniqueness, transcript-version compatibility with the 2011 design, and probe summarization rules remain open. A module passes the 80% display gate only when the matched probes are unique within the 78 locked genes; such a pass does not authorize adding GSE56081 to the default summary.

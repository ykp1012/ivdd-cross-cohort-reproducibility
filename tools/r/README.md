# R Environment Specification

The S7-S9 exploratory random-effects syntheses require R 4.4.1 and the
packages listed in requirements.tsv. Run the following command from the
repository root after installing R:

~~~powershell
Rscript .\tools\r\check_required_packages.R
~~~

The checker reports missing or version-mismatched packages and does not modify
the system library. Install the exact recorded versions in an isolated R
library before running the R scripts in scripts/.

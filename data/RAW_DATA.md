# Raw-library provenance

The study reports a 1,012-entry South African Natural Compounds Database (SANCDB) input set, of which 995 structures were parsed and 267 unique warhead-bearing compounds were retained before β-accessibility triage.

The original third-party 1,012-entry raw library file is **not contained in this release archive**. The repository therefore does not claim that `michael_filter.py` can reproduce the 1,012 -> 995 step without the exact source snapshot. To reproduce that stage, retrieve the same SANCDB snapshot used by the authors (Hatherley et al., 2015; Diallo et al., 2021), preserve the downloaded file unmodified, record its SHA-256 checksum and retrieval date, and run `michael_filter.py`.

The processed 267-hit CSV, a derived 267-hit SDF, and the 185-compound accessible CSV are included so that downstream open-source triage and analysis can be rerun.

# General transcriptomics workflow in Python

Reproducible scaffold for sample validation, differential-expression summaries, pathway-level interpretation, and downstream modeling. A dependency-free reference path is provided for continuous verification; the original PyDESeq2/GSEA/scikit-learn research script is retained under `src/legacy/`.

## Quick start

```bash
python workflow.py --input examples/input.csv --config config/workflow.json --output output
python -m unittest discover -s tests -v
```

The synthetic example writes `results.csv`, a checksum-bearing manifest, and a summary figure. See `docs/ENVIRONMENT.md` before running the advanced implementation.

Production analyses must lock genome/reference releases, define contrasts before modeling, keep feature selection inside resampling, correct for multiplicity, and avoid treating technical replicates as independent samples.

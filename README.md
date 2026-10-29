# pvalue-adjuster

**Correct p-values for multiple comparisons with the standard methods (Bonferroni, Sidak, Holm,
Hochberg, Benjamini-Hochberg, Benjamini-Yekutieli, Storey q-values), straight from a CSV.**

Part of the *Research Automation Toolkit*: 100 open tools for daily research work.

## Why

Testing many outcomes, genes, voxels, items or subgroups inflates false positives. Reviewers ask
"did you correct for multiple comparisons?", and the answer should name the method and be
reproducible. Spreadsheet implementations of step-up procedures are notoriously error-prone; these
implementations reproduce R's `p.adjust` exactly.

## Methods

| Method | Controls | Notes |
|--------|----------|-------|
| `bonferroni` | FWER | simplest, conservative |
| `sidak` | FWER | exact for independent tests |
| `holm` | FWER | uniformly more powerful than Bonferroni, no assumptions |
| `hochberg` | FWER | step-up, assumes independence/positive dependence |
| `bh` (`fdr`) | FDR | Benjamini-Hochberg, standard in omics |
| `by` | FDR | Benjamini-Yekutieli, arbitrary dependence |
| `--qvalue` | FDR | Storey q-values with pi0 estimate (lambda = 0.5) |

## Usage

```bash
pip install .
pvalue-adjuster 0.01 0.04 0.03 0.005 0.2 -m holm -m bh
pvalue-adjuster --csv results.csv -c p_value -m bh --qvalue -o results_adjusted.csv
```

## Testing

```bash
python -m unittest discover -s tests -t .
```

## License

MIT (c) JANIN A APURBA

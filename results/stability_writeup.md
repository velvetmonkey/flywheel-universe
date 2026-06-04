# Stability Experiment Writeup

Verdict: **does not hold**.

This is a preliminary, competitive, hardware-aware stability run. It is not a validation or proof.

## Degradation Slopes

```text
                       method     mean     std  count
             hybrid_frobenius  0.01406 0.01789      2
hybrid_support_random_weights  0.01328 0.01900      2
              random_budgeted -0.01789 0.01082      2
                       static -0.03275 0.00538      2
              static_budgeted -0.00676 0.00205      2
```

## Cut Quality Mean +/- Std

```text
  sigma                        method      mean      std  count
0.00000              hybrid_frobenius 718.00000 16.97056      2
0.00000 hybrid_support_random_weights 718.00000 16.97056      2
0.00000               random_budgeted 745.50000  3.53553      2
0.00000                        static 766.00000  0.00000      2
0.00000               static_budgeted 746.00000 11.31371      2
0.25000              hybrid_frobenius 728.00000 12.72792      2
0.25000 hybrid_support_random_weights 728.00000 12.72792      2
0.25000               random_budgeted 749.00000  2.82843      2
0.25000                        static 766.50000  3.53553      2
0.25000               static_budgeted 744.50000  2.12132      2
0.50000              hybrid_frobenius 718.00000  4.24264      2
0.50000 hybrid_support_random_weights 713.00000 11.31371      2
0.50000               random_budgeted 747.50000  7.77817      2
0.50000                        static 759.50000 10.60660      2
0.50000               static_budgeted 747.00000  2.82843      2
1.00000              hybrid_frobenius 731.00000  5.65685      2
1.00000 hybrid_support_random_weights 731.00000  5.65685      2
1.00000               random_budgeted 733.50000 10.60660      2
1.00000                        static 742.50000  6.36396      2
1.00000               static_budgeted 740.50000  9.19239      2
```

## Normalized Degradation Mean +/- Std

```text
  sigma                        method    mean     std  count
0.00000              hybrid_frobenius 1.00000 0.00000      2
0.00000 hybrid_support_random_weights 1.00000 0.00000      2
0.00000               random_budgeted 1.00000 0.00000      2
0.00000                        static 1.00000 0.00000      2
0.00000               static_budgeted 1.00000 0.00000      2
0.25000              hybrid_frobenius 1.01400 0.00624      2
0.25000 hybrid_support_random_weights 1.01400 0.00624      2
0.25000               random_budgeted 1.00470 0.00097      2
0.25000                        static 1.00065 0.00462      2
0.25000               static_budgeted 0.99813 0.01798      2
0.50000              hybrid_frobenius 1.00035 0.02955      2
0.50000 hybrid_support_random_weights 0.99350 0.03924      2
0.50000               random_budgeted 1.00267 0.00568      2
0.50000                        static 0.99151 0.01385      2
0.50000               static_budgeted 1.00143 0.01140      2
1.00000              hybrid_frobenius 1.01830 0.01619      2
1.00000 hybrid_support_random_weights 1.01830 0.01619      2
1.00000               random_budgeted 0.98388 0.00956      2
1.00000                        static 0.96932 0.00831      2
1.00000               static_budgeted 0.99265 0.00273      2
```

## G1 Half-Mean-Degree Falsification

```text
           method        mean      std  count
           greedy 11361.00000 51.02940      3
hebbian_frobenius 11321.00000 44.30576      3
 hybrid_frobenius 11321.00000 44.30576      3
  random_budgeted 11488.00000 21.93171      3
  static_budgeted 11422.00000 56.50664      3
```

## Discipline Notes

- Headline comparison is `hybrid_frobenius` vs `static_budgeted`; raw `static` is context.
- The plotted degradation ratio normalizes each method to its own sigma=0 paired baseline.
- G1 SDP was not run; the ceiling reference is known best cut 11624.
- Results are reported as mean +/- std, not best-of.

## Runtime Caps

- heterogeneity rows 40 of staged diagnostic target 300; headline tables use 40 rows from complete four-sigma paired cells
- full robust 10x10 heterogeneity matrix (2000 rows) was not run
- G14 skipped because it was not present locally

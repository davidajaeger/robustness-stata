# robustness (Stata)

Range tests for equality and equivalence across specifications, from saved bootstrap draws. This is the Stata implementation of the method in Jaeger (2026), "Robustness? Range Tests for Equality and Equivalence Across Specifications." It is the counterpart of the [R `robustness`](https://github.com/davidajaeger/robustness) package: given the same bootstrap draws, the two produce the same statistics.

## What it does

Applied work routinely presents several specifications of the same coefficient and declares the results "robust" when the estimates look similar. The range of the estimates is the implicit object of that claim, but its joint sampling distribution is rarely reported. `robustness` computes, for each comparison set:

- **`R*(.95)` (the minimum equivalence bound):** the smallest tolerance within which the specifications can be certified equivalent at the 5 percent level, in the units of the coefficient. It is the `(1 - alpha)` quantile of the uncentred bootstrap range, which is the paper's definition.
- **`R*(.50)`:** the median of the bootstrap range, a point estimate of the range.
- **`p_R` (the range-based equality test):** the bootstrap p-value for the null that the specifications share a common probability limit.
- the robustness ratio `R*(.95)/|theta_bar|`.

The Wald statistic, its bootstrap p-value, and the Wald-based bound (`W*`) are also computed and returned in `r()`.

The command is a **Stage-2 post-processor: it computes statistics from draws and does not generate them.** It reads its inputs from disk and leaves the data in memory untouched, so it can be run at any point in a session.

## Installation

```
net install robustness, from("https://raw.githubusercontent.com/davidajaeger/robustness-stata/main")
```

To update an existing install to the latest version, add `replace`, then clear the cached program and Mata engine:

```
net install robustness, from("https://raw.githubusercontent.com/davidajaeger/robustness-stata/main") replace
discard
```

(After the paper is published the package will also be available via `ssc install robustness`.)

## Syntax

```
robustness using DRAWSFILE, meta(METAFILE) comps(COMPSFILE) [alpha(numlist) maxdrop(#)]
```

It reads three files, all produced by your bootstrap-generation step:

- **Draws file** (`using`): one row per bootstrap replication, with variables `rep coef1 se1 coef2 se2 ...`, one coef/se pair per specification in column order. Optional `n1 n2 ...` carry per-spec bootstrap sample sizes. The `se` columns are checked for presence but are not used by any statistic, which reads the `coef` columns only. Draws must be **uncentred**; all recentring happens inside the command.
- **Metadata file** (`meta()`): one row per specification, variables `k label theta se`, with `theta` the full-sample estimate. `k` is **required**: it is the specification index, must list each specification exactly once over `1..K`, and the command sorts by it, so `k` fixes which metadata row maps to `coef`*k* / `se`*k* in the draws file. Optional `n` carries the full-sample size.
- **Comparisons file** (`comps()`): one row per comparison, variables `comp_name` and `comp_cols`, where `comp_cols` is a space-separated list of 1-indexed column numbers. A specification may appear at most once in a comparison (duplicates are rejected). `comp_name` is used as an identifier in the output and in `r()`; any label is accepted, but one that is not a valid Stata name is converted with `strtoname` (for example, "All specs" becomes `All_specs`), and the conversion is reported.

`alpha()` sets the significance levels for the equivalence bounds (default `0.50 0.05`); `maxdrop()` sets the maximum percentage of incomplete replications tolerated before the command aborts (default `1`).

## Output

Two panels matching the paper's reporting convention:

- **Panel A** — one row per specification: full-sample estimate, standard error, full-sample n, bootstrap-average n.
- **Panel B** — one row per comparison: `K`, mean estimate, observed range, `R*(.50)`, `R*(.95)`, `p_R`, and the robustness ratio.

Full results, including the Wald statistics and any extra alphas, are returned in `r(table)`, `r(specs)`, and `r(extra)`. The equivalence-bound columns are named `Rstar_50`, `Rstar_95` (and `Wstar_50`, `Wstar_95`), where the suffix is the quantile level `1 - alpha`. See `help robustness`.

## Example

```
robustness using bsdraws.dta, meta(bsdraws_meta.dta) comps(bsdraws_comps.dta)
```

## Generating the draws

The command consumes draws; it does not produce them. The [replication package](https://github.com/davidajaeger/robustness-replications) for Jaeger (2026) contains the generation scripts for each application in the paper. The single requirement the command cannot verify is that the **same resampled units were used for every specification on each replication**. Resampling independently per specification destroys the joint distribution and produces wrong `p_R` and `R*` with no warning. Enforce this in the generation step.

## Citing

The software implements the method; please cite the paper.

> Jaeger, David A. (2026). Robustness? Range Tests for Equality and Equivalence Across Specifications.

## License

MIT.

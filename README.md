# robustness (Stata)

Range tests for equality and equivalence across specifications, from saved bootstrap draws. This is the Stata implementation of the method in Jaeger (2026), "Robustness? Range Tests for Equality and Equivalence Across Specifications." It is the counterpart of the [R `robustness`](https://github.com/davidajaeger/robustness) package: given the same bootstrap draws, the two produce the same statistics.

## What it does

Applied work routinely presents several specifications of the same coefficient and declares the results "robust" when the estimates look similar. The range of the estimates is the implicit object of that claim, but its joint sampling distribution is rarely reported. `robustness` computes, for each comparison set:

- **`R*(.95)` (the minimum equivalence bound):** the smallest tolerance within which the specifications can be certified equivalent at the 5 percent level, in the units of the coefficient. It is the `(1 - alpha)` quantile of the uncentred bootstrap range, which is the paper's definition.
- **`R*(.50)`:** the median of the uncentred bootstrap range distribution.
- **`p_R` (the range-based equality test):** the bootstrap p-value for the null that the specifications share a common probability limit.
- the robustness ratio `R*(.95)/|theta_bar|`.

The Wald statistic, its bootstrap p-value, and the Wald-based bound (`W*`) are also computed and returned in `r()`. The Wald requires a full-rank contrast covariance. If distinct specifications are collinear in the bootstrap draws (for example one is a fixed shift or a linear combination of others) that covariance is singular, and `W`, `p_W`, and `W*` are returned as missing with a warning rather than computed with a generalized inverse. The range statistics (`R`, `p_R`, `R*`) do not use the contrast covariance and are unaffected. Duplicate specification references within a comparison (such as `comp_cols = "1 1 2"`) are a separate case, rejected when the comparisons file is read.

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
robustness using DRAWSFILE, meta(METAFILE) comps(COMPSFILE) [alpha(numlist) maxdrop(#) saving(filename[, replace])]
```

It reads three files, all produced by your bootstrap-generation step:

- **Draws file** (`using`): one row per bootstrap replication. The `coef1 coef2 ...` columns are required, one per specification in column order. Optional `rep`, per-spec `se1 se2 ...`, and per-spec `n1 n2 ...` columns may be present; the `se` columns are ignored (no statistic reads them — the Wald uses the bootstrap covariance of the `coef` draws), and `n`, if present for every spec, drives the average-n reporting in Panel A. Draws must be **uncentred**; all recentring happens inside the command.
- **Metadata file** (`meta()`): one row per specification, variables `k label theta se`, with `theta` the full-sample estimate. `k` is **required**: it is the specification index, must list each specification exactly once over `1..K`, and the command sorts by it, so `k` fixes which metadata row maps to `coef`*k* / `se`*k* in the draws file. Optional `n` carries the full-sample size.
- **Comparisons file** (`comps()`): one row per comparison, variables `comp_name` and `comp_cols`, where `comp_cols` is a space-separated list of 1-indexed column numbers. A specification may appear at most once in a comparison (duplicates are rejected). `comp_name` is used as an identifier in the output and in `r()`; any label is accepted, but one that is not a valid Stata name is converted with `strtoname` (for example, "All specs" becomes `All_specs`), and the conversion is reported.

`alpha()` sets the significance levels for the equivalence bounds (default `0.50 0.05`); `maxdrop()` sets the maximum percentage of incomplete replications tolerated before the command aborts (default `1`).

`saving(filename[, replace])` writes the per-replication bootstrap statistics — the distributions the reported summaries collapse to scalars — to a `.dta` for plotting. The file is long, one row per comparison-by-draw, with `comparison draw range_unc range_rc wald_unc wald_rc`. The `(1 - alpha)` quantile of `range_unc` is `R*`; `p_R` is the Monte Carlo p-value `(1 + #{range_rc >= observed range})/(B + 1)`. The data in memory are left untouched.

## Output

Two panels matching the paper's reporting convention:

- **Panel A** — one row per specification: full-sample estimate, standard error, full-sample n, bootstrap-average n.
- **Panel B** — one row per comparison: `K`, mean estimate, observed range, `R*(.50)`, `R*(.95)`, `p_R`, and the robustness ratio.

Full results, including the Wald statistics and any extra alphas, are returned in `r(table)`, `r(specs)`, and `r(extra)`. The equivalence-bound columns are named `Rstar_50`, `Rstar_95` (and `Wstar_50`, `Wstar_95`), where the suffix is the quantile level `1 - alpha`. See `help robustness`.

## Example

```
robustness using bsdraws.dta, meta(bsdraws_meta.dta) comps(bsdraws_comps.dta)
```

Save and plot the bootstrap range distribution — the object the "robustness" claim implicitly invokes — with the equivalence bounds marked:

```
robustness using bsdraws.dta, meta(bsdraws_meta.dta) comps(bsdraws_comps.dta) saving(rdist.dta, replace)
matrix T = r(table)
use rdist.dta, clear
histogram range_unc if comparison=="main", xline(`=T["main","Rstar_50"]' `=T["main","Rstar_95"]')
```

The example reads the bounds from `r(table)` rather than recomputing them from the saved draws on purpose: `R*` is a type-1 order statistic (no interpolation), so a line placed with `summarize, detail` percentiles, which interpolate, can sit a fraction off the reported bound.

## Generating the draws

The command consumes draws; it does not produce them. The [replication package](https://github.com/davidajaeger/robustness-replications) for Jaeger (2026) contains the generation scripts for each application in the paper. The single requirement the command cannot verify is that the **same resampled units were used for every specification on each replication**. Resampling independently per specification destroys the joint distribution and produces wrong `p_R` and `R*` with no warning. Enforce this in the generation step.

## Citing

The software implements the method; please cite the paper.

> Jaeger, David A. (2026). Robustness? Range Tests for Equality and Equivalence Across Specifications.

## License

MIT.

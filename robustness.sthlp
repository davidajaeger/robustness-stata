{smcl}
{* *! version 1.4.1  Jaeger (2026)}{...}
{viewerjumpto "Syntax" "robustness##syntax"}{...}
{viewerjumpto "Description" "robustness##description"}{...}
{viewerjumpto "Options" "robustness##options"}{...}
{viewerjumpto "Input files" "robustness##inputs"}{...}
{viewerjumpto "Output" "robustness##output"}{...}
{viewerjumpto "Remarks" "robustness##remarks"}{...}
{viewerjumpto "Examples" "robustness##examples"}{...}
{viewerjumpto "Stored results" "robustness##results"}{...}
{viewerjumpto "References" "robustness##references"}{...}
{title:Title}

{phang}
{bf:robustness} {hline 2} Range tests for equality and equivalence across
specifications, from saved bootstrap draws

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:robustness}
{cmd:using} {it:drawsfile}{cmd:,}
{opth m:eta(filename)}
{opth c:omps(filename)}
[{opth a:lpha(numlist)}
{opt maxd:rop(#)}
{opt sav:ing(filename[, replace])}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt :{opth m:eta(filename)}}metadata file with full-sample estimates{p_end}
{synopt :{opth c:omps(filename)}}comparisons file defining the spec groups{p_end}

{syntab:Optional}
{synopt :{opth a:lpha(numlist)}}significance levels for equivalence margins.
Default {cmd:alpha(0.50 0.05)}, matching Panel B in the paper. Additional
alphas are computed for every comparison and returned in {cmd:r(extra)} but are
not added to Panel B{p_end}
{synopt :{opt maxd:rop(#)}}maximum percent of incomplete replications before
aborting. Default {cmd:maxdrop(1)}{p_end}
{synopt :{opt sav:ing(filename[, replace])}}save the per-replication bootstrap
statistics to a {cmd:.dta} for plotting{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:robustness} computes the range-test statistics in Jaeger (2026) from
bootstrap draws saved by an application's generation step. Output is two
panels matching the paper's reporting convention:

{phang}
{bf:Panel A: Specification-level estimates.} One row per specification, with
columns for the full-sample point estimate, standard error, full-sample n,
and bootstrap-average n.

{phang}
{bf:Panel B: Comparison-set statistics.} One row per comparison, with columns
for K, mean point estimate, observed range, R*(.50), R*(.95),
bootstrap p-value p_R, and the robustness ratio R*(.95)/|theta_bar|.

{pstd}
Equivalence bounds are computed at alpha = .50 and .05, matching Panel B in
Jaeger (2026). The Wald statistic, its bootstrap p-value, and the
Wald-based equivalence bound are computed and available in {cmd:r()} but are
not printed.

{pstd}
The command reads all of its inputs from disk. It preserves the data in
memory and restores it on exit, so it can be run at any point in a session
without disturbing the user's data.

{marker options}{...}
{title:Options}

{phang}
{opth meta(filename)} names the metadata file. Required. See {help robustness##inputs:Input files}.

{phang}
{opth comps(filename)} names the comparisons file. Required.

{phang}
{opth alpha(numlist)} controls the equivalence-bound computation. Panel B
always reports R*(.50) and R*(.95), matching the reporting convention of
Jaeger (2026); these two columns are not user-tunable. The {opt alpha()}
option lets you compute additional equivalence bounds at other significance
levels. The additional bounds are returned in the matrix {cmd:r(extra)} for
every comparison, with columns {cmd:Rstar_}{it:XX} and {cmd:Wstar_}{it:XX}
for each extra alpha (where {it:XX} is the quantile level 1-alpha as two
digits; for example, alpha=0.10 -> {cmd:Rstar_90}). The additional bounds
are not added to Panel B's printed output. If {opt alpha()} is not supplied,
only the default .50 and .05 are computed. Each value must be strictly
between 0 and 1.

{phang}
{opt maxdrop(#)} sets the maximum percentage of bootstrap replications that
may be incomplete before the command aborts. A replication is incomplete if
any coefficient draw in the comparison is non-finite, typically because a
specification failed to converge on that resample. Incomplete replications
are dropped and counted. If the dropped share exceeds {opt maxdrop()} the
command stops, because the surviving draws may no longer represent the
intended distribution. The default is {cmd:1}.

{pstd}
{opt saving(filename[, replace])} writes the per-replication bootstrap
statistics to {it:filename} (a {cmd:.dta}), the distributions that the
reported summaries collapse to scalars. The file is in long form, one row per
comparison-by-draw, with a string variable {cmd:comparison} and the numeric
variables {cmd:draw}, {cmd:range_unc}, {cmd:range_rc}, {cmd:wald_unc},
{cmd:wald_rc}. The uncentred series ({cmd:range_unc}, {cmd:wald_unc}) are the
ones whose {cmd:(1-alpha)} quantiles are {cmd:R*} and {cmd:W*}; the recentred
series ({cmd:range_rc}, {cmd:wald_rc}) are the ones whose tails at or above the
observed statistic give {cmd:p_R} and {cmd:p_W}. Draws are renumbered 1 to B
within each comparison after dropping incomplete replications, so the same
{cmd:draw} number does not correspond to the same resample across comparisons.
Specify {cmd:replace} to overwrite an existing file. The data in memory are
left untouched. See the examples below for plotting.

{pstd}
When marking {cmd:R*} on a plot of the saved draws, read the bound from
{cmd:r(table)} rather than recomputing a quantile from the saved series.
{cmd:R*} is a type-1 order statistic (no interpolation), whereas
{cmd:summarize, detail} and most quantile routines interpolate, so a line
computed that way can sit a fraction off the reported bound. The plotting
examples below take the bound from {cmd:r(table)} for exactly this reason.

{marker inputs}{...}
{title:Input files}

{pstd}
{bf:Draws file} ({cmd:using}). One observation per bootstrap replication.
The {cmd:coef1} {cmd:coef2} ... columns are required, one per specification
in column order. Optional {cmd:rep} and per-spec {cmd:se1} {cmd:se2} ...
columns may be present but are ignored: no statistic reads them (the Wald
uses the bootstrap covariance of the {cmd:coef} draws). Optional {cmd:n1}
{cmd:n2} ... carry per-spec bootstrap sample sizes; if present for every
spec, the command reports the bootstrap-average n per spec in Panel A. The
coefficient draws must be raw, that is uncentred. All recentring happens
inside the command.

{pstd}
{bf:Metadata file} ({opt meta()}). One observation per specification.
Variables {cmd:k}, {cmd:label}, {cmd:theta}, {cmd:se}, where {cmd:k} is the
specification index and {cmd:theta} is the full-sample point estimate.
{cmd:k} is required and must list each specification exactly once over
1 to K. The command sorts by {cmd:k}, so the rows need not be supplied in
order, but {cmd:k} fixes the mapping to the draws file: specification
{cmd:k} pairs with {cmd:coef}{it:k} and {cmd:se}{it:k}. Optional {cmd:n}
carries the full-sample sample size, shown in Panel A if present.

{pstd}
{bf:Comparisons file} ({opt comps()}). One observation per comparison.
Variables {cmd:comp_name} and {cmd:comp_cols}, where {cmd:comp_cols} is a
space-separated list of 1-indexed column numbers. A comparison needs at
least two specifications, every column number must be an integer in the
range 1 to K, and a specification may appear at most once in a comparison
(duplicates are rejected). {cmd:comp_name} is used as an identifier in the
output table and in {cmd:r()}; any label is accepted, but one that is not a
valid Stata name is converted with {help strtoname:strtoname()} (for example
"All specs" becomes {cmd:All_specs}), and the command reports the conversion.

{marker output}{...}
{title:Output}

{pstd}
Panel B's robustness ratio R*(.95)/|theta_bar| is the heuristic ratio
defined in Jaeger (2026), Section 4. When |theta_bar| is close to zero,
interpret it with caution: the ratio can be large even when R*(.95) is
economically small. In such cases, judge the bound directly in coefficient
units rather than as a ratio.

{pstd}
A note is printed below Panel B if any comparison has bootstrap-average
sample sizes that differ across its specifications. The per-spec n values
appear in Panel A.

{marker remarks}{...}
{title:Remarks}

{pstd}
The single requirement the command cannot verify is that the same resampled
units were used for all specifications on each replication. Resampling
independently per specification destroys the joint distribution across
specifications and produces wrong p_R and R* with no error and no
warning. The guarantee must be enforced in the generation step. The
generation examples accompanying Jaeger (2026) enforce and teach it.

{pstd}
The Wald statistics require a full-rank contrast covariance. When the
specifications are collinear in the bootstrap draws (for example one is a fixed
shift or a linear combination of others) that covariance is singular, and the
command reports {cmd:W}, {cmd:p_W}, and {cmd:W*} as missing with a warning
rather than inverting it with a generalized inverse. The range statistics
({cmd:R}, {cmd:p_R}, {cmd:R*}) do not use the contrast covariance and are
reported normally. Exact duplicate columns are a separate case and are rejected
when the draws are read.

{marker examples}{...}
{title:Examples}

{pstd}Run all comparisons at the default alphas{p_end}
{phang2}{cmd:. robustness using bsdraws.dta, meta(bsdraws_meta.dta) comps(bsdraws_comps.dta)}{p_end}

{pstd}Request additional alphas. The bounds are computed for all
comparisons and returned in {cmd:r(extra)}.{p_end}
{phang2}{cmd:. robustness using bsdraws.dta, meta(bsdraws_meta.dta) comps(bsdraws_comps.dta) alpha(0.50 0.10 0.05)}{p_end}

{pstd}Save the per-replication statistics, then plot the equivalence-bound
distribution for one comparison with the median and 95th-percentile bounds
marked. The bounds are read from {cmd:r(table)}, which the command leaves
behind.{p_end}
{phang2}{cmd:. robustness using bsdraws.dta, meta(bsdraws_meta.dta) comps(bsdraws_comps.dta) saving(rdist.dta, replace)}{p_end}
{phang2}{cmd:. matrix T = r(table)}{p_end}
{phang2}{cmd:. use rdist.dta, clear}{p_end}
{phang2}{cmd:. histogram range_unc if comparison=="main", xline(`=T["main","Rstar_50"]' `=T["main","Rstar_95"]')}{p_end}

{pstd}Or visualize the equality test: the recentred range distribution with
the observed range marked. The mass at or beyond the line is {cmd:p_R}.{p_end}
{phang2}{cmd:. histogram range_rc if comparison=="main", xline(`=T["main","R"]')}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:robustness} stores the following in {cmd:r()}.

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Scalars}{p_end}
{synopt:{cmd:r(nspecs)}}number of specifications K{p_end}
{synopt:{cmd:r(ncomps)}}number of comparisons{p_end}
{synopt:{cmd:r(B)}}number of bootstrap replications{p_end}

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Macros}{p_end}
{synopt:{cmd:r(comparisons)}}names of the comparisons computed{p_end}

{synoptset 16 tabbed}{...}
{p2col 5 16 20 2: Matrices}{p_end}
{synopt:{cmd:r(specs)}}Nspecs x 4 matrix of Panel A data. Rows are the spec
labels from the metadata; columns are {cmd:theta}, {cmd:se}, {cmd:n_full},
{cmd:n_boot}{p_end}
{synopt:{cmd:r(table)}}Ncomps x 12 matrix of Panel B data. Rows are the
comparison names; columns are {cmd:theta_bar}, {cmd:R}, {cmd:p_R}, {cmd:W},
{cmd:p_W}, {cmd:Rstar_50}, {cmd:Rstar_95}, {cmd:Wstar_50},
{cmd:Wstar_95}, {cmd:ratio}, {cmd:K}, {cmd:dropped}{p_end}
{synopt:{cmd:r(extra)}}Ncomps x (2 * n_extras) matrix, present only when
{cmd:alpha()} requests significance levels beyond .50 and .05. Rows are the
comparison names; columns are {cmd:Rstar_}{it:XX} and
{cmd:Wstar_}{it:XX} for each extra alpha, where {it:XX} is the quantile level
1-alpha as two digits (alpha=0.10 -> {cmd:Rstar_90}){p_end}

{marker references}{...}
{title:References}

{phang}
Jaeger, D. A. 2026. Robustness? Range Tests for Equality and Equivalence
Across Specifications. Working paper.

{title:Author}

{pstd}David A. Jaeger, University of St Andrews.{p_end}

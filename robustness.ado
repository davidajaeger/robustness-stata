*! version 1.4.1  Jaeger (2026)  Range Tests for Equality and Equivalence
*! robustness -- compute robustness statistics from saved bootstrap draws
program define robustness, rclass
    version 16.1

    /*=======================================================================
      robustness -- canonical Stage 2 post-processor for the robustness
      statistics in Jaeger (2026), "Robustness? Range Tests for Equality
      and Equivalence Across Specifications."

      Reads three files produced by an application's bootstrap-generation
      step and computes, for each comparison, the equality test (R, p_R,
      W, p_W) and the equivalence bounds (R* and W*) at one or
      more significance levels.

      The command reads its inputs from disk and does not touch the data in
      memory. The single requirement it cannot verify is that the same
      resampled units were used for all specifications on each replication.
      That guarantee lives in the generation step.

      Syntax
      ------
      robustness using DRAWSFILE, Meta(string) Comps(string)
              [ Alpha(numlist) MAXDrop(real 1) SAVing(string) ]

        using DRAWSFILE   B-row file of raw (uncentered) bootstrap draws. The
                          coef1..coefK columns are required, one per spec in
                          column order. Optional rep, per-spec se1 se2 ... and
                          per-spec sample sizes n1 n2 ... may be included; the
                          se columns are ignored, and if n is present for every
                          spec the command reports average bootstrap n per spec
                          in Panel A. None of these enter any statistic.
        meta()            K-row metadata file, variables k label theta se;
                          optional n column carries the full-sample sample
                          size and is shown in Panel A if present.
        comps()           comparisons file, variables comp_name comp_cols
        alpha()           significance levels for equivalence margins.
                          Default 0.50 0.05, matching Panel B in the paper.
                          Additional alphas are computed for all comparisons
                          and returned in r(extra) as a matrix; they are
                          not added to Panel B's printed output.
        maxdrop()         max percent of incomplete reps before aborting
                          (default 1)
        saving()          saving(filename [, replace]) writes the per-
                          replication bootstrap statistics to a .dta, long,
                          one row per comparison-by-draw, with variables
                          comparison draw range_unc range_rc wald_unc wald_rc.
                          These are the distributions the summaries collapse:
                          the (1-alpha) quantile of range_unc is R*, the share
                          of range_rc at or above the observed range is p_R.
                          Intended for plotting. The data in memory are left
                          untouched.
    =======================================================================*/

    syntax using/, Meta(string) Comps(string) ///
        [ Alpha(numlist >0 <1 sort) MAXDrop(real 1) SAVing(string) ]

    * The Panel B layout in Jaeger (2026) reports R* at alpha = .50 and .05
    * (that is, R*_{.50} and R*_{.95}). Those two are always computed and
    * shown. Users may request additional alphas via alpha(); the additional
    * bounds are computed for every comparison and returned in r(extra), not
    * added to Panel B.
    if "`alpha'" == "" local alpha "0.50 0.05"
    local extra_alphas ""
    foreach a of local alpha {
        if `a' != 0.50 & `a' != 0.05 {
            local extra_alphas "`extra_alphas' `a'"
        }
    }
    local extra_alphas = strtrim("`extra_alphas'")
    if `maxdrop' < 0 | `maxdrop' > 100 {
        di as error "maxdrop() must be between 0 and 100."
        exit 198
    }

    * Optional saving(filename [, replace]): write the per-replication bootstrap
    * statistics (the distributions whose quantiles give R* and whose tails give
    * p_R) to a .dta for plotting. Parse the filename and the lone permitted
    * suboption here, before the data are touched.
    local svfile ""
    local svreplace ""
    if `"`saving'"' != "" {
        gettoken svfile svtail : saving, parse(",")
        local svfile = strtrim(`"`svfile'"')
        local svtail = strtrim(`"`svtail'"')
        if `"`svtail'"' != "" {
            gettoken comma svtail : svtail, parse(",")
            local svtail = strtrim(`"`svtail'"')
            if `"`svtail'"' != "replace" {
                di as error "saving() allows only the optional suboption 'replace'."
                exit 198
            }
            local svreplace "replace"
        }
    }

    * Confirm all three files exist before touching the user's data.
    capture confirm file "`using'"
    if _rc {
        di as error "draws file not found: `using'"
        exit 601
    }
    capture confirm file "`meta'"
    if _rc {
        di as error "meta() file not found: `meta'"
        exit 601
    }
    capture confirm file "`comps'"
    if _rc {
        di as error "comps() file not found: `comps'"
        exit 601
    }

    * Preserve the user's data. Everything below reads from files with plain
    * use, and restore returns the user's data untouched on exit, including
    * on error. This avoids any dependence on frame local-scoping rules.
    preserve

    * --------------------------------------------------------------------
    * 1. Metadata: validate, enforce spec ordering via k, build theta_hat
    * --------------------------------------------------------------------
    quietly use "`meta'", clear

    * Required variables. k is the specification index and is enforced: it
    * binds metadata row r to the draws-file column coef`r'.
    foreach v in k label theta se {
        capture confirm variable `v'
        if _rc {
            di as error "Metadata file is missing required variable '`v''."
            di as error "Required variables: k label theta se."
            restore
            exit 111
        }
    }
    qui count
    local nspecs = r(N)
    if `nspecs' < 2 {
        di as error "Need at least 2 specifications, found `nspecs'."
        restore
        exit 2000
    }

    * k must be a permutation of 1..nspecs: numeric, nonmissing, integer, no
    * duplicates, with min 1 and max nspecs. For a file of nspecs rows those
    * conditions force k to be exactly {1,...,nspecs}. We then sort by k so
    * that row position equals specification number and the position-based
    * read below is correct.
    capture confirm numeric variable k
    if _rc {
        di as error "Metadata variable 'k' must be numeric (the specification index 1..`nspecs')."
        restore
        exit 109
    }
    qui count if missing(k) | k != floor(k)
    if r(N) {
        di as error "Metadata 'k' must be a nonmissing integer in every row."
        restore
        exit 198
    }
    tempvar ktag
    qui duplicates tag k, generate(`ktag')
    qui count if `ktag' > 0
    if r(N) {
        di as error "Metadata 'k' has duplicate values; each specification must appear exactly once."
        restore
        exit 198
    }
    qui summarize k, meanonly
    if r(min) != 1 | r(max) != `nspecs' {
        di as error "Metadata 'k' must span 1..`nspecs' with no gaps (found min `=r(min)', max `=r(max)')."
        restore
        exit 198
    }
    sort k

    * Full-sample n per spec is optional in the metadata. Shown if present.
    capture confirm variable n
    local meta_has_n = (_rc == 0)
    di _n as text "Metadata: `meta'"
    di "  Number of specifications: `nspecs'"
    matrix theta_hat   = J(`nspecs', 1, .)
    matrix se_hat      = J(`nspecs', 1, .)
    matrix n_full_hat  = J(`nspecs', 1, .)
    forvalues k = 1/`nspecs' {
        local lab = label[`k']
        matrix theta_hat[`k', 1] = theta[`k']
        matrix se_hat[`k', 1]    = se[`k']
        if `meta_has_n' matrix n_full_hat[`k', 1] = n[`k']
        * Each spec label gets its own local. The engine reads label_k for
        * k = 1..nspecs. Avoids fragile string-parsing on labels that may
        * themselves contain spaces or other delimiters.
        local label_`k' "`lab'"
    }

    * --------------------------------------------------------------------
    * 2. Comparisons: validate column references
    * --------------------------------------------------------------------
    quietly use "`comps'", clear
    foreach v in comp_name comp_cols {
        capture confirm variable `v'
        if _rc {
            di as error "Comparisons file is missing required variable '`v''."
            di as error "Required variables: comp_name comp_cols."
            restore
            exit 111
        }
    }
    qui count
    local ncomps = r(N)
    if `ncomps' < 1 {
        di as error "No comparisons defined."
        restore
        exit 2000
    }
    di _n as text "Comparisons: `comps'"
    di "  Number of comparisons: `ncomps'"

    local comparisons ""
    local used_safe ""
    forvalues c = 1/`ncomps' {
        local cname = comp_name[`c']
        local ccols = comp_cols[`c']
        if "`cname'" == "" {
            di as error "Comparison `c' has an empty name."
            restore
            exit 2000
        }
        if "`ccols'" == "" {
            di as error "Comparison '`cname'' has no columns."
            restore
            exit 2000
        }
        local m = 0
        foreach col of local ccols {
            capture confirm integer number `col'
            if _rc {
                di as error "Comparison '`cname'' references non-integer column '`col''."
                restore
                exit 2000
            }
            if `col' < 1 | `col' > `nspecs' {
                di as error "Comparison '`cname'' references column `col', outside 1..`nspecs'."
                restore
                exit 2000
            }
            local m = `m' + 1
        }

        * Reject duplicate columns. Each specification may appear at most once;
        * a repeat would overstate K and make the bootstrap covariance
        * singular. list dups returns the repeated tokens; uniq trims the
        * message to one mention of each. Token-safe: "1" and "12" never match.
        local cdups : list dups ccols
        if "`cdups'" != "" {
            local cdups : list uniq cdups
            di as error "Comparison '`cname'' lists specification(s) `cdups' more than once; each spec may appear at most once."
            restore
            exit 198
        }

        if `m' < 2 {
            di as error "Comparison '`cname'' has fewer than 2 specifications."
            restore
            exit 2000
        }
        * Comparison names are used as identifiers: as local-macro keys here
        * and as matrix row names in r(table). They must therefore be valid
        * Stata names. strtoname() converts any label (spaces, hyphens, and so
        * on) into one; collisions are disambiguated with a numeric suffix. The
        * safe name is what appears in the printed table and in r(); a note is
        * printed whenever it differs from the label the user supplied.
        local safe = strtoname("`cname'")
        local j = 1
        local pos : list posof "`safe'" in used_safe
        while `pos' > 0 {
            local j = `j' + 1
            local safe = strtoname("`cname'") + "_`j'"
            local pos : list posof "`safe'" in used_safe
        }
        local used_safe "`used_safe' `safe'"
        if "`safe'" != "`cname'" {
            di as text "  (comparison '`cname'' recorded as '`safe'')"
        }
        local comp_`safe' "`ccols'"
        local comparisons "`comparisons' `safe'"
        di "  `safe': specs `ccols'"
    }
    local comparisons = strtrim("`comparisons'")

    * --------------------------------------------------------------------
    * 3. Draws: verify the column layout, then run the engine
    * --------------------------------------------------------------------
    quietly use "`using'", clear
    local coef_vars ""
    local n_vars ""
    local has_n 1
    forvalues k = 1/`nspecs' {
        capture confirm variable coef`k'
        if _rc {
            di as error "Draws file is missing variable 'coef`k''."
            di as error "Expected coef1..coef`nspecs' to match the `nspecs' specs in metadata."
            restore
            exit 111
        }
        local coef_vars "`coef_vars' coef`k'"
        * Per-spec se`k' and sample size n`k' are optional and are not read by
        * any statistic (the Wald uses the bootstrap covariance of the coef
        * draws). They are ignored if present. If every spec has n`k', the
        * command reports average n per spec; otherwise n reporting is skipped.
        capture confirm variable n`k'
        if _rc local has_n 0
        else   local n_vars "`n_vars' n`k'"
    }
    if !`has_n' local n_vars ""
    qui count
    local B = r(N)
    if `B' < `nspecs' {
        di as error "Only `B' draws for `nspecs' specs. Need at least `nspecs'."
        restore
        exit 2001
    }
    di _n as text "Draws: `using'"
    di "  Bootstrap replications: `B'"

    * Clear any staging matrices left over from a previous run, so the
    * conditional surfacing of __rob_extra later doesn't pick up stale data.
    capture matrix drop __rob_specs
    capture matrix drop __rob_table
    capture matrix drop __rob_extra

    * The draws are the data in memory now, so st_view in the engine reads
    * them directly. theta_hat is a Stata matrix, visible to Mata. The
    * locals coef_vars, comparisons, comp_*, alpha, maxdrop, nspecs are all
    * in this program's scope.
    mata: _robustness_engine()

    * If requested, write the per-replication statistics to a .dta. The draws
    * are still in memory here, so _rob_save_draws() copies what it needs, then
    * clears and rebuilds memory as the output dataset; the restore below brings
    * the user's data back regardless.
    if `"`svfile'"' != "" {
        mata: _rob_save_draws()
        quietly save `"`svfile'"', `svreplace'
        di _n as text "Per-replication statistics saved to " as result `"`svfile'"' ///
            as text " (`=_N' rows: comparison, draw, range_unc, range_rc, wald_unc, wald_rc)."
    }

    * Restore the user's data.
    restore

    * Run-level returns: scalars and macro.
    return local comparisons "`comparisons'"
    return scalar ncomps = `ncomps'
    return scalar nspecs = `nspecs'
    return scalar B      = `B'

    * Matrix returns: the engine staged these as named Stata matrices.
    * Surface them through return matrix, then drop the staging matrices.
    return matrix specs = __rob_specs
    return matrix table = __rob_table
    capture confirm matrix __rob_extra
    if _rc == 0 {
        return matrix extra = __rob_extra
    }
end


* ==========================================================================
* Mata engine. Defined once, persists for the session.
* ==========================================================================
version 16.1
mata:
mata set matastrict on

// Type-1 sample quantile: the inverse of the empirical c.d.f., i.e. the
// smallest order statistic x_(j) with at least a share p of the data <= x_(j),
// which is x_(ceil(n*p)). Matches R quantile(type=1) and is the paper's
// definition of R*_{1-alpha}. Drops missing before computing.
real scalar _rob_quantile(real colvector v, real scalar p) {
    real colvector sv
    real scalar    n, idx
    sv  = sort(select(v, v :< .), 1)
    n   = rows(sv)
    if (n == 0) return(.)
    idx = ceil(n * p)
    if (idx < 1) idx = 1
    if (idx > n) idx = n
    return(sv[idx])
}

// Grand-mean contrast matrix, (K-1) x K.
real matrix _rob_contrast(real scalar K) {
    real matrix R
    real scalar i
    R = J(K-1, K, -1/K)
    for (i=1; i<=K-1; i++) R[i,i] = R[i,i] + 1
    return(R)
}

// Parse space-separated integers to a row vector.
real rowvector _rob_str2cols(string scalar s) {
    string rowvector tok
    real rowvector v
    real scalar i
    tok = tokens(s)
    v   = J(1, cols(tok), .)
    for (i=1; i<=cols(tok); i++) v[1,i] = strtoreal(tok[i])
    return(v)
}

// True if the symmetric PSD contrast covariance is numerically rank deficient,
// so the Wald is undefined. A relative eigenvalue tolerance is used rather than
// rank(), because a default rank tolerance can miss structured collinearity
// (e.g. one specification a constant shift of another) whose true zero
// eigenvalue floating-point rounding lifts to a tiny positive value. The gap is
// wide in practice: full-rank cases sit near rcond ~ 1e-1, deficient ones near
// 1e-16.
real scalar _rob_rank_deficient(real matrix RVR) {
    real rowvector ev
    ev = symeigenvalues(RVR)
    if (max(ev) <= 0) return(1)
    return(min(ev) <= 1e-12 * max(ev))
}

// Per-comparison statistics. Reads the uncentred draws two ways: uncentred
// for the equivalence bounds, recentred for the equality p-values. Returns a
// 12-element column vector in the canonical r(table) order:
//
//   1  theta_bar   mean of the estimates in the comparison
//   2  R           observed range, max(theta) - min(theta)
//   3  p_R         range equality p-value (recentred bootstrap)
//   4  W           observed Wald statistic
//   5  p_W         Wald equality p-value (recentred bootstrap)
//   6  Rstar_50    R*_{.50}, median of the uncentred bootstrap range
//   7  Rstar_95    R*_{.95}, the minimum equivalence bound
//   8  Wstar_50    W*_{.50}, sqrt of the .50 uncentred Wald quantile
//   9  Wstar_95    W*_{.95}, Wald-scale equivalence bound
//   10 ratio       Rstar_95 / |theta_bar| (missing if theta_bar == 0)
//   11 K           number of specifications in the comparison
//   12 dropped     incomplete bootstrap replications for this comparison
//
// Aborts if too few complete replications remain or the dropped share
// exceeds maxdrop.
real colvector _rob_compute(real matrix theta, real matrix DRAWS,
                            real rowvector cols, real scalar maxdrop,
                            string scalar label) {
    real scalar    B, K, B_orig, dropped, dropshare, theta_bar, W_obs, R_obs
    real scalar    p_R, p_W, ratio, Wstar_50, Wstar_95
    real matrix    th, D, Rmat, RVR, RVRinv, Bd, Bdc, Dc
    real colvector R_unc, W_unc, R_rc, W_rc, out

    th = theta[cols', 1]            // K x 1 full-sample estimates
    D  = DRAWS[., cols]             // B x K uncentred bootstrap draws
    K  = cols(D)

    // Keep only replications complete across these specs; count and police.
    B_orig    = rows(D)
    D         = select(D, rowmissing(D) :== 0)
    B         = rows(D)
    dropped   = B_orig - B
    dropshare = dropped / B_orig
    if (B < K) {
        printf("\n%s: only %g complete reps of %g, need at least %g.\n",
               label, B, B_orig, K)
        printf("  ABORTING. Too few usable replications.\n")
        exit(error(2001))
    }
    if (dropshare > maxdrop) {
        printf("\n%s: %g of %g reps incomplete (%5.2f%%), exceeds %5.2f%% limit.\n",
               label, dropped, B_orig, 100*dropshare, 100*maxdrop)
        printf("  ABORTING. The surviving draws may not represent the intended distribution.\n")
        printf("  Check the generation step. A common cause is a specification failing\n")
        printf("  to converge on many resamples.\n")
        exit(error(2001))
    }

    // Grand-mean contrast and the bootstrap covariance of the K contrasts.
    Rmat = _rob_contrast(K)
    RVR  = Rmat * variance(D) * Rmat'

    theta_bar = mean(th)

    // RANGE statistics. These never use the contrast covariance, so they are
    // always defined.
    R_obs = max(th) - min(th)
    R_unc = rowmax(D) :- rowmin(D)
    // Recentred draws: subtract each spec's deviation from the cross-spec mean
    // (theta_hat_k - theta_bar), imposing Delta = 0. The common +theta_bar
    // shift cancels in the range and in the contrast, so it does not affect the
    // p-value; it is kept to mirror the paper's definition.
    Dc    = D :- th' :+ theta_bar
    R_rc  = rowmax(Dc) :- rowmin(Dc)
    p_R   = (1 + sum(R_rc :>= R_obs)) / (B + 1)

    // WALD statistics. They require a full-rank contrast covariance. Duplicate
    // columns are rejected at parse time, so a singular RVR here means the
    // specifications are genuinely collinear in the bootstrap draws. Rather
    // than invert with a generalized inverse (which would return a degenerate
    // Wald), report W, p_W, and W* as missing and keep the range results.
    if (_rob_rank_deficient(RVR)) {
        W_obs = .; p_W = .; Wstar_50 = .; Wstar_95 = .
        printf("\n%s: contrast covariance is rank deficient (collinear specifications).\n", label)
        printf("  Wald statistics (W, p_W, W*) set to missing; range statistics (R, p_R, R*) are unaffected.\n")
    }
    else {
        RVRinv   = invsym(RVR)
        W_obs    = ((Rmat*th)' * RVRinv * (Rmat*th))[1,1]
        Bd       = Rmat * D'
        W_unc    = colsum(Bd :* (RVRinv * Bd))'
        Bdc      = Rmat * Dc'
        W_rc     = colsum(Bdc :* (RVRinv * Bdc))'
        // Monte Carlo p-value (1 + #)/(B + 1): the observed statistic joins its
        // own reference set, so it is bounded away from zero and uniform under
        // the null by exchangeability (Davison and Hinkley 1997).
        p_W      = (1 + sum(W_rc :>= W_obs)) / (B + 1)
        Wstar_50 = sqrt(_rob_quantile(W_unc, 0.50))
        Wstar_95 = sqrt(_rob_quantile(W_unc, 0.95))
    }

    // Robustness ratio uses R*_{.95} (range-based, always defined).
    if (abs(theta_bar) > 0) ratio = _rob_quantile(R_unc, 0.95) / abs(theta_bar)
    else                    ratio = .

    out      = J(12, 1, .)
    out[ 1]  = theta_bar
    out[ 2]  = R_obs
    out[ 3]  = p_R
    out[ 4]  = W_obs
    out[ 5]  = p_W
    out[ 6]  = _rob_quantile(R_unc, 0.50)          // R*_{.50}
    out[ 7]  = _rob_quantile(R_unc, 0.95)          // R*_{.95}
    out[ 8]  = Wstar_50                            // W*_{.50}
    out[ 9]  = Wstar_95                            // W*_{.95}
    out[10]  = ratio
    out[11]  = K
    out[12]  = dropped
    return(out)
}

// Build the per-replication statistics dataset for saving(). For every
// comparison it recomputes the four bootstrap series that the summary path
// produces internally -- the uncentred range and Wald (whose quantiles are
// R* and W*) and the recentred range and Wald (whose tails give p_R and p_W)
// -- and stacks them long, one row per (comparison, complete draw). The
// arithmetic mirrors _rob_compute exactly; this routine only exposes the
// per-draw series rather than collapsing them to quantiles.
//
// Called while the draws are still in memory. It copies the coef columns with
// st_data (not a view), then clears memory and rebuilds it as the output
// dataset; the caller saves it and restores the user's data.
void _rob_save_draws() {
    real matrix      theta, DRAWS, D, Dc, Rmat, RVR, RVRinv, Bd, Bdc, M, block
    real colvector   keeprows, ru, rr, wu, wr, th
    real rowvector   cols
    real scalar      nc, c, K, Bc, theta_bar, mx, i
    string rowvector coefvars, compnames
    string colvector names, nblock

    theta     = st_matrix("theta_hat")
    coefvars  = tokens(st_local("coef_vars"))
    DRAWS     = st_data(., coefvars)
    compnames = tokens(st_local("comparisons"))
    nc        = cols(compnames)

    M     = J(0, 5, .)
    names = J(0, 1, "")

    for (c = 1; c <= nc; c++) {
        cols = _rob_str2cols(st_local("comp_" + compnames[c]))
        th   = theta[cols', 1]

        D        = DRAWS[., cols]
        keeprows = rowmissing(D) :== 0
        D        = select(D, keeprows)
        Bc       = rows(D)
        if (Bc == 0) continue
        K        = cols(D)

        Rmat      = _rob_contrast(K)
        RVR       = Rmat * variance(D) * Rmat'
        theta_bar = mean(th)

        ru  = rowmax(D) :- rowmin(D)
        Dc  = D :- th' :+ theta_bar
        rr  = rowmax(Dc) :- rowmin(Dc)

        // Wald columns only when the contrast covariance is full rank; under
        // rank deficiency they are written missing, matching W/p_W/W*.
        if (!_rob_rank_deficient(RVR)) {
            RVRinv = invsym(RVR)
            Bd     = Rmat * D'
            wu     = colsum(Bd :* (RVRinv * Bd))'
            Bdc    = Rmat * Dc'
            wr     = colsum(Bdc :* (RVRinv * Bdc))'
        }
        else {
            wu = J(Bc, 1, .)
            wr = J(Bc, 1, .)
        }

        block  = ((1::Bc), ru, rr, wu, wr)
        M      = M \ block
        nblock = J(Bc, 1, compnames[c])
        names  = names \ nblock
    }

    // Longest comparison name, for the string variable width.
    mx = 1
    for (i = 1; i <= rows(names); i++) mx = max((mx, strlen(names[i])))

    stata("clear")
    if (rows(M) == 0) return
    st_addobs(rows(M))
    (void) st_addvar("str" + strofreal(mx), "comparison")
    (void) st_addvar("double", ("draw", "range_unc", "range_rc", "wald_unc", "wald_rc"))
    st_sstore(., "comparison", names)
    st_store(., ("draw", "range_unc", "range_rc", "wald_unc", "wald_rc"), M)
}

// Compute average per-spec sample size across bootstrap reps. Returns a
// column vector of length K (avg n for each spec in the comparison).
// Descriptive only; does not enter any statistic.
real colvector _rob_avg_n_comp(real matrix N, real rowvector cols) {
    real matrix    Nc
    real colvector avg, col_k
    real scalar    K, j
    Nc = N[., cols]
    K  = cols(Nc)
    avg = J(K, 1, .)
    for (j=1; j<=K; j++) {
        col_k = select(Nc[., j], Nc[., j] :< .)
        if (rows(col_k) > 0) avg[j] = mean(col_k)
    }
    return(avg)
}

// Right-justify a string in a field of given width. Used to align headers
// directly above the column positions produced by printf format specs.
string scalar _rob_rpad(string scalar s, real scalar w) {
    real scalar pad, i
    string scalar out
    pad = w - strlen(s)
    if (pad <= 0) return(s)
    out = ""
    for (i=1; i<=pad; i++) out = out + " "
    return(out + s)
}

// Left-justify a string in a field of given width.
string scalar _rob_lpad(string scalar s, real scalar w) {
    real scalar pad, i
    string scalar out
    pad = w - strlen(s)
    out = s
    if (pad <= 0) return(out)
    for (i=1; i<=pad; i++) out = out + " "
    return(out)
}

// Print Panel A: per-specification estimates.
// Columns: spec label, theta_hat, SE, n (full sample), avg boot n.
// Last two are shown if available; if neither is, just theta and SE.
// Headers are right-justified above the data positions produced by the
// corresponding %f format specs, so they align column-by-column.
void _rob_print_panel_a(string rowvector labels, real matrix theta_hat,
                        real matrix se_hat, real matrix n_full,
                        real colvector n_boot_all, real scalar has_n_full,
                        real scalar has_n_boot, real scalar divwidth) {
    real scalar K, k, maxlab, labw
    string scalar header

    // Column widths matching the data-row format specs.
    real scalar w_theta, w_se, w_n
    w_theta = 12         // %12.5f
    w_se    = 8          // %8.5f
    w_n     = 10         // %10.0f

    K = rows(theta_hat)

    // Width the spec label column to the longest label.
    maxlab = 4   // minimum "Spec"
    for (k=1; k<=K; k++) {
        if (strlen(labels[k]) > maxlab) maxlab = strlen(labels[k])
    }
    labw = maxlab + 2

    printf("\n%s\n", "-"*divwidth)
    printf("  Panel A: Specification-level estimates\n")
    printf("%s\n", "-"*divwidth)

    // Header row: "Spec" left-justified in the label column, then each
    // numeric header right-justified to align with its data column.
    header = "  " + _rob_lpad("Spec", labw)
    header = header + _rob_rpad("theta_hat", w_theta)
    header = header + "  " + _rob_rpad("SE", w_se)
    if (has_n_full) header = header + "    " + _rob_rpad("n_full",     w_n)
    if (has_n_boot) header = header + "    " + _rob_rpad("avg_n_boot", w_n)
    printf("%s\n", header)

    for (k=1; k<=K; k++) {
        printf("  %s", _rob_lpad(labels[k], labw))
        printf("%12.5f  %8.5f", theta_hat[k, 1], se_hat[k, 1])
        if (has_n_full) printf("    %10.0f", n_full[k, 1])
        if (has_n_boot) printf("    %10.0f", n_boot_all[k])
        printf("\n")
    }
}

// Print Panel B: comparison-set statistics.
// Reads from the canonical results matrix (12 rows, columns indexed by
// the column order documented in _rob_compute). The printed columns are:
//   Comparison set, K, theta_bar, R(theta), R*(.50), R*(.95), p_R, Rob. ratio
// Headers right-justified to align with the data positions produced by the
// %f format specs.
void _rob_print_panel_b(string rowvector cnames, real matrix results,
                        real scalar maxlab, real scalar divwidth) {
    real scalar    nc, c, labw, K_c
    real scalar    theta_bar, R_obs, dR50, dR05, pR, ratio
    string scalar  cn, ratio_str, header

    // Column widths matching the data-row format specs.
    real scalar w_K, w_theta, w_R, w_d50, w_d05, w_pR, w_ratio
    w_K     = 3          // %3.0f
    w_theta = 9          // %9.5f
    w_R     = 9          // %9.5f
    w_d50   = 9          // %9.5f
    w_d05   = 9          // %9.5f
    w_pR    = 7          // %7.4f
    w_ratio = 10         // %10.4f

    nc   = cols(cnames)
    labw = maxlab + 2
    if (labw < 16) labw = 16    // min header "Comparison set"

    printf("\n%s\n", "-"*divwidth)
    printf("  Panel B: Comparison-set statistics\n")
    printf("%s\n", "-"*divwidth)

    // Header row.
    header = "  " + _rob_lpad("Comparison set", labw)
    header = header + " "    + _rob_rpad("K",          w_K)
    header = header + "   "  + _rob_rpad("theta_bar",  w_theta)
    header = header + "   "  + _rob_rpad("R(theta)",   w_R)
    header = header + "   "  + _rob_rpad("R*(.50)",    w_d50)
    header = header + "   "  + _rob_rpad("R*(.95)",    w_d05)
    header = header + "  "   + _rob_rpad("p_R",        w_pR)
    header = header + "       " + _rob_rpad("Rob. ratio", w_ratio)
    printf("%s\n", header)

    // Data rows. Column indices match the order in _rob_compute:
    //   1 theta_bar, 2 R, 3 p_R, 6 Rstar_50, 7 Rstar_95, 10 ratio, 11 K
    for (c=1; c<=nc; c++) {
        cn        = cnames[c]
        theta_bar = results[ 1, c]
        R_obs     = results[ 2, c]
        pR        = results[ 3, c]
        dR50      = results[ 6, c]
        dR05      = results[ 7, c]
        ratio     = results[10, c]
        K_c       = results[11, c]

        if (ratio < .) ratio_str = sprintf("%10.4f", ratio)
        else           ratio_str = "         ."

        printf("  %s", _rob_lpad(cn, labw))
        printf(" %3.0f   %9.5f   %9.5f   %9.5f   %9.5f  %7.4f       %s\n",
               K_c, theta_bar, R_obs, dR50, dR05, pR, ratio_str)
    }
}

// Print a one-line note below Panel B about a comparison's dropped reps.
// B (the total bootstrap reps) is passed in; the per-comparison results
// matrix no longer carries it since it's a run-level constant.
void _rob_print_drop_notes(string rowvector cnames, real matrix results,
                           real scalar B) {
    real scalar nc, c, dropped
    nc = cols(cnames)
    for (c=1; c<=nc; c++) {
        dropped = results[12, c]
        if (dropped > 0) {
            printf("  Note: comparison '%s' dropped %g of %g reps as incomplete.\n",
                   cnames[c], dropped, B)
        }
    }
}

// Print sample-variation note if any comparison's specs differ in avg n.
// We compute the per-comparison avg-n range here and print the worst case.
void _rob_print_n_note(real matrix N, string rowvector cnames,
                       string scalar comp_prefix) {
    real scalar    nc, c, varying, range_min, range_max
    real colvector avg
    string scalar  compspec
    string rowvector tok
    real rowvector cols
    real scalar    i

    nc = cols(cnames)
    varying = 0
    for (c=1; c<=nc; c++) {
        compspec = st_local(comp_prefix + cnames[c])
        tok = tokens(compspec)
        cols = J(1, cols(tok), .)
        for (i=1; i<=cols(tok); i++) cols[1,i] = strtoreal(tok[i])
        avg = _rob_avg_n_comp(N, cols)
        range_min = min(avg)
        range_max = max(avg)
        if (range_max > range_min) {
            varying = 1
            break
        }
    }
    if (varying) {
        printf("  Note: average bootstrap sample size differs across specifications in at least one comparison.\n")
        printf("        See Panel A for per-spec sample sizes.\n")
    }
}

// Extra-alpha bounds for every comparison. Returns an Ncomps x (2*n_extras)
// matrix, columns alternating Rstar, Wstar for each extra alpha. Uses the
// uncentred range and Wald, the same objects as the main path.
real matrix _rob_compute_extras(real matrix theta, real matrix DRAWS,
                                string scalar comparisons,
                                string scalar comp_prefix,
                                string scalar extralist) {
    real scalar      nc, c, na, a, K, alpha, wald_ok
    real matrix      D, Rmat, RVR, RVRinv, Bd, out_extra
    real colvector   R_unc, W_unc
    real rowvector   cols
    string rowvector compnames, alphas

    compnames = tokens(comparisons)
    nc        = cols(compnames)
    alphas    = tokens(extralist)
    na        = cols(alphas)
    out_extra = J(nc, 2 * na, .)

    for (c=1; c<=nc; c++) {
        cols   = _rob_str2cols(st_local(comp_prefix + compnames[c]))
        D      = select(DRAWS[., cols], rowmissing(DRAWS[., cols]) :== 0)
        K      = cols(D)
        Rmat   = _rob_contrast(K)
        RVR    = Rmat * variance(D) * Rmat'
        R_unc  = rowmax(D) :- rowmin(D)
        // Wald bounds only when the contrast covariance is full rank; otherwise
        // the Wstar columns stay missing (out_extra is initialized to .).
        wald_ok = !_rob_rank_deficient(RVR)
        if (wald_ok) {
            RVRinv = invsym(RVR)
            Bd     = Rmat * D'
            W_unc  = colsum(Bd :* (RVRinv * Bd))'
        }
        for (a=1; a<=na; a++) {
            alpha = strtoreal(alphas[a])
            out_extra[c, 2*(a-1) + 1] = _rob_quantile(R_unc, 1 - alpha)
            if (wald_ok) {
                out_extra[c, 2*(a-1) + 2] = sqrt(_rob_quantile(W_unc, 1 - alpha))
            }
        }
    }
    return(out_extra)
}

// Column names for r(extra): Rstar_XX, Wstar_XX, where XX is the quantile
// level 1-alpha as two digits (alpha=0.10 -> Rstar_90, Wstar_90).
string rowvector _rob_extras_colnames(string scalar extralist) {
    real scalar      na, a, tag, alpha
    string rowvector alphas, names
    string scalar    suffix
    alphas = tokens(extralist)
    na     = cols(alphas)
    names  = J(1, 2 * na, "")
    for (a=1; a<=na; a++) {
        alpha  = strtoreal(alphas[a])
        tag    = round((1 - alpha) * 100)
        suffix = (tag < 10 ? "0" : "") + strofreal(tag)
        names[2*(a-1) + 1] = "Rstar_" + suffix
        names[2*(a-1) + 2] = "Wstar_" + suffix
    }
    return(names)
}

// Driver: pulls locals set by the ado, runs every comparison.
void _robustness_engine() {
    real matrix      DRAWS, theta, se_hat, n_full, N, results
    string scalar    alphalist, compname, compspec, nvars
    string rowvector comps, labels
    real scalar      c, maxdrop, has_n_boot, has_n_full, nspecs, B, k, maxlab
    real colvector   n_boot_avg, comp_stats, col_k

    st_view(DRAWS=., ., tokens(st_local("coef_vars")))
    theta     = st_matrix("theta_hat")
    se_hat    = st_matrix("se_hat")
    n_full    = st_matrix("n_full_hat")
    alphalist = st_local("alpha")
    maxdrop   = strtoreal(st_local("maxdrop")) / 100
    nspecs    = rows(theta)
    B         = rows(DRAWS)

    // Optional per-spec bootstrap sample sizes. Descriptive only.
    nvars = st_local("n_vars")
    has_n_boot = (strtrim(nvars) != "")
    if (has_n_boot) st_view(N=., ., tokens(nvars))
    else            N = J(0, 0, .)

    // Full-sample n is available if any value is nonmissing.
    has_n_full = (sum(n_full :< .) > 0)

    // Parse labels: each is in a separate local label_k.
    labels = J(1, nspecs, "")
    for (k=1; k<=nspecs; k++) {
        labels[k] = st_local("label_" + strofreal(k))
    }

    // Compute per-spec bootstrap-mean n (descriptive only) across all reps.
    n_boot_avg = J(nspecs, 1, .)
    if (has_n_boot) {
        for (k=1; k<=nspecs; k++) {
            col_k = select(N[., k], N[., k] :< .)
            if (rows(col_k) > 0) n_boot_avg[k] = mean(col_k)
        }
    }

    // Compute label-column widths and row widths for both panels, then
    // use the wider as the divider for both (uniform within application).
    real scalar maxlab_a, maxlab_b, labw_a, labw_b, w_a, w_b, divwidth
    real matrix specs_mat, extras_mat
    string rowvector specs_colnames, table_colnames, extras_colnames

    maxlab_a = 4   // "Spec"
    for (k=1; k<=nspecs; k++) {
        if (strlen(labels[k]) > maxlab_a) maxlab_a = strlen(labels[k])
    }
    labw_a = maxlab_a + 2
    // Panel A data row width: 2 ("  ") + labw_a + 12 (theta_hat) + 2 + 8 (SE)
    //                         + 14 (n_full, if shown) + 14 (n_boot, if shown)
    w_a = 2 + labw_a + 22 + (has_n_full ? 14 : 0) + (has_n_boot ? 14 : 0)

    // Header banner.
    printf("\n%s\n", "="*70)
    printf("  Robustness Statistics\n")
    printf("  Jaeger (2026), Range Tests for Equality and Equivalence\n")
    printf("                 Across Specifications\n")
    printf("  alpha = %s   B = %s   K = %s   max drop = %s%%\n",
           alphalist, strofreal(B), strofreal(nspecs),
           st_local("maxdrop"))
    printf("%s\n", "="*70)

    // Compute every comparison. The canonical results matrix is 12 rows by
    // Ncomps columns internally (for efficient column access during
    // printing). Column order matches _rob_compute's documented layout.
    comps = tokens(st_local("comparisons"))
    results = J(12, cols(comps), .)
    maxlab_b = strlen("Comparison set")  // min header width
    for (c=1; c<=cols(comps); c++) {
        compname = comps[1, c]
        if (strlen(compname) > maxlab_b) maxlab_b = strlen(compname)
        compspec = st_local("comp_" + compname)
        comp_stats = _rob_compute(theta, DRAWS, _rob_str2cols(compspec),
                                  maxdrop, compname)
        for (k=1; k<=12; k++) results[k, c] = comp_stats[k]
    }
    labw_b = maxlab_b + 2
    if (labw_b < 16) labw_b = 16
    // Panel B data row width: 2 ("  ") + labw_b + format spec width 77.
    w_b = 2 + labw_b + 77

    // Unified divider: the wider of the two panels.
    divwidth = max((w_a, w_b))

    // Panel A printing.
    _rob_print_panel_a(labels, theta, se_hat, n_full, n_boot_avg,
                       has_n_full, has_n_boot, divwidth)

    // Panel B printing.
    _rob_print_panel_b(comps, results, maxlab_b, divwidth)
    _rob_print_drop_notes(comps, results, B)
    if (has_n_boot) _rob_print_n_note(N, comps, "comp_")
    printf("  Note: the robustness ratio is R*(.95)/|theta_bar|. When |theta_bar| is close to zero,\n")
    printf("        interpret it with caution.\n")

    // --- Stage matrices for the ado to return via return matrix ---
    //
    // We can't write to r() directly here because the rclass ado will
    // clear and rebuild r() on the return statements after this engine
    // exits. So we stage everything as named Stata matrices and let the
    // ado use 'return matrix' to publish them.

    // __rob_specs: Panel A data. Nspecs x 4, rows = spec labels,
    // cols = (theta, se, n_full, n_boot).
    specs_mat = (theta, se_hat, n_full, n_boot_avg)
    st_matrix("__rob_specs", specs_mat)
    st_matrixrowstripe("__rob_specs", (J(nspecs, 1, ""), labels'))
    specs_colnames = ("theta", "se", "n_full", "n_boot")
    st_matrixcolstripe("__rob_specs", (J(4, 1, ""), specs_colnames'))

    // __rob_table: Panel B data, transposed to comparison-per-row.
    st_matrix("__rob_table", results')
    st_matrixrowstripe("__rob_table", (J(cols(comps), 1, ""), comps'))
    table_colnames = ("theta_bar", "R", "p_R", "W", "p_W",
                      "Rstar_50", "Rstar_95", "Wstar_50", "Wstar_95",
                      "ratio", "K", "dropped")
    st_matrixcolstripe("__rob_table", (J(12, 1, ""), table_colnames'))

    // __rob_extra: only if extras requested.
    if (st_local("extra_alphas") != "") {
        extras_mat      = _rob_compute_extras(theta, DRAWS,
                                              st_local("comparisons"), "comp_",
                                              st_local("extra_alphas"))
        extras_colnames = _rob_extras_colnames(st_local("extra_alphas"))
        st_matrix("__rob_extra", extras_mat)
        st_matrixrowstripe("__rob_extra", (J(cols(comps), 1, ""), comps'))
        st_matrixcolstripe("__rob_extra",
                           (J(cols(extras_colnames), 1, ""), extras_colnames'))
        printf("\n  Additional bounds at alpha = %s in r(extra).\n",
               st_local("extra_alphas"))
        printf("  Type 'matrix list r(extra)' to display.\n")
    }

    printf("\n  Results returned as r(specs) and r(table).\n")
    printf("  Type 'matrix list r(table)' to display.\n")

    printf("\n%s\n", "="*70)
    printf("  Done.\n")
    printf("%s\n", "="*70)
}
end

*! version 1.0.0  Jaeger (2026)  Range Tests for Equality and Equivalence
*! robustness -- compute robustness statistics from saved bootstrap draws
program define robustness, rclass
    version 16.1

    /*=======================================================================
      robustness -- canonical Stage 2 post-processor for the robustness
      statistics in Jaeger (2026), "Robustness? Range Tests for Equality
      and Equivalence Across Specifications."

      Reads three files produced by an application's bootstrap-generation
      step and computes, for each comparison, the equality test (R, p_R,
      W, p_W) and the equivalence margins (delta*_R and delta*_W) at one or
      more significance levels.

      The command reads its inputs from disk and does not touch the data in
      memory. The single requirement it cannot verify is that the same
      resampled units were used for all specifications on each replication.
      That guarantee lives in the generation step.

      Syntax
      ------
      robustness using DRAWSFILE, Meta(string) Comps(string)
              [ Alpha(numlist) MAXDrop(real 1) ]

        using DRAWSFILE   B-row file of raw bootstrap draws, with variables
                          rep coef1 se1 coef2 se2 ... (uncentred draws)
                          optional per-spec sample sizes n1 n2 ... may be
                          included; if present for every spec, the command
                          reports average bootstrap n per spec in Panel A.
                          They are descriptive and do not enter any statistic.
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
    =======================================================================*/

    syntax using/, Meta(string) Comps(string) ///
        [ Alpha(numlist >0 <1 sort) MAXDrop(real 1) ]

    * The Panel B layout in Jaeger (2026) reports delta* at alpha = .50 and
    * .05. Those two are always computed and shown. Users may request
    * additional alphas via alpha(); the additional bounds are computed for
    * the last comparison and stored in r() but not added to Panel B.
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
    * 1. Metadata: validate, build theta_hat
    * --------------------------------------------------------------------
    quietly use "`meta'", clear
    foreach v in label theta se {
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
        if `m' < 2 {
            di as error "Comparison '`cname'' has fewer than 2 specifications."
            restore
            exit 2000
        }
        local comp_`cname' "`ccols'"
        local comparisons "`comparisons' `cname'"
        di "  `cname': specs `ccols'"
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
        capture confirm variable se`k'
        if _rc {
            di as error "Draws file is missing variable 'se`k''."
            restore
            exit 111
        }
        local coef_vars "`coef_vars' coef`k'"
        * Per-spec sample size n`k' is optional. It is descriptive only and
        * does not enter any statistic. If every spec has it, the command
        * reports average n per spec. If any is absent, n reporting is
        * skipped and the coef/se statistics are unaffected.
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

// Quantile, R type 7 (h = p*(n-1)+1). Drops missing before computing.
real scalar _rob_quantile(real colvector v, real scalar p) {
    real colvector sv
    real scalar    n, h
    sv = select(v, v :< .)
    n  = rows(sv)
    if (n == 0) return(.)
    sv = sort(sv, 1)
    h  = p * (n - 1) + 1
    if (h <= 1) return(sv[1])
    if (h >= n) return(sv[n])
    return(sv[floor(h)] + (h - floor(h)) * (sv[ceil(h)] - sv[floor(h)]))
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

// Compute per-comparison statistics. Pure function: no side effects, no
// scratch scalars. Returns a 12-element column vector in the canonical
// column order used by r(table):
//
//   1  theta_bar    (mean of theta_hat across comparison specs)
//   2  R            (observed range)
//   3  p_R          (range-test bootstrap p-value)
//   4  W            (observed Wald statistic)
//   5  p_W          (Wald-test bootstrap p-value)
//   6  delta_R_50   (range-based equivalence bound, alpha=.50)
//   7  delta_R_05   (range-based equivalence bound, alpha=.05)
//   8  delta_W_50   (Wald-based equivalence bound, alpha=.50)
//   9  delta_W_05   (Wald-based equivalence bound, alpha=.05)
//   10 ratio        (delta_R_05 / |theta_bar|, missing if theta_bar==0)
//   11 K            (number of specs in the comparison)
//   12 dropped      (number of incomplete bootstrap reps for this comparison)
//
// Aborts on too many incomplete reps; max_drop enforces the bound.
real colvector _rob_compute(real matrix theta, real matrix DRAWS,
                            real rowvector cols, real scalar maxdrop,
                            string scalar label) {
    real scalar    B, K, i, B_orig, dropped, dropshare
    real matrix    th, D, Vhat, Rmat, RVR, RVRinv, d, d_b, d_b_rc
    real colvector W_boot, R_boot, W_boot_rc, R_boot_rc, _ok, out
    real scalar    W_obs, R_obs, theta_bar, pW, pR, ratio
    real scalar    q_R_50, q_R_05, q_W_50, q_W_05
    real scalar    delta_R_50, delta_R_05, delta_W_50, delta_W_05

    th = theta[cols', 1]
    D  = DRAWS[., cols]
    K  = cols(D)

    B_orig    = rows(D)
    _ok       = rowmissing(D) :== 0
    D         = select(D, _ok)
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

    Vhat   = variance(D)
    Rmat   = _rob_contrast(K)
    RVR    = Rmat * Vhat * Rmat'
    RVRinv = invsym(RVR)

    d         = Rmat * th
    W_obs     = (d' * RVRinv * d)[1,1]
    R_obs     = max(th) - min(th)
    theta_bar = mean(th)

    W_boot    = J(B, 1, .)
    R_boot    = J(B, 1, .)
    W_boot_rc = J(B, 1, .)
    R_boot_rc = J(B, 1, .)
    for (i=1; i<=B; i++) {
        d_b          = Rmat * D[i,.]'
        W_boot[i]    = (d_b' * RVRinv * d_b)[1,1]
        R_boot[i]    = max(D[i,.]) - min(D[i,.])
        d_b_rc       = Rmat * (D[i,.]' :- th :+ theta_bar)
        W_boot_rc[i] = (d_b_rc' * RVRinv * d_b_rc)[1,1]
        R_boot_rc[i] = max(D[i,.] :- th' :+ theta_bar) ///
                      - min(D[i,.] :- th' :+ theta_bar)
    }
    pW = mean(W_boot_rc :>= W_obs)
    pR = mean(R_boot_rc :>= R_obs)

    q_R_50 = _rob_quantile(R_boot :- R_obs, 0.50)
    q_R_05 = _rob_quantile(R_boot :- R_obs, 0.95)
    q_W_50 = _rob_quantile(W_boot :- W_obs, 0.50)
    q_W_05 = _rob_quantile(W_boot :- W_obs, 0.95)

    delta_R_50 = R_obs + q_R_50
    delta_R_05 = R_obs + q_R_05
    delta_W_50 = sqrt(max((W_obs + q_W_50, 0)))
    delta_W_05 = sqrt(max((W_obs + q_W_05, 0)))

    if (abs(theta_bar) > 0) ratio = delta_R_05 / abs(theta_bar)
    else                    ratio = .

    out = J(12, 1, .)
    out[ 1] = theta_bar
    out[ 2] = R_obs
    out[ 3] = pR
    out[ 4] = W_obs
    out[ 5] = pW
    out[ 6] = delta_R_50
    out[ 7] = delta_R_05
    out[ 8] = delta_W_50
    out[ 9] = delta_W_05
    out[10] = ratio
    out[11] = K
    out[12] = dropped
    return(out)
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
//   Comparison set, K, theta_bar, R(theta), d*(.50), d*(.05), p_R, Rob. ratio
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
    header = header + "   "  + _rob_rpad("d*(.50)",    w_d50)
    header = header + "   "  + _rob_rpad("d*(.05)",    w_d05)
    header = header + "  "   + _rob_rpad("p_R",        w_pR)
    header = header + "       " + _rob_rpad("Rob. ratio", w_ratio)
    printf("%s\n", header)

    // Data rows. Column indices match the order in _rob_compute:
    //   1 theta_bar, 2 R, 3 p_R, 6 delta_R_50, 7 delta_R_05, 10 ratio, 11 K
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
        else           ratio_str = "       inf"

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

// Compute extra-alpha bounds for the given comparison. Recomputes the
// bootstrap range and Wald distributions (cheap), then stores each
// alpha's delta*_R and delta*_W in Stata scalars __rob_delta_R_aXX and
// __rob_delta_W_aXX where XX is the two-digit alpha (e.g. alpha=.10 -> a10).
// Compute extra-alpha bounds for all comparisons. Returns a Ncomps x
// (2 * n_extras) matrix. Columns are arranged as alternating
// delta_R_aXX, delta_W_aXX pairs for each extra alpha in extralist order.
// Column names are built separately by _rob_extras_colnames.
real matrix _rob_compute_extras(real matrix theta, real matrix DRAWS,
                                string scalar comparisons, string scalar comp_prefix,
                                string scalar extralist) {
    real scalar      nc, c, na, a, K, B, i, alpha, q_R, q_W
    real scalar      W_obs, R_obs
    real matrix      th, D, Vhat, Rmat, RVR, RVRinv, d, d_b, out
    real colvector   W_boot, R_boot, _ok
    real rowvector   cols
    string rowvector compnames, alphas

    compnames = tokens(comparisons)
    nc        = cols(compnames)
    alphas    = tokens(extralist)
    na        = cols(alphas)
    out       = J(nc, 2 * na, .)

    for (c=1; c<=nc; c++) {
        cols = _rob_str2cols(st_local(comp_prefix + compnames[c]))
        th = theta[cols', 1]
        D  = DRAWS[., cols]
        K  = cols(D)
        _ok = rowmissing(D) :== 0
        D   = select(D, _ok)
        B   = rows(D)

        Vhat   = variance(D)
        Rmat   = _rob_contrast(K)
        RVR    = Rmat * Vhat * Rmat'
        RVRinv = invsym(RVR)

        d     = Rmat * th
        W_obs = (d' * RVRinv * d)[1,1]
        R_obs = max(th) - min(th)

        W_boot = J(B, 1, .)
        R_boot = J(B, 1, .)
        for (i=1; i<=B; i++) {
            d_b       = Rmat * D[i,.]'
            W_boot[i] = (d_b' * RVRinv * d_b)[1,1]
            R_boot[i] = max(D[i,.]) - min(D[i,.])
        }

        for (a=1; a<=na; a++) {
            alpha = strtoreal(alphas[a])
            q_R   = _rob_quantile(R_boot :- R_obs, 1-alpha)
            q_W   = _rob_quantile(W_boot :- W_obs, 1-alpha)
            out[c, 2*(a-1) + 1] = R_obs + q_R
            out[c, 2*(a-1) + 2] = sqrt(max((W_obs + q_W, 0)))
        }
    }

    return(out)
}

// Build column names for r(extra) from the extra-alphas list. Returns a
// row vector of length 2 * n_extras: delta_R_aXX, delta_W_aXX alternating,
// where XX is the two-digit alpha (alpha=.10 -> "a10").
string rowvector _rob_extras_colnames(string scalar extralist) {
    real scalar      na, a, atag, alpha
    string rowvector alphas, names
    string scalar    suffix

    alphas = tokens(extralist)
    na     = cols(alphas)
    names  = J(1, 2 * na, "")
    for (a=1; a<=na; a++) {
        alpha = strtoreal(alphas[a])
        atag  = round(alpha * 100)
        if (atag < 10) suffix = "a0" + strofreal(atag)
        else           suffix = "a"  + strofreal(atag)
        names[2*(a-1) + 1] = "delta_R_" + suffix
        names[2*(a-1) + 2] = "delta_W_" + suffix
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
    printf("  Note: the robustness ratio is delta*(.05)/|theta_bar|. When |theta_bar| is close to zero,\n")
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
                      "delta_R_50", "delta_R_05", "delta_W_50", "delta_W_05",
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

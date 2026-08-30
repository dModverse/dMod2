// parvec subsetting and concatenation
//
// A `parvec` is a named numeric carrying `deriv` [p x theta], optionally
// `deriv2` [p x theta x theta], and `fixed` (the names whose deriv row is
// absent). Subsetting and concatenating it is pure bookkeeping, but the
// composition protocol does it a few hundred times per objective evaluation on
// vectors of a dozen elements, where the interpreter costs more than the work.
//
// Both entry points reproduce the R implementations in R/parClass.R exactly,
// including which cases keep the incoming array rather than copying it, and
// return R_NilValue for anything they do not handle so the R side can fall
// back. Name lookup goes through R's global string cache: two names are equal
// exactly when their CHARSXP pointers are.

#include <Rcpp.h>
#include <cstring>
#include <vector>

using namespace Rcpp;

namespace {

// Same strings in the same order. CHARSXP identity is string identity.
inline bool same_names(SEXP a, SEXP b) {
  if (a == b) return true;
  if (a == R_NilValue || b == R_NilValue) return false;
  if (TYPEOF(a) != STRSXP || TYPEOF(b) != STRSXP) return false;
  const R_xlen_t n = XLENGTH(a);
  if (n != XLENGTH(b)) return false;
  for (R_xlen_t i = 0; i < n; ++i)
    if (STRING_ELT(a, i) != STRING_ELT(b, i)) return false;
  return true;
}

// Row names of a matrix or the first dimnames entry of an array.
inline SEXP dim_names(SEXP a, int which) {
  SEXP dn = Rf_getAttrib(a, R_DimNamesSymbol);
  if (dn == R_NilValue || TYPEOF(dn) != VECSXP || XLENGTH(dn) <= which)
    return R_NilValue;
  return VECTOR_ELT(dn, which);
}

// Positions of `want` in `have`, 0-based, -1 when absent.
inline void match_pos(SEXP want, SEXP have, std::vector<int>& out) {
  const R_xlen_t n = XLENGTH(want);
  out.assign((size_t) n, -1);
  if (have == R_NilValue) return;
  SEXP m = PROTECT(Rf_match(have, want, NA_INTEGER));
  const int* mi = INTEGER(m);
  for (R_xlen_t i = 0; i < n; ++i)
    out[(size_t) i] = (mi[i] == NA_INTEGER) ? -1 : mi[i] - 1;
  UNPROTECT(1);
}

// Rows `keep` of a numeric matrix, dimnames carried along.
SEXP subset_rows_mat(SEXP m, const std::vector<int>& keep, SEXP keep_names) {
  const int nr = Rf_nrows(m), nc = Rf_ncols(m);
  const int k  = (int) keep.size();
  SEXP out = PROTECT(Rf_allocMatrix(REALSXP, k, nc));
  const double* src = REAL(m);
  double* dst = REAL(out);
  for (int j = 0; j < nc; ++j)
    for (int i = 0; i < k; ++i)
      dst[i + (R_xlen_t) j * k] = src[keep[(size_t) i] + (R_xlen_t) j * nr];
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 0, keep_names);
  SET_VECTOR_ELT(dn, 1, dim_names(m, 1));
  Rf_setAttrib(out, R_DimNamesSymbol, dn);
  UNPROTECT(2);
  return out;
}

// Rows `keep` of a [p x q x r] array.
SEXP subset_rows_arr(SEXP a, const std::vector<int>& keep, SEXP keep_names) {
  SEXP d = Rf_getAttrib(a, R_DimSymbol);
  const int nr = INTEGER(d)[0], n2 = INTEGER(d)[1], n3 = INTEGER(d)[2];
  const int k = (int) keep.size();
  SEXP out = PROTECT(Rf_alloc3DArray(REALSXP, k, n2, n3));
  const double* src = REAL(a);
  double* dst = REAL(out);
  for (int c = 0; c < n3; ++c)
    for (int j = 0; j < n2; ++j) {
      const R_xlen_t so = (R_xlen_t) nr * (j + (R_xlen_t) n2 * c);
      const R_xlen_t doff = (R_xlen_t) k * (j + (R_xlen_t) n2 * c);
      for (int i = 0; i < k; ++i) dst[doff + i] = src[so + keep[(size_t) i]];
    }
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(dn, 0, keep_names);
  SET_VECTOR_ELT(dn, 1, dim_names(a, 1));
  SET_VECTOR_ELT(dn, 2, dim_names(a, 2));
  Rf_setAttrib(out, R_DimNamesSymbol, dn);
  UNPROTECT(2);
  return out;
}

inline bool is_3d_real(SEXP a) {
  if (a == R_NilValue || TYPEOF(a) != REALSXP) return false;
  SEXP d = Rf_getAttrib(a, R_DimSymbol);
  return d != R_NilValue && XLENGTH(d) == 3;
}

// The names in `nms` that carry no deriv row.
SEXP fixed_names(SEXP nms, SEXP rn) {
  std::vector<int> pos;
  match_pos(nms, rn, pos);
  int nmiss = 0;
  for (size_t i = 0; i < pos.size(); ++i) if (pos[i] < 0) ++nmiss;
  if (nmiss == 0) return R_NilValue;
  SEXP out = PROTECT(Rf_allocVector(STRSXP, nmiss));
  int w = 0;
  for (size_t i = 0; i < pos.size(); ++i)
    if (pos[i] < 0) SET_STRING_ELT(out, w++, STRING_ELT(nms, (R_xlen_t) i));
  UNPROTECT(1);
  return out;
}

SEXP parvec_class() {
  SEXP cl = PROTECT(Rf_allocVector(STRSXP, 2));
  SET_STRING_ELT(cl, 0, Rf_mkChar("parvec"));
  SET_STRING_ELT(cl, 1, Rf_mkChar("numeric"));
  UNPROTECT(1);
  return cl;
}

// deriv / deriv2 for a target name set, or R_NilValue when no row survives.
// Returns the input array untouched when its rows already are `nms` in order.
SEXP rows_for(SEXP arr, SEXP nms, bool three_d) {
  if (arr == R_NilValue) return R_NilValue;
  SEXP rn = three_d ? dim_names(arr, 0) : dim_names(arr, 0);
  if (same_names(rn, nms)) return arr;
  std::vector<int> pos;
  match_pos(nms, rn, pos);
  std::vector<int> keep;
  std::vector<R_xlen_t> src;
  keep.reserve(pos.size());
  src.reserve(pos.size());
  for (size_t i = 0; i < pos.size(); ++i)
    if (pos[i] >= 0) { keep.push_back(pos[i]); src.push_back((R_xlen_t) i); }
  if (keep.empty()) return R_NilValue;
  SEXP kn = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t) keep.size()));
  for (size_t i = 0; i < src.size(); ++i)
    SET_STRING_ELT(kn, (R_xlen_t) i, STRING_ELT(nms, src[i]));
  SEXP out = three_d ? subset_rows_arr(arr, keep, kn)
                     : subset_rows_mat(arr, keep, kn);
  UNPROTECT(1);
  return out;
}

}  // namespace

// Attach the derivative attributes of a parvec to an already-subset numeric.
// `out` carries the values and names; `deriv`/`deriv2` are the originals, whose
// rows are restricted to those names. Returns R_NilValue for shapes it does not
// handle so `[.parvec` can fall back. Index semantics stay in R: this never
// looks at how `out` was selected.
// [[Rcpp::export]]
SEXP parvec_attach(SEXP out, SEXP deriv, SEXP deriv2) {
  if (TYPEOF(out) != REALSXP) return R_NilValue;
  SEXP nms = Rf_getAttrib(out, R_NamesSymbol);
  if (nms == R_NilValue && XLENGTH(out) > 0) return R_NilValue;
  if (deriv != R_NilValue && (TYPEOF(deriv) != REALSXP || !Rf_isMatrix(deriv)))
    return R_NilValue;
  if (deriv2 != R_NilValue && !is_3d_real(deriv2)) return R_NilValue;

  const R_xlen_t n = XLENGTH(out);
  // Always a fresh object: `out` is bound in the caller, so writing attributes
  // into it would be visible there. (Reference-count macros are not API.)
  SEXP res = PROTECT(Rf_shallow_duplicate(out));
  SEXP d1 = PROTECT(rows_for(deriv, nms, false));
  SEXP d2 = PROTECT(rows_for(deriv2, nms, true));
  Rf_setAttrib(res, Rf_install("deriv"), d1);
  Rf_setAttrib(res, Rf_install("deriv2"), d2);
  if (d1 != R_NilValue && Rf_nrows(d1) < (int) n) {
    SEXP fx = PROTECT(fixed_names(nms, dim_names(d1, 0)));
    Rf_setAttrib(res, Rf_install("fixed"), fx);
    UNPROTECT(1);
  } else {
    Rf_setAttrib(res, Rf_install("fixed"), R_NilValue);
  }
  SEXP cl = PROTECT(parvec_class());
  Rf_setAttrib(res, R_ClassSymbol, cl);
  UNPROTECT(4);
  return res;
}

// Concatenation of parvecs. Returns R_NilValue for anything it does not handle.
// [[Rcpp::export]]
SEXP parvec_concat(SEXP lst) {
  if (TYPEOF(lst) != VECSXP) return R_NilValue;
  const R_xlen_t m = XLENGTH(lst);
  if (m == 0) return R_NilValue;

  R_xlen_t n = 0;
  bool any_deriv = false, any_deriv2 = false;
  int theta = -1;
  SEXP theta_names = R_NilValue, theta2_names = R_NilValue,
       theta3_names = R_NilValue;
  for (R_xlen_t k = 0; k < m; ++k) {
    SEXP e = VECTOR_ELT(lst, k);
    if (TYPEOF(e) != REALSXP) return R_NilValue;
    if (Rf_getAttrib(e, R_NamesSymbol) == R_NilValue && XLENGTH(e) > 0)
      return R_NilValue;
    n += XLENGTH(e);
    SEXP d = Rf_getAttrib(e, Rf_install("deriv"));
    if (d != R_NilValue) {
      if (TYPEOF(d) != REALSXP || !Rf_isMatrix(d)) return R_NilValue;
      // Every block must already live in one theta basis; the R path relies on
      // rbind() for that and would bind mismatched columns silently.
      if (!any_deriv) theta_names = dim_names(d, 1);
      if (theta < 0) theta = Rf_ncols(d);
      else if (Rf_ncols(d) != theta) return R_NilValue;
      any_deriv = true;
    }
    SEXP d2 = Rf_getAttrib(e, Rf_install("deriv2"));
    if (d2 != R_NilValue) {
      if (!is_3d_real(d2)) return R_NilValue;
      SEXP dd = Rf_getAttrib(d2, R_DimSymbol);
      if (INTEGER(dd)[1] != INTEGER(dd)[2]) return R_NilValue;
      // The R path takes the deriv2 row labels from each block's own names, so
      // only handle blocks where those already agree.
      if (!same_names(dim_names(d2, 0), Rf_getAttrib(e, R_NamesSymbol)))
        return R_NilValue;
      if (!any_deriv2) {
        theta2_names = dim_names(d2, 1); theta3_names = dim_names(d2, 2);
      }
      if (theta < 0) theta = INTEGER(dd)[1];
      else if (INTEGER(dd)[1] != theta) return R_NilValue;
      any_deriv2 = true;
    }
  }

  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, n));
  double* ov = REAL(out);
  R_xlen_t w = 0;
  for (R_xlen_t k = 0; k < m; ++k) {
    SEXP e = VECTOR_ELT(lst, k);
    SEXP en = Rf_getAttrib(e, R_NamesSymbol);
    const R_xlen_t ne = XLENGTH(e);
    const double* ev = REAL(e);
    for (R_xlen_t i = 0; i < ne; ++i, ++w) {
      ov[w] = ev[i];
      SET_STRING_ELT(nms, w, STRING_ELT(en, i));
    }
  }
  {
    SEXP chk = PROTECT(Rf_match(nms, nms, NA_INTEGER));
    const int* ci = INTEGER(chk);
    for (R_xlen_t i = 0; i < n; ++i)
      if (ci[i] != (int) i + 1) { UNPROTECT(3); return R_NilValue; }
    UNPROTECT(1);
  }
  Rf_setAttrib(out, R_NamesSymbol, nms);

  int nprot = 2;
  if (any_deriv) {
    // Rows of the blocks that carry a deriv, stacked in order -- rbind().
    R_xlen_t nr = 0;
    for (R_xlen_t k = 0; k < m; ++k) {
      SEXP d = Rf_getAttrib(VECTOR_ELT(lst, k), Rf_install("deriv"));
      if (d != R_NilValue) nr += Rf_nrows(d);
    }
    SEXP J = PROTECT(Rf_allocMatrix(REALSXP, (int) nr, theta)); ++nprot;
    SEXP rn = PROTECT(Rf_allocVector(STRSXP, nr)); ++nprot;
    double* jv = REAL(J);
    R_xlen_t r0 = 0;
    for (R_xlen_t k = 0; k < m; ++k) {
      SEXP d = Rf_getAttrib(VECTOR_ELT(lst, k), Rf_install("deriv"));
      if (d == R_NilValue) continue;
      const int dr = Rf_nrows(d);
      const double* dv = REAL(d);
      SEXP drn = dim_names(d, 0);
      for (int j = 0; j < theta; ++j)
        for (int i = 0; i < dr; ++i)
          jv[(r0 + i) + (R_xlen_t) j * nr] = dv[i + (R_xlen_t) j * dr];
      for (int i = 0; i < dr; ++i)
        SET_STRING_ELT(rn, r0 + i,
                       drn == R_NilValue ? R_BlankString : STRING_ELT(drn, i));
      r0 += dr;
    }
    SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2)); ++nprot;
    SET_VECTOR_ELT(dn, 0, rn);
    SET_VECTOR_ELT(dn, 1, theta_names);
    Rf_setAttrib(J, R_DimNamesSymbol, dn);
    Rf_setAttrib(out, Rf_install("deriv"), J);

    if (any_deriv2) {
      // Blocks without a deriv2 contribute a zero slab, as the R path does.
      SEXP H = PROTECT(Rf_alloc3DArray(REALSXP, (int) n, theta, theta)); ++nprot;
      double* hv = REAL(H);
      memset(hv, 0, sizeof(double) * (size_t) n * theta * theta);
      R_xlen_t o0 = 0;
      for (R_xlen_t k = 0; k < m; ++k) {
        SEXP e = VECTOR_ELT(lst, k);
        const R_xlen_t ne = XLENGTH(e);
        SEXP d2 = Rf_getAttrib(e, Rf_install("deriv2"));
        if (d2 != R_NilValue) {
          const int p = INTEGER(Rf_getAttrib(d2, R_DimSymbol))[0];
          const double* sv = REAL(d2);
          for (int c = 0; c < theta; ++c)
            for (int j = 0; j < theta; ++j)
              for (int i = 0; i < p; ++i)
                hv[(o0 + i) + n * (j + (R_xlen_t) theta * c)] =
                    sv[i + (R_xlen_t) p * (j + (R_xlen_t) theta * c)];
        }
        o0 += ne;
      }
      SEXP dn2 = PROTECT(Rf_allocVector(VECSXP, 3)); ++nprot;
      SET_VECTOR_ELT(dn2, 0, nms);
      SET_VECTOR_ELT(dn2, 1, theta2_names);
      SET_VECTOR_ELT(dn2, 2, theta3_names);
      Rf_setAttrib(H, R_DimNamesSymbol, dn2);
      Rf_setAttrib(out, Rf_install("deriv2"), H);
    }
    if (nr < n) {
      SEXP fx = PROTECT(fixed_names(nms, rn)); ++nprot;
      Rf_setAttrib(out, Rf_install("fixed"), fx);
    }
  }

  SEXP cl = PROTECT(parvec_class()); ++nprot;
  Rf_setAttrib(out, R_ClassSymbol, cl);
  UNPROTECT(nprot);
  return out;
}

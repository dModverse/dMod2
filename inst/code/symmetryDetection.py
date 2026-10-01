# Structural identifiability of ODE models. Self-contained module driven from
# R/symmetryDetection.R via reticulate. Exact engines below carry these internal
# tags; the R-side `method` maps onto them ("polynomial" -> "liesym"):
#
#   "observability": local identifiability from the rank of the
#     observability-identifiability matrix O = d/dz [g, L_f g, ...]. For a
#     rational model O is built by a Taylor-mode construction over GF(p) at a
#     rational point (no symbolic O), the rank certified across several primes;
#     this scales to large, deep systems. Non-rational observables route to a
#     symbolic build whose rank is certified at several generic points, with an
#     optional closed-form nullspace reconstructed only in the relevant
#     variables.
#   "liesym": the polynomial Lie-symmetry ansatz of Merkt et al. 2015
#     (PRE 92, 012920); the determining system X(g)=0, [f,X]=0 (and the Lie
#     chain X(L_f^k g)=0 and the steady-state field [f_ss,X]=0) is solved in
#     exact rational arithmetic (modular GF(p) + CRT + rational reconstruction,
#     validated, sympy fallback; exact=False selects a floating-point path) and
#     every generator is verified symbolically.
#   "scaling": scaling (toric) symmetries from the integer kernel of the
#     monomial-exponent conditions; pure integer linear algebra.
#
# Observation functions may be non-rational with rational gradient (log10);
# gradients are taken before any polynomial is formed. A clean parse_expr
# symbol table keeps state and parameter names unmangled, and an optional
# symengine backend accelerates the symbolic build.

import io
import os
import sys
import math
import tokenize
import random
from fractions import Fraction

import numpy as np
import sympy as spy
from sympy.parsing.sympy_parser import parse_expr
from sympy.polys import galoistools as _gft
from sympy.polys.domains import ZZ as _ZZ

try:
    import symengine as _seng
    _HAVE_SYMENGINE = True
except Exception:
    _HAVE_SYMENGINE = False

# module state set per run
_EXACT = True
_BACKEND = "sympy"
_warned_symengine = False

# primes below 1518500213, msolve's bound for FGLM
_PRIMES = [1518500183, 1518500173, 1518500171, 1518500143,
           1518500141, 1518500131, 1518500101, 1518500077]


# ---- utilities: parsing, symbol tables, modular nullspace & rational reconstruction --

def _select_backend(name):
    global _BACKEND, _warned_symengine
    if name == "symengine" and not _HAVE_SYMENGINE:
        if not _warned_symengine:
            sys.stdout.write("symengine not available; using sympy.\n")
            sys.stdout.flush()
            _warned_symengine = True
        _BACKEND = "sympy"
    else:
        _BACKEND = name


def _diff(expr, var):
    """Differentiate with the selected backend, return a sympy expression."""
    if _BACKEND == "symengine" and _HAVE_SYMENGINE:
        try:
            d = _seng.diff(_seng.sympify(expr), _seng.sympify(var))
            return d._sympy_()
        except Exception:
            return spy.diff(expr, var)
    return spy.diff(expr, var)


###########################################################################
#####################     input parsing (clean)     ######################
###########################################################################

def _as_list(x):
    if x is None:
        return []
    if isinstance(x, str):
        return [x]
    return list(x)


def _clean(line):
    return line.replace('"', '').replace(',', '').replace('^', '**').strip()


def _build_symbol_table(all_lines):
    """Map every non-function identifier to a sympy Symbol, so names like
    E, I, S, O, Q, N are taken as model variables, not sympy constants."""
    syms, funcs = set(), set()
    for raw in all_lines:
        line = _clean(raw)
        if not line:
            continue
        try:
            toks = list(tokenize.generate_tokens(io.StringIO(line).readline))
        except (tokenize.TokenError, IndentationError):
            toks = []
        for i, t in enumerate(toks):
            if t.type == tokenize.NAME:
                nxt = toks[i + 1] if i + 1 < len(toks) else None
                if nxt is not None and nxt.type == tokenize.OP and nxt.string == '(':
                    funcs.add(t.string)
                else:
                    syms.add(t.string)
    syms -= funcs
    return {s: spy.Symbol(s) for s in syms}


def _function_aliases():
    """Functions used in dMod that sympy does not provide by name; mapped to
    differentiable sympy expressions so their gradients are rational."""
    x = spy.Symbol('_x_')
    return {
        'log10': spy.Lambda(x, spy.log(x) / spy.log(10)),
        'log2': spy.Lambda(x, spy.log(x) / spy.log(2)),
        'exp10': spy.Lambda(x, spy.Integer(10) ** x),
    }


def _make_parse(local_dict):
    def parse(rhs):
        return parse_expr(_clean(rhs), local_dict=local_dict, evaluate=True)
    return parse


def _make_local_parse(all_lines):
    """Symbol table (model names as symbols, not sympy constants) with the dMod
    function aliases, plus a parser bound to it."""
    local = _build_symbol_table(all_lines)
    local.update(_function_aliases())
    return local, _make_parse(local)


def _read_equations(lines, parse):
    variables, functions = [], []
    for raw in lines:
        line = _clean(raw)
        if '=' not in line:
            continue
        lhs, rhs = line.split('=', 1)
        variables.append(spy.Symbol(lhs.strip()))
        functions.append(parse(rhs))
    allsyms = set()
    for f in functions:
        allsyms |= set(spy.sympify(f).free_symbols)
    parameters = sorted(allsyms - set(variables), key=spy.default_sort_key)
    return variables, functions, parameters


###########################################################################
#####################     exact polynomial class     #####################
###########################################################################

class Apoly:
    """Sparse multivariate polynomial. Coefficients are kept exact (sympy
    Rational) when _EXACT, else float64 (the floating-point path)."""

    def __init__(self, expr, variables, rs):
        self.vars = variables
        self.rs = rs
        if expr is None:
            self.coefs = []
            self.exps = []
            return
        poly = spy.Poly(spy.sympify(expr), variables).as_dict()
        if rs is None:
            self.coefs = list(poly.values())
        else:
            coefsTmp = list(poly.values())
            self.coefs = [None] * len(coefsTmp)
            for i in range(len(coefsTmp)):
                if _EXACT:
                    row = np.empty(len(rs), dtype=object)
                    for j, r in enumerate(rs):
                        row[j] = spy.diff(coefsTmp[i], r) if coefsTmp[i].has(r) \
                            else spy.Integer(0)
                else:
                    row = np.zeros(len(rs))
                    for j, r in enumerate(rs):
                        if coefsTmp[i].has(r):
                            row[j] = float(spy.diff(coefsTmp[i], r))
                self.coefs[i] = row
        self.exps = [np.array(e) for e in poly.keys()]

    def __repr__(self):
        return str(self.coefs) + '\n' + str(self.exps)

    def getCopy(self):
        newPoly = Apoly(None, self.vars, self.rs)
        newPoly.coefs = [c.copy() if hasattr(c, "copy") else c for c in self.coefs]
        newPoly.exps = [e.copy() for e in self.exps]
        return newPoly

    def add(self, otherPoly):
        for i in range(len(otherPoly.exps)):
            for j in range(len(self.exps)):
                if np.array_equal(otherPoly.exps[i], self.exps[j]):
                    self.coefs[j] = self.coefs[j] + otherPoly.coefs[i]
                    if not np.any(self.coefs[j]):
                        self.coefs.pop(j)
                        self.exps.pop(j)
                    break
            else:
                self.coefs.append(otherPoly.coefs[i])
                self.exps.append(otherPoly.exps[i])

    def sub(self, otherPoly):
        for i in range(len(otherPoly.exps)):
            for j in range(len(self.exps)):
                if np.array_equal(otherPoly.exps[i], self.exps[j]):
                    self.coefs[j] = self.coefs[j] - otherPoly.coefs[i]
                    if not np.any(self.coefs[j]):
                        self.coefs.pop(j)
                        self.exps.pop(j)
                    break
            else:
                self.coefs.append(-1 * otherPoly.coefs[i])
                self.exps.append(otherPoly.exps[i])

    def mul(self, otherPoly):
        newPoly = Apoly(None, self.vars, self.rs)
        n = len(self.coefs) * len(otherPoly.coefs)
        newPoly.coefs = [0] * n
        newPoly.exps = [0] * n
        k = 0
        for i in range(len(otherPoly.exps)):
            for j in range(len(self.exps)):
                # works only because at most one factor carries rs
                newPoly.coefs[k] = otherPoly.coefs[i] * self.coefs[j]
                newPoly.exps[k] = otherPoly.exps[i] + self.exps[j]
                k += 1
        i = 0
        while i < len(newPoly.coefs):
            j = i + 1
            while j < len(newPoly.coefs):
                if np.array_equal(newPoly.exps[i], newPoly.exps[j]):
                    newPoly.exps.pop(j)
                    newPoly.coefs[i] = newPoly.coefs[i] + newPoly.coefs.pop(j)
                else:
                    j += 1
            i += 1
        return newPoly

    def diff(self, j):
        newPoly = self.getCopy()
        i = 0
        while i < len(newPoly.exps):
            if newPoly.exps[i][j] != 0:
                newPoly.coefs[i] = newPoly.coefs[i] * int(newPoly.exps[i][j])
                newPoly.exps[i][j] -= 1
                i += 1
            else:
                newPoly.coefs.pop(i)
                newPoly.exps.pop(i)
        return newPoly

    def as_expr(self):
        expr = 0
        for i in range(len(self.coefs)):
            fact = 1
            for j in range(len(self.vars)):
                fact = fact * self.vars[j] ** int(self.exps[i][j])
            if self.rs is None:
                expr += self.coefs[i] * fact
            else:
                coef = 0
                for j in range(len(self.rs)):
                    coef += self.rs[j] * self.coefs[i][j]
                expr += coef * fact
        return spy.nsimplify(expr)


###########################################################################
#####################     infinitesimal ansatz     #######################
###########################################################################

# ---- polynomial (liesym) Lie-symmetry engine: ansatz, determining system, verify -----

def _sym(name):
    return spy.Symbol(name)


def giveDegree(vars, i, p, summand, poly, num, k, rs):
    if i == len(vars) - 1:
        rs.append(_sym('r_' + str(vars[k]) + '_' + str(num)))
        poly += rs[-1] * summand * vars[i] ** p
        return poly, num + 1
    else:
        for j in range(p + 1):
            poly, num = giveDegree(vars, i + 1, p - j, summand * vars[i] ** j,
                                   poly, num, k, rs)
    return poly, num


def makeAnsatz(ansatz, allVariables, m, q, pMax, fixed):
    n = len(allVariables)
    rs = []
    infis = []

    if ansatz == 'uni':
        for k in range(n):
            infis.append(spy.sympify(0))
            if allVariables[k] in fixed:
                continue
            for p in range(pMax + 1):
                rs.append(_sym('r_' + str(allVariables[k]) + '_' + str(p)))
                infis[-1] += rs[-1] * allVariables[k] ** p
        diffInfis = [[0] * n]
        for i in range(n):
            diffInfis[0][i] = spy.diff(infis[i], allVariables[i])

    elif ansatz == 'par':
        for k in range(n):
            infis.append(spy.sympify(0))
            if allVariables[k] in fixed:
                continue
            num = 0
            for p in range(pMax + 1):
                vari = list(allVariables[m + q:])
                if k < (m + q):
                    vari.append(allVariables[k])
                    kp = len(vari) - 1
                else:
                    kp = k - (m + q)
                degree, num = giveDegree(vari, 0, p, 1, 0, num, kp, rs)
                infis[-1] += degree
        diffInfis = [[0] * n]
        for i in range(n):
            diffInfis[0][i] = spy.diff(infis[i], allVariables[i])

    elif ansatz == 'multi':
        for k in range(n):
            infis.append(spy.sympify(0))
            if allVariables[k] in fixed:
                continue
            num = 0
            for p in range(pMax + 1):
                if k < m:
                    vari = list(allVariables[:m]) + list(allVariables[m + q:])
                    kp = k
                elif k < m + q:
                    vari = list(allVariables[:])
                    kp = k
                else:
                    vari = list(allVariables[m + q:])
                    kp = k - (m + q)
                degree, num = giveDegree(vari, 0, p, 1, 0, num, kp, rs)
                infis[-1] += degree
        diffInfis = [0] * n
        for i in range(n):
            diffInfis[i] = [0] * n
        for i in range(n):
            for j in range(n):
                diffInfis[i][j] = spy.diff(infis[i], allVariables[j])
    else:
        raise UserWarning("ansatz must be one of 'uni', 'par', 'multi'")

    return infis, diffInfis, rs


def transformInfisToPoly(infis, diffInfis, allVariables, rs):
    n = len(allVariables)
    k = len(diffInfis)
    infisPoly = [Apoly(infis[i], allVariables, rs) for i in range(n)]
    diffInfisPoly = [[0] * n for _ in range(k)]
    for a in range(k):
        for i in range(n):
            diffInfisPoly[a][i] = Apoly(diffInfis[a][i], allVariables, rs)
    return infisPoly, diffInfisPoly


###########################################################################
#####################     determining equations     ######################
###########################################################################

def _quotient_derivatives(numerators, denominators, allVariables):
    m = len(numerators)
    n = len(allVariables)
    derivativesNum = [[None] * n for _ in range(m)]
    for k in range(m):
        for l in range(n):
            d = Apoly(None, allVariables, None)
            d.add(numerators[k].diff(l).mul(denominators[k]))
            d.sub(numerators[k].mul(denominators[k].diff(l)))
            derivativesNum[k][l] = d
    return derivativesNum


def doEquation(k, numerators, denominators, derivativesNum, infis, diffInfis,
               allVariables, rs, ansatz):
    n = len(allVariables)
    m = len(numerators)
    polynomial = Apoly(None, allVariables, rs)
    if ansatz in ('uni', 'par'):
        polynomial.add(diffInfis[0][k].mul(denominators[k]).mul(numerators[k]))
        for i in range(n):
            polynomial.sub(infis[i].mul(derivativesNum[k][i]))
    else:  # multi
        for j in range(m):
            summand = diffInfis[k][j].mul(denominators[k]).mul(numerators[j])
            for l in range(m):
                if l != j:
                    summand = summand.mul(denominators[l])
            polynomial.add(summand)
        for i in range(n):
            summand = infis[i].mul(derivativesNum[k][i])
            for l in range(m):
                if l != k:
                    summand = summand.mul(denominators[l])
            polynomial.sub(summand)
    return list(polynomial.coefs)


def _exp_det_rows(fields, obsExprs, infisSym, diffInfisSym, allVariables, rs,
                  ansatz, m):
    """Determining rows of a model with exponentials: each canonical exponential
    W = r^(tau/L) (_exp_leaf_canon) is an extra polynomial variable, algebraically
    independent of the coordinates, and enters derivatives by the chain rule. The
    ansatz stays in the coordinates."""
    taken = {str(v) for v in allVariables}
    for e in list(fields) + list(obsExprs):
        taken |= {str(x) for x in spy.sympify(e).free_symbols}
    ctx = _ExpCtx(taken)
    lc = _exp_leaf_canon([spy.sympify(e) for e in list(fields) + list(obsExprs)], ctx)
    if 'why' in lc:
        raise UserWarning(lc['why'])
    F, G = lc['exprs'][:len(fields)], lc['exprs'][len(fields):]
    atoms = lc['atoms']
    extra = set()
    for e in lc['exprs'] + [w['tau'] for w in atoms]:
        extra |= e.free_symbols
    extVars = list(allVariables) + sorted(extra - set(allVariables), key=str)

    def D(e, z):
        out = spy.diff(e, z)
        for w in atoms:
            dW = spy.diff(e, w['W'])
            if dW != 0:
                out += dW * w['W'] * ctx.lnOf(w['r']) * D(w['tau'] / w['L'], z)
        return out

    def coefRows(expr):
        num = spy.expand(spy.fraction(spy.together(expr))[0])
        return list(Apoly(num, extVars, rs).coefs) if num != 0 else []

    n = len(allVariables)
    rows = []
    for k in range(m):
        if ansatz in ('uni', 'par'):
            lhs = diffInfisSym[0][k] * F[k]
        else:
            lhs = sum(diffInfisSym[k][j] * F[j] for j in range(m))
        rhs = sum(infisSym[i] * D(F[k], allVariables[i])
                  for i in range(n) if infisSym[i] != 0)
        rows.extend(coefRows(lhs - rhs))
    for g in G:
        rows.extend(coefRows(sum(infisSym[l] * D(g, allVariables[l])
                                 for l in range(n) if infisSym[l] != 0)))
    back = {sm: spy.log(pr) for pr, sm in ctx.lnSym.items()}
    back[ctx.ec] = spy.exp(spy.Rational(1, getattr(ctx, 'L0', 1)))
    return rows, back


def _obs_rows(obsExpr, infisSym, allVariables, rs):
    # X(g) = sum_l xi_l dg/dvar_l = 0. Only the gradient of g enters, so g may
    # be non-rational (e.g. log10) as long as every dg/dvar is rational; the
    # numerator of together(X(g)) is then polynomial in the variables.
    n = len(allVariables)
    Xg = spy.Integer(0)
    for l in range(n):
        d = _diff(obsExpr, allVariables[l])
        if d != 0:
            Xg += infisSym[l] * d
    num, _ = spy.fraction(spy.together(Xg))
    num = spy.expand(num)
    try:
        poly = Apoly(num, allVariables, rs)
    except spy.PolynomialError:
        raise UserWarning(
            "Observable gradient is not rational (e.g. exp/sqrt mixed with "
            "other terms). The liesym engine needs rational gradients; strip "
            "an invertible outer function g = phi(h) to its argument h, or use "
            "method='observability'.")
    return list(poly.coefs)


###########################################################################
#####################     exact / float solvers     ######################
###########################################################################

def _rref_mod_p(A, p):
    """Vectorized int64 Gauss-Jordan over GF(p); returns (R, pivots). Eliminates
    only rows nonzero in the pivot column (the scaling matrix is sparse) by a rank-1
    update; products < p^2 < 2^62 fit int64."""
    A = (np.asarray(A, dtype=np.int64) % p)
    nrows, ncols = A.shape
    pivots = []
    r = 0
    for c in range(ncols):
        nz = np.nonzero(A[r:, c])[0]
        if nz.size == 0:
            continue
        piv = r + int(nz[0])
        if piv != r:
            A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        # the pivot row is zero left of c, so the updates start at c
        A[r, c:] = (A[r, c:] * inv) % p
        col = A[:, c].copy()
        col[r] = 0
        nzr = np.nonzero(col)[0]
        if nzr.size:
            A[nzr, c:] = (A[nzr, c:] - np.outer(col[nzr], A[r, c:])) % p
        pivots.append(c)
        r += 1
        if r == nrows:
            break
    return A, pivots


def _rational_reconstruct(a, m):
    """Recover n/d with a == n*d^{-1} (mod m), |n|,|d| bounded by sqrt(m/2)."""
    a %= m
    if a == 0:
        return spy.Integer(0)
    bound = math.isqrt(m // 2)
    r0, r1 = m, a
    s0, s1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound:
        return None
    return spy.Rational(r1, s1)


def symRatReconBig(residues, primes):
    """Arbitrary-precision per-row CRT + rational reconstruction, free of the u128
    cap (4 primes) of the C++ symRatRecon. Each row of `residues` is one coefficient
    across `primes`; returns {'num', 'den'} as decimal strings, den '0' where no lift
    exists. For log-carrying coefficients whose height exceeds the 4-prime bound."""
    from sympy.ntheory.modular import crt
    mods = [int(p) for p in primes]
    M = 1
    for m in mods:
        M *= m
    rows = [[int(x) for x in row] for row in np.asarray(residues, dtype=object)]
    num, den = [], []
    for row in rows:
        res = [r % m for r, m in zip(row, mods)]
        x, _ = crt(mods, res)
        val = _rational_reconstruct(int(x), int(M))
        if val is None:
            num.append("0"); den.append("0")
        else:
            val = spy.Rational(val)
            num.append(str(int(val.p))); den.append(str(int(val.q)))
    return {'num': num, 'den': den}


def _crt_nullspace(reduceModP, ncols):
    """Multi-prime GF(p) + CRT + rational-reconstruction nullspace, shared by the
    two exact nullspace routines below. `reduceModP(p)` yields the matrix mod p.
    Primes whose pivot set differs from the first (unlucky rank drop) are skipped;
    the first four agreeing ones are lifted to Q. Returns the unvalidated basis as
    sympy column vectors (callers validate), or None if no residue lifts."""
    ref_pivots = None
    free = None
    residues = {}
    mods = []
    # the first four primes reduce in parallel threads (the eliminations are numpy
    # array operations, which release the GIL); further primes only on a pivot clash
    first = list(_PRIMES[:4])
    try:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=len(first)) as ex:
            pre = dict(zip(first, ex.map(lambda q: _rref_mod_p(reduceModP(q), q), first)))
    except Exception:
        pre = {}
    for p in _PRIMES:
        R, pivots = pre[p] if p in pre else _rref_mod_p(reduceModP(p), p)
        if ref_pivots is None:
            ref_pivots = pivots
            pivset = set(pivots)
            free = [c for c in range(ncols) if c not in pivset]
        elif pivots != ref_pivots:
            continue
        mods.append(p)
        for ki in range(len(ref_pivots)):
            for f in free:
                residues.setdefault((ki, f), []).append(int(R[ki, f]) % p)
        if len(mods) >= 4:
            break
    if not mods:
        return None

    from sympy.ntheory.modular import crt
    exact = {}
    for key, res in residues.items():
        x, Mmod = crt(mods, res)
        val = _rational_reconstruct(int(x), int(Mmod))
        if val is None:
            return None
        exact[key] = val

    basis = []
    for f in free:
        v = [spy.Integer(0)] * ncols
        v[f] = spy.Integer(1)
        for ki, c in enumerate(ref_pivots):
            v[c] = -exact[(ki, f)]
        basis.append(spy.Matrix(v))
    return basis


def _modular_nullspace(M):
    """Exact nullspace of sympy Matrix M via _crt_nullspace, validated symbolically.
    Returns list of sympy column vectors, or None if reconstruction/validation
    fails (caller falls back to sympy)."""
    nrows, ncols = M.shape
    if nrows == 0:
        return [spy.Matrix([1 if j == c else 0 for j in range(ncols)])
                for c in range(ncols)]
    Aint = []
    for i in range(nrows):
        row = [spy.Rational(M[i, j]) for j in range(ncols)]
        L = 1
        for r in row:
            L = spy.ilcm(L, r.q)
        Aint.append([int(r * L) for r in row])

    basis = _crt_nullspace(lambda p: [[x % p for x in row] for row in Aint], ncols)
    if basis is None:
        return None
    for v in basis:
        if not (M * v).is_zero_matrix:
            return None
    return basis


def _modular_nullspace_int(rowsInt, ncols):
    """Integer-matrix analogue of _modular_nullspace (scaling determining matrix),
    in numpy int with no per-entry sympy. Returns sympy column vectors, or None if
    reconstruction or validation fails."""
    A0 = np.asarray(rowsInt, dtype=np.int64)
    basis = _crt_nullspace(lambda p: A0 % p, ncols)
    if basis is None:
        return None

    # exact validation over the integers: clear denominators per basis vector, then
    # A0 @ W == 0 in one object-dtype matmul (exact big-int, no overflow)
    if basis:
        cols = []
        for v in basis:
            L = 1
            for e in v:
                L = spy.ilcm(L, spy.Rational(e).q)
            cols.append([int(spy.Rational(e) * L) for e in v])
        W = np.array(cols, dtype=object).T
        if np.any(A0.astype(object).dot(W) != 0):
            return None
    return basis


def exactIntKernel(rows, ncols):
    """Primitive integer basis of the kernel of an integer matrix, as plain int
    lists (R-friendly): exactNullspace + the _scaling_gens content clearing."""
    basis = exactNullspace([[int(x) for x in r] for r in rows], int(ncols),
                           integer=True)
    out = []
    for v in basis:
        L = 1
        for e in v:
            L = spy.ilcm(L, spy.Rational(e).q)
        w = [int(spy.Rational(e) * L) for e in v]
        G = 0
        for e in w:
            G = spy.igcd(G, abs(e))
        if G > 1:
            w = [e // G for e in w]
        out.append(w)
    return out


def exactNullspace(rows, ncols, integer=False):
    """Nullspace basis as sympy column vectors. integer=True (integer rows) tries
    the numpy-int GF(p) + CRT path first, falling back to sympy."""
    if not rows:
        return [spy.Matrix([1 if j == c else 0 for j in range(ncols)])
                for c in range(ncols)]
    if integer:
        basis = _modular_nullspace_int(rows, ncols)
        if basis is not None:
            return basis
    M = spy.Matrix([[spy.sympify(x) for x in row] for row in rows])
    # the modular solver needs a rational matrix; transcendental constants
    # (e.g. log(10) from a log10 observable) route to the exact sympy nullspace
    if all(bool(e.is_rational) for e in M):
        basis = _modular_nullspace(M)
        if basis is not None:
            return basis
    return M.nullspace()


# ---- floating-point path (exact=False) ----

def legacyNullspace(rows, ncols):
    # Float SVD-based nullspace, valid for rectangular and underdetermined systems.
    from scipy.linalg import null_space
    if not rows:
        base = np.eye(ncols)
    else:
        A = np.array([[float(x) for x in row] for row in rows])
        base = null_space(A)
    return [spy.Matrix([spy.nsimplify(base[i, l], rational=True)
                        for i in range(ncols)])
            for l in range(base.shape[1])]


###########################################################################
#####################     post-processing / report     ###################
###########################################################################

def checkForCommonFactor(infisTmp, allVariables, m):
    factors = []
    for i in range(len(allVariables)):
        if infisTmp[i] != 0:
            fac = spy.factor(infisTmp[i])
            if fac.is_Add:
                factors = [infisTmp[i]]
            elif fac.is_Mul:
                factors = list(fac.args)
            else:
                factors = [fac]
            break
    i = 0
    while i < len(factors):
        if factors[i].is_number:
            factors.pop(i)
        elif factors[i] in allVariables[:m]:
            factors.pop(i)
        elif factors[i].is_Add:
            factors.pop(i)
        elif factors[i].is_Pow:
            if not factors[i].args[0].is_Add:
                factors[i] = factors[i].args[0]
                i += 1
            else:
                factors.pop(i)
        else:
            i += 1
    for i in range(1, len(infisTmp)):
        if infisTmp[i] == 0:
            continue
        fac = spy.factor(infisTmp[i])
        if fac.is_Mul:
            factorsTmp = list(fac.args)
        else:
            factorsTmp = [fac]
        j = 0
        while j < len(factors):
            k = 0
            while k < len(factorsTmp):
                if factorsTmp[k].is_number:
                    factorsTmp.pop(k)
                elif factorsTmp[k] in allVariables[:m]:
                    factorsTmp.pop(k)
                elif factorsTmp[k].is_Add:
                    factorsTmp.pop(k)
                elif factorsTmp[k].is_Pow:
                    if not factorsTmp[k].args[0].is_Add:
                        factorsTmp[k] = factorsTmp[k].args[0]
                        k += 1
                    else:
                        factorsTmp.pop(k)
                else:
                    k += 1
            if factors[j] in factorsTmp:
                j += 1
            else:
                factors.pop(j)
        if len(factors) != 0:
            continue
        else:
            break
    return len(factors) != 0


def buildTransformation(infis, allVariables):
    n = len(allVariables)
    epsilon = spy.Symbol('epsilon')
    transformations = [0] * n
    tType = [False] * 6
    for i in range(n):
        if infis[i] == 0:
            transformations[i] = allVariables[i]
        else:
            poly = spy.Poly(infis[i], allVariables).as_dict()
            monomials = list(poly.keys())
            coefs = list(poly.values())
            if len(monomials) == 1:
                p = None
                broke = False
                for j in range(n):
                    if monomials[0][j] != 0:
                        if j == i and p is None:
                            p = monomials[0][i]
                        elif p is None and monomials[0][j] == 1:
                            p = -1 - j
                        else:
                            transformations[i] = '-?-'
                            tType[0] = True
                            broke = True
                            break
                if not broke:
                    if p is None:
                        transformations[i] = allVariables[i] + epsilon * coefs[0]
                        tType[2] = True
                    elif p <= 0:
                        transformations[i] = allVariables[i] + epsilon * coefs[0] * allVariables[-p - 1]
                        tType[5] = True
                    elif p == 1:
                        transformations[i] = spy.exp(epsilon * coefs[0]) * allVariables[i]
                        tType[1] = True
                    else:
                        transformations[i] = spy.simplify(
                            allVariables[i] / (1 - (p - 1) * epsilon * allVariables[i] ** (p - 1)) ** (spy.sympify(1) / (p - 1)))
                        if p == 2:
                            tType[3] = True
                        else:
                            tType[4] = True
            else:
                transformations[i] = '-?-'
                tType[0] = True

    labels = ['unknown', 'scaling', 'translation', 'MM-like', 'p>2', 'gen. translation']
    parts = [labels[i] for i in range(6) if tType[i]]
    return transformations, 'Type: ' + ', '.join(parts)


def _canon_int_vector(prim, gens):
    """Scale a polynomial vector by a rational constant to canonical integer-primitive
    form (clear denominators, divide by content, first nonzero LC positive), so
    gauge-equivalent representatives coincide for a fixed `gens` order."""
    coeffs = []
    for e in prim:
        if e == 0:
            continue
        for c in spy.Poly(e, *gens).coeffs():
            coeffs.append(spy.Rational(c))
    if not coeffs:
        return prim
    L = spy.Integer(1)
    for c in coeffs:
        L = spy.ilcm(L, c.q)
    scaled = [spy.expand(e * L) for e in prim]
    g = 0
    for e in scaled:
        if e == 0:
            continue
        for c in spy.Poly(e, *gens).coeffs():
            g = spy.igcd(g, int(c))
    if g and g != 1:
        scaled = [spy.expand(e / g) for e in scaled]
    for e in scaled:
        if e != 0:
            if spy.Poly(e, *gens).LC() < 0:
                scaled = [spy.expand(-x) for x in scaled]
            break
    return scaled


def _factor_map(e):
    """{irreducible factor: multiplicity} of a nonzero polynomial expression and its
    rational content, without expanding a product that is already factored."""
    c, fl = spy.factor_list(e)
    out = {}
    for f, m in fl:
        # a factor fixed up to sign: keep one orientation, move the sign to c
        if spy.Poly(f).LC() < 0:
            f = -f
            if m % 2:
                c = -c
        out[f] = out.get(f, 0) + int(m)
    return c, out


def _classify_factored(infis, gens):
    """classifyDirection for wide directions: lcm of the denominators and gcd of the
    numerators from factor lists (exact over Q), the canonical representative kept
    as a product. None when a component is not rational."""
    try:
        nz = [i for i, e in enumerate(infis) if e != 0]
        fr = [spy.fraction(spy.together(infis[i])) for i in nz]
        num = [_factor_map(a) for a, _ in fr]
        den = [_factor_map(b) for _, b in fr]
        # lcm of the denominators, gcd of the numerators times lcm/den
        L = {}
        for _, dm in den:
            for f, m in dm.items():
                L[f] = max(L.get(f, 0), m)
        cleared = []
        for (cn, nm), (cd, dm) in zip(num, den):
            fm = dict(nm)
            for f, m in L.items():
                fm[f] = fm.get(f, 0) + m - dm.get(f, 0)
            cleared.append((cn / cd, {f: m for f, m in fm.items() if m}))
        common = None
        for _, fm in cleared:
            common = dict(fm) if common is None else \
                {f: min(m, fm.get(f, 0)) for f, m in common.items() if fm.get(f, 0)}
        coeffs = [c for c, _ in cleared]
        # integer-primitive: clear the rational coefficients, divide by their gcd
        dl = 1
        for c in coeffs:
            dl = spy.ilcm(dl, spy.fraction(spy.nsimplify(c))[1])
        ints = [int(c * dl) for c in coeffs]
        g = 0
        for v in ints:
            g = spy.igcd(g, v)
        g = g or 1
        if ints[0] < 0:
            g = -g
        prim = {}
        degs = []
        for i, (c, fm) in zip(nz, cleared):
            rest = {f: m - common.get(f, 0) for f, m in fm.items() if m - common.get(f, 0)}
            e = spy.Integer(int(c * dl) // g)
            d = 0
            for f, m in rest.items():
                e = e * f**m
                d += m * spy.Poly(f, *gens).total_degree()
            prim[i] = e
            degs.append(d)
        isScaling = all(degs[k] == 1 and prim[i].free_symbols == {gens[i]}
                        and spy.Poly(prim[i], gens[i]).total_degree() == 1
                        for k, i in enumerate(nz))
        comp_out = {str(gens[i]): str(prim[i]) for i in nz}
        if isScaling:
            weights = {str(gens[i]): int(spy.Poly(prim[i], gens[i]).LC()) for i in nz}
            return {'type': 'scaling', 'weights': weights,
                    'components': comp_out, 'degree': 1}
        return {'type': 'general', 'weights': None, 'components': comp_out,
                'degree': int(max(degs)) if degs else 0}
    except Exception:
        return None


def classifyDirection(vector):
    """Classify one non-identifiability direction and return its canonical
    poly-primitive differential generator.

    `vector` is {coordinate: component_str}, the components eta_i of
    X = sum_i eta_i d/dz_i (a scaling as {z_i: "w_i*z_i"}). X is defined up to a
    factor h(z); the gauge is fixed by clearing the common denominator, dividing by
    the polynomial content and normalising to integer-primitive form.

      scaling : every eta_i = w_i * z_i. The orbit is a ray, so any one coordinate
                can be fixed at any value.
      general : anything else. The orbit is curved, so fixing a coordinate is not
                free; the direction is removed by reparametrising onto invariants.

    Returns {'type', 'weights', 'components', 'degree'}: integer scaling weights
    (None unless scaling), the canonical generator {z_i: str(eta_i)} and its total
    degree (-1 when not polynomial)."""
    items = [(str(k), str(v)) for k, v in vector.items()]
    fallback = {'type': 'general', 'weights': None,
                'components': {k: v for k, v in items}, 'degree': -1}
    if not items:
        return {'type': 'general', 'weights': None, 'components': {}, 'degree': 0}
    # one moved coordinate: z*d/dz spans the same direction
    moved = [(k, v) for k, v in items if v.strip() != '0']
    if len(moved) == 1:
        z = moved[0][0]
        return {'type': 'scaling', 'weights': {z: 1}, 'components': {z: z}, 'degree': 1}

    # symbol table so names like E, I, S, N are model symbols, not sympy constants
    _, parse = _make_local_parse([k + " = " + v for k, v in items])
    comps = {}
    try:
        for k, v in items:
            comps[spy.Symbol(k)] = parse(v)
    except Exception:
        return fallback

    varset = set(comps.keys())
    for e in comps.values():
        varset |= set(spy.sympify(e).free_symbols)
    gens = sorted(varset, key=spy.default_sort_key)
    n = len(gens)
    infis = [comps.get(z, spy.Integer(0)) for z in gens]

    # a wide direction (many symbols) goes through factor lists: the multivariate gcd
    # of the expanded numerators costs minutes on an entry like
    # ksec*(Km + R1 + ... + R30)**2, its factor list is read off the product
    if n > 12:
        wide = _classify_factored(infis, gens)
        if wide is not None:
            return wide

    # canonical poly-primitive representative (gauge fix): clear the common
    # denominator, then divide by the polynomial content of the numerators
    try:
        dens = [spy.fraction(spy.together(e))[1] for e in infis]
        L = spy.Integer(1)
        for d in dens:
            L = spy.lcm(L, d)
        nums = [spy.expand(spy.cancel(e * L)) for e in infis]
        for e in nums:
            spy.Poly(e, *gens)                 # raises if not polynomial in gens
    except Exception:
        return fallback
    g = spy.Integer(0)
    for e in nums:
        g = spy.gcd(g, e)
    if g == 0:
        g = spy.Integer(1)
    prim = [spy.expand(spy.cancel(e / g)) for e in nums]
    prim = _canon_int_vector(prim, gens)

    # read the class off the canonical representative
    isScaling, maxdeg = True, 0
    for i, e in enumerate(prim):
        if e == 0:
            continue
        P = spy.Poly(e, *gens)
        maxdeg = max(maxdeg, P.total_degree())
        monoms = list(P.as_dict().keys())
        if not (len(monoms) == 1 and sum(monoms[0]) == 1 and monoms[0][i] == 1):
            isScaling = False

    comp_out = {str(gens[i]): str(prim[i]) for i in range(n) if prim[i] != 0}
    if isScaling:
        weights = {}
        for i in range(n):
            if prim[i] != 0:
                key = tuple(1 if j == i else 0 for j in range(n))
                weights[str(gens[i])] = int(spy.Poly(prim[i], *gens).as_dict()[key])
        return {'type': 'scaling', 'weights': weights,
                'components': comp_out, 'degree': 1}
    return {'type': 'general', 'weights': None,
            'components': comp_out, 'degree': int(maxdeg)}


def certifyInSpan(vector, generators):
    """Is the direction `vector` {coord: xi_i} a nonzero constant linear combination
    of the Lie-symmetry `generators`, i.e. an exact polynomial Lie point symmetry?
    Decided exactly by a span test on evaluations at several integer points.
    Returns {'certified': bool}."""
    items = [(str(k), str(v)) for k, v in vector.items()]
    G = [{str(k): str(v) for k, v in g.items()} for g in generators]
    if not G or not items:
        return {'certified': False}
    all_lines = [k + " = " + v for k, v in items] + \
                [k + " = " + v for g in G for k, v in g.items()]
    _, parse = _make_local_parse(all_lines)
    coords = sorted(set(k for k, _ in items) | set(k for g in G for k in g.keys()))
    try:
        xi = {k: parse(v) for k, v in items}
        Gp = [{k: parse(v) for k, v in g.items()} for g in G]
    except Exception:
        return {'certified': False}
    m = len(Gp)
    Arows, brows = [], []
    step, pts = 0, 0
    while pts < m + 3 and step < 400:
        step += 1
        subs = {spy.Symbol(coords[i]): spy.Integer(2 + (step * 13 + i * 29) % 89)
                for i in range(len(coords))}
        good, Ar, br = True, [], []
        for c in coords:
            try:
                bv = spy.Rational(xi.get(c, spy.Integer(0)).subs(subs))
                gv = [spy.Rational(Gp[j].get(c, spy.Integer(0)).subs(subs))
                      for j in range(m)]
            except Exception:
                good = False
                break
            Ar.append(gv)
            br.append(bv)
        if good:
            Arows.extend(Ar)
            brows.extend(br)
            pts += 1
    if not Arows:
        return {'certified': False}
    A = spy.Matrix(Arows)
    b = spy.Matrix(brows)
    if all(x == 0 for x in b):
        return {'certified': False}
    return {'certified': A.rank() == A.row_join(b).rank()}


def printTransformations(infisAll, allVariables, verified):
    n = len(allVariables)
    print('\n\n' + str(len(infisAll)) + ' transformation(s) found:')
    for l in range(len(infisAll)):
        for i in range(n):
            infisAll[l][i] = spy.nsimplify(infisAll[l][i])
        trans, ttype = buildTransformation(infisAll[l], allVariables)
        flag = ''
        if verified is not None:
            flag = '  [verified]' if verified[l] else '  [UNVERIFIED]'
        print('-' * 60)
        print('#' + str(l + 1) + ': ' + ttype + flag)
        for i in range(n):
            if infisAll[l][i] != 0:
                print('  {0:>14s} : {1:s}  ->  {2:s}'.format(
                    str(allVariables[i]), str(infisAll[l][i]), str(trans[i])))


###########################################################################
#####################     verification     ###############################
###########################################################################

def _is_zero(expr):
    if expr == 0:
        return True
    e = spy.together(spy.sympify(expr))
    num, _ = spy.fraction(e)
    return spy.simplify(spy.expand(num)) == 0


def _verify_generator(infis, allVariables, diffEquations, obsExprs):
    n = len(allVariables)
    # observation / Lie-chain invariance: X(o) = 0
    for o in obsExprs:
        s = 0
        for l in range(n):
            if infis[l] != 0:
                s += infis[l] * spy.diff(o, allVariables[l])
        if not _is_zero(s):
            return False
    # flow conditions [f, X] = 0; f has components only for dynamic states
    # (index < len(fields)), zero for parameters and inputs
    fields = diffEquations
    mf = len(fields)
    for k in range(mf):
        fk = fields[k]
        bracket = 0
        for l in range(n):
            fl = fields[l] if l < mf else 0
            if fl != 0 and infis[k] != 0:
                bracket += fl * spy.diff(infis[k], allVariables[l])
            if infis[l] != 0 and fk != 0:
                bracket -= infis[l] * spy.diff(fk, allVariables[l])
        if not _is_zero(bracket):
            return False
    return True


###########################################################################
#####################     core driver     ################################
###########################################################################

def _rationalize(functions, allVariables):
    """Split each expression into numerator/denominator Apoly pairs."""
    nums, dens = [], []
    for fexpr in functions:
        rational = spy.together(spy.sympify(fexpr))
        nums.append(Apoly(spy.numer(rational), allVariables, None))
        dens.append(Apoly(spy.denom(rational), allVariables, None))
    return nums, dens


def symmetryDetection(allVariables, diffEquations, obsFunctions,
                      ansatz='uni', pMax=2, inputs=(), fixed=(), lieOrder=0,
                      allTrafos=False, verify=True):
    n = len(allVariables)
    m = len(diffEquations)
    q = len(inputs)

    sys.stdout.write('Preparing equations...')
    sys.stdout.flush()

    infisSym, diffInfis, rs = makeAnsatz(ansatz, allVariables, m, q, pMax, list(fixed))
    diffInfisSym = diffInfis
    infis, diffInfis = transformInfisToPoly(infisSym, diffInfis, allVariables, rs)

    diffEquations = [_hyp_to_exp(f) for f in diffEquations]
    obsFunctions = [_hyp_to_exp(g) for g in obsFunctions]
    expMode = any(_has_exp(e) for e in list(diffEquations) + list(obsFunctions))
    if not expMode:
        numerators, denominators = _rationalize(diffEquations, allVariables)
        derivativesNum = _quotient_derivatives(numerators, denominators, allVariables)

    # observation invariance plus the Lie-derivative chain X(L_f^k g) = 0
    fieldExprs = [spy.together(spy.sympify(f)) for f in diffEquations]
    obsExprs = []
    for g in obsFunctions:
        expr = spy.sympify(g)
        chain = [expr]
        for _ in range(int(lieOrder)):
            cur = chain[-1]
            lie = 0
            for i in range(m):
                d = _diff(cur, allVariables[i])
                if d != 0:
                    lie += fieldExprs[i] * d
            chain.append(spy.together(lie))
        obsExprs.extend(chain)
    h = len(obsExprs)

    sys.stdout.write('done\nBuilding system...')
    sys.stdout.flush()

    rows, expBack = [], {}
    if expMode:
        rows, expBack = _exp_det_rows(fieldExprs[:m], obsExprs, infisSym, diffInfisSym,
                                      allVariables, rs, ansatz, m)
    else:
        for k in range(m):
            rows.extend(doEquation(k, numerators, denominators, derivativesNum,
                                   infis, diffInfis, allVariables, rs, ansatz))
        for k in range(h):
            rows.extend(_obs_rows(obsExprs[k], infisSym, allVariables, rs))

    ncols = len(rs)
    sys.stdout.write('done\nSolving system of size %dx%d (%s)...'
                     % (len(rows), ncols, 'exact' if _EXACT else 'float'))
    sys.stdout.flush()

    if _EXACT:
        basis = exactNullspace(rows, ncols)
    else:
        basis = legacyNullspace(rows, ncols)

    sys.stdout.write('done\nProcessing results...\n')
    sys.stdout.flush()

    infisAll = []
    for v in basis:
        infisTmp = [0] * n
        for i in range(n):
            poly = infis[i].getCopy()
            poly.rs = [v[j] for j in range(ncols)]
            infisTmp[i] = (_tidy_logs(poly.as_expr().xreplace(expBack))
                           if expBack else poly.as_expr())
        if allTrafos or not checkForCommonFactor(infisTmp, allVariables, m):
            infisAll.append(infisTmp)

    verified = None
    if verify:
        verified = [_verify_generator(infisTmp, allVariables, fieldExprs, obsExprs)
                    for infisTmp in infisAll]

    printTransformations(infisAll, allVariables, verified)

    # structured return
    result = []
    for l, infisTmp in enumerate(infisAll):
        trans, ttype = buildTransformation(infisTmp, allVariables)
        nz = [i for i in range(n) if infisTmp[i] != 0]
        result.append({
            'infinitesimals': {str(allVariables[i]): str(infisTmp[i]) for i in nz},
            'transformation': {str(allVariables[i]): str(trans[i]) for i in nz},
            'type': ttype,
            'verified': (None if verified is None else bool(verified[l])),
        })
    return result


# raised when an expression is not a rational function of the coordinates
class _NotRational(Exception):
    pass

###########################################################################
#####################   observability tape compiler   ####################
###########################################################################

# Opcodes for the flat tape consumed by the C++ kernel (src/symmetry_kernel.cpp).
# Integer powers are expanded into multiplications, so the kernel needs only
# these four operations.
_OP_CONST, _OP_ADD, _OP_MUL, _OP_INV = 0, 1, 2, 3


# ---- modular observability: tape emission and steady-state solves --------------------

def _is_rational_expr(e):
    """True if e is a rational function of its symbols (no transcendental
    function, no non-integer power)."""
    e = spy.sympify(e)
    for a in spy.preorder_traversal(e):
        if a.is_Function:
            return False
        if a.is_Pow and not a.exp.is_Integer:
            return False
    return True


def _emit_tape_shared(fexpr, gexpr, slotOf, base):
    """Emit a straight-line tape over an explicit slot map of states and leaves;
    instruction i writes slot base + i. States get their own slots, seeded from an
    initial condition (multi-condition compiler)."""
    op, a, b, cnum, cden = [], [], [], [], []
    memo = {}

    def emit(opcode, aa, bb, num=0, den=1):
        op.append(int(opcode))
        a.append(int(aa))
        b.append(int(bb))
        cnum.append(str(num))
        cden.append(str(den))
        return base + len(op) - 1

    def build(e):
        s = memo.get(e)
        if s is not None:
            return s
        if e.is_Symbol:
            nm = str(e)
            if nm not in slotOf:
                raise _NotRational()
            memo[e] = slotOf[nm]
            return slotOf[nm]
        if e.is_Integer:
            s = emit(_OP_CONST, 0, 0, int(e), 1)
        elif e.is_Rational:
            s = emit(_OP_CONST, 0, 0, int(e.p), int(e.q))
        elif e.is_Float:
            r = spy.nsimplify(e, rational=True)
            if not r.is_Rational:
                r = spy.Rational(e)
            s = emit(_OP_CONST, 0, 0, int(r.p), int(r.q))
        elif e.is_Add:
            args = list(e.args)
            s = build(args[0])
            for t in args[1:]:
                s = emit(_OP_ADD, s, build(t))
        elif e.is_Mul:
            args = list(e.args)
            s = build(args[0])
            for t in args[1:]:
                s = emit(_OP_MUL, s, build(t))
        elif e.is_Pow and e.exp.is_Integer:
            n = int(e.exp)
            baseS = build(e.base)
            if n == 0:
                s = emit(_OP_CONST, 0, 0, 1, 1)
            else:
                s = baseS
                for _ in range(abs(n) - 1):
                    s = emit(_OP_MUL, s, baseS)
                if n < 0:
                    s = emit(_OP_INV, s, 0)
        else:
            raise _NotRational()
        memo[e] = s
        return s

    repl, red = spy.cse(list(fexpr) + list(gexpr))
    for sym, sub in repl:
        memo[sym] = build(sub)
    outslots = [build(e) for e in red]
    return op, a, b, cnum, cden, outslots


_ssModularCache = {}
# compiled t0 event value/derivative term lists keyed by the event tuple and the
# parameter order
_eventCompileCache = {}

# test seam: force the per-point solve
_SS_FORCE_POINT = False

# force the numeric per-point steady-state seed; a test seam to cross-check it
# against the compiled steady-state IC tape
_FORCE_CONSTRAINT_SEED = False


def setSteadyStateForcePointSolve(on=True):
    """Force the per-point steady-state solve; returns the previous setting."""
    global _SS_FORCE_POINT
    prev = _SS_FORCE_POINT
    _SS_FORCE_POINT = bool(on)
    return prev


def setForceConstraintSeed(on=True):
    """Force the numeric per-point steady-state seed and return the previous setting;
    a test seam to cross-check it against the compiled steady-state IC tape."""
    global _FORCE_CONSTRAINT_SEED
    prev = _FORCE_CONSTRAINT_SEED
    _FORCE_CONSTRAINT_SEED = bool(on)
    return prev


def _modp_rational(val, p):
    """Map a sympy Rational (or integer-valued expression) to its residue mod p."""
    r = spy.Rational(val)
    return (int(r.p) % p) * pow(int(r.q) % p, p - 2, p) % p


def _poly_terms(expr, gens):
    """Compile a rational expression in `gens` to (numerator, denominator) term
    lists [(exponent tuple, (p, q))] with rational coefficients as integer pairs.
    Prime-independent, so the compile is cached and only the per-prime reduction
    (in _eval_terms) repeats."""
    num, den = spy.fraction(spy.together(expr))
    pos = {g: k for k, g in enumerate(gens)}
    def terms(e):
        # the Poly over the generators present only; a dense one over all of them
        # recurses through every level
        used = [g for g in gens if g in e.free_symbols]
        if not used:
            r = spy.Rational(e)
            return [(tuple([0] * len(gens)), (int(r.p), int(r.q)))]
        P = spy.Poly(e, *used)
        out = []
        for monom, coef in P.terms():
            r = spy.Rational(coef)
            full = [0] * len(gens)
            for g, m in zip(used, monom):
                full[pos[g]] = int(m)
            out.append((tuple(full), (int(r.p), int(r.q))))
        return out
    return terms(num), terms(den)


def _bipoly(expr, stateGens, paramGens):
    """Compile a polynomial in (states, params) to a list of (state-monomial
    exponent tuple, param-coefficient term list), the latter from _poly_terms, so
    the per-point state-polynomial is rebuilt by evaluating the parameter
    coefficients mod p."""
    P = spy.Poly(expr, *stateGens)
    out = []
    for sm, coef in P.terms():
        ct = _poly_terms(coef, paramGens)
        out.append((tuple(int(m) for m in sm), ct))
    return out


def _solve_mod(A, B, p):
    """Solve A X = B over GF(p) by Gauss-Jordan (A: n x n, B: n x k, integer lists
    already reduced mod p). Returns X as a list of rows, or None if A is singular."""
    n = len(A)
    if n == 0:
        return [[] for _ in range(0)]
    k = len(B[0]) if B and B[0] else 0
    M = [list(A[i]) + list(B[i]) for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if M[r][col] % p), None)
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        inv = pow(M[col][col] % p, p - 2, p)
        M[col] = [(x * inv) % p for x in M[col]]
        for r in range(n):
            if r != col and M[r][col] % p:
                f = M[r][col] % p
                M[r] = [(M[r][c] - f * M[col][c]) % p for c in range(n + k)]
    return [M[i][n:] for i in range(n)]


_invModCache = {}


def _inv_mod(a, p):
    """Modular inverse of a mod p, memoised per (a, p): term-list denominators are
    a small fixed set. inv(0)=0 as with Fermat."""
    a %= p
    if a <= 1:
        return a
    key = (a, p)
    v = _invModCache.get(key)
    if v is None:
        v = pow(a, p - 2, p)
        _invModCache[key] = v
    return v


def _eval_poly_terms(terms, ptvals, p):
    """Value of one term list (from _poly_terms) at the integer point `ptvals`
    (gens order, already mod p) as a residue mod p."""
    acc = 0
    for monom, (cn, cd) in terms:
        c = (cn % p) * _inv_mod(cd, p) % p          # cd is a compile-time constant
        for k, e in enumerate(monom):
            if e:
                c = c * pow(ptvals[k], e, p) % p
        acc = (acc + c) % p
    return acc


def _eval_terms(numden, ptvals, p):
    """Residue num/den mod p of a (num, den) term list at the point `ptvals`."""
    nv = _eval_poly_terms(numden[0], ptvals, p)
    dv = _eval_poly_terms(numden[1], ptvals, p)
    if dv == 1:                                     # polynomial entry: no inverse
        return nv
    return nv * pow(dv % p, p - 2, p) % p


def _eval_terms_guarded(numden, ptvals, p):
    """Residue num/den mod p, or None when the denominator vanishes mod p (a Fermat
    inverse of 0 would silently return 0)."""
    dv = _eval_poly_terms(numden[1], ptvals, p)
    if dv % p == 0:
        return None
    nv = _eval_poly_terms(numden[0], ptvals, p)
    if dv == 1:
        return nv
    return nv * pow(dv, p - 2, p) % p


def _reduce_modp(e, p):
    """Reduce the integer coefficients of a polynomial expression modulo p."""
    e = spy.sympify(e)
    syms = sorted(e.free_symbols, key=str)
    if not syms:
        return spy.Integer(int(e) % p)
    return spy.Poly(e, *syms, modulus=p).as_expr()


def evalRationalMod(expr, names, vals, q):
    """Value of the rational expression `expr` at the integer point names -> vals,
    reduced modulo the prime q. Returns None if the denominator vanishes mod q.
    Used to certify a reconstructed identifiability direction at a fresh prime."""
    q = int(q)
    if not isinstance(names, (list, tuple)): names = [names]
    if not isinstance(vals, (list, tuple)): vals = [vals]
    local, parse = _make_local_parse([str(expr)] + [str(n) for n in names])
    e = spy.together(spy.sympify(parse(str(expr))).subs(
        {spy.Symbol(str(n)): spy.Integer(int(v))
         for n, v in zip(names, vals)}))
    num, den = spy.fraction(e)
    num, den = int(num), int(den) % q
    if den == 0:
        return None
    return int((num % q) * pow(den, q - 2, q) % q)


_evalBatchCache = {}


def evalRationalModBatch(exprs, names, points, q):
    """Values of the rational expressions `exprs` (strings) at each integer point
    (rows of `points`, columns aligned with `names`), reduced modulo the prime q.
    The lambdified expressions are memoised on (exprs, names), since compiling
    dominates evaluation in the sampling loops; evaluation is exact (Fraction).
    Returns one list per point; -1 marks a vanishing denominator or a failure."""
    q = int(q)
    exprs = [str(e) for e in _as_list(exprs)]
    names = [str(n) for n in _as_list(names)]
    points = _as_list(points)
    if points and not isinstance(points[0], (list, tuple)):
        points = [points]
    key = (tuple(exprs), tuple(names))
    fns = _evalBatchCache.get(key)
    if fns is None:
        local, parse = _make_local_parse(exprs + names)
        argsyms = [local.get(n, spy.Symbol(n)) for n in names]
        fns = []
        for ex in exprs:
            try:
                fns.append(spy.lambdify(argsyms, parse(ex), modules="math"))
            except Exception:
                fns.append(None)
        if len(_evalBatchCache) > 512:
            _evalBatchCache.clear()
        _evalBatchCache[key] = fns
    out = []
    for pt in points:
        vals = [Fraction(int(v)) for v in pt]
        row = []
        for f in fns:
            v = None
            if f is not None:
                try:
                    v = Fraction(f(*vals))
                except Exception:
                    v = None
            if v is None or v.denominator % q == 0:
                row.append(-1)
            else:
                row.append(int(v.numerator % q) *
                           pow(v.denominator % q, q - 2, q) % q)
        out.append(row)
    return out


def evalRationalBatch(exprs, names, points):
    """Exact rational values of `exprs` at each integer point as decimal 'num'/'den'
    strings (den '0' marks a failure), for cofactor-matrix sampling. Shares the
    cache of evalRationalModBatch."""
    exprs = [str(e) for e in _as_list(exprs)]
    names = [str(n) for n in _as_list(names)]
    points = _as_list(points)
    if points and not isinstance(points[0], (list, tuple)):
        points = [points]
    key = (tuple(exprs), tuple(names))
    fns = _evalBatchCache.get(key)
    if fns is None:
        local, parse = _make_local_parse(exprs + names)
        argsyms = [local.get(n, spy.Symbol(n)) for n in names]
        fns = []
        for ex in exprs:
            try:
                fns.append(spy.lambdify(argsyms, parse(ex), modules="math"))
            except Exception:
                fns.append(None)
        if len(_evalBatchCache) > 512:
            _evalBatchCache.clear()
        _evalBatchCache[key] = fns
    nums, dens = [], []
    for pt in points:
        vals = [Fraction(int(v)) for v in pt]
        nrow, drow = [], []
        for f in fns:
            v = None
            if f is not None:
                try:
                    v = Fraction(f(*vals))
                except Exception:
                    v = None
            if v is None:
                nrow.append("0"); drow.append("0")
            else:
                nrow.append(str(v.numerator)); drow.append(str(v.denominator))
        nums.append(nrow); dens.append(drow)
    return {'num': nums, 'den': dens}


def _parse_generators(gens, extra_lines):
    """Shared symbol table + parsed xi dicts for a list of generator dicts."""
    gens = _as_list(gens)
    all_lines = list(extra_lines) + \
        [str(v) for g in gens for v in g.values()] + \
        [str(kk) for g in gens for kk in g.keys()]
    local, parse = _make_local_parse(all_lines)
    xis = [{str(kk): parse(str(v)) for kk, v in g.items()} for g in gens]
    return xis, parse


def moduleReduceReplay(rows, cols, schedule):
    """Exact sympy replay of a module-reduction pivot schedule, one call for the
    whole block. `rows`: per-generator dicts {col: expr-string} (missing = 0)
    over the coordinate list `cols`; `schedule`: 1-based [row, col] pivot pairs.
    Returns {'rows': [{col: string} nonzero entries], 'combo': [[strings]]}, the
    combination coefficients tracking each reduced row over the originals."""
    cols = [str(c) for c in _as_list(cols)]
    rows = _as_list(rows)
    k = len(rows)
    all_lines = [str(v) for r in rows for v in r.values()] + cols
    local, parse = _make_local_parse(all_lines)
    E = [{c: (parse(str(r[c])) if c in r and str(r[c]) != "0"
              else spy.Integer(0)) for c in cols} for r in rows]
    C = [[spy.Integer(1 if i == j else 0) for j in range(k)] for i in range(k)]
    for st in _as_list(schedule):
        r = int(st[0]) - 1
        cc = cols[int(st[1]) - 1]
        piv = E[r][cc]
        if piv == 0:
            continue
        for j in range(k):
            if j == r or E[j][cc] == 0:
                continue
            m = spy.cancel(E[j][cc] / piv)
            for c2 in cols:
                E[j][c2] = spy.cancel(E[j][c2] - m * E[r][c2])
            for l in range(k):
                C[j][l] = spy.cancel(C[j][l] - m * C[r][l])
    return {'rows': [{c: str(E[j][c]) for c in cols if E[j][c] != 0}
                     for j in range(k)],
            'combo': [[str(x) for x in row] for row in C]}


def darbouxCofactors(candidates, generators):
    """Batch Darboux division test with proportionality dedup. For each candidate
    polynomial string, X(P)/P must be polynomial for EVERY generator (dicts
    {coord: xi-string}); duplicates proportional to an earlier keeper are
    dropped. Returns {'Ps': [strings], 'lams': [[cofactor strings per gen]]}."""
    cands = [str(c) for c in _as_list(candidates)]
    xis, parse = _parse_generators(generators, cands)
    Ps, lams, kept = [], [], []
    for cs in cands:
        try:
            P = parse(cs)
        except Exception:
            continue
        if getattr(P, "is_number", False):
            continue
        lam = []
        ok = True
        for xi in xis:
            XP = spy.Integer(0)
            for v, comp in xi.items():
                XP = XP + comp * spy.diff(P, spy.Symbol(v))
            q = spy.cancel(spy.expand(XP) / P)
            if str(spy.fraction(q)[1]) != "1":
                ok = False
                break
            lam.append(str(q))
        if not ok:
            continue
        if any(getattr(spy.cancel(P / Q), "is_number", False) for Q in kept):
            continue
        kept.append(P)
        Ps.append(str(P))
        lams.append(lam)
    return {'Ps': Ps, 'lams': lams}


def verifyInvariants(invariants, generators):
    """Batch exact verification X(I) = 0: for each invariant string, the Lie
    derivative along every generator must cancel to literal 0 (a proof, not a
    sample). Returns a list of booleans."""
    invs = [str(i) for i in _as_list(invariants)]
    xis, parse = _parse_generators(generators, invs)
    out = []
    def _real(e):
        return e.subs({t: spy.Symbol(t.name, real=True) for t in e.free_symbols})

    for s in invs:
        try:
            # real symbols: dMod coordinates are real, and the atan/Abs forms of
            # the quadrature stage only cancel under that assumption
            I = _real(parse(s))
            ok = True
            for xi in xis:
                XI = spy.Integer(0)
                for v, comp in xi.items():
                    XI = XI + _real(comp) * spy.diff(I, spy.Symbol(v, real=True))
                XI = spy.cancel(spy.together(XI))
                if str(XI) != "0" and spy.simplify(XI) != 0:
                    # cancel() suffices for rational/exp invariants; the simplify
                    # fallback covers atan/Abs forms from the quadrature stage
                    ok = False
                    break
            out.append(ok)
        except Exception:
            out.append(False)
    return out


def integratingFactorIntegral(Mstr, xi1, xi2, z1, z2):
    """First integral by quadrature from a 2D integrating factor:
    dI = M*(xi2 dz1 - xi1 dz2), so I = int M*xi2 dz1 + phi(z2) with phi fixed by
    the z2-part. Returns the integral string, or None when the quadrature stays
    unevaluated or the exactness check X(I) = 0 fails."""
    local, parse = _make_local_parse([str(Mstr), str(xi1), str(xi2),
                                      str(z1), str(z2)])
    try:
        M = parse(str(Mstr))
        f1 = parse(str(xi1))
        f2 = parse(str(xi2))
        # real symbols throughout: ratint then prefers real atan/log forms, and
        # a leftover complex-logarithm pair can be realified by re()
        real = {s: spy.Symbol(s.name, real=True)
                for s in (M.free_symbols | f1.free_symbols | f2.free_symbols)}
        M = M.subs(real)
        f1 = f1.subs(real)
        f2 = f2.subs(real)
        s1 = spy.Symbol(str(z1), real=True)
        s2 = spy.Symbol(str(z2), real=True)
        F = spy.integrate(M * f2, s1)
        rest = spy.cancel(spy.together(-M * f1 - spy.diff(F, s2)))
        phi = spy.integrate(rest, s2)
        I = F + phi
        if I.has(spy.Integral):
            return None
        if I.has(spy.I):
            # the imaginary remainder of a real first integral is constant
            I = spy.simplify(spy.re(spy.expand_complex(I)))
            if I.has(spy.I):
                return None
        XI = spy.together(f1 * spy.diff(I, s1) + f2 * spy.diff(I, s2))
        if spy.simplify(XI) != 0:
            return None
        return str(I)
    except Exception:
        return None


_LINSOLVE_OPS_CAP = 2000  # bail to the numeric solver past this elimination size


def _linear_solution(polys, solveStates):
    """Solve polys = 0 for `solveStates` by generic linear elimination, each state a
    rational function of the parameters. Returns {state: expr}, or None when a
    coupled residual remains or an expression exceeds `_LINSOLVE_OPS_CAP`. A point
    where a pivot denominator vanishes mod p is caught at evaluation."""
    remSet = set(solveStates)
    remP = [spy.sympify(pl) for pl in polys]
    elim = []
    progress = True
    while progress and remSet:
        progress = False
        for idx in range(len(remP)):
            pl = remP[idx]
            if pl is None or not (pl.free_symbols & remSet):
                continue
            if pl.count_ops() > _LINSOLVE_OPS_CAP:
                return None
            picked = None
            for v in pl.free_symbols & remSet:
                pv = spy.Poly(pl, v)
                if pv.degree() != 1:
                    continue
                a = pv.nth(1)
                if (a.free_symbols & remSet) or a.is_zero:
                    continue
                picked = (v, spy.cancel(-pv.nth(0) / a))
                break
            if picked is None:
                continue
            v, expr = picked
            elim.append((v, expr))
            remP[idx] = None
            for j in range(len(remP)):
                if remP[j] is not None and v in remP[j].free_symbols:
                    sub = remP[j].subs(v, expr)
                    if sub.count_ops() > _LINSOLVE_OPS_CAP:
                        return None
                    remP[j] = spy.fraction(spy.together(sub))[0]
            remSet.discard(v)
            progress = True
    if remSet:
        return None
    sol = {}
    for v, expr in reversed(elim):
        sol[v] = spy.cancel(expr.subs(sol))
    return sol


def _compile_linear_plan(polys, solveStates, paramSyms):
    """Compile the linear-elimination solution of polys = 0 to per-state
    prime-independent (num, den) term lists in `solveStates` order. Returns
    (True, perStateTerms) when every state eliminates linearly, else (False, None)."""
    sol = _linear_solution(polys, solveStates)
    if sol is None:
        return False, None
    return True, [(str(s), _poly_terms(sol[s], paramSyms)) for s in solveStates]


# ---- fast numeric state solve (dict polynomials mod p, no sympy in the hot loop) ----
# Linear elimination on dict polynomials; a coupled core goes to _solve_states_reduced
# (up to two states) or to msolve.

def _eval_bipoly_dict(bip, paramvals, p):
    """The state polynomial at numeric parameter values mod p, as a dict
    {state-exponent-tuple: coeff}."""
    d = {}
    for sm, ct in bip:
        c = _eval_terms(ct, paramvals, p)
        if c:
            d[sm] = (d.get(sm, 0) + c) % p
    return {m: c for m, c in d.items() if c}


def _dp_mul(a, b, nv, p):
    """Product of two dict polynomials mod p."""
    out = {}
    for ma, ca in a.items():
        for mb, cb in b.items():
            key = tuple(ma[k] + mb[k] for k in range(nv))
            v = (out.get(key, 0) + ca * cb) % p
            if v:
                out[key] = v
            elif key in out:
                del out[key]
    return out


def _dp_subst(d, i, expr, nv, p):
    """Substitute var i -> expr (a dict poly not containing var i) into d, mod p."""
    maxe = max((m[i] for m in d), default=0)
    pows = [{tuple([0] * nv): 1}]
    for _ in range(maxe):
        pows.append(_dp_mul(pows[-1], expr, nv, p))
    out = {}
    for m, c in d.items():
        term = pows[m[i]]
        base = tuple(0 if k == i else m[k] for k in range(nv))
        for mm, cc in term.items():
            key = tuple(base[k] + mm[k] for k in range(nv))
            v = (out.get(key, 0) + c * cc) % p
            if v:
                out[key] = v
            elif key in out:
                del out[key]
    return out


def _solve_states_fast(dps, solveStates, p):
    """Interior point of f = 0 over GF(p) from the per-point state polynomials in dict
    form. Linear elimination by dict-polynomial substitution mod p, then numeric
    back-substitution; a coupled residual goes to _solve_states_reduced (core of at most
    two states) or _solve_states_msolve. Returns (sol {Symbol: Integer}, None) or
    (None, fail-dict)."""
    nv = len(solveStates)
    dps0 = dps
    dps = [dict(d) for d in dps]
    remIdx = set(range(nv))
    elim = []                          # (var_index, expr dict poly over remaining vars)
    progress = True
    while progress and remIdx:
        progress = False
        for eqi in range(len(dps)):
            d = dps[eqi]
            if d is None or not d:
                continue
            present = sorted((i for i in remIdx if any(m[i] for m in d)),
                             key=lambda i: str(solveStates[i]))   # match the symbolic order
            if not present:
                continue
            picked = None
            for i in present:
                if max(m[i] for m in d) != 1:          # degree in var i must be 1
                    continue
                c1 = 0
                ok = True
                for m in d:                            # coeff of var_i^1 must be const
                    if m[i] == 1:
                        if any(m[k] for k in range(nv) if k != i):
                            ok = False
                            break
                        c1 = d[m]
                if not ok or c1 % p == 0:
                    continue
                inv = pow(c1 % p, p - 2, p)
                c0 = {m: c for m, c in d.items() if m[i] == 0}   # var_i^0 part
                expr = {m: (c * (-inv)) % p for m, c in c0.items()}
                expr = {m: c for m, c in expr.items() if c}
                picked = (i, expr)
                break
            if picked is None:
                continue
            i, expr = picked
            elim.append((i, expr))
            dps[eqi] = None
            for j in range(len(dps)):
                if dps[j] and any(m[i] for m in dps[j]):
                    dps[j] = _dp_subst(dps[j], i, expr, nv, p)
            remIdx.discard(i)
            progress = True

    sol = {}
    coupledIdx = [i for i in range(nv) if i in remIdx]
    if coupledIdx:
        red = _solve_states_reduced(dps0, solveStates, p)
        if red is not None:
            return red
        return _solve_states_msolve(dps0, solveStates, p)

    valById = {i: int(sol[solveStates[i]]) % p for i in coupledIdx}
    for i, expr in reversed(elim):     # resolve eliminated vars numerically, latest first
        acc = 0
        for m, c in expr.items():
            t = c
            for k in range(nv):
                if m[k]:
                    t = t * pow(valById.get(k, 0), m[k], p) % p
            acc = (acc + t) % p
        valById[i] = acc
        sol[solveStates[i]] = spy.Integer(acc)
    return sol, None


# ---- reduced coupled solve: substitution to a small core, resultant, univariate roots ----
# A state that enters a single balance, linearly, is an output of the network: it is
# peeled and solved last. The rest is eliminated by linear pivots, rational where the
# pivot coefficient depends on other states (x = -b/a, the other balances multiplied by
# a^deg). A core of at most two states is solved by a resultant and root finding over
# GF(p); every candidate is back-substituted and checked against f = 0.

_REDUCE_TERMS_CAP = 2000    # past this polynomial size msolve is faster than the reduction


def _dp_add(a, b, p):
    """Sum of two dict polynomials mod p."""
    out = dict(a)
    for m, c in b.items():
        v = (out.get(m, 0) + c) % p
        if v:
            out[m] = v
        elif m in out:
            del out[m]
    return out


def _dp_coeffs(d, i):
    """Coefficients {k: dict poly free of var i} of d as a polynomial in var i."""
    cs = {}
    for m, c in d.items():
        cs.setdefault(m[i], {})[m[:i] + (0,) + m[i + 1:]] = c
    return cs


def _dp_subst_rat(d, i, a, b, nv, p):
    """d at var i = -b/a, times a^deg_i(d): a dict polynomial free of var i."""
    cs = _dp_coeffs(d, i)
    deg = max(cs)
    one = {tuple([0] * nv): 1}
    nb = {m: (-c) % p for m, c in b.items()}
    powA, powB = [one], [one]
    for _ in range(deg):
        powA.append(_dp_mul(powA[-1], a, nv, p))
        powB.append(_dp_mul(powB[-1], nb, nv, p))
    out = {}
    for k, ck in cs.items():
        out = _dp_add(out, _dp_mul(_dp_mul(ck, powB[k], nv, p), powA[deg - k], nv, p), p)
    return out


def _dp_value(d, vals, p):
    """Value of a dict polynomial at the point `vals` (residues mod p)."""
    acc = 0
    for m, c in d.items():
        t = c
        for k, e in enumerate(m):
            if e:
                t = t * pow(vals[k], e, p) % p
        acc = (acc + t) % p
    return acc


def _reduce_system(dps, nv, p):
    """Peel outputs and eliminate linear pivots. Returns (plan, sinks, resid, core):
    plan and sinks as (var, a, b) with var = -b/a, in elimination order; resid the
    remaining balances over the core states. None past _REDUCE_TERMS_CAP."""
    dps = [dict(d) for d in dps]
    alive = [j for j in range(len(dps)) if dps[j]]
    rem = set(range(nv))
    sinks = []
    changed = True
    while changed:
        changed = False
        for i in sorted(rem):
            js = [j for j in alive if any(m[i] for m in dps[j])]
            if len(js) == 1 and max(m[i] for m in dps[js[0]]) == 1:
                cs = _dp_coeffs(dps[js[0]], i)
                sinks.append((i, cs[1], cs.get(0, {})))
                alive.remove(js[0])
                rem.discard(i)
                changed = True
                break
    plan = []
    while True:
        # max degree of each remaining state per balance, and in how many it occurs
        degOf = {}
        for j in alive:
            dg = {}
            for m in dps[j]:
                for i in rem:
                    if m[i] > dg.get(i, 0):
                        dg[i] = m[i]
            degOf[j] = dg
        occ = {i: sum(1 for j in alive if i in degOf[j]) for i in rem}
        # constant pivot coefficients first, then a state's own balance, then the
        # pivot touching the fewest balances with the lowest coefficient degree
        best = None
        for j in alive:
            for i, dg in degOf[j].items():
                if dg != 1:
                    continue
                cs = _dp_coeffs(dps[j], i)
                a = cs[1]
                key = (0 if all(not any(m) for m in a) else 1, 0 if i == j else 1,
                       max(sum(m) for m in a) + occ[i] - 1,
                       len(a) + len(cs.get(0, {})), i, j)
                if best is None or key < best[0]:
                    best = (key, j, i, a, cs.get(0, {}))
        if best is None:
            break
        _, j, i, a, b = best
        plan.append((i, a, b))
        alive.remove(j)
        rem.discard(i)
        for k in alive:
            if any(m[i] for m in dps[k]):
                dps[k] = _dp_subst_rat(dps[k], i, a, b, nv, p)
                if len(dps[k]) > _REDUCE_TERMS_CAP:
                    return None
    return plan, sinks, [dps[j] for j in alive if dps[j]], sorted(rem)


def _gf_roots(f, p):
    """Distinct nonzero roots in GF(p), ascending, of the dense univariate polynomial
    f (highest degree first)."""
    f = _gft.gf_strip([int(c) % p for c in f])
    if len(f) <= 1:
        return []
    f = _gft.gf_monic(f, p, _ZZ)[1]
    h = _gft.gf_sub(_gft.gf_pow_mod([1, 0], p, f, p, _ZZ), [1, 0], p, _ZZ)
    g = _gft.gf_gcd(f, h, p, _ZZ)
    if len(g) <= 1:
        return []
    lin = [g] if len(g) == 2 else _gft.gf_edf_zassenhaus(g, 1, p, _ZZ)
    return sorted(r for r in ((-fac[1]) % p for fac in lin) if r)


def _gf_resultant(f, g, p):
    """Resultant of two dense univariate polynomials over GF(p) by Euclid."""
    f, g = _gft.gf_strip(f), _gft.gf_strip(g)
    if not f or not g:
        return 0
    res = 1
    while True:
        df, dg = len(f) - 1, len(g) - 1
        if dg == 0:
            return res * pow(g[0], df, p) % p
        r = _gft.gf_rem(f, g, p, _ZZ)
        if not r:
            return 0
        if df % 2 and dg % 2:
            res = -res
        res = res * pow(g[0], df - (len(r) - 1), p) % p
        f, g = g, r


def _gf_interpolate(xs, ys, p):
    """Dense polynomial (highest degree first) through the points (xs, ys) over GF(p),
    by Newton's divided differences."""
    n = len(xs)
    c = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) * pow((xs[i] - xs[i - j]) % p, p - 2, p) % p
    poly = [c[-1]]
    for i in range(n - 2, -1, -1):
        poly = _gft.gf_add(_gft.gf_mul(poly, [1, (-xs[i]) % p], p, _ZZ), [c[i]], p, _ZZ)
    return _gft.gf_strip(poly)


def _bi_at_u(P, a, p):
    """Dense polynomial in v of the bivariate dict polynomial P{(du, dv): c} at u = a."""
    dv = max(k[1] for k in P)
    out = [0] * (dv + 1)
    for (i, j), c in P.items():
        out[dv - j] = (out[dv - j] + c * pow(a, i, p)) % p
    return _gft.gf_strip(out)


def _core_candidates(resid, core, p):
    """Candidate values (tuples aligned with `core`) of the core states, or None when
    the core has more than two states."""
    if not core:
        return [()] if not resid else []
    if len(core) > 2 or len(resid) != len(core):
        return None
    if len(core) == 1:
        i = core[0]
        g = None
        for d in resid:
            q = [0] * (max(m[i] for m in d) + 1)
            for m, c in d.items():
                q[len(q) - 1 - m[i]] = c
            q = _gft.gf_strip(q)
            g = q if g is None else _gft.gf_gcd(g, q, p, _ZZ)
        return [(r,) for r in _gf_roots(g, p)]
    # two states (u, v): the resultant in v by evaluation at u = a and interpolation.
    # A common factor g of the balances (a cleared denominator) is divided out at each
    # a through the gcd in v; the values stay polynomial in a (times a power of lc_v g)
    iu, iv = core
    P1, P2 = [{(m[iu], m[iv]): c for m, c in d.items()} for d in resid]
    if max(k[1] for k in P1) < 1 or max(k[1] for k in P2) < 1:
        return None
    bound = max(sum(k) for k in P1) * max(sum(k) for k in P2)
    rng = random.Random(p)
    dg = min(len(_gft.gf_gcd(_bi_at_u(P1, a, p), _bi_at_u(P2, a, p), p, _ZZ)) - 1
             for a in (rng.randrange(1, p), rng.randrange(1, p)))
    n1, n2 = max(k[1] for k in P1), max(k[1] for k in P2)
    xs, ys = [], []
    a = 0
    while len(xs) <= bound:
        a += 1
        f1, f2 = _bi_at_u(P1, a, p), _bi_at_u(P2, a, p)
        if len(f1) - 1 != n1 or len(f2) - 1 != n2:
            continue
        G = _gft.gf_gcd(f1, f2, p, _ZZ)
        if len(G) - 1 > dg:
            r = 0                      # a common root on top of g: a root of the resultant
        else:
            r = _gf_resultant(_gft.gf_quo(f1, G, p, _ZZ), _gft.gf_quo(f2, G, p, _ZZ), p)
        xs.append(a)
        ys.append(r)
    res = _gf_interpolate(xs, ys, p)
    if not res:
        return None
    out = []
    for a in _gf_roots(res, p):
        G = _gft.gf_gcd(_bi_at_u(P1, a, p), _bi_at_u(P2, a, p), p, _ZZ)
        out.extend((a, b) for b in _gf_roots(G, p))
    return out


def _solve_states_reduced(dps, solveStates, p):
    """Interior point of f = 0 over GF(p) by _reduce_system and a resultant on the
    core. Returns (sol, None), (None, fail-dict) when no candidate is an interior
    point, or None when the core is too large for this solve."""
    nv = len(solveStates)
    red = _reduce_system(dps, nv, p)
    if red is None:
        return None
    plan, sinks, resid, core = red
    try:
        cands = _core_candidates(resid, core, p)
    except Exception:
        return None
    if cands is None:
        return None
    steps = list(reversed(plan)) + list(reversed(sinks))
    for cp in cands:
        vals = [0] * nv
        for i, x in zip(core, cp):
            vals[i] = x
        ok = True
        for i, a, b in steps:
            av = _dp_value(a, vals, p)
            if av == 0:
                ok = False
                break
            vals[i] = (-_dp_value(b, vals, p)) * pow(av, p - 2, p) % p
            if vals[i] == 0:
                ok = False
                break
        if ok and all(_dp_value(d, vals, p) == 0 for d in dps):
            return {solveStates[i]: spy.Integer(vals[i]) for i in range(nv)}, None
    return None, {'ok': False, 'why': 'no consistent interior point'}


# ---- coupled solve by msolve ----------------------------------------------------------
# Constant pivots first (rational pivots add spurious components); the rest goes to
# msolve, saturated by w * prod(x) - 1, as x_i = -p_i(t) / d(t) over the roots of h(t).

_MSOLVE_TIMEOUT = 600
_MSOLVE_PATH = ""        # set by R


class MsolveMissing(RuntimeError):
    pass


def _msolve_binary():
    """Path of the msolve executable."""
    path = _MSOLVE_PATH or os.environ.get("DMOD_MSOLVE", "")
    if not path or not os.path.isfile(path):
        raise MsolveMissing(
            "symmetryDetection(equilibrate = TRUE): this steady state needs msolve; "
            "see dMod2::install_libs(\"msolve\").")
    return path


def _const_reduce(dps, nv, p):
    """Eliminate pivots x_i = -b/a0 with constant a0. Returns (plan, rest), plan as
    (i, b, 1/a0) in elimination order."""
    dps = [dict(d) for d in dps if d]
    plan = []
    while True:
        pick = None
        for j, d in enumerate(dps):
            for i in range(nv):
                if not any(m[i] for m in d) or max(m[i] for m in d) != 1:
                    continue
                cs = _dp_coeffs(d, i)
                a = cs[1]
                if len(a) == 1 and not any(next(iter(a))):
                    pick = (j, i, cs.get(0, {}), pow(next(iter(a.values())) % p, p - 2, p))
                    break
            if pick:
                break
        if pick is None:
            return plan, dps
        j, i, b, inv = pick
        dps.pop(j)
        plan.append((i, b, inv))
        a1 = {tuple([0] * nv): 1}
        b1 = {m: c * inv % p for m, c in b.items()}      # x_i = -(b * inv)
        dps = [(_dp_subst_rat(d, i, a1, b1, nv, p) if any(m[i] for m in d) else d) for d in dps]
        dps = [d for d in dps if d]


def _msolve_points(polys, vars_idx, nv, p):
    """Interior points of `polys` over GF(p). Returns a list of {index: value} or a
    fail-dict."""
    import ast, subprocess, tempfile
    names = {i: "x%d" % i for i in vars_idx}

    def pstr(d):
        out = []
        for m, c in d.items():
            mono = "*".join(names[i] + ("^%d" % m[i] if m[i] > 1 else "") for i in vars_idx if m[i])
            out.append("%d*%s" % (c % p, mono) if mono else "%d" % (c % p))
        return "+".join(out) or "0"

    eqs = [pstr(d) for d in polys] + ["w*" + "*".join(names[i] for i in vars_idx) + "-1"]
    vs = [names[i] for i in vars_idx] + ["w"]
    with tempfile.TemporaryDirectory() as td:
        fi, fo = os.path.join(td, "in.ms"), os.path.join(td, "out.ms")
        with open(fi, "w") as fh:
            fh.write(", ".join(vs) + "\n" + str(p) + "\n" + ",\n".join(eqs) + "\n")
        r = subprocess.run([_msolve_binary(), "-f", fi, "-o", fo, "-P", "1", "-t", "1"],
                           capture_output=True, text=True, timeout=_MSOLVE_TIMEOUT)
        txt = open(fo).read() if os.path.exists(fo) else ""
    if r.returncode != 0 or "too large" in (r.stdout + r.stderr) or not txt.strip():
        return {'ok': False, 'why': 'msolve failed: %s' % (r.stdout + r.stderr).strip()[-200:]}
    data = ast.literal_eval(txt.strip().rstrip(":"))
    if data[0] != 0:
        return {'ok': False, 'why': 'positive-dimensional steady state'}
    _, nvar, deg, vnames, lin, par = data[1]
    if deg <= 0 or par[0] != 1:
        return {'ok': False, 'why': 'no consistent interior point'}
    elim, den, params = par[1][0][1], par[1][1][1], par[1][2]
    lin = [int(c) % p for c in lin]
    last = len(vnames) - 1
    if lin[last] == 0:
        return {'ok': False, 'why': 'msolve: unexpected linear form'}
    idxOf = {names[i]: i for i in vars_idx}

    def ev(coeffs, x):          # coefficients lowest degree first
        acc = 0
        for c in reversed(coeffs):
            acc = (acc * x + c) % p
        return acc

    out = []
    for t in _gf_roots(list(reversed([int(c) % p for c in elim])), p):
        dv = ev(den, t)
        if dv == 0:
            continue
        dinv = pow(dv, p - 2, p)
        vals = {}
        for k, q in enumerate(params):
            vals[vnames[k]] = (-ev(q[0][1], t)) * dinv % p
        # last variable from the linear form
        rest = (t - sum(lin[k] * vals[vnames[k]] for k in range(last))) % p
        vals[vnames[last]] = rest * pow(lin[last], p - 2, p) % p
        out.append({idxOf[nm]: v for nm, v in vals.items() if nm in idxOf})
    return out


def _solve_states_msolve(dps, solveStates, p):
    """Interior point of f = 0 over GF(p), checked against f."""
    nv = len(solveStates)
    plan, rest = _const_reduce(dps, nv, p)
    # outputs: a state in a single remaining balance, linearly, is solved after msolve
    sinks = []
    while True:
        occ = {}
        for j, d in enumerate(rest):
            for i in {i for m in d for i in range(nv) if m[i]}:
                occ.setdefault(i, []).append(j)
        pick = next(((i, js[0]) for i, js in sorted(occ.items()) if len(js) == 1 and
                     max(m[i] for m in rest[js[0]]) == 1 and len(occ) > 1), None)
        if pick is None:
            break
        i, j = pick
        cs = _dp_coeffs(rest[j], i)
        sinks.append((i, cs[1], cs.get(0, {})))
        rest = rest[:j] + rest[j + 1:]
    varsLeft = sorted({i for d in rest for m in d for i in range(nv) if m[i]})
    if rest and not varsLeft:
        return None, {'ok': False, 'why': 'no consistent interior point'}
    cands = _msolve_points(rest, varsLeft, nv, p) if rest else [{}]
    if isinstance(cands, dict):
        return None, cands
    for cand in cands:
        vals = [0] * nv
        for i, v in cand.items():
            vals[i] = v
        ok = True
        for i, a, b in reversed(sinks):
            av = _dp_value(a, vals, p)
            if av == 0:
                ok = False
                break
            vals[i] = (-_dp_value(b, vals, p)) * pow(av, p - 2, p) % p
        if not ok:
            continue
        for i, b, inv in reversed(plan):
            vals[i] = (-_dp_value(b, vals, p)) * inv % p
        if all(vals) and all(_dp_value(d, vals, p) == 0 for d in dps):
            return {solveStates[i]: spy.Integer(vals[i]) for i in range(nv)}, None
    return None, {'ok': False, 'why': 'no consistent interior point'}


def _compile_t0events(events, paramNames):
    """Compile each t0 event's value and its parameter derivatives to
    prime-independent term lists over the parameters, so a point evaluates by
    integer arithmetic. Keyed by the event tuple and parameter order."""
    key = (tuple((str(e['var']), str(e['value']), str(e['method'])) for e in events),
           tuple(paramNames))
    compiled = _eventCompileCache.get(key)
    if compiled is not None:
        return compiled
    paramSyms = [spy.Symbol(nm) for nm in paramNames]
    eloc = _build_symbol_table([str(e['value']) for e in events] + list(paramNames))
    eloc.update(_function_aliases())
    eparse = _make_parse(eloc)
    compiled = []
    for evt in events:
        val = spy.sympify(eparse(str(evt['value'])))
        duals = {}
        for nm in paramNames:
            d = spy.diff(val, spy.Symbol(nm))
            duals[nm] = None if d == 0 else _poly_terms(d, paramSyms)
        compiled.append({'var': str(evt['var']), 'method': str(evt['method']),
                         'valT': _poly_terms(val, paramSyms), 'duals': duals})
    _eventCompileCache[key] = compiled
    return compiled


def _ss_compile(model, stateNames, paramNames, forcings, heldStateNames=()):
    """Prime-independent compile of f = 0 for the modular steady-state solve,
    memoised in _ssModularCache. Inputs are the already-normalised model/state/param
    lists and the forcings set. Returns the cached tuple (paramSyms, solveStates,
    polys, Jx, Jt, gens, JxTerms, JtTerms, polyBi, genericLinear, linTerms, pointPlan).

    `heldStateNames` are conserved-moiety pivot states with a free resting value
    (reduceCQ = FALSE). f = 0 is rank-deficient by one equation per moiety, so these
    states act as parameters of the point solve and their dependent equations are
    dropped, leaving a square system (`pointPlan`, None when nothing is held).
    `solveStates` and Jx/Jt still cover all non-forcing states, so the joint
    determining system leaves the pivot directions free."""
    key = (tuple(model), tuple(stateNames), tuple(paramNames), tuple(sorted(forcings)),
           tuple(sorted(heldStateNames)))
    cached = _ssModularCache.get(key)
    if cached is not None:
        return cached
    local, parse = _make_local_parse(model + list(paramNames))
    rhsByName = {}
    for raw in model:
        line = _clean(raw)
        if '=' not in line:
            continue
        lhs, rhs = line.split('=', 1)
        rhsByName[lhs.strip()] = parse(rhs)
    stateSyms = [spy.Symbol(nm) for nm in stateNames]
    paramSyms = [spy.Symbol(nm) for nm in paramNames]
    solveStates = [s for s in stateSyms if str(s) not in forcings]
    forcingSubs = {spy.Symbol(nm): spy.Integer(0) for nm in forcings}
    fSym = [spy.sympify(rhsByName[str(s)]).subs(forcingSubs) for s in solveStates]
    polys = [spy.expand(spy.fraction(spy.together(fi))[0]) for fi in fSym]
    Jx = spy.Matrix([[spy.diff(fi, sj) for sj in solveStates] for fi in fSym])
    Jt = {str(th): spy.Matrix([[spy.diff(fi, th)] for fi in fSym])
          for th in paramSyms}
    gens = list(paramSyms) + list(solveStates)
    nSc = len(solveStates)
    JxTerms = [[_poly_terms(Jx[i, j], gens) for j in range(nSc)]
               for i in range(nSc)]
    JtTerms = {str(th): [_poly_terms(Jt[str(th)][i], gens) for i in range(nSc)]
               for th in paramSyms}
    polyBi = [_bipoly(pl, list(solveStates), list(paramSyms)) for pl in polys]
    genericLinear, linTerms = _compile_linear_plan(polys, solveStates, paramSyms)
    # held states are a moiety pivot set (R: .cq_pivot_decomposition), so their own
    # equations are exactly the dependent rows of f = 0
    heldSet = set(heldStateNames)
    pointPlan = None
    if heldSet:
        heldSyms = [s for s in solveStates if str(s) in heldSet]
        pointStates = [s for s in solveStates if str(s) not in heldSet]
        pointParams = list(paramSyms) + heldSyms
        keepPolys = [polys[i] for i, s in enumerate(solveStates) if str(s) not in heldSet]
        pointPolyBi = [_bipoly(pl, list(pointStates), list(pointParams)) for pl in keepPolys]
        pointGenericLinear, pointLinTerms = _compile_linear_plan(
            keepPolys, pointStates, pointParams)
        pointPlan = (pointStates, pointParams, heldSyms, pointPolyBi,
                     pointGenericLinear, pointLinTerms)
    cached = (paramSyms, solveStates, polys, Jx, Jt, gens, JxTerms, JtTerms,
              polyBi, genericLinear, linTerms, pointPlan)
    _ssModularCache[key] = cached
    return cached


def solveSteadyStateModular(model, stateNames, paramNames, paramVals, prime,
                            forcings=None, backend='sympy', t0events=None,
                            recast=None, lVals=None, jointMode=False,
                            heldStates=None):
    """Numeric point on the interior component of f = 0 over GF(prime) with its
    implicit-function-theorem parameter sensitivities.

    `model` is a list of "X = rhs" lines; `stateNames`/`paramNames` the state and
    parameter order; `paramVals` a {name: int} map of residues mod `prime`;
    `forcings` are held at 0; `t0events` are dose-style events composed onto the
    point after the solve. Returns {'ok': True, 'xstar': [...], 'dx': {param:
    [...]}} (state-ordered) or {'ok': False, 'why': ...}."""
    _select_backend(backend)
    model = _as_list(model)
    forcings = set(_as_list(forcings))
    stateNames = _as_list(stateNames)
    paramNames = _as_list(paramNames)
    p = int(prime)

    heldStates = {str(k): int(v) % p for k, v in (heldStates or {}).items()}
    cached = _ss_compile(model, stateNames, paramNames, forcings,
                         tuple(sorted(heldStates.keys())))
    (paramSyms, solveStates, polys, Jx, Jt, gens, JxTerms, JtTerms, polyBi,
     genericLinear, linTerms, pointPlan) = cached

    paramvals = [int(paramVals.get(str(th), 0)) % p for th in paramSyms]

    if pointPlan is not None:
        # held pivot states act as parameters; solve the square reduced system
        (pointStates, pointParams, heldSyms, pointPolyBi,
         pointGenericLinear, pointLinTerms) = pointPlan
        ppvals = [(heldStates[str(th)] if str(th) in heldStates
                   else int(paramVals.get(str(th), 0)) % p) for th in pointParams]
        sub = None
        if pointGenericLinear and not _SS_FORCE_POINT:
            sub = {}
            for nm, terms in pointLinTerms:
                val = _eval_terms_guarded(terms, ppvals, p)
                if val is None:
                    sub = None
                    break
                sub[spy.Symbol(nm)] = spy.Integer(val)
        if sub is None:
            dps = [_eval_bipoly_dict(bip, ppvals, p) for bip in pointPolyBi]
            sub, fail = _solve_states_fast(dps, list(pointStates), p)
            if fail is not None:
                return fail
        sol = {}
        for s in solveStates:
            sol[s] = (spy.Integer(heldStates[str(s)]) if str(s) in heldStates
                      else spy.Integer(int(sub[s]) % p))
    else:
        # fast path: each state is a compiled rational function of the parameters,
        # evaluated in integer arithmetic; a vanishing pivot denominator mod p (None)
        # routes the point to the symbolic solve
        sol = None
        if genericLinear and not _SS_FORCE_POINT:
            sol = {}
            for nm, terms in linTerms:
                val = _eval_terms_guarded(terms, paramvals, p)
                if val is None:
                    sol = None
                    break
                sol[spy.Symbol(nm)] = spy.Integer(val)
        if sol is None:
            dps = [_eval_bipoly_dict(bip, paramvals, p) for bip in polyBi]
            sol, fail = _solve_states_fast(dps, list(solveStates), p)
            if fail is not None:
                return fail

    nS = len(solveStates)
    valBy = {str(s): int(sol[s]) % p for s in solveStates}
    sIdx = {str(s): i for i, s in enumerate(solveStates)}
    # point in gens order (params then solve-states), already reduced mod p
    ptvals = [int(paramVals.get(str(th), 0)) % p for th in paramSyms] + \
             [valBy[str(s)] for s in solveStates]
    JtBy = {str(th): [_eval_terms(JtTerms[str(th)][i], ptvals, p) for i in range(nS)]
            for th in paramSyms}
    Jxeff = [[_eval_terms(JxTerms[i][j], ptvals, p) for j in range(nS)]
             for i in range(nS)]

    # resting Jacobian df_rest for the joint determining system (tangency
    # df_rest . xi = 0 in uneliminated (x, theta) coordinates), snapshot before the
    # recast folding below mutates Jxeff/JtBy; R stacks these rows
    if pointPlan is not None:
        # reduced tangency: non-pivot rows and state columns; pivot state columns
        # become parameter columns, as the pivot's resting value is free
        heldIdx = [i for i, s in enumerate(solveStates) if str(s) in heldStates]
        keepIdx = [i for i, s in enumerate(solveStates) if str(s) not in heldStates]
        dfJx = [[Jxeff[i][j] for j in keepIdx] for i in keepIdx]
        dfStateCols = [str(solveStates[j]) for j in keepIdx]
        dfJt = {str(th): [JtBy[str(th)][i] for i in keepIdx] for th in paramSyms}
        for hj in heldIdx:
            dfJt[str(solveStates[hj])] = [Jxeff[i][hj] for i in keepIdx]
        dfParamCols = [str(th) for th in paramSyms] + \
                      [str(solveStates[hj]) for hj in heldIdx]
    else:
        dfJx = [list(row) for row in Jxeff]
        dfJt = {nm: list(col) for nm, col in JtBy.items()}
        dfStateCols = [str(sst) for sst in solveStates]
        dfParamCols = [str(th) for th in paramSyms]

    if jointMode:
        # joint mode needs only valBy and df_rest; returning here also skips the
        # (occasionally singular) dual solve, which would reject a good point
        return {'ok': True, 'stateNames': list(stateNames), 'valBy': dict(valBy),
                'dfJx': dfJx, 'dfJt': dfJt, 'dfStateCols': dfStateCols,
                'dfParamCols': dfParamCols}

    # power/Hill recast. A normal entry holds E = base^exp generic and folds its chain
    # rule into the base and exponent columns. An inverted entry solves for E (linear
    # in the balance), holds base and L = log(base) generic, and takes their duals
    # from dE by the inverse chain rule. base0, E0, L0 are resting values.
    recastOut = []
    if recast:
        lVals = lVals or {}
        for rc in recast:
            E, L, base, exp = rc['E'], rc['L'], rc['base'], rc['exp']
            inverted = bool(rc.get('inverted'))
            E0 = valBy[E] if E in valBy else int(paramVals.get(E, 0)) % p
            L0 = int(lVals.get(L, 0)) % p
            base0 = valBy[base] if base in valBy else int(paramVals.get(base, 0)) % p
            expv = int(paramVals.get(exp, 0)) % p
            if not inverted:
                binv = pow(base0, p - 2, p) if base0 % p else 0
                dEbase = expv * E0 % p * binv % p
                dEexp = E0 * L0 % p
                jE = JtBy.get(E, [0] * nS)
                if base in sIdx:
                    bc = sIdx[base]
                    for i in range(nS):
                        Jxeff[i][bc] = (Jxeff[i][bc] + jE[i] * dEbase) % p
                elif base in JtBy:
                    JtBy[base] = [(JtBy[base][i] + jE[i] * dEbase) % p for i in range(nS)]
                if exp in JtBy:
                    JtBy[exp] = [(JtBy[exp][i] + jE[i] * dEexp) % p for i in range(nS)]
            recastOut.append({'E': E, 'L': L, 'base': base, 'exp': exp,
                              'inverted': inverted, 'E0': E0, 'L0': L0,
                              'base0': base0, 'expv': expv})

    # solve Jxeff * dx = -Jt for every parameter column at once over GF(p)
    params_l = [str(th) for th in paramSyms]
    B = [[(-JtBy[params_l[j]][i]) % p for j in range(len(params_l))]
         for i in range(nS)]
    X = _solve_mod(Jxeff, B, p)
    if X is None:
        return {'ok': False, 'why': 'singular jacobian mod p'}
    dxBy = {params_l[j]: [X[i][j] % p for i in range(nS)]
            for j in range(len(params_l))}

    # assemble in the full state order; held states stay at 0 with no duals
    xstar = [valBy.get(nm, 0) for nm in stateNames]
    dx = {}
    for nm in paramNames:
        col = dxBy.get(nm, [0] * nS)
        dx[nm] = [col[sIdx[s]] if s in sIdx else 0 for s in stateNames]

    # parameter-duals for the recast icSeed rows. Each entry seeds its generic
    # partner (gen) and L: a normal entry derives them from dx[base], an inverted
    # one from dx[E]. gen is E for a normal entry and base for an inverted one.
    for rc in recastOut:
        E, L, base, exp = rc['E'], rc['L'], rc['base'], rc['exp']
        E0, L0, base0, expv = rc['E0'], rc['L0'], rc['base0'], rc['expv']
        gD = {}
        lD = {}
        if not rc['inverted']:
            bi = sIdx.get(base)
            binv = pow(base0, p - 2, p) if base0 % p else 0
            dEbase = expv * E0 % p * binv % p
            dEexp = E0 * L0 % p
            for nm in paramNames:
                dxb = dxBy[nm][bi] if (bi is not None and nm in dxBy) else 0
                gD[nm] = (dEbase * dxb + (dEexp if nm == exp else 0)) % p
                lD[nm] = binv * dxb % p
            rc['gen'], rc['gen0'] = E, E0
        else:
            ei = sIdx.get(E)
            xe = pow(expv, p - 2, p) if expv % p else 0
            einv = pow(E0, p - 2, p) if E0 % p else 0
            for nm in paramNames:
                dE = dxBy[nm][ei] if (ei is not None and nm in dxBy) else 0
                gD[nm] = (base0 * xe % p * einv % p * dE
                          - (base0 * L0 % p * xe if nm == exp else 0)) % p
                lD[nm] = (xe * einv % p * dE - (L0 * xe if nm == exp else 0)) % p
            rc['gen'], rc['gen0'] = base, base0
        rc['genDual'] = gD
        rc['lDual'] = lD
        for k in ('E', 'base', 'exp', 'E0', 'base0', 'expv'):
            rc.pop(k, None)

    # compose the t0 events onto the seed (value mod p and its parameter-duals)
    events = list(t0events) if t0events else []
    if events:
        idxOfState = {nm: i for i, nm in enumerate(stateNames)}
        for evt in _compile_t0events(events, list(paramNames)):
            X = evt['var']
            if X not in idxOfState:
                continue
            i = idxOfState[X]
            v0 = _eval_terms(evt['valT'], paramvals, p)
            vdu = {nm: (_eval_terms(evt['duals'][nm], paramvals, p)
                        if evt['duals'][nm] is not None else 0)
                   for nm in paramNames}
            meth = evt['method']
            if meth == 'replace':
                xstar[i] = v0
                for nm in paramNames:
                    dx[nm][i] = vdu[nm]
            elif meth == 'add':
                xstar[i] = (xstar[i] + v0) % p
                for nm in paramNames:
                    dx[nm][i] = (dx[nm][i] + vdu[nm]) % p
            elif meth == 'multiply':
                old = xstar[i]
                oldd = {nm: dx[nm][i] for nm in paramNames}
                xstar[i] = old * v0 % p
                for nm in paramNames:
                    dx[nm][i] = (old * vdu[nm] + v0 * oldd[nm]) % p
    return {'ok': True, 'xstar': xstar, 'dx': dx, 'stateNames': list(stateNames),
            'recast': recastOut,
            'dfJx': dfJx, 'dfJt': dfJt, 'dfStateCols': dfStateCols,
            'dfParamCols': dfParamCols, 'valBy': dict(valBy)}


# ---- forward steady-state solve (choose the resting states, solve for rates) ----------
# Inverts the backward solve: choose resting states and free parameters, solve f = 0
# for a turnover subset of rates. Mass-action rates enter f linearly, so this solves
# at every prime, and a direction depending on x* (e.g. a Hill term C3^n)
# reconstructs as a rational in (theta, states) on a slice shared across primes.

def _forward_rate_pick(rhsByName, solveStates, paramNames, forcings, keepFree=None):
    """Match each state to a rate entering its balance linearly, preferring a negative
    turnover term (rate times a monomial containing the state) found in few balances.
    `keepFree` (the direction's support) is never solved for. Most-constrained states go
    first. Returns one rate name per solve state, or None if no matching exists."""
    keepFree = set(keepFree or [])
    # never solve for a recast coordinate (E, L): E multiplies a turnover rate, which
    # would make f bilinear in the rates
    paramset = {spy.Symbol(pn) for pn in paramNames
                if pn not in keepFree and not pn.startswith('_E_') and not pn.startswith('_L_')}
    forc = {spy.Symbol(nm) for nm in forcings}
    cand, appear = {}, {}
    for s in solveStates:
        f = spy.sympify(rhsByName[s]); sSym = spy.Symbol(s); lst = []
        for r in sorted(f.free_symbols & paramset, key=str):
            if r in forc:
                continue
            try:
                pv = spy.Poly(f, r)
            except spy.PolynomialError:
                continue
            if pv.degree() != 1 or r in pv.nth(1).free_symbols:
                continue
            coef = pv.nth(1)
            turnover = sSym in coef.free_symbols
            neg = coef.could_extract_minus_sign()
            lst.append((0 if (turnover and neg) else 1 if turnover else 2, str(r)))
            appear[str(r)] = appear.get(str(r), 0) + 1
        cand[s] = lst
    order = sorted(solveStates, key=lambda s: len(cand[s]))
    used, assign = set(), {}
    for s in order:
        opts = sorted((sc, appear.get(r, 99), r) for sc, r in cand[s] if r not in used)
        if not opts:
            return _rate_matching(cand, appear, order)
        assign[s] = opts[0][2]; used.add(opts[0][2])
    return [assign[s] for s in solveStates]


def _rate_matching(cand, appear, order):
    """A complete matching of states to candidate rates by augmenting paths, each state
    trying its rates in preference order; None if there is none."""
    prefs = {s: [r for _, _, r in sorted((sc, appear.get(r, 99), r) for sc, r in cand[s])]
             for s in order}
    owner = {}
    def augment(s, seen):
        for r in prefs[s]:
            if r in seen:
                continue
            seen.add(r)
            if r not in owner or augment(owner[r], seen):
                owner[r] = s
                return True
        return False
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 10 * len(order) + 100))
    for s in order:
        if not augment(s, set()):
            return None
    byState = {s: r for r, s in owner.items()}
    return [byState[s] for s in cand]


def _forward_plan(rhsByName, solveStates, paramNames, forcings, keepFree=None):
    """Forward plan when balances share rates: a balance without a rate of its own is
    solved for its own state, before the rates. Returns (stateSolve, solveRates), the
    state-solved balances in dependency order and the rates, or None. A balance whose
    rate coefficients depend on the others' at a random point is state-solved too."""
    pre = []
    for _ in range(len(solveStates)):
        r = _forward_plan_once(rhsByName, solveStates, paramNames, forcings, keepFree, pre)
        if r is None or r[0] != 'dependent':
            return r
        new = [u for u in r[1] if u not in pre]
        if not new:
            return None
        pre = pre + new[:1]
    return None


def _forward_plan_once(rhsByName, solveStates, paramNames, forcings, keepFree=None,
                       preBs=()):
    keepFree = set(keepFree or [])
    paramset = {spy.Symbol(pn) for pn in paramNames
                if pn not in keepFree and not pn.startswith('_E_') and not pn.startswith('_L_')}
    forc = {spy.Symbol(nm) for nm in forcings}
    polys = {s: spy.sympify(rhsByName[s]) for s in solveStates}
    params = {s: polys[s].free_symbols & paramset for s in solveStates}
    cand = {}
    for s in solveStates:
        lst = []
        for r in sorted(params[s] - forc, key=str):
            try:
                pv = spy.Poly(polys[s], r)
            except spy.PolynomialError:
                continue
            if pv.degree() == 1 and r not in pv.nth(1).free_symbols:
                lst.append(str(r))
        cand[s] = lst
    stateSet = set(solveStates)
    # parameters sharing a monomial; no two of them are solve rates
    pnames = {str(r) for r in paramset}
    partners = {}
    for s in solveStates:
        for mono in spy.Add.make_args(spy.expand(polys[s])):
            ps = [str(x) for x in mono.free_symbols if str(x) in pnames]
            for a in ps:
                partners.setdefault(a, set()).update(x for x in ps if x != a)
    banned = set()
    bs = list(preBs)                          # balances solved for their own state
    for u in bs:
        try:
            pu = spy.Poly(polys[u], spy.Symbol(u))
        except spy.PolynomialError:
            return None
        if pu.degree() != 1 or spy.Symbol(u) in pu.nth(1).free_symbols:
            return None
    while True:
        forbidden = {str(r) for u in bs for r in params[u]} | banned
        rest = [s for s in solveStates if s not in bs]
        match = {}                            # rate -> balance, by augmenting paths
        def augment(s, seen):
            for r in cand[s]:
                if r in forbidden or r in seen:
                    continue
                seen.add(r)
                if r not in match or augment(match[r], seen):
                    match[r] = s
                    return True
            return False
        for s in sorted(rest, key=lambda s: len(cand[s])):
            augment(s, set())
        # ban the solve rate with the most partners among the chosen, then rematch
        chosen = set(match)
        clash = sorted((r for r in chosen if partners.get(r, set()) & chosen),
                       key=lambda r: (-len(partners.get(r, ())), r))
        if clash:
            banned.add(clash[0])
            continue
        matched = set(match.values())
        unmatched = [s for s in rest if s not in matched]
        if not unmatched:
            break
        for u in unmatched:
            try:
                pu = spy.Poly(polys[u], spy.Symbol(u))
            except spy.PolynomialError:
                return None
            if pu.degree() != 1 or spy.Symbol(u) in pu.nth(1).free_symbols:
                return None
        bs.extend(unmatched)
    # dependency order: a state-solved balance reads the other states it contains
    deps = {u: {str(x) for x in polys[u].free_symbols if str(x) in stateSet} & set(bs) - {u}
            for u in bs}
    order = []
    while len(order) < len(bs):
        ready = [u for u in bs if u not in order and deps[u] <= set(order)]
        if not ready:
            return None
        order.extend(sorted(ready))
    # rates as a column basis of the coefficient matrix at a random point, turnover
    # rates first, so the rate system is nonsingular
    rest = [s for s in solveStates if s not in bs]
    forbidden = {str(r) for u in bs for r in params[u]}
    Q = 2147483629
    rng = random.Random(len(rest))
    syms = set().union(*(polys[s].free_symbols for s in rest)) if rest else set()
    pt = {x: rng.randrange(2, Q - 1) for x in syms}
    def score(r):
        best = 2
        for s in rest:
            if r not in cand[s]:
                continue
            coef = spy.Poly(polys[s], spy.Symbol(r)).nth(1)
            if spy.Symbol(s) in coef.free_symbols:
                best = min(best, 0 if coef.could_extract_minus_sign() else 1)
        return best
    cols = sorted({r for s in rest for r in cand[s]} - forbidden,
                  key=lambda r: (score(r), sum(r in cand[s] for s in rest), r))
    colVec = {}
    def vecOf(r):
        if r not in colVec:
            colVec[r] = [int(spy.Poly(polys[s], spy.Symbol(r)).nth(1).xreplace(pt)) % Q
                         if r in cand[s] else 0 for s in rest]
        return colVec[r]
    def pickBasis(colOrder):
        basis, chosen = [], []                # echelon rows (pivot, vector) mod Q
        for r in colOrder:
            if partners.get(r, set()) & set(chosen):
                continue
            v = list(vecOf(r))
            for piv, b in basis:
                if v[piv]:
                    fct = v[piv]
                    v = [(x - fct * y) % Q for x, y in zip(v, b)]
            piv = next((i for i, x in enumerate(v) if x), None)
            if piv is None:
                continue
            inv = pow(v[piv], Q - 2, Q)
            basis.append((piv, [x * inv % Q for x in v]))
            chosen.append(r)
            if len(chosen) == len(rest):
                break
        return chosen
    chosen = pickBasis(cols)
    # the partner rule can block the greedy order: shuffles within each preference class
    keyOf = {r: score(r) for r in cols}
    for att in range(40):
        if len(chosen) == len(rest):
            break
        rs = random.Random(att)
        shuffled = sorted(cols, key=lambda r: (keyOf[r], rs.random()))
        cand2 = pickBasis(shuffled)
        if len(cand2) > len(chosen):
            chosen = cand2
    if len(chosen) < len(rest):
        # the balances whose coefficient rows depend on the earlier ones
        rows, dep = [], []
        for s in rest:
            v = [int(spy.Poly(polys[s], spy.Symbol(r)).nth(1).xreplace(pt)) % Q
                 if r in cand[s] else 0 for r in cols]
            for piv, b in rows:
                if v[piv]:
                    fct = v[piv]
                    v = [(x - fct * y) % Q for x, y in zip(v, b)]
            piv = next((i for i, x in enumerate(v) if x), None)
            if piv is None:
                dep.append(s)
                continue
            inv = pow(v[piv], Q - 2, Q)
            rows.append((piv, [x * inv % Q for x in v]))
        return ('dependent', dep)
    return order, chosen


_forwardCache = {}


def _forward_compile(model, stateNames, paramNames, forcings, solveRates, stateSolve=()):
    """Cached prime-independent compile of the forward solve: each solve rate's coefficient
    and the rate-free constant per balance as term lists over rgens = (non-solve params) +
    solve states, plus the steady-state Jacobian term lists. A balance in `stateSolve` is
    compiled as its own state's coefficient and remainder instead. Returns {'bad': reason}
    if the balances are not linear in the chosen rates."""
    stateSolve = tuple(stateSolve)
    key = (tuple(model), tuple(stateNames), tuple(paramNames), tuple(sorted(forcings)),
           tuple(solveRates), stateSolve)
    c = _forwardCache.get(key)
    if c is not None:
        return c
    (paramSyms, solveStates, polys, Jx, Jt, gens0, JxTerms, JtTerms,
     polyBi, genLin, linTerms, _pp) = _ss_compile(model, stateNames, paramNames, forcings)
    solveSet = set(solveRates); rSyms = [spy.Symbol(r) for r in solveRates]
    rgens = [th for th in paramSyms if str(th) not in solveSet] + list(solveStates)
    rgenset = set(rgens)
    idxOf = {str(st): i for i, st in enumerate(solveStates)}
    stateT = []
    for u in stateSolve:
        pu = spy.Poly(spy.sympify(polys[idxOf[u]]), spy.Symbol(u))
        c1, c0 = pu.nth(1).as_expr(), pu.nth(0).as_expr()
        if (c1.free_symbols | c0.free_symbols) - rgenset:
            return {'bad': 'state-solved balance %s holds a solve rate' % u}
        stateT.append((u, _poly_terms(c1, rgens), _poly_terms(c0, rgens)))
    rateRows = [i for i in range(len(solveStates)) if str(solveStates[i]) not in set(stateSolve)]
    coefT, constT = [], []
    for i in rateRows:
        fexp = spy.sympify(polys[i]); row = []
        for r in rSyms:
            dexp = spy.diff(fexp, r)
            if dexp.free_symbols - rgenset:
                return {'bad': 'not linear in the chosen rates at %s' % str(solveStates[i])}
            row.append(_poly_terms(dexp, rgens))
        cexp = fexp.subs({r: spy.Integer(0) for r in rSyms})
        if cexp.free_symbols - rgenset:
            return {'bad': 'unsubstituted symbol in balance %s' % str(solveStates[i])}
        coefT.append(row); constT.append(_poly_terms(cexp, rgens))
    c = {'paramNames': [str(s) for s in paramSyms], 'solveStates': [str(s) for s in solveStates],
         'rgens': [str(g) for g in rgens], 'coefT': coefT, 'constT': constT, 'stateT': stateT,
         'JxTerms': JxTerms, 'JtTerms': JtTerms, 'gens': [str(g) for g in gens0]}
    _forwardCache[key] = c
    return c


def solveForwardModular(model, stateNames, paramNames, stateVals, paramVals, prime,
                        forcings=None, solveRates=None, keepFree=None, backend='sympy',
                        stateSolve=None):
    """Solve f = 0 over GF(prime) for a turnover subset of rates, given chosen resting
    `stateVals` and `paramVals` (residues; forcings held at 0). `solveRates` (one per
    non-forcing state outside `stateSolve`) is auto-picked avoiding `keepFree` if None;
    the balances in `stateSolve` are solved for their own state first. Returns the
    joint-mode payload of the backward solve plus {'rates', 'solveRates', 'stateSolve'},
    or {'ok': False, 'why'}."""
    _select_backend(backend)
    p = int(prime)
    model = _as_list(model); forcings = set(_as_list(forcings))
    stateNames = _as_list(stateNames); paramNames = _as_list(paramNames)
    svd = {str(k): int(v) % p for k, v in dict(stateVals).items()}
    pvd = {str(k): int(v) % p for k, v in dict(paramVals).items()}
    stateSolve = [] if stateSolve is None else _as_list(stateSolve)
    if solveRates is None:
        (_ps, _ss, _polys, *_r) = _ss_compile(model, stateNames, paramNames, forcings)
        rhsByName = {str(_ss[i]): _polys[i] for i in range(len(_ss))}
        ssNames = [str(s) for s in _ss]
        solveRates = _forward_rate_pick(rhsByName, ssNames, paramNames, forcings, keepFree=keepFree)
        stateSolve = []
        if solveRates is None:
            plan = _forward_plan(rhsByName, ssNames, paramNames, forcings, keepFree=keepFree)
            if plan is None:
                return {'ok': False, 'why': 'no complete forward rate matching'}
            stateSolve, solveRates = plan
    solveRates = _as_list(solveRates)
    c = _forward_compile(model, stateNames, paramNames, forcings, solveRates, stateSolve)
    if 'bad' in c:
        return {'ok': False, 'why': c['bad']}
    solveStates = c['solveStates']; nS = len(solveStates); nR = len(c['coefT'])
    if len(solveRates) != nR:
        return {'ok': False, 'why': 'rate/state count mismatch (%d rates, %d balances)'
                % (len(solveRates), nR)}
    def cvfree(nm):
        return svd[nm] if nm in svd else (pvd[nm] if nm in pvd else 0)
    rgIdx = {nm: k for k, nm in enumerate(c['rgens'])}
    rgv = [cvfree(nm) for nm in c['rgens']]
    # state-solved balances first, each x = -c0/c1 at the values known so far
    for u, c1T, c0T in c['stateT']:
        c1 = _eval_terms(c1T, rgv, p)
        if c1 == 0:
            return {'ok': False, 'why': 'vanishing state pivot mod p at %s' % u}
        xu = (-_eval_terms(c0T, rgv, p)) * pow(c1, p - 2, p) % p
        if xu == 0:
            return {'ok': False, 'why': 'zero resting state %s' % u}
        svd[u] = xu
        rgv[rgIdx[u]] = xu
    A = [[_eval_terms(c['coefT'][i][j], rgv, p) for j in range(nR)] for i in range(nR)]
    b = [[(-_eval_terms(c['constT'][i], rgv, p)) % p] for i in range(nR)]
    X = _solve_mod(A, b, p)
    if X is None:
        return {'ok': False, 'why': 'singular forward system mod p'}
    ratesDict = {solveRates[j]: int(X[j][0]) % p for j in range(len(solveRates))}
    def cvfull(nm):
        if nm in ratesDict: return ratesDict[nm]
        return svd[nm] if nm in svd else (pvd[nm] if nm in pvd else 0)
    ptvals = [cvfull(nm) for nm in c['gens']]
    dfJx = [[_eval_terms(c['JxTerms'][i][j], ptvals, p) for j in range(nS)] for i in range(nS)]
    dfJt = {th: [_eval_terms(c['JtTerms'][th][i], ptvals, p) for i in range(nS)]
            for th in c['paramNames']}
    return {'ok': True, 'rates': ratesDict, 'solveRates': list(solveRates),
            'stateSolve': list(stateSolve),
            'solveStates': list(solveStates), 'valBy': {s: cvfull(s) for s in solveStates},
            'dfJx': dfJx, 'dfJt': dfJt, 'dfStateCols': list(solveStates),
            'dfParamCols': list(c['paramNames'])}


_HYPERBOLIC = (spy.sinh, spy.cosh, spy.tanh, spy.coth, spy.sech, spy.csch)


def _hyp_to_exp(e):
    """Hyperbolic functions written in exponentials."""
    e = spy.sympify(e)
    if not e.has(*_HYPERBOLIC):
        return e
    return e.replace(lambda x: isinstance(x, _HYPERBOLIC), lambda x: x.rewrite(spy.exp))


def _has_exp(e):
    """True if e contains an exponential or the number E."""
    e = spy.sympify(e)
    return e.has(spy.E) or any(_exp_atom(a) is not None
                               for a in e.atoms(spy.exp, spy.Pow))


def _strip_log_obs(g):
    """An observable kappa*(sum a_i log(h_i)) + c, with a_i rational and kappa 1 or
    1/log(b) for a numeric b, is a strictly monotone function of exp(g/kappa) =
    prod h_i^a_i * exp(c/kappa) and carries the same information; that expression is
    returned. Anything else comes back unchanged."""
    g = spy.sympify(g)
    logs = [a for a in g.atoms(spy.log) if not a.args[0].is_number]
    if not logs:
        return g
    consts = [spy.Integer(1)] + [a for a in g.atoms(spy.log) if a.args[0].is_number]
    dums = [spy.Dummy() for _ in logs]
    gd = g.xreplace(dict(zip(logs, dums)))
    for kap in consts:
        g2 = spy.expand(gd * kap)
        coeffs = [g2.coeff(d) for d in dums]
        rest = spy.expand(g2 - sum(c * d for c, d in zip(coeffs, dums)))
        if rest.has(*dums):
            continue
        if not all(c.is_Rational for c in coeffs):
            continue
        if any(not a.args[0].is_number for a in rest.atoms(spy.log)):
            continue
        out = spy.Integer(1)
        for c, L in zip(coeffs, logs):
            out = out * L.args[0] ** c
        for t in spy.Add.make_args(rest):
            nl = [a for a in t.atoms(spy.log) if a.args[0].is_number]
            out = out * (nl[0].args[0] ** spy.cancel(t / nl[0]) if len(nl) == 1
                         else spy.exp(t))
        if not out.atoms(spy.log):
            return out
    return g


def _exp_atom(e):
    """(base, exponent) of an exponential atom with a numeric base, else None."""
    if isinstance(e, spy.exp):
        return spy.E, e.args[0]
    if e.is_Pow and e.base.is_number and e.base.is_positive and not e.exp.is_number:
        return e.base, e.exp
    return None


def _detect_log_params(exprs, states):
    """Parameters that occur only in exponents of numeric-base powers, base^(c*theta +
    ...) with c rational and one base per parameter. Returns {theta: base}. A constant
    term in the exponent must leave a rational factor base^k."""
    baseOf, bad = {}, set()
    for e in exprs:
        for at in spy.sympify(e).atoms(spy.exp, spy.Pow):
            be = _exp_atom(at)
            if be is None:
                continue
            base, ex = be
            terms = spy.expand(ex).as_coefficients_dict()
            k = terms.pop(spy.Integer(1), spy.Integer(0))
            ok = k == 0 or (base != spy.E and base.is_Rational and k.is_Integer)
            for t, cf in terms.items():
                if not (t.is_Symbol and cf.is_Rational and t not in states) or not ok:
                    bad |= ex.free_symbols
                    break
                if baseOf.get(t, base) != base:
                    bad.add(t)
                baseOf[t] = base
    return {t: b for t, b in baseOf.items() if t not in bad}


def _log_param_sub(lp):
    """Rewrite base^(sum c_i theta_i + k) as base^k * prod X_i^c_i, X_i = base^theta_i."""
    def repl(at):
        be = _exp_atom(at)
        if be is None:
            return at
        base, ex = be
        terms = spy.expand(ex).as_coefficients_dict()
        k = terms.pop(spy.Integer(1), spy.Integer(0))
        if not terms or any(t not in lp or lp[t][1] != base for t in terms):
            return at
        out = base ** k if k != 0 else spy.Integer(1)
        for t, cf in terms.items():
            out = out * lp[t][0] ** cf
        return out

    def sub(e):
        e = spy.sympify(e)
        if not lp:
            return e
        return e.replace(lambda x: _exp_atom(x) is not None, repl)
    return sub


def _log_params_only_in_atoms(baseOf, exprs):
    """Drop every candidate that still occurs after the rewrite: it enters somewhere
    other than an exponent and cannot be traded for base^theta."""
    baseOf = dict(baseOf)
    while baseOf:
        sub = _log_param_sub({t: (spy.Dummy(), b) for t, b in baseOf.items()})
        left = set()
        for e in exprs:
            left |= sub(e).free_symbols
        drop = [t for t in baseOf if t in left]
        if not drop:
            break
        for t in drop:
            del baseOf[t]
    return baseOf


def signChart(f, gs, positive):
    """abs() and sign() resolved with the declared signs (abs(v) = v for positive
    v). `f` maps state -> rhs, `gs` is a list of observable dicts. Returns {'f',
    'g'} as strings, or {'why': ...} for max, min, a step (not analytic) or an
    abs() of undecided sign."""
    fl = dict(f)
    gl = [dict(g) for g in gs]
    lines = list(fl.values()) + [v for g in gl for v in g.values()] + [str(k) for k in fl]
    local, parse = _make_local_parse(lines)
    local = dict(local)
    local.update({'abs': spy.Abs, 'sign': spy.sign})
    parse = _make_parse(local)
    pos = None if positive is True else {str(x) for x in (positive or [])}
    names = {str(sym) for sym in local.values() if isinstance(sym, spy.Symbol)}
    toPos = {spy.Symbol(n): spy.Symbol(n, positive=True)
             for n in names if pos is None or n in pos}
    back = {v: k for k, v in toPos.items()}

    def one(k, rhs):
        if any(t in rhs for t in ('max(', 'min(', 'pmax(', 'pmin(', 'ifelse(',
                                  'Heaviside(')):
            raise ValueError('%s: max, min and steps are not analytic; analyse each '
                             'regime on its own, e.g. through `conditions`' % k)
        e = parse(rhs).xreplace(toPos)
        if e.atoms(spy.Abs, spy.sign):
            raise ValueError('%s: the sign inside abs() or sign() is not decided by '
                             'the coordinates declared `positive`' % k)
        return str(e.xreplace(back))
    try:
        Fn = {k: one(k, str(v)) for k, v in fl.items()}
        Gn = [{k: one(k, str(v)) for k, v in g.items()} for g in gl]
    except ValueError as err:
        return {'why': str(err)}
    return {'f': Fn, 'g': Gn}


def _log_numbers(e, consts):
    """Rewrite log of a positive rational as sum e_p*lognum_p over its prime factors,
    collecting the constants lognum_p in `consts`. Exact: log 4 = 2 log 2 is kept, and
    logs of distinct primes are independent (Schanuel)."""
    rep = {}
    for a in e.atoms(spy.log):
        n = a.args[0]
        if not n.is_number:
            continue
        n = spy.nsimplify(n)
        if not (n.is_Rational and n > 0):
            raise ValueError('log(%s) of a number that is not a positive rational' % n)
        out = spy.Integer(0)
        for q, sgn in ((n.p, 1), (n.q, -1)):
            for pr, ex in spy.factorint(q).items():
                sym = spy.Symbol('lognum_%d' % pr)
                consts.add(str(sym))
                out += sgn * ex * sym
        rep[a] = out
    return e.xreplace(rep)


def logArgChart(f, gs, positive, taken, ics=None):
    """Log chart for positive arguments of log() and fractional powers: log(a + b*v),
    a and b free of states and v, a + b*v positive, becomes L with v = (exp(L) - a)/b;
    log(v) and v^(p/q) are a = 0, b = 1. A state v becomes L with rhs b*f_v*exp(-L)
    and initial value log(a + b*v0). `f` maps state -> rhs, `gs` and `ics` are lists
    of dicts, `positive` True or a list of names. Returns None when nothing qualifies,
    else {'f', 'g', 'ic', 'consts', 'map': [{'v', 'L', 'a', 'b'}]}, or {'why': ...}."""
    fl = dict(f)
    gl = [dict(g) for g in gs]
    il = [dict(ic) for ic in (ics or [])]
    lines = (list(fl.values()) + [v for g in gl for v in g.values()] +
             [v for ic in il for v in ic.values()] +
             [str(k) for k in fl] + [str(t) for t in taken])
    local, parse = _make_local_parse(lines)
    F = {k: parse(str(v)) for k, v in fl.items()}
    G = [{k: parse(str(v)) for k, v in g.items()} for g in gl]
    I = [{k: parse(str(v)) for k, v in ic.items()} for ic in il]
    pos = None if positive is True else {str(x) for x in (positive or [])}
    states = {spy.Symbol(k) for k in fl}
    one, zero = spy.Integer(1), spy.Integer(0)

    def isPos(sym):
        return pos is None or str(sym) in pos

    def positiveOn(e):
        # sign certificate on the declared domain
        sub = {x: spy.Symbol(str(x), positive=True) for x in e.free_symbols if isPos(x)}
        return bool(e.xreplace(sub).is_positive)

    chart = {}                               # v -> (a, b)

    def claim(v, a, b):
        if v in chart and (spy.simplify(chart[v][0] - a) != 0 or
                           spy.simplify(chart[v][1] - b) != 0):
            raise ValueError('%s enters two different logarithms, which no single '
                             'chart makes rational' % v)
        chart[v] = (a, b)

    def scan(e):
        e = e.xreplace({at: spy.log(spy.factor(spy.together(at.args[0])))
                        for at in e.atoms(spy.log)})
        e = spy.expand_log(e, force=True)
        for at in e.atoms(spy.log):
            arg = at.args[0]
            if arg.is_number:
                continue
            if arg.is_Symbol:
                if isPos(arg):
                    claim(arg, zero, one)
                continue
            # affine in one symbol: a state first, then a parameter
            for v in sorted(arg.free_symbols, key=lambda x: (x not in states, str(x))):
                b = spy.diff(arg, v)
                if b.free_symbols & (states | {v}):
                    continue
                a = spy.expand(arg - b * v)
                if a.free_symbols & states or not positiveOn(arg):
                    continue
                claim(v, a, b)
                break
        for pw in e.atoms(spy.Pow):
            if pw.base.is_Symbol and pw.exp.is_Rational and not pw.exp.is_Integer:
                if isPos(pw.base):
                    claim(pw.base, zero, one)

    try:
        for e in F.values():
            scan(e)
        for g in G:
            for e in g.values():
                scan(_strip_log_obs(e))
        # a free-exponent base with a given initial value has no leaf for the recast
        given = {spy.Symbol(k) for ic in I for k in ic}
        for e in list(F.values()) + [e for g in G for e in g.values()]:
            for pw in e.atoms(spy.Pow):
                if pw.base in given and not pw.exp.is_number and isPos(pw.base):
                    claim(pw.base, zero, one)
        # an initial value of a charted state enters as the logarithm of its argument
        while True:
            n0 = len(chart)
            for ic in I:
                for k, e in ic.items():
                    ks = spy.Symbol(k)
                    if ks in chart:
                        a, b = chart[ks]
                        scan(spy.log(a + b * e))
            if len(chart) == n0:
                break
    except ValueError as err:
        return {'why': str(err)}
    if not chart:
        return None

    used = {str(t) for t in taken} | {str(k) for k in fl}
    L = {}
    for v in sorted(chart, key=str):
        nm = 'log_%s' % v
        while nm in used:
            nm += '_'
        used.add(nm)
        L[v] = spy.Symbol(nm, real=True)
    sub0 = {v: (spy.exp(L[v]) - chart[v][0]) / chart[v][1] for v in chart}
    # a and b may hold other charted symbols
    sub = dict(sub0)
    for _ in chart:
        sub = {v: e.xreplace(sub) for v, e in sub0.items()}
    consts = set()

    def tidy(e):
        e = e.xreplace(sub)
        e = e.xreplace({at: spy.log(spy.cancel(at.args[0])) for at in e.atoms(spy.log)})
        return _log_numbers(spy.powsimp(spy.expand_log(e, force=True)), consts)
    Fn = {}
    for k, e in F.items():
        ks = spy.Symbol(k)
        if ks in chart:
            Fn[str(L[ks])] = str(spy.powsimp(spy.expand(
                tidy(chart[ks][1] * e) * spy.exp(-L[ks]))))
        else:
            Fn[k] = str(tidy(e))
    Gn = [{k: str(tidy(e)) for k, e in g.items()} for g in G]
    In = []
    for ic in I:
        d = {}
        for k, e in ic.items():
            ks = spy.Symbol(k)
            if ks in chart:
                a, b = chart[ks]
                le = tidy(spy.log(a + b * e))
                if any(not at.args[0].is_number for at in le.atoms(spy.log)):
                    return {'why': 'the initial value %s = %s does not split into '
                            'logarithms of positive coordinates' % (k, e)}
                try:
                    d[str(L[ks])] = str(_log_numbers(le, consts))
                except ValueError as err:
                    return {'why': str(err)}
            else:
                d[k] = str(tidy(e))
        In.append(d)
    return {'f': Fn, 'g': Gn, 'ic': In, 'consts': sorted(consts),
            'map': [{'v': str(v), 'L': str(L[v]), 'a': str(chart[v][0]),
                     'b': str(chart[v][1])} for v in sorted(chart, key=str)]}


def logArgEvent(var, value, method, maps):
    """One event in the log chart: v -> (exp(L) - a)/b in its value, and an event
    on a charted state acts on L: a replacement by log(a + b*value), and for
    a = 0 a multiplication as the addition of log(value). Returns {'var',
    'value', 'method', 'consts'} as strings or {'why': ...}."""
    maps = list(maps)
    lines = [str(value)] + [str(m[k]) for m in maps for k in ('v', 'L', 'a', 'b')]
    local, parse = _make_local_parse(lines)
    ent = {str(m['v']): (spy.Symbol(str(m['L']), real=True), parse(str(m['a'])),
                         parse(str(m['b']))) for m in maps}
    sub = {parse(v): (spy.exp(Lv) - a) / b for v, (Lv, a, b) in ent.items()}
    e = parse(str(value)).subs(sub)
    if str(var) not in ent:
        return {'var': str(var), 'value': str(e), 'method': str(method), 'consts': []}
    Lv, a, b = ent[str(var)]
    if method == 'add' and e == 0:
        return {'var': str(Lv), 'value': '0', 'method': 'add', 'consts': []}
    if method == 'replace':
        le = spy.log(a + b * e)
        kind = 'replace'
    elif method == 'multiply' and a == 0:
        le = spy.log(e)
        kind = 'add'
    else:
        return {'why': 'a %s event on %s, which is analysed as %s, is not supported'
                % (method, var, 'log(%s)' % (a + b * parse(str(var))))}
    le = spy.powsimp(spy.expand_log(le, force=True))
    if any(not at.args[0].is_number for at in le.atoms(spy.log)):
        return {'why': 'the event value %s for %s does not split into logarithms of '
                'positive coordinates' % (value, var)}
    consts = set()
    try:
        le = _log_numbers(le, consts)
    except ValueError as err:
        return {'why': str(err)}
    return {'var': str(Lv), 'value': str(le), 'method': kind, 'consts': sorted(consts)}


def logArgBackVec(vec, maps, transform=False):
    """A generator (or, with `transform`, a finite transformation) of the log chart
    back in the user's symbols. L -> log(a + b*v) everywhere; for a component
    eta_v = ((a + b*v)*eta_L - eta(a) - v*eta(b))/b with eta(a) = grad(a) . eta,
    for a transformation v = (exp(T_L) - a(T))/b(T)."""
    maps = list(maps)
    lines = [str(x) for x in vec.values()] + [str(k) for k in vec] + \
        [str(m[k]) for m in maps for k in ('v', 'L', 'a', 'b')]
    local, parse = _make_local_parse(lines)
    V = {str(k): parse(str(x)) for k, x in vec.items()}
    ent = [(parse(str(m['v'])), str(m['L']), parse(str(m['a'])), parse(str(m['b'])))
           for m in maps]
    Lnames = {Ln for _, Ln, _, _ in ent}
    vp = {v: spy.Symbol(str(v), positive=True) for v, _, _, _ in ent}
    toL = {parse(Ln): spy.log(a + b * v) for v, Ln, a, b in ent}
    out = {k: e.subs(toL) for k, e in V.items() if k not in Lnames}
    # charts whose a, b hold another charted symbol come after it
    order = sorted(ent, key=lambda t: len((t[2].free_symbols | t[3].free_symbols) &
                                          set(vp)))
    for v, Ln, a, b in order:
        eL = V.get(Ln, spy.Integer(0)).subs(toL)
        if transform:
            T = {parse(k): e for k, e in out.items()}
            aT, bT = a.xreplace(T), b.xreplace(T)
            if eL == 0 and aT == a and bT == b:
                continue
            nv = (spy.exp(eL if Ln in V else spy.log(a + b * v)) - aT) / bT
        else:
            grad = lambda h: sum((spy.diff(h, parse(k)) * e for k, e in out.items()),
                                 spy.Integer(0))
            nv = ((a + b * v) * eL - grad(a) - v * grad(b)) / b
        out[str(v)] = nv
    res = {}
    for k, e in out.items():
        e = spy.powsimp(spy.expand_log(e.xreplace(vp), force=True))
        e = spy.factor_terms(spy.cancel(e)).xreplace({s: v for v, s in vp.items()})
        if e != 0 or k in vec:
            res[k] = str(e)
    return res


def scalingWeights(comps):
    """{name: component} -> {name: weight} when every component is a number times
    its own coordinate, else None."""
    local, parse = _make_local_parse([str(v) for v in comps.values()] +
                                     [str(k) for k in comps])
    out = {}
    for k, v in comps.items():
        w = spy.cancel(parse(str(v)) / parse(str(k)))
        if not w.is_Rational:
            return None
        out[str(k)] = str(w)
    return out


def logChart(gens, coords):
    """Generators in the chart X = b^theta for every coordinate theta that enters
    them only through b^(c*theta), or whose own component carries 1/log(b): there
    they are rational. `gens` is a list of {name: expression} (None for a support-only
    direction). Returns {'map': [...], 'gens': [...]} or None when there is no such
    coordinate or the chart does not make every component rational."""
    exprsOf = []
    names = [str(c) for c in coords]
    for g in gens:
        if g is None:
            exprsOf.append(None)
            continue
        local, parse = _make_local_parse([str(v) for v in g.values()] + names)
        exprsOf.append({str(k): spy.sympify(parse(str(v))) for k, v in g.items()})
    allE = [e for g in exprsOf if g for e in g.values()]
    coordSyms = {spy.Symbol(n) for n in names}
    baseOf = _detect_log_params(allE, set())
    baseOf = {t: b for t, b in baseOf.items() if t in coordSyms}
    # a coordinate whose own component carries 1/log(b) and that enters nowhere else
    for g in exprsOf:
        if not g:
            continue
        for k, v in g.items():
            t = spy.Symbol(k)
            nl = {a.args[0] for a in v.atoms(spy.log) if a.args[0].is_number}
            if t in baseOf or len(nl) != 1:
                continue
            if all(t not in e.free_symbols for e in allE):
                baseOf[t] = nl.pop()
    baseOf = _log_params_only_in_atoms(baseOf, allE)
    if not baseOf:
        return None
    taken = set(names)
    lp = {}
    for t, b in sorted(baseOf.items(), key=lambda kv: str(kv[0])):
        nm = 'exp_%s' % t
        while nm in taken:
            nm += '_'
        taken.add(nm)
        lp[t] = (spy.Symbol(nm), b)
    sub = _log_param_sub(lp)
    out = []
    for g in exprsOf:
        if g is None:
            out.append(None)
            continue
        h = {}
        for k, v in g.items():
            t = spy.Symbol(k)
            if t in lp:
                X, b = lp[t]
                h[str(X)] = sub(X * spy.log(b) * v)
            else:
                h[k] = sub(v)
        h = {k: spy.cancel(v) for k, v in h.items()}
        # a translation of theta is X*log(b) d/dX: the constant factor goes
        if any(v.atoms(spy.log) for v in h.values()):
            for b in {b for _, b in lp.values()}:
                h2 = {k: spy.cancel(v / spy.log(b)) for k, v in h.items()}
                if not any(v.atoms(spy.log) for v in h2.values()):
                    h = h2
                    break
        for v in h.values():
            if v.atoms(spy.log, spy.exp) or any(t in v.free_symbols for t in lp):
                return None
        out.append({k: str(v) for k, v in h.items()})
    return {'map': [{'theta': str(t), 'X': str(X), 'base': 'E' if b == spy.E else str(b)}
                    for t, (X, b) in lp.items()],
            'gens': out}


def expChart(gens, coords, positive):
    """Generators in the chart L = log(v) for every positive coordinate v that enters
    them through log(v), or that a scaling with symbolic weights moves:
    eta_L = eta_v/v, v = exp(L). `gens` as in logChart,
    `positive` True or a list of names. Returns {'map': [{'v', 'L'}], 'gens': [...]}
    or None when no coordinate qualifies or a component stays transcendental."""
    names = [str(c) for c in coords]
    pos = None if positive is True else {str(x) for x in (positive or [])}
    exprsOf = []
    for g in gens:
        if g is None:
            exprsOf.append(None)
            continue
        local, parse = _make_local_parse([str(v) for v in g.values()] + names)
        exprsOf.append({str(k): spy.expand_log(spy.sympify(parse(str(v))), force=True)
                        for k, v in g.items()})
    coordSyms = {spy.Symbol(n) for n in names}
    cand = set()
    for g in exprsOf:
        for e in (g or {}).values():
            for a in e.atoms(spy.log):
                if a.args[0] in coordSyms and (pos is None or str(a.args[0]) in pos):
                    cand.add(a.args[0])
    # a scaling with symbolic weights w translates log(z) at the rate w
    for g in exprsOf:
        if not g:
            continue
        supp = {spy.Symbol(k) for k in g}
        w = {k: spy.cancel(e / spy.Symbol(k)) for k, e in g.items()}
        if any(ww.free_symbols & supp or ww.atoms(spy.log, spy.exp) for ww in w.values()):
            continue
        if all(ww.is_number for ww in w.values()):
            continue
        if all(pos is None or str(z) in pos for z in supp):
            cand |= supp
    if not cand:
        return None
    taken = set(names)
    L = {}
    for v in sorted(cand, key=str):
        nm = 'log_%s' % v
        while nm in taken:
            nm += '_'
        taken.add(nm)
        L[v] = spy.Symbol(nm, real=True)
    sub = {v: spy.exp(L[v]) for v in cand}
    out = []
    for g in exprsOf:
        if g is None:
            out.append(None)
            continue
        h = {}
        for k, e in g.items():
            ks = spy.Symbol(k)
            e2 = spy.expand_log(e.subs(sub), force=True)
            if ks in cand:
                h[str(L[ks])] = spy.cancel(spy.powsimp(e2 * spy.exp(-L[ks])))
            else:
                h[k] = spy.cancel(spy.powsimp(e2))
        for v in h.values():
            if v.atoms(spy.log, spy.exp):
                return None
        out.append({k: str(v) for k, v in h.items()})
    return {'map': [{'v': str(v), 'L': str(L[v])} for v in sorted(cand, key=str)],
            'gens': out}


def expChartBack(expr, vNames, lNames, solveFor=None):
    """An expression of the chart L = log(v) back in v; with `solveFor` naming an L
    the expression is its value and exp() of it, the value of v, is returned."""
    asList = lambda v: list(v) if isinstance(v, (list, tuple)) else [v]
    vNames, lNames = asList(vNames), asList(lNames)
    local, parse = _make_local_parse([str(expr)] + [str(x) for x in vNames + lNames])
    e = parse(str(expr))
    vp = {str(l): spy.Symbol(str(v), positive=True) for v, l in zip(vNames, lNames)}
    e = e.subs({parse(l): spy.log(v) for l, v in vp.items()})
    if solveFor is not None and str(solveFor) in vp:
        e = spy.exp(e)
    e = spy.powsimp(spy.expand_log(e, force=True))
    return str(e.subs({v: spy.Symbol(str(v)) for v in vp.values()}))


def logChartBack(expr, xNames, thetas, bases, solveFor=None):
    """An expression of the chart X = b^theta back in theta. With `solveFor` naming
    an X, the expression is the value of that X and log_b of it is returned, the
    value of theta."""
    asList = lambda v: list(v) if isinstance(v, (list, tuple)) else [v]
    xNames, thetas, bases = asList(xNames), asList(thetas), asList(bases)
    local, parse = _make_local_parse([str(expr)] + [str(x) for x in xNames + thetas])
    e = parse(str(expr))
    subs = {}
    logb = None
    for X, th, b in zip(xNames, thetas, bases):
        bb = spy.E if str(b) == 'E' else spy.sympify(str(b))
        subs[parse(str(X))] = bb ** parse(str(th))
        if solveFor is not None and str(X) == str(solveFor):
            logb = spy.log(bb)
    e = e.subs(subs)
    if logb is not None:
        e = spy.cancel(spy.expand_log(spy.log(e), force=True) / logb)
    return str(e)


def logParamBacksub(expr, comp, xNames, thetas, bases):
    """A reported component in the rational coordinates X = base^theta back in theta:
    X -> base^theta everywhere, and the component of X itself divided by X*log(base)
    when `comp` names an X."""
    asList = lambda v: list(v) if isinstance(v, (list, tuple)) else [v]
    xNames, thetas, bases = asList(xNames), asList(thetas), asList(bases)
    local, parse = _make_local_parse([str(expr)] + [str(x) for x in xNames + thetas])
    e = parse(str(expr))
    subs = {}
    for X, th, b in zip(xNames, thetas, bases):
        bb = spy.E if str(b) == 'E' else spy.sympify(str(b))
        Xs = parse(str(X))
        if str(comp) == str(X):
            e = e / (Xs * spy.log(bb))
        subs[Xs] = bb ** parse(str(th))
    return str(spy.powsimp(spy.cancel(e.subs(subs))))


def _detect_power_atoms(perCond, states=()):
    """Find base^exp terms with a non-numeric exponent. base must be a single
    symbol and exp = c*n (c rational, n a symbol other than a state); returns the
    unique (base, n) pairs, or None if a power is outside this form (then it stays
    non-rational). A numeric base is an exponential, see _apply_exp_recast."""
    pairs = set()
    states = set(states)
    for (f_c, g_c, ic_c, f_ss) in perCond:
        for e in list(f_c) + list(g_c) + list(ic_c.values()) + list(f_ss):
            for pw in spy.sympify(e).atoms(spy.Pow):
                if pw.exp.is_number or _exp_atom(pw) is not None:
                    continue
                if not pw.base.is_Symbol:
                    return None
                c, rest = pw.exp.as_coeff_Mul()
                if not (rest.is_Symbol and c.is_rational) or rest in states:
                    return None
                pairs.add((pw.base, rest))
    return sorted(pairs, key=lambda t: (str(t[0]), str(t[1])))


def _apply_power_recast(pairs, S, perCond):
    """Recast each base^(c*n) as E^c with E = base^n a new state, and add a
    companion L = log(base) per base. E' = n*E*base'/base, L' = base'/base when
    base is a state (0 when it is a parameter). E, L are appended to the state
    list; f_ss keeps only the real states (E is held generic in the solve)."""
    recast = {}
    Lof = {}
    for (base, n) in pairs:
        recast[(base, n)] = spy.Symbol('_E_%s_%s' % (base, n))
        if base not in Lof:
            Lof[base] = spy.Symbol('_L_%s' % base)

    def sub(expr):
        def repl(pw):
            if pw.exp.is_number or not pw.base.is_Symbol:
                return pw
            c, rest = pw.exp.as_coeff_Mul()
            E = recast.get((pw.base, rest))
            return E ** c if E is not None else pw
        return spy.sympify(expr).replace(lambda x: x.is_Pow, repl)

    extra = [recast[(b, n)] for (b, n) in pairs] + \
            [Lof[b] for b in sorted(Lof, key=str)]
    newPerCond = []
    for (f_c, g_c, ic_c, f_ss) in perCond:
        f_r = [sub(e) for e in f_c]
        g_r = [sub(e) for e in g_c]
        ic_r = {k: sub(v) for k, v in ic_c.items()}
        ss_r = [sub(e) for e in f_ss]
        rhsOf = {str(X): f_r[i] for i, X in enumerate(S)}
        for (b, n) in pairs:
            E = recast[(b, n)]
            brhs = rhsOf.get(str(b), spy.Integer(0))
            ic_r[str(E)] = E
            f_r.append(n * E * brhs / b if brhs != 0 else spy.Integer(0))
        for b in sorted(Lof, key=str):
            L = Lof[b]
            brhs = rhsOf.get(str(b), spy.Integer(0))
            ic_r[str(L)] = L
            f_r.append(brhs / b if brhs != 0 else spy.Integer(0))
        newPerCond.append((f_r, g_r, ic_r, ss_r))
    meta = [{'E': str(recast[(b, n)]), 'L': str(Lof[b]),
             'base': str(b), 'exp': str(n)} for (b, n) in pairs]
    return list(S) + extra, newPerCond, meta


def recastBacksub(expr, eNames, lNames, bases, exps):
    """Substitute recast coordinates back, E -> base**exp and L -> log(base), and
    cancel. The name vectors are aligned per recast atom; the shared symbol table
    keeps names like E, I, N as symbols."""
    asList = lambda v: list(v) if isinstance(v, (list, tuple)) else [v]
    eNames, lNames, bases, exps = (asList(eNames), asList(lNames),
                                   asList(bases), asList(exps))
    local, parse = _make_local_parse(
        [str(expr)] + [str(x) for x in eNames + lNames + bases + exps])
    e = parse(str(expr))
    subs = {}
    for E, L, base, exp in zip(eNames, lNames, bases, exps):
        b = parse(str(base))
        subs[parse(str(E))] = b ** parse(str(exp))
        subs[parse(str(L))] = spy.log(b)
    return str(spy.cancel(e.subs(subs)))


# ---- exponentials of states: auxiliary states and generic exponential leaves ---------
#
# b^u with u in the states is carried by an auxiliary state X = r^phi (r the canonical
# base, phi = t/D a term of u), X' = log(r) X phi'. Its initial value is rational in
# generic leaves W = r^(tau/L), one per term tau of the initial exponents, tied to the
# other leaves by dW = W log(r) d(tau/L) (stacked onto the codistribution in R).
# Exact at a generic point: exponentials of Q-linearly independent terms are
# algebraically independent over the rational functions (Ax 1971). e^c and
# log(prime) are fixed leaves.

def _canon_exp_base(b):
    """(r, k) with b = r^k, r = E or a positive rational that is not a perfect power;
    None for any other base."""
    if b == spy.E:
        return spy.E, 1
    b = spy.nsimplify(b, rational=True)
    if not (b.is_Rational and b.is_positive) or b == 1:
        return None
    sign = 1
    if b < 1:
        b, sign = 1 / b, -1
    fp, fq = spy.factorint(b.p), spy.factorint(b.q)
    g = 0
    for e in list(fp.values()) + list(fq.values()):
        g = math.gcd(g, int(e))
    num = den = 1
    for pr, e in fp.items():
        num *= pr ** (e // g)
    for pr, e in fq.items():
        den *= pr ** (e // g)
    return spy.Rational(num, den), sign * g


def _exp_terms(b, u):
    """b^u = r^(c0 + sum q*t): (r, c0, {t: q}) with q rational, or None."""
    rk = _canon_exp_base(b)
    if rk is None:
        return None
    r, k = rk
    u = spy.expand(spy.nsimplify(k * u, rational=True))
    terms = dict(u.as_coefficients_dict())
    c0 = terms.pop(spy.Integer(1), spy.Integer(0))
    if not (c0.is_Rational and all(q.is_Rational for q in terms.values())):
        return None
    return r, c0, terms


def _log_basis(r):
    """log(r) as {1: 1} for r = E, else {prime: exponent}."""
    if r == spy.E:
        return {1: 1}
    d = dict(spy.factorint(r.p))
    for pr, e in spy.factorint(r.q).items():
        d[pr] = d.get(pr, 0) - e
    return d


def _terms_independent(atoms):
    """True if 1 and the exponents tau*log(r) of the atoms [(r, tau)] are linearly
    independent over Q, with 1 and the logs of the primes independent. Checked at
    random integer points."""
    if not atoms:
        return True
    comps = [_log_basis(r) for r, _ in atoms]
    basis = sorted({b for d in comps for b in d}, key=str)
    taus = [t for _, t in atoms]
    syms = sorted(set().union(*[t.free_symbols for t in taus]), key=str)
    m, nb = len(atoms), len(basis)
    rng = np.random.default_rng(7)
    rows, npts = [], 0
    for _ in range(4 * (m + 3)):
        pt = {s: spy.Integer(int(v)) for s, v in zip(syms, rng.integers(2, 997, len(syms)))}
        vals = [t.xreplace(pt) for t in taus]
        vals = [v if v.is_Rational else spy.nsimplify(v) for v in vals]
        if not all(v.is_Rational for v in vals):
            continue
        for bi, b in enumerate(basis):
            rows.append([comps[k].get(b, 0) * vals[k] for k in range(m)] +
                        [1 if j == bi else 0 for j in range(nb)])
        npts += 1
        if npts >= m + 3:
            break
    if npts < m + 1:
        return False
    # full rank over GF(p) implies full rank over Q
    for p in (2147483629, 2147483587):
        if any(spy.Rational(v).q % p == 0 for r in rows for v in r):
            continue
        modRows = [[_modp_rational(v, p) for v in r] for r in rows]
        if _rank_mod(modRows, p) == m + nb:
            return True
    return spy.Matrix(rows).rank() == m + nb


class _ExpCtx:
    """Fresh symbols and constant leaves shared by one exponential rewrite: e^(1/L0)
    as `ec` (a power of it until L0 is known), log(p) as one leaf per prime p."""
    def __init__(self, taken):
        self.taken = set(taken)
        self.ec = self.fresh('_e_', positive=True)
        self.lnSym = {}

    def fresh(self, stem, **kw):
        nm = stem
        while nm in self.taken:
            nm += '_'
        self.taken.add(nm)
        return spy.Symbol(nm, **kw)

    def lnOf(self, r):
        """log(r) through one leaf per prime."""
        out = spy.Integer(0)
        for pr, e in _log_basis(r).items():
            if pr == 1:
                out += e
                continue
            if pr not in self.lnSym:
                self.lnSym[pr] = self.fresh('_ln%d_' % pr, positive=True)
            out += e * self.lnSym[pr]
        return out

    @staticmethod
    def rpow(r, e):
        return spy.exp(e) if r == spy.E else r ** e

    def const(self, r, c0):
        """r^c0, None when it is irrational and not a power of e."""
        if c0 == 0:
            return spy.Integer(1)
        if r == spy.E:
            return self.ec ** c0
        return r ** c0 if c0.is_Integer else None

    def numbers(self, e):
        """E as the leaf ec, log(b) of a number as a multiple of a log leaf."""
        e = spy.sympify(e).xreplace({spy.E: self.ec})
        rep = {}
        for L in e.atoms(spy.log):
            if L.args[0].is_number:
                rk = _canon_exp_base(L.args[0])
                if rk is not None:
                    rep[L] = rk[1] * self.lnOf(rk[0])
        return e.xreplace(rep) if rep else e


def _exp_leaf_canon(exprs, ctx):
    """Every exponential in `exprs`, innermost first, as a product of generic leaves
    W = r^(tau/L), one per term tau, and powers of ctx.ec. Returns {'exprs', 'atoms'}
    with atoms [{'W', 'r', 'tau', 'L'}], after ec -> ec^L0; or {'why': ...}."""
    exprs = [ctx.numbers(e) for e in exprs]
    atoms = []
    while True:
        found = {}
        for e in exprs:
            for at in e.atoms(spy.exp, spy.Pow):
                if at in found or _exp_atom(at) is None:
                    continue
                if any(_exp_atom(a) is not None
                       for a in _exp_atom(at)[1].atoms(spy.exp, spy.Pow)):
                    continue
                dec = _exp_terms(*_exp_atom(at))
                if dec is None or ctx.const(dec[0], dec[1]) is None:
                    return {'why': 'the exponential %s is not supported by '
                            'symEngine = "modular"; try symEngine = "symbolic"' % at}
                found[at] = dec
        if not found:
            break
        den = {}
        for r, c0, terms in found.values():
            for t, q in terms.items():
                den[(r, t)] = spy.ilcm(den.get((r, t), 1), q.q)
        keyW = {}
        for (r, t), L in sorted(den.items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]))):
            W = ctx.fresh('_ew%d_' % (len(atoms) + 1), positive=True)
            keyW[(r, t)] = W
            atoms.append({'W': W, 'r': r, 'tau': t, 'L': L})
        repl = {}
        for at, (r, c0, terms) in found.items():
            v = ctx.const(r, c0)
            for t, q in terms.items():
                v = v * keyW[(r, t)] ** int(q * den[(r, t)])
            repl[at] = v
        exprs = [e.xreplace(repl) for e in exprs]
    # e^c as integer powers of the leaf e^(1/L0)
    L0 = 1
    for e in exprs + [w['tau'] for w in atoms]:
        for pw in e.atoms(spy.Pow):
            if pw.base == ctx.ec:
                L0 = spy.ilcm(L0, spy.Rational(pw.exp).q)
    ecSub = {ctx.ec: ctx.ec ** L0}
    exprs = [e.xreplace(ecSub) for e in exprs]
    for w in atoms:
        w['tau'] = w['tau'].xreplace(ecSub)
    ctx.L0 = L0
    if not _terms_independent([(w['r'], w['tau']) for w in atoms]):
        return {'why': 'exponents that are linearly dependent over the rationals '
                'are not supported'}
    return {'exprs': exprs, 'atoms': atoms}


def _exp_back(atoms, ctx, used):
    """Back-substitution map {name: value string} for the leaves and constants."""
    back = {str(w['W']): str(ctx.rpow(w['r'], w['tau'] / w['L'])) for w in atoms}
    if ctx.ec in used:
        back[str(ctx.ec)] = str(spy.exp(spy.Rational(1, ctx.L0)))
    for r, s in ctx.lnSym.items():
        back[str(s)] = 'log(%s)' % r
    return back


def _apply_exp_recast(S, perCond, evPer, tmPer, taken):
    """Replace every exponential in the per-condition model by auxiliary states and
    generic exponential leaves (see above). `evPer` holds each condition's later
    events as {'var', 'value', 'method'} with sympy values, `tmPer` its segment start
    time or None. Returns a dict with the extended S, perCond, evPer and tmPer, the
    relation list [(W, z, expr)], the fixed constant leaves, the back-substitution map
    and the codimension added; or {'why': ...} for an unsupported form."""
    ctx = _ExpCtx(taken)
    rpow, const = ctx.rpow, ctx.const
    S = list(S)
    K = len(perCond)
    f = [[ctx.numbers(e) for e in p[0]] for p in perCond]
    g = [[ctx.numbers(e) for e in p[1]] for p in perCond]
    ic = [{k: ctx.numbers(v) for k, v in p[2].items()} for p in perCond]
    fss = [list(p[3]) for p in perCond]
    ev = [[dict(e, value=ctx.numbers(e['value'])) for e in evs] for evs in evPer]
    tm = [None if t is None else ctx.numbers(t) for t in tmPer]

    # ---- dynamics: auxiliary states, innermost exponential first
    aux = []                         # {'X', 'r', 'phi'}
    nReal = len(S)
    while True:
        states = set(S)

        def stateAtom(a):
            be = _exp_atom(a)
            return be is not None and bool(be[1].free_symbols & states)
        found = {}
        for c in range(K):
            for e in f[c] + g[c]:
                for at in e.atoms(spy.exp, spy.Pow):
                    if at in found or not stateAtom(at):
                        continue
                    if any(stateAtom(a) for a in _exp_atom(at)[1].atoms(spy.exp, spy.Pow)):
                        continue
                    dec = _exp_terms(*_exp_atom(at))
                    if dec is None or const(dec[0], dec[1]) is None:
                        return {'why': 'the exponential %s is not supported by '
                            'symEngine = "modular"; try symEngine = "symbolic"' % at}
                    found[at] = dec
        if not found:
            break
        den = {}
        for r, c0, terms in found.values():
            for t, q in terms.items():
                if t.free_symbols & states:
                    den[(r, t)] = spy.ilcm(den.get((r, t), 1), q.q)
        keyX = {}
        for (r, t), D in sorted(den.items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]))):
            X = ctx.fresh('_ex%d_' % (len(aux) + 1))
            keyX[(r, t)] = X
            aux.append({'X': X, 'r': r, 'phi': t / D})
        repl = {}
        for at, (r, c0, terms) in found.items():
            v = const(r, c0)
            for t, q in terms.items():
                if (r, t) in keyX:
                    v = v * keyX[(r, t)] ** int(q * den[(r, t)])
                else:
                    v = v * rpow(r, q * t)
            repl[at] = v
        f = [[e.xreplace(repl) for e in fc] for fc in f]
        g = [[e.xreplace(repl) for e in gc] for gc in g]
        newAux = aux[len(aux) - len(keyX):]
        for c in range(K):
            rhsOf = dict(zip(S, f[c]))
            for a in newAux:
                phi = a['phi']
                dphi = sum(spy.diff(phi, s) * rhsOf[s] for s in phi.free_symbols & states)
                f[c].append(spy.expand(a['X'] * ctx.lnOf(a['r']) * dphi))
                icSub = {s: ic[c][str(s)] for s in phi.free_symbols & states}
                ic[c][str(a['X'])] = rpow(a['r'], phi.xreplace(icSub))
                fss[c].append(spy.Integer(0))
        S += [a['X'] for a in newAux]

    # an auxiliary state no observable depends on stays at 1 in that condition
    idx = {X: i for i, X in enumerate(S)}
    for c in range(K):
        used = set()
        front = set()
        for e in f[c][:nReal] + g[c]:
            front |= e.free_symbols
        while front:
            nxt = set()
            for a in aux:
                if a['X'] in front and a['X'] not in used:
                    used.add(a['X'])
                    nxt |= f[c][idx[a['X']]].free_symbols
            front = nxt - used
        for a in aux:
            if a['X'] not in used:
                f[c][idx[a['X']]] = spy.Integer(0)
                ic[c][str(a['X'])] = spy.Integer(1)

    # later events on a state inside an exponent move the auxiliary states with it
    Sset = set(S)
    for c in range(K):
        out = []
        for e in ev[c]:
            out.append(e)
            queue = [e]
            while queue:
                cur = queue.pop(0)
                for a in aux:
                    phi = a['phi']
                    vs = [s for s in phi.free_symbols if str(s) == str(cur['var'])]
                    if not vs or f[c][idx[a['X']]] == 0:
                        continue
                    vsym = vs[0]
                    others = (phi.free_symbols & Sset) - {vsym}
                    d1 = spy.diff(phi, vsym)
                    if cur['method'] == 'replace' and not others:
                        new = {'var': str(a['X']), 'method': 'replace',
                               'value': rpow(a['r'], phi.xreplace({vsym: cur['value']}))}
                    elif (cur['method'] == 'add' and spy.diff(d1, vsym) == 0 and
                          not (d1.free_symbols & Sset)):
                        new = {'var': str(a['X']), 'method': 'multiply',
                               'value': rpow(a['r'], cur['value'] * d1)}
                    else:
                        return {'why': 'an event (%s) on %s, which enters an exponent, '
                                'is not supported' % (cur['method'], cur['var'])}
                    out.append(new)
                    queue.append(new)
        ev[c] = out

    # ---- leaf level: generic exponential leaves W
    refs = ([('f', c, i) for c in range(K) for i in range(len(f[c]))] +
            [('g', c, i) for c in range(K) for i in range(len(g[c]))] +
            [('ic', c, k) for c in range(K) for k in ic[c]] +
            [('ev', c, i) for c in range(K) for i in range(len(ev[c]))] +
            [('tm', c, None) for c in range(K) if tm[c] is not None])
    store = {'f': f, 'g': g, 'ic': ic}
    vals = []
    for kind, c, i in refs:
        vals.append(ev[c][i]['value'] if kind == 'ev' else tm[c] if kind == 'tm'
                    else store[kind][c][i])
    lc = _exp_leaf_canon(vals, ctx)
    if 'why' in lc:
        return lc
    for (kind, c, i), v in zip(refs, lc['exprs']):
        if kind == 'ev':
            ev[c][i]['value'] = v
        elif kind == 'tm':
            tm[c] = v
        else:
            store[kind][c][i] = v
    wAtoms = lc['atoms']

    constSyms = set(ctx.lnSym.values())
    used = set()
    for v in lc['exprs'] + [w['tau'] for w in wAtoms]:
        used |= v.free_symbols
    if ctx.ec in used:
        constSyms.add(ctx.ec)
    rel = []
    for w in wAtoms:
        phi = w['tau'] / w['L']
        for z in sorted(phi.free_symbols - constSyms, key=str):
            rel.append((str(w['W']), str(z),
                        spy.together(w['W'] * ctx.lnOf(w['r']) * spy.diff(phi, z))))
    perCond = [(f[c], g[c], ic[c], fss[c]) for c in range(K)]
    return {'S': S, 'perCond': perCond, 'evPer': ev, 'tmPer': tm, 'rel': rel,
            'consts': sorted(str(s) for s in constSyms),
            'atoms': [str(w['W']) for w in wAtoms], 'back': _exp_back(wAtoms, ctx, used),
            'nAux': len(aux), 'codim': len(aux) + len(wAtoms)}


# ---- trigonometric functions: half-angle states and generic half-angle leaves --------
#
# sin, cos and tan of u = c0 + k*pi/2 + sum q*t are rational in T = tan(phi/2), phi =
# t/D a term of u, through e^(i*phi) = (1 + i*T)/(1 - i*T). A term in the states is
# carried by an auxiliary state T' = (1 + T^2)/2 phi', a term in the leaves by a
# generic leaf V with dV = (1 + V^2)/2 d(phi) (stacked like an exponential leaf), a
# rational c0 by the fixed leaf tan(1/(2*M0)). Exact at a generic point by Ax (1971)
# applied to e^(i*phi); real and imaginary exponents are independent of each other.

_TRIG = (spy.sin, spy.cos, spy.tan)
_TRIG_RECIPROCAL = {spy.cot: lambda a: spy.cos(a) / spy.sin(a),
                    spy.sec: lambda a: 1 / spy.cos(a),
                    spy.csc: lambda a: 1 / spy.sin(a)}


def _has_trig(e):
    e = spy.sympify(e)
    return e.has(*_TRIG) or e.has(*_TRIG_RECIPROCAL)


def _trig_basic(e, piSym):
    """cot, sec and csc through sin and cos, a model symbol pi as the number."""
    e = spy.sympify(e)
    if piSym is not None:
        e = e.xreplace({piSym: spy.pi})
    return e.replace(lambda x: type(x) in _TRIG_RECIPROCAL,
                     lambda x: _TRIG_RECIPROCAL[type(x)](x.args[0]))


def _trig_terms(fn, u):
    """fn(u) with u = c0 + cp*pi + sum q*t as (c0, cp, {t: q}), all rational; tan in
    the doubled argument, tan(u) = sin(2u)/(1 + cos(2u)). None for another form."""
    if fn is spy.tan:
        u = 2 * u
    u = spy.expand(spy.nsimplify(u, rational=True))
    terms = dict(u.as_coefficients_dict())
    c0 = terms.pop(spy.Integer(1), spy.Integer(0))
    cp = terms.pop(spy.pi, spy.Integer(0))
    if not (c0.is_Rational and cp.is_Rational and (2 * cp).is_Integer and
            all(q.is_Rational for q in terms.values())):
        return None
    # an inner trig atom left in a term is one the later leaf pass resolves
    if any(_has_exp(t) or not _is_rational_expr(
            t.replace(lambda x: isinstance(x, _TRIG), lambda x: spy.Dummy())
             .xreplace({spy.pi: 1})) for t in terms):
        return None
    return c0, cp, terms


def _trig_value(fn, factors, k, rest):
    """fn(u) for e^(i*v) = i^k e^(i*rest) prod ((1 + i*T)/(1 - i*T))^m over factors
    [(T, m)], v = u or 2u for tan; cos and sin of rest stay as they are."""
    I = spy.I
    reps, A, n = {}, spy.Integer(1), spy.Integer(1)
    for T, m in factors:
        d = spy.Dummy(real=True)
        reps[d] = T
        A *= (1 + I * d) ** m if m > 0 else (1 - I * d) ** (-m)
        n *= (1 + T ** 2) ** abs(m)
    Z = I ** (int(k) % 4) * spy.expand(A ** 2)
    if rest != 0:
        C, Sn = spy.Dummy(real=True), spy.Dummy(real=True)
        reps[C], reps[Sn] = spy.cos(rest), spy.sin(rest)
        Z = Z * (C + I * Sn)
    re, im = spy.expand(Z).as_real_imag()
    re, im = re.xreplace(reps), im.xreplace(reps)
    if fn is spy.tan:
        return spy.cancel(im / (n + re))
    return spy.cancel((im if fn is spy.sin else re) / n)


def _trig_find(exprs, inScope, why):
    """Innermost trig atoms of `exprs` whose argument `inScope` accepts, with their
    decomposition; {'why': ...} for an unsupported argument."""
    found = {}
    for e in exprs:
        for at in e.atoms(*_TRIG):
            if at in found or not inScope(at):
                continue
            if any(inScope(a) for a in at.args[0].atoms(*_TRIG)):
                continue
            dec = _trig_terms(type(at), at.args[0])
            if dec is None:
                return {'why': why % at}
            found[at] = dec
    return found


def _trig_den(found, keep):
    """Per term t that `keep` takes the least D with q*D integer over all atoms."""
    den = {}
    for c0, cp, terms in found.values():
        for t, q in terms.items():
            if keep(t):
                den[t] = spy.ilcm(den.get(t, 1), q.q)
    return den


def _apply_trig_recast(S, perCond, evPer, tmPer, taken):
    """Replace every sin, cos and tan (and cot, sec, csc) in the per-condition model by
    auxiliary half-angle states and generic half-angle leaves (see above). Arguments
    and return as _apply_exp_recast, or {'why': ...} for an unsupported form."""
    ctx = _ExpCtx(taken)
    pis = ctx.fresh('_pi_', positive=True)
    piSym = next((s for p in perCond
                  for e in list(p[0]) + list(p[1]) + list(p[2].values())
                  for s in spy.sympify(e).free_symbols if str(s) == 'pi'), None)
    if piSym is not None and piSym in set(S):
        piSym = None
    tb = lambda e: _trig_basic(e, piSym)
    S = list(S)
    K = len(perCond)
    f = [[tb(e) for e in p[0]] for p in perCond]
    g = [[tb(e) for e in p[1]] for p in perCond]
    ic = [{k: tb(v) for k, v in p[2].items()} for p in perCond]
    fss = [[tb(e).xreplace({spy.pi: pis}) for e in p[3]] for p in perCond]
    ev = [[dict(e, value=tb(e['value'])) for e in evs] for evs in evPer]
    tm = [None if t is None else tb(t) for t in tmPer]
    unsupported = ('the argument of %s is not supported by symEngine = "modular"; '
                   'try symEngine = "symbolic"')

    # ---- dynamics: auxiliary half-angle states, innermost first
    aux = []                         # {'X', 'phi'}
    while True:
        states = set(S)
        stateAtom = lambda a: bool(a.args[0].free_symbols & states)
        found = _trig_find([e for c in range(K) for e in f[c] + g[c]], stateAtom,
                           unsupported)
        if 'why' in found:
            return found
        if not found:
            break
        den = _trig_den(found, lambda t: bool(t.free_symbols & states))
        keyX = {}
        for t, D in sorted(den.items(), key=lambda kv: str(kv[0])):
            X = ctx.fresh('_tx%d_' % (len(aux) + 1))
            keyX[t] = X
            aux.append({'X': X, 'phi': t / D})
        repl = {}
        for at, (c0, cp, terms) in found.items():
            factors = [(keyX[t], int(q * den[t]))
                       for t, q in terms.items() if t in keyX]
            rest = c0 + sum((q * t for t, q in terms.items() if t not in keyX),
                            spy.Integer(0))
            repl[at] = _trig_value(type(at), factors, 2 * cp, rest)
        f = [[e.xreplace(repl) for e in fc] for fc in f]
        g = [[e.xreplace(repl) for e in gc] for gc in g]
        newAux = aux[len(aux) - len(keyX):]
        for c in range(K):
            rhsOf = dict(zip(S, f[c]))
            for a in newAux:
                phi = a['phi']
                dphi = sum(spy.diff(phi, s) * rhsOf[s]
                           for s in phi.free_symbols & states)
                f[c].append(spy.expand((1 + a['X'] ** 2) * dphi / 2))
                icSub = {s: ic[c][str(s)] for s in phi.free_symbols & states}
                ic[c][str(a['X'])] = spy.tan(phi.xreplace(icSub) / 2)
                fss[c].append(spy.Integer(0))
        S += [a['X'] for a in newAux]

    # an auxiliary state no observable depends on stays at 0 in that condition
    nReal = len(S) - len(aux)
    idx = {X: i for i, X in enumerate(S)}
    for c in range(K):
        used, front = set(), set()
        for e in f[c][:nReal] + g[c]:
            front |= e.free_symbols
        while front:
            nxt = set()
            for a in aux:
                if a['X'] in front and a['X'] not in used:
                    used.add(a['X'])
                    nxt |= f[c][idx[a['X']]].free_symbols
            front = nxt - used
        for a in aux:
            if a['X'] not in used:
                f[c][idx[a['X']]] = spy.Integer(0)
                ic[c][str(a['X'])] = spy.Integer(0)

    # a later replace event on a state inside an argument resets the auxiliary state
    Sset = set(S)
    for c in range(K):
        out = []
        for e in ev[c]:
            out.append(e)
            for a in aux:
                phi = a['phi']
                vs = [s for s in phi.free_symbols if str(s) == str(e['var'])]
                if not vs or f[c][idx[a['X']]] == 0:
                    continue
                if e['method'] != 'replace' or (phi.free_symbols & Sset) - {vs[0]}:
                    return {'why': 'an event (%s) on %s, which enters the argument of '
                            'a trigonometric function, is not supported'
                            % (e['method'], e['var'])}
                out.append({'var': str(a['X']), 'method': 'replace',
                            'value': spy.tan(phi.xreplace({vs[0]: e['value']}) / 2)})
        ev[c] = out

    # ---- leaf level: generic half-angle leaves V and the constant leaves
    refs = ([('f', c, i) for c in range(K) for i in range(len(f[c]))] +
            [('g', c, i) for c in range(K) for i in range(len(g[c]))] +
            [('ic', c, k) for c in range(K) for k in ic[c]] +
            [('ev', c, i) for c in range(K) for i in range(len(ev[c]))] +
            [('tm', c, None) for c in range(K) if tm[c] is not None])
    store = {'f': f, 'g': g, 'ic': ic}
    vals = [ev[c][i]['value'] if kind == 'ev' else tm[c] if kind == 'tm'
            else store[kind][c][i] for kind, c, i in refs]
    # rational offsets c0 as powers of e^(i/M0), M0 over all arguments
    atoms, M0 = [], 1
    for e in vals:
        for at in e.atoms(*_TRIG):
            dec = _trig_terms(type(at), at.args[0])
            if dec is not None:
                M0 = spy.ilcm(M0, dec[0].q)
    tc = ctx.fresh('_tc_', positive=True)
    while True:
        found = _trig_find(vals, lambda a: True, unsupported)
        if 'why' in found:
            return found
        if not found:
            break
        den = _trig_den(found, lambda t: True)
        keyV = {}
        for t, D in sorted(den.items(), key=lambda kv: str(kv[0])):
            V = ctx.fresh('_tw%d_' % (len(atoms) + 1))
            keyV[t] = V
            atoms.append({'V': V, 'phi': t / D})
        repl = {}
        for at, (c0, cp, terms) in found.items():
            factors = [(keyV[t], int(q * den[t])) for t, q in terms.items()]
            if c0 != 0:
                if not (c0 * M0).is_Integer:
                    return {'why': unsupported % at}
                factors.append((tc, int(c0 * M0)))
            repl[at] = _trig_value(type(at), factors, 2 * cp, spy.Integer(0))
        vals = [e.xreplace(repl) for e in vals]
    if not _terms_independent([(spy.E, a['phi'].xreplace({spy.pi: spy.Symbol('_pi_')}))
                               for a in atoms]):
        return {'why': 'arguments of trigonometric functions that are linearly '
                'dependent over the rationals are not supported'}

    vals = [e.xreplace({spy.pi: pis}) for e in vals]
    for a in atoms:
        a['phi'] = a['phi'].xreplace({spy.pi: pis})
    for (kind, c, i), v in zip(refs, vals):
        if kind == 'ev':
            ev[c][i]['value'] = v
        elif kind == 'tm':
            tm[c] = v
        else:
            store[kind][c][i] = v

    used = set()
    for v in vals + [a['phi'] for a in atoms]:
        used |= v.free_symbols
    consts = {pis, tc} & used
    rel = []
    for a in atoms:
        for z in sorted(a['phi'].free_symbols - consts, key=str):
            rel.append((str(a['V']), str(z),
                        spy.together((1 + a['V'] ** 2) * spy.diff(a['phi'], z) / 2)))
    back = {str(a['V']): str(spy.tan(a['phi'].xreplace({pis: spy.pi}) / 2))
            for a in atoms}
    if pis in consts:
        back[str(pis)] = 'pi'
    if tc in consts:
        back[str(tc)] = str(spy.tan(spy.Rational(1, 2 * M0)))
    perCond = [(f[c], g[c], ic[c], fss[c]) for c in range(K)]
    return {'S': S, 'perCond': perCond, 'evPer': ev, 'tmPer': tm, 'rel': rel,
            'consts': sorted(str(s) for s in consts),
            'atoms': [str(a['V']) for a in atoms], 'back': back,
            'nAux': len(aux), 'codim': len(aux) + len(atoms)}


def _merge_recast(a, b):
    """Two recast results as one: b applied after a."""
    if a is None or b is None:
        return b if a is None else a
    back = dict(a['back'])
    back.update(b['back'])
    return dict(b, rel=a['rel'] + b['rel'], atoms=a['atoms'] + b['atoms'], back=back,
                consts=sorted(set(a['consts']) | set(b['consts'])),
                nAux=a['nAux'] + b['nAux'], codim=a['codim'] + b['codim'])


def _tidy_logs(e):
    """Common factors pulled out and each numeric sum of logs as one log."""
    e = spy.factor_terms(spy.cancel(e))
    return e.replace(lambda x: x.is_Add and x.is_number and x.has(spy.log),
                     lambda x: spy.logcombine(x, force=True))


def dropMonomialContent(vector, skip=None):
    """A direction {coordinate: component} with denominators cleared, divided by the
    power product of symbols every component carries and signed so that the first
    coordinate has a positive component: the same direction. Components in `skip`
    are scaled along but do not set the factor."""
    skip = set(_as_list(skip)) if skip is not None else set()
    items = [(str(k), str(v)) for k, v in dict(vector).items()]
    local, parse = _make_local_parse([v for _, v in items])
    es = [spy.cancel(parse(v)) for _, v in items]
    L = spy.Integer(1)
    for (k, _), e in zip(items, es):
        if k not in skip:
            L = spy.lcm(L, spy.fraction(e)[1])
    es = [spy.factor_terms(spy.cancel(e * L)) for e in es]
    first = [e for (k, _), e in sorted(zip(items, es)) if e != 0 and k not in skip]
    if first and first[0].could_extract_minus_sign():
        es = [-e for e in es]
    common = None
    for (k, _), e in zip(items, es):
        if k in skip:
            continue
        pw = {}
        num, den = spy.fraction(e)
        for sgn, part in ((1, num), (-1, den)):
            for f in spy.Mul.make_args(part):
                b, x = f.as_base_exp()
                if b.is_Symbol and x.is_Integer:
                    pw[b] = pw.get(b, 0) + sgn * int(x)
        common = pw if common is None else {
            b: (min(x, common[b]) if x > 0 else max(x, common[b]))
            for b, x in pw.items() if b in common and (x > 0) == (common[b] > 0)}
    m = spy.Integer(1)
    for b, x in (common or {}).items():
        m *= b ** x
    return {k: str(spy.cancel(e / m)) for (k, _), e in zip(items, es)}


def expBacksub(expr, names, values):
    """A reported component in the exponential leaves back in the model symbols:
    each name is replaced by its value until none is left."""
    asList = lambda v: list(v) if isinstance(v, (list, tuple)) else [v]
    names, values = asList(names), asList(values)
    local, parse = _make_local_parse([str(expr)] + [str(v) for v in values] +
                                     [str(n) for n in names])
    e = parse(str(expr))
    sub = {parse(str(n)): parse(str(v)) for n, v in zip(names, values)}
    for _ in range(len(sub) + 1):
        e2 = e.xreplace(sub)
        if e2 == e:
            break
        e = e2
    return str(_tidy_logs(spy.powsimp(spy.cancel(e))))


def trigFullAngle(vector, skip=None):
    """A direction {coordinate: component} over half-angle leaves V = tan(u/2) in
    sin(u) and cos(u): denominators cleared and common factors divided out in V, then
    divided by (1 + V^2)^k with the least k that makes every component a polynomial
    in sin(u) and cos(u). Components in `skip` are scaled along but do not set the
    factor."""
    skip = set(_as_list(skip)) if skip is not None else set()
    items = [(str(k), str(v)) for k, v in dict(vector).items()]
    local, parse = _make_local_parse([v for _, v in items])
    es = [spy.sympify(parse(v)) for _, v in items]
    halves = sorted(set().union(*[e.atoms(spy.tan) for e in es]), key=str)
    if not halves:
        return dict(items)
    Vs = [spy.Dummy('V') for _ in halves]
    es = [spy.cancel(e.xreplace(dict(zip(halves, Vs)))) for e in es]
    keep = [k not in skip for k, _ in items]
    L = spy.Integer(1)
    for e, kp in zip(es, keep):
        if kp:
            L = spy.lcm(L, spy.fraction(e)[1])
    es = [spy.cancel(e * L) for e in es]
    G = spy.Integer(0)
    for e, kp in zip(es, keep):
        if kp and e != 0:
            G = spy.gcd(G, e)
    if G != 0:
        es = [spy.cancel(e / G) for e in es]
    # sin(u) = 2V/(1 + V^2), cos(u) = (1 - V^2)/(1 + V^2), 1 + V^2 = 2/(1 + cos(u))
    rels, sub = [], {}
    for h, V in zip(halves, Vs):
        deg = max((spy.degree(e, V) for e, kp in zip(es, keep) if kp and e != 0),
                  default=0)
        es = [e / (1 + V ** 2) ** ((deg + 1) // 2) for e in es]
        sn, cs = spy.sin(2 * h.args[0]), spy.cos(2 * h.args[0])
        sub[V] = sn / (1 + cs)
        rels.append((sn, cs))
    out = {}
    for (k, _), e in zip(items, es):
        v = spy.cancel(spy.together(e.xreplace(sub)))
        num, den = spy.fraction(v)
        for sn, cs in rels:
            num = spy.expand(num).subs(sn ** 2, 1 - cs ** 2)
            den = spy.expand(den).subs(sn ** 2, 1 - cs ** 2)
        v = spy.factor(spy.cancel(num / den))
        if spy.count_ops(v) < 200:
            v = spy.factor(spy.trigsimp(v))
        out[k] = str(v)
    return out


# ---- observability tape compiler (multi-condition, shared coordinate space) ----------

def _compose_t0_events(icMap, t0evs, pval, subsMap):
    """Apply a segment's t0 events to the initial-state map in place: the jet starts
    at x0 = E(x_ss), the f = 0 constraint stays on the pre-event x_ss. Events on
    non-coordinate states are ignored."""
    for e in t0evs or []:
        Xv = spy.Symbol(str(e['var']))
        if Xv not in icMap:
            continue
        val = spy.sympify(pval(e['value'])).subs(subsMap)
        meth = str(e['method'])
        if meth == 'replace':
            icMap[Xv] = val
        elif meth == 'add':
            icMap[Xv] = icMap[Xv] + val
        elif meth == 'multiply':
            icMap[Xv] = icMap[Xv] * val
    return icMap


def jetGenerator(conds, support, anchor, maxOrder=8, maxChars=200000,
                 timeLimit=60.0):
    """A direction with support S in closed form from the jets.

    Every Lie derivative L_f^k h evaluated at the initial values is invariant under
    the direction X, so X_S is orthogonal to the S-gradient of every jet. With
    |S| - 1 jets whose S-gradients are independent, X_S is their generalised cross
    product: X_i = (-1)^i det(M without column i), M the gradient rows. For |S| = 2
    this is X = (d_b Phi, -d_a Phi), the Hamiltonian field of one jet Phi.

    `conds` holds per condition {'f': ['x = rhs', ...], 'g': [h, ...],
    'ic': ['x = x0', ...]} over the leaves, or a chain {'segments': [...]} whose later
    segments carry their boundary events ('ev') instead of initial values: at order 0
    in every gap length the state does not move across a gap, so a later segment
    starts from the earlier initial values with its events applied, and its jets
    there are invariant like any other (each gap monomial of a row annihilates X).
    A segment whose boundary time depends on the coordinates ends the chain. The
    first segments are tried alone first (fewer, smaller jets); the later ones join
    only when those fall short of the rank.
    Jets are taken order by order, lowest first, rows kept greedily while they raise
    the rank at a random rational point. Returns {'ok': True, 'vector': {z: entry}}
    normalised to 1 at `anchor`, or {'ok': False, 'reason': ...}. The caller
    verifies the result."""
    import time as _time
    t0 = _time.time()
    primary, extra = _jet_expand_chains(conds)
    res = _jet_core(primary, support, anchor, maxOrder, maxChars, timeLimit, t0)
    if not res['ok'] and extra and str(res.get('reason', '')).startswith('rank'):
        res = _jet_core(primary + extra, support, anchor, maxOrder, maxChars,
                        timeLimit, t0)
    return res


def _jet_core(conds, support, anchor, maxOrder, maxChars, timeLimit, t0):
    """jetGenerator() on a fixed list of jet conditions, within timeLimit from t0."""
    import time as _time
    support = [str(s) for s in _as_list(support)]
    anchor = str(anchor)
    m = len(support)
    if m < 2 or anchor not in support:
        return {'ok': False, 'reason': 'support'}
    lines = []
    for c in conds:
        lines += list(_as_list(c['f'])) + list(_as_list(c['ic'])) + \
            [str(x) for x in _as_list(c['g'])]
    lines += support
    local, parse = _make_local_parse(lines)
    zs = [local.get(s, spy.Symbol(s)) for s in support]
    rng = random.Random(20260927)
    point = {}

    def num(e):
        for sym in e.free_symbols:
            if sym not in point:
                point[sym] = spy.Integer(rng.randint(2, 997))
        return e.xreplace(point)

    parsed = []
    for c in conds:
        fl = [l.split(' = ', 1) for l in _as_list(c['f'])]
        states = [local.get(a.strip(), spy.Symbol(a.strip())) for a, _ in fl]
        rhs = [spy.sympify(parse(b)) for _, b in fl]
        icd = {}
        for l in _as_list(c['ic']):
            a, b = l.split(' = ', 1)
            icd[local.get(a.strip(), spy.Symbol(a.strip()))] = spy.sympify(parse(b))
        gs = [spy.sympify(parse(str(x))) for x in _as_list(c['g'])]
        parsed.append((states, rhs, icd, gs))
    rows, numrows = [], []
    rank = 0
    jets = [list(p[3]) for p in parsed]
    for order in range(maxOrder + 1):
        for ci, (states, rhs, icd, gs) in enumerate(parsed):
            for h in jets[ci]:
                if _time.time() - t0 > timeLimit:
                    return {'ok': False, 'reason': 'time limit'}
                he = h.xreplace(icd) if icd else h
                row = [spy.diff(he, z) for z in zs]
                if all(r == 0 for r in row):
                    continue
                # a minor of rows this large does not simplify in bounded time
                if sum(len(str(r)) for r in row) > maxChars // 20:
                    continue
                nr = [num(r) for r in row]
                if spy.Matrix(numrows + [nr]).rank() > rank:
                    rows.append(row); numrows.append(nr); rank += 1
                    if rank == m - 1:
                        break
            if rank == m - 1:
                break
        if rank == m - 1:
            break
        if order == maxOrder:
            break
        for ci, (states, rhs, icd, gs) in enumerate(parsed):
            new = []
            for h in jets[ci]:
                d = sum((spy.diff(h, states[i]) * rhs[i] for i in range(len(states))),
                        spy.Integer(0))
                if len(str(d)) > maxChars:
                    return {'ok': False, 'reason': 'jet too large'}
                new.append(d)
            jets[ci] = new
    if rank < m - 1:
        return {'ok': False, 'reason': 'rank %d of %d within order %d' % (rank, m - 1, maxOrder)}
    M = spy.Matrix(rows)
    comps = []
    for i in range(m):
        sub = M[:, [j for j in range(m) if j != i]]
        # the minor by cofactor expansion, kept as a product/sum tree (no expansion)
        comps.append((-1) ** i * (sub.det(method='berkowitz') if m > 2 else sub[0, 0]))
        if _time.time() - t0 > timeLimit:
            return {'ok': False, 'reason': 'time limit'}
    ia = support.index(anchor)
    if comps[ia] == 0:
        return {'ok': False, 'reason': 'anchor component vanishes'}
    vec = {}
    for i in range(m):
        e = spy.together(comps[i] / comps[ia])
        # cancel() and factor() only while they stay cheap; the caller verifies the
        # entry modulo a fresh prime either way
        if len(str(e)) <= maxChars // 40:
            e = spy.factor(spy.cancel(e))
        if len(str(e)) > maxChars:
            return {'ok': False, 'reason': 'closed form too large'}
        if e != 0:
            vec[support[i]] = str(e)
        if _time.time() - t0 > timeLimit:
            return {'ok': False, 'reason': 'time limit'}
    return {'ok': True, 'vector': vec, 'order': int(order)}


def _jet_expand_chains(conds):
    """Chains of segments as independent jet conditions at gap order 0: segment k starts
    from segment k-1's initial values with its own boundary events applied. Returns the
    first segments and the later ones as two lists."""
    out, extra = [], []
    for c in _as_list(conds):
        c = dict(c)
        if 'segments' not in c:
            out.append(c)
            continue
        segs = [dict(sg) for sg in _as_list(c['segments'])]
        if not segs:
            continue
        ic = [str(l) for l in _as_list(segs[0]['ic'])]
        out.append({'f': segs[0]['f'], 'g': segs[0]['g'], 'ic': ic})
        for sg in segs[1:]:
            if not sg.get('timeFixed', True):
                break
            icd = dict(l.split(' = ', 1) for l in ic)
            for ev in _as_list(sg.get('ev', [])):
                ev = dict(ev)
                X, v, how = str(ev['var']), str(ev['value']), str(ev['method'])
                old = icd.get(X, X)
                icd[X] = ('(%s)' % v if how == 'replace' else
                          '(%s) + (%s)' % (old, v) if how == 'add' else
                          '(%s)*(%s)' % (old, v))
            ic = ['%s = %s' % (X, e) for X, e in icd.items()]
            extra.append({'f': sg['f'], 'g': sg['g'], 'ic': ic})
    return out, extra


def _lie_reach(f_c, g_c, S, ic_c, evs, f_ss=None):
    """Structural first Lie order at which each symbol can enter the output jets of one
    segment, a lower bound for the order its column first becomes nonzero: 0 for the
    symbols of the observables, else one more than the least order of a state whose
    right-hand side holds it (L_f^{k+1} h = sum_i dL_f^k h/dx_i f_i). An initial value
    or an event value enters with its state; at an implicit steady state (f_ss) every
    parameter of the balance may enter with the earliest state. Returns {name: order}."""
    import collections
    Sstr = [str(X) for X in S]
    rhs = {}
    for i, X in enumerate(Sstr):
        try:
            rhs[X] = {str(s) for s in spy.sympify(f_c[i]).free_symbols}
        except Exception:
            rhs[X] = set()
    dist = {}
    queue = collections.deque()
    for e in g_c:
        for s in spy.sympify(e).free_symbols:
            if str(s) not in dist:
                dist[str(s)] = 0
                queue.append(str(s))
    while queue:
        u = queue.popleft()
        for s in rhs.get(u, ()):
            if s not in dist:
                dist[s] = dist[u] + 1
                queue.append(s)
    out = dict(dist)

    def lower(name, k):
        if k < out.get(name, 10 ** 9):
            out[name] = k

    for X, e in dict(ic_c).items():
        if str(X) in dist:
            for s in spy.sympify(e).free_symbols:
                lower(str(s), dist[str(X)])
    for ev in evs or []:
        X = str(ev['var'])
        if X in dist:
            for s in spy.sympify(ev['value']).free_symbols:
                lower(str(s), dist[X])
    if f_ss is not None:
        reached = [dist[X] for X in Sstr if X in dist]
        if reached:
            k0 = min(reached)
            for e in f_ss:
                for s in spy.sympify(e).free_symbols:
                    lower(str(s), k0)
    return out


def compileObservabilityTapeMulti(model, observation, conditionSubs, conditionIC0,
                                  fixed=None, parameters=None, backend='sympy',
                                  equilibrate=False, forcings=None,
                                  segEquilibrate=None, conditionEvents=None,
                                  conditionT0Events=None, jointSteadyState=False,
                                  jointFixedStates=None, heldStateParams=None,
                                  conditionObs=None, conditionTimes=None,
                                  keepCoords=None):
    """Compile one observability tape per experimental condition over a shared
    coordinate space, for the multi-condition observability path.

    `conditionObs` optionally replaces `observation` per condition, so an observable
    measured in some conditions enters only their rows. `conditionSubs` holds one
    {symbol: replacement} map per condition: a number bakes the symbol, a symbol renames
    it (e.g. a knockdown rate). `conditionIC0` holds one {state: expression} map per
    condition with the t0+ initial values composed in R (steady state, `initial`, t0
    events); a state without entry starts free.

    Returns the per-condition tapes (each with an IC tape seeding the initial values and
    their duals), the shared leaf/state layout and the dual-carrying coordinates z (free
    initial values + free parameters, minus `fixed`). A non-rational right-hand side,
    observable or initial value returns {'ok': False, 'nonrational': ...}."""
    _select_backend(backend)
    model = _as_list(model)
    observation = _as_list(observation)
    conditionSubs = conditionSubs or []
    conditionIC0 = conditionIC0 or []
    forcings = set(_as_list(forcings))
    K = len(conditionSubs)
    # one flag per (condition, segment): equilibrated first segments are seeded from
    # the steady state, all others from an IC tape (possibly over carry coordinates)
    if segEquilibrate is None:
        segEq = [bool(equilibrate)] * K
    else:
        segEq = [bool(x) for x in list(segEquilibrate)]
        segEq += [bool(equilibrate)] * (K - len(segEq))

    # missing/None entries of `conditionObs` fall back to `observation`
    conditionObs = list(conditionObs) if conditionObs else []
    obsPerCond = [_as_list(conditionObs[c]) if (c < len(conditionObs) and
                                                conditionObs[c] is not None)
                  else observation for c in range(K)]

    extra = []
    for d in conditionSubs:
        for k, v in dict(d).items():
            extra += [str(k), str(v)]
    for d in conditionIC0:
        for k, v in dict(d).items():
            extra += [str(k), str(v)]
    all_lines = (model + observation + [l for o in obsPerCond for l in o] +
                 _as_list(parameters) + extra)
    # symbol collection is order-insensitive and conditions repeat the same
    # trafo/IC strings, so tokenize each distinct line once
    all_lines = list(dict.fromkeys(all_lines))
    local, parse = _make_local_parse(all_lines)

    variables, diffEquations, _ = _read_equations(model, parse)
    obsRead = [_read_equations(o, parse) for o in obsPerCond]
    obsVarsPer = [r[0] for r in obsRead]
    fixedNames = set(str(s) for s in
                     [local.get(nm, spy.Symbol(nm)) for nm in _as_list(fixed)])
    # held-variable moieties (equilibrate + reduceCQ = FALSE): each pivot state's resting
    # value is a shared initial-value parameter (heldMap), so in joint mode the moiety
    # freedom lands on a reported parameter, not a per-condition state column
    heldStateParams = dict(heldStateParams or {})
    heldMap = {str(k): spy.Symbol(str(v)) for k, v in heldStateParams.items()}
    heldParamSyms = [spy.Symbol(str(v)) for v in heldStateParams.values()]

    isConst = [spy.sympify(e) == 0 for e in diffEquations]
    S = [variables[i] for i in range(len(variables)) if not isConst[i]]
    Srhs = [spy.sympify(diffEquations[i]) for i in range(len(variables))
            if not isConst[i]]
    nS = len(S)

    parseCache = {}

    def pval(s):
        key = str(s)
        e = parseCache.get(key)
        if e is None:
            e = spy.sympify(parse(key))
            parseCache[key] = e
        return e

    # substitution memoised on the transitively relevant subset of the map (conditions
    # repeat the same replacements); exact, since other items can never fire and sympy
    # orders dict items canonically per item
    subsCache = {}

    def subsMemo(e, smap):
        e = spy.sympify(e)
        rel, frontier = {}, set(e.free_symbols)
        while True:
            add = {k: v for k, v in smap.items() if k not in rel and k in frontier}
            if not add:
                break
            rel.update(add)
            for v in add.values():
                frontier |= v.free_symbols
        if not rel:
            return e
        key = (e, tuple(sorted(rel.items(),
                               key=lambda kv: spy.default_sort_key(kv[0]))))
        r = subsCache.get(key)
        if r is None:
            r = e.subs(rel)
            subsCache[key] = r
        return r

    perCond = []
    subsPer = []            # each condition's substitution map, aligned with perCond
    for c in range(K):
        subsMap = {}
        for k, v in dict(conditionSubs[c]).items():
            subsMap[local.get(str(k), spy.Symbol(str(k)))] = pval(v)
        ic0 = dict(conditionIC0[c]) if c < len(conditionIC0) else {}
        f_c = [subsMemo(e, subsMap) for e in Srhs]
        g_c = [subsMemo(e, subsMap) for e in obsRead[c][1]]
        # resting-state model for the equilibrate solve: forcings at 0, no events,
        # non-forcing per-condition substitutions baked in
        forcZero = {local.get(nm, spy.Symbol(nm)): spy.Integer(0) for nm in forcings}
        subsMapNF = {k: v for k, v in subsMap.items() if str(k) not in forcings}
        f_ss = [subsMemo(subsMemo(e, subsMapNF), forcZero) for e in Srhs]
        ic_c = {}
        for X in S:
            e = pval(ic0[str(X)]) if str(X) in ic0 else X
            ic_c[str(X)] = subsMemo(e, subsMap)
        perCond.append((f_c, g_c, ic_c, f_ss))
        subsPer.append(subsMap)

    # a logarithmic observable carries the information of its argument; hyperbolic
    # functions are exponentials
    perCond = [([_hyp_to_exp(e) for e in f_c],
                [_hyp_to_exp(_strip_log_obs(e)) for e in g_c],
                {k: _hyp_to_exp(v) for k, v in ic_c.items()},
                [_hyp_to_exp(e) for e in f_ss])
               for (f_c, g_c, ic_c, f_ss) in perCond]

    # log-parametrised parameters: one that occurs only as base^(c*theta) is replaced
    # by the rational coordinate X = base^theta, and the report maps X back to theta
    conditionEvents = conditionEvents or []
    conditionTimes = conditionTimes or []
    evExprs = []
    for c in range(K):
        for e in (list(conditionEvents[c]) if c < len(conditionEvents) else []):
            evExprs.append(subsMemo(pval(e['value']), subsPer[c]))
        if c < len(conditionTimes) and conditionTimes[c] is not None:
            evExprs.append(subsMemo(pval(conditionTimes[c]), subsPer[c]))
    lpExprs = [e for (f_c, g_c, ic_c, f_ss) in perCond
               for e in list(f_c) + list(g_c) + list(ic_c.values()) + list(f_ss)]
    lpBase = _detect_log_params(lpExprs + evExprs, set(S))
    lpBase = _log_params_only_in_atoms(lpBase, lpExprs + evExprs)
    taken = set()
    for e in lpExprs + evExprs:
        taken |= {str(x) for x in spy.sympify(e).free_symbols}
    lp = {}
    for t, b in lpBase.items():
        nm = 'exp_%s' % t
        while nm in taken:
            nm += '_'
        taken.add(nm)
        lp[t] = (spy.Symbol(nm), b)
    lpSub = _log_param_sub(lp)
    if lp:
        perCond = [([lpSub(e) for e in f_c], [lpSub(e) for e in g_c],
                    {k: lpSub(v) for k, v in ic_c.items()}, [lpSub(e) for e in f_ss])
                   for (f_c, g_c, ic_c, f_ss) in perCond]
        for t, (X, b) in lp.items():
            if str(t) in fixedNames:
                fixedNames.add(str(X))
    logParams = [{'X': str(X), 'theta': str(t), 'base': 'E' if b == spy.E else str(b)}
                 for t, (X, b) in sorted(lp.items(), key=lambda kv: str(kv[0]))]

    # each condition's later events and segment start time, in the leaves
    evPer = [[{'var': str(e['var']), 'method': str(e['method']),
               'value': lpSub(subsMemo(pval(e['value']), subsPer[c]))}
              for e in (list(conditionEvents[c]) if c < len(conditionEvents) else [])]
             for c in range(K)]
    tmPer = [lpSub(subsMemo(pval(conditionTimes[c]), subsPer[c]))
             if c < len(conditionTimes) and conditionTimes[c] is not None else None
             for c in range(K)]

    # power/Hill recast: base^exp (exp a parameter) becomes a state E with
    # E' = exp*E*base'/base plus L = log(base), tied to (base, exp) downstream, so f
    # stays rational (sound: base, log(base), base^exp are algebraically independent).
    # E is held generic in the f = 0 solve, or a free-initial-value leaf when transient.
    nReal = nS
    powerRecast = []
    invSolveName = {}
    pairs = _detect_power_atoms(perCond, S)
    if pairs is None:
        return {'ok': False, 'nonrational':
                ['unsupported power form: base must be a symbol and exponent c*param']}
    if pairs:
        S, perCond, powerRecast = _apply_power_recast(pairs, S, perCond)
        nS = len(S)
        if equilibrate:
            # a base with a linear turnover term is solved directly; otherwise its
            # balance is linear in E and it is "inverted": E is solved, base and L
            # stay generic
            realSet = {str(X) for X in S[:nReal]}
            bareSyms = set()
            for (_f, _g, _ic, ss_r) in perCond:
                for e in ss_r[:nReal]:
                    bareSyms |= spy.fraction(spy.together(e))[0].free_symbols
            for rc in powerRecast:
                rc['inverted'] = (rc['base'] in realSet
                                  and spy.Symbol(rc['base']) not in bareSyms)
            byBase = {}
            for rc in powerRecast:
                byBase.setdefault(rc['base'], []).append(rc)
            for base, group in byBase.items():
                if group[0]['inverted'] and len(group) > 1:
                    return {'ok': False, 'nonrational':
                            ['free exponent base %s carries multiple independent '
                             'powers with no linear turnover term' % base]}
            invSolveName = {rc['base']: rc['E']
                            for rc in powerRecast if rc['inverted']}
        else:
            # transient: no steady-state solve, so nothing is inverted
            for rc in powerRecast:
                rc['inverted'] = False

    def anyIn(pred):
        return (any(pred(e) for (f_c, g_c, ic_c, f_ss) in perCond
                    for e in list(f_c) + list(g_c) + list(ic_c.values())) or
                any(pred(e['value']) for evs in evPer for e in evs) or
                any(t is not None and pred(t) for t in tmPer))

    def takenNames():
        taken = {str(x) for x in S} | set(fixedNames) | {str(X) for X, _ in lp.values()}
        for (f_c, g_c, ic_c, f_ss) in perCond:
            for e in list(f_c) + list(g_c) + list(ic_c.values()) + list(f_ss):
                taken |= {str(x) for x in spy.sympify(e).free_symbols}
        for e in lpExprs + evExprs:
            taken |= {str(x) for x in spy.sympify(e).free_symbols}
        return taken

    # trigonometric functions: half-angle states and leaves (_apply_trig_recast), before
    # the exponentials so that an exponent may contain them
    scalS, scalPerCond = list(S), perCond
    trigRecast = None
    if anyIn(_has_trig):
        if equilibrate:
            return {'ok': False, 'why': 'equilibrate = TRUE does not support '
                    'trigonometric functions; give the steady state through `trafo` '
                    'or start from free initial values'}
        trigRecast = _apply_trig_recast(S, perCond, evPer, tmPer, takenNames())
        if 'why' in trigRecast:
            return {'ok': False, 'why': trigRecast['why']}
        S, perCond = trigRecast['S'], trigRecast['perCond']
        evPer, tmPer = trigRecast['evPer'], trigRecast['tmPer']
        nS = len(S)
        fixedNames |= set(trigRecast['consts'])

    # exponentials: auxiliary states and generic exponential leaves (_apply_exp_recast)
    expRecast = None
    if anyIn(_has_exp):
        if equilibrate:
            return {'ok': False, 'why': 'equilibrate = TRUE does not support exp() '
                    'or b^x of a state, or of a parameter that also enters elsewhere; '
                    'give the steady state through `trafo` or start from free '
                    'initial values'}
        expRecast = _apply_exp_recast(S, perCond, evPer, tmPer, takenNames())
        if 'why' in expRecast:
            return {'ok': False, 'why': expRecast['why']}
        S, perCond = expRecast['S'], expRecast['perCond']
        evPer, tmPer = expRecast['evPer'], expRecast['tmPer']
        nS = len(S)
        fixedNames |= set(expRecast['consts'])
    if trigRecast is not None:
        expRecast = _merge_recast(trigRecast, expRecast)

    nonrational = []
    ratOK = set()            # deduped exprs recur across conditions; check each once

    def _rat(e):
        if e in ratOK:
            return True
        ok = _is_rational_expr(e)
        if ok:
            ratOK.add(e)
        return ok

    for c, (f_c, g_c, ic_c, f_ss) in enumerate(perCond):
        for X, e in zip(S, f_c):
            if not _rat(e):
                nonrational.append('d%s/dt = %s' % (X, e))
        for y, e in zip(obsVarsPer[c], g_c):
            if not _rat(e):
                nonrational.append('%s = %s' % (y, e))
        for X in S:
            if not _rat(ic_c[str(X)]):
                nonrational.append('%s(0) = %s' % (X, ic_c[str(X)]))
    if nonrational:
        return {'ok': False, 'nonrational': nonrational}

    # a state has a free initial-value leaf when its initial value still contains the
    # state symbol (possibly under a dose); equilibrate-seeded states have none
    freeState = {str(X): False for X in S}
    for c, (f_c, g_c, ic_c, f_ss) in enumerate(perCond):
        if segEq[c]:
            # joint mode: states stay free coordinates (seeded to x* by R, constraint
            # as stacked df rows), keeping directions low-degree in (x, theta)
            if jointSteadyState:
                jfs = set(jointFixedStates or [])
                for X in S:
                    if (str(X) not in forcings and str(X) not in jfs
                            and str(X) not in heldMap):
                        freeState[str(X)] = True
            continue  # equilibrate-seeded: no state is a free coordinate (unless joint)
        for X in S:
            if X in spy.sympify(ic_c[str(X)]).free_symbols:
                freeState[str(X)] = True

    # boundary state-dose values (e.g. a second-dose parameter) may appear only in
    # an event, never in f/g/ic; collect their symbols so they become coordinates
    conditionEvents = conditionEvents or []
    eventVals = []
    for c in range(K):
        evs = list(conditionEvents[c]) if c < len(conditionEvents) else []
        eventVals.append([lpSub(spy.sympify(pval(e['value']))) for e in evs])

    # f_ss is included so parameters the perturbed dynamics drop still count
    paramset = set()
    seenPS = set()           # deduped exprs recur across conditions; walk each once
    for (f_c, g_c, ic_c, f_ss) in perCond:
        for e in list(f_c) + list(g_c) + list(ic_c.values()) + list(f_ss):
            if e in seenPS:
                continue
            seenPS.add(e)
            paramset |= set(spy.sympify(e).free_symbols)
    for vs in eventVals:
        for e in vs:
            paramset |= set(e.free_symbols)
    for evs in evPer:
        for e in evs:
            paramset |= set(spy.sympify(e['value']).free_symbols)
    if expRecast is not None:
        byName = {str(x): x for x in paramset}
        for rl in expRecast['rel']:
            paramset |= set(rl[2].free_symbols)
            paramset.add(byName.get(rl[1], spy.Symbol(rl[1])))
    # a segment's left boundary may sit at a time given in the parameters
    for c, tm in enumerate(conditionTimes):
        if tm is not None and c < K:
            paramset |= set(lpSub(pval(tm)).free_symbols)
    paramset -= set(S)
    paramset |= set(heldParamSyms)   # held-variable initial-value parameters
    # coordinates a chart absorbed stay coordinates
    paramset |= {local.get(str(n), spy.Symbol(str(n))) for n in _as_list(keepCoords)}
    params = sorted(paramset, key=spy.default_sort_key)
    # a base seen only under a free exponent (Km in C^n/(Km^n+C^n)) vanished into its
    # recast atom; re-add it so the recast relation ties it and the Km co-scaling is
    # reported as an exact scaling
    if powerRecast:
        known = set(str(s) for s in S) | set(str(s) for s in params) | set(forcings)
        for rc in powerRecast:
            if rc['base'] not in known:
                params.append(spy.Symbol(rc['base'])); known.add(rc['base'])
        params = sorted(params, key=spy.default_sort_key)

    freeStates = [X for X in S if freeState[str(X)]]
    leafNames = [str(X) for X in freeStates] + [str(s) for s in params]
    leafSlot = {nm: i for i, nm in enumerate(leafNames)}
    nLeaves = len(leafNames)
    slotOf = dict(leafSlot)
    for i in range(nS):
        slotOf[str(S[i])] = nLeaves + i
    base = nLeaves + nS

    tapes = []
    reachCache = {}  # structural Lie reach per (model, initial values, events)
    emitCache = {}   # (f_c, g_c, obs names) recur across segments/conditions
    ssSolCache = {}  # linear resting-state solutions, keyed on the f_ss tuple
    icCache = {}     # emitted IC tapes, keyed on the seed-expression tuple
    for c, (f_c, g_c, ic_c, f_ss) in enumerate(perCond):
        subsMap = subsPer[c]
        ekey = (tuple(f_c), tuple(g_c), tuple(str(v) for v in obsVarsPer[c]))
        cached = emitCache.get(ekey)
        if cached is None:
            try:
                emitted = _emit_tape_shared(f_c, g_c, slotOf, base)
            except _NotRational:
                return {'ok': False}
            # this segment's substituted model as lines for the scaling peel,
            # exponentials kept
            f0, g0 = scalPerCond[c][0], scalPerCond[c][1]
            eNum = {spy.E: spy.exp(1, evaluate=False)}
            mLines = ['%s = %s' % (str(scalS[i]), spy.sympify(f0[i]).xreplace(eNum))
                      for i in range(len(scalS))]
            oLines = ['%s = %s' % (str(obsVarsPer[c][j]), spy.sympify(g0[j]).xreplace(eNum))
                      for j in range(len(g0))]
            cached = (emitted, mLines, oLines)
            emitCache[ekey] = cached
        (op, a, b, cnum, cden, outslots), mLines, oLines = cached
        fOut = outslots[:nS]
        # structural first Lie order of every symbol in this segment's jets; the R
        # saturation counts no flat step before it (a readout behind a transit chain is
        # flat at rank 0 until the chain has filled)
        reachKey = (ekey, tuple(sorted((str(k), str(v)) for k, v in ic_c.items())),
                    tuple((str(e['var']), str(e['value'])) for e in evPer[c]),
                    bool(segEq[c]))
        reach = reachCache.get(reachKey)
        if reach is None:
            reach = _lie_reach(f_c, g_c, S, ic_c, evPer[c],
                               f_ss if segEq[c] else None)
            reachCache[reachKey] = reach
        tape = {
            'op': op, 'a': a, 'b': b, 'cnum': cnum, 'cden': cden,
            'stateSlots': [nLeaves + i for i in range(nS)],
            'fOut': fOut,
            'gOut': outslots[nS:nS + len(g_c)],
            'icLeaf': [-1] * nS, 'icNum': ['0'] * nS, 'icDen': ['1'] * nS,
            'modelLines': mLines,
            'obsLines': oLines,
            # initial values, boundary events and start time, for verifyScalings()
            'verIC': ['%s = %s' % (str(X), str(ic_c[str(X)])) for X in S if str(X) in ic_c],
            'verEv': [{'var': str(e['var']), 'method': str(e['method']),
                       'value': str(e['value'])} for e in evPer[c]],
            'verT0': [{'var': str(e['var']), 'method': str(e['method']),
                       'value': str(e['value'])} for e in
                      (conditionT0Events[c] if conditionT0Events and
                       c < len(conditionT0Events) else [])],
            'verTime': '' if tmPer[c] is None else str(tmPer[c]),
        }
        if segEq[c] and jointSteadyState:
            # joint mode: each non-forcing state is an identity leaf that R seeds to
            # x*; the resting model is kept so R can solve x* and read [Jx|Jt]
            jfs = set(jointFixedStates or [])
            # held pivots seed from their parameter (heldMap), forcings and forced-zero
            # states from 0
            icMap = {X: (heldMap[str(X)] if str(X) in heldMap
                         else X if (str(X) not in forcings and str(X) not in jfs)
                         else spy.Integer(0)) for X in S}
            # t0 events act on the pre-event coordinate x_ss (seeded by R from valBy);
            # the IC tape carries x0 = E(x_ss) with its chain-rule duals
            _compose_t0_events(icMap, conditionT0Events[c]
                               if conditionT0Events and c < len(conditionT0Events)
                               else [], pval, subsMap)
            try:
                icOp, icA, icB, icCnum, icCden, icOut = _emit_tape_shared(
                    [icMap[X] for X in S], [], leafSlot, nLeaves)
            except _NotRational:
                return {'ok': False}
            tape.update({'icOp': icOp, 'icA': icA, 'icB': icB, 'icCnum': icCnum,
                         'icCden': icCden, 'icOut': icOut})
            tape['constraintModel'] = [
                '%s = %s' % (invSolveName.get(str(S[i]), str(S[i])),
                             spy.sympify(f_ss[i])) for i in range(nReal)]
        elif segEq[c]:
            # equilibrate-seeded first segment: a generically linear resting state (no
            # recast) is solved symbolically and emitted as an IC tape; otherwise R
            # seeds it numerically per point with its IFT duals (icSeed)
            ssIC = None
            if not powerRecast and not _FORCE_CONSTRAINT_SEED:
                sskey = tuple(f_ss)
                if sskey in ssSolCache:
                    ssSol = ssSolCache[sskey]
                else:
                    solveS = [S[i] for i in range(nReal)
                              if str(S[i]) not in forcings]
                    fssPolys = [spy.expand(spy.fraction(spy.together(
                                    spy.sympify(f_ss[i])))[0])
                                for i in range(nReal) if str(S[i]) not in forcings]
                    ssSol = _linear_solution(fssPolys, solveS)
                    ssSolCache[sskey] = ssSol
                if ssSol is not None:
                    # recast coordinates (beyond nReal) stay generic; forcings and
                    # dead states stay 0
                    icMap = {X: (ssSol[X] if X in ssSol
                                 else X if i >= nReal else spy.Integer(0))
                             for i, X in enumerate(S)}
                    # t0 events compose onto the resting state
                    _compose_t0_events(icMap, conditionT0Events[c]
                                       if conditionT0Events and
                                       c < len(conditionT0Events) else [],
                                       pval, subsMap)
                    try:
                        icOp, icA, icB, icCnum, icCden, icOut = _emit_tape_shared(
                            [icMap[X] for X in S], [], leafSlot, nLeaves)
                        ssIC = {'icOp': icOp, 'icA': icA, 'icB': icB,
                                'icCnum': icCnum, 'icCden': icCden, 'icOut': icOut}
                    except _NotRational:
                        ssIC = None
            if ssIC is not None:
                tape.update(ssIC)
            else:
                tape['constraintModel'] = [
                    '%s = %s' % (invSolveName.get(str(S[i]), str(S[i])),
                                 spy.sympify(f_ss[i])) for i in range(nReal)]
        else:
            # IC tape: seed each state from its initial-value expression in the
            # leaves (free initial values, carry coordinates, doses, resets).
            ickey = tuple(ic_c[str(X)] for X in S)
            icEmitted = icCache.get(ickey)
            if icEmitted is None:
                try:
                    icEmitted = _emit_tape_shared(list(ickey), [], leafSlot,
                                                  nLeaves)
                except _NotRational:
                    return {'ok': False}
                icCache[ickey] = icEmitted
            icOp, icA, icB, icCnum, icCden, icOut = icEmitted
            tape.update({'icOp': icOp, 'icA': icA, 'icB': icB,
                         'icCnum': icCnum, 'icCden': icCden, 'icOut': icOut})
            # the symbolic jets of this segment, for jetGenerator(): dynamics,
            # observables and initial values over the leaves
            tape['jetF'] = ['%s = %s' % (str(S[i]), str(f_c[i])) for i in range(nS)]
            tape['jetG'] = [str(e) for e in g_c]
            tape['jetIC'] = ['%s = %s' % (str(X), str(ic_c[str(X)])) for X in S]
            # a later segment's boundary events, and whether its start time is free of
            # the coordinates (the jets of a moving boundary need the time shift too)
            tape['jetEv'] = [{'var': str(e['var']), 'method': str(e['method']),
                              'value': str(e['value'])} for e in evPer[c]]
            tape['jetTimeFixed'] = bool(tmPer[c] is None or
                                        not spy.sympify(tmPer[c]).free_symbols)
        # events at this segment's left boundary, applied by the kernel to the
        # propagated state; their values as an order-0 tape over the leaves
        evs = evPer[c]
        if evs:
            idxOfState = {str(X): i for i, X in enumerate(S)}
            methodCode = {'replace': 0, 'add': 1, 'multiply': 2}
            keep = [e for e in evs if str(e['var']) in idxOfState]
            if keep:
                vals = [e['value'] for e in keep]
                try:
                    evOp, evA, evB, evCnum, evCden, evO = _emit_tape_shared(
                        [], vals, leafSlot, nLeaves)
                except _NotRational:
                    return {'ok': False}
                tape.update({
                    'evOp': evOp, 'evA': evA, 'evB': evB, 'evCnum': evCnum,
                    'evCden': evCden, 'evOut': evO,
                    'evVarIdx': [idxOfState[str(e['var'])] for e in keep],
                    'evMethod': [methodCode[str(e['method'])] for e in keep]})
        # the segment's left boundary time as an order-0 tape: its value sets the gap
        # lengths, its duals carry a time that depends on the coordinates
        if tmPer[c] is not None:
            tv = tmPer[c]
            try:
                tmOp, tmA, tmB, tmCnum, tmCden, tmO = _emit_tape_shared(
                    [], [tv], leafSlot, nLeaves)
            except _NotRational:
                return {'ok': False}
            tape.update({'tmOp': tmOp, 'tmA': tmA, 'tmB': tmB, 'tmCnum': tmCnum,
                         'tmCden': tmCden, 'tmOut': tmO})
        tape['_reach'] = reach
        tapes.append(tape)

    zStateNames = [str(X) for X in freeStates if str(X) not in fixedNames]
    zParamNames = [str(s) for s in params if str(s) not in fixedNames]
    znames = zStateNames + zParamNames
    zSlots = [leafSlot[nm] for nm in znames]
    # per tape, the first Lie order at which each coordinate can enter (-1: never)
    for t in tapes:
        r = t.pop('_reach', None) or {}
        t['lieReach'] = [int(r.get(nm, -1)) for nm in znames]

    out = {
        'ok': True,
        'tapes': tapes,
        'nLeaves': nLeaves,
        'nStates': nS,
        'zSlots': zSlots,
        'znames': znames,
        'leafNames': leafNames,
        'logParams': logParams,
    }
    if expRecast is not None:
        # relation rows dW = W log(r) d(tau/L): one tape output per (W, leaf) pair
        rel = expRecast['rel']
        try:
            rOp, rA, rB, rCnum, rCden, rOut = _emit_tape_shared(
                [rl[2] for rl in rel], [], leafSlot, nLeaves)
        except _NotRational:
            return {'ok': False, 'why': 'an exponent that is not rational is not '
                    'supported by symEngine = "modular"; try symEngine = "symbolic"'}
        out['expRelation'] = {'W': [rl[0] for rl in rel], 'z': [rl[1] for rl in rel],
                              'op': rOp, 'a': rA, 'b': rB, 'cnum': rCnum,
                              'cden': rCden, 'out': rOut}
        out['expAtoms'] = expRecast['atoms']
        out['expBack'] = {'names': list(expRecast['back'].keys()),
                          'values': list(expRecast['back'].values())}
        out['expCodim'] = expRecast['codim']
    if equilibrate:
        out['equilibrate'] = True
        out['stateNames'] = [str(X) for X in S]
        out['paramNames'] = [str(s) for s in params]
        out['forcings'] = sorted(forcings)
        out['realStateNames'] = [invSolveName.get(str(X), str(X))
                                 for X in S[:nReal]]
        out['powerRecast'] = powerRecast
        out['heldStateParams'] = heldStateParams
        if jointSteadyState:
            out['jointSteadyState'] = True
            # state z-columns (the free states that entered znames), for the R joint
            # branch to seed on-manifold and stack the df constraint rows
            out['zStateNames'] = zStateNames
    elif powerRecast:
        # transient free exponent: E and L are free-initial-value leaves; R stacks the
        # recast relation rows and reports in real states + parameters
        out['recastTransient'] = True
        out['stateNames'] = [str(X) for X in S]
        out['paramNames'] = [str(s) for s in params]
        out['forcings'] = sorted(forcings)
        out['realStateNames'] = [str(X) for X in S[:nReal]]
        out['powerRecast'] = powerRecast
        # the recast atom coordinates (E and L), auxiliary to the physical report
        out['recastAtomNames'] = ([rc['E'] for rc in powerRecast] +
                                  sorted(set(rc['L'] for rc in powerRecast)))
    return out


###########################################################################
#####################     scaling symmetries     ########################
###########################################################################

# ---- scaling (toric) symmetries from the integer kernel ------------------------------

def _poly_monomials(expr, zvars):
    """Numerator and denominator monomial-exponent lists of expr over zvars
    (other symbols are treated as weight-zero coefficients). Exponents may be
    rational, as in sqrt(x). Raises if expr is not a ratio of such sums.
    Read term by term: a dense Poly nests one level per coordinate."""
    idx = {z: i for i, z in enumerate(zvars)}
    zset = set(zvars)
    zero = spy.Integer(0)

    zpos = {str(z): i for i, z in enumerate(zvars)}

    def mons_se(poly):
        # symengine expands a product of sums tens of times faster than sympy; the
        # terms are read the same way, a coordinate inside anything but a power of
        # itself is refused
        import symengine as _se
        ep = _se.expand(_se.sympify(poly))
        terms = ep.args if isinstance(ep, _se.Add) else (ep,)
        out = {}
        for t in terms:
            if t == 0:
                continue
            a = [zero] * len(zvars)
            for k, v in t.as_powers_dict().items():
                ks = str(k)
                inExp = any(str(fs) in zpos for fs in getattr(v, 'free_symbols', ()))
                if ks in zpos and not inExp:
                    a[zpos[ks]] = spy.Rational(str(v))
                elif inExp or any(str(fs) in zpos
                                  for fs in getattr(k, 'free_symbols', ())):
                    # a coordinate in an exponent (exp(x) is E**x here) or inside a
                    # function
                    raise ValueError('not a monomial in the coordinates: %s' % t)
            out[tuple(a)] = None
        return list(out)

    def mons(poly):
        try:
            return mons_se(poly)
        except ValueError:
            raise
        except Exception:
            pass
        out = {}
        for t in spy.Add.make_args(spy.expand(poly)):
            if t == 0:
                continue
            a = [zero] * len(zvars)
            for k, v in t.as_powers_dict().items():
                if k in idx:
                    a[idx[k]] = spy.Rational(v)
                elif k.free_symbols & zset:
                    raise ValueError('not a monomial in the coordinates: %s' % t)
            out[tuple(a)] = None
        return list(out)
    p, q = spy.fraction(spy.together(spy.sympify(expr)))
    return mons(p), mons(q)


def _exp_split(expr, logs=False):
    """expr with every exponential atom (exp(u), or b^u with a numeric base), every
    trigonometric function of u, and with `logs` every log(u), replaced by a fresh
    symbol, innermost first, and the list of the u. Each is invariant exactly when u
    is: log(lambda^w*u) = log(u) + w*log(lambda)."""
    exps = []

    def rec(e):
        if not e.args:
            return e
        be = _exp_atom(e)
        if be is not None:
            exps.append(rec(be[1]))
            return spy.Dummy('exp')
        if isinstance(e, _TRIG) or type(e) in _TRIG_RECIPROCAL:
            exps.append(rec(e.args[0]))
            return spy.Dummy('trig')
        if logs and isinstance(e, spy.log):
            exps.append(rec(e.args[0]))
            return spy.Dummy('log')
        return e.func(*[rec(a) for a in e.args])
    return rec(_hyp_to_exp(expr)), exps


def _scaling_rows(diffEquations, obsFunctions, m, zvars, interOffset, logs=False,
                  seen=None):
    """Sparse {col: coeff} monomial-exponent rows of one (f, g) system for the scaling
    kernel: weight columns 0..nz-1 over zvars, intermediate columns from interOffset.
    Returns (rows, ninter, skipped). An exponential has weight zero and its exponent
    is an invariant, like an observable; with `logs` so is the argument of a log(),
    off for the power recast, which keeps log(base) as a coordinate. An (expression,
    target) pair already in `seen` is skipped."""
    nz = len(zvars)
    exprs = []          # (numer monomials, denom monomials, target weight vector)
    skipped = 0
    split = [_exp_split(e, logs) for e in list(obsFunctions) + list(diffEquations[:m])]
    obsFunctions = ([s[0] for s in split[:len(obsFunctions)]] +
                    [u for s in split for u in s[1]])
    diffEquations = [s[0] for s in split[len(split) - m:]] if m else []
    seen = set() if seen is None else seen

    def fresh(e, t):
        key = (e, t)
        if key in seen:
            return False
        seen.add(key)
        return True

    for g in obsFunctions:
        if not fresh(g, -1):
            continue
        try:
            pmon, qmon = _poly_monomials(g, zvars)
        except Exception:
            skipped += 1
            continue
        exprs.append((pmon, qmon, [0] * nz))
    for i in range(m):
        if not fresh(diffEquations[i], i):
            continue
        try:
            pmon, qmon = _poly_monomials(diffEquations[i], zvars)
        except Exception:
            skipped += 1
            continue
        tvec = [0] * nz
        tvec[i] = 1     # state i is column i of zvars; x_i has weight c_i
        exprs.append((pmon, qmon, tvec))

    def diffRow(a, b, t=None):
        # (a - b - t) . c = 0 over the weight columns, times the lcm of the denominators
        d = [x - y for x, y in zip(a, b)]
        if t is not None:
            d = [x - int(y) for x, y in zip(d, t)]
        nzj = [(j, x) for j, x in enumerate(d) if x != 0]
        if not nzj:
            return None
        D = 1
        for _, x in nzj:
            q = int(getattr(x, 'q', 1))
            if q != 1:
                D = D * q // math.gcd(D, q)
        return {j: int(x * D) for j, x in nzj}

    # every numerator monomial has the weight of the first, every denominator monomial
    # that of the first denominator one, and the two differ by the target weight: the
    # rows of the intermediate-column form with the intermediates eliminated
    rows, keys = [], set()
    for pmon, qmon, tvec in exprs:
        cand = [diffRow(a, pmon[0]) for a in pmon[1:]] + \
               [diffRow(b, qmon[0]) for b in qmon[1:]]
        if pmon and qmon:
            cand.append(diffRow(pmon[0], qmon[0], tvec))
        for row in cand:
            if row is None:
                continue
            k = tuple(sorted(row.items()))
            if k not in keys:
                keys.add(k)
                rows.append(row)
    return rows, 0, skipped


def _materialize_rows(rows, ncols):
    """Dense matrix (list of lists) from sparse {col: coeff} rows."""
    out = []
    for r in rows:
        dense = [0] * ncols
        for c, v in r.items():
            dense[c] = v
        out.append(dense)
    return out


def _scaling_gens(rows, ncols, nz):
    """Reduced integer generators (primitive, over the first nz weight columns)
    of the scaling kernel given the determining rows."""
    basis = exactNullspace(rows, ncols, integer=True) if rows else []
    cvecs = [[v[j] for j in range(nz)] for v in basis
             if any(v[j] != 0 for j in range(nz))]
    gens = []
    if cvecs:
        red, _ = spy.Matrix(cvecs).rref()
        for r in range(red.rows):
            row = [red[r, j] for j in range(nz)]
            if all(x == 0 for x in row):
                continue
            L = 1
            for x in row:
                L = spy.ilcm(L, spy.Rational(x).q)
            ints = [int(spy.Rational(x) * L) for x in row]
            d = 0
            for x in ints:
                d = spy.igcd(d, x)
            if d:
                ints = [x // d for x in ints]
            gens.append(ints)
    return gens


def _weighted_degree(e, w):
    """Weighted degree of `e` with symbol weights `w` (name -> number, default 0), or None
    when `e` is not weighted-homogeneous."""
    if e.is_Number or e.is_NumberSymbol:
        return 0
    if e.is_Symbol:
        return w.get(str(e), 0)
    if e.is_Add:
        ds = [_weighted_degree(a, w) for a in e.args]
        return ds[0] if None not in ds and len(set(ds)) == 1 else None
    if e.is_Mul:
        tot = 0
        for a in e.args:
            d = _weighted_degree(a, w)
            if d is None:
                return None
            tot += d
        return tot
    if e.is_Pow:
        b, x = e.args
        db = _weighted_degree(b, w)
        if x.is_Number:
            return None if db is None else db * x
        return 0 if db == 0 and _weighted_degree(x, w) == 0 else None
    if e.args:
        return 0 if all(_weighted_degree(a, w) == 0 for a in e.args) else None
    return None


def verifyScalings(scalings, tapes):
    """Exact check that each scaling (name -> weight) leaves every tape's model, observables,
    initial values, events and start time invariant: a right-hand side has the weight of
    its state, an observable weight 0, a replaced or added value the weight of its target,
    a factor and a time weight 0. Returns one bool per scaling."""
    lines = []
    for t in tapes:
        lines += _as_list(t.get('modelLines', [])) + _as_list(t.get('obsLines', []))
        lines += _as_list(t.get('verIC', []))
        lines += ['_ = %s' % e['value'] for e in _as_list(t.get('verEv', [])) +
                  _as_list(t.get('verT0', []))]
        if t.get('verTime'):
            lines.append('_ = %s' % t['verTime'])
    local, parse = _make_local_parse(lines)

    def split(line):
        lhs, rhs = _clean(line).split('=', 1)
        return lhs.strip(), spy.sympify(parse(rhs))

    checks = []                       # (expression, target: state name, None for 0)
    for t in tapes:
        for l in _as_list(t.get('modelLines', [])) + _as_list(t.get('verIC', [])):
            if '=' in l:
                x, e = split(l)
                checks.append((e, x))
        for l in _as_list(t.get('obsLines', [])):
            if '=' in l:
                checks.append((split(l)[1], None))
        for ev in _as_list(t.get('verEv', [])) + _as_list(t.get('verT0', [])):
            e = spy.sympify(parse(str(ev['value'])))
            checks.append((e, None if str(ev['method']) == 'multiply' else str(ev['var'])))
        if t.get('verTime'):
            checks.append((spy.sympify(parse(t['verTime'])), None))
    out = []
    for sc in scalings:
        w = {str(k): spy.Rational(str(v)) for k, v in dict(sc).items()}
        ok = True
        for e, target in checks:
            if e == 0:                # zero has every weight
                continue
            d = _weighted_degree(e, w)
            if d is None or d != (w.get(target, 0) if target else 0):
                ok = False
                if os.environ.get('DMOD_SYM_LIEDIAG'):
                    print('[liediag] scaling fails at %s -> %s: degree %s' %
                          (target, str(e)[:120], d), flush=True)
                break
        out.append(ok)
    return out


# ---- certificate for general directions: an invariant distribution ------------------
# Fields X_j over (states, parameters) with X_j(h) = 0 and [F, X_j] in their span annihilate
# every Lie derivative of h; the flow of F keeps the span. Checked at random points mod p.

class _DistCheck:
    """Rational expressions over the coordinates, evaluated exactly mod p."""

    def __init__(self, lines):
        self.local, self.parse = _make_local_parse(list(lines))
        self.fns = {}

    def expr(self, text):
        return spy.sympify(self.parse(str(text)))

    def value(self, e, syms, vals, p, skey=None):
        """e at the point vals (Fractions aligned with syms) mod p; None on a pole or a
        non-rational value. `skey` names the symbol list in the cache key, which saves
        hashing a long list of symbols on every call."""
        key = (e, tuple(syms) if skey is None else skey)
        fn = self.fns.get(key)
        if fn is None:
            # decimals as exact rationals, rationals as Fraction: no float enters
            ex = e.xreplace({x: spy.Rational(str(x)) for x in e.atoms(spy.Float)})
            frac = spy.Function('Fraction')
            ex = ex.xreplace({x: frac(x.p, x.q) for x in ex.atoms(spy.Rational)
                              if not x.is_Integer})
            fn = spy.lambdify(syms, ex, modules=[{'Fraction': Fraction}])
            self.fns[key] = fn
        try:
            v = fn(*vals)
        except ZeroDivisionError:
            return None
        if isinstance(v, int):
            return v % p
        if isinstance(v, Fraction):
            if v.denominator % p == 0:
                return None
            return v.numerator % p * pow(v.denominator % p, p - 2, p) % p
        return None


def _span_coeffs(B, v, p):
    """Coefficients of v over the vectors B (lists of residues) mod p, earlier vectors
    first and dependent ones at zero, or None if v is not in their span. An incremental
    echelon that carries each row's combination of B; p < 2^31 keeps products in int64."""
    n = len(B)
    rows = []
    for k, b in enumerate(B):
        x = np.asarray(b, dtype=np.int64) % p
        comb = np.zeros(n, dtype=np.int64)
        comb[k] = 1
        for piv, r, cm in rows:
            f = int(x[piv])
            if f:
                x = (x - f * r) % p
                comb = (comb - f * cm) % p
        nz = np.flatnonzero(x)
        if nz.size == 0:
            continue
        piv = int(nz[0])
        inv = pow(int(x[piv]), p - 2, p)
        rows.append((piv, x * inv % p, comb * inv % p))
    y = np.asarray(v, dtype=np.int64) % p
    coef = np.zeros(n, dtype=np.int64)
    for piv, r, cm in rows:
        f = int(y[piv])
        if f:
            y = (y - f * r) % p
            coef = (coef + f * cm) % p
    if y.any():
        return None
    return [int(c) for c in coef]


def _rank_mod(rows, p):
    """Rank of a list of rows (lists of ints) over GF(p)."""
    m = [list(r) for r in rows]
    rank, ncol = 0, (len(m[0]) if m else 0)
    for c in range(ncol):
        piv = next((i for i in range(rank, len(m)) if m[i][c] % p), None)
        if piv is None:
            continue
        m[rank], m[piv] = m[piv], m[rank]
        inv = pow(m[rank][c] % p, p - 2, p)
        m[rank] = [x * inv % p for x in m[rank]]
        for i in range(len(m)):
            if i != rank and m[i][c] % p:
                fct = m[i][c]
                m[i] = [(a - fct * b) % p for a, b in zip(m[i], m[rank])]
        rank += 1
    return rank


def _in_span(basis, v, p):
    return _rank_mod(basis + [v], p) == _rank_mod(basis, p)


def _coeffs(basis, v, p):
    """Coefficients c with sum c_i basis_i = v over GF(p) for independent rows `basis`,
    or None if v is not in their span."""
    n = len(basis)
    if n == 0:
        return [] if not any(x % p for x in v) else None
    m = len(v)
    # columns: basis rows as unknowns; rows: coordinates
    A = [[basis[k][r] % p for k in range(n)] + [v[r] % p] for r in range(m)]
    piv, row = [], 0
    for c in range(n):
        pr = next((r for r in range(row, m) if A[r][c]), None)
        if pr is None:
            continue
        A[row], A[pr] = A[pr], A[row]
        inv = pow(A[row][c], p - 2, p)
        A[row] = [x * inv % p for x in A[row]]
        for r in range(m):
            if r != row and A[r][c]:
                fct = A[r][c]
                A[r] = [(a - fct * b) % p for a, b in zip(A[r], A[row])]
        piv.append(c)
        row += 1
    if any(A[r][n] for r in range(row, m)):
        return None
    c = [0] * n
    for r, col in enumerate(piv):
        c[col] = A[r][n]
    return c


class _DualMod:
    """a + b*eps over GF(p) with eps^2 = 0, for directional derivatives of lambdified
    rational expressions; constants are ints or Fractions."""
    __slots__ = ('a', 'b', 'p')

    def __init__(self, a, b, p):
        self.a, self.b, self.p = a % p, b % p, p

    def _c(self, o):
        if isinstance(o, _DualMod):
            return o
        if isinstance(o, Fraction):
            if o.denominator % self.p == 0:
                raise ZeroDivisionError
            return _DualMod(o.numerator % self.p * pow(o.denominator % self.p, self.p - 2,
                                                       self.p), 0, self.p)
        if isinstance(o, int):
            return _DualMod(o, 0, self.p)
        raise TypeError('unsupported operand')

    def __add__(self, o):
        o = self._c(o)
        return _DualMod(self.a + o.a, self.b + o.b, self.p)

    __radd__ = __add__

    def __neg__(self):
        return _DualMod(-self.a, -self.b, self.p)

    def __sub__(self, o):
        o = self._c(o)
        return _DualMod(self.a - o.a, self.b - o.b, self.p)

    def __rsub__(self, o):
        return self._c(o) - self

    def __mul__(self, o):
        o = self._c(o)
        return _DualMod(self.a * o.a, self.a * o.b + self.b * o.a, self.p)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._c(o)
        if o.a % self.p == 0:
            raise ZeroDivisionError
        inv = pow(o.a, self.p - 2, self.p)
        q = self.a * inv % self.p
        return _DualMod(q, (self.b - q * o.b) * inv, self.p)

    def __rtruediv__(self, o):
        return self._c(o) / self

    def __pow__(self, n):
        if isinstance(n, Fraction) and n.denominator == 1:
            n = n.numerator
        if not isinstance(n, int):
            raise TypeError('non-integer power')
        if n < 0:
            return _DualMod(1, 0, self.p) / (self ** (-n))
        if n == 0:
            return _DualMod(1, 0, self.p)
        an = pow(self.a, n - 1, self.p)
        return _DualMod(an * self.a, n * an * self.b, self.p)

    def __pos__(self):
        return self


def verifyDistribution(regimes, gens, events, lifts, p, npts=2, seed=1):
    """Certificate that every direction of `lifts` lies in the kernel of all orders.

    `regimes`: list of {'f': ["X = rhs"], 'g': ["y = h"]}; `gens[r]`: fields of regime r,
    each a {coordinate: expression}; `lifts`: one per chain, {'regime': r, 'point':
    {coordinate: int}, 'vectors': [{coordinate: int}], 'times': [expression]}, the start
    state and the kernel directions lifted onto it; `events`: the joins between
    consecutive segments in chain order, {'from', 'to', 'ev': [{var, method, value}],
    'lift'}. Per chain the fields that the directions reach are closed under the
    brackets with each regime's field and carried through the joins; (A) and (B) are
    checked on them. Returns {'ok': bool, 'why': str}."""
    p, npts = int(p), int(npts)
    lines = []
    for rg in regimes:
        lines += list(_as_list(rg['f'])) + list(_as_list(rg['g']))
    for gl in gens:
        for gg in gl:
            lines += ['_ = %s' % v for v in dict(gg).values()]
    for ev in events:
        lines += ['_ = %s' % e['value'] for e in _as_list(ev['ev'])]
    ck = _DistCheck(lines)
    rng = random.Random(seed)
    states = []
    for rg in regimes:
        for l in _as_list(rg['f']):
            nm = _clean(l).split('=', 1)[0].strip()
            if nm not in states:
                states.append(nm)
    F, H = [], []
    for rg in regimes:
        rhs = {}
        for l in _as_list(rg['f']):
            lhs, e = _clean(l).split('=', 1)
            rhs[lhs.strip()] = ck.expr(e)
        F.append(rhs)
        H.append([ck.expr(_clean(l).split('=', 1)[1]) for l in _as_list(rg['g'])])
    Xm = [[{str(k): ck.expr(v) for k, v in dict(gg).items()} for gg in gl] for gl in gens]
    syms = set()
    for rhs in F:
        for e in rhs.values():
            syms |= e.free_symbols
    for hs in H:
        for e in hs:
            syms |= e.free_symbols
    for gl in Xm:
        for gg in gl:
            for e in gg.values():
                syms |= e.free_symbols
    for ev in events:
        for e in _as_list(ev['ev']):
            syms |= ck.expr(e['value']).free_symbols
    for lf in lifts:
        syms |= {spy.Symbol(str(k)) for k in dict(lf['point'])}
        for v in _as_list(lf['vectors']):
            syms |= {spy.Symbol(str(k)) for k in dict(v)}
    coords = states + sorted({str(s) for s in syms} - set(states))
    inModel = set()
    for rhs in F:
        for e in rhs.values():
            inModel |= {str(x) for x in e.free_symbols}
    for hs in H:
        for e in hs:
            inModel |= {str(x) for x in e.free_symbols}
    free = [nm for nm in coords if nm not in states and nm not in inModel]
    symList = [ck.local.get(nm, spy.Symbol(nm)) for nm in coords]
    zero = spy.Integer(0)
    U = [{nm: spy.Integer(1)} for nm in free]

    def comp(field, nm):
        return field.get(nm, zero)

    # regimes share most right-hand sides, so the partials are memoised per expression
    dcache = {}

    def dpart(e, s):
        k = (e, s)
        d = dcache.get(k)
        if d is None:
            d = dcache[k] = spy.diff(e, s)
        return d

    def derivDir(e, field):
        return sum((dpart(e, s) * comp(field, str(s)) for s in e.free_symbols
                    if str(s) in field), zero)

    def vec(exprs, vals):
        out = []
        for e in exprs:
            v = ck.value(e, symList, vals, p, skey='vd')
            if v is None:
                return None
            out.append(v)
        return out

    posOf = {nm: k for k, nm in enumerate(coords)}

    def fieldAt(x, vals):
        out = [0] * len(coords)
        for nm, e in x.items():
            if nm not in posOf or e == zero:
                continue
            v = ck.value(e, symList, vals, p, skey='vd')
            if v is None:
                return None
            out[posOf[nm]] = v
        return out

    def randPoint():
        return [Fraction(rng.randrange(2, p - 1)) for _ in coords]

    bracketCache = {}

    def bracket(r, j):
        key = (r, j)
        if key not in bracketCache:
            x, Fr = Xm[r][j], {nm: F[r].get(nm, zero) for nm in states}
            bracketCache[key] = [derivDir(comp(Fr, nm), x) - derivDir(comp(x, nm), Fr)
                                 if nm in Fr else -derivDir(comp(x, nm), Fr)
                                 for nm in coords]
        return bracketCache[key]

    def fieldsAt(r, vals):
        """Values of the model fields of regime r and of the unit fields at a point."""
        out = [fieldAt(x, vals) for x in Xm[r]]
        if any(o is None for o in out):
            return None
        return out + [[1 if nm == u else 0 for nm in coords] for u in free]

    def support(r, vals, v):
        """Model fields of regime r in the decomposition of v at the point, or None."""
        B = fieldsAt(r, vals)
        if B is None:
            return None
        c = _span_coeffs(B, v, p)
        if c is None:
            return None
        nx = len(Xm[r])
        return {k for k in range(nx) if c[k]}

    # lambdified vector functions over the coordinates, evaluated on duals: the
    # derivative of an expression list along a direction without symbolic brackets
    lamCache = {}

    def lam(key, exprs):
        fn = lamCache.get(key)
        if fn is None:
            ex = []
            frac = spy.Function('Fraction')
            for e in exprs:
                e = spy.sympify(e)
                e = e.xreplace({x: spy.Rational(str(x)) for x in e.atoms(spy.Float)})
                e = e.xreplace({x: frac(x.p, x.q) for x in e.atoms(spy.Rational)
                                if not x.is_Integer})
                ex.append(e)
            fn = lamCache[key] = spy.lambdify(symList, ex, modules=[{'Fraction': Fraction}])
        return fn

    def along(key, exprs, vals, dirv):
        """Values and derivatives of exprs at vals along dirv (residues), or None."""
        args = [_DualMod(int(v.numerator) % p * pow(int(v.denominator) % p, p - 2, p),
                         d, p) for v, d in zip(vals, dirv)]
        try:
            out = lam(key, exprs)(*args)
        except Exception:
            return None
        res = []
        for o in out:
            if isinstance(o, _DualMod):
                res.append((o.a, o.b))
            elif isinstance(o, Fraction):
                if o.denominator % p == 0:
                    return None
                res.append((o.numerator % p * pow(o.denominator % p, p - 2, p), 0))
            elif isinstance(o, int):
                res.append((o % p, 0))
            else:
                return None
        return res

    # the rational right-hand sides of a regime go through duals; one with a function
    # (exp) is differentiated symbolically, where the terms along a field cancel
    ratOf = {}

    def ratIdx(r):
        if r not in ratOf:
            ratOf[r] = [k for k, nm in enumerate(states)
                        if not F[r].get(nm, zero).atoms(spy.Function)]
        return ratOf[r]

    def bracketAt(r, j, vals):
        """[X_j, F^(r)] at vals, or None on a pole: DF.X - DX.F over the coordinates."""
        X = Xm[r][j]
        x = fieldAt(X, vals)
        if x is None:
            return None
        nS = len(states)
        ri = ratIdx(r)
        df = [0] * nS
        fx = along(('F', r), [F[r].get(states[k], zero) for k in ri], vals, x)
        if fx is None:
            return None
        for k, (_, d) in zip(ri, fx):
            df[k] = d
        supp = set(X)
        for k in set(range(nS)) - set(ri):
            e = F[r].get(states[k], zero)
            if not supp & {str(t) for t in e.free_symbols}:
                continue
            v = ck.value(derivDir(e, X), symList, vals, p, skey='vd')
            if v is None:
                return None
            df[k] = v
        Fr = {nm: F[r].get(nm, zero) for nm in states}
        out = []
        for k, nm in enumerate(coords):
            cx = comp(X, nm)
            dx = 0
            if cx != zero:
                dx = ck.value(derivDir(cx, Fr), symList, vals, p, skey='vd')
                if dx is None:
                    return None
            out.append(((df[k] if k < nS else 0) - dx) % p)
        return out

    def close(r, S):
        """Close the field set S of regime r under brackets with F^(r); check (A), (B)."""
        S, todo = set(S), list(S)
        while todo:
            j = todo.pop()
            for _ in range(npts):
                vals = randPoint()
                x = fieldAt(Xm[r][j], vals)
                hr = [h for h in H[r] if not h.atoms(spy.Function)]
                hs = [h for h in H[r] if h.atoms(spy.Function)]
                ha = along(('H', r), hr, vals, x) if x is not None else None
                hb = vec([derivDir(h, Xm[r][j]) for h in hs], vals) if hs else []
                if ha is None or hb is None or any(d for _, d in ha) or any(hb):
                    return None, 'regime %d: a field changes the output' % (r + 1)
                bv = bracketAt(r, j, vals)
                if bv is None:
                    return None, 'regime %d: a bracket has a pole' % (r + 1)
                sup = support(r, vals, bv)
                if sup is None:
                    return None, ('regime %d: the fields are not invariant under the '
                                  'dynamics' % (r + 1))
                for k in sup - S:
                    S.add(k)
                    todo.append(k)
        return S, ''

    def image(ev_list, vals, x):
        """D E applied to the field or vector x at the point, and E of the point."""
        valOf = dict(zip(coords, vals))
        img = dict(valOf)
        xd = dict(zip(coords, x))
        for e in ev_list:
            ve = ck.expr(e['value'])
            v = ck.value(ve, symList, vals, p, skey='vd')
            if v is None:
                return None, None
            dv = 0
            for s_ in ve.free_symbols:
                d = ck.value(spy.diff(ve, s_), symList, vals, p, skey='vd')
                if d is None:
                    return None, None
                dv = (dv + d * x[coords.index(str(s_))]) % p if str(s_) in coords else dv
            s0 = int(valOf[e['var']]) % p
            img[e['var']] = Fraction({'replace': v, 'add': (s0 + v) % p,
                                      'multiply': s0 * v % p}[e['method']])
            xd[e['var']] = {'replace': dv, 'add': (xd[e['var']] + dv) % p,
                            'multiply': (xd[e['var']] * v + s0 * dv) % p}[e['method']]
        return [xd[nm] for nm in coords], [img[nm] for nm in coords]

    joinsOf = {}
    for ev in events:
        if ev.get('lift') is not None:
            joinsOf.setdefault(int(ev['lift']), []).append(ev)

    for li, lf in enumerate(lifts):
        r = int(lf['regime'])
        pt = {str(k): int(v) for k, v in dict(lf['point']).items()}
        vals0 = [Fraction(pt.get(nm, 0)) for nm in coords]
        taus = [ck.expr(t) for t in _as_list(lf.get('times', []))]
        vecs = [{str(k): int(x) % p for k, x in dict(v).items()}
                for v in _as_list(lf['vectors'])]
        # (D) the lifted directions in the span at the start, and the fields they reach
        S = set()
        for vd in vecs:
            sup = support(r, vals0, [vd.get(nm, 0) for nm in coords])
            if sup is None:
                return {'ok': False, 'why': 'a kernel direction is not in the span of the fields'}
            S |= sup
            for tau in taus:
                dv = [ck.value(spy.diff(tau, s_), symList, vals0, p, skey='vd') for s_ in tau.free_symbols]
                if any(x is None for x in dv) or sum(
                        x * vd.get(str(sy), 0) for x, sy in zip(dv, tau.free_symbols)) % p:
                    return {'ok': False, 'why': 'a kernel direction moves an event time'}
        S, why = close(r, S)
        if S is None:
            return {'ok': False, 'why': why}
        us = [[vd.get(nm, 0) if nm in free else 0 for nm in coords] for vd in vecs]
        # (C) the joins of the chain in order
        for ev in joinsOf.get(li, []):
            r2 = int(ev['to'])
            evs = _as_list(ev['ev'])
            S2 = set()
            for _ in range(npts):
                vals = randPoint()
                items = [fieldAt(Xm[r][j], vals) for j in S] + us
                for x in items:
                    if x is None:
                        return {'ok': False, 'why': 'a field has a pole before an event'}
                    xi, zi = image(evs, vals, x)
                    if xi is None:
                        return {'ok': False, 'why': 'an event value has a pole'}
                    sup = support(r2, zi, xi)
                    if sup is None:
                        return {'ok': False, 'why': 'an event leaves the span of the fields'}
                    S2 |= sup
            S, why = close(r2, S2)
            if S is None:
                return {'ok': False, 'why': why}
            r = r2
    return {'ok': True, 'why': ''}


def extendFields(fields, icLines, eventGroups, states):
    """Components of each field on the parameters that are the value of an initial value
    or an event, from compatibility with that map. `eventGroups` holds the events of one
    time each, applied together as E. x0 = X gives eta_X = eta_x(x0); an event x -> x + v,
    v, x v on a state s gives eta_v = eta_s(E x) - eta_s(x), eta_s(E x),
    (eta_s(E x) - v eta_s(x)) / x_s. A component that depends on the states means no
    extension, and the field is dropped. Returns the extended fields."""
    icLines = _as_list(icLines)
    groups = [_as_list(g) for g in _as_list(eventGroups)]
    lines = list(icLines) + ['_ = %s' % e['value'] for g in groups for e in g]
    for fld in _as_list(fields):
        lines += ['_ = %s' % v for v in dict(fld).values()]
    local, parse = _make_local_parse(lines + ['_ = %s' % st for st in _as_list(states)])
    stateSet = set(_as_list(states))
    sym = lambda nm: local.get(nm, spy.Symbol(nm))
    stateSyms = {sym(nm) for nm in stateSet}
    out = []
    for fld in _as_list(fields):
        eta = {str(k): spy.sympify(parse(str(v))) for k, v in dict(fld).items()}
        ext = {}
        for l in icLines:
            lhs, rhs = _clean(l).split('=', 1)
            x, e = lhs.strip(), spy.sympify(parse(rhs))
            if e.is_Symbol and str(e) not in stateSet:
                ext[str(e)] = eta.get(x, spy.Integer(0)).subs(sym(x), e)
        for g in groups:
            vals = {str(ev['var']): spy.sympify(parse(str(ev['value']))) for ev in g}
            E = {}
            for ev in g:
                s_, v = str(ev['var']), vals[str(ev['var'])]
                xs = sym(s_)
                E[xs] = {'add': xs + v, 'replace': v, 'multiply': xs * v}[ev['method']]
            for ev in g:
                s_, v = str(ev['var']), vals[str(ev['var'])]
                if not (v.is_Symbol and str(v) not in stateSet):
                    continue
                xs = sym(s_)
                es = eta.get(s_, spy.Integer(0))
                eE = es.xreplace(E) if E else es
                ext[str(v)] = {'add': eE - es, 'replace': eE,
                               'multiply': (eE - v * es) / xs}[ev['method']]
        ok = True
        for k, val in ext.items():
            val = spy.simplify(val)
            if val.free_symbols & stateSyms:
                ok = False
                break
            if val != 0:
                eta[k] = val
        if ok:
            out.append({k: str(v) for k, v in eta.items()})
    return out


def _ser_mul(a, b, p, T):
    out = [0] * (T + 1)
    for i, x in enumerate(a):
        if x:
            for j in range(T + 1 - i):
                if b[j]:
                    out[i + j] = (out[i + j] + x * b[j]) % p
    return out


def _ser_inv(a, p, T):
    """Inverse of a truncated power series with a unit constant term, or None."""
    if a[0] % p == 0:
        return None
    b = [0] * (T + 1)
    b[0] = pow(a[0] % p, p - 2, p)
    for k in range(1, T + 1):
        acc = sum(a[i] * b[k - i] for i in range(1, k + 1)) % p
        b[k] = (-acc) * b[0] % p
    return b


def _ser_terms(terms, gser, p, T, pw):
    """A term list at series values of the generators, truncated at eps^T; `pw` caches
    the powers of each generator."""
    out = [0] * (T + 1)
    for expo, (cn, cd) in terms:
        c = cn % p * pow(cd % p, p - 2, p) % p
        if not c:
            continue
        acc = [c] + [0] * T
        for g, e in enumerate(expo):
            if e:
                key = (g, e)
                if key not in pw:
                    s_ = [1] + [0] * T
                    for _ in range(e):
                        s_ = _ser_mul(s_, gser[g], p, T)
                    pw[key] = s_
                acc = _ser_mul(acc, pw[key], p, T)
        out = [(x + y) % p for x, y in zip(out, acc)]
    return out


def _ser_value(numden, gser, p, T, pw):
    num = _ser_terms(numden[0], gser, p, T, pw)
    den = _ser_terms(numden[1], gser, p, T, pw)
    if den == [1] + [0] * T:
        return num
    inv = _ser_inv(den, p, T)
    return None if inv is None else _ser_mul(num, inv, p, T)


_contCache = {}


def continueRestingState(model, stateNames, paramNames, paramVals, dirs, restVals, prime,
                         order, forcings=None):
    """Resting state of `model` as a power series in eps along the parameters
    paramVals + eps dirs, continued from restVals, a resting state at eps = 0, by
    Newton steps: x_k = -J0^-1 [eps^k] f(x_<k). Returns {'ok', 'valBy': {state:
    [x_0, ..., x_T]}, 'dfJx': [[series]], 'dfJt': {param: [series per balance]},
    'dfStateCols', 'dfParamCols'} or {'ok': False, 'why'}."""
    p, T = int(prime), int(order)
    model = _as_list(model)
    forcings = set(_as_list(forcings))
    stateNames = _as_list(stateNames)
    paramNames = _as_list(paramNames)
    (paramSyms, solveStates, polys, Jx, Jt, gens, JxTerms, JtTerms, polyBi, genLin,
     linTerms, pointPlan) = _ss_compile(model, stateNames, paramNames, forcings)
    key = (tuple(model), tuple(stateNames), tuple(paramNames), tuple(sorted(forcings)))
    cached = _contCache.get(key)
    if cached is None:
        fT = [_poly_terms(pl, gens) for pl in polys]
        # Newton linearises the cleared numerators, not f
        JnT = [[_poly_terms(spy.diff(pl, s_), gens) for s_ in solveStates] for pl in polys]
        cached = (fT, JnT)
        _contCache[key] = cached
    fT, JnT = cached
    pv = {str(k): int(v) % p for k, v in dict(paramVals).items()}
    dv = {str(k): int(v) % p for k, v in dict(dirs).items()}
    rv = {str(k): int(v) % p for k, v in dict(restVals).items()}
    nS = len(solveStates)
    pser = [[pv.get(str(th), 0)] + [dv.get(str(th), 0)] + [0] * (T - 1) if T >= 1
            else [pv.get(str(th), 0)] for th in paramSyms]
    xser = [[rv.get(str(s_), 0)] + [0] * T for s_ in solveStates]

    def f_series():
        pw = {}
        out = []
        for t in fT:
            v = _ser_value(t, pser + xser, p, T, pw)
            if v is None:
                return None
            out.append(v)
        return out

    f0 = f_series()
    if f0 is None or any(v[0] for v in f0):
        return {'ok': False, 'why': 'the start is not a resting state at eps = 0'}
    g0 = [pv.get(str(th), 0) for th in paramSyms] + [rv.get(str(s_), 0) for s_ in solveStates]
    J0 = [[_eval_terms(JnT[i][j], g0, p) for j in range(nS)] for i in range(nS)]
    for k in range(1, T + 1):
        fs = f_series()
        if fs is None:
            return {'ok': False, 'why': 'a pole along the continuation'}
        X = _solve_mod(J0, [[(-fs[i][k]) % p] for i in range(nS)], p)
        if X is None:
            return {'ok': False, 'why': 'a singular resting Jacobian'}
        for i in range(nS):
            xser[i][k] = X[i][0] % p
    pw = {}
    gser = pser + xser

    def ser_of(terms):
        v = _ser_value(terms, gser, p, T, pw)
        return v if v is not None else [0] * (T + 1)

    dfJx = [[ser_of(JxTerms[i][j]) for j in range(nS)] for i in range(nS)]
    dfJt = {str(th): [ser_of(JtTerms[str(th)][i]) for i in range(nS)] for th in paramSyms}
    return {'ok': True, 'valBy': {str(solveStates[i]): xser[i] for i in range(nS)},
            'dfJx': dfJx, 'dfJt': dfJt,
            'dfStateCols': [str(x) for x in solveStates],
            'dfParamCols': [str(x) for x in paramSyms]}


_LIFT_CACHE = {}


def liftStart(icLines, t0events, point, vector, p, joint=False):
    """Start state of a chain and a direction lifted onto it: the initial values `icLines`
    ("X = e", unless `joint`, where the resting state and its direction are in `point`
    and `vector`), then the t0 events. Returns {'point', 'vector'} over states and
    parameters, or None on a pole."""
    p = int(p)
    icLines = _as_list(icLines)
    evs = _as_list(t0events)
    pt = {str(k): int(v) % p for k, v in dict(point).items()}
    vv = {str(k): int(v) % p for k, v in dict(vector).items()}
    names = sorted(pt)
    # parsed lines, lambdified values and derivatives are reused across points and primes
    key = (tuple(icLines), tuple((e['var'], str(e['value']), e['method']) for e in evs),
           tuple(names))
    ck = _LIFT_CACHE.get(key)
    if ck is None:
        if len(_LIFT_CACHE) > 256:
            _LIFT_CACHE.clear()
        ck = _DistCheck(icLines + ['_ = %s' % e['value'] for e in evs] +
                        ['_ = %s' % k for k in names])
        ck.parsed, ck.derivs = {}, {}
        _LIFT_CACHE[key] = ck
    syms = [ck.local.get(n, spy.Symbol(n)) for n in names]
    vals = [Fraction(pt[n]) for n in names]

    def parsed(text):
        e = ck.parsed.get(text)
        if e is None:
            e = ck.parsed[text] = ck.expr(text)
        return e

    skey = tuple(names)

    def val_and_dir(e):
        v = ck.value(e, syms, vals, p, skey=skey)
        ds = ck.derivs.get(e)
        if ds is None:
            ds = ck.derivs[e] = [(str(s), spy.diff(e, s)) for s in e.free_symbols]
        dv = 0
        for s, de in ds:
            d = ck.value(de, syms, vals, p, skey=skey)
            if d is None:
                return None, None
            dv = (dv + d * vv.get(s, 0)) % p
        return v, dv

    if not joint:
        for l in icLines:
            lhs, e = _clean(l).split('=', 1)
            x, dx = val_and_dir(parsed(e))
            if x is None:
                return None
            pt[lhs.strip()], vv[lhs.strip()] = x, dx
    for e in evs:
        v, dv = val_and_dir(parsed(str(e['value'])))
        if v is None:
            return None
        s0, d0 = pt.get(e['var'], 0), vv.get(e['var'], 0)
        if e['method'] == 'replace':
            pt[e['var']], vv[e['var']] = v, dv
        elif e['method'] == 'add':
            pt[e['var']], vv[e['var']] = (s0 + v) % p, (d0 + dv) % p
        else:
            pt[e['var']], vv[e['var']] = s0 * v % p, (d0 * v + s0 * dv) % p
    return {'point': pt, 'vector': vv}


def _scaling_nonid(gens, znames, nz):
    """Scaling entries (support, integer vector, type) from generators."""
    nonId = []
    for c in gens:
        comp = {znames[j]: c[j] for j in range(nz) if c[j] != 0}
        nonId.append({'support': sorted(comp.keys()),
                      'vector': {k: str(v) for k, v in comp.items()},
                      'type': 'scaling'})
    return nonId


def _scaling_gens_recast(gens, znames, nz, recast):
    """Impose c_E = exp * c_base on the integer scaling lattice span(gens) and return
    the physical generators, weights possibly symbolic in the exponent. No integer
    weight satisfies c_E = nhill * c_base, so the plain kernel misses Hill scalings;
    intersecting over Q(exp) recovers them exactly (e.g. xi_kinh = -nhill * kinh).
    Holds for inverted and normal recasts alike."""
    if not gens:
        return []
    G = spy.Matrix(gens)                      # rows = basis vectors, cols = coords
    idx = {nm: j for j, nm in enumerate(znames)}
    relRows = []
    for rc in recast:
        E, base, exp = str(rc['E']), str(rc['base']), str(rc['exp'])
        if E not in idx or base not in idx:
            continue
        expSym = spy.sympify(exp)
        jE, jB = idx[E], idx[base]
        # (G[:,E] - exp*G[:,base]) . alpha = 0 over the basis coefficients alpha
        relRows.append([G[i, jE] - expSym * G[i, jB] for i in range(G.rows)])
    if not relRows:
        return [[spy.Integer(x) for x in g] for g in gens]
    alphas = spy.Matrix(relRows).nullspace()   # over Q(exp)
    out = []
    for a in alphas:
        v = [spy.together(sum(a[i] * G[i, j] for i in range(G.rows)))
             for j in range(nz)]
        # L = log(base) shifts by c_base*log(lam) under the scaling; the kernel above
        # leaves it at 0, the joint nullspace carries the shift
        for rc in recast:
            base, L = str(rc['base']), str(rc['L'])
            if base in idx and L in idx:
                v[idx[L]] = v[idx[base]]
        # clear denominators so the weights are polynomial in the exponents
        dens = [spy.denom(x) for x in v if x != 0]
        Lden = spy.Integer(1)
        for d in dens:
            Lden = spy.lcm(Lden, d)
        v = [spy.expand(x * Lden) for x in v]
        if any(x != 0 for x in v):
            out.append(v)
    return out


def scalingSymmetries(allVariables, diffEquations, obsFunctions, m, params,
                      fixed=(), verbose=True):
    """Exact scaling symmetries z_i -> lam^{c_i} z_i from the integer kernel of the
    monomial-exponent conditions (observables invariant, f_i scaling like x_i); no
    expression swell. Symbols in `fixed` get weight zero."""
    zvars = list(allVariables[:m]) + list(params)
    nz = len(zvars)
    znames = [str(s) for s in zvars]
    fixedset = set(str(s) for s in fixed)

    sparse, ninter, skipped = _scaling_rows(diffEquations, obsFunctions, m, zvars, nz,
                                            logs=True)
    ncols = nz + ninter
    rows = _materialize_rows(sparse, ncols)
    for j in range(nz):                      # a fixed coordinate does not scale
        if znames[j] in fixedset:
            row = [0] * ncols
            row[j] = 1
            rows.append(row)

    nonId = _scaling_nonid(_scaling_gens(rows, ncols, nz), znames, nz)

    if verbose:
        print('-' * 60)
        print('%d scaling symmetry/ies (exact integer kernel)%s:'
              % (len(nonId), '' if not skipped else
                 ', %d non-polynomial term(s) skipped' % skipped))
        for d in nonId:
            print('  ' + ', '.join('%s^(%s)' % (k, d['vector'][k])
                                   for k in d['support']))

    return {
        'method': 'scaling',
        'count': len(nonId),
        'nonIdentifiable': nonId,
        'coordinates': znames,
    }


def scalingSymmetriesMulti(perCondModel, perCondObs, inputs=None, fixed=None,
                           recast=None, logs=False, extraModel=None, extraObs=None):
    """Scaling symmetries common to every condition: the integer kernel of all
    conditions' monomial-exponent rows stacked over a shared weight space (own
    intermediate columns each), i.e. the intersection of the per-condition lattices.
    Coordinates are states plus parameters; `inputs` and `fixed` do not scale.
    `extraModel` ("X = e", e of the weight of state X) and `extraObs` ("_ = e", e of
    weight 0) add constraints per condition, such as initial values and events."""
    perCondModel = [_as_list(m) for m in perCondModel]
    perCondObs = [_as_list(o) for o in perCondObs]
    extraModel = [_as_list(m) for m in (extraModel or [[] for _ in perCondModel])]
    extraObs = [_as_list(o) for o in (extraObs or [[] for _ in perCondModel])]
    inputset = set(_as_list(inputs))
    fixedset = set(_as_list(fixed))
    K = len(perCondModel)
    if K == 0:
        return {'method': 'scaling', 'count': 0, 'nonIdentifiable': [],
                'coordinates': []}

    all_lines = [l for lines in perCondModel + perCondObs + extraModel + extraObs
                 for l in lines]
    local, parse = _make_local_parse(all_lines)

    stateNames = [_clean(l).split('=', 1)[0].strip() for l in perCondModel[0]]
    stateSyms = [local.get(nm, spy.Symbol(nm)) for nm in stateNames]

    perF, perG, paramset = [], [], set()
    for c in range(K):
        rhs = {}
        for l in perCondModel[c]:
            lhs, expr = _clean(l).split('=', 1)
            rhs[lhs.strip()] = spy.sympify(parse(expr))
        f = [rhs[nm] for nm in stateNames]
        g = [spy.sympify(parse(_clean(l).split('=', 1)[1])) for l in perCondObs[c]]
        perF.append(f)
        perG.append(g)
        for e in f + g:
            paramset |= set(spy.sympify(e).free_symbols)
    perXF, perXG = [], []
    for c in range(K):
        xf = []
        for l in extraModel[c]:
            lhs, expr = _clean(l).split('=', 1)
            if lhs.strip() in stateNames:
                e = spy.sympify(parse(expr))
                xf.append((stateNames.index(lhs.strip()), e))
                paramset |= set(e.free_symbols)
        xg = [spy.sympify(parse(_clean(l).split('=', 1)[1])) for l in extraObs[c]]
        for e in xg:
            paramset |= set(e.free_symbols)
        perXF.append(xf)
        perXG.append(xg)

    paramset -= set(stateSyms)
    paramset -= {local.get(nm, spy.Symbol(nm)) for nm in inputset}
    params = sorted(paramset, key=spy.default_sort_key)
    # re-add free-exponent bases the recast eliminated (Km in C^n/(Km^n+C^n)), so
    # _scaling_gens_recast can impose c_E = exp*c_base on them
    for rc in _as_list(recast):
        b = local.get(str(rc['base']), spy.Symbol(str(rc['base'])))
        if (b not in set(stateSyms) and b not in params
                and str(b) not in inputset and str(b) not in fixedset):
            params.append(b)
    params = sorted(params, key=spy.default_sort_key)
    zvars = stateSyms + params
    nz = len(zvars)
    znames = [str(s) for s in zvars]
    m = len(stateSyms)

    rows, interOffset, skipped, seen = [], nz, 0, set()
    for c in range(K):
        sparse, ninter, sk = _scaling_rows(perF[c], perG[c], m, zvars, interOffset,
                                           logs, seen)
        rows.extend(sparse)
        interOffset += ninter
        skipped += sk
        for i, e in perXF[c]:
            fx = [spy.Integer(0)] * m
            fx[i] = e
            sparse, ninter, sk = _scaling_rows(fx, [], m, zvars, interOffset, logs, seen)
            rows.extend(sparse)
            interOffset += ninter
            skipped += sk
        if perXG[c]:
            sparse, ninter, sk = _scaling_rows([spy.Integer(0)] * m, perXG[c], m, zvars,
                                               interOffset, logs, seen)
            rows.extend(sparse)
            interOffset += ninter
            skipped += sk
    if skipped:
        print('scalingSymmetriesMulti: %d non-polynomial term(s) skipped '
              '(a scaling they would forbid may be over-reported)' % skipped)
    ncols = interOffset
    dense = _materialize_rows(rows, ncols)
    for j in range(nz):
        if znames[j] in fixedset:
            row = [0] * ncols
            row[j] = 1
            dense.append(row)

    gens = _scaling_gens(dense, ncols, nz)
    # the integer kernel treats each recast E as free; recover the Hill scalings
    recast = _as_list(recast)
    if recast:
        physical = _scaling_gens_recast(gens, znames, nz, recast)
        nonId = _scaling_nonid(physical, znames, nz)
    else:
        nonId = _scaling_nonid(gens, znames, nz)
    return {'method': 'scaling', 'count': len(nonId), 'nonIdentifiable': nonId,
            'coordinates': znames}


###########################################################################
#####################  pure-symbolic observability  ######################
###########################################################################

# ---- pure-symbolic observability cross-check -----------------------------------------

def _lie_deriv(h, states, rhs):
    """Lie derivative L_f h = sum_i (dh/dx_i) f_i; parameters are constant."""
    return sum(spy.diff(h, states[i]) * rhs[i] for i in range(len(states)))


def observabilitySympyMulti(model, observation, conditionSubs=None, conditionIC0=None,
                            fixed=None, parameters=None, inputs=None, backend='sympy',
                            conditionObs=None):
    """Multi-condition pure-symbolic observability, the exact cross-check of the
    modular multi-condition engine (same coordinate and substitution semantics as
    compileObservabilityTapeMulti).

    Condition k substitutes `conditionSubs[k]` into f and g, starts the jet at
    `conditionIC0[k]` and may replace `observation` by `conditionObs[k]`. A direction
    is non-identifiable iff it lies in the nullspace of the observability matrices
    stacked over one shared coordinate space. Single segment only and no equilibrate;
    give a steady state through the initial values (`trafo`)."""
    _select_backend(backend)
    asL = lambda v: list(v) if isinstance(v, (list, tuple)) else ([] if v is None else [v])
    model = [str(l) for l in asL(model)]
    observation = [str(l) for l in asL(observation)]
    conditionSubs = [dict(c) for c in asL(conditionSubs)] or [{}]
    conditionIC0 = [dict(c) for c in asL(conditionIC0)]
    K = len(conditionSubs)
    forcings = set(str(x) for x in asL(inputs))
    conditionObs = list(conditionObs) if conditionObs else []
    obsPerCond = [[str(l) for l in asL(conditionObs[c])]
                  if (c < len(conditionObs) and conditionObs[c] is not None)
                  else observation for c in range(K)]

    all_lines = model + observation + asL(parameters)
    for o in obsPerCond:
        all_lines += o
    for c in conditionSubs + conditionIC0:
        for k, v in c.items():
            all_lines += [str(k), str(v)]
    local, parse = _make_local_parse(all_lines)

    variables, diffEquations, _ = _read_equations(model, parse)
    obsFunPer = [_read_equations(o, parse)[1] for o in obsPerCond]
    fixedNames = set(str(s) for s in
                     [local.get(nm, spy.Symbol(nm)) for nm in asL(fixed)])

    # constant (u' = 0) states are dropped; R has baked their per-condition value
    # into conditionSubs, so they enter f/g as substituted constants
    isConst = [spy.sympify(e) == 0 for e in diffEquations]
    S = [variables[i] for i in range(len(variables)) if not isConst[i]]
    Srhs = [spy.sympify(diffEquations[i]) for i in range(len(variables))
            if not isConst[i]]

    def pval(s):
        return spy.sympify(parse(str(s)))

    # per-condition substituted dynamics/observation and initial state
    perCond = []
    for c in range(K):
        subsMap = {local.get(str(k), spy.Symbol(str(k))): pval(v)
                   for k, v in conditionSubs[c].items()}
        ic0 = conditionIC0[c] if c < len(conditionIC0) else {}
        f_c = [e.subs(subsMap) for e in Srhs]
        g_c = [spy.sympify(e).subs(subsMap) for e in obsFunPer[c]]
        ic_c = {}
        for X in S:
            e = pval(ic0[str(X)]) if str(X) in ic0 else X
            ic_c[X] = spy.sympify(e).subs(subsMap)
        perCond.append((f_c, g_c, ic_c))

    # a state is a free coordinate iff its symbol survives in its initial value in
    # some condition; parameters are all remaining free symbols; both minus `fixed`
    freeStates = [X for X in S
                  if any(X in ic_c[X].free_symbols for (_, _, ic_c) in perCond)]
    paramset = set()
    for (f_c, g_c, ic_c) in perCond:
        for e in list(f_c) + list(g_c) + list(ic_c.values()):
            paramset |= set(spy.sympify(e).free_symbols)
    paramset -= set(S)
    params = sorted(paramset, key=spy.default_sort_key)

    z = [X for X in freeStates if str(X) not in fixedNames] + \
        [s for s in params if str(s) not in fixedNames]
    znames = [str(s) for s in z]
    nz = len(z)
    if nz == 0:
        return {'ok': True, 'rank': 0, 'dim': 0, 'lieOrder': 0,
                'nonIdentifiable': [], 'identifiable': True, 'coordinates': []}

    # per Lie order, each condition's rows d/dz[L_{f_c}^order g_c at ic_c]; stop at
    # full rank or after two orders without rank gain, counted only past the
    # structural first order (a readout behind a transit chain is flat until it fills)
    reachPer = [_lie_reach(f_c, g_c, S, ic_c, []) for (f_c, g_c, ic_c) in perCond]
    firstOrd = [min([r[nm] for r in reachPer if nm in r] or [-1]) for nm in znames]
    lo0 = max([k for k in firstOrd if k >= 0] or [0])
    rows = []
    jets = [list(g_c) for (f_c, g_c, ic_c) in perCond]     # order 0: the observables
    prev, flat, order = -1, 0, 0
    while True:
        for ci, (f_c, g_c, ic_c) in enumerate(perCond):
            for h in jets[ci]:
                he = spy.sympify(h).subs(ic_c)             # evaluate at x(0) = ic_c
                rows.append([spy.diff(he, zj) for zj in z])
        rank = spy.Matrix(rows).rank()
        if rank == prev:
            flat = flat + 1 if order > lo0 else flat
        else:
            flat = 0
        if rank >= nz or (flat >= 2 and order >= 1) or order >= nz + lo0:
            break
        prev = rank
        order += 1
        for ci, (f_c, g_c, ic_c) in enumerate(perCond):
            jets[ci] = [spy.expand(_lie_deriv(h, S, f_c)) for h in jets[ci]]

    M = spy.Matrix(rows)
    rank = M.rank()
    nonId = []
    if rank < nz:
        for vec in M.nullspace():
            v = [spy.cancel(vec[i]) for i in range(nz)]
            nz_i = [i for i in range(nz) if v[i] != 0]
            nonId.append({'support': sorted(znames[i] for i in nz_i),
                          'vector': {znames[i]: str(v[i]) for i in nz_i},
                          'type': 'general', 'closedForm': True})
    return {'ok': True, 'rank': int(rank), 'dim': int(nz), 'lieOrder': int(order),
            'nonIdentifiable': nonId, 'identifiable': bool(rank >= nz),
            'coordinates': znames}


###########################################################################
#####################     R entry point     ##############################
###########################################################################

# ---- R-facing entry point and I/O ----------------------------------------------------

def symmetryDetectiondMod(model, observation,
                          ansatz='uni', pMax=2, inputs=None, fixed=None,
                          allTrafos=False, lieOrder=0, exact=True,
                          verify=True, backend='sympy', parameters=None,
                          method='liesym'):
    global _EXACT
    _EXACT = bool(exact)
    _select_backend(backend)

    model = _as_list(model)
    observation = _as_list(observation)
    inputNames = _as_list(inputs)
    fixedNames = _as_list(fixed)
    extraParams = _as_list(parameters)

    all_lines = model + observation + inputNames + fixedNames + extraParams
    local, parse = _make_local_parse(all_lines)

    sys.stdout.write('\nReading input...')
    sys.stdout.flush()

    variables, diffEquations, params = _read_equations(model, parse)

    _, obsFunctions, params = _read_observation(
        observation, variables, params, parse)

    # explicit parameters
    for pn in extraParams:
        s = local.get(pn, spy.Symbol(pn))
        if s not in params and s not in variables:
            params.append(s)

    inputSyms = [local.get(nm, spy.Symbol(nm)) for nm in inputNames]
    fixedSyms = [local.get(nm, spy.Symbol(nm)) for nm in fixedNames]
    for s in inputSyms:
        if s in params:
            params.remove(s)

    # a declared input that is also a state (a known switch) is a constant, not
    # a dynamic variable: drop its equation so it appears once, as an input
    keep = [i for i, v in enumerate(variables) if v not in inputSyms]
    variables = [variables[i] for i in keep]
    diffEquations = [diffEquations[i] for i in keep]

    allVariables = list(variables) + list(inputSyms) + list(params)

    sys.stdout.write('done\n')
    sys.stdout.flush()

    if str(method) == 'scaling':
        return scalingSymmetries(
            allVariables, diffEquations, obsFunctions, len(variables), params,
            fixed=fixedSyms)

    return symmetryDetection(
        allVariables, diffEquations, obsFunctions, ansatz=ansatz, pMax=int(pMax),
        inputs=inputSyms, fixed=fixedSyms, lieOrder=int(lieOrder),
        allTrafos=bool(allTrafos), verify=bool(verify))


def _read_observation(observation, variables, parameters, parse):
    obsVars, obsFunctions, obsParameters = _read_equations(observation, parse)
    for var in variables:
        if var in obsParameters:
            obsParameters.remove(var)
    for par in parameters:
        if par in obsParameters:
            obsParameters.remove(par)
    return obsVars, obsFunctions, parameters + obsParameters

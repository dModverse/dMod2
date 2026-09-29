# AlyssaPetit version 1.4
# Use with python 3.x
#
# A new core on the same interface as v1.3. Every flux is expanded into atoms,
# products and quotients of positive symbols with a positive coefficient, and
# every balance becomes a linear form over these atoms with exact rational
# coefficients (the stoichiometry). An unknown, a state or a rate constant,
# is solved from a combination of balances found by a linear program: all atoms
# that contain the unknown on one side, only positive atoms on the other. So
#   - a solution is a ratio of positive sums by construction, nothing is ever
#     expanded to prove it,
#   - differences cancel exactly between balances (moiety sums, recycling
#     loops, a ligand that leaves only through its receptors), because
#     identical atoms carry one coefficient,
#   - solutions stay lazy and small: they reference each other, the LP keeps
#     the references acyclic, and one resolution happens at output time.
# Where the greedy choice of unknowns runs into a dead end, the states whose
# solutions block it are left free and the search restarts.

import os
import sys
import csv
import re
import time
import random
from fractions import Fraction

import numpy
import sympy
from sympy.parsing.sympy_parser import parse_expr as _parse_expr_sympy
import re as _re_pe
# Model names are symbols: every identifier that is not called is bound to a Symbol,
# so species such as Ci, Si, E, S, Q or gamma do not resolve to sympy objects
# (Symbol*Ci raised "unsupported operand type(s)" for a state named Ci).
_PE_IDENT = _re_pe.compile(r'\b([A-Za-z_][A-Za-z0-9_]*)\b(?!\s*\()')


def parse_expr(s, local_dict=None, **kwargs):
    s = str(s)
    ld = {nm: sympy.Symbol(nm) for nm in set(_PE_IDENT.findall(s))}
    if local_dict:
        ld.update(local_dict)
    return _parse_expr_sympy(s, local_dict=ld, **kwargs)
from scipy.optimize import linprog

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import AlyssaPetit_ver1_3 as _v13

# Output: 0 = result only, 1 = progress and result (default), 2 = full trace.
_VERBOSE=1

def _trace(msg):
    if _VERBOSE>=2:
        print(msg, flush=True)

_BLANK=[False]

def _say(msg):
    # progress; consecutive blank lines collapse into one
    if _VERBOSE<1 or (msg=='' and _BLANK[0]):
        return
    _BLANK[0]=(msg=='')
    print(msg, flush=True)

# Notes are collected during the run and printed once, after the summary.
_NOTES=[]
_STUCK=[None]

def _note(msg):
    _NOTES.append(msg)

def _print_notes():
    for msg in _NOTES:
        print('  note: '+msg, flush=True)


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Model
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

def _read_model(filename, injections):
    # States, fluxes and stoichiometry from the csv written by steadyStates().
    # Forcings are held at 0: substituted into the fluxes, their rows dropped.
    with open(filename) as fh:
        rows=list(csv.reader(fh, delimiter=','))
    states=[parse_expr(s) for s in rows[0][2:]]
    forc={parse_expr(s) for s in injections}
    fluxes=[]
    stoich=[]
    for row in rows[1:]:
        f=parse_expr(row[1].replace('^', '**')).xreplace({s: 0 for s in forc})
        col=[sympy.Rational(v) if v!='' else sympy.S.Zero for v in row[2:]]
        if f==0:
            continue
        fluxes.append(f)
        stoich.append(col)
    keep=[i for i, s in enumerate(states) if s not in forc]
    states=[states[i] for i in keep]
    stoich=[[col[i] for i in keep] for col in stoich]
    return states, fluxes, stoich

def _zero_states(states, fluxes, stoich, volumes):
    # States that are 0 at every steady state: a set whose mass only leaks
    # (v1.3's linear program, in amounts). They are removed and substituted by 0.
    zero=[]
    while True:
        if not states or not fluxes:
            break
        SM=sympy.Matrix(len(states), len(fluxes),
                        lambda i, j: stoich[j][i])
        cluster=_v13.FindSinkCluster(SM, F=fluxes, X=states, volumes=volumes)
        if not cluster:
            break
        zs=[states[i] for i in cluster]
        zero+=zs
        rep={s: 0 for s in zs}
        keep=[i for i in range(len(states)) if i not in cluster]
        nf, ns=[], []
        for f, col in zip(fluxes, stoich):
            f=f.xreplace(rep)
            if f==0:
                continue
            nf.append(f)
            ns.append([col[i] for i in keep])
        states=[states[i] for i in keep]
        fluxes, stoich=nf, ns
    return states, fluxes, stoich, zero

def _atoms(states, fluxes, stoich):
    # Expand every flux into atoms; the balance of state i is
    # sum_a B[i][a] * atom_a with exact rational B.
    index={}
    atoms=[]
    B=[dict() for _ in states]
    for f, col in zip(fluxes, stoich):
        for t in sympy.Add.make_args(sympy.expand(f)):
            coeff, mono=t.as_coeff_Mul()
            if not coeff.is_positive:
                raise ValueError('flux term '+str(t)+' is not positive')
            if mono not in index:
                index[mono]=len(atoms)
                atoms.append(mono)
            a=index[mono]
            for i, s in enumerate(col):
                if s!=0:
                    v=B[i].get(a, Fraction(0))+Fraction(int(s.p), int(s.q))*Fraction(int(coeff.p), int(coeff.q))
                    if v==0:
                        B[i].pop(a, None)
                    else:
                        B[i][a]=v
    return atoms, B

def _linear_factor(atom, s):
    # True if s is a plain factor of the atom's numerator, degree 1.
    n, d=sympy.fraction(atom)
    if d.has(s):
        return False
    return n.as_powers_dict().get(s, 0)==1

def _candidates(states, atoms, neglect):
    # Symbols linear in at least one atom: states and rate constants. Per
    # symbol, the atoms linear in it (usable on its side of an equation) and
    # the atoms it sits in otherwise (a denominator, a square), which the
    # equation for it must leave out.
    occ={}
    for a, at in enumerate(atoms):
        for s in at.free_symbols:
            occ.setdefault(s, []).append(a)
    cands={}
    nonlin={}
    held={s for s in occ if str(s) in neglect}
    for s, idx in occ.items():
        if s in held:
            continue
        lin=[a for a in idx if _linear_factor(atoms[a], s)]
        # a rate constant whose fluxes all carry a neglected symbol would be
        # solved as a ratio over it, and that symbol may be 0
        if lin and all(atoms[a].free_symbols & held for a in lin) and s not in states:
            continue
        if lin:
            cands[s]=lin
            nonlin[s]={a for a in idx if a not in lin}
    return cands, nonlin


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Exact linear algebra
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

def _row_basis(B, n_atoms):
    # Independent balances, exact: Gaussian elimination over Fractions.
    rows=[dict(r) for r in B if r]
    basis=[]
    pivots=[]
    for r in rows:
        r=dict(r)
        for (p, b) in zip(pivots, basis):
            if p in r:
                f=r[p]/b[p]
                for k, v in b.items():
                    nv=r.get(k, Fraction(0))-f*v
                    if nv==0:
                        r.pop(k, None)
                    else:
                        r[k]=nv
        if r:
            p=min(r)
            basis.append(r)
            pivots.append(p)
    return basis

class _Span:
    # Incrementally maintained span of used combination vectors (coordinates
    # over the row basis), to keep the chosen equations independent.
    def __init__(self, dim):
        self.dim=dim
        self.vecs=[]
        self.piv=[]
    def reduce(self, w):
        w=list(w)
        for p, v in zip(self.piv, self.vecs):
            if w[p]!=0:
                f=w[p]/v[p]
                w=[a-f*b for a, b in zip(w, v)]
        return w
    def independent(self, w):
        return any(x!=0 for x in self.reduce(w))
    def add(self, w):
        r=self.reduce(w)
        p=next(i for i, x in enumerate(r) if x!=0)
        self.vecs.append(r)
        self.piv.append(p)
    def complement(self):
        # Coordinates not yet pivots: a vector with a component there is new.
        used=set(self.piv)
        return [i for i in range(self.dim) if i not in used]


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# One unknown by linear program
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

def _rationalize(x):
    for D in (10**3, 10**6, 10**9):
        yield [Fraction(v).limit_denominator(D) for v in x]

def _combine(w, basis):
    c={}
    for wi, r in zip(w, basis):
        if wi==0:
            continue
        for a, v in r.items():
            nv=c.get(a, Fraction(0))+wi*v
            if nv==0:
                c.pop(a, None)
            else:
                c[a]=nv
    return c

def _own_row(u, Au, B_rows, row_of, forced, span, basis):
    # A state's own balance, if it already isolates the state: the natural
    # direct solve, preferred over any combination the LP might find.
    i=row_of.get(u)
    if i is None:
        return None
    c=dict(B_rows[i])
    if not c:
        return None
    Au=set(Au)
    if any(a in forced for a in c):
        return None
    neg=sum(v for a, v in c.items() if a in Au)
    if any(v>0 for a, v in c.items() if a in Au) or neg>=0:
        return None
    if any(v<0 for a, v in c.items() if a not in Au):
        return None
    if not any(v>0 for a, v in c.items() if a not in Au):
        return None
    w=_coords(c, basis)
    if w is None or not span.independent(w):
        return None
    return w, c

def _sign_of_sum(e):
    # +1 / -1 if every term of the expanded sum has that sign, 0 if mixed or zero.
    terms=sympy.Add.make_args(sympy.expand(e))
    if e==0:
        return 0
    signs={1 if t.as_coeff_Mul()[0].is_positive else -1 for t in terms}
    return signs.pop() if len(signs)==1 else 0

def _compact(e):
    # A polynomial written with its most frequent symbol factored out, recursively:
    # no expansion, and a sum of positive terms stays one.
    f=sympy.factor(e)
    if not f.is_Add:
        return f
    terms=list(sympy.Add.make_args(sympy.expand(e)))
    count={}
    for t in terms:
        for x in t.free_symbols:
            count[x]=count.get(x, 0)+1
    best=max(sorted(count, key=str), key=lambda x: count[x], default=None)
    if best is None or count[best]<2:
        return sympy.Add(*terms)
    inner=[t for t in terms if t.has(best)]
    rest=[t for t in terms if not t.has(best)]
    return best*_compact(sympy.expand(sympy.Add(*inner)/best))+_compact(sympy.Add(*rest))

def _own_quadratic(u, atoms, B_rows, row_of, forced, span, basis):
    """A state's own balance, quadratic in the state once its denominators are
    cleared: a*u^2 = (P - N)*u + c with a, c, P, N positive sums. Returns
    (w, root) or None.

    The positive root is written 2*c/(sqrt((P - N)^2 + 4*a*c) + N - P), free of
    subtractions outside the square when P = 0, and (P + sqrt(P^2 + 4*a*c))/(2*a)
    when N = 0.
    """
    i=row_of.get(u)
    if i is None:
        return None
    c=B_rows[i]
    if not c or any(a in forced for a in c):
        return None
    w=_coords(c, basis)
    if w is None or not span.independent(w):
        return None
    expr=sympy.Add(*[sympy.Rational(v.numerator, v.denominator)*atoms[a] for a, v in c.items()])
    # the denominators are products of positive sums: the numerator keeps the sign
    num, _=sympy.fraction(sympy.together(expr))
    try:
        poly=sympy.Poly(num, u)
    except sympy.PolynomialError:
        return None
    if poly.degree()!=2 or any(s.has(u) for s in poly.free_symbols_in_domain):
        return None
    a2, a1, a0=poly.all_coeffs()
    s2, s0=_sign_of_sum(a2), _sign_of_sum(a0)
    if s2==0 or s0==0 or s2==s0:
        return None
    if s2>0:
        a2, a1, a0=-a2, -a1, -a0
    A, C=_compact(-a2), _compact(a0)
    terms=sympy.Add.make_args(sympy.expand(a1))
    P=_compact(sympy.Add(*[t for t in terms if t.as_coeff_Mul()[0].is_positive]))
    N=_compact(-sympy.Add(*[t for t in terms if not t.as_coeff_Mul()[0].is_positive]))
    if P==0 and N==0:
        return w, sympy.sqrt(C/A)
    root=sympy.sqrt((P-N)**2+4*A*C)
    if N==0:
        return w, (P+root)/(2*A)
    return w, 2*C/(root+N-P)

def _coords(c, basis):
    # Coordinates of the linear form c in the row basis (exact), or None.
    c=dict(c)
    w=[Fraction(0)]*len(basis)
    for k, b in enumerate(basis):
        p=min(b)
        if p in c:
            f=c[p]/b[p]
            w[k]=f
            for a, v in b.items():
                nv=c.get(a, Fraction(0))-f*v
                if nv==0:
                    c.pop(a, None)
                else:
                    c[a]=nv
    return None if c else w

def _solve_unknown(u, Au, basis, Rnp, forced, span):
    """Combination w of the balances that isolates u, or None.

    c = w^T R: atoms of u (Au) carry c <= 0 summing to -1, every other atom
    c >= 0, atoms in `forced` c = 0 (they would close a reference cycle).
    Minimises the positive side, which keeps the equation sparse. The float
    optimum is rationalised and re-checked exactly; it must be independent of
    the equations already used.
    """
    r, A=Rnp.shape
    Au=set(Au)
    forced=set(forced)
    others=[a for a in range(A) if a not in Au and a not in forced]
    if not others:
        return None
    A_ub=numpy.vstack([Rnp[:, sorted(Au)].T, -Rnp[:, others].T])
    b_ub=numpy.zeros(A_ub.shape[0])
    free_u=sorted(Au-set(forced))
    if not free_u:
        return None
    eq=[Rnp[:, free_u].sum(axis=1)]
    beq=[-1.0]
    for a in forced:
        eq.append(Rnp[:, a])
        beq.append(0.0)
    obj=Rnp[:, others].sum(axis=1)
    extra=[None]+[(i, s) for i in span.complement() for s in (1, -1)]
    for ex in extra:
        Aeq=numpy.array(eq)
        Beq=numpy.array(beq)
        Aub, Bub=A_ub, b_ub
        if ex is not None:
            # force a component outside the span of the used equations
            i, s=ex
            red=numpy.zeros(r)
            red[i]=-s
            Aub=numpy.vstack([A_ub, red])
            Bub=numpy.append(b_ub, -1e-4)
        res=linprog(obj, A_ub=Aub, b_ub=Bub, A_eq=Aeq, b_eq=Beq,
                    bounds=[(None, None)]*r, method='highs')
        if not res.success:
            continue
        for w in _rationalize(res.x):
            c=_combine(w, basis)
            if any(c.get(a, 0)>0 for a in Au) or sum(c.get(a, 0) for a in Au)>=0:
                continue
            if any(c.get(a, 0)<0 for a in c if a not in Au):
                continue
            if any(c.get(a, 0)!=0 for a in forced):
                continue
            if not any(v>0 for a, v in c.items() if a not in Au):
                continue
            if not span.independent(w):
                continue
            return w, c
    return None

def _formula(u, c, Au, atoms):
    # scale to coprime integers: the LP normalisation is arbitrary
    from math import gcd
    L=1
    for v in c.values():
        L=L*v.denominator//gcd(L, v.denominator)
    G=0
    for v in c.values():
        G=gcd(G, int(v*L))
    c={a: v*L/G for a, v in c.items()}
    num=sympy.Add(*[sympy.Rational(v.numerator, v.denominator)*atoms[a]
                    for a, v in c.items() if a not in Au and v>0])
    den=sympy.Add(*[sympy.Rational(-v.numerator, v.denominator)*(atoms[a]/u)
                    for a, v in c.items() if a in Au and v<0])
    return num/den


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Search
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

def _reaches(start, target, deps):
    # Is `target` reachable from `start` along solution references?
    seen=set()
    stack=[start]
    while stack:
        x=stack.pop()
        if x==target:
            return True
        if x in seen:
            continue
        seen.add(x)
        stack.extend(deps.get(x, ()))
    return False

def _group_unknown(order, cands, nonlin, atoms, atom_syms, basis, Rnp, span,
                   sol, deps, B_rows, row_of, state_set):
    # One side of an unsolved state's balance, every atom of it carrying its
    # own unsolved rate constant. The fluxes of that side share their sum by
    # new free flux ratios, flux_g = r_g * flux_g0 (v1.3's r_*), and the sum
    # fixes k_g0.
    for y in order:
        if y in sol or y not in row_of:
            continue
        row=B_rows[row_of[y]]
        for sign in (-1, 1):
            side=[a for a, v in row.items() if (v>0)==(sign>0)]
            if len(side)<2:
                continue
            params=[]
            for a in side:
                ps=[p for p in atom_syms[a] if p in cands and p not in state_set
                    and p not in sol and a in cands[p]]
                if not ps:
                    params=None
                    break
                params.append(min(ps, key=lambda p: (len(cands[p]), str(p))))
            if params is None:
                continue
            G=list(dict.fromkeys(params))
            if len(G)<2:
                continue
            Au=set()
            forced=set()
            for p in G:
                Au|=set(cands[p])
                forced|=nonlin[p]
                bad={v for v in sol if _reaches(v, p, deps)}
                forced|={a for a in range(len(atoms)) if atom_syms[a] & bad}
            res=_solve_unknown(None, Au, basis, Rnp, forced, span)
            if res is None:
                continue
            w, c=res
            # flux of each rate constant in this equation, per unit rate constant;
            # one whose flux cancels between the balances is not part of it
            phi={p: sympy.Add(*[sympy.Rational(-c[a].numerator, c[a].denominator)*atoms[a]/p
                                for a in cands[p] if c.get(a, 0)<0]) for p in G}
            G=[p for p in G if phi[p]!=0]
            if not G:
                continue
            g0=G[0]
            rho={p: sympy.Symbol('r_'+str(y)+'_'+str(i)) for i, p in enumerate(G) if i>0}
            return g0, [(p, rho[p]) for p in G[1:]], w, c, Au, phi
    return None

def _greedy(order, cands, nonlin, atoms, basis, Rnp, r, check_time, B_rows, row_of,
            quadratic=(), roots=None):
    """Pick r unknowns, in `order`, each the first one an LP can isolate.

    States in `quadratic` may take the positive root of their own balance; they
    are tried before the first rate constant, so a root is preferred to a pivot;
    the states that took one are appended to `roots`.
    Returns (solutions, None) or (None, blockers): the solved unknowns whose
    references alone kept a remaining candidate from an equation.
    """
    span=_Span(len(basis))
    sol={}
    deps={}
    solved=[]
    ratios=[]
    state_set=set(row_of)
    unknowns=set(cands)|state_set
    atom_syms=[at.free_symbols for at in atoms]

    def root_of_some_state():
        for y in quadratic:
            if y in sol:
                continue
            bad={v for v in sol if _reaches(v, y, deps)}
            forced={a for a in range(len(atoms)) if atom_syms[a] & bad}
            res=_own_quadratic(y, atoms, B_rows, row_of, forced, span, basis)
            if res is not None:
                return y, res
        return None

    while len(solved)<r:
        if check_time():
            return None, None
        picked=None
        quad=None
        tried_quad=not quadratic
        for u in order:
            if u in sol:
                continue
            if not tried_quad and u not in state_set:
                tried_quad=True
                quad=root_of_some_state()
                if quad is not None:
                    break
            # a solved symbol whose solution leads back to u must not appear
            bad={v for v in sol if _reaches(v, u, deps)}
            forced={a for a in range(len(atoms)) if atom_syms[a] & bad} | nonlin[u]
            res=_own_row(u, cands[u], B_rows, row_of, forced, span, basis)
            if res is None:
                res=_solve_unknown(u, cands[u], basis, Rnp, forced, span)
            if res is not None:
                picked=(u, res)
                break
        if picked is None and quad is None and not tried_quad:
            quad=root_of_some_state()
        if quad is not None:
            u, (w, f)=quad
            span.add(w)
            sol[u]=f
            deps[u]={s for s in f.free_symbols if s in unknowns}
            solved.append(u)
            if roots is not None:
                roots.append(u)
            continue
        if picked is None:
            # no single unknown: let one side of a balance share a scale
            grp=_group_unknown(order, cands, nonlin, atoms, atom_syms, basis, Rnp,
                               span, sol, deps, B_rows, row_of, state_set)
            if grp is not None:
                g0, others, w, c, Au, phi=grp
                num=sympy.Add(*[sympy.Rational(v.numerator, v.denominator)*atoms[a]
                                for a, v in c.items() if a not in Au and v>0])
                span.add(w)
                # flux_g0 * (1 + sum r) = num; flux_g = r_g * flux_g0
                sol[g0]=num/(phi[g0]*(1+sympy.Add(*[rho for _, rho in others])))
                deps[g0]={s for s in sol[g0].free_symbols if s in unknowns}
                solved.append(g0)
                for p, rho in others:
                    sol[p]=rho*g0*phi[g0]/phi[p]
                    deps[p]={s for s in sol[p].free_symbols if s in unknowns}
                    ratios.append((str(rho), str(p), str(g0)))
                continue
            # the smallest set of solutions whose references alone keep some
            # candidate from an equation: leave exactly those for later
            blockers=None
            for u in order:
                if u in sol:
                    continue
                bad={v for v in sol if _reaches(v, u, deps)}
                if not bad or (blockers is not None and len(bad)>=len(blockers)):
                    continue
                if _solve_unknown(u, cands[u], basis, Rnp, nonlin[u], span) is not None:
                    blockers=bad
            blockers=blockers or set()
            # the balances the used ones do not span, named by their state
            state_of={i: y for y, i in row_of.items()}
            rest=_Span(len(basis))
            for w in span.vecs:
                rest.add(w)
            left=[]
            for i, row in enumerate(B_rows):
                w=_coords(row, basis) if row else None
                if w is not None and rest.independent(w):
                    rest.add(w)
                    left.append(str(state_of[i]))
            _STUCK[0]=(len(solved), left, [str(u) for u in solved])
            return None, blockers
        u, (w, c)=picked
        f=_formula(u, c, set(cands[u]), atoms)
        span.add(w)
        sol[u]=f
        deps[u]={s for s in f.free_symbols if s in unknowns}
        solved.append(u)
    extra=[(u, f) for u, f in sol.items() if u not in solved]
    return ([(u, sol[u]) for u in solved]+extra, ratios), None


# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# Entry point, same interface as v1.3
# –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

def Alyssa(filename,
           injections=[],
           givenCQs=[],
           neglect=[],
           sparsifyLevel=0,
           outputFormat='R',
           testSteady='fast',
           walltime=0,
           simplify=True,
           solveQuadratic=False,
           positive=True,
           branches=False,
           priority=[],
           verbose=True,
           volumes={}):
    global _VERBOSE
    _VERBOSE=2 if verbose=='full' else int(bool(verbose))
    del _NOTES[:]
    _BLANK[0]=False
    _STUCK[0]=None
    t0=time.time()
    check_time=lambda: walltime>0 and time.time()-t0>walltime
    filename=str(filename)
    injections=[str(s) for s in injections]
    neglect={str(s) for s in neglect}
    priority=[str(s) for s in priority]
    if positive is not True:
        raise ValueError('version 1.4 requires positive = TRUE')
    if outputFormat!='R':
        raise ValueError('version 1.4 writes outputFormat = "R" only')

    states, fluxes, stoich=_read_model(filename, injections)
    ODE_states=list(states)
    ODE=[sympy.Add(*[col[i]*f for f, col in zip(fluxes, stoich)]) for i in range(len(states))]
    states, fluxes, stoich, zero=_zero_states(states, fluxes, stoich,
                                              _v13.ParseVolumes(volumes))
    atoms, B=_atoms(states, fluxes, stoich)
    basis=_row_basis(B, len(atoms))
    r=len(basis)
    Rnp=numpy.array([[float(row.get(a, 0)) for a in range(len(atoms))] for row in basis])
    cands_all, nonlin=_candidates(states, atoms, neglect)
    # a given conserved quantity keeps one of its states free, the first by
    # default; if that state's balance is left, the next one is tried
    cq_groups=[]
    for cq in givenCQs:
        lhs=str(cq).split('=')[0]
        g=[sympy.Symbol(nm) for nm in re.findall(r'[A-Za-z_][A-Za-z0-9_]*', lhs)]
        g=[x for x in dict.fromkeys(g) if x in states]
        if g:
            cq_groups.append(g)
    cq_pick=[0]*len(cq_groups)

    def keep_free():
        free=[]
        for g, k in zip(cq_groups, cq_pick):
            rest=[x for x in g if x not in free]
            if k<len(rest):
                free.append(rest[k])
        return free
    cq_free=keep_free()
    cands={x: v for x, v in cands_all.items() if x not in cq_free}
    _say('steadyStates: '+str(len(states))+' states'+
         (', '+str(len(zero))+' zero a priori' if zero else '')+
         ', '+str(r)+' independent balances over '+str(len(atoms))+' flux terms')
    _say('')

    known={str(s) for s in cands_all}
    unknown=[p for p in priority if p not in known]
    if unknown:
        _note('priority ignores '+', '.join(unknown)+' (no state or rate constant)')
    state_set=set(states)
    prio_rank={p: i for i, p in enumerate(priority)}

    # A state whose flux terms occur in no other balance (an mRNA, a protein)
    # is safe to solve first. One that hands mass on (a precursor, a chain, a
    # receptor pool) would put its rate constants into its solution, and a
    # later balance may need exactly those: such states come after the rate
    # constants.
    atom_rows={}
    for i, row in enumerate(B):
        for a in row:
            atom_rows.setdefault(a, set()).add(i)
    isolated={s for i, s in enumerate(states) if all(atom_rows[a]=={i} for a in B[i])}

    def ordering(defer):
        # user priority first, then isolated states, rate constants, coupled
        # states, and last whatever blocked a former run
        def key(s):
            tier=0 if s in isolated else (2 if s in state_set else 1)
            return (0 if str(s) in prio_rank else 1, prio_rank.get(str(s), 0),
                    s in defer, tier, str(s))
        return sorted(cands, key=key)

    ratios=[]
    solution=None
    blockers=None
    no_root=set()
    attempt=0
    for rotation in range(1+sum(len(g) for g in cq_groups)):
        defer=set()
        no_root=set()
        prev=set()
        for _ in range(4*len(cands)+len(states)+1):
            attempt+=1
            order=ordering(defer)
            quadratic=[]
            if solveQuadratic:
                quadratic=[s for s in order if s in state_set]+\
                          sorted((s for s in states if s not in cands and str(s) not in neglect
                                  and s not in cq_free), key=str)
                quadratic=[s for s in quadratic if s not in no_root]
            roots=[]
            result, blockers=_greedy(order, cands, nonlin, atoms, basis, Rnp, r, check_time,
                                     B, {s: i for i, s in enumerate(states)}, quadratic, roots)
            solution=None if result is None else result[0]
            if solution is not None or blockers is None:
                ratios=result[1] if result is not None else []
                break
            n_used, left, used=_STUCK[0]
            head='  attempt '+str(attempt)+': '+str(n_used)+'/'+str(r)+' balances, left '+', '.join(left)
            # full trace: which unknowns this attempt solved differently
            gained=[u for u in used if u not in prev]
            lost=sorted(prev-set(used))
            if prev and (gained or lost):
                _trace('    '+' '.join(['+'+u for u in gained]+['-'+u for u in lost]))
            prev=set(used)
            # a kept state of a conserved quantity is changed before anything else
            left_syms={sympy.Symbol(x) for x in left}
            if any(cq_pick[i]+1<len(g) and set(g) & left_syms for i, g in enumerate(cq_groups)):
                _say(head)
                break
            # a root spends a balance a pivot would not: those that blocked go first
            if roots:
                banned=[u for u in roots if u in blockers] or roots[-1:]
                no_root|=set(banned)
                _say(head+'; next: '+', '.join(str(u) for u in banned)+' not as a root')
                continue
            new=blockers-defer
            if not new:
                _say(head)
                break
            defer|=new
            _say(head+'; next: '+', '.join(sorted(str(s) for s in new))+' last')
        if solution is not None or blockers is None or check_time():
            break
        # a conserved quantity whose kept state's balance is left moves on
        left_syms={sympy.Symbol(x) for x in _STUCK[0][1]}
        moved=next((i for i, g in enumerate(cq_groups)
                    if cq_pick[i]+1<len(g) and set(g) & left_syms), None)
        if moved is None:
            break
        old=cq_free[moved]
        cq_pick[moved]+=1
        cq_free=keep_free()
        cands={x: v for x, v in cands_all.items() if x not in cq_free}
        _say('  -> conserved quantity '+' + '.join(str(x) for x in cq_groups[moved])+': '+
             str(cq_free[moved])+' free instead of '+str(old))
        _say('')
    if no_root and solution is not None:
        _note('no root for '+', '.join(sorted(str(s) for s in no_root))+
              ': it would take a balance another unknown needs')
    if attempt>1:
        _say('')
    if solution is None:
        if check_time():
            print('No positive steady state found: walltime exceeded.', flush=True)
        elif _STUCK[0] is not None:
            left=_STUCK[0][1]
            print('No positive steady state found: the balance'+('s' if len(left)>1 else '')+
                  ' of '+', '.join(left)+' '+('are' if len(left)>1 else 'is')+
                  ' left without a positive unknown.', flush=True)
            print('  neglect, priority'+('' if solveQuadratic else ' or solveQuadratic = TRUE')+
                  ' change which unknowns are solved.', flush=True)
        else:
            print('No positive steady state found.', flush=True)
        _print_notes()
        return 0

    n_state=sum(1 for u, _ in solution if u in state_set)
    if ratios:
        _say('  '+str(len(ratios))+' new ratio parameter(s): '+
             ', '.join(rho for rho, _, _ in ratios))
    _say('  '+str(len(solution))+' unknowns solved: '+str(n_state)+' states, '+
         str(len(solution)-n_state)+' rate constants')

    # lazy equations, definitions before uses
    solved_names={str(u) for u, _ in solution}
    eqOut=[str(u)+' = '+str(f) for u, f in solution]
    eqOut+=[str(s)+' = '+str(s) for s in states if str(s) not in solved_names]
    eqOut=_v13._topo_sort_eqs(eqOut)

    if simplify:
        full=isinstance(simplify, str) and simplify.lower()=='full'
        _say('  simplifying '+str(len(solution))+' expressions ...')
        out=[]
        for eq in eqOut:
            ls, rs=eq.split(' = ', 1)
            # a root is built compact; simplifying it expands the discriminant
            if ls!=rs and 'sqrt' not in rs:
                if check_time():
                    _note('walltime exceeded while simplifying, the rest is left as it is')
                    simplify=False
                else:
                    rs=_v13._finalSimplify(rs, full)
            out.append(ls+' = '+rs)
        eqOut=out

    if testSteady in ('fast', 'modp', 'exact', 'T'):
        bad=_v13._steady_test_fast([str(o) for o in ODE], eqOut, zero)
        if bad:
            _note('not steady: '+', '.join('d'+str(ODE_states[i])+'/dt' for i in bad))
        test_status='FAILED' if bad else 'passed (mod p)'
    else:
        test_status='skipped'

    # resolve: a dMod trafo substitutes all entries at once
    substituted=set()
    for i in range(len(eqOut)):
        ls, rs=eqOut[i].split(' = ', 1)
        if ls==rs:
            continue
        pat=re.compile(r'\b'+re.escape(ls)+r'\b')
        for j in range(i+1, len(eqOut)):
            ls2, rs2=eqOut[j].split(' = ', 1)
            new=pat.sub('('+rs+')', rs2)
            if new!=rs2:
                eqOut[j]=ls2+' = '+new
                substituted.add(j)
    # reprinted by sympy, which drops the redundant brackets; factor_terms only
    # pulls out common factors, cancel() would expand the nested sums
    if simplify:
        for j in sorted(substituted):
            ls2, rs2=eqOut[j].split(' = ', 1)
            eqOut[j]=ls2+' = '+str(sympy.factor_terms(parse_expr(rs2)))

    zero_names=[str(s) for s in zero]+[s for s in injections if s not in {str(z) for z in zero}]
    ret=[s+'=0' for s in zero_names]+[eq.replace(' = ', '=', 1) for eq in eqOut]

    if _VERBOSE>=1:
        print('I obtained the following equations:\n', flush=True)
        for s_ in zero_names:
            print('\t'+s_+' = 0\n', flush=True)
        for eq in eqOut:
            ls, rs=eq.split(' = ', 1)
            print('\t'+ls+' = "'+rs+'",\n', flush=True)
        print('Number of Species:  '+str(len(ODE_states)+len(injections)), flush=True)
        print('Number of Equations:  '+str(len(eqOut)+len(zero_names)), flush=True)
        print('Number of new introduced variables:  '+str(len(ratios)), flush=True)
        for rho, p, g0 in ratios:
            print('\t'+rho+' = flux('+p+') / flux('+g0+') > 0', flush=True)

    pivots=[str(u) for u, _ in solution if u not in state_set]
    free=[str(s) for s in states if str(s) not in solved_names]
    _say('')
    print('Steady state: '+str(len(solution))+' expressions, '+str(len(zero_names))+
          ' states at 0, test '+test_status+', '+str(round(time.time()-t0))+' s', flush=True)
    if pivots:
        print('  rate constants solved for: '+', '.join(pivots), flush=True)
    if free:
        print('  free states: '+', '.join(free)+
              (' ('+', '.join(str(s) for s in cq_free)+' by conserved quantities)' if cq_free else ''),
              flush=True)
    if ratios:
        print('  new flux ratios: '+', '.join(rho for rho, _, _ in ratios), flush=True)
    _print_notes()
    return ret

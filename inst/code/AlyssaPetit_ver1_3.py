# AlyssaPetit version 1.3
# Use with python 3.x
#
# Based on v1.2. Changes:
#   - Lazy solves: solutions are recorded, never substituted into F or earlier
#     equations; one textual resolution at output time (topological order).
#     Removes the expression blowup that made cancel()/GCD dominate runtime.
#   - Solved states enter the dependency graph as aliases, so cycle detection
#     sees dependencies through recorded solutions.
#   - Lock guard: a direct solve whose locked rate constants would strand a
#     cycle state is skipped up front instead of found by rollback later.
#   - Affinity test by polynomial degree instead of cancel(); rejections
#     cached while F is unchanged; single-row SM*F products.
#   - Steady-state test evaluates mod p by environment extension instead of
#     symbolic substitution.

import numpy
import sympy
from sympy import Matrix, simplify as _simplify, expand, solve, cancel, posify, factor, fraction, oo
from numpy import shape, zeros, concatenate
from sympy.parsing.sympy_parser import parse_expr
from sympy.matrices import *
from sympy.matrices import matrix_multiply_elementwise
from scipy.optimize import linprog
import csv
import re
import time
import random
from random import shuffle

def LCS(s1, s2):
    m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]

def SolveSymbLES(A,b):
    dim=shape(A)[0]
    Asave=A[:]
    Asave=Matrix(dim, dim, Asave)
    #printmatrix(Asave)
    #print(b)
    determinant=Asave.det()
    if(determinant==0):
        #print('Determinant of LCL-calculation is zero! Try to specify LCLs yourself!')
        return([])
    result=[]
    for i in range(dim):
        A=Matrix(dim,dim,Asave)
        A.col_del(i)
        A=A.col_insert(i,b)
        result.append(_simplify(A.det()/determinant))
    
    return(result)

def CutStringListatSymbol(liste, symbol):
    out=[]    
    for el in liste:
        if(symbol in el):
            add=el.split(symbol)
        else:
            add=[el]
        out=out+add
    return(out)

def FillwithRanNum(M):
    dimx=len(M.row(0))
    dimy=len(M.col(0))
    ranM=zeros(dimy, dimx)
    parlist=[]
    ranlist=[]
    for i in M[:]:
        if(i!=0):
            if(str(i)[0]=='-'):
                parlist.append(str(i)[1:])
            else:
                parlist.append(str(i))
    parlist=list(set(parlist))
    for symbol in [' - ', ' + ', '*', '/', '(',')']:
        parlist=CutStringListatSymbol(parlist,symbol)
    parlist=list(set(parlist))
    temp=[]    
    for i in parlist:
        if(i!=''):
            if(not is_number(i)):
                temp.append(i)
                ranlist.append(random.random())
    parlist=temp
    for i in range(dimy):
        for j in range(dimx):
            ranM[i,j]=M[i,j]
            if(ranM[i,j]!=0):
                for p in range(len(parlist)):
                   ranM[i,j]=ranM[i,j].subs(parse_expr(parlist[p]),ranlist[p])
    return(ranM)

def FindLinDep(M, tol=1e-12):
    ranM=FillwithRanNum(M)
    Q,R=numpy.linalg.qr(ranM)
    for i in range(shape(R)[0]):
        for j in range(shape(R)[1]):
            if(abs(R[i,j]) < tol):
                R[i,j]=0.0
                
    LinDepList=[]
    for i in range(shape(R)[0]):
        if(R[i][i]==0):
            LinDepList.append(i)
    
    return(LinDepList)

def FindLCL(M, X):
    LCL=[]    
    LinDepList=FindLinDep(M)
    i=0
    counter=0
    deleted_rows=[]
    states=Matrix(X[:])
    while(LinDepList!=[]):
        i=LinDepList[0]
        testM=FillwithRanNum(M)
        rowliste=list(numpy.nonzero(testM[:,i])[0])
        colliste=[i]        
        for z in range(i):
            for k in rowliste:        
                for j in range(i):
                    jliste=list(numpy.nonzero(testM[:,j])[0])
                    if(k in jliste):
                        rowliste=rowliste+jliste
                        colliste=colliste+[j]
            rowliste=list(set(rowliste))
            colliste=list(set(colliste))
        rowliste.sort()
        colliste.sort()
        colliste.pop()        
        rowlisteTry=rowliste[0:(len(colliste))]
        vec=SolveSymbLES(M[rowlisteTry,colliste],M[rowlisteTry,i])
        shufflecounter=0
        while(vec==[] and shufflecounter < 100):
            shuffle(rowliste)
            shufflecounter=shufflecounter+1
            rowlisteTry=rowliste[0:(len(colliste))]
            vec=SolveSymbLES(M[rowlisteTry,colliste],M[rowlisteTry,i])
        if(shufflecounter==100):
            print('Problems while finding conserved quantities!',flush=True)
            return([],0)
        counter=counter+1
        try:
            mat=[states[l] for l in colliste]
            test=parse_expr('0')
            for v in range(0,len(vec)):
                test=test-parse_expr(str(vec[v]))*parse_expr(str(mat[v]))
        except:
            return([],0)
        partStr=str(test)+' + '+str(states[i])
        partStr=partStr.split(' + ')
        partStr2=[]
        for index in range(len(partStr)):
            partStr2=partStr2+partStr[index].split('-')
        partStr=partStr2
        if(len(partStr) > 1):        
            CLString=LCS(str(partStr[0]),str(partStr[1]))
            for ps in range(2,len(partStr)):
                CLString=LCS(CLString,str(partStr[ps]))
        else:
            CLString=str(partStr[0])
        if(CLString==''):
            CLString=str(counter)
        LCL.append(str(test)+' + '+str(states[i])+' = '+'total'+CLString)
        M.col_del(i)
        states.row_del(i)
        deleted_rows.append(i+counter-1)
        LinDepList=FindLinDep(M)
    return(LCL, deleted_rows)

def printmatrix(M):    
    lengths=[]
    for i in range(len(M.row(0))):
        lengths.append(0)
        for j in range(len(M.col(0))):
            lengths[i]=max(lengths[i],len(str(M.col(i)[j])))          
    string=''.ljust(5)
    string2=''.ljust(5)
    for j in range(len(M.row(0))):
        string=string+(str(j)).ljust(lengths[j]+2)
        for k in range(lengths[j]+2):        
            string2=string2+('-')        
    print(string)
    print(string2)
    for i in range(len(M.col(0))):
        string=str(i).ljust(4) + '['
        for j in range(len(M.row(0))):
            if(j==len(M.row(0))-1):
                string=string+str(M.row(i)[j]).ljust(lengths[j])
            else:
                string=string+(str(M.row(i)[j])+', ').ljust(lengths[j]+2)        
        print(string+']',flush=True)    
    return()
    
def printgraph(G):
    for el in G:
        print(el+': '+str(G[el]),flush==True)
    return()
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def FindSinkCluster(SM, eps=1e-8, M=1e4):
    # Structural mass-balance test: find c >= 0, c != 0 such that
    # c^T * SM <= 0 componentwise with at least one strictly negative entry.
    # Support of such c is a subset of states whose total mass monotonically
    # leaks out of the system -> all of them must be zero in steady state.
    # Generalises checkNegRows (which only catches single-row sinks).
    n_states=SM.rows
    n_flux=SM.cols
    if n_states==0 or n_flux==0:
        return []
    SM_np=numpy.array(SM.tolist(), dtype=float)
    A_ub=SM_np.T
    b_ub=numpy.zeros(n_flux)
    c_obj=SM_np.sum(axis=1)
    for i in range(n_states):
        bounds=[(0.0, M)]*n_states
        bounds[i]=(1.0, 1.0)
        res=linprog(c=c_obj, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs')
        if res.success and res.fun is not None and res.fun<-eps:
            return [j for j in range(n_states) if res.x[j]>eps]
    return []

def _zero_out_state(row_idx, SM, F, X, zeroStates):
    xrow=X[row_idx]
    zeroStates.append(xrow)
    counter=0
    for i in range(len(F)):
        if F[i-counter].subs(xrow, 1)!=F[i-counter]:
            if F[i-counter].subs(xrow, 0)==0:
                F.row_del(i-counter)
                SM.col_del(i-counter)
                counter=counter+1
            else:
                F[i-counter]=F[i-counter].subs(xrow, 0)
    X.row_del(row_idx)
    SM.row_del(row_idx)

def _flip_sign(cls):
    # Sign class under multiplication by -1.
    return {'+': '-', '-': '+'}.get(cls, cls)

def _mono_sign_definite(mono, positive_syms):
    # Monomial sign is set by its coefficient only if every symbolic factor is
    # known non-negative: declared positive, or raised to an even power.
    for base, exp in mono.as_powers_dict().items():
        if not base.free_symbols or base in positive_syms:
            continue
        if getattr(exp, 'is_even', False):
            continue
        return False
    return True

def _sign_class(expr, positive_syms=None):
    # Sign of a polynomial under positive_syms (None = all symbols positive):
    # '+'/'-' if all monomials share a provable sign, '0' if zero, else '+/-'.
    expr=sympy.expand(expr)
    if expr==0:
        return '0'
    has_pos=False
    has_neg=False
    for t in sympy.Add.make_args(expr):
        coeff, mono=t.as_coeff_Mul()
        if positive_syms is not None and not _mono_sign_definite(mono, positive_syms):
            return '+/-'
        if coeff.is_negative:
            has_neg=True
        elif coeff.is_positive:
            has_pos=True
        else:
            return '+/-'
    if has_pos and has_neg:
        return '+/-'
    if has_pos:
        return '+'
    if has_neg:
        return '-'
    return '0'

def _rational_sign_class(expr, positive_syms=None):
    # Sign of a rational function: cancel, then combine numerator/denominator
    # signs. Returns '+'/'-'/'+/-'/'0'.
    expr=cancel(expr)
    numer, denom=sympy.fraction(expr)
    ns=_sign_class(numer, positive_syms)
    ds=_sign_class(denom, positive_syms)
    if ns=='0':
        return '0'
    if ns=='+/-' or ds=='+/-':
        return '+/-'
    if (ns=='+' and ds=='+') or (ns=='-' and ds=='-'):
        return '+'
    if (ns=='+' and ds=='-') or (ns=='-' and ds=='+'):
        return '-'
    return '+/-'

def _normalizeSign(expr):
    # factor() normalises the sign of each irreducible factor, which can park a
    # -1 in front of a quotient whose value is positive.
    numer, denom=sympy.fraction(expr)
    if _sign_class(denom)=='-':
        return (-numer)/(-denom)
    return expr


def _exclusiveFluxPivots(SM, F, fluxpars, index, neglect):
    """Rate constants that can absorb row `index`'s ODE without touching another.

    The parameter must factor exactly one flux, whose column is confined to this
    row and which occurs in no other flux. Influx sides come first: they give the
    "production = consumption" form, which is a pure sum.
    """
    cands=[]
    for k in range(SM.cols):
        if SM[index,k]==0:
            continue
        fp=fluxpars[k]
        if fp is None or str(fp) in neglect:
            continue
        if any(SM[r,k]!=0 for r in range(SM.rows) if r!=index):
            continue
        if any(fp in F[j].free_symbols for j in range(len(F)) if j!=k):
            continue
        cands.append((0 if SM[index,k]>0 else 1, fp))
    cands.sort(key=lambda c: c[0])
    return [fp for _, fp in cands]


def _reportSignIndefinite(node, sol, sol_cls, positive_syms, solved,
                          blocked, tried, solveQuadratic):
    """Report a steady state that no pivot could make manifestly non-negative."""
    numer, denom=sympy.fraction(cancel(sol))
    print("",flush=True)
    print("    ======================================================",flush=True)
    print("    STEADY STATE IS NOT MANIFESTLY NON-NEGATIVE",flush=True)
    print("    ======================================================",flush=True)
    print(f"    State: {node}     sign class: {sol_cls}"
          f"  (numerator {_sign_class(numer, positive_syms)},"
          f" denominator {_sign_class(denom, positive_syms)})",flush=True)
    print("",flush=True)
    print(f"      {node} = {sol}",flush=True)
    print("",flush=True)
    negs=[t for t in sympy.Add.make_args(sympy.expand(numer))
          if t.as_coeff_Mul()[0].is_negative]
    if negs:
        print("    Negative contributions in the numerator:",flush=True)
        for t in negs[:6]:
            via=sorted(str(sy) for sy in t.free_symbols if str(sy) in solved)
            tag="   <-- via already-solved "+", ".join(via) if via else ""
            print(f"      {t}{tag}",flush=True)
        if len(negs)>6:
            print(f"      ... and {len(negs)-6} more",flush=True)
        print("",flush=True)
    print("    Why:",flush=True)
    print(f"      The direct positive-solve pass could not certify {node}, so it",flush=True)
    print("      fell through to here. Upstream substitutions then left the",flush=True)
    print(f"      numerator a difference rather than a sum of production terms --",flush=True)
    print(f"      typically because a {node}-proportional consumption term was",flush=True)
    print(f"      replaced by a constant sink. {node} > 0 is therefore a",flush=True)
    print("      constraint on the parameters, not an identity, and the",flush=True)
    print("      expression changes sign over parameter space.",flush=True)
    print("",flush=True)
    print("    What you can do:",flush=True)
    _n=[0]
    def nextItem():
        _n[0]+=1
        return str(_n[0])
    if tried:
        print(f"      {nextItem()}. These rate pivots were tried and are"
              f" sign-indefinite too: {tried}",flush=True)
    if blocked:
        print(f"      {nextItem()}. Remove from 'neglect' to free a rate pivot:"
              f" {blocked}",flush=True)
    if not solveQuadratic:
        print(f"      {nextItem()}. solveQuadratic = TRUE -- admits closed-form"
              " positive roots at",flush=True)
        print("         the price of sqrt terms.",flush=True)
    print(f"      {nextItem()}. priority = ... to steer which state or rate is"
          " resolved first.",flush=True)
    print(f"      {nextItem()}. positive = FALSE to stop requiring sign-definiteness."
          " Note this",flush=True)
    print("         also changes which pivots get certified, so the whole",flush=True)
    print("         resolution order may differ.",flush=True)
    print("    ======================================================",flush=True)


def _try_positive_direct_solve(SM, F, X, positive_syms=None, neglect=(),
                               rejected=None, aliases=None):
    # Find a state y whose own ODE is (a) linear in y and (b) has a
    # structurally positive closed-form steady-state solution y = -In/Out
    # under positive_syms (None = all symbols positive). Returns
    # (row_idx, y, sol) or None. Used BEFORE cycle breaking to resolve the
    # parts of the network that don't actually need an r_* helper -- this
    # avoids leaking negative r_* contributions into downstream states whose
    # own ODE would otherwise have been positive by construction.
    # States in `neglect` are skipped: their ODE must be spent on a flux
    # pivot so they stay free parameters of the trafo.
    #
    # `rejected` (a set of state names) caches refusals: F does not change
    # between solves (lazy scheme), so a refused row stays refused until the
    # next cycle break. `aliases` are the recorded solutions' references; a
    # state reached by its own equation through a solved state is a cycle in
    # disguise and left to cycle breaking.
    X_names={str(X[j]) for j in range(len(X))}
    for i in range(len(X)):
        y=X[i]
        nm=str(y)
        if nm in neglect or (rejected is not None and nm in rejected):
            continue
        eq=sympy.S.Zero
        row=SM.row(i)
        for k in range(SM.cols):
            if row[k]!=0:
                eq=eq+row[k]*F[k]
        if eq==0:
            if rejected is not None: rejected.add(nm)
            continue
        if aliases:
            eq_names={str(s) for s in eq.free_symbols}
            if nm in _alias_closure(eq_names & set(aliases), aliases, X_names):
                if rejected is not None: rejected.add(nm)
                continue
        # Affinity in y by polynomial degree -- no GCD. together() only
        # collects; a y-bearing denominator means the row is genuinely
        # rational in y and a direct solve does not apply.
        num, den=sympy.fraction(sympy.together(eq))
        if den.has(y):
            if rejected is not None: rejected.add(nm)
            continue
        try:
            poly=sympy.Poly(num, y)
        except sympy.PolynomialError:
            if rejected is not None: rejected.add(nm)
            continue
        if poly.degree()!=1:
            if rejected is not None: rejected.add(nm)
            continue
        Out_n, In_n=poly.all_coeffs()
        In_cls=_rational_sign_class(In_n/den, positive_syms)
        Out_cls=_rational_sign_class(Out_n/den, positive_syms)
        if In_cls=='0' or Out_cls=='0':
            if rejected is not None: rejected.add(nm)
            continue
        if (In_cls=='+' and Out_cls=='-') or (In_cls=='-' and Out_cls=='+'):
            sol=cancel(-In_n/Out_n)
            return (i, y, sol)
        if rejected is not None: rejected.add(nm)
    return None

def _simplify_with_sqrt(expr):
    # Decompose a rational expression with at most one sqrt subterm into
    # (P + Q*sqrt(D))/R, factor each piece independently, and recombine.
    # SymPy's generic simplify/factor treats sqrt(.) as an opaque atom and
    # therefore can't reduce expressions of this shape -- it just leaves the
    # numerator as a giant sum mixing rational and irrational monomials.
    # The decomposition below catches the structure the quadratic solver
    # actually produces and reliably compacts it.
    expr = cancel(expr)
    sqrts = [s for s in expr.atoms(sympy.Pow) if s.exp == sympy.S.Half]
    if not sqrts:
        n, d = fraction(expr)
        return factor(n)/factor(d)
    if len(sqrts) != 1:
        return expr
    sq = sqrts[0]
    try:
        poly = sympy.Poly(expr, sq)
    except sympy.PolynomialError:
        return expr
    if poly.degree() > 1:
        return expr
    coeffs = poly.all_coeffs()
    if len(coeffs) == 1:
        n, d = fraction(coeffs[0])
        return factor(n)/factor(d)
    Q_part, P_part = coeffs[0], coeffs[1]
    P_n, P_d = fraction(cancel(P_part))
    Q_n, Q_d = fraction(cancel(Q_part))
    R = sympy.lcm(P_d, Q_d)
    P_scaled = cancel(P_n * R / P_d)
    Q_scaled = cancel(Q_n * R / Q_d)
    D_factored = factor(sq.args[0])
    sq_factored = sympy.sqrt(D_factored)
    return (factor(P_scaled) + factor(Q_scaled)*sq_factored) / factor(R)

def _try_positive_quadratic_solve(SM, F, X, positive_syms=None, branches=False,
                                  neglect=(), rejected=None, aliases=None):
    # Resolve a state y whose ODE is quadratic in y (a*y^2 + b*y + c = 0,
    # normalised to sign(a)=+) as a closed-form root, avoiding sympy.solve().
    # sign(c)=- gives the unique positive root (-b + sqrt(disc))/(2a); with
    # branches=True, sign(c)=+ & sign(b)=- gives two positive roots, emitted
    # with a selector branch_y in {-1,+1}. Other patterns pivot instead.
    # Returns (row_idx, y, sol, branch) or None. States in `neglect` are
    # skipped; `rejected` and `aliases` as in _try_positive_direct_solve.
    X_names={str(X[j]) for j in range(len(X))}
    for i in range(len(X)):
        y=X[i]
        nm=str(y)
        if nm in neglect or (rejected is not None and nm in rejected):
            continue
        eq=sympy.S.Zero
        row=SM.row(i)
        for k in range(SM.cols):
            if row[k]!=0:
                eq=eq+row[k]*F[k]
        if eq==0:
            if rejected is not None: rejected.add(nm)
            continue
        if aliases:
            eq_names={str(s) for s in eq.free_symbols}
            if nm in _alias_closure(eq_names & set(aliases), aliases, X_names):
                if rejected is not None: rejected.add(nm)
                continue
        # SS denominators are positive sums, so zeros of eq and its numerator coincide.
        num, _den = sympy.fraction(sympy.together(eq))
        try:
            poly = sympy.Poly(num, y)
        except sympy.PolynomialError:
            if rejected is not None: rejected.add(nm)
            continue
        if poly.degree() != 2:
            if rejected is not None: rejected.add(nm)
            continue
        a, b, c = poly.all_coeffs()
        a_cls=_rational_sign_class(a, positive_syms)
        c_cls=_rational_sign_class(c, positive_syms)
        if a_cls=='-':   # normalise to sign(a)=+ (roots unchanged)
            a, b, c = -a, -b, -c
            a_cls, c_cls = '+', _flip_sign(c_cls)
        if a_cls!='+':
            if rejected is not None: rejected.add(nm)
            continue
        disc = b**2 - 4*a*c
        if c_cls=='-':
            return (i, y, cancel((-b + sympy.sqrt(disc)) / (2*a)), None)
        if branches and c_cls=='+' and _rational_sign_class(b, positive_syms)=='-':
            branch=sympy.Symbol('branch_'+str(y))
            return (i, y, cancel((-b + branch*sympy.sqrt(disc)) / (2*a)), branch)
        if rejected is not None: rejected.add(nm)
    return None

def checkNegRows(M):
    NegRows=[]
    if((M==Matrix(0,0,[])) | (M==Matrix(0,1,[])) | (M==Matrix(1,0,[]))):
        return(NegRows)
    else:        
        for i in range(len(M.col(0))):
            foundPos=False
            for j in range(len(M.row(i))):
                if(M[i,j]>0):
                    foundPos=True
            if(foundPos==False):
                NegRows.append(i)    
        return(NegRows)
    
def checkPosRows(M):
    PosRows=[]
    if((M==Matrix(0,0,[])) | (M==Matrix(0,1,[])) | (M==Matrix(1,0,[]))):
        return(PosRows)
    else: 
        for i in range(len(M.col(0))):
            foundNeg=False
            for j in range(len(M.row(i))):
                if(M[i,j]<0):
                    foundNeg=True
            if(foundNeg==False):
                PosRows.append(i)    
        return(PosRows)             

def _alias_closure(names, aliases, keep):
    # Expand `names` through `aliases` ({solved name -> referenced state
    # names}) and return the reachable subset of `keep` (the unsolved
    # states). A visited set guards against reference cycles.
    out=set()
    stack=list(names)
    seen=set()
    while stack:
        nm=stack.pop()
        if nm in seen:
            continue
        seen.add(nm)
        if nm in keep:
            out.add(nm)
        if nm in aliases:
            stack.extend(aliases[nm])
    return out

def _states_in_cycles(graph):
    # Names on at least one directed cycle: members of a nontrivial strongly
    # connected component, or self-looped. Iterative Tarjan.
    index={}
    low={}
    onstack=set()
    stack=[]
    result=set()
    counter=[0]
    for root in graph:
        if root in index:
            continue
        work=[(root, iter(graph.get(root, ())))]
        index[root]=low[root]=counter[0]; counter[0]+=1
        stack.append(root); onstack.add(root)
        while work:
            node, it=work[-1]
            advanced=False
            for nxt in it:
                if nxt not in graph:
                    continue
                if nxt not in index:
                    index[nxt]=low[nxt]=counter[0]; counter[0]+=1
                    stack.append(nxt); onstack.add(nxt)
                    work.append((nxt, iter(graph.get(nxt, ()))))
                    advanced=True
                    break
                elif nxt in onstack:
                    low[node]=min(low[node], index[nxt])
            if advanced:
                continue
            work.pop()
            if work:
                parent=work[-1][0]
                low[parent]=min(low[parent], low[node])
            if low[node]==index[node]:
                scc=[]
                while True:
                    w=stack.pop(); onstack.discard(w)
                    scc.append(w)
                    if w==node:
                        break
                if len(scc)>1:
                    result.update(scc)
                elif scc[0] in graph.get(scc[0], ()):
                    result.add(scc[0])
    return result

def DetermineGraphStructure(SM, F, X, neglect, aliases=None):
    # `aliases` maps an already-solved state to the state names its recorded
    # solution references; a dependency through a solved state is still a
    # dependency, so cycle detection stays sound under lazy solves.
    graph={}
    nrows=SM.rows
    ncols=SM.cols
    aliases=aliases or {}
    X_names={str(X[j]) for j in range(len(X))}
    # Precompute free_symbols per flux and nonzero columns per row
    F_syms=[F[k].free_symbols for k in range(ncols)]
    row_nz=[[k for k in range(ncols) if SM[i,k]!=0] for i in range(nrows)]
    for i in range(nrows):
        liste=[]
        xi=X[i]
        # ODE[i] = sum_k SM[i,k]*F[k], so symbols in ODE[i] are the union
        # of symbols in F[k] for all k where SM[i,k] != 0
        ode_syms=set()
        for k in row_nz[i]:
            ode_syms|=F_syms[k]
        ode_names={str(s) for s in ode_syms}
        via_alias=_alias_closure(ode_names & set(aliases), aliases, X_names)
        for j in range(len(X)):
            xj=X[j]
            nm_j=str(xj)
            if(xj in ode_syms):
                if(j==i):
                    # Self-dependence: check if nonlinear (degree >= 2 in any flux)
                    nonlinear=False
                    for k in row_nz[i]:
                        if(xi in F_syms[k]):
                            if(xi in F[k].diff(xi).free_symbols):
                                nonlinear=True
                                break
                    if(nonlinear or nm_j in via_alias):
                        liste.append(nm_j)
                else:
                    liste.append(nm_j)
            elif(nm_j in via_alias):
                # Reached through a solved state's recorded solution.
                liste.append(nm_j)
            else:
                if(j==i):
                    liste.append(nm_j)
        graph[str(X[i])]=liste
    for el in neglect:
        if(parse_expr(el) in X):
            if not el in graph:
                graph[el]=[el]
            else:
                if(el not in graph[el]):
                    graph[el].append(el)
    return(graph)

def find_all_simple_cycles(graph, max_cycles=20000):
    """Enumerate all simple directed cycles in `graph`.

    `graph` is a dict {node_name: [successor_names]}. Each returned cycle is
    a list of node names (not closed). A Johnson-style DFS that restricts
    paths to the start node's index-subgraph, so every cycle is emitted
    exactly once.

    Caps at `max_cycles` to guard against pathological enumerations. For the
    AlyssaPetit pipeline the practical bound is the graph size after
    cycle-breaking (e.g. TGFb with 24 states).
    """
    cycles=[]
    nodes=list(graph.keys())
    index={n: i for i, n in enumerate(nodes)}
    def backtrack(start, current, visited, path):
        for nxt in graph.get(current, []):
            if nxt == start:
                cycles.append(list(path))
                if len(cycles) >= max_cycles:
                    return True
            elif nxt in index and index[nxt] > index[start] and nxt not in visited:
                visited.add(nxt)
                path.append(nxt)
                if backtrack(start, nxt, visited, path):
                    return True
                path.pop()
                visited.remove(nxt)
        return False
    for s in nodes:
        if backtrack(s, s, {s}, [s]):
            print(f'   Warning: cycle enumeration capped at {max_cycles}.', flush=True)
            break
    return cycles


def _side_type(flux_col_counts):
    """Classify one flux side using the global column-nonzero counts.

    Mirrors `getType` from Severin Bang's Julia port (helperFunctions.jl:947).
    `flux_col_counts[i]` is the number of ODEs in which flux `i` occurs
    (i.e. the count of nonzero entries in its SM column). Returns:

      1 = single flux on this side, appears only in this ODE (safe).
      2 = multiple fluxes on this side, each only in its own ODE (safe).
      3 = otherwise -- resolving this side may create new cycles.
      None = no fluxes on this side (caller marks dontUseThisSide=1).
    """
    if len(flux_col_counts) == 0:
        return None
    if len(flux_col_counts) == 1:
        return 1 if flux_col_counts[0] == 1 else 3
    if all(c == 1 for c in flux_col_counts):
        return 2
    return 3


def printPriorityTable(rows, top=5):
    """Log the head of the priority table -- the candidates that were in
    contention for the cycle break about to happen, in ranked order."""
    print(f'   Priority table (top {min(top, len(rows))} of {len(rows)}):',flush=True)
    for rank, r in enumerate(rows[:top], start=1):
        side='CQ ' if r['isOutflux'] is None else ('out' if r['isOutflux'] else 'in ')
        fps=','.join(str(fp) for fp in r['fluxPars']) or '-'
        print(f"     {rank}. {r['species']:<18s} {side} prio={r['userRank']} "
              f"type={r['type']} risk={r['propagationRisk']} len={r['fluxLength']} "
              f"cyc={r['NoCycleOccur']} rhs={r['OccInRhs']}"
              + ('  UNUSABLE' if r['dontUseThisSide'] else '')
              + f"  [{fps}]",flush=True)


def genPriorityTable(SM, F, fluxpars, X, LCLs, SSgraph, neglect, locked_fluxpars=None,
                     priority=()):
    """Global priority table for cycle breaking / state resolution.

    Idea and sort criteria ported from Severin Bang's Julia reimplementation
    of AlyssaPetit (`helperFunctions.jl::genPriorityTable`). Instead of
    picking the best pair from the FIRST found cycle (as `GetBestPair` does),
    this builds a global table over ALL remaining states x both flux sides,
    then sorts by six keys so the truly best resolution is always on top:

        (dontUseThisSide asc, OccInCycles desc, type asc,
         fluxLength asc, NoCycleOccur desc, OccInRhs asc)

    Each row represents one resolvable state/side combination and carries:
      - spIndex, species   : which state this row is about
      - isOutflux          : True = outflux side, False = influx side,
                             None = CQ-type (no flux side involved)
      - type               : 0 CQ / 1 safe-single / 2 safe-multi / 3 may
                             create new cycles (see _side_type)
      - fluxLength         : number of fluxes on the chosen side
      - fluxParIDs, fluxPars: per-side flux indices and their flux parameters
      - dontUseThisSide    : 1 if any flux parameter is in `neglect` or the
                             side has no flux -- the side is unusable
      - userRank           : position of the first symbol of `priority`
                             matching this row (its species or one of its
                             flux parameters); len(priority) if none match
      - OccInCycles        : 1 if the state is in any simple cycle, else 0
      - NoCycleOccur       : number of simple cycles containing the state
      - OccInRhs           : number of OTHER ODE equations whose RHS
                             contains this species

    `priority` is a user-supplied sequence of symbol names (states and/or
    rate parameters). It outranks every structural key except
    `dontUseThisSide` and `OccInCycles`: the solver still refuses unusable
    sides and still only breaks states that actually sit in a cycle, but
    among those the user's order decides. Naming a state promotes breaking
    that state's row; naming a rate parameter promotes the row whose pivot
    side contains it.

    Returns `(rows, cycles)` where rows is the sorted list of candidate
    dicts and cycles is the list of simple cycles used to compute the
    OccInCycles / NoCycleOccur columns.
    """
    n = len(X)
    X_names = [str(X[i]) for i in range(n)]
    name_to_index = {nm: i for i, nm in enumerate(X_names)}

    # Simple-cycle enumeration on the steady-state dependency graph.
    cycles = find_all_simple_cycles(SSgraph)
    noCycleOccur = [0]*n
    for cyc in cycles:
        for nm in cyc:
            if nm in name_to_index:
                noCycleOccur[name_to_index[nm]] += 1
    occInCycles = [1 if c > 0 else 0 for c in noCycleOccur]

    # occInRhs[i] = number of ODE rows that mention species i on their RHS,
    # derived directly from SSgraph to avoid re-traversing symbolic flux
    # expressions (matches giveFluxesAndSteadyStates in the Julia port).
    occInRhs = [0]*n
    for _src, deps in SSgraph.items():
        for dep in deps:
            if dep in name_to_index:
                occInRhs[name_to_index[dep]] += 1

    # Global column-nonzero counts for all fluxes (how many ODEs each flux
    # participates in). Used by _side_type to distinguish safe from cycle-
    # creating resolutions.
    col_counts = [CountNZE(SM.col(c)) for c in range(SM.cols)]

    neglect_set = set(str(x) for x in neglect)
    # Flux parameters that were consumed by a direct-positive-solve earlier
    # and therefore must not be reparameterized by cycle-breaking -- doing so
    # would create a self-referential substitution (the fp's new expression
    # is derived from an F already containing a state solution that itself
    # contains fp). See the "Direct positive-solve pass" block in Alyssa()
    # for why.
    locked_set = set(str(x) for x in (locked_fluxpars or set()))
    prio_index = {str(nm): i for i, nm in enumerate(priority)}
    prio_none = len(prio_index)

    def _user_rank(*names):
        return min((prio_index[nm] for nm in names if nm in prio_index),
                   default=prio_none)

    _lcl_lhs = [parse_expr(lcl.split(' = ')[0]) for lcl in LCLs]

    def _state_in_active_cq(nm):
        parsed = parse_expr(nm)
        for lhs in _lcl_lhs:
            if lhs.subs(parsed, 1) != lhs:
                return True
        return False

    rows = []
    for i in range(n):
        nm = X_names[i]

        # CQ-type row: the state is in an active conservation law -- we can
        # drop its ODE entirely and recover the state via the CQ later.
        # Only emit when the state is actually in a cycle, otherwise we'd
        # eagerly consume LCLs that aren't blocking progress.
        if occInCycles[i] and _state_in_active_cq(nm):
            rows.append({
                'spIndex': i, 'species': nm, 'isOutflux': None,
                'type': 0, 'fluxLength': 0,
                'fluxParIDs': [], 'fluxPars': [],
                'dontUseThisSide': 0,
                'OccInCycles': occInCycles[i],
                'NoCycleOccur': noCycleOccur[i],
                'OccInRhs': occInRhs[i],
                'propagationRisk': 0,
                'userRank': _user_rank(nm),
            })
            continue

        row_SM = SM.row(i)
        for is_out in (False, True):
            flux_ids = []
            for c in range(SM.cols):
                v = row_SM[c]
                if is_out and v < 0:
                    flux_ids.append(c)
                elif (not is_out) and v > 0:
                    flux_ids.append(c)
            flux_pars = [fluxpars[c] for c in flux_ids]
            fc_counts = [col_counts[c] for c in flux_ids]
            typ = _side_type(fc_counts)
            dont = 0
            if typ is None:
                typ = 3
                dont = 1
            # A flux without its own rate constant (an irreducible sum) offers
            # nothing to solve for, so the whole side is unusable as a pivot.
            if any(fp is None for fp in flux_pars):
                dont = 1
            if any(str(fp) in neglect_set for fp in flux_pars):
                dont = 1
            # A locked flux disqualifies the WHOLE side, not just itself:
            # `type`/`fluxLength` are what certify a pivot as positive, and a
            # type-1 pivot solves the state's FULL ODE, so a hidden second flux
            # turns the solve into influx - other_outflux.
            if any(fp is not None and str(fp) in locked_set for fp in flux_pars):
                dont = 1
            # Propagation risk (Python-only extension, see sort comment):
            # sum over pivot-side flux columns of (col_count - 1), i.e.
            # the number of OTHER ODEs whose flux expressions get rewired
            # by this pivot's trafo. A zero-risk pivot (all col_counts==1)
            # is either type 1 or type 2 -- by construction it cannot
            # introduce sign conflicts downstream because each substituted
            # flux parameter only appears in the pivot's own ODE.
            prop_risk = sum(max(c - 1, 0) for c in fc_counts)
            rows.append({
                'spIndex': i, 'species': nm, 'isOutflux': is_out,
                'type': typ, 'fluxLength': len(flux_ids),
                'fluxParIDs': flux_ids, 'fluxPars': flux_pars,
                'dontUseThisSide': dont,
                'OccInCycles': occInCycles[i],
                'NoCycleOccur': noCycleOccur[i],
                'OccInRhs': occInRhs[i],
                'propagationRisk': prop_risk,
                'userRank': _user_rank(nm, *(str(fp) for fp in flux_pars)),
            })

    # Sort key extensions over Severin Bang's 6-key Julia ranking
    # (helperFunctions.jl::genPriorityTable). The Julia ranking has:
    #   (dontUseThisSide, -OccInCycles, type, fluxLength,
    #    -NoCycleOccur, OccInRhs).
    # We insert two extra keys to preserve dMod's manifestly-positive
    # rational output contract, which Julia's pipeline does not need:
    #
    #  (a) propagationRisk -- sum of (col_count - 1) over the chosen
    #      side's flux columns, i.e. how many OTHER ODEs get rewired by
    #      the trafo. A zero-risk side is exactly type 1 or type 2, in
    #      which case every substituted flux parameter appears only in
    #      the pivot's ODE and no sign mixing can occur downstream. Type
    #      3 sides always have propagationRisk >= 1; among them we
    #      prefer the smallest risk.
    #
    #  (b) userRank -- the caller's explicit `priority` order, inserted
    #      directly below OccInCycles so it decides among all rows that are
    #      legitimately breakable, but can never force an unusable side or
    #      a state outside every cycle.
    #
    #  (c) rowop_penalty -- a binary flag that demotes type-3 anz==1
    #      pivots. Those take the direct-solve + SM row-manipulation
    #      path, which INSERTS new flux terms (with inherited negative
    #      stoichiometry) into other ODEs. That is strictly worse than
    #      the type-3 anz>1 ratio-parameter trafo, which only REPLACES
    #      existing flux terms with positive rational substitutes.
    #
    # These keys are Python-side additions; neither appears in Severin
    # Bang's Julia port (the Julia output format and solver semantics
    # don't carry a positivity contract, so the original 6-key ranking
    # is sufficient there).
    def _rowop_penalty(r):
        return 1 if (r['type'] == 3 and r['fluxLength'] == 1) else 0
    rows.sort(key=lambda r: (
        r['dontUseThisSide'],
        -r['OccInCycles'],
        r['userRank'],
        r['type'],
        r['propagationRisk'],
        _rowop_penalty(r),
        r['fluxLength'],
        -r['NoCycleOccur'],
        r['OccInRhs'],
    ))
    return rows, cycles


def GetNegFluxParameters(SM, fluxpars, X, node):
    row=list(X).index(parse_expr(node))
    liste=[]
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]<0):
            liste.append(fluxpars[i])
    return(liste)

def GetPosFluxParameters(SM, fluxpars, X, node):
    row=list(X).index(parse_expr(node))
    liste=[]    
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]>0):
            liste.append(fluxpars[i])
    return(liste)
        
def GetType(node, fp, fluxpars, LCLs):
    for LCL in LCLs:
        ls=parse_expr(LCL.split(' = ')[0])
        if(ls.subs(parse_expr(node),1)!=ls):
            return(0)
    if(GetAppearances(fp, fluxpars)==1):
        if(GetDimension(node)==1):
            return(1)
        else:
            return(2)
    else:
        return(3)
        
def GetFluxParameter(flux, X):
    """Rate constant the whole flux is proportional to, or None.

    A flux parameter is a top-level factor that kills the flux when set to
    zero, so that `flux = fluxpar * prefactor` -- the form every trafo body
    below relies on. A flux that is an irreducible sum has no such factor and
    yields None.
    """
    if(flux.args==()):
        return(flux)
    for el in flux.args:
        if(el not in X and not is_number(str(el)) and flux.subs(el, 0)==0):
            return(el)
    return(None)

def GetAppearances(fp, fluxpars, SM):
    anz=0
    cols = [i for i, x in enumerate(fluxpars) if x == fp]
    #col=list(fluxpars).index(fp)
    for i in cols:
        for j in range(len(SM.col(i))):
            if(SM.col(i)[j]!=0):
                anz=anz+1
    return(anz)

def GetDimension(node, X, SM, getSign=False):
    row=list(X).index(parse_expr(node))
    anzminus=0
    anzappearminus=0
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]<0):
            anzappearminus=anzappearminus+CountNZE(SM.col(i))
            anzminus=anzminus+1
    anzplus=0
    anzappearplus=0
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]>0):
            anzappearplus=anzappearplus+CountNZE(SM.col(i))
            anzplus=anzplus+1
    if(not getSign):
        return(min(anzminus, anzplus))
    else:
        if(anzminus<anzplus or (anzminus==anzplus and anzappearminus<anzappearplus)):
            return(anzminus, "minus")
        else:
            return(anzplus, "plus")
            
def GetOutfluxes(node, X, SM, F, fluxpars):
    row=list(X).index(parse_expr(node))
    outsum=0
    out=[]
    fps=[]
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]<0):
            outsum=outsum-SM.row(row)[i]*F[i]
            out.append(-SM.row(row)[i]*F[i])
            fps.append(fluxpars[i])
    return(out, outsum, fps)

def GetInfluxes(node, X, SM, F, fluxpars):
    row=list(X).index(parse_expr(node))
    outsum=0
    out=[]
    fps=[]
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]>0):
            outsum=outsum+SM.row(row)[i]*F[i]
            out.append(SM.row(row)[i]*F[i])
            fps.append(fluxpars[i])
    return(out, outsum, fps)

def _resolvePivotSum(sumOpp, fluxes, fps, nenner, state):
    """Opposite-side flux sum with the pivot relations substituted into itself.

    Cycle removal distributes the opposite-side sum S over the pivot fluxes,
    fp_j * prefactor_j = S * r_j / nenner. S can contain one of those fp_j: an
    earlier removal replaces the column of the flux it pivoted on by the sum it
    was solved against, so a recycling loop (X -> X_int -> X) puts X's own
    internalisation rate back onto X's influx side. S is then an equation, not
    a value. It is linear in the fp_j, so substitute and solve once -- without
    this the emitted trafo defines fp_0 in terms of itself.
    """
    fps_sym=[parse_expr(str(fp)) for fp in fps]
    S_expr=sympy.sympify(sumOpp)
    shared=S_expr.free_symbols & set(fps_sym)
    if(not shared):
        return(sumOpp)
    S=sympy.Dummy('S_opp')
    repl={}
    for j, fp in enumerate(fps_sym):
        weight=1 if j==0 else parse_expr('r_'+state+'_'+str(j))
        repl[fp]=S*weight/(nenner*cancel(fluxes[j]/fp))
    root=solve(S_expr.subs(repl, simultaneous=True)-S, S)
    shared_str=sorted(str(s) for s in shared)
    if(len(root)!=1):
        print('   WARNING: '+str(shared_str)+' feed(s) back into the flux sum '
              'balancing '+str(state)+' and could not be resolved -- the '
              'equations stay recurrent. Please report this bug!',flush=True)
        return(sumOpp)
    print('   Recycled pivot rate constant(s) '+str(shared_str)+
          ' resolved in the flux sum balancing '+str(state),flush=True)
    return(cancel(root[0]))

def FindNodeToSolve(graph, SM=None, F=None, X=None):
    leaves=[el for el in graph if graph[el]==[]]
    if not leaves:
        return(None)
    if len(leaves)==1 or SM is None or F is None or X is None:
        return(leaves[0])
    # Among all solvable leaf nodes, pick the one whose ODE row has the
    # fewest symbolic operations. Solving a "small" equation first keeps
    # the recorded solutions compact. Only leaf rows are assembled -- no
    # full SM*F product.
    best=None
    best_cost=float('inf')
    for el in leaves:
        idx=list(X).index(parse_expr(el))
        expr=sympy.S.Zero
        for k in range(SM.cols):
            if SM[idx,k]!=0:
                expr=expr+SM[idx,k]*F[k]
        cost=expr.count_ops()
        if cost<best_cost:
            best_cost=cost
            best=el
    return(best)

def CountNZE(V):
    counter=0
    for v in V:
        if(v!=0):
            counter=counter+1
    return(counter)

# testSteady='fast' checks f_i(x_ss) == 0 by Schwartz-Zippel: evaluate each
# residual at random points in GF(p); a non-zero rational function vanishes
# there with probability <= deg/p, so a few points give near-certainty cheaply.
# Primes are == 3 mod 4 so a field sqrt is d**((p+1)//4) for a quadratic residue.
_MODP_PRIMES=(2147483647, 1000000007)

class _ResampleModp(Exception):
    pass
class _UnsupportedModp(Exception):
    pass

def _collect_exponent_syms(expr):
    # Symbols that appear inside a non-numeric exponent (e.g. Hill coefficients
    # in C3**nhill). These are assigned small integers so the exponent stays a
    # concrete integer during modular evaluation.
    syms=set()
    for pw in expr.atoms(sympy.Pow):
        if not pw.exp.is_number:
            syms|=pw.exp.free_symbols
    return syms

def _eval_exponent(e, env_int):
    if e.is_number:
        return int(e)
    return int(e.subs(env_int))

def _eval_modp(expr, env, p):
    # Evaluate expr at the point env (symbol -> int) in GF(p). Raises
    # _ResampleModp on a zero denominator / non-residue sqrt, _UnsupportedModp
    # on a node we do not model (caller falls back to the symbolic test).
    if expr.is_Integer:
        return int(expr) % p
    if expr.is_Rational:
        return (int(expr.p) * pow(int(expr.q), -1, p)) % p
    if expr.is_Symbol:
        if expr not in env:
            # Used before it is defined -- the equations are not resolvable in
            # the given order (a self-referential entry). Hand the residual to
            # the symbolic test instead of dying with a KeyError.
            raise _UnsupportedModp()
        return env[expr] % p
    if expr.is_Add:
        return sum(_eval_modp(a, env, p) for a in expr.args) % p
    if expr.is_Mul:
        r=1
        for a in expr.args:
            r=r*_eval_modp(a, env, p) % p
        return r
    if expr.is_Pow:
        base, e = expr.args
        if e.is_Rational and e.q == 2:
            d=_eval_modp(base, env, p)
            if d == 0:
                if int(e.p) > 0:
                    return 0
                raise _ResampleModp()
            if pow(d, (p-1)//2, p) != 1:
                raise _ResampleModp()
            s=pow(d, (p+1)//4, p)
            return pow(s, int(e.p) % (p-1), p)
        n=_eval_exponent(e, env)
        bv=_eval_modp(base, env, p)
        if n == 0:
            return 1
        if bv == 0:
            if n > 0:
                return 0
            raise _ResampleModp()
        return pow(bv, n % (p-1), p)
    if expr.is_number:
        v=sympy.nsimplify(expr)
        if v.is_Rational:
            return (int(v.p) * pow(int(v.q), -1, p)) % p
    raise _UnsupportedModp()

def _steady_test_fast(ODE, eqOut, zeroStates, trials=3, primes=_MODP_PRIMES):
    # Schwartz-Zippel steady-state check WITHOUT symbolic substitution: the
    # solved symbols are evaluated in dependency order (eqOut is
    # topologically sorted, definitions first) into the same environment,
    # then every residual is evaluated at that point. Returns the sorted
    # list of indices of ODEs with a provably nonzero residual, or None if
    # a node is unsupported / no valid sample point was found (caller falls
    # back to the exact symbolic test).
    pairs=[]
    for eq in eqOut:
        ls, rs=eq.split(' = ', 1)
        if ls==rs:
            continue
        pairs.append((parse_expr(ls), parse_expr(rs)))
    odes=[sympy.sympify(o) for o in ODE]
    all_exprs=[r for _, r in pairs]+odes
    exp_syms=set()
    for e in all_exprs:
        exp_syms|=_collect_exponent_syms(e)
    defined={l for l, _ in pairs}
    zero_syms=set(zeroStates)
    base=set()
    for e in all_exprs:
        base|=e.free_symbols
    base-=defined
    bad=set()
    for p in primes:
        good=0
        attempts=0
        while good < trials and attempts < trials*30:
            attempts+=1
            env={}
            for s in base:
                name=str(s)
                if s in zero_syms:
                    env[s]=0
                elif name.startswith('branch_'):
                    env[s]=random.choice([1, p-1])
                elif s in exp_syms:
                    env[s]=random.randint(2, 5)
                else:
                    env[s]=random.randint(1, p-1)
            try:
                for l, r in pairs:
                    env[l]=_eval_modp(r, env, p)
                vals=[_eval_modp(o, env, p) for o in odes]
            except _ResampleModp:
                continue
            except _UnsupportedModp:
                return None
            good+=1
            for i, v in enumerate(vals):
                if v % p != 0:
                    bad.add(i)
        if good == 0:
            return None
    return sorted(bad)

def _topo_sort_eqs(eqs):
    # Order 'lhs = rhs' entries so every definition precedes its uses -- the
    # textual resolution substitutes entry i into entries j > i, so a
    # reference must point backwards. Free-parameter entries (lhs == rhs)
    # define nothing. Stable: among ready entries the earliest wins. Falls
    # back to the given order if the entries reference each other cyclically.
    n=len(eqs)
    parts=[eq.split(' = ', 1) for eq in eqs]
    patterns=[re.compile(r'\b'+re.escape(ls)+r'\b') if ls!=rs else None
              for ls, rs in parts]
    for i, (ls, rs) in enumerate(parts):
        if patterns[i] is not None and patterns[i].search(rs):
            print('   WARNING: '+ls+' is defined in terms of itself -- no order '
                  'resolves it. Please report this bug!',flush=True)
    deps=[set() for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i==j or patterns[j] is None:
                continue
            if patterns[j].search(parts[i][1]):
                deps[i].add(j)
    order=[]
    emitted=set()
    while len(order) < n:
        pick=next((i for i in range(n)
                   if i not in emitted and deps[i] <= emitted), None)
        if pick is None:
            print('   Warning: steady-state equations reference each other '
                  'cyclically; output may stay recurrent.', flush=True)
            return eqs
        emitted.add(pick)
        order.append(pick)
    return [eqs[i] for i in order]


def _finalSimplify(rs, full):
    """One right-hand side through the final-simplification pipeline.

    `full` selects cancel + posify + simplify(ratio=oo) + factor(num)/factor(den);
    otherwise a single sympy.simplify. sqrt takes the (P + Q*sqrt(D))/R route
    either way.
    """
    e_raw=parse_expr(rs)
    if(e_raw.has(sympy.sqrt)):
        e_pos, reps = posify(e_raw)
        e_pos = _simplify_with_sqrt(e_pos)
        return str(e_pos.subs(reps))
    if(full):
        e_pos, reps=posify(cancel(e_raw))
        e_pos=_simplify(e_pos, ratio=oo, rational=True)
        n, d=fraction(e_pos)
        return str(_normalizeSign((factor(n)/factor(d)).subs(reps)))
    return str(_normalizeSign(_simplify(e_raw)))


def Alyssa(filename,
          injections=[],
          givenCQs=[],
          neglect=[],
          sparsifyLevel = 2,
          outputFormat='R',
          testSteady='T',
          walltime=0,
          simplify=True,
          solveQuadratic=False,
          positive=True,
          branches=False,
          priority=[]):
    filename=str(filename)
    _start_time = time.time()
    def _check_walltime():
        if walltime > 0 and time.time() - _start_time > walltime:
            return True
        return False
    # Positive symbols for the sign-based solvers: True = all (None sentinel),
    # False = none, or an explicit list of names.
    if positive is True or positive is False:
        positive_syms = None if positive else set()
    elif isinstance(positive, str):
        positive_syms = {parse_expr(positive)}
    elif isinstance(positive, (list, tuple, set)):
        positive_syms = set(parse_expr(str(s)) for s in positive)
    else:
        positive_syms = None if bool(positive) else set()
    # Sign-definiteness is only decidable when every symbol is known positive.
    _enforcePositive = (positive_syms is None)
    # States certified positive by a lazy solve join positive_syms (rolled
    # back with the snapshot below).
    _added_positive=set()
    # Names the solver must not spend on a state expression / flux pivot.
    # Kept as a set of plain strings: `str(sym) in neglect` is the only test
    # anywhere in the pipeline.
    neglect=set(str(n) for n in neglect)
    # Caller-supplied resolution order for the priority table (states and/or
    # rate parameters, most-preferred first). Validated against the model
    # once X and fluxpars exist -- see the check below the flux-parameter
    # extraction.
    priority=[str(pr) for pr in priority]
    file=csv.reader(open(filename), delimiter=',')
    print('Reading csv-file ...',flush=True)
    L=[]
    nrrow=0
    nrcol=0
    for row in file:
        nrrow=nrrow+1
        nrcol=len(row)
        L.append(row)
        
    nrspecies=nrcol-2
    
##### Remove injections  
    counter=0
    for i in range(1,len(L)):
        if(L[i-counter][1] in injections):
            L.remove(L[i-counter])
            counter=counter+1       
    
##### Define flux vector F
    # Track fluxes zeroed by an injection substitution; their injected states
    # are treated as exogenous (see the zero-flux handling below).
    F=[]
    flux_zeroed_by_injection=[]

    for i in range(1,len(L)):
        F.append(L[i][1])
        F[i-1]=F[i-1].replace('^','**')
        F[i-1]=parse_expr(F[i-1])
        zeroed_by_inj=False
        for inj in injections:
            new_expr=F[i-1].subs(parse_expr(inj),0)
            if new_expr!=F[i-1] and new_expr==0:
                zeroed_by_inj=True
            F[i-1]=new_expr
        flux_zeroed_by_injection.append(zeroed_by_inj)
    F=Matrix(F)
    #print(F)
##### Define state vector X
    X=[]
    X=L[0][2:]
    for i in range(len(X)):
        X[i]=parse_expr(X[i])               
    X=Matrix(X)
    #print(X)
    Xo=X.copy()
        
##### Define stoichiometry matrix SM
    SM=[]
    for i in range(len(L)-1):
    	SM.append(L[i+1][2:])        
    for i in range(len(SM)):
    	for j in range(len(SM[0])):
    		if (SM[i][j]==''):
    			SM[i][j]='0'
    		SM[i][j]=parse_expr(SM[i][j])    
    SM=Matrix(SM)
    SM=SM.T
    SMorig=SM.copy()

    
##### Check for zero fluxes
    # A flux zeroed by an injection (forcing) drops only the injected state rows
    # (the forcing is exogenous); co-reactants and products stay so NegRows /
    # PosRows / FindSinkCluster can still flag e.g. a complex left without a
    # source. Any other zero flux drops only its column.
    injection_syms=set()
    for inj in injections:
        try:
            injection_syms.add(parse_expr(inj))
        except Exception:
            pass
    icounter=0
    jcounter=0
    n_flux_orig=len(F)
    for i in range(n_flux_orig):
        idx=i-icounter
        if F[idx]==0:
            from_injection=flux_zeroed_by_injection[i] if i<len(flux_zeroed_by_injection) else False
            F.row_del(idx)
            if from_injection:
                col=SM.col(idx)
                rows_to_drop=[r for r in range(SM.rows) if col[r]!=0 and X[r] in injection_syms]
                for r in sorted(rows_to_drop, reverse=True):
                    X.row_del(r)
                    SM.row_del(r)
                    SMorig.row_del(r)
            SM.col_del(idx)
            SMorig.col_del(idx)
            icounter=icounter+1

    print('Removed '+str(icounter)+' fluxes that are a priori zero!',flush=True)
    #printmatrix(SM)
    #print(F)
    #print(X)
    #print(UsedRC)
#####Check if some species are zero and remove them from the system
    zeroStates=[]
    NegRows=checkNegRows(SM)
    PosRows=checkPosRows(SM)
    while True:
        while((NegRows!=[]) | (PosRows!=[])):
            if(NegRows!=[]):
                row=NegRows[0]
                zeroStates.append(X[row])
                counter=0
                for i in range(len(F)):
                    if(F[i-counter].subs(X[row],1)!=F[i-counter] and F[i-counter].subs(X[row],0)==0):
                        F.row_del(i-counter)
                        SM.col_del(i-counter)
                        counter=counter+1
                    else:
                        if(F[i-counter].subs(X[row],1)!=F[i-counter] and F[i-counter].subs(X[row],0)!=0):
                            F[i-counter]=F[i-counter].subs(X[row],0)
                X.row_del(row)
                SM.row_del(row)
            else:
                row=PosRows[0]
                zeroFluxes=[]
                for j in range(len(SM.row(row))):
                    if(SM.row(row)[j]!=0):
                        zeroFluxes.append(F[j])
                for k in zeroFluxes:
                    StateinFlux=[]
                    for state in X:
                        if(k.subs(state,1)!=k):
                            StateinFlux.append(state)
                    if(len(StateinFlux)==1):
                        zeroStates.append(StateinFlux[0])
                        row=list(X).index(StateinFlux[0])
                        counter=0
                        for i in range(len(F)):
                            if(F[i-counter].subs(X[row],1)!=F[i-counter]):
                                if(F[i-counter].subs(X[row],0)==0):
                                    F.row_del(i-counter)
                                    SM.col_del(i-counter)
                                else:
                                    F[i-counter]=F[i-counter].subs(X[row],0)
                                counter=counter+1
            NegRows=checkNegRows(SM)
            PosRows=checkPosRows(SM)
        # Row-local checks exhausted. Try structural sink-cluster detection:
        # find a subset C of states with c^T*SM <= 0 (and strictly negative in
        # some column) -- the combined mass of C leaks monotonically, so every
        # state in C must be zero in steady state. Classic example: TGFb alone
        # looks OK (binding/dissociation balance) but {TGFb, R1_TGFb} is a
        # sink because the complex degrades.
        sink=FindSinkCluster(SM)
        if not sink:
            break
        for row_idx in sorted(sink, reverse=True):
            _zero_out_state(row_idx, SM, F, X, zeroStates)
        NegRows=checkNegRows(SM)
        PosRows=checkPosRows(SM)
    #printmatrix(SM)
    #print(F)
    #print(X)
    nrspecies=nrspecies-len(zeroStates)
    if(nrspecies==0):
        print('All states are zero!',flush=True)
        return(0)
    else:
        if(zeroStates==[]):
            print('No states found that are a priori zero!',flush=True)
        else:
            print('These states are zero:',flush=True)
            for state in zeroStates:
                print('\t'+str(state),flush=True)
    
    nrspecies=nrspecies+len(zeroStates)

##### Identify linearities, bilinearities and multilinearities
    SMF=SM*F
    SMF_exp=[expand(e) for e in SMF]

    BLList=[]
    MLList=[]
    for i in range(len(SMF)):
        LHS=str(SMF_exp[i])
        LHS=LHS.replace(' ','')
        LHS=LHS.replace('-','+')
        LHS=LHS.replace('**2','tothepowerof2')
        LHS=LHS.replace('**3','tothepowerof3')
        exprList=LHS.split('+')
        for expr in exprList:
            VarList=expr.split('*')
            counter=0
            factors=[]
            for j in range(len(X)):
                anz=0
                if(str(X[j]) in VarList):
                    anz=1
                    factors.append(X[j])
                if((str(X[j])+'tothepowerof2') in VarList):
                    anz=2 
                    factors.append(X[j])
                    factors.append(X[j])
                if((str(X[j])+'tothepowerof3') in VarList):
                    anz=3
                    factors.append(X[j])
                    factors.append(X[j])
                    factors.append(X[j])
                counter=counter+anz
            if(counter==2):
                string=''            
                for l in range(len(factors)):
                    if(l==len(factors)-1):
                        string=string+str(factors[l])
                    else:
                        string=string+str(factors[l])+'*'
                if(not(string in BLList)):
                    BLList.append(string)
            if(counter>2):
                string=''            
                for l in range(len(factors)):
                    if(l==len(factors)-1):
                        string=string+str(factors[l])
                    else:
                        string=string+str(factors[l])+'*'
                if(not(string in MLList)):
                    MLList.append(string)
        
    COPlusLIPlusBL=[]
    for i in range(len(SMF)):
        COPlusLIPlusBL.append(SMF[i])
        for j in range(len(MLList)):
            ToSubs=SMF_exp[i].coeff(MLList[j])
            COPlusLIPlusBL[i]=expand(COPlusLIPlusBL[i]-ToSubs*parse_expr(MLList[j]))
            
    COPlusLI=[]
    for i in range(len(COPlusLIPlusBL)):
        COPlusLI.append(COPlusLIPlusBL[i])
        for j in range(len(BLList)):
            ToSubs=expand((COPlusLIPlusBL)[i]).coeff(BLList[j])
            COPlusLI[i]=expand(COPlusLI[i]-ToSubs*parse_expr(BLList[j]))
    
##### C*X contains linear terms
    C=zeros(len(COPlusLI),len(X))  
    for i in range(len(COPlusLI)):
    	for j in range(len(X)):
    		C[i*len(X)+j]=expand((COPlusLI)[i]).coeff(X[j])
        
##### ML contains multilinearities
    ML=expand(Matrix(SMF)-Matrix(COPlusLIPlusBL))
##### BL contains bilinearities
    BL=expand(Matrix(COPlusLIPlusBL)-Matrix(COPlusLI))    
#### CM is coefficient matrix of linearities
    CM=C        
#####CMBL gives coefficient matrix of bilinearities
    CMBL=[]
    if(BLList!=[]):
        for i in range(len(BLList)):
            CVBL=[]
            for k in range(len(BL)):
                CVBL.append(BL[k].coeff(BLList[i]))
            CMBL.append(CVBL)            
    else:
        CVBL=[]
        for k in range(len(BL)):
            CVBL.append(0)
        CMBL.append(CVBL)
    
    CMBL=Matrix(CMBL).T 
    
#####CMML gives coefficient matrix of multilinearities
#####Summarize multilinearities and bilinearities 
    if(MLList!=[]):
        CMML=[]
        for i in range(len(MLList)):
            CVML=[]
            for k in range(len(ML)):
                CVML.append(expand(ML[k]).coeff(MLList[i]))
            CMML.append(CVML)    
        CMML=Matrix(CMML).T  
        BLList=BLList+MLList
        CMBL=Matrix(concatenate((CMBL,CMML),axis=1))
      
    for i in range(len(BLList)):
        BLList[i]=parse_expr(BLList[i])
       
    if(BLList!=[]):    
        CMbig=Matrix(concatenate((CM,CMBL),axis=1))
    else:
        CMbig=Matrix(CM)      

#### Split purely additive fluxes
    # fluxpars must stay column-aligned and the trafo bodies need
    # flux = fluxpar * prefactor, which a sum like (k_basal + k_ind*C3) has no
    # factor for. Give each summand its own column; SM*F is unchanged. Summands
    # that could carry a minus sign are left alone -- that would put a negative
    # flux in a positive-stoichiometry column.
    col=0
    nsplit=0
    while col<len(F):
        if(GetFluxParameter(F[col], X) is not None):
            col=col+1
            continue
        summands=sympy.Add.make_args(F[col])
        if(len(summands)<2 or
           any(t.could_extract_minus_sign() for t in summands) or
           any(GetFluxParameter(t, X) is None for t in summands)):
            col=col+1
            continue
        for t in summands:
            SM=SM.col_insert(SM.cols, SM.col(col))
            F=F.row_insert(len(F), Matrix(1,1,[t]))
        SM.col_del(col)
        F.row_del(col)
        nsplit=nsplit+1
    if(nsplit>0):
        print('Split '+str(nsplit)+' additive flux(es) into one column per summand.',flush=True)

#### Save ODE equations for testing solutions at the end    
    print('Rank of SM is '+str(SM.rank()) + '!',flush=True)
    SMorig=SM.copy()
    ODE=SMorig*F
#### Get Flux Parameters
    # None for an irreducible sum the split could not take apart; the entry is
    # kept for column alignment and genPriorityTable refuses to pivot on it.
    fluxpars=[GetFluxParameter(flux, X) for flux in F]

    if priority:
        known={str(x) for x in Xo} | {str(fp) for fp in fluxpars if fp is not None}
        unknown=[pr for pr in priority if pr not in known]
        if unknown:
            print('Warning: priority entries match no state or rate parameter '
                  'and are ignored: '+str(unknown),flush=True)
        print('Priority order: '+str(priority),flush=True)

#### Find conserved quantities
    # Computed before sparsification so we can protect CQ-involved state rows
    # during sparsification (see below). FindLCL consumes CMbig, built from
    # the un-sparsified SM*F decomposition above -- independent of sparsify.
    if(givenCQs==[]):
        print('\nFinding conserved quantities ...',flush=True)
        LCLs, rowsToDel=FindLCL(CMbig.transpose(), X)
    else:
        print('\nI took the given conserved quantities!',flush=True)
        LCLs=givenCQs
    LCLs_original=list(LCLs)
    if(LCLs!=[]):
        print(LCLs,flush=True)
    else:
        print('System has no conserved quantities!',flush=True)

##### Sparsification disabled (see Severin Bang's Julia reimplementation)
    # v1.2 used to call Sparsify() here to combine rows of the stoichiometric
    # matrix (rank-preserving) to shorten each ODE. That was removed when we
    # adopted the priority-table-based selection from the Julia port: the
    # priority table's `type` classification (see genPriorityTable) relies on
    # the untouched in/out-flux sign structure of each row, and sparsify can
    # destroy that invariant for states not protected by the CQ-exemption. The
    # sparsifyLevel parameter is kept in the signature for R-side
    # backward-compatibility and is silently ignored.
    if sparsifyLevel != 0:
        print('Note: sparsifyLevel='+str(sparsifyLevel)+' is ignored -- sparsification '
              'was removed in favour of the priority-table cycle-breaking heuristic.',
              flush=True)
#### Define graph structure
    print('\nDefine graph structure ...\n',flush=True)
    
    SSgraph=DetermineGraphStructure(SM, F, X, neglect)    
    #printgraph(SSgraph)
    #print(fluxpars)
#### Priority-table cycle breaking (ported from Severin Bang's Julia
#### reimplementation, helperFunctions.jl::genPriorityTable).
####
#### Each iteration builds a global priority table over all remaining states
#### and both flux sides, sorts by (dontUseThisSide, OccInCycles desc, type,
#### fluxLength, NoCycleOccur desc, OccInRhs), and resolves the top-ranked
#### candidate. This replaces the old FindCycle+GetBestPair pair, which
#### picked the best pair from the FIRST found cycle only -- the global
#### ranking avoids locally-good-but-globally-poor choices and is also what
#### makes the removal of the sparsify step tolerable (the table's `type`
#### classification implicitly penalises sides whose resolution would create
#### new cycles, so the algorithm no longer depends on an a-priori compacted
#### stoichiometry).
    # Track flux parameters that must not be reparameterized by cycle
    # breaking (because an earlier direct-positive-solve already baked them
    # into a state's final expression -- reparameterizing them now would
    # create self-referential eqOut entries, see _try_positive_direct_solve).
    locked_fluxpars=set()
    fluxpar_str_to_sym={str(fp): fp for fp in fluxpars if fp is not None}
    # Lazy-solve bookkeeping: recorded solutions are NOT substituted into F
    # or earlier equations. solved_refs[name] holds the state names a
    # recorded solution references, feeding the alias-aware graph.
    solved_refs={}
    Xo_names={str(s) for s in Xo}
    priority_rows, cycles_list = genPriorityTable(SM, F, fluxpars, X, LCLs, SSgraph, neglect, locked_fluxpars,
                                                  priority)
#### Remove cycles step by step
    gesnew=0
    # (name, numerator flux, reference flux) per r_*: a trafo splits the
    # opposite side's total in the proportions 1 : r_1 : r_2 : ...
    newvars=[]
    eqOut=[]

    def _apply_state_solve(row_idx, y, sol, label):
        # Lazy: record the solution, lock its rate constants, remember which
        # states it references, and drop the row. F and earlier eqOut entries
        # stay untouched -- the one textual resolution at output time
        # substitutes everything, so every affine/sign test in between runs
        # on small expressions. Locks are transitive already: a referenced
        # solved state locked its own rate constants when it was recorded.
        eqOut.append(str(y)+' = '+str(sol))
        print('   Solved '+label+': '+str(y),flush=True)
        sol_sym_names={str(s) for s in sol.free_symbols}
        for nm, sym in fluxpar_str_to_sym.items():
            if nm in sol_sym_names:
                locked_fluxpars.add(sym)
        solved_refs[str(y)]=sol_sym_names & Xo_names
        if positive_syms is not None:
            positive_syms.add(y)
            _added_positive.add(y)
        X.row_del(row_idx)
        SM.row_del(row_idx)

    def _in_active_cq(nm):
        parsed=parse_expr(nm)
        for lcl in LCLs:
            lhs=parse_expr(lcl.split(' = ')[0])
            if lhs.subs(parsed, 1)!=lhs:
                return True
        return False

    def _usable_side_exists(i, locked_names):
        # A pivot side of row i is usable if it has fluxes and none of its
        # rate constants is missing, neglected, or locked. Mirrors the
        # dontUseThisSide rule of genPriorityTable.
        row=SM.row(i)
        for want_pos in (True, False):
            ids=[c for c in range(SM.cols)
                 if (row[c]>0 if want_pos else row[c]<0)]
            if ids and all(fluxpars[c] is not None
                           and str(fluxpars[c]) not in neglect
                           and str(fluxpars[c]) not in locked_names
                           for c in ids):
                return True
        return False

    def _strands_cycle_state(new_lock_names, solving_name, cyc_states):
        # Would adding new_lock_names to the locked set leave a cycle state
        # without any usable pivot side (and no rescuing conservation law)?
        # Pure integer/graph work -- this is the up-front replacement for
        # discovering the deadlock via rollback after the fact.
        locked_now={str(fp) for fp in locked_fluxpars}
        prospective=locked_now | set(new_lock_names)
        name_to_row={str(X[i]): i for i in range(len(X))}
        for nm in cyc_states:
            if nm==solving_name or nm not in name_to_row:
                continue
            if _in_active_cq(nm):
                continue
            i=name_to_row[nm]
            if _usable_side_exists(i, locked_now) and \
               not _usable_side_exists(i, prospective):
                return nm
        return None

    def _drain_positive_solves():
        # Solve states with a certified-positive closed form until none is
        # left. Linear is tried first (cheap, no sqrt), quadratic only when
        # solveQuadratic is on. F never changes here (lazy), so rejections
        # are cached for the whole sweep; each accepted solve is vetoed if
        # its locks would strand a cycle state.
        nonlocal SSgraph
        if not allow_positive_solves:
            return False
        rejected_lin=set()
        rejected_quad=set()
        cyc_states=_states_in_cycles(SSgraph)
        progress=False
        while True:
            cand=_try_positive_direct_solve(SM, F, X, positive_syms, neglect,
                                            rejected_lin, solved_refs)
            label='directly'
            if cand is None and solveQuadratic:
                quad=_try_positive_quadratic_solve(SM, F, X, positive_syms,
                                                   branches, neglect,
                                                   rejected_quad, solved_refs)
                if quad is not None:
                    row_idx, y, sol, branch=quad
                    cand=(row_idx, y, sol)
                    label='quadratically' + (' (2 branches)'
                                             if branch is not None else '')
            if cand is None:
                break
            row_idx, y, sol=cand
            new_locks={str(s) for s in sol.free_symbols} & set(fluxpar_str_to_sym)
            victim=_strands_cycle_state(new_locks, str(y), cyc_states)
            if victim is not None:
                print('   Skipping direct solve of '+str(y)+': locking its '
                      'rate constants would leave '+victim+' without a '
                      'usable pivot side.',flush=True)
                rejected_lin.add(str(y))
                rejected_quad.add(str(y))
                continue
            _apply_state_solve(row_idx, y, sol, label)
            progress=True
            SSgraph=DetermineGraphStructure(SM, F, X, neglect, solved_refs)
            cyc_states=_states_in_cycles(SSgraph)
        return progress

    # The lock guard in _drain_positive_solves refuses solves that would
    # strand a cycle state, but it reasons on the current graph -- keep the
    # snapshot/rollback as a safety net for deadlocks it cannot foresee.
    allow_positive_solves=True
    retried_without_positive_solves=False
    _snapshot=(SM.copy(), F.copy(), X.copy(), list(fluxpars), list(LCLs),
               counter, gesnew, list(newvars))

    if cycles_list:
        print('\nDirect positive-solve pass ...',flush=True)
        if _drain_positive_solves():
            priority_rows, cycles_list = genPriorityTable(SM, F, fluxpars, X, LCLs, SSgraph, neglect, locked_fluxpars,
                                                          priority)
        else:
            print('   (nothing resolvable positively up-front)',flush=True)

    while cycles_list:
        # Re-check positive-solves before every cycle-break: a break's trafo
        # rewrites fluxes, which may have turned a previously-unsolvable
        # state into a linear-in-itself or quadratic-with-positive-root one.
        if _drain_positive_solves():
            priority_rows, cycles_list = genPriorityTable(SM, F, fluxpars, X, LCLs, SSgraph, neglect, locked_fluxpars,
                                                          priority)
            if not cycles_list:
                break
        row = priority_rows[0]
        minType  = row['type']
        state2Rem = row['species']
        index    = row['spIndex']
        isOutflux = row['isOutflux']
        anz      = row['fluxLength']
        # Map priority-row's side selection onto the (sign, signChanged)
        # flags the existing type-2/3 bodies already consume -- they branch
        # on `(sign=="minus") XOR signChanged`, so setting signChanged=False
        # and sign from isOutflux gives the body the right branch directly.
        sign = "minus" if isOutflux else ("plus" if isOutflux is not None else None)
        signChanged = False
        fp2Rem = row['fluxPars'][0] if (minType == 1 and row['fluxPars']) else None

        print('Removing cycle '+str(counter),flush=True)
        printPriorityTable(priority_rows)
        # Unresolvable candidate: no usable side for the top priority row.
        if row['dontUseThisSide'] and minType != 0:
            # Withholding only the blocking solves does find resolutions with
            # no ratio parameter at all, but they cost two orders of magnitude in
            # expression size (TGFb, Smad2: 869 -> 30164 chars) and stop
            # simplifying. Hence all-or-nothing.
            if locked_fluxpars and not retried_without_positive_solves:
                allow_positive_solves=False
                retried_without_positive_solves=True
                print('   Every side is blocked by a flux parameter locked by the '
                      'positive-solve',flush=True)
                print('   pass -- rolling back and retrying without it.',flush=True)
                SM, F, X = (_snapshot[0].copy(), _snapshot[1].copy(),
                            _snapshot[2].copy())
                fluxpars=list(_snapshot[3])
                LCLs=list(_snapshot[4])
                counter, gesnew = _snapshot[5], _snapshot[6]
                newvars=list(_snapshot[7])
                eqOut=[]
                locked_fluxpars=set()
                solved_refs.clear()
                for s in _added_positive:
                    positive_syms.discard(s)
                _added_positive.clear()
                SSgraph=DetermineGraphStructure(SM, F, X, neglect)
                priority_rows, cycles_list = genPriorityTable(SM, F, fluxpars, X, LCLs, SSgraph, neglect,
                                                              locked_fluxpars, priority)
                continue
            # Take any simple cycle the state participates in for the
            # diagnostic. Falls back to a single-node "cycle" if the table
            # picked a non-cycle state (shouldn't happen when cycles_list is
            # non-empty but keeps the error path robust).
            diag_cycle = next((c for c in cycles_list if state2Rem in c),
                              [state2Rem])
            unique_nodes=list(dict.fromkeys(diag_cycle))
            print("",flush=True)
            print("    ======================================================",flush=True)
            print("    CYCLE CANNOT BE REMOVED",flush=True)
            print("    ======================================================",flush=True)
            print(f"    Cycle: {unique_nodes}",flush=True)
            print("",flush=True)
            # --- Per-node diagnosis ---
            print("    State of each node:",flush=True)
            for node in unique_nodes:
                dim, sign = GetDimension(node, X, SM, True)
                negfps = GetNegFluxParameters(SM, fluxpars, X, node)
                posfps = GetPosFluxParameters(SM, fluxpars, X, node)
                neg_in_neglect = [str(fp) for fp in negfps if str(fp) in neglect]
                pos_in_neglect = [str(fp) for fp in posfps if str(fp) in neglect]
                print(f"      {node}:",flush=True)
                if(len(posfps)==0 and len(negfps)==0):
                    print(f"        0 influxes, 0 outfluxes",flush=True)
                    print(f"        All fluxes were absorbed by previous cycle removals.",flush=True)
                    print(f"        The solver has no free parameter left to solve for {node}.",flush=True)
                elif(len(posfps)==0):
                    print(f"        0 influxes, {len(negfps)} outflux(es): {[str(fp) for fp in negfps]}",flush=True)
                    print(f"        No production term -> cannot balance in/outfluxes.",flush=True)
                elif(len(negfps)==0):
                    print(f"        {len(posfps)} influx(es): {[str(fp) for fp in posfps]}, 0 outfluxes",flush=True)
                    print(f"        No degradation/consumption term -> cannot balance in/outfluxes.",flush=True)
                else:
                    blocked = neg_in_neglect + pos_in_neglect
                    if(blocked):
                        print(f"        influxes: {[str(fp) for fp in posfps]}, outfluxes: {[str(fp) for fp in negfps]}",flush=True)
                        print(f"        Blocked by neglect: {blocked}",flush=True)
                    else:
                        print(f"        influxes: {[str(fp) for fp in posfps]}, outfluxes: {[str(fp) for fp in negfps]}",flush=True)
            # --- Conserved quantities context ---
            print("",flush=True)
            print("    Conserved quantities (original):",flush=True)
            if(LCLs_original):
                for lcl in LCLs_original:
                    used = "(available)" if lcl in LCLs else "(already used)"
                    # check if any node from cycle appears in this CQ
                    involves_cycle = False
                    ls=parse_expr(lcl.split(' = ')[0])
                    for node in unique_nodes:
                        if(ls.subs(parse_expr(node),1)!=ls):
                            involves_cycle = True
                    tag = " <-- involves " + ", ".join(unique_nodes) if involves_cycle else ""
                    print(f"      {lcl}  {used}{tag}",flush=True)
            else:
                print("      (none detected)",flush=True)
            # --- Actionable suggestions ---
            print("",flush=True)
            print("    What you can do:",flush=True)
            has_absorbed = any(
                len(GetPosFluxParameters(SM, fluxpars, X, n))==0 and
                len(GetNegFluxParameters(SM, fluxpars, X, n))==0
                for n in unique_nodes
            )
            has_blocked = any(
                any(str(fp) in neglect for fp in
                    GetPosFluxParameters(SM, fluxpars, X, n) +
                    GetNegFluxParameters(SM, fluxpars, X, n))
                for n in unique_nodes
            )
            neglected_nodes = [n for n in unique_nodes if n in neglect]
            nodes_str = ", ".join(unique_nodes)
            _item=[0]
            def nextItem():
                _item[0]+=1
                return str(_item[0])
            if(has_absorbed):
                print(f"      {nextItem()}. Supply a conserved quantity (givenCQs) that includes {nodes_str}.",flush=True)
                print(f"         This lets the solver express {nodes_str} in terms of other states",flush=True)
                print(f"         and a total-amount parameter, without needing flux parameters.",flush=True)
                print(f"         Example: givenCQs = c(\"{unique_nodes[0]} + ... = total{unique_nodes[0]}\")",flush=True)
            if(has_blocked):
                print(f"      {nextItem()}. Remove blocked parameters from 'neglect'.",flush=True)
            if(neglected_nodes):
                print(f"      {nextItem()}. Remove {neglected_nodes} from 'neglect'.",flush=True)
                print(f"         A neglected state must spend its ODE on a rate-parameter pivot,",flush=True)
                print(f"         and no usable flux side is left for it.",flush=True)
            print(f"      {nextItem()}. Review the model reactions involving {nodes_str}:",flush=True)
            print(f"         Does {nodes_str} participate in enough independent reactions?",flush=True)
            print(f"         A state needs both production and consumption to be solvable.",flush=True)
            # Extra note when the model has CQ-involved bilinear coupling that
            # the sign-preserving sparsify protection cannot resolve -- this is
            # the failure mode of models like TGFb (pSmad2/pSmad3/Smad4 linked
            # by k_form*pSmad2*pSmad3*Smad4 -> C3).
            cq_related = any(
                any(str(parse_expr(lcl.split(' = ')[0]).subs(parse_expr(n),1))
                    != lcl.split(' = ')[0] for lcl in LCLs_original)
                for n in unique_nodes
            )
            if cq_related:
                print(f"      {nextItem()}. This model has bilinear coupling between conservation-law states.",flush=True)
                print(f"         No flux-parameter pivot keeps the steady state a manifestly",flush=True)
                print(f"         positive rational function here. Try solveQuadratic = TRUE, which",flush=True)
                print(f"         admits closed-form positive roots at the price of sqrt terms, or",flush=True)
                print(f"         steer the pivot choice with 'priority' / 'neglect'.",flush=True)
            print("    ======================================================",flush=True)
            return(0)
        if(minType==0):
            for LCL in LCLs:
                ls=parse_expr(LCL.split(' = ')[0])
                if(ls.subs(parse_expr(state2Rem),1)!=ls):
                    LCL2Rem=LCL
            LCLs.remove(LCL2Rem)
            eqOut.append(state2Rem+' = '+state2Rem)
            print('   '+str(state2Rem)+' --> '+'Done by CQ',flush=True)
        if(minType==1):
            eq=sympy.S.Zero
            for k in range(SM.cols):
                if SM[index,k]!=0:
                    eq=eq+SM[index,k]*F[k]
            sol=solve(eq, fp2Rem, simplify=False)[0]
            eqOut.append(str(fp2Rem)+' = '+str(sol))
            print('   '+str(state2Rem)+' --> '+str(fp2Rem),flush=True)
        if(minType==2):
            negs, sumnegs, negfps=GetOutfluxes(state2Rem, X, SM, F, fluxpars)
            poss, sumposs, posfps=GetInfluxes(state2Rem, X, SM, F, fluxpars)
            if(anz==1):
                print("Error in Type Determination. Please report this bug!",flush=True)
                return(0)
            else:
                nenner=1
                for j in range(anz):
                    if(j>0):
                        nenner=nenner+parse_expr('r_'+state2Rem+'_'+str(j))
                trafoList=[]
                if((sign=="minus" and not signChanged) or (sign=="plus" and signChanged)):
                    # A pivot rate constant can sit on the opposite side as
                    # well (recycling loops): resolve the sum against the pivot
                    # relations before distributing it, see _resolvePivotSum.
                    sumposs=_resolvePivotSum(sumposs, negs, negfps, nenner, state2Rem)
                    for j in range(len(negs)):
                        flux=negs[j]
                        fp=negfps[j]
                        prefactor=flux/fp
                        if(j==0):
                            trafoList.append(str(fp)+' = ('+str(sumposs)+')*1/('+str(nenner)+')*1/('+str(prefactor)+')')
                        else:
                            gesnew=gesnew+1
                            newvars.append(('r_'+state2Rem+'_'+str(j), str(fp), str(negfps[0])))
                            trafoList.append(str(fp)+' = ('+str(sumposs)+')*'+'r_'+state2Rem+'_'+str(j)+'/('+str(nenner)+')*1/('+str(prefactor)+')')                        
                    print('   '+str(state2Rem)+' --> '+str(negfps),flush=True)
                    
                else:
                    sumnegs=_resolvePivotSum(sumnegs, poss, posfps, nenner, state2Rem)
                    for j in range(len(poss)):
                        flux=poss[j]
                        fp=posfps[j]
                        prefactor=flux/fp
                        if(j==0):
                            trafoList.append(str(fp)+' = ('+str(sumnegs)+')*1/('+str(nenner)+')*1/('+str(prefactor)+')')
                        else:
                            gesnew=gesnew+1
                            newvars.append(('r_'+state2Rem+'_'+str(j), str(fp), str(posfps[0])))
                            trafoList.append(str(fp)+' = ('+str(sumnegs)+')*'+'r_'+state2Rem+'_'+str(j)+'/('+str(nenner)+')*1/('+str(prefactor)+')')
                    print('   '+str(state2Rem)+' --> '+str(posfps),flush=True)
                for eq in trafoList:
                    eqOut.append(eq)
        if(minType==3):
            negs, sumnegs, negfps=GetOutfluxes(state2Rem, X, SM, F, fluxpars)
            poss, sumposs, posfps=GetInfluxes(state2Rem, X, SM, F, fluxpars)
            if(anz==1):
                if((sign=="minus" and not signChanged) or (sign=="plus" and signChanged)):
                    fp2Rem=negfps[0]
                    flux=negs[0]
                else:
                    fp2Rem=posfps[0]
                    flux=poss[0]
                eq=sympy.S.Zero
                for k in range(SM.cols):
                    if SM[index,k]!=0:
                        eq=eq+SM[index,k]*F[k]
                sol=solve(eq, fp2Rem, simplify=False)[0]
                eqOut.append(str(fp2Rem)+' = '+str(sol))
                FsearchFlux = matrix_multiply_elementwise(abs(SM[index,:]),F.T)
                colindex=list(FsearchFlux).index(flux)
                for row2repl in range(len(SM.col(0))):
                    if(SM[row2repl,colindex]!=0 and row2repl!=index):
                        SM=SM.row_insert(row2repl,SM.row(row2repl)-(SM[row2repl,colindex]/SM[index,colindex])*SM.row(index))
                        SM.row_del(row2repl+1)
                #print('HELP',flush=True)
            else:
                nenner=1
                for j in range(anz):
                    if(j>0):
                        nenner=nenner+parse_expr('r_'+state2Rem+'_'+str(j))
                trafoList=[]
                # Each substituted flux parameter on the chosen side becomes
                # sum_opp * r_j / nenner (r_0=1, r_j for j>0). After the
                # substitution, the RAW flux F[j] (i.e. flux/fp_j rescaled)
                # equals sum_opp * r_j / (nenner * |SM[pivot,j]|) -- the
                # stoichiometric multiplier is absorbed by the prefactor.
                # When we split the pivot-side column into len(opp_fps) new
                # columns (one per opposite-side flux), each new raw flux
                # must therefore carry r_j / |SM[pivot,j]|, not just r_j.
                # The pre-priority-table code dropped BOTH the r_j weight
                # (masked when only +-1 stoichiometries were involved) and
                # the stoichiometry division (masked further because
                # sparsify tended to produce +-1 rows). Both factors are
                # needed for the pivot and non-pivot rows to stay
                # consistent.
                if((sign=="minus" and not signChanged) or (sign=="plus" and signChanged)):
                    # A pivot rate constant can sit on the opposite side as
                    # well (recycling loops): resolve the sum against the pivot
                    # relations before distributing it, see _resolvePivotSum.
                    sumposs=_resolvePivotSum(sumposs, negs, negfps, nenner, state2Rem)
                    for j in range(len(negs)):
                        flux=negs[j]
                        fp=negfps[j]
                        prefactor=flux/fp
                        if(j==0):
                            trafoList.append(str(fp)+' = ('+str(sumposs)+')*1/('+str(nenner)+')*1/('+str(prefactor)+')')
                            r_weight = parse_expr('1')
                        else:
                            gesnew=gesnew+1
                            newvars.append(('r_'+state2Rem+'_'+str(j), str(fp), str(negfps[0])))
                            trafoList.append(str(fp)+' = ('+str(sumposs)+')*'+'r_'+state2Rem+'_'+str(j)+'/('+str(nenner)+')*1/('+str(prefactor)+')')
                            r_weight = parse_expr('r_'+state2Rem+'_'+str(j))

                        FsearchFlux = matrix_multiply_elementwise(abs(SM[index,:]),F.T)
                        colindex=list(FsearchFlux).index(flux)
                        s_abs = abs(SM[index, colindex])
                        for k in range(len(posfps)):
                            SM=SM.col_insert(len(SM.row(0)),SM.col(colindex))
                            F=F.row_insert(len(F),Matrix(1,1,[cancel(poss[k]*r_weight/(nenner*s_abs))]))
                            fluxpars.append(posfps[k])
                        SM.col_del(colindex)
                        F.row_del(colindex)
                        fluxpars.__delitem__(colindex)
                    print('   '+str(state2Rem)+' --> '+str(negfps),flush=True)

                else:
                    sumnegs=_resolvePivotSum(sumnegs, poss, posfps, nenner, state2Rem)
                    for j in range(len(poss)):
                        flux=poss[j]
                        fp=posfps[j]
                        prefactor=flux/fp
                        if(j==0):
                            trafoList.append(str(fp)+' = ('+str(sumnegs)+')*1/('+str(nenner)+')*1/('+str(prefactor)+')')
                            r_weight = parse_expr('1')
                        else:
                            gesnew=gesnew+1
                            newvars.append(('r_'+state2Rem+'_'+str(j), str(fp), str(posfps[0])))
                            trafoList.append(str(fp)+' = ('+str(sumnegs)+')*'+'r_'+state2Rem+'_'+str(j)+'/('+str(nenner)+')*1/('+str(prefactor)+')')
                            r_weight = parse_expr('r_'+state2Rem+'_'+str(j))
                        FsearchFlux = matrix_multiply_elementwise(abs(SM[index,:]),F.T)
                        colindex=list(FsearchFlux).index(flux)
                        s_abs = abs(SM[index, colindex])
                        for k in range(len(negfps)):
                            SM=SM.col_insert(len(SM.row(0)),SM.col(colindex))
                            F=F.row_insert(len(F),Matrix(1,1,[cancel(negs[k]*r_weight/(nenner*s_abs))]))
                            fluxpars.append(negfps[k])
                        SM.col_del(colindex)
                        F.row_del(colindex)
                        fluxpars.__delitem__(colindex)
                    print('   '+str(state2Rem)+' --> '+str(posfps),flush=True)
                for eq in trafoList:
                    eqOut.append(eq)
        X.row_del(index)
        SM.row_del(index)
        SSgraph=DetermineGraphStructure(SM, F, X, neglect, solved_refs)
        priority_rows, cycles_list = genPriorityTable(SM, F, fluxpars, X, LCLs, SSgraph, neglect, locked_fluxpars,
                                                      priority)
        counter=counter+1
    print('There is no cycle in the system!\n',flush=True)
    
#### Solve remaining equations
    eqOut.reverse()
    print('Solving remaining equations ...\n',flush=True)
    while len(X)>0:
        # Rebuild with aliases each round: a recorded solution's references
        # count as dependencies, so leaf order stays a valid resolution
        # order even though nothing is substituted into F.
        SSgraph=DetermineGraphStructure(SM, F, X, neglect, solved_refs)
        node=FindNodeToSolve(SSgraph, SM, F, X)
        if node is None:
            print('   WARNING: no solvable leaf left -- aborting solve phase.',flush=True)
            break
        xi=parse_expr(node)
        index=list(X).index(xi)
        # After cycle-breaking, this row is guaranteed linear in xi.
        # Compute just the row (avoids a full SM*F).
        expr=sympy.S.Zero
        for k in range(SM.cols):
            if SM[index,k]!=0:
                expr=expr+SM[index,k]*F[k]
        # Direct linear solve: expr = In + Out*xi  =>  xi = -In/Out.
        # No sympy.solve(), no simplify -- cancel keeps the result compact.
        In=expr.subs(xi, 0)
        Out=expr.diff(xi)
        if cancel(Out)==0:
            # The ODE has collapsed into a 0=0 identity after upstream
            # substitutions -- typically because the state is constrained by
            # a conservation law whose other members have already been
            # resolved. Keep `xi` as a free parameter.
            sol=xi
        else:
            sol=cancel(-In/Out)
        # Unlike _try_positive_direct_solve this phase has no sign gate: an
        # upstream substitution can turn an xi-proportional sink into a constant
        # one, leaving `sol` a difference. The ODE is linear in xi, so -In/Out is
        # the unique root -- the only way out is to spend the row on one of the
        # state's own rate constants, which keeps xi free and the pivot a sum.
        sol_cls=_rational_sign_class(sol, positive_syms) if _enforcePositive else '+'
        pivoted=False
        tried=[]
        if node in neglect or sol_cls not in ('+','0'):
            for fp in _exclusiveFluxPivots(SM, F, fluxpars, index, neglect):
                psol=solve(expr, fp, simplify=False)
                if not psol:
                    continue
                psol=cancel(psol[0])
                if _enforcePositive and \
                   _rational_sign_class(psol, positive_syms) not in ('+','0'):
                    tried.append(str(fp))
                    continue
                eqOut.insert(0,str(fp)+' = '+str(psol))
                print(f'   Solved {node} --> {fp}'
                      +('  (neglect)' if node in neglect else '  (positivity)'),
                      flush=True)
                pivoted=True
                break
        if not pivoted:
            if node in neglect:
                print(f"   WARNING: {node} is in 'neglect' but no usable rate "
                      f"pivot is left -- solving for {node} itself.",flush=True)
            if sol_cls not in ('+','0'):
                _reportSignIndefinite(
                    node, sol, sol_cls, positive_syms,
                    {e.split(' = ',1)[0] for e in eqOut},
                    sorted({str(bfp) for bfp in
                            _exclusiveFluxPivots(SM, F, fluxpars, index, ())
                            if str(bfp) in neglect}),
                    tried, solveQuadratic)
                return(0)
            eqOut.insert(0,node+' = '+str(sol))
            if sol is not xi:
                solved_refs[node]={str(s) for s in sol.free_symbols} & Xo_names
                if positive_syms is not None:
                    positive_syms.add(xi)
            print(f'   Solved {node}',flush=True)
        X.row_del(index)
        SM.row_del(index)

#### Order equations: definitions before uses
    # Solutions are recorded lazily (in terms of each other), so the
    # resolution below and the mod-p test both need a topological order.
    eqOut=_topo_sort_eqs(eqOut)

#### Optional final simplification (applied once per expression, after the loop)
    # simplify:
    #   True   -> sympy.simplify once per expression (default, heuristic)
    #   'full' -> aggressive pipeline: cancel + posify + simplify(ratio=oo)
    #            + separate factor of numerator and denominator. Slower,
    #            but often gives noticeably more compact formulas.
    #   False  -> skip entirely (expressions stay in cancel() normal form)
    _full=(isinstance(simplify, str) and simplify.lower()=='full')
    if(simplify):
        print(('Full-simplifying' if _full else 'Simplifying')+
              ' final expressions ...',flush=True)
        _simplify_aborted=False
        for i in range(len(eqOut)):
            if _check_walltime():
                # Budget exhausted -- leave the remaining expressions in
                # their cancelled-but-not-simplified form. They're still
                # correct (cancel() kept them in rational normal form),
                # just bulkier than they could be.
                print(f'   Walltime exceeded after simplifying {i}/{len(eqOut)} '
                      f'expressions -- leaving the rest un-simplified.',flush=True)
                _simplify_aborted=True
                break
            ls, rs = eqOut[i].split(' = ', 1)
            eqOut[i]=ls+' = '+_finalSimplify(rs, _full)
            print(f'   Simplified {ls}',flush=True)
    
#### Test Solution
    # testSteady: 'fast' = probabilistic GF(p), 'exact' = symbolic, else skip.
    if testSteady in ('fast', 'modp', 'exact', 'T'):
        fast = testSteady in ('fast', 'modp')
        print('Testing Steady State'+(' (probabilistic mod p)' if fast else '')+'...\n',flush=True)
        NonSteady=False
        bad=None
        if fast:
            # Environment extension over the topologically sorted equations;
            # no symbolic substitution into the ODEs at all.
            bad=_steady_test_fast(ODE, eqOut, zeroStates)
            if bad:
                for i in bad:
                    print('   Equation '+str(ODE[i]),flush=True)
                    print('   is nonzero at a random point in GF(p)',flush=True)
                NonSteady=True
            elif bad is None:
                print('   (falling back to the exact symbolic test)',flush=True)
        if not fast or bad is None:
            # Exact test: substitute dependents before their dependencies
            # (reverse topological order), so every reference resolves.
            subs_pairs=[(parse_expr(ls), parse_expr(rs))
                        for ls, rs in (eq.split(' = ', 1) for eq in reversed(eqOut))
                        if ls!=rs]
            for i in range(len(ODE)):
                expr=parse_expr(str(ODE[i]))
                for zs in zeroStates:
                    expr=expr.subs(zs, 0)
                for lsym, rexpr in subs_pairs:
                    expr=expr.subs(lsym, rexpr)
                expr=_simplify(expr)
                if(expr!=0):
                    print('   Equation '+str(ODE[i]),flush=True)
                    print('   results:'+str(expr),flush=True)
                    NonSteady=True
        if(NonSteady):
            print('Solution is wrong!\n',flush=True)
        elif fast:
            print('Solution is correct (almost surely, mod p)!\n',flush=True)
        else:
            print('Solution is correct!\n',flush=True)
    else:
        print('Skipping the Testing of Steady State...\n',flush=True)
    
#### Print Equations
    # Echo every forcing as = 0: injections are held at 0 and dropped as
    # exogenous, so add any not already zeroed structurally.
    zero_names={str(s) for s in zeroStates}
    for inj in injections:
        if inj not in zero_names:
            zeroStates.append(parse_expr(inj))
            zero_names.add(inj)
    print('I obtained the following equations:\n',flush=True)
    if(outputFormat=='M'):
        eqOutReturn=[]
        for state in zeroStates:
            print('\tinit_'+str(state)+'  "0"'+'\n',flush=True)
            eqOutReturn.append('init_'+str(state)+'  "0"')
        for i in range(len(eqOut)):
            ls, rs = eqOut[i].split('=')
            ls=parse_expr(ls)
            rs=parse_expr(rs)
            for j in range(i,len(eqOut)):
                ls2, rs2 = eqOut[j].split('=')
                rs2=parse_expr(rs2)
                rs2=rs2.subs(ls,rs)
                eqOut[j]=str(ls2)+'='+str(rs2)
            for state in Xo:
                ls=ls.subs(state, parse_expr('init_'+str(state)))
                rs=rs.subs(state, parse_expr('init_'+str(state)))
            eqOut[i]=str(ls)+'  "'+str(rs)+'"'
                            
        for i in range(len(eqOut)):
            eqOut[i]=eqOut[i].replace('**','^')
                    
        for eq in eqOut:
            print('\t'+eq+'\n',flush=True)
            eqOutReturn.append(eq)            
        
    else:
        # Resolve the recurrence before printing: a dMod trafo substitutes all
        # entries at once ('M' above does the same inline). Textually, not via
        # sympy -- subs() re-flattens the nested rational and throws away the
        # factored form. \b works as a boundary because '_' is a word character,
        # so 'R1' matches neither 'R1mRNA' nor 'R1_R2_TGFb'.
        substituted=set()
        for i in range(len(eqOut)):
            ls, rs = eqOut[i].split(' = ', 1)
            if(ls==rs):
                continue
            pattern=re.compile(r'\b'+re.escape(ls)+r'\b')
            for j in range(i+1, len(eqOut)):
                ls2, rs2 = eqOut[j].split(' = ', 1)
                resolved=pattern.sub('('+rs+')', rs2)
                if(resolved!=rs2):
                    eqOut[j]=ls2+' = '+resolved
                    substituted.add(j)
        # Substituting nests two separately simplified expressions, so a factor
        # can straddle the new boundary. Re-run what changed.
        if(simplify and substituted):
            print(('Full-simplifying' if _full else 'Simplifying')+
                  ' resolved expressions ...',flush=True)
            for j in sorted(substituted):
                if _check_walltime():
                    print('   Walltime exceeded -- leaving the remaining resolved '
                          'expressions un-simplified.',flush=True)
                    break
                ls2, rs2 = eqOut[j].split(' = ', 1)
                eqOut[j]=ls2+' = '+_finalSimplify(rs2, _full)
                print(f'   Simplified {ls2}',flush=True)
        eqOutReturn=[]
        for state in zeroStates:
            print('\t'+str(state)+' = 0'+'\n',flush=True)
            eqOutReturn.append(str(state)+'=0')
        for eq in eqOut:
            ls, rs = eq.split(' = ')
            print('\t'+ls+' = "'+rs+'",'+'\n',flush=True)
            eqOutReturn.append(ls+'='+rs)
    print('Number of Species:  '+str(nrspecies),flush=True)
    print('Number of Equations:  '+str(len(eqOut)+len(zeroStates)),flush=True)
    print('Number of new introduced variables:  '+str(gesnew),flush=True)
    for nm, fp, ref in newvars:
        print('\t'+nm+' = flux('+fp+') / flux('+ref+') > 0',flush=True)
    return(eqOutReturn)

// The steady states of Pimpl() in C++: the roots of its system by cppDE's pseudo-transient
// continuation, from the start the parameters give, with their derivatives by the implicit
// function theorem. A root depends on the parameters, the start and the controls only.

#ifndef DMOD2_PIMPL_KERNEL_H
#define DMOD2_PIMPL_KERNEL_H

#include <Rinternals.h>

#include <string>
#include <vector>

namespace dmod {
namespace pimpl {

struct State;

// x exp(C' lambda) on C x = total for the n positive states xa; C is k x n, column-major.
void onManifold(double* xa, int n, const double* C, int k, const double* total);

// The steady state of `st` at its parameters `pv` from the start `guess` (the dependent states),
// with derivatives of order `order`.
struct Request {
  State* st = nullptr;
  std::vector<double> pv, guess;
  int order = 0;
};

// The root x [ndep], D [ndep, npar] and D2 [ndep, npar, npar], where the order asked for them
// and the implicit function theorem gives them.
struct Root {
  std::vector<double> x, D, D2;
  bool hasD = false, hasD2 = false;
};

// The roots of `req`, repeated values from memory, the others in `threads` threads, the same for
// every thread count. False with `error` where one is not found; `warnings` collects what R
// reports.
bool roots(const std::vector<Request>& req, std::vector<Root>& out, int threads, std::string& error,
           std::vector<std::string>& warnings);

// The state behind an external pointer of pimpl_state(), nullptr where it is gone.
State* stateOf(SEXP xp);

}  // namespace pimpl
}  // namespace dmod

#endif

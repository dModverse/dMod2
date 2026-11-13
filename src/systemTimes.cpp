#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif
#include <Rcpp.h>

// Idle and total processor time of the machine in seconds since boot, summed
// over the cores; NA outside Windows.
// [[Rcpp::export(name = ".systemTimes")]]
Rcpp::NumericVector systemTimes() {
#ifdef _WIN32
  FILETIME idle, kernel, user;
  if (GetSystemTimes(&idle, &kernel, &user)) {
    auto sec = [](const FILETIME& t) {
      return (static_cast<double>(t.dwHighDateTime) * 4294967296.0 + t.dwLowDateTime) * 1e-7;
    };
    // kernel time includes the idle time
    return Rcpp::NumericVector::create(Rcpp::_["idle"] = sec(idle),
                                       Rcpp::_["total"] = sec(kernel) + sec(user));
  }
#endif
  return Rcpp::NumericVector::create(Rcpp::_["idle"] = NA_REAL, Rcpp::_["total"] = NA_REAL);
}

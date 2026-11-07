// The heap memory the process freed, returned to the system where the C library allows it.

#include <Rcpp.h>
#ifdef __GLIBC__
#include <malloc.h>
#endif

// [[Rcpp::export]]
void release_heap() {
#ifdef __GLIBC__
  malloc_trim(0);
#endif
}

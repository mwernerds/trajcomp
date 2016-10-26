#include<Rcpp.h>
using namespace Rcpp;
// Following Simon: http://permalink.gmane.org/gmane.comp.lang.r.devel/27627
// This should check for interrupts without longjumping away from destructors


static void chkIntFn(void *dummy) {
  R_CheckUserInterrupt();
}

// this will call the above in a top-level context so it won't longjmp-out of your context
bool checkInterrupt() {
  return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}



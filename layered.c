#include "grid/cartesian1D.h"
#include "saint-venant.h"

scalar lamda[];
scalar etta[];

int main() {
  X0 = 0.;
  L0 = 21.;
  G = 9.81;
  N = 256;

 nu = 0.01;
  lambda_b = lamda;

  nl = 2;  run();
  nl = 5;  run();
  nl = 15; run();
}

h[right]   = dirichlet(0.6);
etta[right] = dirichlet(0.6);

scalar hc[];

event init (i = 0) {
  foreach() {
    zb[] = max(0., 0.2*(1. - 1./sq(5.75/2.)*sq(x - 10.)));
    hc[] = h[]  = 0.6 - zb[];
  }

 event ("friction");

 for (vector u in ul) {
    u.n[left] = dirichlet(h[left] ? 1./h[left] : 0.);
    u.n[right] = neumann(0.);
  }
}

event friction (i++) {
  foreach() {
    double U = 0.;
    int l = 0;
    for (vector u in ul)
      U += u.x[]*layer[l++];
    double S = 25., k = G/(sq(S)*pow(h[],1./3.))*fabs(U);
    lamda[] = k > 0. ? nu/k : 0.;
  }
  boundary ({lamda});
}

event logfile (t += 0.1; i <= 100000) {
  double dh = change (h, hc);
  if (i > 0 && dh < 1e-4)
    return 1;
}

event output (t = end) {
  char name[80];
  sprintf (name, "end-%d", nl);
  FILE * fp = nl == 15 ? stderr : fopen (name, "w");
  foreach() {
    fprintf (fp, "%g %g %g\n", x, etta[], zb[]);
    if (nl == 15) {
      double z = zb[];
      int l = 0;
      printf ("%g %g %g\n", x, z, u.x[]);
      for (vector u in ul) {
	z += layer[l++]*h[];
	printf ("%g %g %g\n", x, z, u.x[]);
      }
      printf ("\n");
    }
  }
}

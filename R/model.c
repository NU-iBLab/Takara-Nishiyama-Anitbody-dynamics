#include <R.h>
/* a trick to keep up with the parameters */
static double parms[7];

#define M33 parms[0]
#define mu1 parms[1]
#define mu2 parms[2]
#define m parms[3]
#define K parms[4]
#define delay3 parms[5]
#define vac3 parms[6]

/* initializers */
void initparms( void (* odeparms)(int *, double *) ) {
    int N = 7;
    odeparms(&N,parms);
}

/* names for states and derivatives */
#define M var[0]

#define dMdt vardot[0]

void derivs( int *neq, double *t, double *var, double *vardot, double *varout, int *ip ) {

    if( ip[0]<1 ) {
        error("nout should be at least 1");
    }
    double dg=0.693;
    double G0=50;

    if( (*t) <= vac3+delay3 ){
        dMdt = -mu1*M;
    }
    if( (*t) > vac3+delay3 ){
        double g = G0*exp(-dg*(*t-vac3));
        dMdt = M33*(pow(g,m))/(pow(g,m)+pow(K,m))- mu2*M;
    }
    
}

typedef struct cond_t {
        int n;               /* No. timesteps in conductance profile */
        double dt;           /* Integration time step (msecs) */
        double *g;           /* Conductance profile for prescribed dt (mu a) */
        double tol;          /* Determines off-threshold for synapse */
        double tau;          /* Rise time for synapse (msecs) */
        double gmax;         /* Peak conductance per synapse (mu) */
} cond;

typedef struct synapse_t
{
        cond *cp;            /* Pointer to profile of synaptic conductance */
        int ind;             /* Indicator of status of synapse */
        int max;             /* Maximum value of status counter */
        double gold;         /* Previous synaptic conductance */
        double vsyn;         /* Reversal potential of synapse (mV) */
} synapse;

void cgs(int, sparse_mat *, double *, double *, double *, double),
     FreeConductanceProfile( cond *),
     UpdateSynapticConductance( synapse *syn);
cond *ConductanceProfile( double, double, double, double);
synapse *SynapticProfile( int, int, double, double, double, cond * );
double c_in(synapse *, long int);

/* Declaration of HH coefficient functions */
double alfa_h( double );
double alfa_m( double );
double alfa_n( double );
double beta_h( double );
double beta_m( double );
double beta_n( double );

double temperature_correction;


 /********************************************
        Computes a conductance profile
  ********************************************/
cond *ConductanceProfile( double dt,   /* Integration time step */
                          double tau,  /* Synaptic time constant */
                          double tol,  /* Determines off-threshold for synapse */
                          double gmax  /* Maximum conductance (mS) */ )
{
    int k;
    cond *con;
    double tmp, told, tnew;

    con = (cond *) malloc( sizeof(cond) );
    con->dt = dt;
    con->tau = tau;
    con->tol = tol;
    con->gmax = gmax;

/* Iterate to find duration of pulse */
    tmp = 1.0-log(tol);
    tnew = tmp;
    do {
        told = tnew;
        tnew = tmp+log(told);
    } while ( fabs(tnew-told)>5.e-11 );
    con->n = ((int) ceil(tau*tnew/dt));
    con->g = (double *) malloc( (con->n)*sizeof(double) );
    con->g[0] = 0.0;
    for ( k=1 ; k<(con->n) ; k++ ) {
        tmp = dt*((double) k)/tau;
        con->g[k] = gmax*tmp*exp(1.0-tmp);
    }
    return con;
}


 /********************************************
     Function to free a conductance profile
  ********************************************/
void FreeConductanceProfile( cond *con)
{
    free(con->g);
    free(con);
    return;
}


 /*******************************************************
    Function to initialise and update synaptic activity
  *******************************************************/
void UpdateSynapticConductance( synapse *syn)
{
    extern unsigned long int ix, iy, iz;
    static double RealRate;
    static int start=1;
    int j, k;
    double deadtime, dt, interval;

/*  Step 1. - Initialise synapse */
    if ( start ) {
        if ( syn->cp ) {
            dt = (syn->cp->dt);
            syn->max = syn->cp->n;
            syn->gold = 0.0;
            deadtime = ((double) syn->cp->n)*dt;
            RealRate = (1000.0/RATE)-deadtime;
            if ( RealRate <= 0.0 ) {
                printf("\nImpossible firing rate requested\n");
                fflush(stdout);
                return NULL;
            }
            interval = -RATE*log(ran( &ix, &iy, &iz))/dt;
            if ( fmod(interval,1.0) <= 0.5 ) {
               syn->ind = -((int) floor(interval));
            } else {
               syn->ind = -((int) ceil(interval));
            }
        } else {
            printf("\nSynapse has no conductance profile!\n");
            return NULL;
        }
        start =0;
    } else {
        (syn->ind)++;
        if ( syn->ind == syn->max ) {
            dt = (syn->cp->dt);
            syn->gold = 0.0;
            interval = -RealRate*log(ran( &ix, &iy, &iz))/dt;
            if ( fmod(interval,1.0) <= 0.5 ) {
               syn->ind = -((int) floor(interval));
            } else {
               syn->ind = -((int) ceil(interval));
            }
        }
    }
    return;
}



 /********************************************************************
                    ALPHA for ACTIVATION OF SODIUM
  *******************************************************************/
double alfa_m( double volt )
{
    extern double temperature_correction;
    double tmp;

    tmp = 2.5-0.1*volt;
    if ( fabs(tmp)<0.001 ) {
        tmp = 1.0/(((tmp/24.0+1.0/6.0)*tmp+0.5)*tmp+1.0);
    } else {
        tmp = tmp/(exp(tmp)-1.0);
    }
    return tmp*temperature_correction;
}


 /********************************************************************
                    BETA for ACTIVATION OF SODIUM
  ********************************************************************/
double beta_m( double volt )
{
    extern double temperature_correction;
    double tmp;

    tmp = volt/18.0;
    return 4.0*temperature_correction*exp(-tmp);
}


 /********************************************************************
                    ALPHA for INACTIVATION OF SODIUM
  ********************************************************************/
double alfa_h( double volt )
{
    extern double temperature_correction;
    double tmp;

    tmp = 0.05*volt;
    return 0.07*temperature_correction*exp(-tmp);
}


 /********************************************************************
                    BETA for INACTIVATION OF SODIUM
  ********************************************************************/
double beta_h( double volt )
{
    extern double temperature_correction;
    double tmp;

    tmp = 3.0-0.1*volt;
    return temperature_correction/(exp(tmp)+1.0);
}


 /********************************************************************
                    ALPHA for ACTIVATION OF POTASSIUM
  ********************************************************************/
double alfa_n( double volt )
{
    extern double temperature_correction;
    double tmp;

    tmp = 1.0-0.1*volt;
    if ( fabs(tmp)<0.001 ) {
        tmp = 0.1/(((tmp/24.0+1.0/6.0)*tmp+0.5)*tmp+1.0);
    } else {
        tmp = 0.1*tmp/(exp(tmp)-1.0);
    }
    return tmp*temperature_correction;
}


 /********************************************************************
                    BETA for ACTIVATION OF POTASSIUM
  ********************************************************************/
double beta_n( double volt )
{
    extern double temperature_correction;
    double tmp;

    tmp = 0.0125*volt;
    return 0.125*temperature_correction*exp(-tmp);
}



void cgs(int getmem, sparse_mat *a, double *b, double *x0, double *x, double tol)
{
    long int i, k, n, finish;
    static int start=1;
    double rho_old, rho_new, alpha, beta, sigma, err;
    static double *p, *q, *u, *v, *r, *rb, *y;

/* Step 1 - Check memory status */
    n = a->n;
    if ( start ) {
        r = (double *) malloc( n*sizeof(double) );
        rb = (double *) malloc( n*sizeof(double) );
        p = (double *) malloc( n*sizeof(double) );
        q = (double *) malloc( n*sizeof(double) );
        u = (double *) malloc( n*sizeof(double) );
        v = (double *) malloc( n*sizeof(double) );
        y = (double *) malloc( n*sizeof(double) );
        start = 0;
    }

/* Step 2 - Initialise residual, p[ ] and q[ ] */
    mat_vec_mult( a, x0, r);
    for ( rho_old=0.0,i=0 ; i<n ; i++ ) {
        r[i] = b[i]-r[i];
        rho_old += r[i]*r[i];
        rb[i] = r[i];
        p[i] = 0.0;
        q[i] = 0.0;
    }
    if ( rho_old<tol*((double) n) ) {
        for ( i=0 ; i<n ; i++ ) x[i] = x0[i];
        return;
    }
    rho_old = 1.0;
    finish = 0;

/* The main loop */
    while ( !finish ) {

/* Compute scale parameter for solution update */
        for ( rho_new=0.0,i=0 ; i<n ; i++ ) rho_new += r[i]*rb[i];
        beta = rho_new/rho_old;

/* Update u[ ] and p[ ] */
        for ( i=0 ; i<n ; i++ ) {
            u[i] = r[i]+beta*q[i];
            p[i] = u[i]+beta*(q[i]+beta*p[i]);
        }

/* Update v[ ] and compute sigma */
        mat_vec_mult( a, p, v);
        for ( sigma=0.0,i=0 ; i<n ; i++ ) sigma += rb[i]*v[i];

/* Compute alpha and update q[ ], v[ ] and x[ ] */
        alpha = rho_new/sigma;
        for ( i=0 ; i<n ; i++ ) {
            q[i] = u[i]-alpha*v[i];
            v[i] = alpha*(u[i]+q[i]);
            x[i] += v[i];
        }

 /* Update r[ ] and estimate error */
        mat_vec_mult( a, v, y);
        for ( err=0.0,i=0 ; i<n ; i++ ) {
            r[i] -= y[i];
            err += r[i]*r[i];
        }
        rho_old = rho_new;
        if ( err<tol*((double) n) ) finish = 1;
    }

/* Check memory status */
    if ( getmem<=0 ) start = 1;
    if ( start ) {
        free(r);
        free(rb);
        free(p);
        free(q);
        free(u);
        free(v);
        free(y);
    }
    return;
}

/**********************************************************
                Returns Gaussian deviate
  **********************************************************/
double normal( double mean, double sigma)
{
    static int start=1;
    long int n;
    static double g1, g2;
    double v1, v2, w, ran(long int *);

    if ( start ) {
        n = 1;
        do {
            v1 = 2.0*ran(&n)-1.0;
            v2 = 2.0*ran(&n)-1.0;
            w = v1*v1+v2*v2;
        } while ( w==0.0 || w>=1.0 );
        w = log(w)/w;
        w = sqrt(-w-w);
        g1 = v1*w;
        g2 = v2*w;
        start = !start;
        return (mean+sigma*g1);
    } else {
        start = !start;
        return (mean+sigma*g2);
    }
}

 /**********************************************************
         Order integers in ascending order by heapsort
  **********************************************************/
void heapsort( long int n, long int *x)
{
    int finish;
    long int i, ir, j, k, tmp;

    if ( n<2 ) return;
    k = n/2;
    ir = n-1;
    finish = 0;
    while ( !finish ) {
        if ( k>0 ) {
            tmp = x[--k];
        } else {
            tmp = x[ir];
            x[ir] = x[0];
            if ( --ir==0 ) {
                x[0] = tmp;
                finish = 1;
            }
        }
        i = k;
        j = 2*k+1;
        while ( j<=ir ) {
            if ( j<ir && x[j]<x[j+1] ) j++;
            if ( tmp<x[j] ) {
                x[i] = x[j];
                i = j;
                j = 2*j+1;
            } else {
                j = ir+1;
            }
        }
        x[i] = tmp;
    }
    return;
}

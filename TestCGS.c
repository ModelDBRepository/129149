#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

 /***************************************************************
    Function to analyse the numerical error of Generalised
    Compartmental Model. A single input is simulated and
    the exact solution determined using the Equivalent Cable
  ***************************************************************/

typedef struct SparseMatrix_t
{
        double *a;
        int *col;
        int *StartRow;
        int n;
        struct SparseMatrix_t *l;
        struct SparseMatrix_t *u;

} SparseMatrix;


double  ran(unsigned int *, unsigned int *, unsigned int *);
void solve( int, double *, double *, double *);

void    Matrix_Vector_Multiply( SparseMatrix *, double *, double *),
        Matrix_Malloc( SparseMatrix *, int, int),
        Matrix_Free( SparseMatrix *),
        cgs1( int *, SparseMatrix *, double *, double *, double),
        cgs( int *, SparseMatrix *, double *, double *, double);

/* Global definitions */
#define            CGS    1.0e-26   /* Tolerance used in CGS algorithm */
#define          NODES    50
#define          NSEED    2         /* Seed for random number generator */

/* Global Variables */
unsigned int ix, iy, iz;

int main( void )
{
    extern unsigned int ix, iy, iz;
    SparseMatrix a;
    extern double pi, dt;
    int k, j, jj, getmem;
    double *vec, *phi, *soln, vold, vnew, *vlam;
    void srand( unsigned int);

/*  Initialise random number generator */
    srand( ((unsigned int) NSEED) );
    ix = rand( );
    iy = rand( );
    iz = rand( );
    soln = (double *) malloc( 4*sizeof(double) );
    vlam = (double *) malloc( 4*sizeof(double) );
    phi = (double *) malloc( 4*sizeof(double) );
    vec = (double *) malloc( 4*sizeof(double) );

/*  Phase 1. - Allocate sparse matrices for synaptic activity */
    Matrix_Malloc( &a, 4, 10 );
    for ( j=0 ; j<3 ; j++ ) {
        jj = 2*j;
        a.StartRow[j] = jj;
        a.a[jj] = 1.0;
        a.col[jj] = j;
        a.a[jj+1] = -1.0;
        a.col[jj+1] = j+1;
    }
    a.StartRow[3] = 6;
    vnew = 0.0;
    vlam[0] = vnew;
    for ( j=0 ; j<3 ; j++ ) {
        jj = 6+j;
        vold = vnew;
        vnew += 0.25;
        vlam[j+1] = vnew;
        a.a[jj] = vnew-vold;
        a.col[jj] = j;
    }
    jj = 9;
    a.a[jj] = 1.0-vnew;
    a.col[jj] = 3;
    a.StartRow[4] = 10;

    for ( j=0 ; j<4 ; j++ ) soln[j] = ran( &ix, &iy, &iz);
    Matrix_Vector_Multiply( &a, soln, phi);
    for ( j=0 ; j<4 ; j++ ) vec[j] = 0.0;
    getmem = 1;
    cgs1( &getmem, &a, phi, vec, CGS);
    for ( k=0 ; k<4 ; k++ ) {
        printf("\nTrue %12.6lf \t Numerical %12.6lf", soln[k], vec[k]);
    }


    solve( 3, vlam, phi, vec);
    printf("\n\n");
    for ( k=0 ; k<4 ; k++ ) {
        printf("\nTrue %12.6lf \t Numerical %12.6lf", soln[k], vec[k]);
    }

    return 0;
}


 /**********************************************************
     Multiplies sparse matrix a[ ][ ] with vector v[ ]
  **********************************************************/
void Matrix_Vector_Multiply( SparseMatrix *a, double *v , double *b)
{
    int i, j, k, n;

    n = a->n;
    for ( j=0 ; j<n ; j++) {
        k = a->StartRow[j+1];
        for( b[j]=0.0,i=(a->StartRow[j]) ; i<k ; i++ ) {
            b[j] += (a->a[i])*v[a->col[i]];
        }
    }
    return;
}


 /***********************************************
        Allocate memory to a sparse matrix
  ***********************************************/
void Matrix_Malloc( SparseMatrix *a, int n, int w)
{
    a->a = (double *) malloc( w*sizeof(double) );
    a->col = (int *) malloc( w*sizeof(int) );
    a->StartRow = (int *) malloc( (n+1)*sizeof(int) );
    a->n = n;
    a->l = malloc(sizeof(SparseMatrix));
    a->u = malloc(sizeof(SparseMatrix));
    a->l->a = (double *) malloc( (2*n-1)*sizeof(double) );
    a->l->col = (int *) malloc( (2*n-1)*sizeof(int) );
    a->l->StartRow = (int *) malloc( (n+1)*sizeof(int) );
    a->l->n = n;
    a->u->a = (double *) malloc( (2*n-1)*sizeof(double) );
    a->u->col = (int *) malloc( (2*n-1)*sizeof(int) );
    a->u->StartRow = (int *) malloc( (n+1)*sizeof(int) );
    a->u->n = n;
    return;
}


 /**********************************************
     De-allocates memory of a sparse matrix
  **********************************************/
void Matrix_Free( SparseMatrix *a)
{
    free(a->a);
    free(a->col);
    free(a->StartRow);
    free(a);
}



 /************************************************************
         Function returns primitive uniform random number.
  ************************************************************/
double ran(unsigned int *ix, unsigned int *iy, unsigned int *iz)
{
    double tmp;

/*  1st item of modular arithmetic  */
    *ix = (171*(*ix))%30269;
/*  2nd item of modular arithmetic  */
    *iy = (172*(*iy))%30307;
/*  3rd item of modular arithmetic  */
    *iz = (170*(*iz))%30323;
/*  Generate random number in (0,1) */
    tmp = ((double) (*ix))/30269.0+((double) (*iy))/30307.0
          +((double) (*iz))/30323.0;
    return fmod(tmp,1.0);
}


void cgs1(int *getmem, SparseMatrix *a, double *b, double *x, double tol)
{
    long int i, k, n, repeat;
    static int start=1;
    double rho_old, rho_new, alpha, beta, sigma, err;
    static double *p, *q, *r, *u, *v, *rb, *y;

/* Step 1 - Check memory status */
    printf("\nEntering 1 ");
    n = a->n;
    if ( start ) {
        *getmem = 1;
        start = 0;
    }
    if ( *getmem ) {
        if ( p ) free(p);
        p = (double *) malloc( n*sizeof(double) );
        if ( q ) free(q);
        q = (double *) malloc( n*sizeof(double) );
        if ( r ) free(r);
        r = (double *) malloc( n*sizeof(double) );
        if ( u ) free(u);
        u = (double *) malloc( n*sizeof(double) );
        if ( v ) free(v);
        v = (double *) malloc( n*sizeof(double) );
        if ( rb ) free(rb);
        rb = (double *) malloc( n*sizeof(double) );
        if ( y ) free(y);
        y = (double *) malloc( n*sizeof(double) );
        *getmem = 0;
    }

/* Step 2 - Initialise residual, p[ ] and q[ ] */
    Matrix_Vector_Multiply( a, x, r);
    for ( rho_old=0.0,i=0 ; i<n ; i++ ) {
        r[i] = b[i]-r[i];
        rho_old += r[i]*r[i];
        rb[i] = r[i];
        p[i] = 0.0;
        q[i] = 0.0;
    }
    if ( rho_old<tol*((double) n) ) return;


/* The main loop */
    rho_old = 1.0;
    do {

/* Compute scale parameter for solution update */
        for ( rho_new=0.0,i=0 ; i<n ; i++ ) rho_new += r[i]*rb[i];
        beta = rho_new/rho_old;

/* Update u[ ] and p[ ] */
        for ( i=0 ; i<n ; i++ ) {
            u[i] = r[i]+beta*q[i];
            p[i] = u[i]+beta*(q[i]+beta*p[i]);
        }

/* Update v[ ] and compute sigma */
        Matrix_Vector_Multiply( a, p, v);
        for ( sigma=0.0,i=0 ; i<n ; i++ ) sigma += rb[i]*v[i];

/* Compute alpha and update q[ ], v[ ] and x[ ] */
        if ( sigma==0.0 ) {
            printf(" Trouble ");
            for (i=0 ; i<n ; i++ ) {
                printf("\n%20.16lf",v[i]);
                getchar( );
            }
        }
        alpha = rho_new/sigma;
        for ( i=0 ; i<n ; i++ ) {
            q[i] = u[i]-alpha*v[i];
            v[i] = alpha*(u[i]+q[i]);
            x[i] += v[i];
        }

 /* Update r[ ] and estimate error */
        Matrix_Vector_Multiply( a, v, y);
        for ( err=0.0,i=0 ; i<n ; i++ ) {
            r[i] -= y[i];
            err += r[i]*r[i];
        }
        rho_old = rho_new;
        repeat = ( err > tol*((double) n) );
    } while ( repeat );
    printf(" Leaving\n");
    return;
}


void cgs(int *getmem, SparseMatrix *a, double *b, double *x, double tol)
{
    int i, k, n, finish;
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
    Matrix_Vector_Multiply( a, x, r);
    for ( i=0 ; i<n ; i++ ) {
        r[i] = b[i]-r[i];
        rb[i] = r[i];
        p[i] = 0.0;
        q[i] = 0.0;
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
        Matrix_Vector_Multiply( a, p, v);
        for ( sigma=0.0,i=0 ; i<n ; i++ ) sigma += rb[i]*v[i];

/* Compute alpha and update q[ ], v[ ] and x[ ] */
        alpha = rho_new/sigma;
        for ( i=0 ; i<n ; i++ ) {
            q[i] = u[i]-alpha*v[i];
            v[i] = alpha*(u[i]+q[i]);
            x[i] += v[i];
        }

 /* Update r[ ] and estimate error */
        Matrix_Vector_Multiply( a, v, y);
        for ( err=0.0,i=0 ; i<n ; i++ ) {
            r[i] -= y[i];
            err += r[i]*r[i];
        }
        rho_old = rho_new;
        if ( err<tol ) finish = 1;
    }

/* Check memory status */
    if ( *getmem<=0 ) start = 1;
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


void solve( int nsyn, double *vlam, double *b, double *soln)
{
    double sum;
    int k;

/*  Step 1. - Forward substitution phase */
    for ( sum=0.0,k=0 ; k<nsyn ; k++ ) {
        soln[k] = b[k];
        sum += vlam[k+1]*soln[k];
    }
    soln[nsyn] = b[nsyn]-sum;

/*  Step 2. - Backward substitution phase */
    for ( k=nsyn-1 ; k>=0 ; k-- ) soln[k] += soln[k+1];
    return;
}

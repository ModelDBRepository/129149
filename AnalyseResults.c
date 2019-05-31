#include <stdio.h>
#include <stdlib.h>
#include <math.h>

 /***************************************************************
   Function to generate exact solutions and locations of stimuli
  ***************************************************************/

/* Global definitions */
#define     INPUT1    "ExactSolution.dat"   /* Exact Solutions */
#define     INPUT2    "NRes500.dat"         /* NEURON Solutions */
#define     INPUT3    "NewRes500.dat"       /* New Solutions */
#define     INPUT4    "Old500.dat"          /* Conventional Solutions */
#define     OUTPUT    "Ratio.dat"           /* Output results */
#define      NODES     500

/* Parameters for solutions */
#define       NSIM     2000  /* Number of simulations */
#define       NREC     10    /* Number of data */

int main( void )
{
    int j, k;
    double **os, **es, **ns, **cs;
    double *n_mu_err, *o_mu_err, *c_mu_err, *n_sd_err, *o_sd_err, *c_sd_err;
    double tmp, dim, nsum, osum, csum;
    FILE *fp;

/*  STEP 1. - Get data */
    os = (double **) malloc( NSIM*sizeof(double *) );
    es = (double **) malloc( NSIM*sizeof(double *) );
    ns = (double **) malloc( NSIM*sizeof(double *) );
    cs = (double **) malloc( NSIM*sizeof(double *) );
    for ( k=0 ; k<NSIM ; k++ ) {
        os[k] = (double *) malloc( NREC*sizeof(double) );
        es[k] = (double *) malloc( NREC*sizeof(double) );
        ns[k] = (double *) malloc( NREC*sizeof(double) );
        cs[k] = (double *) malloc( NREC*sizeof(double) );
    }
    fp = fopen(INPUT1,"r");
    for ( k=0 ; k<NSIM ; k++ ) {
        for ( j=0 ; j<NREC ; j++ ) fscanf(fp,"%lf", &es[k][j]);
    }
    fclose(fp);

    fp = fopen(INPUT2,"r");
    for ( k=0 ; k<NSIM ; k++ ) {
        for ( j=0 ; j<NREC ; j++ ) fscanf(fp,"%lf", &ns[k][j]);
    }
    fclose(fp);

    fp = fopen(INPUT3,"r");
    for ( k=0 ; k<NSIM ; k++ ) {
        for ( j=0 ; j<NREC ; j++ ) fscanf(fp,"%lf", &os[k][j]);
    }
    fclose(fp);

    fp = fopen(INPUT4,"r");
    for ( k=0 ; k<NSIM ; k++ ) {
        for ( j=0 ; j<NREC ; j++ ) fscanf(fp,"%lf", &cs[k][j]);
    }
    fclose(fp);

/*  STEP 2. - Compute statistics of results */
    n_mu_err = (double *) malloc( NREC*sizeof(double) );
    n_sd_err = (double *) malloc( NREC*sizeof(double) );
    o_mu_err = (double *) malloc( NREC*sizeof(double) );
    o_sd_err = (double *) malloc( NREC*sizeof(double) );
    c_mu_err = (double *) malloc( NREC*sizeof(double) );
    c_sd_err = (double *) malloc( NREC*sizeof(double) );

    for ( j=0 ; j<NREC ; j++ ) {
        n_mu_err[j] = 0.0;
        n_sd_err[j] = 0.0;
        o_mu_err[j] = 0.0;
        o_sd_err[j] = 0.0;
        c_mu_err[j] = 0.0;
        c_sd_err[j] = 0.0;
    }
    for ( k=0 ; k<NSIM ; k++ ) {
        for ( j=0 ; j<NREC ; j++ ) {
            tmp = fabs(es[k][j]-ns[k][j]);
            tmp /= fabs(es[k][j]);
            n_mu_err[j] += tmp;
            n_sd_err[j] += tmp*tmp;

            tmp = fabs(es[k][j]-os[k][j]);
            tmp /= fabs(es[k][j]);
            o_mu_err[j] += tmp;
            o_sd_err[j] += tmp*tmp;

            tmp = fabs(es[k][j]-cs[k][j]);
            tmp /= fabs(cs[k][j]);
            c_mu_err[j] += tmp;
            c_sd_err[j] += tmp*tmp;
        }
    }
    dim = ((double) NSIM);
    for ( j=0 ; j<NREC ; j++ ) {
        n_mu_err[j] /= dim;
        n_sd_err[j] /= dim;
        n_sd_err[j] = sqrt(n_sd_err[j]-n_mu_err[j]*n_mu_err[j]);
        o_mu_err[j] /= dim;
        o_sd_err[j] /= dim;
        o_sd_err[j] = sqrt(o_sd_err[j]-o_mu_err[j]*o_mu_err[j]);
        c_mu_err[j] /= dim;
        c_sd_err[j] /= dim;
        c_sd_err[j] = sqrt(c_sd_err[j]-c_mu_err[j]*c_mu_err[j]);
    }

/*  STEP 3 - Output results */
    fp = fopen(OUTPUT,"w");
    fprintf(fp,"%s\n\n","NEURON RESULTS");
    for ( j=0 ; j<NREC ; j++ ) fprintf(fp,"%12.8lf", n_mu_err[j] );
    fprintf(fp,"\n");
    for ( j=0 ; j<NREC ; j++ ) fprintf(fp,"%12.8lf", log10(n_mu_err[j]) );
    fprintf(fp,"\n\n");
    for ( j=0 ; j<NREC ; j++ ) fprintf(fp,"%12.8lf", n_sd_err[j] );
    fprintf(fp,"\n");
    for ( j=0 ; j<NREC ; j++ ) fprintf(fp,"%12.8lf", log10(n_sd_err[j]) );
    fprintf(fp,"\n\n\n");
    fprintf(fp,"%s\n\n","OUR RESULTS");
    for ( j=0 ; j<NREC ; j++ ) fprintf(fp,"%12.8lf", o_mu_err[j] );
    fprintf(fp,"\n");
    for ( j=0 ; j<NREC ; j++ ) fprintf(fp,"%12.8lf", log10(o_mu_err[j]) );
    fprintf(fp,"\n\n");
    for ( j=0 ; j<NREC ; j++ ) fprintf(fp,"%12.8lf", o_sd_err[j] );
    fprintf(fp,"\n");
    for ( j=0 ; j<NREC ; j++ ) fprintf(fp,"%12.8lf", log10(o_sd_err[j]) );
    fprintf(fp,"\n\n\n");
    fprintf(fp,"%s\n\n","CONVENTIONAL RESULTS");
    for ( j=0 ; j<NREC ; j++ ) fprintf(fp,"%12.8lf", c_mu_err[j] );
    fprintf(fp,"\n");
    for ( j=0 ; j<NREC ; j++ ) fprintf(fp,"%12.8lf", log10(c_mu_err[j]) );
    fprintf(fp,"\n\n");
    for ( j=0 ; j<NREC ; j++ ) fprintf(fp,"%12.8lf", c_sd_err[j] );
    fprintf(fp,"\n");
    for ( j=0 ; j<NREC ; j++ ) fprintf(fp,"%12.8lf", log10(c_sd_err[j]) );
    fclose(fp);

/*  STEP 3 - Output results
    fp = fopen(OUTPUT,"a");
    fprintf(fp,"%s \t %d \n\n","RESULTS", NODES);
    for ( nsum=osum=csum=0.0,j=0 ; j<NREC ; j++ ) {
        nsum += n_mu_err[j];
        osum += o_mu_err[j];
        csum += c_mu_err[j];
    }
    fprintf(fp,"%12.8lf \t %12.8lf\n", nsum/osum, csum/osum);
    for ( nsum=osum=csum=0.0,j=0 ; j<NREC ; j++ ) {
        nsum += n_sd_err[j];
        osum += o_sd_err[j];
        csum += c_sd_err[j];
    }
    fprintf(fp,"%12.8lf \t %12.8lf\n\n", nsum/osum, csum/osum);

    fclose(fp);*/
    return 0;
}

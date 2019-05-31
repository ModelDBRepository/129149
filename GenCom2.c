#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct sparse_mat_t
{
        double *a;
        int *col;
        int *start_row;
        int n;
        struct sparse_mat_t *l;
        struct sparse_mat_t *u;

} sparse_mat;


typedef struct syn_cond_t {
        int n;               /* No. timesteps in the pulse */
        int id;              /* Identifier of synapse type */
        double dt;           /* Integration time step (msecs) */
        double *g;           /* Conductance profile for prescribed dt (mS/cm^2) */
        double tol;          /* Threshold activity for synapse */
        double tau;          /* Rise time for synapse (msecs) */
        double gmax;         /* Maximum conductance per synapse (mS/cm^2) */
} syn_cond;


typedef struct contact_t
{
        int id;                     /* Identifies contact type */

        double xc;                  /* X coordinate of contact */
        double yc;                  /* Y coordinate of contact */
        double zc;                  /* Z coordinate of contact */
        double dc;                  /* Dendritic diameter at contact */

        double xp;                  /* Projected X coordinate */
        double yp;                  /* Projected Y coordinate */
        double zp;                  /* Projected Z coordinate */

        double sd;                  /* Shortest distance to dendrite (micron) */
        double pl;                  /* Measurement of physical length (micron) */
        double amp;                 /* Strength of contact */

        int xn;                     /* Node of contact */
        syn_cond *cp;               /* Profile of synaptic conductance */
        int hist;                   /* Stores location of firing cycle */
        double vsyn;                /* Reversal potential for synapse (mV) */
        double deadtime;            /* Inactivation time after firing (msec) */
        struct contact_t *next;     /* Address of next contact */
} contact;


typedef struct soma_t
{
/*  Physical properties of soma */
        int nobs;                   /* No. of observations in somal specification */

        double *x;                  /* X-coords of defining point */
        double *y;                  /* Y-coords of defining point */
        double *z;                  /* Z-coords of defining point */
        double *d;                  /* Diameter of soma at point (x,y,z) */
        double p_len;               /* Length of soma */

/*  Static biophysical properties of soma */
        double cs;                  /* Somal membrane capacitance (mu F/cm^2) */
        double ga;                  /* Intracellular conductance of soma (mS/cm) */
        double gs;                  /* Membrane conductance of soma (mS/cm^2) */

/*  Contact information */
        contact *conlist;           /* List of contacts */
        int ncon;
} soma;


typedef struct branch_t
{
/*  Connectivity of branch */
        struct branch_t *parent;    /* Pointer to parent branch */
        struct branch_t *child;     /* Pointer to child branch */
        struct branch_t *peer;      /* Pointer to a peer branch */

/*  Physical properties of branch */
        int nobs;                   /* No. of observations in branch specification */
        double *x;                  /* X-coordinate of defining point */
        double *y;                  /* Y-coordinate of defining point */
        double *z;                  /* Z-coordinate of defining point */
        double *d;                  /* Diameter of dendrite */
        double p_len;               /* Length of branch */
        double *pl;                 /* Measurement of physical length (micron) */

/*  Biophysical properties of branch */
        double cm;                  /* Dendritic membrane capacitance (mu F/cm^2) */
        double ga;                  /* Intracellular conductance of dendrite (mS/cm) */
        double gm;                  /* Membrane conductance of dendrite (mS/cm^2) */

/*  Node information for spatial representation */
        int nodes;                  /* Total number nodes spanning branch */
        int junct;                  /* Junction node of the branch */
        int first;                  /* Internal node connected to junction */

/*  Contact information */
        contact *conlist;           /* List of contacts */
        int ncon;
} branch;


typedef struct dendrite_t
{
        branch *root;               /* Pointer to root branch of dendrite */
        double p_len;               /* Total length of dendrite */
} dendrite;


typedef struct neuron_t
{
        int ndend;                  /* Number of dendrites */
        dendrite *dendlist;         /* Pointer to an array of dendrites */
        soma *s;                    /* Soma structure */
} neuron;



/* Function type declarations */
neuron *Load_Sampled_Neuron(char *);
neuron *Load_Test_Neuron(char *);
void Destroy_Sampled_Neuron(neuron *);
void init_branch( branch *, int, double, double, double);
void BuildContactInfo(contact *, branch *, branch **);
void remove_branch( branch **, branch *);
void build_dendrite( branch **, branch *);
void clean_dendrite( branch *);
void destroy_dendrite( branch *);
int count_branches( branch *, branch *);
double branch_length( branch *, branch *);
void enumerate_nodes( branch *, int *);
void mat_vec_mult( sparse_mat *, double *, double *);
void mat_malloc( sparse_mat *, int, int);
void mat_free( sparse_mat *);
void Generate_Dendrite(branch *, int *);
void init_test_branch( branch *, int , double , double , double );
void build_test_dendrite( branch **, branch *);
void LU_Factor(sparse_mat *, int *);
void LU_Solve( sparse_mat *, double *, double *);
void Assign_Current( branch *, double *, double );
void count_branch ( branch *, int * );
void count_dendritic_length ( branch *, double *);
void phys_lengths( branch *, double * );
syn_cond *cond_profile( double , double , double , double );
void free_cond_profile( syn_cond *);
void Initialise_Contact( branch *, syn_cond *);
double ran(unsigned long int *, unsigned long int *, unsigned long int *);
void InsertContacts( branch *);


/*  Declaration of HH coefficient functions */
double alfa_h( double );
double alfa_m( double );
double alfa_n( double );
double beta_h( double );
double beta_m( double );
double beta_n( double );

/* Global definitions */
#define        CELSIUS    18.5    /* Celsius temperature of neuron */
#define             CS    1.0
#define             GA    14.286
#define             CM    1.0
#define             GM    0.091
#define         OUTPUT    "contactinfo.dat"
#define           TEND    1000.0
#define             NT    100
#define             DT    0.1
#define            SIN    0.0e-3
#define             RS    0.002
#define         NNODES    100
#define         SENTRY    50.0
#define          NSEED    10

/* Synaptic Properties */

#define          GMAX          1.e-1  /* Maximum synaptic conductance (mS) */
#define           TAU           0.50  /* Rise-time for injected current (m sec) */
#define          VSYN           55.0  /* Synaptic reversal potential (mV) */
#define         DTIME            0.0  /* Dead time between synaptic activation (m sec) */
#define          RATE           30.0  /* Average number of spikes per second */

/* Global Variables */
sparse_mat lhs, rhs;
unsigned long int ix, iy, iz;
double veq;

int main( int argc, char **argv )
{
    int k, id, start, nodes, nc, i, in, FirstNode, first;
    int counter, nb, nspk;
    extern unsigned long int ix, iy, iz;
    double *v, *x, *y, max, as, gama, xold, xnew, vs, frac,
           arg, sum, tmp, tmp1, tmp2, dt, tnow, tout, hnew,
           hold, nnew, nold, mold, mnew, vnew, vold, jold,
           jnew, neq, meq, heq, jeq, gval, pi,
           cval, length, h, vmid, *LhsStore, *RhsStore;
    double g_na = 120.0, g_k = 36.0, g_l = 0.3, v_na = 55.0,
           v_k = -72.0, v_l = -49.416, fact=25.0;
    void srand( unsigned int);
    neuron *n;
    syn_cond *cond;
    extern sparse_mat lhs, rhs;
    extern double veq;
    branch *bnow;
    FILE *fp;

/*  Initialise random number generator */
    srand( ((unsigned int) NSEED) );
    ix = rand();
    iy = rand();
    iz = rand();

/*  Compute ancillary variables */
    pi = 4.0*atan(1.0);
    as = 4.0*pi*RS*RS;

/*  Compute Equilibrium Potential */
    vold = -62.0;
    mold = alfa_m(vold)/(alfa_m(vold)+beta_m(vold));
    hold = alfa_h(vold)/(alfa_h(vold)+beta_h(vold));
    nold = alfa_n(vold)/(alfa_n(vold)+beta_n(vold));
    jold = g_na*pow(mold,3)*hold*(vold-v_na)+g_k*pow(nold,4)*(vold-v_k)+g_l*(vold-v_l);

    vnew = -58.0;
    mnew = alfa_m(vnew)/(alfa_m(vnew)+beta_m(vnew));
    hnew = alfa_h(vnew)/(alfa_h(vnew)+beta_h(vnew));
    nnew = alfa_n(vnew)/(alfa_n(vnew)+beta_n(vnew));
    jnew = g_na*pow(mnew,3)*hnew*(vnew-v_na)+g_k*pow(nnew,4)*(vnew-v_k)+g_l*(vnew-v_l);

    if ( jold*jnew > 0.0 ) {
        printf(" No zero found \n");
        return(0);
    } else {
        while ( fabs(vold-vnew) > 5.e-7 ) {
            veq = 0.5*(vold+vnew);
            meq = alfa_m(veq)/(alfa_m(veq)+beta_m(veq));
            heq = alfa_h(veq)/(alfa_h(veq)+beta_h(veq));
            neq = alfa_n(veq)/(alfa_n(veq)+beta_n(veq));
            jeq = g_na*pow(meq,3)*heq*(veq-v_na)+g_k*pow(neq,4)*(veq-v_k)+g_l*(veq-v_l);
            if ( jeq*jold > 0.0 ) {
                vold = veq;
            } else {
                vnew = veq;
            }
        }
    }
    vnew = veq;

    mnew = alfa_m(vnew)/(alfa_m(vnew)+beta_m(vnew));
    hnew = alfa_h(vnew)/(alfa_h(vnew)+beta_h(vnew));
    nnew = alfa_n(vnew)/(alfa_n(vnew)+beta_n(vnew));

    v_na -= veq;
    v_k -= veq;
    v_l -= veq;

/*  Load sampled neuron */
    if ( argc != 2 && argc != 3 ) {
        printf("\n Invoke program with load <input>\n");
        return(1);
    } else if ( argc == 2 ) {
        n = Load_Sampled_Neuron( argv[1] );
        if ( !n ) {
            printf("\n Failed to find sampled neuron\n");
            return(1);
        }

        for(i=0; i<n->ndend; i++) InsertContacts( n->dendlist[i].root);

    } else if ( strcmp(argv[1], "test") == 0  || strcmp(argv[1], "Test") == 0 ) {
        n = Load_Test_Neuron( argv[2] );

        for(i=0; i<n->ndend; i++) InsertContacts( n->dendlist[i].root);

        length = 0.0;
        nb = 0;

        for(i=0; i<n->ndend; i++) count_branch( n->dendlist[i].root, &nb );
        for(i=0; i<n->ndend; i++) count_dendritic_length( n->dendlist[i].root, &length );

        h = length/((double) NNODES - nb);
        for(i=0; i<n->ndend; i++) phys_lengths( n->dendlist[i].root, &h );

        if ( !n ) {
            printf("\n Failed to find sampled neuron\n");
            return(1);
        }
    } else {
        printf("\n Failed to find sampled neuron\n");
        return(1);
    }

/*  Enumerate Nodes */

    FirstNode = 0;
    for ( k=0 ; k<n->ndend ; k++ ) enumerate_nodes( n->dendlist[k].root, &FirstNode );
    for ( k=0 ; k<n->ndend ; k++ ) n->dendlist[k].root->junct = FirstNode;
    printf("Number of nodes is %d\n", FirstNode+1);

/* Construct Sparse Matrices */

    counter = 0;
    nodes = FirstNode+1;
    mat_malloc( &lhs, nodes, 3*nodes-2 );
    mat_malloc( &rhs, nodes, 3*nodes-2 );
    lhs.n = rhs.n = nodes;
    lhs.a[3*nodes-3] = rhs.a[3*nodes-3] = 0.0;
    lhs.start_row[0] = rhs.start_row[0] = 0;

    for ( k=0 ; k<n->ndend ; k++ ) {
        bnow = n->dendlist[k].root;
        Generate_Dendrite( bnow, &counter);
    }

    for ( k=0 ; k<n->ndend ; k++ ) {
        bnow = n->dendlist[k].root;
        lhs.a[counter] = 0.5*bnow->d[1]*bnow->pl[1];
        rhs.a[counter] = -(bnow->d[0])*(bnow->d[1])/bnow->pl[1];
        lhs.col[counter] = rhs.col[counter] = bnow->first;
        lhs.a[3*nodes-3] += 1.5*(bnow->d[0])*bnow->pl[1];
        rhs.a[3*nodes-3] += (bnow->d[0])*(bnow->d[1])/bnow->pl[1];
        counter++;
    }
    lhs.col[3*nodes-3] = rhs.col[3*nodes-3] = nodes-1;
    lhs.start_row[nodes] = rhs.start_row[nodes] = 3*nodes-2;

    dt = 1.0/((double) NT);

    /* Constructing Conductance Profile */
    cond = cond_profile( dt, TAU, 5.e-5, GMAX);
    for ( k=0 ; k<n->ndend ; k++ ) {
        bnow = n->dendlist[k].root;
        Initialise_Contact( bnow, cond);
    }

    for ( i=0 ; i<3*nodes-2 ; i++ ) {
        rhs.a[i] *= GA;
        rhs.a[i] += GM*lhs.a[i];
        lhs.a[i] *= CM;
        rhs.a[i] *= 0.5*dt;
        lhs.a[i] += rhs.a[i];
        rhs.a[i] = lhs.a[i] - 2.0*rhs.a[i];
    }
    lhs.a[3*nodes-3] += (4.0*as/pi)*CS;
    rhs.a[3*nodes-3] += (4.0*as/pi)*CS;

/*  Store Lhs and Rhs matrices */
    LhsStore = (double *) malloc( (3*nodes-2)*sizeof(double) );
    RhsStore = (double *) malloc( (3*nodes-2)*sizeof(double) );
    for ( i=0 ; i<3*nodes-2 ; i++ ) {
        LhsStore[i] = lhs.a[i];
        RhsStore[i] = rhs.a[i];
    }

    v = (double *) malloc( (nodes)*sizeof(double) );
    x = (double *) malloc( (nodes)*sizeof(double) );
    y = (double *) malloc( (nodes)*sizeof(double) );

    for ( i=0 ; i<nodes ; i++ ) v[i] = 0.0;

/*  Initialise temporal integration */

    tnow = 0.0;
    tout = DT;
    vold = vnew = vmid = 0.0;
    first = 1;
    nspk = 0;

    while ( tnow < TEND ) {
        tnow += dt;
        vs = veq+v[nodes-1];

        mold = mnew;
        nold = nnew;
        hold = hnew;

        tmp1 = dt*alfa_m(vs);
        tmp2 = dt*beta_m(vs);
        tmp = tmp1+tmp2;
        mnew = mold+(tmp1-mold*tmp)/(1.0+0.5*tmp);
        tmp1 = dt*alfa_n(vs);
        tmp2 = dt*beta_n(vs);
        tmp = tmp1+tmp2;
        nnew = nold+(tmp1-nold*tmp)/(1.0+0.5*tmp);
        tmp1 = dt*alfa_h(vs);
        tmp2 = dt*beta_h(vs);
        tmp = tmp1+tmp2;
        hnew = hold+(tmp1-hold*tmp)/(1.0+0.5*tmp);

        gval = g_na*hnew*pow(mnew,3);
        cval = v_na*gval;
        tmp = g_k*pow(nnew,4);
        gval += tmp;
        cval += tmp*v_k;
        gval += g_l;
        cval += g_l*v_l;

        tmp = as*4.0*dt/pi;
        gval *= 0.5*tmp*fact;
        cval *= tmp*fact;

        for ( i=0 ; i<nodes ; i++ ) x[i] = 0.0;
        tmp = 2.0*dt/pi;
        for ( i=0 ; i<n->ndend ; i++ ) Assign_Current(n->dendlist[i].root, x, tmp);
        lhs.a[3*nodes-3] += gval;
        rhs.a[3*nodes-3] -= gval;

        LU_Factor(&lhs, &first);

        mat_vec_mult(&rhs,v,y);
        for ( i=0 ; i<nodes ; i++ ) x[i] += y[i];
        x[nodes-1] -= 4.0*dt*SIN/pi;
        x[nodes-1] += cval;

        LU_Solve( &lhs, v, x );

        vnew = v[nodes-1];

        if ( vmid > SENTRY ) {
            if ( vold < vmid && vnew < vmid ) {
                tmp1 = tnow + 0.5*dt*(vnew-vold)/(2.0*vmid-vold-vnew);
                tmp2 = vmid + 0.125*pow(vnew-vold,2)/(2.0*vmid-vold-vnew);
                if ( nspk == 0 ) fp = fopen("out.res", "w");
                if ( nspk != 0 ) fp = fopen("out.res", "a");
                nspk++;
                fprintf(fp,"Spike %d of %lf mv at time %lf ms \n", nspk, tmp2, tmp1);
                fclose(fp);
            }
        }
        vold = vmid;
        vmid = vnew;

        if ( tnow > tout ) {
            printf("\rReached time %5.1lf ms\t", tout);
            printf("\nNumerical Voltage %12.6lf mV\n",v[nodes-1]);
            tout += DT;
        }

        for ( i=0 ; i<3*nodes-2 ; i++ ) {
            lhs.a[i] = LhsStore[i];
            rhs.a[i] = RhsStore[i];
        }
    }

    Destroy_Sampled_Neuron( n );
    return(0);
}


 /************************************************************
                 Function to assign current
 **************************************************************/
void Assign_Current(branch *bnow, double *x, double fac )
{
    int i, j, node;
    contact *con;
    double InterSpikeInterval, tmp, gnow;
    extern sparse_mat lhs, rhs;
    extern double veq;

    if (bnow->child != NULL ) Assign_Current(bnow->child, x, fac );
    if (bnow->peer != NULL ) Assign_Current(bnow->peer, x, fac );

    con = bnow->conlist;
    while ( con != NULL ) {
        if ( con->hist < 0 ) {
            (con->hist)++;
        } else if ( con->hist == (con->cp->n)-1 ) {
            InterSpikeInterval = DTIME-(1000.0/RATE-DTIME)*log(ran(&ix, &iy, &iz));
            InterSpikeInterval /= (con->cp->dt);
            if ( fmod(InterSpikeInterval,1.0) > 0.5 ) {
                con->hist = -((int) ceil(InterSpikeInterval));
            } else {
                con->hist = -((int) floor(InterSpikeInterval));
            }
        } else {
            node = con->xn;
            j = lhs.start_row[node];
            while ( node != lhs.col[j] ) j++;
            gnow = con->cp->g[con->hist];
            rhs.a[j] -= fac*gnow;
            x[node] = gnow;
            (con->hist)++;
            gnow = con->cp->g[con->hist];
            lhs.a[j] += fac*gnow;
            x[node] = fac*(x[node]+gnow)*(con->vsyn);
        }
        con = con->next;
    }
    return;
}


 /******************************************************
           Function to initialise a contact
  ******************************************************/
void Initialise_Contact( branch *bnow, syn_cond *cond)
{
    int i, k;
    extern unsigned long int ix, iy, iz;
    extern double veq;
    double InterSpikeInterval, rate=RATE, dtime=DTIME;
    branch *btmp;
    contact *con;

/*  Step 1 - Recurse to the end of the dendrite */
    if ( bnow->child != NULL ) Initialise_Contact( bnow->child, cond);
    if ( bnow->peer != NULL ) Initialise_Contact( bnow->peer, cond);

/*  Step 2 - Fill in matrix entries for terminal points */
    con = bnow->conlist;
    while ( con != NULL ) {
        con->cp = cond;
        con->deadtime = dtime;
        con->vsyn = VSYN-veq;
        if ( rate*dtime >= 1000.0 ) {
            printf("\nImpossible rate requested - Deadtime ignored!");
            InterSpikeInterval = -(1000.0/rate)*log(ran(&ix, &iy, &iz));
        } else {
            InterSpikeInterval = dtime-(1000.0/rate-dtime)*log(ran(&ix, &iy, &iz));
        }
        InterSpikeInterval /= cond->dt;
        if ( fmod(InterSpikeInterval,1.0) > 0.5 ) {
            con->hist = -((int) ceil(InterSpikeInterval));
        } else {
            con->hist = -((int) floor(InterSpikeInterval));
        }
        con = con->next;
    }
    return;
}


/* Gets current profile */
syn_cond *cond_profile( double dt,   /* Integration time step */
                        double tau,  /* Synapse time constant */
                        double tol,  /* Significance indicator */
                        double gmax  /* Max input conductance (mS/cm^2) */)
{
    int i;
    syn_cond *out;
    double tmp, told, tnew;

    out = (syn_cond *) malloc( sizeof(syn_cond) );
    out->dt = dt;
    out->tau = tau;
    out->tol = tol;
    out->gmax = gmax;

/* Iterate to find duration of conductance pulse */
    tmp = 1.0-log(tol);
    tnew = tmp;
    do {
        told = tnew;
        tnew = tmp+log(told);
    } while ( fabs(tnew-told)>5.e-11 );
    out->n = ((int) ceil(tau*tnew/dt) );
    out->g = (double *) malloc( (out->n)*sizeof(double) );
    for ( i=0 ; i<(out->n) ; i++ ) {
        tmp = (dt*((double) i)/tau);
        out->g[i] = gmax*tmp*exp(1.0-tmp);
    }
    return out;
}


/* Free current profile */
void free_cond_profile( syn_cond *cond)
{
    free(cond->g);
    free(cond);
    return;
}



 /******************************************************
        Function to constuct sparse matrices
  ******************************************************/

void Generate_Dendrite( branch *bnow, int *counter)
{
    int i, k;
    extern sparse_mat lhs, rhs;
    branch *btmp;
    contact *con;
    double SumL, SumR;

/* Step 1 - Recurse to the end of the dendrite */

    if ( bnow->child != NULL ) Generate_Dendrite( bnow->child, counter);
    if ( bnow->peer != NULL ) Generate_Dendrite( bnow->peer, counter);

    for ( k=bnow->first-bnow->nobs+2,i=bnow->nobs-1 ; i>0 ; i--,k++ ) {

/* Step 2 - Fill in matrix entries for terminal points */

        if ( bnow->child == NULL && i == bnow->nobs - 1 ) {

            lhs.a[*counter] = 1.5*(bnow->d[i])*(bnow->pl[i]-bnow->pl[i-1]);
            rhs.a[*counter] = (bnow->d[i-1]*bnow->d[i])/(bnow->pl[i]-bnow->pl[i-1]);
            lhs.col[*counter] = rhs.col[*counter] = k;
            (*counter)++;

            lhs.a[*counter] = 0.5*(bnow->d[i-1])*(bnow->pl[i] - bnow->pl[i-1]);
            rhs.a[*counter] = -(bnow->d[i-1]*bnow->d[i])/(bnow->pl[i] - bnow->pl[i-1]);

            if ( k == bnow->first ) {
                lhs.col[*counter] = rhs.col[*counter] = bnow->junct;
            } else {
                lhs.col[*counter] = rhs.col[*counter] = k + 1;
            }

            (*counter)++;
            lhs.start_row[k+1] = rhs.start_row[k+1] = *counter;

/* Step 3 - Fill in matrix entries for branch points */

        } else if ( bnow->child != NULL && i == bnow->nobs - 1 ) {

            btmp = bnow->child;
            SumR = SumL = 0.0;

            while ( btmp != NULL ) {

                lhs.a[*counter] = 0.5*btmp->d[1]*btmp->pl[1];
                rhs.a[*counter] = -(btmp->d[0]*btmp->d[1])/btmp->pl[1];
                lhs.col[*counter] = rhs.col[*counter] = btmp->first;
                (*counter)++;
                SumL += 1.5*btmp->d[0]*btmp->pl[1];
                SumR += (btmp->d[0]*btmp->d[1])/btmp->pl[1];
                btmp = btmp->peer;
            }

            lhs.a[*counter] = SumL+1.5*bnow->d[i]*(bnow->pl[i]-bnow->pl[i-1]);
            rhs.a[*counter] = SumR+(bnow->d[i-1]*bnow->d[i])/(bnow->pl[i]-bnow->pl[i-1]);
            lhs.col[*counter] = rhs.col[*counter] = k;
            (*counter)++;

            lhs.a[*counter] = 0.5*(bnow->d[i-1])*(bnow->pl[i]-bnow->pl[i-1]);
            rhs.a[*counter] = -(bnow->d[i-1]*bnow->d[i])/(bnow->pl[i] - bnow->pl[i-1]);

            if ( k == bnow->first ) {
                lhs.col[*counter] = rhs.col[*counter] = bnow->junct;
            } else {
                lhs.col[*counter] = rhs.col[*counter] = k + 1;
            }

            (*counter)++;
            lhs.start_row[k+1] = rhs.start_row[k+1] = *counter;

        } else {

    /* Step 4 - Fill in matrix entries for internal point */

            lhs.a[*counter] = 0.5*(bnow->d[i+1])*(bnow->pl[i+1] - bnow->pl[i]);
            rhs.a[*counter] = -(bnow->d[i]*bnow->d[i+1])/(bnow->pl[i+1] - bnow->pl[i]);
            lhs.col[*counter] = rhs.col[*counter] = k - 1;
            (*counter)++;

            lhs.a[*counter] = 1.5*(bnow->d[i])*(bnow->pl[i+1] - bnow->pl[i-1]);
            rhs.a[*counter] = (bnow->d[i-1]*bnow->d[i])/(bnow->pl[i] - bnow->pl[i-1])
                            + (bnow->d[i]*bnow->d[i+1])/(bnow->pl[i+1] - bnow->pl[i]);
            lhs.col[*counter] = rhs.col[*counter] = k;
            (*counter)++;

            lhs.a[*counter] = 0.5*(bnow->d[i-1])*(bnow->pl[i] - bnow->pl[i-1]);
            rhs.a[*counter] = -(bnow->d[i-1]*bnow->d[i])/(bnow->pl[i] - bnow->pl[i-1]);
            lhs.col[*counter] = rhs.col[*counter] = k + 1;

            if ( k == bnow->first ) {
                lhs.col[*counter] = rhs.col[*counter] = bnow->junct;
            } else {
                lhs.col[*counter] = rhs.col[*counter] = k + 1;
            }

            (*counter)++;

            lhs.start_row[k+1] = rhs.start_row[k+1] = *counter;
        }
    }

    /*  Step 5 - Assign enumerated node to contact */

    con = bnow->conlist;
    while ( con != NULL ) {
        if ( con->xn == 0 ) {
            con->xn = bnow->junct;
        } else {
            con->xn = bnow->first+1-con->xn;
        }
        con = con->next;
    }
    return;
}


neuron *Load_Sampled_Neuron(char *filename)
{
    int j, k, ncon, n, id, connected, ignored;
    double tmp, piby2, xold, xnew, yold, ynew, zold, znew, diam,
           xl, xr, dl, dr, px, py, pz, min;
    neuron *cell;
    soma *s;
    contact *oldcon, *newcon, *firstcon;
    branch *bold, *bnew, *first_branch, *bopt;
    char temp[100];
    FILE *input;

/*  STEP 1. - Open neuron data file */
    printf("\nOpening file %s\n",filename);
    if ( (input=fopen(filename,"r"))==NULL ) return NULL;

/*  STEP 2. - Get memory for neuron structure */
    cell = (neuron *) malloc( sizeof(neuron) );
    s = cell->s = (soma *) malloc( sizeof(soma) );

/*  STEP 3, - Get information about soma */
    fscanf(input,"%s",temp);
    fscanf(input,"%d",&n);
    printf("Identified a %s defined by %d data\n", temp, n);

/*  STEP 4. - Initialise soma structure */
    s->x = (double *) malloc( n*sizeof(double) );
    s->y = (double *) malloc( n*sizeof(double) );
    s->z = (double *) malloc( n*sizeof(double) );
    s->d = (double *) malloc( n*sizeof(double) );
    s->ga = GA;
    s->cs = CS;
    s->conlist = NULL;
    s->ncon = 0;
    printf("%s initialised\n", temp);

/*  STEP 5. - Get soma morphological data */
    fscanf(input,"%lf %lf %lf %lf",&xold, &yold, &zold, &diam);
    s->x[0] = xold;
    s->y[0] = yold;
    s->z[0] = zold;
    s->d[0] = diam;
    for ( k=j=1,s->p_len=0.0 ; k<n ; k++ ) {
        fscanf(input,"%lf %lf %lf %lf",&xnew, &ynew, &znew, &diam);
        tmp = pow(xnew-xold,2)+pow(ynew-yold,2)+pow(znew-zold,2);
        if ( tmp > 0.01 ) {
            s->x[j] = xold = xnew;
            s->y[j] = yold = ynew;
            s->z[j] = zold = znew;
            s->d[j] = diam;
            j++;
            s->p_len += sqrt(tmp);
        }
    }
    s->nobs = j;
    printf("Finished soma\n");

/*  STEP 6. - Get branch and contact data */
    firstcon = oldcon = NULL;
    bold = NULL;
    while  ( fscanf(input,"%s", temp)!=EOF ) {
        if ( strcmp(temp, "Branch") == 0 || strcmp(temp, "branch") == 0 ) {
            fscanf(input, "%d", &n);
            printf("Found a branch defined by %d nodes\n", n);
            bnew = (branch *) malloc( sizeof(branch) );
            if ( bold != NULL) {
                bold->child = bnew;
            } else {
                first_branch = bnew;
            }
            bnew->parent = bold;
            bnew->peer = NULL;
            bnew->child = NULL;

/*  STEP 7. - Initialise branch */
            init_branch( bnew, n, CM, GM, GA);
            fscanf(input,"%lf %lf %lf %lf", &xold, &yold, &zold, &diam);
            bnew->x[0] = xold;
            bnew->y[0] = yold;
            bnew->z[0] = zold;
            bnew->d[0] = diam;
            for ( bnew->p_len=bnew->pl[0]=0.0,k=j=1 ; k<n ; k++ ) {
                fscanf(input,"%lf %lf %lf %lf", &xnew, &ynew, &znew, &diam);
                tmp = pow(xnew-xold,2)+pow(ynew-yold,2)+pow(znew-zold,2);
                if ( tmp > 0.01 ) {
                    bnew->p_len += sqrt(tmp);
                    bnew->pl[j] = bnew->p_len;
                    bnew->x[j] = xold = xnew;
                    bnew->y[j] = yold = ynew;
                    bnew->z[j] = zold = znew;
                    bnew->d[j] = diam;
                    j++;
                }
            }
            bnew->nobs = j;
            bold = bnew;
        } else if ( strcmp(temp, "Marker") == 0 || strcmp(temp, "marker") == 0 ) {

/*  STEP 8. - Initialise marker */
            fscanf(input, "%d %d", &id, &n);
            printf("Found and initialised %d branch contacts of type %d\n", n, id);
            for ( k=0 ; k<n ; k++ ) {
                newcon = (contact *) malloc( sizeof(contact) );
                newcon->next = NULL;
                newcon->sd = NULL;
                newcon->id = id;
                newcon->cp = NULL;
                if ( oldcon != NULL ) {
                    oldcon->next = newcon;
                } else {
                    firstcon = newcon;
                }
                fscanf(input,"%lf %lf %lf %lf", &newcon->xc, &newcon->yc, &newcon->zc, &newcon->dc );
                oldcon = newcon;
            }
        } else {
            printf("Unknown block type %s!\n", temp);
            return NULL;
        }
    }
    fclose(input);

/*  STEP 9. - Associate contacts with branches and soma */
    ignored = 0;
    while ( firstcon ) {
        bold = first_branch;
        bopt = NULL;
        while ( bold ) {
            BuildContactInfo( firstcon, bold, &bopt);
            bold = bold->child;
        }
        newcon = firstcon->next;
        if ( firstcon->sd > 4.0 ) {

/*  STEP 10a. - Check for proximity to soma */
            px = firstcon->xc;
            py = firstcon->yc;
            pz = firstcon->zc;

/*  STEP 1. - First stage is different from others */
            xnew = s->x[0]; ynew = s->y[0]; znew = s->z[0];
            firstcon->sd = min = sqrt(pow(xnew-px,2)+pow(ynew-py,2)+pow(znew-pz,2))-(s->d[0]);

/*  STEP 2. - Second stage compares points and projected points */
            for ( k=1 ; k<s->nobs ; k++ ) {
                xnew = s->x[k]; ynew = s->y[k]; znew = s->z[k];
                min = sqrt(pow(xnew-px,2)+pow(ynew-py,2)+pow(znew-pz,2))-(s->d[k]);
                if ( min < firstcon->sd ) firstcon->sd = min;
            }
            if ( firstcon->sd < 4.0 ) {
                oldcon = s->conlist;
                if ( oldcon ) {
                    while ( oldcon->next ) oldcon = oldcon->next;
                    oldcon->next = firstcon;
                } else {
                    s->conlist = firstcon;
                }
                firstcon->next = NULL;
                s->ncon++;
            } else {
                free(firstcon);
                ignored++;
            }
        } else {

/*  STEP 3. - Add contact to conlist of optimal branch */
            oldcon = bopt->conlist;
            if ( oldcon ) {
                while ( oldcon->next ) oldcon = oldcon->next;
                oldcon->next = firstcon;
            } else {
                bopt->conlist = firstcon;
            }
            firstcon->next = NULL;
            bopt->ncon++;
        }
        firstcon = newcon;
    }
    printf("\nIgnored %d contacts", ignored);

/*  STEP 11. - Count dendritic branches at soma */
    bold = first_branch;
    n = 0;
    while ( bold ) {
        bnew = first_branch;
        do {
            k = bnew->nobs-1;
            tmp = pow(bold->x[0]-bnew->x[k],2)+
                  pow(bold->y[0]-bnew->y[k],2)+
                  pow(bold->z[0]-bnew->z[k],2);
            connected = ( tmp < 0.01 );
            bnew = bnew->child;
        } while ( bnew && !connected );
        if ( !connected ) n++;
        bold = bold->child;
    }
    cell->ndend = n;
    printf("\n\nTree contains %d individual dendrite(s) ...\n", n);

/* STEP 12. - Identify somal dendrites but extract nothing */
    cell->dendlist = (dendrite *) malloc( (cell->ndend)*sizeof(dendrite) );
    bold = first_branch;
    n = 0;
    while ( n < cell->ndend ) {
        bnew = first_branch;
        do {
            k = bnew->nobs-1;
            tmp = pow(bold->x[0]-bnew->x[k],2)+
                  pow(bold->y[0]-bnew->y[k],2)+
                  pow(bold->z[0]-bnew->z[k],2);
            connected = ( tmp < 0.01 );
            bnew = bnew->child;
        } while ( bnew && !connected );
        if ( !connected ) {
            cell->dendlist[n].root = bold;
            n++;
        }
        bold = bold->child;
    }

/*  STEP 13. - Extract root of each dendrite from dendrite list */
    for ( k=0 ; k<cell->ndend ; k++ ) {
        bold = cell->dendlist[k].root;
        remove_branch( &first_branch, bold);
    }

/*  STEP 14. - Build each dendrite from its root branch */
    for ( k=0 ; k<cell->ndend ; k++ ) {
        build_dendrite( &first_branch, cell->dendlist[k].root);
        clean_dendrite( cell->dendlist[k].root);
    }
    if ( first_branch != NULL ) printf("\nWarning: Unconnected branch segments still exist\n");
    return cell;
}



/****************************************************************

                Function to initialise a BRANCH

****************************************************************/

void init_branch( branch *b, int n, double cm, double gm, double ga)
{

/*  Allocate memory for spatial orientation of dendrite */
    b->x = (double *) malloc( n*sizeof(double) );
    b->y = (double *) malloc( n*sizeof(double) );
    b->z = (double *) malloc( n*sizeof(double) );

/*  Allocate memory for branch geometry */
    b->d = (double *) malloc( n*sizeof(double) );
    b->pl = (double *) malloc( n*sizeof(double) );

/*  Set parameter values */
    b->cm = cm;
    b->gm = gm;
    b->ga = ga;

/*  Initialise node information */
    b->nodes = 0;
    b->first = 0;
    b->junct = 0;

/*  Initialise contact information */
    b->conlist = NULL;
    b->ncon = 0;
    return;
}

/**************************************************************

                Function to initialise a TEST BRANCH

***************************************************************/

void init_test_branch( branch *b, int n, double cm, double gm,double ga)
{

/*  Allocate memory for spatial orientation of dendrite */
    b->x = (double *) malloc( 2*sizeof(double) );
    b->y = (double *) malloc( 2*sizeof(double) );
    b->z = (double *) malloc( 2*sizeof(double) );

/*  Allocate memory for branch geometry */
    b->d = (double *) malloc( n*sizeof(double) );
    b->pl = (double *) malloc( n*sizeof(double) );

/*  Set parameter values */
    b->cm = cm;
    b->gm = gm;
    b->ga = ga;

/*  Initialise node information */
    b->nodes = 0;
    b->first = 0;
    b->junct = 0;

/*  Initialise contact information */
    b->conlist = NULL;
    b->ncon = 0;
    return;
}

/***************************************************************

             Function to build CONTACT information

***************************************************************/

void BuildContactInfo(contact *con, branch *b, branch **bopt)
{
    int k;
    double px, py, pz, tmp, xold, xnew, yold, ynew, zold, znew,
           numer, denom, xmin, ymin, zmin, min;

    px = con->xc;
    py = con->yc;
    pz = con->zc;

/*  STEP 1. - First stage is different from others */
    xnew = b->x[0]; ynew = b->y[0]; znew = b->z[0];
    min = sqrt(pow(xnew-px,2)+pow(ynew-py,2)+pow(znew-pz,2));
    if ( con->sd == NULL || ( con->sd != NULL && min < con->sd ) ) {
        con->sd = min;
        con->xp = xnew; con->yp = ynew; con->zp = znew;
        con->pl = 0.0;
        *bopt = b;
    }

/*  STEP 2. - Second stage compares points and projected points */
    for ( k=1 ; k<b->nobs ; k++ ) {

        xold = xnew; yold = ynew; zold = znew;
        xnew = b->x[k]; ynew = b->y[k]; znew = b->z[k];
        numer = (xnew-xold)*(px-xold)+(ynew-yold)*(py-yold)+(znew-zold)*(pz-zold);
        denom = pow(xnew-xold,2)+pow(ynew-yold,2)+pow(znew-zold,2);

/*  STEP 2a. - Project onto branch */
        if ( 0.0 <= numer && numer <= denom ) {
            tmp = numer/denom;
            xmin = (1.0-tmp)*xold+tmp*xnew;
            ymin = (1.0-tmp)*yold+tmp*ynew;
            zmin = (1.0-tmp)*zold+tmp*znew;
            min = sqrt(pow(xmin-px,2)+pow(ymin-py,2)+pow(zmin-pz,2));
            if ( !(con->sd) || ( con->sd && min < con->sd ) ) {
                con->sd = min;
                con->xp = xmin; con->yp = ymin; con->zp = zmin;
                con->pl = (1.0-tmp)*b->pl[k-1]+tmp*b->pl[k];
                *bopt = b;
            }
        }

/*  STEP 2b. - Check proximity to points of branch */
        min = sqrt(pow(xnew-px,2)+pow(ynew-py,2)+pow(znew-pz,2));
        if ( !(con->sd) || ( con->sd && min < con->sd ) ) {
            con->sd = min;
            con->xp = xnew; con->yp = ynew; con->zp = znew;
            con->pl = b->pl[k];
            *bopt = b;
        }
    }
    return;
}

/***************************************************************

          Function to remove a branch from a branch list

****************************************************************/

void remove_branch(branch **head, branch *b)
{
    if ( *head == NULL || b == NULL ) return;
    if ( *head == b ) {
        *head = b->child;
        if ( *head != NULL )  (*head)->parent = NULL;
    } else {
        b->parent->child = b->child;
        if ( b->child != NULL ) b->child->parent = b->parent;
    }
    b->parent = NULL;
    b->child = NULL;
    return;
}

/***************************************************************

          Function to build a dendrite from its root

***************************************************************/

void build_dendrite( branch **head, branch *root)
{
    int k;
    double tmp;
    branch *bnow, *bnext, *btmp;

    bnow = *head;
    while ( bnow != NULL ) {

/* Store bnow's child in case it's corrupted */
        bnext = bnow->child;

/* Search if proximal end of bnow is connected to distal end of root */
        k = (root->nobs)-1;
        tmp = pow(bnow->x[0]-root->x[k],2)+pow(bnow->y[0]-root->y[k],2)+
              pow(bnow->z[0]-root->z[k],2);
        if ( tmp <= 0.01 ) {

/* Take bnow out of the branch list */
            remove_branch( head, bnow);

/* Connect bnow to the root as the child or a peer of the child.
   Initialise childs' children and peers to NULL as default */
            bnow->child = NULL;
            bnow->peer = NULL;
            bnow->parent = root;

/* Inform root about its child if it's the first child, or add
   new child to first child's peer list */
            if ( root->child != NULL ) {
                btmp = root->child;
                while ( btmp->peer != NULL ) btmp = btmp->peer;
                btmp->peer = bnow;
            } else {
                root->child = bnow;
            }
        }

/* Initialise bnow to next branch in list */
        bnow = bnext;
    }

/* Iterate through remaining tree */
    if ( root->child ) build_dendrite( head, root->child);
    if ( root->peer ) build_dendrite( head, root->peer);
    return;
}

/***************************************************************

          Function to remove peerless children

***************************************************************/
void clean_dendrite( branch *root)
{
    int k, np, nc, mem, n;
    double tmp, sarea;
    contact *con;
    branch *btmp, *brem;

/*  Iterate through remaining tree */
    if ( root->child != NULL ) clean_dendrite( root->child );
    if ( root->peer != NULL ) clean_dendrite( root->peer );

/*  Extend original parent limb */
    brem = root->child;
    if ( brem != NULL && brem->peer == NULL ) {
        root->child = brem->child;
        if ( brem->child ) {
            brem->child->parent = root;
            btmp = brem->child->peer;
            while ( btmp ) {
                btmp->parent = root;
                btmp = btmp->peer;
            }
        }
//        root->nodes += (brem->nodes)-1;
        np = root->nobs;
        nc = brem->nobs;
        mem = np+nc-1;
        root->nobs = mem;
        root->x = (double *) realloc( (void *)root->x, mem*sizeof(double) );
        for ( k=np ; k<mem ; k++ ) root->x[k] = brem->x[k-np+1];
        root->y = (double *) realloc( (void *)root->y, mem*sizeof(double) );
        for ( k=np ; k<mem ; k++ ) root->y[k] = brem->y[k-np+1];
        root->z = (double *) realloc( (void *)root->z, mem*sizeof(double) );
        for ( k=np ; k<mem ; k++ ) root->z[k] = brem->z[k-np+1];
        root->d = (double *) realloc( (void *)root->d, mem*sizeof(double) );
        for ( k=np ; k<mem ; k++ ) root->d[k] = brem->d[k-np+1];
        root->pl = (double *) realloc( (void *)root->pl, mem*sizeof(double) );
        for ( k=np ; k<mem ; k++ ) root->pl[k] = root->p_len+brem->pl[k-np+1];
        root->ncon += brem->ncon;
        con = root->conlist;
        if ( con ) {
            while ( con->next ) con = con->next;
            con->next = brem->conlist;
        } else {
            root->conlist = brem->conlist;
        }
        brem->conlist = NULL;
        free(brem->x);
        free(brem->y);
        free(brem->z);
        free(brem->d);
        free(brem->pl);
        free ( brem );
    }
    return;
}

/***************************************************************

                Function to destroy a NEURON

****************************************************************/

void Destroy_Sampled_Neuron(neuron *cell)
{
    int i;
    contact *prevcon, *nextcon;

/*  Free Soma */

    if ( cell->s != NULL) {
        if ( cell->s->x ) free ( cell->s->x );
        if ( cell->s->y ) free ( cell->s->y );
        if ( cell->s->z ) free ( cell->s->z );
        if ( cell->s->d ) free ( cell->s->d );
        prevcon = cell->s->conlist;
        while ( prevcon ) {
            nextcon = prevcon->next;
            free ( prevcon );
            prevcon = nextcon;
        }
        free ( cell->s );
    }
    for ( i=0 ; i<cell->ndend ; i++ ) destroy_dendrite( cell->dendlist[i].root );
    free(cell);
    return;
}

/***************************************************************

               Function to destroy DENDRITE

****************************************************************/

void destroy_dendrite( branch *b )
{
    int i;
    contact *prevcon, *nextcon;

    if ( b->child ) destroy_dendrite(b->child);
    if ( b->peer ) destroy_dendrite(b->peer);
    free(b->x);
    free(b->y);
    free(b->z);
    free(b->d);
    free(b->pl);
    prevcon = b->conlist;
    while ( prevcon ) {
        nextcon = prevcon->next;
        free ( prevcon );
        prevcon = nextcon;
    }
    free ( b );
    return;
}

/***************************************************************

             Function to count number of branches

***************************************************************/

int count_branches( branch *bstart, branch *bnow)
{
    static int n;

    if ( bstart == bnow ) n = 0;
    if ( bnow ) {
        if ( bnow->child ) count_branches(bstart, bnow->child);
        if ( bnow->peer ) count_branches(bstart, bnow->peer);
        n++;
    }
    return n;
}

/***************************************************************

             Function to find length of dendrite from
                    current branch to tips.

****************************************************************/

double branch_length( branch *bstart, branch *bnow)
{
    static double length;

    if ( bstart == bnow ) length = 0.0;
    if ( bnow ) {
        if ( bnow->child ) branch_length(bstart, bnow->child);
        if ( bnow->peer ) branch_length(bstart, bnow->peer);
        length += bnow->p_len;
    }
    return length;
}




/**********************************************************

        Function to enumerate the nodes on a dendrite

**********************************************************/

void enumerate_nodes(branch *bnow, int *FirstNode )
{
    branch *btmp;

    if ( (bnow->child) != NULL ) enumerate_nodes( bnow->child, FirstNode );
    if ( (bnow->peer) != NULL ) enumerate_nodes( bnow->peer, FirstNode );

    if ( bnow->child != NULL ) {
        btmp = bnow->child;
        while( btmp != NULL ){
            btmp->junct = *FirstNode;
            btmp = btmp->peer;
        }
    }

    bnow->first = *FirstNode + bnow->nobs - 2;
    *FirstNode += (bnow->nobs)-1;
    return;
}


/*************************************************************

    Multiplies sparse matrix a[ ][ ] with vector v[ ]

*************************************************************/

void mat_vec_mult( sparse_mat *a, double *v , double *b)
{
    int i, j, ntop, n;

    n = a->n;
    for ( j=0 ; j<n ; j++) {
        ntop = a->start_row[j+1];
        for( b[j]=0.0,i=(a->start_row[j]) ; i<ntop ; i++ ) {
            b[j] += (a->a[i])*v[a->col[i]];
        }
    }
    return;
}


/*******************************************************************

    Allocates memory to a sparse matrix - returns 0 if successful

*******************************************************************/

void mat_malloc( sparse_mat *a, int n, int w)
{
    a->a = (double *) malloc( w*sizeof(double) );
    a->col = (int *) malloc( w*sizeof(int) );
    a->start_row = (int *) malloc( (n+1)*sizeof(int) );
    a->n = n;
    a->l = malloc(sizeof(sparse_mat));
    a->u = malloc(sizeof(sparse_mat));
    a->l->a = (double *) malloc( (2*n-1)*sizeof(double) );
    a->l->col = (int *) malloc( (2*n-1)*sizeof(int) );
    a->l->start_row = (int *) malloc( (n+1)*sizeof(int) );
    a->l->n = n;
    a->u->a = (double *) malloc( (2*n-1)*sizeof(double) );
    a->u->col = (int *) malloc( (2*n-1)*sizeof(int) );
    a->u->start_row = (int *) malloc( (n+1)*sizeof(int) );
    a->u->n = n;

    if ( !a->start_row ) exit(1);

}


/* De-allocates memory of a sparse matrix */
void mat_free( sparse_mat *a)
{
    free(a->a);
    free(a->col);
    free(a->start_row);
    free(a);
}



 /****************************************************
        Recursive function to insert test contacts
  ****************************************************/
void InsertContacts( branch *b)
{
    int j, k, xn;
    double old, new;
    contact *con;

    if ( b->child != NULL ) InsertContacts(b->child);
    if ( b->peer != NULL ) InsertContacts(b->peer);

    con = b->conlist;
    while ( con != NULL ) {
        xn = -1;
        for ( k=0 ; k<b->nobs ; k++ ) {
            if ( fabs(b->pl[k]-con->pl) < 5.e-11 ) xn = k;
        }
        if ( xn > -1 ) {
            con->xn = xn;
        } else {
            (b->nobs)++;
            b->pl = (double *) realloc( b->pl, (b->nobs)*sizeof(double) );
            k = 0;
            while ( b->pl[k] < con->pl ) k++;
            new = b->pl[k];
            b->pl[k] = con->pl;
            con->xn = k;
            for ( j=k+1 ; j<b->nobs ; j++ ) {
                old = new;
                new = b->pl[j];
                b->pl[j] = old;
            }
        }
        con = con->next;
    }
    return;
}



 /***************************************************************
            Function To Load A Test Neuron
 ***************************************************************/
neuron *Load_Test_Neuron(char *filename)
{
    int j, k, ncon, n, id, connected, ignored;
    double tmp, piby2, xl, xr, dl, dr, px, py, pz, min, radius, dx;
    neuron *cell;
    contact *newcon;
    branch *bold, *bnew, *first_branch, *bopt;
    char temp[100];
    FILE *input;

/*  STEP 1. - Open neuron data file */
    printf("\nOpening file %s\n",filename);
    if ( (input=fopen(filename,"r"))==NULL ) return NULL;

/*  STEP 2. - Get memory for neuron structure */
    cell = (neuron *) malloc( sizeof(neuron) );

/*  STEP 3. - Get branch and contact data */
    bold = NULL;
    while  ( fscanf(input,"%s", temp)!=EOF ) {
        if ( strcmp(temp, "Branch") == 0 || strcmp(temp, "branch") == 0 ) {
            fscanf(input, "%d", &n);
            printf("Found a branch\n");
            bnew = (branch *) malloc( sizeof(branch) );
            if ( bold != NULL) {
                bold->child = bnew;
            } else {
                first_branch = bnew;
            }
            bnew->parent = bold;
            bnew->peer = NULL;
            bnew->child = NULL;

/*  STEP 4. - Initialise branch */
            init_test_branch( bnew, n, CM, GM, GA);
            fscanf(input,"%lf %lf %lf", &(bnew->x[0]), &(bnew->y[0]), &(bnew->z[0]));
            fscanf(input,"%lf %lf %lf", &(bnew->x[1]), &(bnew->y[1]), &(bnew->z[1]));
            fscanf(input,"%lf %lf %lf %lf", &(bnew->p_len), &radius, &(bnew->ga), &(bnew->gm));
            dx = bnew->p_len / ((double) n-1 );
            for ( j=0; j<n; j++) {
                bnew->pl[j] = dx*((double) j );
                bnew->d[j] = 2.0*radius;
            }
            bnew->nobs = n;
            bold = bnew;
        } else if ( strcmp(temp, "Marker") == 0 || strcmp(temp, "marker") == 0 ) {

/*  STEP 5. - Initialise marker */

            printf("Found and initialised a branch contact\n");
            newcon = (contact *) malloc( sizeof(contact) );
            newcon->next = NULL;
            newcon->sd = NULL;
            newcon->cp = NULL;
            fscanf(input,"%lf %lf", &newcon->pl, &newcon->amp );
            bnew->conlist = newcon;
        }
    }
    fclose(input);

/*  STEP 6. - Count dendritic branches at soma */
    bold = first_branch;
    n = 0;
    while ( bold ) {
        bnew = first_branch;
        do {
            tmp = pow(bold->x[0]-bnew->x[1],2)+
                  pow(bold->y[0]-bnew->y[1],2)+
                  pow(bold->z[0]-bnew->z[1],2);
            connected = ( tmp < 0.01 );
            bnew = bnew->child;
        } while ( bnew && !connected );
        if ( !connected ) n++;
        bold = bold->child;
    }
    cell->ndend = n;
    printf("\n\nTree contains %d individual dendrite(s) ...\n", n);

/* STEP 7. - Identify somal dendrites but extract nothing */
    cell->dendlist = (dendrite *) malloc( (cell->ndend)*sizeof(dendrite) );
    bold = first_branch;
    n = 0;
    while ( n < cell->ndend ) {
        bnew = first_branch;
        do {
            bnew = bnew->child;
        } while ( bnew );
        cell->dendlist[n].root = bold;
        n++;
        bold = bold->child;
    }

/*  STEP 8. - Extract root of each dendrite from dendrite list */
    for ( k=0 ; k<cell->ndend ; k++ ) {
        bold = cell->dendlist[k].root;
        remove_branch( &first_branch, bold);
    }

/*  STEP 9. - Build each test dendrite from its root branch */
    for ( k=0 ; k<cell->ndend ; k++ ) {
        build_test_dendrite( &first_branch, cell->dendlist[k].root);
    }
    if ( first_branch != NULL ) printf("\nWarning: Unconnected branch segments still exist\n");
    return cell;
}

/**************************************************************

        Function to build a test dendrite from its root

***************************************************************/
void build_test_dendrite( branch **head, branch *root)
{
    double tmp;
    branch *bnow, *bnext, *btmp;

    bnow = *head;
    while ( bnow != NULL ) {

/* Store bnow's child in case it's corrupted */
        bnext = bnow->child;

/* Search if proximal end of bnow is connected to distal end of root */
        tmp = pow(bnow->x[0]-root->x[1],2)+pow(bnow->y[0]-root->y[1],2)+
              pow(bnow->z[0]-root->z[1],2);
       if ( tmp <= 0.01 ) {

/* Take bnow out of the branch list */
            remove_branch( head, bnow);

/* Connect bnow to the root as the child or a peer of the child.
   Initialise childs' children and peers to NULL as default */
            bnow->child = NULL;
            bnow->peer = NULL;
            bnow->parent = root;

/* Inform root about its child if it's the first child, or add
   new child to first child's peer list */
            if ( root->child != NULL ) {
                btmp = root->child;
                while ( btmp->peer != NULL ) btmp = btmp->peer;
                btmp->peer = bnow;
            } else {
                root->child = bnow;
            }
        }

/* Initialise bnow to next branch in list */
        bnow = bnext;
    }

/* Iterate through remaining tree */
    if ( root->child ) build_test_dendrite( head, root->child);
    if ( root->peer ) build_test_dendrite( head, root->peer);
    return;
}



 /***************************************************************
            Function To Factorise A Sparse Matrix
 ***************************************************************/
void LU_Factor(sparse_mat *m, int *start)
{
    double tmp, sum;
    int i, j, k, r, n, cl, cu, col, row;

/*  Step 1. - Identify matrix dimension */
    n = m->n;

/*  Step 2. - Fill column vectors for triangular matrices */
    if ( *start ) {
        cl = cu = 0;
        for ( i=k=0 ; i<n ; i++ ) {
            m->l->start_row[i] = cl;
            m->u->start_row[i] = cu;
            while ( m->col[k] < i ) m->l->col[cl++] = m->col[k++];
            m->l->col[cl++] = m->col[k];
            m->u->col[cu++] = m->col[k++];
            while ( k < m->start_row[i+1] ) m->u->col[cu++] = m->col[k++];
        }
        m->l->start_row[n] = cl;
        m->u->start_row[n] = cu;
        *start = 0;
    }

/*  Step 3. - Fill remaining entries of L and U row by row */
    cl = cu = 0;
    for ( i=0 ; i<n ; i++ ) {
        for ( k=m->start_row[i] ; k < m->start_row[i+1] ; k++ ) {
            if ( m->col[k] < i ) {
                sum = m->a[k];
                for ( j=m->l->start_row[i] ; j<cl ; j++ ) {
                    col = m->l->col[j];
                    row = m->u->start_row[col];
                    while ( m->u->col[row] < m->col[k] && row < m->u->start_row[col+1] ) row++;
                    if ( m->u->col[row] == m->col[k] ) sum -= (m->l->a[j])*(m->u->a[row]);
                }
                row = m->u->start_row[m->l->col[cl]];
                while ( m->u->col[row] < m->col[k] ) row++;
                m->l->a[cl++] = sum/(m->u->a[row]);
            } else if ( m->col[k] == i ) {
                m->l->a[cl++] = 1.0;
                sum = m->a[k];
                for ( j=m->l->start_row[i] ; j<m->l->start_row[i+1]-1; j++ ) {
                    col = m->l->col[j];
                    row = m->u->start_row[col];
                    while ( m->u->col[row] < i && row < m->u->start_row[col+1] ) row++;
                    if ( m->u->col[row] == i ) sum -= (m->l->a[j])*(m->u->a[row]);
                }
                m->u->a[cu++] = sum;
            } else {
                sum = m->a[k];
                for ( j=m->l->start_row[i] ; j<m->l->start_row[i+1]-1; j++ ) {
                    col = m->l->col[j];
                    row = m->u->start_row[col];
                    while ( m->u->col[row] < m->col[k] && row < m->u->start_row[col+1] ) row++;
                    if ( m->u->col[row] == m->col[k] ) sum -= (m->l->a[j])*(m->u->a[row]);
                }
                m->u->a[cu++] = sum;
            }
        }
    }
    return;
}


/*******************************************************************

          Function to Solve the matrix problem

*******************************************************************/

void LU_Solve(sparse_mat *m, double *x, double *b )
{
    int i,j;
    double *z;

    z = (double *) malloc( (m->n)*sizeof(double) );

    for ( i=0 ; i<m->n ; i++ ) {
        z[i] = b[i];
        for (j=m->l->start_row[i];j<m->l->start_row[i+1]-1;j++)
        { z[i] -= (m->l->a[j])*(z[(m->l->col[j])]); }
        z[i] /= m->l->a[m->l->start_row[i+1]-1];
    }

    for ( i = (m->n) - 1 ; i>=0 ; i-- ) {
        x[i] = z[i];
        for (j=m->u->start_row[i]+1;j<m->u->start_row[i+1];j++)
        { x[i] -= (m->u->a[j])*(x[m->u->col[j]]); }
        x[i] /= m->u->a[m->u->start_row[i]];
    }
    free(z);
    return;
}


 /********************************************************************
                    ALPHA for ACTIVATION OF SODIUM
  *******************************************************************/
double alfa_m( double volt )
{
    double tmp;
    static double fac;
    static int start=1;

    if ( start ) {
        fac = pow(3.0,0.1*CELSIUS-0.63);
        start = !start;
    }
    tmp = -0.1*(volt+35.0);
    if ( fabs(tmp)<0.001 ) {
        tmp = 1.0/(((tmp/24.0+1.0/6.0)*tmp+0.5)*tmp+1.0);
    } else {
        tmp = tmp/(exp(tmp)-1.0);
    }
    return tmp*fac;
}



 /********************************************************************
                    BETA for ACTIVATION OF SODIUM
  ********************************************************************/
double beta_m( double volt )
{
    double tmp;
    static double fac;
    static int start=1;

    if ( start ) {
        fac = pow(3.0,0.1*CELSIUS-0.63);
        start = !start;
    }
    tmp = (volt+60.0)/18.0;
    return 4.0*fac*exp(-tmp);
}



 /********************************************************************
                    ALPHA for INACTIVATION OF SODIUM
  ********************************************************************/
double alfa_h( double volt )
{
    double tmp;
    static double fac;
    static int start=1;

    if ( start ) {
        fac = pow(3.0,0.1*CELSIUS-0.63);
        start = !start;
    }
    tmp = 0.05*(volt+60.0);
    return 0.07*fac*exp(-tmp);
}



 /********************************************************************
                    BETA for INACTIVATION OF SODIUM
  ********************************************************************/
double beta_h( double volt )
{
    double tmp;
    static double fac;
    static int start=1;

    if ( start ) {
        fac = pow(3.0,0.1*CELSIUS-0.63);
        start = !start;
    }
    tmp = -0.1*(volt+30.0);
    return fac/(exp(tmp)+1.0);
}



 /********************************************************************
                    ALPHA for ACTIVATION OF POTASSIUM
  ********************************************************************/
double alfa_n( double volt )
{
    double tmp;
    static double fac;
    static int start=1;

    if ( start ) {
        fac = pow(3.0,0.1*CELSIUS-0.63);
        start = !start;
    }
    tmp = -0.1*(volt+50.0);
    if ( fabs(tmp)<0.001 ) {
        tmp = 0.1/(((tmp/24.0+1.0/6.0)*tmp+0.5)*tmp+1.0);
    } else {
        tmp = 0.1*tmp/(exp(tmp)-1.0);
    }
    return tmp*fac;
}



 /********************************************************************
                    BETA for ACTIVATION OF POTASSIUM
  ********************************************************************/
double beta_n( double volt )
{
    double tmp;
    static double fac;
    static int start=1;

    if ( start ) {
        fac = pow(3.0,0.1*CELSIUS-0.63);
        start = !start;
    }
    tmp = 0.0125*(volt+60.0);
    return 0.125*fac*exp(-tmp);
}

 /*********************************************************************
                    Function to count number of branches
 *********************************************************************/

void count_branch ( branch *b, int *nb )
{
    (*nb)++;
    if ( b->child != NULL ) count_branch ( b->child, nb );
    if ( b->peer != NULL ) count_branch ( b->peer, nb );
    return;
}

 /*********************************************************************
            Function to count total dendritic length
 *********************************************************************/

void count_dendritic_length ( branch *b, double *length )
{
    (*length) += b->p_len;
    if ( b->child != NULL ) count_dendritic_length ( b->child, length );
    if ( b->peer != NULL ) count_dendritic_length ( b->peer, length );
    return;
}

 /********************************************************************
            Function to re-compute physical lengths
 ********************************************************************/

void phys_lengths( branch *b, double *h )
{
    int i,n,k;
    double tmp, hnow, *dtmp, *ptmp, *ptr;


    n = ((int) ceil((b->p_len)/(*h))) + 1;
    hnow = b->p_len/((double) n-1);
    dtmp = (double *) malloc( n*sizeof(double));
    ptmp = (double *) malloc( n*sizeof(double));
    for ( i=0; i<(n-1); i++ ) ptmp[i] = ((double) i)*hnow;
    ptmp[n-1] = b->p_len;

    dtmp[0] = b->d[0];
    for ( i=1 ; i<(n-1) ; i++ ) {
        k = 0;
        while ( ptmp[i] > b->pl[k] ) k++;
        dtmp[i] = b->d[k-1] + (b->d[k] - b->d[k-1])*(ptmp[i] - b->pl[k-1])/(b->pl[k]-b->pl[k-1]);
    }
    dtmp[n-1] = b->d[b->nobs-1];

    ptr = b->d;
    b->d = dtmp;
    free(ptr);
    ptr = b->pl;
    b->pl = ptmp;
    free(ptr);

    b->nobs = n;

    if (b->child != NULL) phys_lengths( b->child, h);
    if (b->peer != NULL) phys_lengths( b->peer, h);
}


 /**********************************************************************
              Function returns primitive uniform random number.
  **********************************************************************/
double ran(unsigned long int *ix,
           unsigned long int *iy,
           unsigned long int *iz)
{
    double tmp;

/*  1st item of modular arithmetic  */
    *ix = (171*(*ix))%30269;
/*  2nd item of modular arithmetic  */
    *iy = (172*(*iy))%30307;
/*  3rd item of modular arithmetic  */
    *iz = (170*(*iz))%30323;
/*  Generate random number in (0,1)  */
    tmp = ((double) (*ix))/30269.0+((double) (*iy))/30307.0
          +((double) (*iz))/30323.0;
    return fmod(tmp,1.0);
}

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

        int xl;                     /* Left hand node */
        int xr;                     /* Right hand node */
        double frac;                /* Fraction of input to left hand node */

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
int count_terminal_branches( branch *, branch *);
void enumerate_nodes( branch *, int *);
void cgs(int, sparse_mat *, double *, double *, double *, double);
void mat_vec_mult( sparse_mat *, double *, double *);
void mat_malloc( sparse_mat *, int, int);
void mat_free( sparse_mat *);
void Generate_Dendrite(branch *, int *);
void init_test_branch( branch *, int , double , double , double );
void build_test_dendrite( branch **, branch *);
void LU_Factor(sparse_mat *, int *);
void LU_Solve( sparse_mat *, double *, double *);
void Input_Current( branch *);
void Assign_Current( branch *, double *, double );
void count_branch ( branch *, int * );
void count_dendritic_length ( branch *, double *);
void phys_lengths( branch *, double * );

/* Global definitions */
#define             CS    1.0
#define             GS    0.091
#define             GA    14.286
#define             CM    1.0
#define             GM    0.091
#define         OUTPUT    "contactinfo.dat"
#define           TEND    10.0
#define             NT    1000
#define             DT    1.0
#define         NNODES    100

/* Parameters for exact solution */
#define            AMP    1.0e-3
#define            SIN    0.0e-3
#define              T    10.0
#define              M    1000

/* Parameters for TestCell1.d3 (unbranched dendrite) */
// #define             RS    0.002
// #define             RD    0.001
// #define              L    0.100
// #define              X    0.055

/* Parameters for TestCell3.d3 (Y-junction) */
// #define             RS    0.002
// #define             RD    0.000354487542
// #define              L    0.0500427736
// #define              X    0.039163520

/* Parameters for TestCell2.d3 (Y-junction) */
#define             RS    0.002
#define             RD    0.000354487542
#define              L    0.0667236981
#define              X    0.029735419

/* Parameters for TestCell4.d3 (large branched dendrite) */
// #define             RS    0.002
// #define             RD    0.000648741708
// #define              L    0.2256605
// #define              X    0.0135280571

/* Global Variables */
sparse_mat lhs, rhs;

int main( int argc, char **argv )
{
    int k, id, start, nodes, nc, i, in, FirstNode;
    int counter, nb;
    double *v, *x, max, *eta, as, gama, *chi, xold, xnew,
           frac, arg, sum, tmp, vs, pi, dt, tnow, tout,
           length, h;
    neuron *n;
    extern sparse_mat lhs, rhs;
    branch *bnow;
    FILE *fp;

/*  Compute ancillary variables */
    pi = 4.0*atan(1.0);
    frac = X/L;

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
    } else if ( strcmp(argv[1], "test") == 0  || strcmp(argv[1], "Test") == 0 ) {
        n = Load_Test_Neuron( argv[2] );
        nb = 0;
        length = 0.0;

        for(i=0; i<n->ndend; i++) count_branch( n->dendlist[i].root, &nb );
        for(i=0; i<n->ndend; i++) count_dendritic_length( n->dendlist[i].root, &length );

        h = length/((double) NNODES - nb);
        for(i=0; i<n->ndend; i++) phys_lengths( n->dendlist[i].root, &h );

        if ( !n ) {
            printf("\n Failed to find sampled neuron\n");
            return(1);
        }
        as = 4.0*pi*RS*RS;
        gama = 0.5*as/(pi*RD*L);
        eta = (double *) malloc( (M+1)*sizeof(double) );
        chi = (double *) malloc( (M+1)*sizeof(double) );
        eta[0] = GM/CM;
        chi[0] = 0.5*(AMP+SIN)/(pi*CM*RD*L*(1.0+gama));
        chi[0] /= eta[0];
        for ( k=1 ; k<=M ; k++ ) {
            xnew = arg = pi*((double) k );
            do {
                xold = xnew;
                xnew = arg-atan(gama*xold);
            } while ( fabs(xold-xnew) > 5.e-7 );
            tmp = cos(xnew);
            eta[k] = (GM+0.5*RD*GA*xnew*xnew/(L*L))/CM;
            chi[k] = tmp*(AMP*cos(xnew*(1.0-frac))+SIN*tmp)/
                         (pi*CM*RD*L*(1.0+gama*tmp*tmp));
            chi[k] /= eta[k];
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
    getchar( );

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

    for( i=0; i<n->ndend ; i++ ) Input_Current(n->dendlist[i].root);

    fp = fopen("out.res", "w");
    fclose(fp);

    dt = 1.0/((double) NT);
    for ( i=0 ; i<3*nodes-2 ; i++ ) {
        rhs.a[i] *= GA;
        rhs.a[i] += GM*lhs.a[i];
        lhs.a[i] *= CM;
        rhs.a[i] *= 0.5*dt;
        lhs.a[i] += rhs.a[i];
        rhs.a[i] = lhs.a[i] - 2.0*rhs.a[i];
    }

/*  Add capacitive term of soma */

    lhs.a[3*nodes-3] += (4.0*as/pi)*(CS+0.5*GS*dt);
    rhs.a[3*nodes-3] += (4.0*as/pi)*(CS-0.5*GS*dt);

    v = (double *) malloc( (nodes)*sizeof(double) );
    x = (double *) malloc( (nodes)*sizeof(double) );
    for ( i=0 ; i<nodes ; i++ ) v[i] = x[i] = 0.0;
    start = 1;
    LU_Factor(&lhs, &start);

/*  Initialise temporal integration */

    tnow = 0.0;
    tout = DT;
    while ( tnow < TEND ) {
        tnow += dt;
        mat_vec_mult(&rhs,v,x);
        x[nodes-1] -= 4.0*dt*SIN/pi;
        tmp = 4.0*dt/pi;
        for ( i=0 ; i<nodes ; i++ ) v[i] = x[i];
        if ( tnow < T ) {
            for ( i=0 ; i<n->ndend ; i++ ) Assign_Current(n->dendlist[i].root, x, tmp);
        }
        LU_Solve( &lhs, v, x );

        if ( tnow > tout && tnow <= T ) {
            printf("\rReached time %5.1lf ms\t", tout);
            for ( vs=0.0,i=M ; i>=0 ; i-- ) {
                arg = tnow*eta[i];
                if ( arg > 20.0 ) {
                    tmp = 1.0;
                } else {
                    tmp = (1.0-exp(-arg));
                }
                vs -= tmp*chi[i];
            }
            printf("\nNumerical Voltage %12.6lf mV\n",v[nodes-1]);
            printf("Exact Voltage %12.6lf mV\n", vs);
            tout += DT;
            fp = fopen("out.res", "a");
            fprintf(fp,"Error %12.6lf% \n",fabs((v[nodes-1]-vs)/(0.01*vs)));
            fclose(fp);
        }
        if ( tnow > tout && tnow > T ) {
            printf("\rReached time %5.1lf ms\t", tout);
            for ( vs=0.0,i=M ; i>=0 ; i-- ) {
                arg = (tnow-T)*eta[i];
                if ( arg < 20.0 ) {
                    tmp = exp(-arg);
                    arg = tnow*eta[i];
                    if ( arg < 20.0 ) tmp -= exp(-arg);
                    vs -= tmp*chi[i];
                }
            }
            printf("\nNumerical Voltage %12.6lf mV\n",v[nodes-1]);
            printf("Exact Voltage %12.6lf mV\n",vs);
            tout += DT;
            fp = fopen("out.res", "a");
            fprintf(fp,"Error %12.6lf% \n",fabs((v[nodes-1]-vs)/(0.01*vs)));
            fclose(fp);
        }
    }

    Destroy_Sampled_Neuron( n );
    return(0);
}

 /******************************************************
        Function to constuct sparse matrices
  ******************************************************/

void Generate_Dendrite( branch *bnow, int *counter)
{
    int i, k;
    extern sparse_mat lhs, rhs;
    branch *btmp;
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
    s->gs = GS;
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

/****************************************************************

        Function to count number of terminal branches

****************************************************************/

int count_terminal_branches( branch *bstart, branch *bnow)
{
    static int n;

    if ( bstart == bnow ) n = 0;
    if ( bnow ) {
        if ( bnow->child ) count_terminal_branches(bstart, bnow->child);
        if ( bnow->peer ) count_terminal_branches(bstart, bnow->peer);
        if ( !bnow->child ) n++;
    }
    return n;
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

/*  STEP 4. - Initialise branch */
            init_test_branch( bnew, n, CM, GM, GA);
            fscanf(input,"%lf %lf %lf", &(bnew->x[0]), &(bnew->y[0]), &(bnew->z[0]));
            fscanf(input,"%lf %lf %lf", &(bnew->x[1]), &(bnew->y[1]), &(bnew->z[1]));
            fscanf(input,"%lf %lf", &(bnew->p_len), &radius);
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
            tmp = pow(bold->x[0]-bnew->x[1],2)+
                  pow(bold->y[0]-bnew->y[1],2)+
                  pow(bold->z[0]-bnew->z[1],2);
            connected = ( tmp < 0.01 );
            bnew = bnew->child;
        } while ( bnew && !connected );
        if ( !connected ) {
            cell->dendlist[n].root = bold;
            n++;
        }
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



 /*************************************************************
            Function to input current to dendrite
 **************************************************************/
void Input_Current( branch *bnow )
{
    int k;
    double tmp;

    if ( bnow->child != NULL ) Input_Current( bnow->child );
    if ( bnow->peer != NULL ) Input_Current( bnow->peer );

    if ( bnow->conlist != NULL ) {
        if ( bnow->conlist->pl <= bnow->pl[1] ) {
            bnow->conlist->xl = bnow->junct;
            bnow->conlist->xr = bnow->first;
            bnow->conlist->frac = (1.0-(bnow->conlist->pl/bnow->pl[1]));
        } else {
            k = 1;
            while ( bnow->conlist->pl > bnow->pl[k] ) k++;
            bnow->conlist->xl = bnow->first-k+2;
            bnow->conlist->xr = bnow->first-k+1;
            bnow->conlist->frac = (bnow->pl[k]-bnow->conlist->pl)/(bnow->pl[k]-bnow->pl[k-1]);
        }
    }
    return;
}

 /************************************************************
                 Function to assign current
 **************************************************************/
void Assign_Current(branch *bnow, double *x, double fac )
{
    if (bnow->child != NULL ) Assign_Current(bnow->child, x, fac );
    if (bnow->peer != NULL ) Assign_Current(bnow->peer, x, fac );

    if ( bnow->conlist != NULL ) {
        x[bnow->conlist->xl] -= fac*(bnow->conlist->frac)*(bnow->conlist->amp);
        x[bnow->conlist->xr] -= fac*(1.0-(bnow->conlist->frac))*(bnow->conlist->amp);
    }
    return;
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

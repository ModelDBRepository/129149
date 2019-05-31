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
        int id;                     /* Identifies contact type */
        double xp;                  /* Location of contact */

        int np;                     /* Proximal neighbour */
        int nd;                     /* Distal neighbour */
        double wp;                  /* Fraction of input carried by proximal neighbour */
        double wd;                  /* Fraction of input carried by distal neighbour */

        double wpp;                 /* Weight of proximal potential in proximal current */
        double wdp;                 /* Weight of distal potential in proximal current */
        double wpd;                 /* Weight of proximal potential in distal current */
        double wdd;                 /* Weight of distal potential in distal current */

        int npp;                    /* Sparse entry occupied by wpp */
        int ndp;                    /* Sparse entry occupied by wdp */
        int npd;                    /* Sparse entry occupied by wpd */
        int ndd;                    /* Sparse entry occupied by wdd */

        double gold;                /* Previous conductance */
        double gnew;                /* Present conductance */

/*  Properties of synaptic conductance profile */
        cond *SynCond;              /* Address of synaptic conductance profile */
        int LocalSynTime;           /* Time to activation of synapse/time since activation */
        int MaxLocalSynTime;        /* Time of inactivation of synapse */
        double vsyn;                /* Reversal potential of synaptic species */
        struct synapse_t *next;     /* Address of next contact */
} synapse;


typedef struct branch_t
{
/*  Connectivity of branch */
        struct branch_t *parent;    /* Address of parent branch */
        struct branch_t *child;     /* Address of child branch */
        struct branch_t *peer;      /* Addresss of peer branch */

/*  Physical properties of branch */
        int id;                     /* Branch identifier */
        int nc;                     /* Number of compartments specifying branch */
        double xl;                  /* X-coordinate of lefthand endpoint */
        double yl;                  /* Y-coordinate of lefthand endpoint */
        double zl;                  /* Z-coordinate of lefthand endpoint */
        double xr;                  /* X-coordinate of righthand endpoint */
        double yr;                  /* Y-coordinate of righthand endpoint */
        double zr;                  /* Z-coordinate of righthand endpoint */
        double diam;                /* Branch diameter (cm) */
        double plen;                /* Branch length (cm) */
        double hseg;                /* Dendritic segment length (cm) */

/*  Node information for spatial representation */
        int nodes;                  /* Total number nodes spanning branch */
        int junct;                  /* Junction node of the branch */
        int first;                  /* Internal node connected to junction */

/*  Contact information */
        synapse *synlist;           /* Branch synapse */
} branch;


typedef struct dendrite_t
{
        branch *root;               /* Pointer to root branch of dendrite */
        double plen;                /* Length of dendrite */
} dendrite;


typedef struct neuron_t
{
        int ndend;                  /* Number of dendrites */
        dendrite *dendlist;         /* Pointer to an array of dendrites */
} neuron;


/* Function type declarations */
cond    *ConductanceProfile( double,  double, double, double );

int     Count_Branches( branch *, branch *),
        Count_Synapses( branch *, branch *);

double  branch_length( branch *, branch *),
        ran(unsigned int *, unsigned int *, unsigned int *),
        alfa_h( double ),
        alfa_m( double ),
        alfa_n( double ),
        beta_h( double ),
        beta_m( double ),
        beta_n( double );

void    Build_Test_Dendrite( branch **, branch *),
        Remove_Branch( branch **, branch *),
        Assign_Branch_Nodes( branch *, double *),
        Enumerate_Nodes( branch *, int *),
        Generate_Dendrite(branch *, int *),
        Initialise_Synapses( branch *),
        Update_Synapses( branch *),
        Matrix_Vector_Multiply( SparseMatrix *, double *, double *),
        Matrix_Malloc( SparseMatrix *, int, int),
        Matrix_Free( SparseMatrix *),
        cgs( int, SparseMatrix *, double *, double *, double *, double);

/* Global definitions */
#define             CS    1.0
#define             GS    0.091
#define             GA    14.286
#define             CM    1.0
#define             GM    0.091
#define         OUTPUT    "K2Spike500.dat"
#define           TEND    11
#define            TAU    0.5
#define            CGS    1.0e-18   /* Tolerance used in CGS algorithm */
#define           GMAX    1.0e-5
#define           RATE    30.0
#define           VSYN    115.0
#define             NT    1000
#define          NODES    500
#define          NSEED    2         /* Seed for random number generator */
#define        CELSIUS    18.5      /* Celsius temperature of neuron */

/* Parameters for exact solution */
#define           NCON    1000        /* Number of contacts */
#define             RS    0.005

/* Global Variables */
SparseMatrix lhs, rhs;
double pi, *SynCurrent;
unsigned int ix, iy, iz;

int main( int argc, char **argv )
{
    extern unsigned int ix, iy, iz;
    extern SparseMatrix lhs, rhs;
    extern double pi, *SynCurrent;
    int k, j, id, nn, nodes, n, nc, i, in, nstep, maxstep, FirstNode,
        NumberOfSynapses, counter, nb, connected, spk, nspk;
    double *v, *x, max, AreaOfSoma, sum, tmp, vs, dt, len, h, sc,
           gs, interval, dx, CellLength, LocusContact, mval,
           nval, hval, aval, bval, vm0, vm1, vm2, tnow;
    double *StoredLHS, *StoredRHS;
    double v_na=115.0, v_k=-12.0, v_l, g_na=120.0, g_k=36.0, g_l=0.3;
    void srand( unsigned int);
    neuron *cell;
    cond *SynCond;
    synapse *newsyn, *syn;
    branch *bo, *bn, *FirstBranch;
    char word[20];
    FILE *fp;

/*  Initialise random number generator */
    fp = fopen(OUTPUT,"w");
    fclose(fp);
    nspk = spk = 0;
    dt = 1.0/((double) NT);
    srand( ((unsigned int) NSEED) );
    ix = rand( );
    iy = rand( );
    iz = rand( );
    SynCond = ConductanceProfile( dt, TAU, 5.e-7, GMAX );

/*  Load Test Neuron */
    maxstep =  1000*NT*TEND;
    pi = 4.0*atan(1.0);
    AreaOfSoma = 4.0*pi*RS*RS;
    if ( argc != 2 ) {
        printf("\n Invoke program with load <input>\n");
        return 1;
    } else {
        printf("\nOpening file %s\n",argv[1]);
        if ( (fp=fopen(argv[1],"r")) == NULL ) {
            printf("\n Test Neuron file not found");
            return 1;
        }
    }

/*  Get branch data */
    bo = NULL;
    while  ( fscanf(fp,"%s",word) != EOF ) {
        if ( strcmp(word,"Branch") == 0 || strcmp(word,"branch") == 0 ) {
            bn = (branch *) malloc( sizeof(branch) );
            fscanf(fp,"%d", &(bn->id) );
            bn->peer = NULL;
            bn->child = NULL;
            bn->synlist = NULL;
            if ( bo != NULL) {
                bo->child = bn;
            } else {
                FirstBranch = bn;
            }
            bn->parent = bo;
            fscanf(fp,"%lf %lf %lf", &(bn->xl), &(bn->yl), &(bn->zl) );
            fscanf(fp,"%lf %lf %lf", &(bn->xr), &(bn->yr), &(bn->zr) );
            fscanf(fp,"%lf %lf", &(bn->plen), &(bn->diam) );
            bo = bn;
        } else {
            printf("Unrecognised dendritic feature\n");
            return 0;
        }
    }
    fclose(fp);

/*  Compute total length of dendrite */
    CellLength = 0.0;
    bn = FirstBranch;
    while ( bn ) {
       CellLength += bn->plen;
       bn = bn->child;
    }

/*  STEP 1. - Randomly place NCON synapses on branches */
    for ( k=0 ; k<NCON ; k++ ) {
        LocusContact = CellLength*ran( &ix, &iy, &iz);
        bn = FirstBranch;
        len = bn->plen;
        while ( LocusContact > len ) {
            bn = bn->child;
            len += bn->plen;
        }
        newsyn = (synapse *) malloc( sizeof(synapse) );
        newsyn->next = NULL;
        newsyn->SynCond = SynCond;
        newsyn->MaxLocalSynTime = SynCond->n;
        newsyn->vsyn = VSYN;
        newsyn->xp = LocusContact-(len-bn->plen);
        syn = bn->synlist;
        if ( syn ) {
            while ( syn->next ) syn = syn->next;
            syn->next = newsyn;
        } else {
            bn->synlist = newsyn;
        }
    }

/*  STEP 2. - Count root branches */
    bo = FirstBranch;
    n = 0;
    while ( bo ) {
        bn = FirstBranch;
        do {
            tmp = pow(bo->xl-bn->xr,2)+pow(bo->yl-bn->yr,2)+
                  pow(bo->zl-bn->zr,2);
            connected = ( tmp < 0.01 );
            bn = bn->child;
        } while ( bn && !connected );
        if ( !connected ) n++;
        bo = bo->child;
    }

/*  STEP 3. - Identify somal dendrites but extract nothing */
    printf("\nTree contains %d individual dendrite(s) ...\n", n);
    cell = (neuron *) malloc( sizeof(neuron) );
    cell->ndend = n;
    cell->dendlist = (dendrite *) malloc( n*sizeof(dendrite) );
    bo = FirstBranch;
    n = 0;
    while ( n < cell->ndend ) {
        bn = FirstBranch;
        do {
            tmp = pow(bo->xl-bn->xr,2)+pow(bo->yl-bn->yr,2)+
                  pow(bo->zl-bn->zr,2);
            connected = ( tmp < 0.01 );
            bn = bn->child;
        } while ( bn );
        if ( !connected ) cell->dendlist[n++].root = bo;
        bo = bo->child;
    }

/*  STEP 4. - Extract root of each dendrite from dendrite list */
    for ( k=0 ; k<cell->ndend ; k++ ) {
        bo = cell->dendlist[k].root;
        Remove_Branch( &FirstBranch, bo);
    }

/*  STEP 5. - Build each test dendrite from its root branch */
    for ( k=0 ; k<cell->ndend ; k++ ) {
        Build_Test_Dendrite( &FirstBranch, cell->dendlist[k].root );
    }
    if ( FirstBranch != NULL ) printf("\nWarning: Unconnected branch segments still exist\n");

/*  STEP 6. - Count number of synapses on Cell */
    NumberOfSynapses = 0;
    for ( k=0 ; k<cell->ndend ; k++ ) {
        bn = cell->dendlist[k].root;
        NumberOfSynapses += Count_Synapses( cell->dendlist[k].root, bn);
    }
    printf("\nNumber of Synapses %d", NumberOfSynapses);

/*  STEP 7. - Count number of dendritic branches */
    for ( nb=k=0 ; k<cell->ndend ; k++ ) {
        bn = cell->dendlist[k].root;
        nb += Count_Branches( bn, bn);
    }
    h = CellLength/((double) NODES-nb);
    for ( k=0 ; k<cell->ndend ; k++ ) Assign_Branch_Nodes( cell->dendlist[k].root, &h);

/*  STEP 8. - Enumerate Nodes */
    FirstNode = 0;
    for ( k=0 ; k<cell->ndend ; k++ ) Enumerate_Nodes( cell->dendlist[k].root, &FirstNode );
    for ( k=0 ; k<cell->ndend ; k++ ) cell->dendlist[k].root->junct = FirstNode;
    printf("\nNumber of nodes is %d\n", FirstNode+1);

/*  STEP 9. - Construct Sparse Matrices */
    nodes = FirstNode+1;
    Matrix_Malloc( &lhs, nodes, 3*nodes-2 );
    Matrix_Malloc( &rhs, nodes, 3*nodes-2 );
    lhs.StartRow[0] = rhs.StartRow[0] = 0;
    for ( counter=k=0 ; k<cell->ndend ; k++ ) {
        bn = cell->dendlist[k].root;
        Generate_Dendrite( bn, &counter);
    }
    lhs.n = rhs.n = nodes;

/*  STEP 10. - Do somal node */
    lhs.a[3*nodes-3] = rhs.a[3*nodes-3] = 0.0;
    for ( k=0 ; k<cell->ndend ; k++ ) {
        bn = cell->dendlist[k].root;
        lhs.a[counter] = (bn->diam)*(bn->hseg)/6.0;
        rhs.a[counter] = -0.25*pi*pow(bn->diam,2)/(bn->hseg);
        lhs.col[counter] = rhs.col[counter] = bn->first;
        lhs.a[3*nodes-3] += 2.0*pi*(bn->diam)*(bn->hseg)/6.0;
        rhs.a[3*nodes-3] += 0.25*pi*pow(bn->diam,2)/(bn->hseg);
        counter++;
    }
    lhs.col[counter] = rhs.col[counter] = nodes-1;
    lhs.StartRow[nodes] = rhs.StartRow[nodes] = counter+1;

/*  STEP 12. - Fill in properties of synapses */
    for( k=0 ; k<cell->ndend ; k++ ) Initialise_Synapses(cell->dendlist[k].root);
    for ( k=0 ; k<3*nodes-2 ; k++ ) {
        rhs.a[k] = 0.5*dt*(GA*rhs.a[k]+GM*lhs.a[k]);
        rhs.a[k] = CM*lhs.a[k]-rhs.a[k];
        lhs.a[k] = 2.0*CM*lhs.a[k]-rhs.a[k];
    }
    lhs.a[3*nodes-3] += AreaOfSoma*CS;
    rhs.a[3*nodes-3] += AreaOfSoma*CS;

/*  STEP 13. - Construct and load vectors to hold transient information */
    StoredLHS = (double *) malloc( (3*nodes-2)*sizeof(double) );
    StoredRHS = (double *) malloc( (3*nodes-2)*sizeof(double) );
    for ( k=0 ; k<3*nodes-2 ; k++ ) {
        StoredLHS[k] = lhs.a[k];
        StoredRHS[k] = rhs.a[k];
    }

/*  STEP 14. - Compute somal conductances and the leakage potential */
    g_na *= AreaOfSoma;
    g_k *= AreaOfSoma;
    g_l *= AreaOfSoma;
    hval = alfa_h(0.0)/(alfa_h(0.0)+beta_h(0.0));
    mval = alfa_m(0.0)/(alfa_m(0.0)+beta_m(0.0));
    nval = alfa_n(0.0)/(alfa_n(0.0)+beta_n(0.0));
    v_l = g_na*pow(mval,3)*hval*v_na+g_k*pow(nval,4)*v_k;
    v_l = -v_l/g_l;

/*  STEP 15. - Construct and initialise potentials and currents */
    v = (double *) malloc( (nodes)*sizeof(double) );
    x = (double *) malloc( (nodes)*sizeof(double) );
    SynCurrent = (double *) malloc( (nodes)*sizeof(double) );
    for ( k=0 ; k<nodes ; k++ ) v[k] = 0.0;

/*  Initialise temporal integration and integrate forward */
    nstep = 0;
    while ( nstep < maxstep ) {
        nstep++;

/*  Phase 1. - Update HH channel variables */
        vs = v[nodes-1];
        aval = dt*alfa_h(vs);
        bval = dt*beta_h(vs);
        tmp = 0.5*(aval+bval);
        hval = (aval+(1.0-tmp)*hval)/(1.0+tmp);
        aval = dt*alfa_m(vs);
        bval = dt*beta_m(vs);
        tmp = 0.5*(aval+bval);
        mval = (aval+(1.0-tmp)*mval)/(1.0+tmp);
        aval = dt*alfa_n(vs);
        bval = dt*beta_n(vs);
        tmp = 0.5*(aval+bval);
        nval = (aval+(1.0-tmp)*nval)/(1.0+tmp);

/*  Phase 2. - Compute somal conductance and contribution to current */
        gs = g_l;                       /* Leakage conductance */
        sc = g_l*v_l;                   /* Leakage contribution to somal current */
        tmp = g_na*hval*pow(mval,3);    /* Sodium conductance */
        gs += tmp;
        sc += tmp*v_na;                 /* Sodium contribution to somal current */
        tmp = g_k*pow(nval,4);          /* Potassium conductance */
        gs += tmp;
        sc += tmp*v_k;                  /* Potasium contribution to somal current */

/*  Phase 3. - Zero LHS, RHS and SynCurrent */
        for ( k=0 ; k<3*nodes-2 ; k++ ) lhs.a[k] = rhs.a[k] = 0.0;
        for ( k=0 ; k<nodes ; k++ ) SynCurrent[k] = 0.0;

/*  Phase 4. - Update synaptic conductances and input */
        for ( k=0 ; k<cell->ndend ; k++ ) {
            bn = cell->dendlist[k].root;
            Update_Synapses( bn );
        }

/*  Phase 5. - Complete the construction of LHS and RHS matrices */
        for ( k=0 ; k<3*nodes-2 ; k++ ) {
            lhs.a[k] += StoredLHS[k];
            rhs.a[k] += StoredRHS[k];
        }
        gs *= 0.5*dt;
        lhs.a[3*nodes-3] += gs;
        rhs.a[3*nodes-3] -= gs;

/*  Phase 6. - Step potential forward */
        Matrix_Vector_Multiply( &rhs, v, x);
        x[nodes-1] += sc*dt;
        for ( k=0 ; k<nodes ; k++ ) x[k] += SynCurrent[k];
        cgs( 1, &lhs, x, v, v, CGS);

/*  Phase 7. - Test for spikes */
        if ( nstep == 1 ) {
            vm2 = 0.0;
            vm1 = v[nodes-1];
        } else {
            vm0 = v[nodes-1];
            if ( !spk ) {
                spk = ( vm0 > 50.0 && vm1 > vm2 && vm1 > vm0 );
                if ( nstep >= 1000*NT && spk ) {
                    tnow = dt*((double) nstep);
                    tmp = tnow+0.5*dt*(vm2+3.0*vm0-4.0*vm1)/(vm0-2.0*vm1+vm2);
                    nn = ((int) floor(tmp))-1000*NT;
                    if ( fmod(tmp,1.0)>0.5 ) nn++;
                    fp = fopen(OUTPUT,"a");
                    fprintf(fp,"%d\n",nn);
                    fclose(fp);
                }
            }
            vm2 = vm1;
            vm1 = vm0;
        }

/*  Phase 8. - Reset spike flag */
        if ( vs < 0.0 && spk == 1 ) {
            spk = 0;
            nspk++;
        }
        if ( nstep%(500*NT) == 0 ) {
            tnow = dt*((double) nstep/NT);
            printf("\rReached time %5.1lf ms \t Spikes so far %d", tnow, nspk);
        }
    }
    return 0;
}


 /******************************************************
     Function to build a test dendrite from its root
  ******************************************************/
void Build_Test_Dendrite( branch **head, branch *root)
{
    double tmp;
    branch *bnow, *bnext, *btmp;

    bnow = *head;
    while ( bnow != NULL ) {

/*  Store bnow's child in case it's corrupted */
        bnext = bnow->child;

/*  Decide if proximal end of bnow is connected to distal end of root */
        tmp = pow(bnow->xl-root->xr,2)+
              pow(bnow->yl-root->yr,2)+
              pow(bnow->zl-root->zr,2);
        if ( tmp <= 0.01 ) {

/*  Remove bnow from the branch list */
            Remove_Branch( head, bnow);

/*  Connect bnow to the root as the child or a peer of the child.
    Initialise childs' children and peers to NULL as default */
            bnow->child = NULL;
            bnow->peer = NULL;
            bnow->parent = root;

/*  Inform root about its child if it's the first child, or add
    new child to first child's peer list */
            if ( root->child != NULL ) {
                btmp = root->child;
                while ( btmp->peer != NULL ) btmp = btmp->peer;
                btmp->peer = bnow;
            } else {
                root->child = bnow;
            }
        }

/*  Initialise bnow to next branch in list */
        bnow = bnext;
    }

/* Iterate through remaining tree */
    if ( root->child ) Build_Test_Dendrite( head, root->child);
    if ( root->peer ) Build_Test_Dendrite( head, root->peer);
    return;
}


 /*********************************************************
        Function to remove a branch from a branch list
  *********************************************************/
void Remove_Branch(branch **head, branch *b)
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



 /*********************************************
      Function to count synapses on a branch
  *********************************************/
int Count_Synapses( branch *bstart, branch *bnow)
{
    static int n;
    synapse *syn;

    if ( bstart == bnow ) n = 0;
    if ( bnow != NULL ) {
        if ( bnow->child ) Count_Synapses(bstart, bnow->child);
        if ( bnow->peer ) Count_Synapses(bstart, bnow->peer);
        syn = bnow->synlist;
        while ( syn ) {
            n++;
            syn = syn->next;
        }
    }
    return n;
}


 /**********************************************
      Function to count number of branches
  **********************************************/
int Count_Branches( branch *bstart, branch *bnow)
{
    static int n;

    if ( bstart == bnow ) n = 0;
    if ( bnow != NULL ) {
        if ( bnow->child ) Count_Branches(bstart, bnow->child);
        if ( bnow->peer ) Count_Branches(bstart, bnow->peer);
        n++;
    }
    return n;
}


 /*******************************************************
       Function to enumerate the nodes on a dendrite
  *******************************************************/
void Enumerate_Nodes(branch *bnow, int *FirstNode )
{
    branch *btmp;

    if ( bnow->child ) Enumerate_Nodes( bnow->child, FirstNode );
    if ( bnow->peer ) Enumerate_Nodes( bnow->peer, FirstNode );

    if ( bnow->child ) {
        btmp = bnow->child;
        while ( btmp ) {
            btmp->junct = *FirstNode;
            btmp = btmp->peer;
        }
    }
    *FirstNode += bnow->nc;
    bnow->first = *FirstNode-1;
    return;
}


 /***************************************************
        Function to constuct sparse matrices
  ***************************************************/
void Generate_Dendrite( branch *b, int *counter)
{
    int k, CurrentNode, nc;
    extern double pi;
    extern SparseMatrix lhs, rhs;
    branch *btmp;
    double SumL, SumR;

/* Step 1 - Recurse to the end of the dendrite */
    if ( b->child ) Generate_Dendrite( b->child, counter);
    if ( b->peer ) Generate_Dendrite( b->peer, counter);

/* Step 2 - Build matrix entries for distal node of branch */
    nc = b->nc;
    CurrentNode = (b->first)-(nc-1);
    if ( b->child ) {
        btmp = b->child;
        SumR = SumL = 0.0;
        while ( btmp ) {
            lhs.a[*counter] = pi*(btmp->diam)*(btmp->hseg)/6.0;
            rhs.a[*counter] = -0.25*pi*pow(btmp->diam,2)/(btmp->hseg);
            SumL += 2.0*pi*(btmp->diam)*(btmp->hseg)/6.0;
            SumR += 0.25*pi*pow(btmp->diam,2)/(btmp->hseg);
            lhs.col[*counter] = rhs.col[*counter] = btmp->first;
            (*counter)++;
            btmp = btmp->peer;
        }
        lhs.a[*counter] = SumL+2.0*pi*(b->diam)*(b->hseg)/6.0;
        rhs.a[*counter] = SumR+0.25*pi*pow(b->diam,2)/(b->hseg);
        lhs.col[*counter] = rhs.col[*counter] = CurrentNode;
        (*counter)++;
        lhs.a[*counter] = pi*(b->diam)*(b->hseg)/6.0;
        rhs.a[*counter] = -0.25*pi*pow(b->diam,2)/(b->hseg);
        if ( CurrentNode == b->first ) {
            lhs.col[*counter] = rhs.col[*counter] = b->junct;
        } else {
            lhs.col[*counter] = rhs.col[*counter] = CurrentNode+1;
        }
        (*counter)++;
        lhs.StartRow[CurrentNode+1] = rhs.StartRow[CurrentNode+1] = *counter;
    } else {
        lhs.a[*counter] = 2.0*pi*(b->diam)*(b->hseg)/6.0;
        rhs.a[*counter] = 0.25*pi*pow(b->diam,2)/(b->hseg);
        lhs.col[*counter] = rhs.col[*counter] = CurrentNode;
        (*counter)++;
        lhs.a[*counter] = pi*(b->diam)*(b->hseg)/6.0;
        rhs.a[*counter] = -0.25*pi*pow(b->diam,2)/(b->hseg);
        if ( CurrentNode == b->first ) {
            lhs.col[*counter] = rhs.col[*counter] = b->junct;
        } else {
            lhs.col[*counter] = rhs.col[*counter] = CurrentNode+1;
        }
        (*counter)++;
        lhs.StartRow[CurrentNode+1] = rhs.StartRow[CurrentNode+1] = *counter;
    }

/* Step 3 - Build matrix entries for internal nodes of branch */
    for ( k=nc-1 ; k>0 ; k-- ) {
        CurrentNode++;
        lhs.a[*counter] = pi*(b->diam)*(b->hseg)/6.0;
        rhs.a[*counter] = -0.25*pi*pow(b->diam,2)/(b->hseg);
        lhs.col[*counter] = rhs.col[*counter] = CurrentNode-1;
        (*counter)++;
        lhs.a[*counter] = 4.0*pi*(b->diam)*(b->hseg)/6.0;
        rhs.a[*counter] = 0.5*pi*pow(b->diam,2)/(b->hseg);
        lhs.col[*counter] = rhs.col[*counter] = CurrentNode;
        (*counter)++;
        lhs.a[*counter] = pi*(b->diam)*(b->hseg)/6.0;
        rhs.a[*counter] = -0.25*pi*pow(b->diam,2)/(b->hseg);
        if ( CurrentNode == b->first ) {
            lhs.col[*counter] = rhs.col[*counter] = b->junct;
        } else {
            lhs.col[*counter] = rhs.col[*counter] = CurrentNode+1;
        }
        (*counter)++;
        lhs.StartRow[CurrentNode+1] = rhs.StartRow[CurrentNode+1] = *counter;
    }
    return;
}


 /***********************************************
       Function to assign synaptic weights
  ***********************************************/
void Initialise_Synapses( branch *b )
{
    extern SparseMatrix lhs, rhs;
    int k;
    double rat, interval;
    synapse *syn;

    if ( b->child ) Initialise_Synapses( b->child );
    if ( b->peer ) Initialise_Synapses( b->peer );

    syn = b->synlist;
    while ( syn ) {
        rat = (syn->xp)/(b->hseg);
        k = ((int) floor(rat));
        rat = fmod(rat,1.0);

/*  Phase 1. - Set up matrix nodes */
        if ( k == 0 ) {
            syn->np = b->junct;
            syn->nd = b->first;
        } else {
            syn->np = b->first-k+1;
            syn->nd = b->first-k;
        }

/*  Phase 2. - Set up weights */
        syn->wpp = (1.0-rat)*(1.0-rat);     /* Weight of proximal potential in proximal current */
        syn->wdp = rat*(1.0-rat);           /* Weight of distal potential in proximal current */
        syn->wpd = rat*(1.0-rat);           /* Weight of proximal potential in distal current */
        syn->wdd = rat*rat;                 /* Weight of distal potential in distal current */
        syn->wp = (1.0-rat)*(syn->vsyn);    /* Weight of input at proximal node */
        syn->wd = rat*(syn->vsyn);          /* Weight of input at distal node */

/*  Phase 3. - Location of sparse entries in equation arising from distal neighbour */
        k = lhs.StartRow[syn->nd];
        while ( syn->nd != lhs.col[k] ) k++;
        syn->ndd = k;
        while ( syn->np != lhs.col[k] ) k++;
        syn->npd = k;

/*  Phase 4. - Location of sparse entries in equation arising from proximal neighbour */
        k = lhs.StartRow[syn->np];
        while ( syn->nd != lhs.col[k] ) k++;
        syn->ndp = k;
        while ( syn->np != lhs.col[k] ) k++;
        syn->npp = k;

/*  Phase 5. - Set initial conductances and firing times */
        syn->gold = syn->gnew = 0.0;
        interval = -(((double) NT*1000)/RATE)*log(ran(&ix, &iy, &iz));
        if ( fmod(interval,1.0) <= 0.5 ) {
            syn->LocalSynTime = -((int) floor(interval));
        } else {
            syn->LocalSynTime = -((int) ceil(interval));
        }
        syn = syn->next;
    }
    return;
}


 /***********************************************
       Function to update status of synapses
  ***********************************************/
void Update_Synapses( branch *b )
{
    extern SparseMatrix lhs, rhs;
    extern double *SynCurrent;
    double interval, gold, gnew, tmp;
    synapse *syn;

    if ( b->child ) Update_Synapses( b->child );
    if ( b->peer ) Update_Synapses( b->peer );

    syn = b->synlist;
    while ( syn ) {
        gold = syn->gold = syn->gnew;
        (syn->LocalSynTime)++;
        if ( syn->LocalSynTime < 1 ) {
            gnew = syn->gnew = 0.0;
        } else if ( syn->LocalSynTime < syn->MaxLocalSynTime ) {
            gnew = syn->gnew = syn->SynCond->g[syn->LocalSynTime];
        } else {
            gnew = syn->gnew = 0.0;
            interval = -(((double) NT*1000)/RATE)*log(ran(&ix, &iy, &iz));
            if ( fmod(interval,1.0) <= 0.5 ) {
                syn->LocalSynTime = -((int) floor(interval));
            } else {
                syn->LocalSynTime = -((int) ceil(interval));
            }
        }
        if ( gold != 0.0 || gnew != 0.0 ) {
            tmp = 0.5*(syn->SynCond->dt);
            SynCurrent[syn->nd] += (syn->wd)*tmp*(gold+gnew);
            SynCurrent[syn->np] += (syn->wp)*tmp*(gold+gnew);

            lhs.a[syn->ndd] += (syn->wdd)*tmp*gnew;
            rhs.a[syn->ndd] -= (syn->wdd)*tmp*gold;

            lhs.a[syn->npd] += (syn->wpd)*tmp*gnew;
            rhs.a[syn->npd] -= (syn->wpd)*tmp*gold;

            lhs.a[syn->ndp] += (syn->wdp)*tmp*gnew;
            rhs.a[syn->ndp] -= (syn->wdp)*tmp*gold;

            lhs.a[syn->npp] += (syn->wpp)*tmp*gnew;
            rhs.a[syn->npp] -= (syn->wpp)*tmp*gold;
        }
        syn = syn->next;
    }
    return;
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


 /**************************************************
            Function to assign branch nodes
  **************************************************/
void Assign_Branch_Nodes( branch *b, double *h )
{
    int k;
    double hseg;

    b->nc = ((int) ceil((b->plen)/(*h)));
    b->hseg = (b->plen)/((double) b->nc);

    if ( b->child ) Assign_Branch_Nodes( b->child, h);
    if ( b->peer ) Assign_Branch_Nodes( b->peer, h);
    return;
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
    tmp = -0.1*(volt-25.0);
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
    tmp = volt/18.0;
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
    tmp = 0.05*volt;
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
    tmp = -0.1*(volt-30.0);
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
    tmp = -0.1*(volt-10.0);
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
    tmp = 0.0125*volt;
    return 0.125*fac*exp(-tmp);
}


 /********************************************
        Computes a conductance profile
  ********************************************/
cond *ConductanceProfile( double dt,   /* Integration time step */
                          double tau,  /* Synaptic time constant */
                          double tol,  /* Determines off-threshold for synapse */
                          double gmax  /* Maximum conductance (mS) */ )
{
    int k;
    cond *out;
    double tmp, told, tnew;

    out = (cond *) malloc( sizeof(cond) );
    out->dt = dt;
    out->tau = tau;
    out->tol = tol;
    out->gmax = gmax;

/* Iterate to find duration of pulse */
    tmp = 1.0-log(tol);
    tnew = tmp;
    do {
        told = tnew;
        tnew = tmp+log(told);
    } while ( fabs(tnew-told)>5.e-11 );
    out->n = ((int) ceil(tau*tnew/dt));
    out->g = (double *) malloc( (out->n)*sizeof(double) );
    out->g[0] = 0.0;
    for ( k=1 ; k<(out->n) ; k++ ) {
        tmp = dt*((double) k)/tau;
        out->g[k] = gmax*tmp*exp(1.0-tmp);
    }
    return out;
}


void cgs(int getmem, SparseMatrix *a, double *b, double *x0, double *x, double tol)
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
    Matrix_Vector_Multiply( a, x0, r);
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

 /***************************************************************
    Function to analyse the numerical error of Generalised
    Compartmental Model. A single input is simulated and
    the exact solution determined using the Equivalent Cable
  ***************************************************************/

/* Global definitions */
#define             CS    1.0
#define             GS    0.091
#define             GA    14.286
#define             CM    1.0
#define             GM    0.091
#define         OUTPUT    "KenSpike200.dat"
#define           TEND    11        /* Must be an integer */
#define            TAU    0.5
#define            CGS    1.0e-26   /* Tolerance used in CGS algorithm */
#define           GMAX    3.0e-5
#define           RATE    30.0
#define           VSYN    115.0
#define             NT    1000
#define          NODES    200
#define          NSEED    2         /* Seed for random number generator */
#define        CELSIUS    18.5      /* Celsius temperature of neuron */

/* Parameters for exact solution */
#define           NCON    100        /* Number of contacts */
#define             RS    0.005

typedef struct SparseMatrix_t
{
        double *a;
        int *col;
        int *StartRow;
        int n;

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
        double vlam;                /* Location of contact */

        int np;                     /* Proximal neighbour */
        int nd;                     /* Distal neighbour */

/*  Properties of synaptic conductance profile */
        cond *SynCond;              /* Address of synaptic conductance profile */
        int SynTime;                /* Time to activation of synapse/time since activation */
        int MaxStep;                /* Time of inactivation of synapse */
        double vsyn;                /* Reversal potential of synaptic species */
        struct synapse_t *NextSyn;  /* Address of next synapse */
} synapse;


typedef struct segment_t
{
/*  Physical properties of segment */
        double rp;                  /* Proximal radius of segment */
        double rd;                  /* Distal radius of segment */
        double hseg;                /* Dendritic segment length (cm) */
        int nsyn;                   /* Number of synapses on segment */

/*  Synaptic conductances */
        double *gold;               /* Old conductance */
        double *gnew;               /* New conductance */

/*  Auxiliary vectors */
        double *vlam;
        double *diag;
        double *vec;
        double *phi;
        double *phi1;
        double *phi2;
        double *phi3;

/*  Solution vector */
        double *soln;

/*  Auxiliary scalars */
        double OldProx1;
        double OldProx2;
        double OldProx3;
        double NewProx1;
        double NewProx2;
        double NewProx3;
        double OldDist1;
        double OldDist2;
        double OldDist3;
        double NewDist1;
        double NewDist2;
        double NewDist3;

/*  Contact information */
        synapse *SynList;           /* List of segment synapses */
} segment;



typedef struct branch_t
{
/*  Connectivity of branch */
        struct branch_t *parent;    /* Address of parent branch */
        struct branch_t *child;     /* Address of child branch */
        struct branch_t *peer;      /* Addresss of peer branch */

/*  Physical properties of branch */
        int id;                     /* Branch identifier */
        int nseg;                   /* Number of segments specifying branch */
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
        segment *SegList;           /* Pointer to an array of segments */
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

int     Count_Synapses( branch *, branch *);

double  branch_length( branch *, branch *),
        ran(unsigned int *, unsigned int *, unsigned int *),
        alfa_h( double ),
        alfa_m( double ),
        alfa_n( double ),
        beta_h( double ),
        beta_m( double ),
        beta_n( double );

void    Build_Test_Dendrite( branch **, branch *),
        Output_Information( branch *, FILE *),
        Remove_Branch( branch **, branch *),
        Enumerate_Nodes( branch *, int *),
        Generate_Dendrite(branch *, int *),
        Initialise_Synapses( branch *),
        Update_Synapses( branch *),
        solve( int, double *, double *, double *),
        Matrix_Vector_Multiply( SparseMatrix *, double *, double *),
        Matrix_Malloc( SparseMatrix *, int, int),
        Matrix_Free( SparseMatrix *),
        Update_Sparse_Matrices( branch *, int *),
        cgs( int *, SparseMatrix *, double *, double *, double);


/* Global Variables */
SparseMatrix lhs, rhs;
double pi, dt, *NewAmp;
unsigned int ix, iy, iz;

int main( int argc, char **argv )
{
    extern unsigned int ix, iy, iz;
    extern SparseMatrix lhs, rhs;
    extern double pi, dt, *NewAmp;
    int k, j, id, nn, nodes, n, nseg, i, in, nstep, maxstep, FirstNode,
        nsyn, NumberOfBranches, counter, connected, spk, nspk, getmem;
    double *v, *x, *OldAmp, max, AreaOfSoma, sum, tmp, vs, len, h, sc,
           gs, interval, dx, CellLength, LocusContact, mval, nval,
           hval, aval, bval, vm0, vm1, vm2, tnow;
    double *StoredLHS, *StoredRHS;
    double v_na=115.0, v_k=-12.0, v_l, g_na=120.0, g_k=36.0, g_l=0.3;
    void srand( unsigned int);
    neuron *cell;
    cond *SynCond;
    synapse *NewSyn, *syn;
    segment *seg;
    branch *bo, *bn, *FirstBranch;
    char word[20];
    FILE *fp, *fp1;

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

/*  STEP 1A. - Determine cell length and number of branches */
    CellLength = 0.0;
    NumberOfBranches = 0;
    bn = FirstBranch;
    while ( bn ) {
        NumberOfBranches++;
        CellLength += bn->plen;
        bn = bn->child;
    }

/*  STEP 1B. - Determine segment lengths and allocate segments */
    h = CellLength/((double) NODES-NumberOfBranches);
    bn = FirstBranch;
    while ( bn ) {
        nseg = bn->nseg = ((int) ceil((bn->plen)/h));
        bn->hseg = (bn->plen)/((double) bn->nseg);
        bn->SegList = (segment *) malloc( nseg*sizeof(segment) );
        for ( k=0 ; k<nseg ; k++ ) {
            (bn->SegList)[k].hseg = bn->hseg;
            (bn->SegList)[k].rp = 0.5*(bn->diam);
            (bn->SegList)[k].rd = 0.5*(bn->diam);
            (bn->SegList)[k].SynList = NULL;
            (bn->SegList)[k].gold = NULL;
            (bn->SegList)[k].gnew = NULL;
            (bn->SegList)[k].soln = NULL;
            (bn->SegList)[k].nsyn = 0;
        }
        bn = bn->child;
    }

/*  STEP 1C. - Randomly place NCON synapses on branches */
    for ( nsyn=0 ; nsyn<NCON ; nsyn++ ) {
        LocusContact = CellLength*ran( &ix, &iy, &iz);
        bn = FirstBranch;
        len = bn->plen;
        while ( LocusContact > len ) {
            bn = bn->child;
            len += bn->plen;
        }
        NewSyn = (synapse *) malloc( sizeof(synapse) );
        NewSyn->NextSyn = NULL;
        NewSyn->SynCond = SynCond;
        NewSyn->MaxStep = SynCond->n;
        NewSyn->vsyn = VSYN;
        LocusContact -= (len-(bn->plen));
        LocusContact /= bn->hseg;
        NewSyn->vlam = fmod(LocusContact,1.0);
        nseg = ((int) floor(LocusContact));
        (bn->SegList)[nseg].nsyn++;
        syn = (bn->SegList)[nseg].SynList;
        if ( syn ) {
            if ( NewSyn->vlam < syn->vlam ) {
                NewSyn->NextSyn = syn;
                (bn->SegList)[nseg].SynList = NewSyn;
            } else {
                while ( syn->NextSyn && NewSyn->vlam > syn->NextSyn->vlam ) syn = syn->NextSyn;
                if ( syn->NextSyn ) NewSyn->NextSyn = syn->NextSyn;
                syn->NextSyn = NewSyn;
            }
        } else {
            (bn->SegList)[nseg].SynList = NewSyn;
        }
    }

/*  STEP 2A. - Count root branches */
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

/*  STEP 2B. - Identify somal dendrites but extract nothing */
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

/*  STEP 2C. - Extract root of each dendrite from dendrite list */
    for ( k=0 ; k<cell->ndend ; k++ ) {
        bo = cell->dendlist[k].root;
        Remove_Branch( &FirstBranch, bo);
    }

/*  STEP 2D. - Build each test dendrite from its root branch */
    for ( k=0 ; k<cell->ndend ; k++ ) {
        Build_Test_Dendrite( &FirstBranch, cell->dendlist[k].root );
    }
    if ( FirstBranch != NULL ) printf("\nWarning: Unconnected branch segments still exist\n");

/*  STEP 2E. - Count number of synapses on Cell */
    for ( nsyn=k=0 ; k<cell->ndend ; k++ ) {
        bn = cell->dendlist[k].root;
        nsyn += Count_Synapses( cell->dendlist[k].root, bn);
    }
    printf("\nNumber of Synapses %d", nsyn);

/*  STEP 3A. - Enumerate Nodes */
    FirstNode = 0;
    for ( k=0 ; k<cell->ndend ; k++ ) Enumerate_Nodes( cell->dendlist[k].root, &FirstNode );
    for ( k=0 ; k<cell->ndend ; k++ ) {
        cell->dendlist[k].root->junct = FirstNode;
    }
    printf("\nNumber of nodes is %d\n", FirstNode+1);
    getchar( );

/*  STEP 3B. - Construct Sparse Matrices */
    nodes = FirstNode+1;
    Matrix_Malloc( &lhs, nodes, 3*nodes-2 );
    Matrix_Malloc( &rhs, nodes, 3*nodes-2 );
    lhs.StartRow[0] = rhs.StartRow[0] = 0;
    for ( counter=k=0 ; k<cell->ndend ; k++ ) {
        bn = cell->dendlist[k].root;
        Generate_Dendrite( bn, &counter);
    }
    lhs.n = rhs.n = nodes;

/*  STEP 3C. - Do somal node */
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

/*  STEP 3D. - Construct Vectors to hold currents */
    OldAmp = (double *) malloc( nodes*sizeof(double) );
    NewAmp = (double *) malloc( nodes*sizeof(double) );
    for ( k=0 ; k<nodes ; k++ ) NewAmp[k] = 0.0;

/*  STEP 3E. - Fill in properties of synapses */
    for( k=0 ; k<cell->ndend ; k++ ) {
        bn = cell->dendlist[k].root;
        Initialise_Synapses(bn);
    }
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

//    fp1 = fopen("KenOutput.dat","w");
//    for( k=0 ; k<3*nodes-2 ; k++ ) {
//        fprintf(fp1,"%12.10lf\n", lhs.a[k]);
//    }
//    fclose(fp1);
//    printf("\nOutput information");
//    getchar( );

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
    for ( k=0 ; k<nodes ; k++ ) v[k] = 0.0;

/*  Initialise temporal integration and integrate forward */
    nstep = 0;
    getmem = 1;
    while ( nstep < maxstep ) {
        nstep++;
//        if ( nstep%1000 == 0 ) printf("\r Msec %d",nstep/1000);

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

/*  Phase 3. - Zero LHS, RHS and synaptic currents */
        for ( k=0 ; k<3*nodes-2 ; k++ ) lhs.a[k] = rhs.a[k] = 0.0;
        for ( k=0 ; k<nodes ; k++ ) {
            OldAmp[k] = NewAmp[k];
            NewAmp[k] = 0.0;
        }

/*  Phase 4. - Update synaptic conductances and input */
        counter = 0;
        for ( k=0 ; k<cell->ndend ; k++ ) {
            bn = cell->dendlist[k].root;
            Update_Synapses( bn );
            bn = cell->dendlist[k].root;
            Update_Sparse_Matrices( bn, &counter);
        }

/*  Phase 4a. - Update soma */
        tmp = 0.5*dt;
        for ( k=0 ; k<cell->ndend ; k++ ) {
            bn = cell->dendlist[k].root;
            seg = &(bn->SegList)[0];
            lhs.a[counter] = tmp*(seg->NewProx2);
            rhs.a[counter] = -tmp*(seg->OldProx2);
            NewAmp[nodes-1] += seg->NewProx3;
            lhs.a[3*nodes-3] += tmp*(seg->NewProx1);
            rhs.a[3*nodes-3] -= tmp*(seg->OldProx1);
            counter++;
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
        tmp = 0.5*dt;
        for ( k=0 ; k<nodes ; k++ ) x[k] -= tmp*(NewAmp[k]+OldAmp[k]);
        cgs( &getmem, &lhs, x, v, CGS);

/*  Phase 7. - Test for spikes */
        if ( nstep == 1 ) {
            vm2 = 0.0;
            vm1 = v[nodes-1];
        } else {
            vm0 = v[nodes-1];
            if ( !spk ) {
                spk = ( vm0 > 50.0 && vm1 > vm2 && vm1 > vm0 );
                if ( spk ) {
                    if ( nstep >= 1000*NT ) {
                        nspk++;
                        tnow = dt*((double) nstep)-1000.0;
                        tmp = tnow+0.5*dt*(vm2+3.0*vm0-4.0*vm1)/(vm0-2.0*vm1+vm2);
                        nn = ((int) floor(tmp));
                        if ( fmod(tmp,1.0)>0.5 ) nn++;
                        fp = fopen(OUTPUT,"a");
                        fprintf(fp,"%d\n",nn);
                        fclose(fp);
                    } else {
                        printf("\nUnrecorded spike fired\n");
                    }
                }
            }
            vm2 = vm1;
            vm1 = vm0;
        }

/*  Phase 8. - Reset spike flag */
        if ( vs < 0.0 && spk == 1 ) spk = 0;
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
    int k;
    synapse *syn;

    if ( bstart == bnow ) n = 0;
    if ( bnow->child ) Count_Synapses(bstart, bnow->child);
    if ( bnow->peer ) Count_Synapses(bstart, bnow->peer);

    for ( k=0 ; k<bnow->nseg ; k++ ) {
        syn = (bnow->SegList)[k].SynList;
        while ( syn ) {
            n++;
            syn = syn->NextSyn;
        }
    }
    return n;
}


 /*******************************************************
       Function to enumerate the nodes on a dendrite
  *******************************************************/
void Enumerate_Nodes(branch *bnow, int *FirstNode )
{
    int node, k;
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
    *FirstNode += bnow->nseg;
    bnow->first = *FirstNode-1;
    return;
}


 /*******************************************************
                 Function to output information
  *******************************************************/
void Output_Information(branch *b, FILE *fp)
{
    int k;
    double pos;
    synapse *syn;
    segment *seg;

    if ( b->child ) Output_Information( b->child, fp);
    if ( b->peer ) Output_Information( b->peer, fp);

    for ( k=0 ; k<b->nseg ; k++ ) {
        seg = &(b->SegList)[k];
        syn = seg->SynList;
        while ( syn ) {
            pos = syn->vlam+((double) k);
            fprintf(fp," %12.6lf \t %d \n",pos, syn->SynTime);
            syn = syn->NextSyn;
        }
    }
    fprintf(fp,"End of branch\n");
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
    nc = b->nseg;
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
    extern double pi;
    int j, jj, k, nsyn;
    double fac, interval, vold, vnew;
    synapse *syn;
    segment *seg;

    if ( b->child ) Initialise_Synapses( b->child );
    if ( b->peer ) Initialise_Synapses( b->peer );

    for ( k=0 ; k<b->nseg ; k++ ) {

/*  Phase 1. - Allocate storage for synaptic activity */
        nsyn = (b->SegList)[k].nsyn;
        if ( nsyn > 0 ) {
            seg = &(b->SegList)[k];
            seg->vlam = (double *) malloc( (nsyn+1)*sizeof(double) );
            seg->diag = (double *) malloc( nsyn*sizeof(double) );
            seg->hseg = b->hseg;
            fac = (seg->hseg)/(pi*(seg->rp)*(seg->rd)*GA);
            vnew = 0.0;
            syn = seg->SynList;
            for ( j=0 ; j<nsyn ; j++ ) {
                (seg->vlam)[j] = vnew;
                vold = vnew;
                vnew = syn->vlam;
                (seg->diag)[j] = fac*(vnew-vold);
                syn = syn->NextSyn;
            }
            seg->vlam[nsyn] = vnew;

/*  Phase 2. - Allocate phi1, phi2, phi3 and diag */
            seg->phi1 = (double *) malloc( nsyn*sizeof(double) );
            seg->phi2 = (double *) malloc( nsyn*sizeof(double) );
            seg->phi3 = (double *) malloc( nsyn*sizeof(double) );
            seg->soln = (double *) malloc( (nsyn+1)*sizeof(double) );
            seg->vec = (double *) malloc( (nsyn+1)*sizeof(double) );
            seg->phi = (double *) malloc( (nsyn+1)*sizeof(double) );
            syn = seg->SynList;
            for ( j=0 ; j<nsyn ; j++ ) {
                (seg->phi1)[j] = 1.0-(syn->vlam);
                (seg->phi2)[j] = syn->vlam;
                (seg->phi3)[j] = -(syn->vsyn);
                syn = syn->NextSyn;
            }

/*  Phase 3. - Allocate and initialise conductances and firing times */
            seg->gold = (double *) malloc( nsyn*sizeof(double) );
            seg->gnew = (double *) malloc( nsyn*sizeof(double) );
            syn = seg->SynList;
            j = 0;
            while ( syn ) {
                (seg->gold)[j] = (seg->gnew)[j] = 0.0;
                interval = -(((double) NT*1000)/RATE)*log(ran(&ix, &iy, &iz));
                if ( fmod(interval,1.0) <= 0.5 ) {
                    syn->SynTime = -((int) floor(interval));
                } else {
                    syn->SynTime = -((int) ceil(interval));
                }
                syn = syn->NextSyn;
                j++;
            }

/*  Phase 4. - Initialise proximal and distal currents */
            seg->NewProx1 = 0.0;
            seg->NewProx2 = 0.0;
            seg->NewProx3 = 0.0;
            seg->NewDist1 = 0.0;
            seg->NewDist2 = 0.0;
            seg->NewDist3 = 0.0;
        }
    }
    return;
}


 /***********************************************
       Function to update status of synapses
  ***********************************************/
void Update_Synapses( branch *b )
{
    int i, j, k, nsyn;
    double interval, tmp;
    synapse *syn;
    segment *seg;

    if ( b->child ) Update_Synapses( b->child );
    if ( b->peer ) Update_Synapses( b->peer );

    for ( k=0 ; k<b->nseg ; k++ ) {
        seg = &(b->SegList)[k];

/*  Phase 1. - Update synaptic conductances */
        syn = seg->SynList;
        j = 0;
        while ( syn ) {
            (seg->gold)[j] = (seg->gnew)[j];
            (syn->SynTime)++;
            if ( syn->SynTime < 1 ) {
                (seg->gnew)[j] = 0.0;
            } else if ( syn->SynTime < syn->MaxStep ) {
                (seg->gnew)[j] = (syn->SynCond->g)[syn->SynTime];
            } else {
                (seg->gnew)[j] = 0.0;
                interval = -(((double) NT*1000)/RATE)*log(ran(&ix, &iy, &iz));
                if ( fmod(interval,1.0) <= 0.5 ) {
                    syn->SynTime = -((int) floor(interval));
                } else {
                    syn->SynTime = -((int) ceil(interval));
                }
            }
            syn = syn->NextSyn;
            j++;
        }

/*  Phase 2. - Update currents at segment endpoints */
        nsyn = seg->nsyn;
        if ( nsyn > 0 ) {
            seg->OldProx1 = seg->NewProx1;
            seg->OldProx2 = seg->NewProx2;
            seg->OldProx3 = seg->NewProx3;
            seg->OldDist1 = seg->NewDist1;
            seg->OldDist2 = seg->NewDist2;
            seg->OldDist3 = seg->NewDist3;

/*  Phase 3a. - Coefficient of VP */
            for ( j=0 ; j<nsyn ; j++ ) {
                (seg->phi)[j] = (seg->gnew)[j]*(seg->phi1)[j];
            }
            (seg->phi)[nsyn] = 0.0;
            solve( nsyn, seg->vlam, seg->phi, seg->soln);

/*  Phase 3b. - Repeat once */
            for ( j=0 ; j<nsyn ; j++ ) {
                for ( tmp=0.0,i=0 ; i<=j ; i++ ) {
                    tmp += (seg->diag)[i]*(seg->soln)[i];
                }
                (seg->phi)[j] -= tmp*(seg->gnew)[j];
            }
            solve( nsyn, seg->vlam, seg->phi, seg->soln);
            seg->NewProx1 = (seg->soln)[0];
            seg->NewDist1 = (seg->soln)[nsyn];

/*  Phase 4a. - Coefficient of VD */
            for ( j=0 ; j<nsyn ; j++ ) {
                (seg->phi)[j] = (seg->gnew)[j]*(seg->phi2)[j];
            }
            (seg->phi)[nsyn] = 0.0;
            solve( nsyn, seg->vlam, seg->phi, seg->soln);

/*  Phase 4b. - Repeat once */
            for ( j=0 ; j<nsyn ; j++ ) {
                for ( tmp=0.0,i=0 ; i<=j ; i++ ) {
                    tmp += (seg->diag)[i]*(seg->soln)[i];
                }
                (seg->phi)[j] -= tmp*(seg->gnew)[j];
            }
            solve( nsyn, seg->vlam, seg->phi, seg->soln);
            seg->NewProx2 = (seg->soln)[0];
            seg->NewDist2 = (seg->soln)[nsyn];

/*  Phase 5a. - Constant term */
            for ( j=0 ; j<nsyn ; j++ ) {
                (seg->phi)[j] = (seg->gnew)[j]*(seg->phi3)[j];
            }
            (seg->phi)[nsyn] = 0.0;
            solve( nsyn, seg->vlam, seg->phi, seg->soln);

/*  Phase 5b. - Repeat once */
            for ( j=0 ; j<nsyn ; j++ ) {
                for ( tmp=0.0,i=0 ; i<=j ; i++ ) {
                    tmp += (seg->diag)[i]*(seg->soln)[i];
                }
                (seg->phi)[j] -= tmp*(seg->gnew)[j];
            }
            solve( nsyn, seg->vlam, seg->phi, seg->soln);
            seg->NewProx3 = (seg->soln)[0];
            seg->NewDist3 = (seg->soln)[nsyn];
        }
    }
    return;
}


 /***************************************************
        Function to constuct sparse matrices
  ***************************************************/
void Update_Sparse_Matrices( branch *b, int *counter)
{
    int k, CurrentNode, nc;
    extern double pi, dt, *NewAmp;
    extern SparseMatrix lhs, rhs;
    branch *btmp;
    segment NewSeg, OldSeg;
    double SumL, SumR, tmp;

/* Step 1 - Recurse to the end of the dendrite */
    if ( b->child ) Update_Sparse_Matrices( b->child, counter);
    if ( b->peer ) Update_Sparse_Matrices( b->peer, counter);

/* Step 2 - Build matrix entries for distal node of branch */
    nc = b->nseg;
    tmp = 0.5*dt;
    CurrentNode = (b->first)-(nc-1);
    if ( b->child ) {
        btmp = b->child;
        SumR = SumL = 0.0;
        while ( btmp ) {
            NewSeg = (btmp->SegList)[0];
            lhs.a[*counter] = tmp*(NewSeg.NewProx2);
            rhs.a[*counter] = -tmp*(NewSeg.OldProx2);
            NewAmp[CurrentNode] += NewSeg.NewProx3;
            SumL += NewSeg.NewProx1;
            SumR += NewSeg.OldProx1;
            (*counter)++;
            btmp = btmp->peer;
        }
        NewSeg = (b->SegList)[nc-1];
        lhs.a[*counter] = tmp*(SumL-NewSeg.NewDist2);
        rhs.a[*counter] = -tmp*(SumR-NewSeg.OldDist2);
        NewAmp[CurrentNode] -= NewSeg.NewDist3;
        (*counter)++;
        lhs.a[*counter] = -tmp*(NewSeg.NewDist1);
        rhs.a[*counter] = tmp*(NewSeg.OldDist1);
        (*counter)++;
    } else {
        NewSeg = (b->SegList)[nc-1];
        lhs.a[*counter] = -tmp*(NewSeg.NewDist2);
        rhs.a[*counter] = tmp*(NewSeg.OldDist2);
        NewAmp[CurrentNode] = -NewSeg.NewDist3;
        (*counter)++;
        lhs.a[*counter] = -tmp*(NewSeg.NewDist1);
        rhs.a[*counter] = tmp*(NewSeg.OldDist1);
        (*counter)++;
    }

/* Step 3 - Build matrix entries for internal nodes of branch */
    for ( k=nc-1 ; k>0 ; k-- ) {
        CurrentNode++;
        OldSeg = NewSeg;
        NewSeg = (b->SegList)[k-1];
        lhs.a[*counter] = tmp*(OldSeg.NewProx2);
        rhs.a[*counter] = -tmp*(OldSeg.OldProx2);
        (*counter)++;
        lhs.a[*counter] = tmp*(OldSeg.NewProx1-NewSeg.NewDist2);
        rhs.a[*counter] = -tmp*(OldSeg.OldProx1-NewSeg.OldDist2);
        NewAmp[CurrentNode] = OldSeg.NewProx3-NewSeg.NewDist3;
        (*counter)++;
        lhs.a[*counter] = -tmp*(NewSeg.NewDist1);
        rhs.a[*counter] = tmp*(NewSeg.OldDist1);
        (*counter)++;
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


void cgs(int *getmem, SparseMatrix *a, double *b, double *x, double tol)
{
    long int i, k, n, repeat;
    static int start=1;
    double rho_old, rho_new, alpha, beta, sigma, err;
    static double *p, *q, *r, *u, *v, *rb, *y;

/* Step 1 - Check memory status */
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
        if ( sigma == 0.0 ) {
            printf(" Trouble ");
            for ( i=0 ; i<n ; i++ ) {
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
              if ( rho_old == 0.0 ) {
            printf(" Trouble rho_old ");
            for ( i=0 ; i<n ; i++ ) {
                printf("\n%20.16lf",v[i]);
                getchar( );
            }
        }
        repeat = ( err > tol*((double) n) );
    } while ( repeat );
    return;
}

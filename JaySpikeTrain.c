#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


 /***************************************************************
    Traditional compartmental model with large scale synaptic
    input and active soma. Spike generating program
  ***************************************************************/

/* Global definitions */
#define             CS    1.0
#define             GS    0.091
#define             GA    14.286
#define             CM    1.0
#define             GM    0.091
#define         OUTPUT    "JaySpike300.dat"
#define           TEND    11             /* Must be an integer */
#define           GMAX    3.0e-5
#define           RATE    30.0
#define           VSYN    115.0
#define             NT    1000
#define          NODES    300
#define          NSEED    2             /* Seed for random number generator */
#define        CELSIUS    18.5          /* Celsius temperature of neuron */
#define            CGS    1.0e-24       /* Tolerance used in CGS algorithm */
#define            TAU    0.5

/* Parameters for exact solution */
#define           NSYN    100          /* Number of synapses */
#define             RS    0.005


typedef struct SparseMatrix_t
{
        double *a;
        int *col;
        int *StartRow;
        int n;
} SparseMatrix;


typedef struct cond_t
{
        int n;                      /* Number of time steps in conductance profile */
        double dt;                  /* Integration time step (msec) */
        double *g;                  /* Conductance profile for prescribed dt (mu a) */
        double tol;                 /* Determines off threshold for synapse */
        double tau;                 /* Rise time for synapse (msec) */
        double gmax;                /* Peak conductance per synapse */
} cond;


typedef struct synapse_t
{
        int id;                     /* Identifies contact type */
        double xp;                  /* Location of contact */
        int cn;                     /* Central node of compartment*/
        int diag;                   /* Diagonal entry in sparse matrix corresponding to cn */
        double gval;                /* Present conductance */

/* Properties of synaptic conductance profile */
        cond *SynCond;              /* Address of synaptic conductance profile */
        int SynTime;                /* Time to activation of syapse/time since activation */
        int MaxStep;                /* Time to inactivation of synapse */
        double vsyn;                /* Reversal potential of synaptic species */
        struct synapse_t *NextSyn;  /* Address of next synapse */
} synapse;


typedef struct branch_t
{
/*  Connectivity of branch */
        struct branch_t *parent;    /* Pointer to parent branch */
        struct branch_t *child;     /* Pointer to child branch */
        struct branch_t *peer;      /* Pointer to a peer branch */

/*  Physical properties of branch */
        int id;                     /* Branch identifier */
        double xl;                  /* X-coordinate of lefthand endpoint */
        double yl;                  /* Y-coordinate of lefthand endpoint */
        double zl;                  /* Z-coordinate of lefthand endpoint */
        double xr;                  /* X-coordinate of righthand endpoint */
        double yr;                  /* Y-coordinate of righthand endpoint */
        double zr;                  /* Z-coordinate of righthand endpoint */
        double diam;                /* Branch diameter (cm) */
        double plen;                /* Branch length (cm) */
        double hseg;                /* Compartment length (cm) */

/*  Node information for spatial representation */
        int nc;                     /* Number of compartments in branch */
        int junct;                  /* Junction node of the branch */
        int first;                  /* Internal node connected to junction */

/*  Synapse information */
        synapse *synlist;           /* Branch contact */
} branch;


typedef struct dendrite_t
{
        branch *root;               /* Pointer to root branch of dendrite */
        double plen;               /* Total length of dendrite */
} dendrite;


typedef struct neuron_t
{
        int ndend;                  /* Number of dendrites */
        dendrite *dendlist;         /* Pointer to an array of dendrites */
} neuron;


/* Function type declarations */
int     Count_Branches( branch *, branch *),
        Count_Synapses( branch *, branch *);

cond *ConductanceProfile( double, double,  double, double );

double  branch_length( branch *, branch *),
        ran(unsigned int *, unsigned int *, unsigned int *),
        alfa_h( double ),
        alfa_m( double ),
        alfa_n( double ),
        beta_h( double ),
        beta_m( double ),
        beta_n( double );

void    Remove_Branch( branch **, branch *),
        Output_Information( branch *, FILE *),
        Assign_Branch_Nodes( branch *, double *),
        Enumerate_Nodes( branch *, int *),
        Matrix_Vector_Multiply( SparseMatrix *, double *, double *),
        Matrix_Malloc( SparseMatrix *, int, int),
        Matrix_Free( SparseMatrix *),
        Generate_Dendrite(branch *, int *),
        Build_Test_Dendrite( branch **, branch *),
        Initialise_Synapses( branch *),
        Update_Synapses( branch *),
        cgs( int *, SparseMatrix *, double *, double *, double);

/* Global Variables */
SparseMatrix lhs, rhs;
double pi, dt, *NewAmp;
unsigned int ix, iy, iz;

int main( int argc, char **argv )
{
    extern unsigned int ix, iy, iz;
    extern double pi, dt, *NewAmp;
    extern SparseMatrix lhs, rhs;
    int j, k, nodes, n, nn, FirstNode, repeat, maxstep, spk, nspk,
        nstep, nb, nsyn, getmem, connected, counter, NumberOfSynapses;
    double *v, *x, *OldAmp, AreaOfSoma, tmp, h, len, CellLength, tnow,
        gs, aval, bval, hval, nval, mval, vs, sc, vm0, vm1, vm2,
        LocusContact;
    void srand( unsigned int);
    double *StoredLHS, *StoredRHS;
    double v_na=115.0, v_k=-12.0, v_l, g_na=120.0, g_k=36.0, g_l=0.3;
    synapse *newsyn, *syn;
    neuron *cell;
    cond *SynCond;
    branch *bold, *bnew, *FirstBranch;
    char word[100];
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
    }
    printf("\nOpening file %s\n",argv[1]);
    if ( (fp=fopen(argv[1],"r")) == NULL ) {
        printf("\n Test Neuron file does not found");
        return 1;
    }

/*  Get branch data */
    bold = NULL;
    while  ( fscanf(fp,"%s", word) != EOF ) {
        if ( strcmp(word, "Branch") == 0 || strcmp(word, "branch") == 0 ) {
            fscanf(fp, "%d", &nodes);
            printf("Found branch %d \n", nodes);
            bnew = (branch *) malloc( sizeof(branch) );
            bnew->id = nodes;
            bnew->peer = NULL;
            bnew->child = NULL;
            bnew->synlist = NULL;
            if ( bold != NULL) {
                bold->child = bnew;
            } else {
                FirstBranch = bnew;
            }
            bnew->parent = bold;
            fscanf(fp,"%lf %lf %lf", &(bnew->xl), &(bnew->yl), &(bnew->zl));
            fscanf(fp,"%lf %lf %lf", &(bnew->xr), &(bnew->yr), &(bnew->zr));
            fscanf(fp,"%lf %lf", &(bnew->plen), &(bnew->diam));
            bold = bnew;
        } else {
            printf("Unrecognised dendritic feature\n");
            return 0;
        }
    }
    fclose(fp);

/*  Compute total length of dendrite */
    CellLength = 0.0;
    bnew = FirstBranch;
    while ( bnew ) {
       CellLength += bnew->plen;
       bnew = bnew->child;
    }

/*  STEP 1. - Randomly place NSYN synapses on branches */
    for ( k=0 ; k<NSYN ; k++ ) {
        LocusContact = CellLength*ran( &ix, &iy, &iz);
        bnew = FirstBranch;
        len = bnew->plen;
        while ( LocusContact > len ) {
            bnew = bnew->child;
            len += bnew->plen;
        }
        newsyn = (synapse *) malloc( sizeof(synapse) );
        newsyn->NextSyn = NULL;
        newsyn->SynCond = SynCond;
        newsyn->MaxStep = SynCond->n;
        newsyn->vsyn = VSYN;
        newsyn->xp = LocusContact-(len-bnew->plen);
        syn = bnew->synlist;
        if ( syn ) {
            if ( newsyn->xp < syn->xp ) {
                newsyn->NextSyn = syn;
                bnew->synlist = newsyn;
            } else {
                while ( syn->NextSyn && newsyn->xp > syn->NextSyn->xp ) syn = syn->NextSyn;
                if ( syn->NextSyn ) newsyn->NextSyn = syn->NextSyn;
                syn->NextSyn = newsyn;
            }
        } else {
            bnew->synlist = newsyn;
        }
    }

/*  STEP 2. - Count root branches */
    bold = FirstBranch;
    n = 0;
    while ( bold ) {
        bnew = FirstBranch;
        do {
            tmp = pow(bold->xl-bnew->xr,2)+
                  pow(bold->yl-bnew->yr,2)+
                  pow(bold->zl-bnew->zr,2);
            connected = ( tmp < 0.01 );
            bnew = bnew->child;
        } while ( bnew && !connected );
        if ( !connected ) n++;
        bold = bold->child;
    }

/*  STEP 3. - Identify somal dendrites but extract nothing */
    printf("\n\nTree contains %d individual dendrite(s) ...\n", n);
    cell = (neuron *) malloc( sizeof(neuron) );
    cell->ndend = n;
    cell->dendlist = (dendrite *) malloc( n*sizeof(dendrite) );
    bold = FirstBranch;
    n = 0;
    while ( n < cell->ndend ) {
        bnew = FirstBranch;
        do {
            tmp = pow(bold->xl-bnew->xr,2)+
                  pow(bold->yl-bnew->yr,2)+
                  pow(bold->zl-bnew->zr,2);
            connected = ( tmp < 0.01 );
            bnew = bnew->child;
        } while ( bnew && !connected );
        if ( !connected ) {
            cell->dendlist[n].root = bold;
            n++;
        }
        bold = bold->child;
    }

/*  STEP 4. - Extract root of each dendrite from dendrite list */
    for ( k=0 ; k<cell->ndend ; k++ ) {
        bold = cell->dendlist[k].root;
        Remove_Branch( &FirstBranch, bold);
    }

/*  STEP 5. - Build each test dendrite from its root branch */
    for ( k=0 ; k<cell->ndend ; k++ ) {
        Build_Test_Dendrite( &FirstBranch, cell->dendlist[k].root );
    }
    if ( FirstBranch != NULL ) printf("\nWarning: Unconnected branch segments still exist\n");

/*  STEP 6. - Count number of synapses on Cell */
    NumberOfSynapses = 0;
    for ( k=0 ; k<cell->ndend ; k++ ) {
        bnew = cell->dendlist[k].root;
        NumberOfSynapses += Count_Synapses( cell->dendlist[k].root, bnew);
    }
    printf("\nNumber of Synapses %d", NumberOfSynapses);

/*  STEP 7. - Count dendritic segments */
    for ( nb=k=0 ; k<cell->ndend ; k++ ) {
        bnew = cell->dendlist[k].root;
        nb += Count_Branches( bnew, bnew);
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
        bnew = cell->dendlist[k].root;
        Generate_Dendrite( bnew, &counter);
    }
    lhs.n = rhs.n = nodes;

/*  STEP 10. - Do somal node */
    lhs.a[3*nodes-3] = rhs.a[3*nodes-3] = 0.0;
    for ( k=0 ; k<cell->ndend ; k++ ) {
        bnew = cell->dendlist[k].root;
        lhs.a[counter] = 0.0;
        rhs.a[counter] = -0.25*pi*pow(bnew->diam,2)/(bnew->hseg);
        rhs.a[3*nodes-3] += 0.25*pi*pow(bnew->diam,2)/(bnew->hseg);
        lhs.a[3*nodes-3] += 0.5*pi*(bnew->diam)*(bnew->hseg);
        lhs.col[counter] = rhs.col[counter] = bnew->first;
        counter++;
    }
    lhs.col[3*nodes-3] = rhs.col[3*nodes-3] = nodes-1;
    lhs.StartRow[nodes] = rhs.StartRow[nodes] = 3*nodes-2;

/*  STEP 11. - Build and store constant components of left and right hand matrices  */
        for ( k=0 ; k<3*nodes-2 ; k++ ) {
        rhs.a[k] = 0.5*dt*(GA*rhs.a[k]+GM*lhs.a[k]);
        rhs.a[k] = CM*lhs.a[k]-rhs.a[k];
        lhs.a[k] = 2.0*CM*lhs.a[k]-rhs.a[k];
    }

/*  Add capacitive term of soma */
    lhs.a[3*nodes-3] += AreaOfSoma*CS;
    rhs.a[3*nodes-3] += AreaOfSoma*CS;
    StoredLHS = (double *) malloc( (3*nodes-2)*sizeof(double) );
    StoredRHS = (double *) malloc( (3*nodes-2)*sizeof(double) );
    for ( k=0 ; k<3*nodes-2 ; k++ ) {
        StoredLHS[k] = lhs.a[k];
        StoredRHS[k] = rhs.a[k];
    }

//    fp1 = fopen("JayOutput.dat","w");
//    for( k=0 ; k<3*nodes-2 ; k++ ) {
//        fprintf(fp1,"%12.10lf\n", lhs.a[k]);
//    }
//    fclose(fp1);
//    printf("\nOutput information");
//    getchar( );

/*  Step 12. - Initialise synapses, soma and membrane potentials */
    for ( k=0 ; k<cell->ndend ; k++ ) {
        bnew = cell->dendlist[k].root;
        Initialise_Synapses(bnew);
    }
    v = (double *) malloc( (nodes)*sizeof(double) );
    x = (double *) malloc( (nodes)*sizeof(double) );
    OldAmp = (double *) malloc( (nodes)*sizeof(double) );
    NewAmp = (double *) malloc( (nodes)*sizeof(double) );
    for ( k=0 ; k<nodes ; k++ ) v[k] = 0.0;

    g_na *= AreaOfSoma;
    g_k *= AreaOfSoma;
    g_l *= AreaOfSoma;
    hval = alfa_h(0.0)/(alfa_h(0.0)+beta_h(0.0));
    mval = alfa_m(0.0)/(alfa_m(0.0)+beta_m(0.0));
    nval = alfa_n(0.0)/(alfa_n(0.0)+beta_n(0.0));
    v_l = g_na*pow(mval,3)*hval*v_na+g_k*pow(nval,4)*v_k;
    v_l = -v_l/g_l;

/*  Step 13. - Run simulation */
    nstep = 0;
    getmem = 1;
    while ( nstep < maxstep ) {
        nstep++;
//        if ( nstep%1000 == 0 ) printf("\r Msec %d",nstep/1000);

/*  Phase 1. - Update Synaptic activity */
        for ( k=0 ; k<3*nodes-2 ; k++ ) lhs.a[k] = rhs.a[k] = 0.0;
        for ( k=0 ; k<nodes ; k++ ) {
            OldAmp[k] = NewAmp[k];
            NewAmp[k] = 0.0;
        }
        for ( k=0 ; k<cell->ndend ; k++ ) {
            bnew = cell->dendlist[k].root;
            Update_Synapses(bnew);
        }
        for ( k=0 ; k<3*nodes-2 ; k++ ) {
            lhs.a[k] += StoredLHS[k];
            rhs.a[k] += StoredRHS[k];
        }

/*  Phase 2. - Update HH channel variables */
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

/* Phase 3. - Compute somal conductance and contribution to current */
        gs = g_l;                       /* Leakage conductance */
        sc = g_l*v_l;                   /* Leakage contribution to somal current */
        tmp = g_na*hval*pow(mval,3);    /* Sodium conductance */
        gs += tmp;
        sc += tmp*v_na;                 /* Sodium contribution to somal current */
        tmp = g_k*pow(nval,4);          /* Potassium conductance */
        gs += tmp;
        sc += tmp*v_k;                  /* Potasium contribution to somal current */

/*  Phase 4. - Complete somal entry of LHS and RHS matrices */
        gs *= 0.5*dt;
        lhs.a[3*nodes-3] += gs;
        rhs.a[3*nodes-3] -= gs;

/*  Phase 5. - Step potential forward */
        Matrix_Vector_Multiply( &rhs, v, x);
        x[nodes-1] += sc*dt;
        tmp = 0.5*dt;
        for ( k=0 ; k<nodes ; k++ ) x[k] -= tmp*(OldAmp[k]+NewAmp[k]);
        cgs( &getmem, &lhs, x, v, CGS);

/* Phase 6. - Test for spikes */
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

/*  Phase 7. - Reset spike flag */
        if ( vs < 0.0 && spk == 1 ) spk = 0;
        if ( nstep%(500*NT) == 0 ) {
            tnow = dt*((double) nstep/NT);
            printf("\rReached time %5.1lf \t Spikes so far %d", tnow, nspk);
        }
    }
    return 0;
}


 /*******************************************************
                 Function to output information
  *******************************************************/
void Output_Information(branch *b, FILE *fp)
{
    int k;
    double pos;
    synapse *syn;

    if ( b->child ) Output_Information( b->child, fp);
    if ( b->peer ) Output_Information( b->peer, fp);

    syn = b->synlist;
    while ( syn ) {
        pos = (syn->xp)/(b->hseg);
        fprintf(fp," %12.6lf \t %d \n",pos, syn->SynTime);
        syn = syn->NextSyn;
    }
    fprintf(fp,"End of branch\n");
    return;
}


 /******************************************************
        Function to constuct sparse matrices
  ******************************************************/
void Generate_Dendrite( branch *b, int *counter)
{
    int k, nc, CurrentNode;
    extern double pi;
    extern SparseMatrix lhs, rhs;
    branch *btmp;
    double SumL, SumR;

/* Step 1 - Recurse to the end of the dendrite */
    if ( b->child ) Generate_Dendrite( b->child, counter);
    if ( b->peer ) Generate_Dendrite( b->peer, counter);

/* Step 2 - Fill in matrix entries for a branch point node */
    nc = b->nc;
    CurrentNode = (b->first)-(nc-1);
    if ( b->child ) {
        btmp = b->child;
        SumR = SumL = 0.0;
        while ( btmp ) {
            lhs.a[*counter] = 0.0;
            rhs.a[*counter] = -0.25*pi*pow(btmp->diam,2)/(btmp->hseg);
            SumL += 0.5*pi*(btmp->diam)*(btmp->hseg);
            SumR += 0.25*pi*pow(btmp->diam,2)/(btmp->hseg);
            lhs.col[*counter] = rhs.col[*counter] = btmp->first;
            (*counter)++;
            btmp = btmp->peer;
        }
        lhs.a[*counter] = SumL+0.5*pi*(b->diam)*(b->hseg);
        rhs.a[*counter] = SumR+0.25*pi*pow(b->diam,2)/(b->hseg);
        lhs.col[*counter] = rhs.col[*counter] = CurrentNode;
        (*counter)++;
        lhs.a[*counter] = 0.0;
        rhs.a[*counter] = -0.25*pi*pow(b->diam,2)/(b->hseg);
        if ( CurrentNode == b->first ) {
            lhs.col[*counter] = rhs.col[*counter] = b->junct;
        } else {
            lhs.col[*counter] = rhs.col[*counter] = CurrentNode+1;
        }
        (*counter)++;
        lhs.StartRow[CurrentNode+1] = rhs.StartRow[CurrentNode+1] = *counter;
    } else {

/* Step 3 - Fill in matrix entries for a terminal node */
        lhs.a[*counter] = 0.5*pi*(b->diam)*(b->hseg);
        rhs.a[*counter] = 0.25*pi*pow(b->diam,2)/(b->hseg);
        lhs.col[*counter] = rhs.col[*counter] = CurrentNode;
        (*counter)++;
        lhs.a[*counter] = 0.0;
        rhs.a[*counter] = -0.25*pi*pow(b->diam,2)/(b->hseg);
        if ( CurrentNode == b->first ) {
            lhs.col[*counter] = rhs.col[*counter] = b->junct;
        } else {
            lhs.col[*counter] = rhs.col[*counter] = CurrentNode+1;
        }
        (*counter)++;
        lhs.StartRow[CurrentNode+1] = rhs.StartRow[CurrentNode+1] = *counter;
    }

/* Step 4 - Fill in matrix entries for an internal node */
    for ( k=nc-1 ; k>0 ; k-- ) {
        CurrentNode++;
        lhs.a[*counter] = 0.0;
        rhs.a[*counter] = -0.25*pi*pow(b->diam,2)/(b->hseg);
        lhs.col[*counter] = rhs.col[*counter] = CurrentNode-1;
        (*counter)++;
        lhs.a[*counter] = pi*(b->diam)*(b->hseg);
        rhs.a[*counter] = 0.5*pi*pow(b->diam,2)/(b->hseg);
        lhs.col[*counter] = rhs.col[*counter] = CurrentNode;
        (*counter)++;
        lhs.a[*counter] = 0.0;
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


 /*************************************************************
          Function to remove a branch from a branch list
  *************************************************************/
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
      Function to count synapses on branches
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
            syn = syn->NextSyn;
        }
    }
    return n;
}


 /*******************************************************
             Function to count number of branches
  *******************************************************/
int Count_Branches( branch *bstart, branch *bnow)
{
    static int n;

    if ( bstart == bnow ) n = 0;
    if ( bnow ) {
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

    btmp = bnow->child;
    while ( btmp ) {
        btmp->junct = *FirstNode;
        btmp = btmp->peer;
    }
    *FirstNode += bnow->nc;
    bnow->first = *FirstNode-1;
    return;
}


 /*************************************************************
        Function to build a test dendrite from its root
  *************************************************************/
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


 /**********************************************************
           Function to initialist synaptic input to dendrite
  **********************************************************/
void Initialise_Synapses( branch *b)
{
    extern unsigned int ix, iy, iz;
    extern SparseMatrix lhs, rhs;
    int k;
    double ratio, interval;
    synapse *syn;

    if ( b->child ) Initialise_Synapses( b->child );
    if ( b->peer ) Initialise_Synapses( b->peer );

    syn = b->synlist;
    while ( syn ) {
        ratio = (syn->xp)/(b->hseg);
        if ( fmod(ratio,1.0) >= 0.5 ) {
            k = ((int) ceil(ratio));
        } else {
            k = ((int) floor(ratio));
        }
        if ( k == 0 ) {
            syn->cn = b->junct;
        } else {
            syn->cn = b->first-k+1;
        }
        k = lhs.StartRow[syn->cn];
        while ( lhs.col[k] != syn->cn ) k++;
        syn->diag = k;
        syn->gval = 0.0;
        interval = -(((double) NT*1000)/RATE)*log(ran(&ix, &iy, &iz));
        if ( fmod(interval,1.0) <= 0.5 ) {
            syn->SynTime = -((int) floor(interval));
        } else {
            syn->SynTime = -((int) ceil(interval));
        }
        syn = syn->NextSyn;
    }
    return;
}


 /**********************************************************
           Function to update synaptic input to dendrite
  **********************************************************/
void Update_Synapses( branch *b)
{
    extern unsigned int ix, iy, iz;
    extern double dt, *NewAmp;
    extern SparseMatrix lhs, rhs;
    int k;
    double interval, tmp, gnew, gold;
    synapse *syn;

    if ( b->child ) Update_Synapses( b->child);
    if ( b->peer ) Update_Synapses( b->peer);

    syn = b->synlist;
    while ( syn ) {
        gold = syn->gval;
        (syn->SynTime)++;
        if ( syn->SynTime < 1 ) {
            syn->gval = 0.0;
        } else if ( syn->SynTime < syn->MaxStep ) {
            syn->gval = syn->SynCond->g[syn->SynTime];
        } else {
            interval = -(((double) NT*1000)/RATE)*log(ran(&ix, &iy, &iz));
            syn->gval = 0.0;
            if ( fmod(interval,1.0) <= 0.5 ) {
                syn->SynTime = -((int) floor(interval));
            } else {
                syn->SynTime = -((int) ceil(interval));
            }
        }
        gnew = syn->gval;
        tmp = 0.5*dt;
        lhs.a[syn->diag] += tmp*gnew;
        rhs.a[syn->diag] -= tmp*gold;
        NewAmp[syn->cn] -= gnew*(syn->vsyn);
        syn = syn->NextSyn;
    }
    return;
}


 /*********************************************************
            Function to assign branch nodes
  *********************************************************/
void Assign_Branch_Nodes( branch *b, double *h )
{
    b->nc = ((int) ceil((b->plen)/(*h)));     /* nc is the number of compartments in branch */
    b->hseg = (b->plen)/((double) b->nc);     /* hseg is the compartment length */

    if ( b->child ) Assign_Branch_Nodes( b->child, h);
    if ( b->peer ) Assign_Branch_Nodes( b->peer, h);
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
/*  Generate random number in (0,1)  */
    tmp = ((double) (*ix))/30269.0+((double) (*iy))/30307.0
          +((double) (*iz))/30323.0;
    return fmod(tmp,1.0);
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

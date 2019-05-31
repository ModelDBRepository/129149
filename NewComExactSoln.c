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


typedef struct contact_t
{
        int id;                     /* Identifies contact type */

        double xp;                  /* Location of contact */
        double amp;                 /* Strength of contact */

        int xl;                     /* Neighbouring node on left */
        int xr;                     /* Neighbouring node on right */
        double fl;                  /* Fraction of input assigned to LH node */
        double fr;                  /* Fraction of input assigned to RH node */

        struct contact_t *next;     /* Address of next contact */
} contact;


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
        double elen;                /* Branch length (eu) */
        double hseg;                /* Dendritic segment length (cm) */

/*  Node information for spatial representation */
        int nodes;                  /* Total number nodes spanning branch */
        int junct;                  /* Junction node of the branch */
        int first;                  /* Internal node connected to junction */

/*  Contact information */
        contact *conlist;           /* Branch contact */
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
int     Count_Branches( branch *, branch *),
        Count_Contacts( branch *, branch *);

double  branch_length( branch *, branch *),
        ran(unsigned int *, unsigned int *, unsigned int *);

void    Build_Test_Dendrite( branch **, branch *),
        Remove_Branch( branch **, branch *),
        Output_Info( branch *, FILE *),
        Destroy_Test_Neuron( neuron *),
        Destroy_Test_Dendrite( branch *),
        Find_Contacts( branch *, double *, double *, int *),
        Assign_Branch_Nodes( branch *, double *),
        Enumerate_Nodes( branch *, int *),
        Generate_Dendrite(branch *, int *),
        Input_Current( branch *),
        Assign_Current( branch *, double *, double ),
        Matrix_Vector_Multiply( SparseMatrix *, double *, double *),
        Matrix_Malloc( SparseMatrix *, int, int),
        Matrix_Free( SparseMatrix *),
        LU_Factor(SparseMatrix *, int *),
        LU_Solve( SparseMatrix *, double *, double *);

/* Global definitions */
#define             CS    1.0
#define             GS    0.091
#define             GA    14.286
#define             CM    1.0
#define             GM    0.091
#define         OUTPUT    "NewRes40.dat"
#define           INFO    "Info20.dat"
#define           TEND    10.0
#define           NSIM    2000      /* Simulations to be done */
#define             NT    1000
#define             DT    1.0
#define          NODES    40
#define          NSEED    7    /* Seed for random number generator */
#define          FSEED    "GenRan40.ran"  /* History of random number generator */

/* Parameters for exact solution */
#define           NCON    75  /* Number of contacts */
#define            AMP    2.0e-5
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
// #define             RS    0.002
// #define             RD    0.000354487542
// #define              L    0.0667236981
// #define              X    0.010

/* Parameters for TestCell4.d3 (large branched dendrite) */
 #define             RS    0.002
 #define             RD    0.000648741708
 #define              L    0.2256605
 #define              X    0.0135280571

/* Global Variables */
SparseMatrix lhs, rhs;
unsigned int ix, iy, iz;

int main( int argc, char **argv )
{
    extern unsigned int ix, iy, iz;
    int k, j, id, start, begin, nodes, n, nc, i, in, nstep,
        maxstep, ncon, FirstNode, NumberOfInput;
    int counter, nb, nsim, connected;
    double *v, *x, max, *eta, *eval, *cval, AreaOfSoma, gama, *chi,
           xold, xnew, fac, arg, sum, tmp, vs, pi, dt, tnow, len,
           h, CableDiameter, ElectrotonicLength, *amp, *loc, input,
           CableLength, dx, CellLength, LocusContact;
    void srand( unsigned int);
    neuron *cell;
    contact *newcon, *oldcon, *con;
    extern SparseMatrix lhs, rhs;
    branch *bnow, *bold, *bnew, *FirstBranch, *CellFirstBranch;
    char word[20];
    FILE *fp;

/*  Initialise simulation counter */
    nsim = 1;
    maxstep = ((int) 1000.0*T);
    if ( (fp=fopen(FSEED,"r"))!=NULL ) {
        while ( fscanf(fp,"%lu %lu %lu", &ix, &iy, &iz )!=EOF ) nsim++;
        fclose(fp);
    } else {
        srand( ((unsigned int) NSEED) );
        ix = rand( );
        iy = rand( );
        iz = rand( );
    }

/*  Load Test Neuron */
    pi = 4.0*atan(1.0);
    if ( argc != 2 ) {
        printf("\n Invoke program with load <input>\n");
        return 1;
    } else {
        printf("\nOpening file %s\n",argv[1]);
        if ( (fp=fopen(argv[1],"r")) == NULL ) {
            printf("\n Test Neuron file not found");
            return 1;
        }

/*  Get branch data */
        bold = NULL;
        while  ( fscanf(fp,"%s",word) != EOF ) {
            if ( strcmp(word,"Branch") == 0 || strcmp(word,"branch") == 0 ) {
                bnew = (branch *) malloc( sizeof(branch) );
                fscanf(fp,"%d", &(bnew->id) );
                bnew->peer = NULL;
                bnew->child = NULL;
                bnew->conlist = NULL;
                if ( bold != NULL) {
                    bold->child = bnew;
                } else {
                    CellFirstBranch = bnew;
                }
                bnew->parent = bold;
                fscanf(fp,"%lf %lf %lf", &(bnew->xl), &(bnew->yl), &(bnew->zl) );
                fscanf(fp,"%lf %lf %lf", &(bnew->xr), &(bnew->yr), &(bnew->zr) );
                fscanf(fp,"%lf %lf", &(bnew->plen), &(bnew->diam) );
                bnew->elen = (bnew->plen)/sqrt(bnew->diam);
                bold = bnew;
            } else if ( strcmp(word,"Marker") == 0 || strcmp(word,"marker") == 0 ) {
                if ( bold == NULL ) {
                    printf("\nMarker is not assigned to a branch\n");
                    return 0;
                }
                printf("Found and initialised a branch contact\n");
                newcon = (contact *) malloc( sizeof(contact) );
                newcon->next = NULL;
                fscanf(fp,"%lf %lf", &newcon->xp, &newcon->amp );
                if ( bnew->conlist == NULL ) {
                    bnew->conlist = newcon;
                } else {
                    con = bnew->conlist;
                    while ( con->next ) con = con->next;
                    con->next = newcon;
                }
            } else {
                printf("Unrecognised dendritic feature\n");
                return 0;
            }
        }
        fclose(fp);
    }

/*  Compute total length of dendrite */
    CellLength = 0.0;
    bnew = CellFirstBranch;
    while ( bnew ) {
       CellLength += bnew->plen;
       bnew = bnew->child;
    }

/*  Start simulation procedure */
    start = 1;
    while ( nsim <= NSIM ) {
        printf("\r Simulation %d", nsim);

/*  Step 1. - Generate a copy of the branch list */
        bnow = CellFirstBranch;
        bold = NULL;
        while ( bnow ) {
            bnew = (branch *) malloc( sizeof(branch) );
            bnew->id = bnow->id;
            bnew->xl = bnow->xl;
            bnew->yl = bnow->yl;
            bnew->zl = bnow->zl;
            bnew->xr = bnow->xr;
            bnew->yr = bnow->yr;
            bnew->zr = bnow->zr;
            bnew->diam = bnow->diam;
            bnew->plen = bnow->plen;
            bnew->elen = bnow->elen;
            bnew->peer = NULL;
            bnew->child = NULL;
            if ( bold ) {
                bold->child = bnew;
            } else {
                FirstBranch = bnew;
            }
            bnew->parent = bold;
            bold = bnew;

            if ( bnow->conlist ) {
                oldcon = NULL;
                con = bnow->conlist;
                while ( con ) {
                    newcon = (contact *) malloc( sizeof(contact) );
                    newcon->next = NULL;
                    newcon->fl = NULL;
                    newcon->fr = NULL;
                    newcon->amp = con->amp;
                    newcon->id = con->id;
                    newcon->xp = con->xp;
                    newcon->xl = NULL;
                    newcon->xr = NULL;
                    if ( oldcon != NULL ) {
                        oldcon->next = newcon;
                    } else {
                        bnew->conlist = newcon;
                    }
                    oldcon = newcon;
                    con = con->next;
                }
            } else {
                bnew->conlist = NULL;
            }
            bnow = bnow->child;
        }

/*  STEP 1. - Randomly place NCON inputs on branches */
        for ( k=0 ; k<NCON ; k++ ) {
            LocusContact = CellLength*ran( &ix, &iy, &iz);
            bnew = FirstBranch;
            len = bnew->plen;
            while ( LocusContact > len ) {
                bnew = bnew->child;
                len += bnew->plen;
            }
            newcon = (contact *) malloc( sizeof(contact) );
            newcon->next = NULL;
            newcon->fl = NULL;
            newcon->fr = NULL;
            newcon->amp = AMP;
            newcon->xl = NULL;
            newcon->xr = NULL;
            newcon->xp = LocusContact-(len-bnew->plen);
            if ( bnew->conlist ) {
                oldcon = bnew->conlist;
                while ( oldcon->next ) oldcon = oldcon->next;
                oldcon->next = newcon;
            } else {
                bnew->conlist = newcon;
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
//        printf("\nTree contains %d individual dendrite(s) ...\n", n);
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
            } while ( bnew );
            if ( !connected ) cell->dendlist[n++].root = bold;
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

/*  STEP 6A. - Identify scaled soma-to-tip electrotonic length */
        ElectrotonicLength = 0.0;
        bnow = cell->dendlist[0].root;
        while ( bnow != NULL ) {
            ElectrotonicLength += bnow->elen;
            bnow = bnow->child;
        }

/*  STEP 6B. - Identify diameter of equivalent cable */
        CableDiameter = 0.0;
        for ( k=0 ; k<cell->ndend ; k++ ) {
            CableDiameter += pow(cell->dendlist[k].root->diam,1.5);
        }
        CableDiameter = pow(CableDiameter,2.0/3.0);

/*  STEP 6C. - Compute length of equivalent cable */
        CableLength = ElectrotonicLength*sqrt(CableDiameter);

/*  STEP 7A. - Count number on inputs on Cell */
        NumberOfInput = 0;
        for ( k=0 ; k<cell->ndend ; k++ ) {
            bnow = cell->dendlist[k].root;
            NumberOfInput += Count_Contacts( cell->dendlist[k].root, bnow);
        }
        amp = (double *) malloc( NumberOfInput*sizeof(double) );
        loc = (double *) malloc( NumberOfInput*sizeof(double) );

/*  STEP 7B. - Identify position of contacts on Equivalent Cable
               and thence their relative position the true cable */
        for ( ncon=k=0 ; k<cell->ndend ; k++ ) {
            bnow = cell->dendlist[k].root;
            Find_Contacts( bnow, loc, amp, &ncon);
        }
        tmp = sqrt(CableDiameter)/CableLength;
        for ( k=0 ; k<ncon ; k++ ) loc[k] *= tmp;
        if ( ncon == 0 ) {
            printf("\n No contact found - Fatal error!");
            return 1;
        }

/*  STEP 7C. - Activate to output branch properties */
//        fp = fopen(INFO,"w");
//        for ( k=0 ; k<cell->ndend ; k++ ) {
//            bnow = cell->dendlist[k].root;
//            Output_Info( bnow, fp);
//        }
//        fclose(fp);

/*  STEP 8A. - Construct eigenvalues for exact solution */
        if ( start ) {
            AreaOfSoma = 4.0*pi*RS*RS;
            gama = AreaOfSoma/(pi*CableDiameter*CableLength);
            eval = (double *) malloc( (M+1)*sizeof(double) );
            cval = (double *) malloc( (M+1)*sizeof(double) );
            eval[0] = 0.0;
            cval[0] = 1.0;
            for ( k=1 ; k<=M ; k++ ) {
                xnew = arg = pi*((double) k );
                do {
                    xold = xnew;
                    xnew = arg-atan(gama*xold);
                } while ( fabs(xold-xnew) > 5.e-7 );
                eval[k] = xnew;
                cval[k] = cos(xnew);
            }

/*  STEP 8B. - Construct time constants */
            eta = (double *) malloc( (M+1)*sizeof(double) );
            eta[0] = GM/CM;
            for ( k=1 ; k<=M ; k++ ) {
                eta[k] = (GM+0.25*CableDiameter*GA*pow(eval[k]/CableLength,2))/CM;
            }
        }
        chi = (double *) malloc( (M+1)*sizeof(double) );
        fac = pi*CM*CableDiameter*CableLength;
        for ( k=0 ; k<=M ; k++ ) {
            chi[k] = SIN*cval[k];
            for ( j=0 ; j<ncon ; j++ ) {
                chi[k] += amp[j]*cos(eval[k]*(1.0-loc[j]));
            }
            chi[k] *= cval[k];
            chi[k] /= fac*eta[k]*(1.0+gama*cval[k]*cval[k]);
        }
        for ( k=1 ; k<=M ; k++ ) chi[k] *= 2.0;

/*  STEP 9. - Count dendritic branches */
        for ( nb=k=0 ; k<cell->ndend ; k++ ) {
            bnow = cell->dendlist[k].root;
            nb += Count_Branches( bnow, bnow);
        }
        h = CellLength/((double) NODES-nb);
        for ( k=0 ; k<cell->ndend ; k++ ) Assign_Branch_Nodes( cell->dendlist[k].root, &h);

/*  STEP 10. - Enumerate Nodes */
        FirstNode = 0;
        for ( k=0 ; k<cell->ndend ; k++ ) Enumerate_Nodes( cell->dendlist[k].root, &FirstNode );
        for ( k=0 ; k<cell->ndend ; k++ ) cell->dendlist[k].root->junct = FirstNode;
//        printf("Number of nodes is %d\n", FirstNode+1);
//        getchar( );

/*  STEP 11. - Construct Sparse Matrices */
        nodes = FirstNode+1;
        if ( start ) {
            Matrix_Malloc( &lhs, nodes, 3*nodes-2 );
            Matrix_Malloc( &rhs, nodes, 3*nodes-2 );
        }
        lhs.StartRow[0] = rhs.StartRow[0] = 0;
        for ( counter=k=0 ; k<cell->ndend ; k++ ) {
            bnow = cell->dendlist[k].root;
            Generate_Dendrite( bnow, &counter);
        }
        lhs.a[3*nodes-3] = rhs.a[3*nodes-3] = 0.0;
        for ( k=0 ; k<cell->ndend ; k++ ) {
            bnow = cell->dendlist[k].root;
            lhs.a[counter] = (bnow->diam)*(bnow->hseg);
            rhs.a[counter] = -pow(bnow->diam,2)/(bnow->hseg);
            lhs.col[counter] = rhs.col[counter] = bnow->first;
            lhs.a[3*nodes-3] += 2.0*(bnow->diam)*(bnow->hseg);
            rhs.a[3*nodes-3] += pow(bnow->diam,2)/(bnow->hseg);
            counter++;
        }
        if ( counter != 3*(nodes-1) ) {
            printf("\nEnumeration error!\n");
            getchar( );
        } else {
//            printf("Generated numerical representation of dendrite\n");
//            getchar( );
        }
        lhs.col[counter] = rhs.col[counter] = nodes-1;
        lhs.StartRow[nodes] = rhs.StartRow[nodes] = counter+1;

/*  STEP 12. - Construct current input */
        for( k=0 ; k<cell->ndend ; k++ ) Input_Current(cell->dendlist[k].root);
        if ( nsim == 1 ) {
            fp = fopen(OUTPUT, "w");
            fclose(fp);
        }
//        printf("Constructed sparse matrices\n");
//        getchar( );
        dt = 1.0/((double) NT);
        for ( k=0 ; k<3*nodes-2 ; k++ ) {
            rhs.a[k] = pi*dt*(0.125*GA*rhs.a[k]+GM*lhs.a[k]/12.0);
            lhs.a[k] *= pi*CM/6.0;
            lhs.a[k] += rhs.a[k];
            rhs.a[k] = lhs.a[k]-2.0*rhs.a[k];
        }

/*  Add current from the soma */
        lhs.a[3*nodes-3] += AreaOfSoma*(CS+0.5*GS*dt);
        rhs.a[3*nodes-3] += AreaOfSoma*(CS-0.5*GS*dt);

        if ( start ) {
            v = (double *) malloc( (nodes)*sizeof(double) );
            x = (double *) malloc( (nodes)*sizeof(double) );
        }
        for ( k=0 ; k<nodes ; k++ ) v[k] = x[k] = 0.0;
        begin = 1;
        LU_Factor(&lhs, &begin);

/*  Initialise temporal integration */
        nstep = 0;
        while ( nstep < maxstep ) {
            Matrix_Vector_Multiply( &rhs, v, x);
            x[nodes-1] -= dt*SIN;
            for ( k=0 ; k<nodes ; k++ ) v[k] = x[k];
            if ( nstep < maxstep ) {
                for ( k=0 ; k<cell->ndend ; k++ ) Assign_Current(cell->dendlist[k].root, x, dt);
            }
            LU_Solve( &lhs, v, x );
            nstep++;
            if ( nstep%1000 == 0 ) {
//                printf("\rReached time %5.1lf ms\t", tout);
                tnow = dt*((double) nstep);
                for ( vs=0.0,k=M ; k>=0 ; k-- ) {
                    arg = tnow*eta[k];
                    if ( arg > 20.0 ) {
                        tmp = 1.0;
                    } else {
                        tmp = (1.0-exp(-arg));
                    }
                    vs -= tmp*chi[k];
                }
//                printf("\nNumerical Voltage %12.6lf mV\n",v[nodes-1]);
//                printf("Exact Voltage %12.6lf mV\n", vs);
                fp = fopen(OUTPUT,"a");
                fprintf(fp,"%20.15lf",v[nodes-1]);
                fclose(fp);
            }
        }
        fp = fopen(OUTPUT, "a");
        fprintf(fp,"\n");
        fclose(fp);

        free(amp);
        free(loc);
        free(chi);
        Destroy_Test_Neuron( cell );

/*  Update seed file */
        if ( nsim == 1 ) {
            fp = fopen(FSEED,"w");
        } else {
            fp = fopen(FSEED,"a");
        }
        fprintf(fp,"%u \t %u \t %u\n", ix, iy, iz);
        fclose(fp);
        if ( start ) start = 0;
        nsim++;
    }
    return 0;
}


void Output_Info( branch *b, FILE *fp)
{
    if ( b->child ) Output_Info( b->child, fp);
    if ( b->peer ) Output_Info( b->peer, fp);

    fprintf(fp,"%3d \t $3d\n", b->id, b->nc);
    return;
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


 /************************************************
         Function to destroy a TEST NEURON
  ************************************************/
void Destroy_Test_Neuron(neuron *cell)
{
    int k;

    for ( k=0 ; k<cell->ndend ; k++ ) {
        Destroy_Test_Dendrite( cell->dendlist[k].root );
    }
    free(cell);
    return;
}


 /***************************************************
              Function to destroy TEST DENDRITE
  ***************************************************/
void Destroy_Test_Dendrite( branch *b )
{
    int i;
    contact *prevcon, *nextcon;

    if ( b->child ) Destroy_Test_Dendrite(b->child);
    if ( b->peer ) Destroy_Test_Dendrite(b->peer);
    if ( b->conlist ) {
        prevcon = b->conlist;
        do {
            nextcon = prevcon->next;
            free(prevcon);
            prevcon = nextcon;
        } while ( prevcon );
    }
    free(b);
    return;
}


 /*********************************************
      Function to count contacts on branches
  *********************************************/
int Count_Contacts( branch *bstart, branch *bnow)
{
    static int n;
    contact *con;

    if ( bstart == bnow ) n = 0;
    if ( bnow != NULL ) {
        if ( bnow->child ) Count_Contacts(bstart, bnow->child);
        if ( bnow->peer ) Count_Contacts(bstart, bnow->peer);
        con = bnow->conlist;
        while ( con ) {
            n++;
            con = con->next;
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


 /***************************************************
       Function to find contacts on a dendrite
  ***************************************************/
void Find_Contacts( branch *b, double *loc, double *amp, int *ncon)
{
    contact *con;
    branch *btmp;

    if ( b->child ) Find_Contacts( b->child, loc, amp, ncon);
    if ( b->peer ) Find_Contacts( b->peer, loc, amp, ncon);

    con = b->conlist;
    while ( con ) {
        amp[(*ncon)] = con->amp;
        loc[(*ncon)] = (con->xp)/sqrt(b->diam);
        btmp = b->parent;
        while ( btmp ) {
            loc[(*ncon)] += btmp->elen;
            btmp = btmp->parent;
        }
        (*ncon)++;
        con = con->next;
    }
    return;
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
            lhs.a[*counter] = (btmp->diam)*(btmp->hseg);
            rhs.a[*counter] = -pow(btmp->diam,2)/(btmp->hseg);
            lhs.col[*counter] = rhs.col[*counter] = btmp->first;
            SumL += 2.0*(btmp->diam)*(btmp->hseg);
            SumR += pow(btmp->diam,2)/(btmp->hseg);
            (*counter)++;
            btmp = btmp->peer;
        }
        lhs.a[*counter] = SumL+2.0*(b->diam)*(b->hseg);
        rhs.a[*counter] = SumR+pow(b->diam,2)/(b->hseg);
        lhs.col[*counter] = rhs.col[*counter] = CurrentNode;
        (*counter)++;
        lhs.a[*counter] = (b->diam)*(b->hseg);
        rhs.a[*counter] = -pow(b->diam,2)/(b->hseg);
        if ( CurrentNode == b->first ) {
            lhs.col[*counter] = rhs.col[*counter] = b->junct;
        } else {
            lhs.col[*counter] = rhs.col[*counter] = CurrentNode+1;
        }
        (*counter)++;
        lhs.StartRow[CurrentNode+1] = rhs.StartRow[CurrentNode+1] = *counter;
    } else {
        lhs.a[*counter] = 2.0*(b->diam)*(b->hseg);
        rhs.a[*counter] = pow(b->diam,2)/(b->hseg);
        lhs.col[*counter] = rhs.col[*counter] = CurrentNode;
        (*counter)++;
        lhs.a[*counter] = (b->diam)*(b->hseg);
        rhs.a[*counter] = -pow(b->diam,2)/(b->hseg);
        if ( CurrentNode == b->first ) {
            lhs.col[*counter] = rhs.col[*counter] = b->junct;
        } else {
            lhs.col[*counter] = rhs.col[*counter] = CurrentNode+1;
        }
        (*counter)++;
        lhs.StartRow[CurrentNode+1] = rhs.StartRow[CurrentNode+1] = *counter;
    }

/* Step 3 - Build matrix entries for internal node of branch */
    for ( k=nc-1 ; k>0 ; k-- ) {
        CurrentNode++;
        lhs.a[*counter] = (b->diam)*(b->hseg);
        rhs.a[*counter] = -pow(b->diam,2)/(b->hseg);
        lhs.col[*counter] = rhs.col[*counter] = CurrentNode-1;
        (*counter)++;
        lhs.a[*counter] = 4.0*(b->diam)*(b->hseg);
        rhs.a[*counter] = 2.0*pow(b->diam,2)/(b->hseg);
        lhs.col[*counter] = rhs.col[*counter] = CurrentNode;
        (*counter)++;
        lhs.a[*counter] = (b->diam)*(b->hseg);
        rhs.a[*counter] = -pow(b->diam,2)/(b->hseg);
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
       Function to input current to dendrite
  ***********************************************/
void Input_Current( branch *b )
{
    int k;
    double len;
    contact *con;

    if ( b->child ) Input_Current( b->child );
    if ( b->peer ) Input_Current( b->peer );

    con = b->conlist;
    while ( con ) {
        k = 0;
        len = b->hseg;
        while ( con->xp > len ) {
            k++;
            len += b->hseg;
        }
        if ( k == 0 ) {
            con->xl = b->junct;
            con->xr = b->first;
            con->fl = 1.0-(con->xp)/(b->hseg);
            con->fr = (con->xp)/(b->hseg);
        } else {
            con->xl = b->first-k+1;
            con->xr = b->first-k;
            con->fl = (len-con->xp)/(b->hseg);
            con->fr = 1.0-(con->fl);
        }
        con = con->next;
    }
    return;
}


 /****************************************
        Function to assign current
  ****************************************/
void Assign_Current(branch *bnow, double *x, double fac )
{
    contact *con;

    if ( bnow->child ) Assign_Current(bnow->child, x, fac );
    if ( bnow->peer ) Assign_Current(bnow->peer, x, fac );

    con = bnow->conlist;
    while ( con ) {
        x[con->xl] -= fac*(con->fl)*(con->amp);
        x[con->xr] -= fac*(con->fr)*(con->amp);
        con = con->next;
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


 /***********************************************
      Function To Factorise A Sparse Matrix
  ***********************************************/
void LU_Factor(SparseMatrix *m, int *start)
{
    double tmp, sum;
    int i, j, k, r, n, cl, cu, col, row;

/*  Step 1. - Identify matrix dimension */
    n = m->n;

/*  Step 2. - Fill column vectors for triangular matrices */
    if ( *start ) {
        cl = cu = 0;
        for ( i=k=0 ; i<n ; i++ ) {
            m->l->StartRow[i] = cl;
            m->u->StartRow[i] = cu;
            while ( m->col[k] < i ) m->l->col[cl++] = m->col[k++];
            m->l->col[cl++] = m->col[k];
            m->u->col[cu++] = m->col[k++];
            while ( k < m->StartRow[i+1] ) m->u->col[cu++] = m->col[k++];
        }
        m->l->StartRow[n] = cl;
        m->u->StartRow[n] = cu;
        *start = 0;
    }

/*  Step 3. - Fill remaining entries of L and U row by row */
    cl = cu = 0;
    for ( i=0 ; i<n ; i++ ) {
        for ( k=m->StartRow[i] ; k < m->StartRow[i+1] ; k++ ) {
            if ( m->col[k] < i ) {
                sum = m->a[k];
                for ( j=m->l->StartRow[i] ; j<cl ; j++ ) {
                    col = m->l->col[j];
                    row = m->u->StartRow[col];
                    while ( m->u->col[row] < m->col[k] && row < m->u->StartRow[col+1] ) row++;
                    if ( m->u->col[row] == m->col[k] ) sum -= (m->l->a[j])*(m->u->a[row]);
                }
                row = m->u->StartRow[m->l->col[cl]];
                while ( m->u->col[row] < m->col[k] ) row++;
                m->l->a[cl++] = sum/(m->u->a[row]);
            } else if ( m->col[k] == i ) {
                m->l->a[cl++] = 1.0;
                sum = m->a[k];
                for ( j=m->l->StartRow[i] ; j<m->l->StartRow[i+1]-1; j++ ) {
                    col = m->l->col[j];
                    row = m->u->StartRow[col];
                    while ( m->u->col[row] < i && row < m->u->StartRow[col+1] ) row++;
                    if ( m->u->col[row] == i ) sum -= (m->l->a[j])*(m->u->a[row]);
                }
                m->u->a[cu++] = sum;
            } else {
                sum = m->a[k];
                for ( j=m->l->StartRow[i] ; j<m->l->StartRow[i+1]-1; j++ ) {
                    col = m->l->col[j];
                    row = m->u->StartRow[col];
                    while ( m->u->col[row] < m->col[k] && row < m->u->StartRow[col+1] ) row++;
                    if ( m->u->col[row] == m->col[k] ) sum -= (m->l->a[j])*(m->u->a[row]);
                }
                m->u->a[cu++] = sum;
            }
        }
    }
    return;
}


 /**************************************************
          Function to Solve the matrix problem
  **************************************************/
void LU_Solve(SparseMatrix *m, double *x, double *b )
{
    int i,j;
    double *z;

    z = (double *) malloc( (m->n)*sizeof(double) );

    for ( i=0 ; i<m->n ; i++ ) {
        z[i] = b[i];
        for (j=m->l->StartRow[i];j<m->l->StartRow[i+1]-1;j++)
        { z[i] -= (m->l->a[j])*(z[(m->l->col[j])]); }
        z[i] /= m->l->a[m->l->StartRow[i+1]-1];
    }

    for ( i = (m->n) - 1 ; i>=0 ; i-- ) {
        x[i] = z[i];
        for (j=m->u->StartRow[i]+1;j<m->u->StartRow[i+1];j++)
        { x[i] -= (m->u->a[j])*(x[m->u->col[j]]); }
        x[i] /= m->u->a[m->u->StartRow[i]];
    }
    free(z);
    return;
}


 /**************************************************
            Function to assign branch nodes
  **************************************************/
void Assign_Branch_Nodes( branch *b, double *h )
{
    int nc, k;
    double hseg;

    nc = ((int) ceil((b->plen)/(*h)));
    b->hseg = (b->plen)/((double) nc);
    b->nc = nc;

    if ( b->child != NULL) Assign_Branch_Nodes( b->child, h);
    if ( b->peer != NULL) Assign_Branch_Nodes( b->peer, h);
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

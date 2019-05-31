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

        int xl;                     /* Left hand node */
        int xr;                     /* Right hand node */
        double frac;                /* Fraction of input to left hand node */

        struct contact_t *next;     /* Address of next contact */
} contact;


typedef struct branch_t
{
/*  Connectivity of branch */
        struct branch_t *parent;    /* Address of parent branch */
        struct branch_t *child;     /* Address of child branch */
        struct branch_t *peer;      /* Addresss of peer branch */

/*  Physical properties of branch */
        int nd;                     /* Number of nodes on branch specification */
        double xl;                  /* X-coordinate of lefthand endpoint */
        double yl;                  /* Y-coordinate of lefthand endpoint */
        double zl;                  /* Z-coordinate of lefthand endpoint */
        double xr;                  /* X-coordinate of righthand endpoint */
        double yr;                  /* Y-coordinate of righthand endpoint */
        double zr;                  /* Z-coordinate of righthand endpoint */
        double diam;                /* Branch diameter (microns) */
        double plen;                /* Branch length (microns) */
        double elen;                /* Branch length (eu) */

        double *d;                  /* Diameter of dendrite (micron) */
        double *x;                  /* Location of nodes on dendrite (micron) */

/*  Node information for spatial representation */
        int nodes;                  /* Total number nodes spanning branch */
        int junct;                  /* Junction node of the branch */
        int first;                  /* Internal node connected to junction */

/*  Contact information */
        contact *conlist;               /* Branch contact */
} branch;


typedef struct dendrite_t
{
        branch *root;               /* Pointer to root branch of dendrite */
        double plen;                /*  length of dendrite */
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
#define         OUTPUT    "contactinfo.dat"
#define           TEND    10.0
#define           NSIM    1      /* Simulations to be done */
#define             NT    1000
#define             DT    1.0
#define          NODES    100
#define          NSEED    7    /* Seed for random number generator */
#define          FSEED    "GenRan.ran"  /* History of random number generator */

/* Parameters for exact solution */
#define           NCON    1  /* Number of contacts */
#define            AMP    0.0e-3
#define            SIN    1.0e-3
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
    int k, j, id, start, begin, nodes, n, nc, i, in,
        ncon, FirstNode, NumberOfInput;
    int counter, nb, nsim, connected;
    double *v, *x, max, *eta, *eval, *cval, AreaOfSoma, gama,
           *chi, xold, xnew, frac, arg, sum, tmp, vs, pi, dt,
           tnow, tout, len, h, CableDiameter, ElectrotonicLength,
           *amp, *loc, input, CableLength, dx, CellLength, LocusContact;
    void srand( unsigned int);
    neuron *cell;
    contact *newcon, *oldcon, *con;
    extern SparseMatrix lhs, rhs;
    branch *bnow, *bold, *bnew, *FirstBranch, *CellFirstBranch;
    char word[20];
    FILE *fp;

/*  Initialise simulation counter */
    nsim = 1;
    start = 1;

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
    if ( argc != 2 ) {
        printf("\n Invoke program with load <input>\n");
        return 1;
    } else {
        printf("\nOpening file %s\n",argv[1]);
        if ( (fp=fopen(argv[1],"r")) == NULL ) {
            printf("\n Test Neuron file does not found");
            return 1;
        }
        pi = 4.0*atan(1.0);

/*  Get branch data */
        bold = NULL;
        while  ( fscanf(fp,"%s",word) != EOF ) {
            if ( strcmp(word,"Branch") == 0 || strcmp(word,"branch") == 0 ) {
                fscanf(fp, "%d", &nodes);
                printf("Found a branch defined by %d nodes\n", nodes);
                bnew = (branch *) malloc( sizeof(branch) );
                bnew->peer = NULL;
                bnew->child = NULL;
                bnew->d = NULL;
                bnew->x = NULL;
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
    while ( nsim <= NSIM ) {

/*  Step 1. - Generate a copy of the branch list */
        bnow = CellFirstBranch;
        bold = NULL;
        while ( bnow != NULL ) {
            bnew = (branch *) malloc( sizeof(branch) );
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
            bnew->d = NULL;
            bnew->x = NULL;
            if ( bold != NULL ) {
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
                    newcon->frac = NULL;
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
            newcon->frac = NULL;
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
        printf("\nTree contains %d individual dendrite(s) ...\n", n);
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

/*  STEP 7B. - Identify position of contacts on Equivalent Cable */
        for ( ncon=k=0 ; k<cell->ndend ; k++ ) {
            bnow = cell->dendlist[k].root;
            Find_Contacts( bnow, loc, amp, &ncon);
        }
        tmp = sqrt(CableDiameter);
        for ( k=0 ; k<ncon ; k++ ) loc[k] *= tmp;
        if ( ncon == 0 ) {
            printf("\n No contact found - Fatal error!");
            return 1;
        }

/*  STEP 8. - Construct exact solution */
        if ( start ) {
            frac = loc[0]/CableLength;
            AreaOfSoma = 4.0*pi*RS*RS;
            gama = AreaOfSoma/(pi*CableDiameter*CableLength);
            eta = (double *) malloc( (M+1)*sizeof(double) );
            chi = (double *) malloc( (M+1)*sizeof(double) );
            eval = (double *) malloc( (M+1)*sizeof(double) );
            cval = (double *) malloc( (M+1)*sizeof(double) );
            eta[0] = GM/CM;
            chi[0] = (AMP+SIN)/(pi*CM*CableDiameter*CableLength*(1.0+gama));
            chi[0] /= eta[0];
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
        }
        for ( k=1 ; k<=M ; k++ ) {
            eta[k] = (GM+0.25*CableDiameter*GA*pow(eval[k]/CableLength,2) )/CM;
            chi[k] = 2.0*cval[k]*(AMP*cos(eval[k]*(1.0-frac))+SIN*cval[k])/
                         (pi*CM*CableDiameter*CableLength*(1.0+gama*cval[k]*cval[k]));
            chi[k] /= eta[k];
        }

/*  STEP 9. - Count dendritic segments */
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
        printf("Number of nodes is %d\n", FirstNode+1);
        getchar( );

/*  STEP 11. - Construct Sparse Matrices */
        nodes = FirstNode+1;
        if ( start ) {
            Matrix_Malloc( &lhs, nodes, 3*nodes-2 );
            Matrix_Malloc( &rhs, nodes, 3*nodes-2 );
        }
        lhs.n = rhs.n = nodes;
        lhs.a[3*nodes-3] = rhs.a[3*nodes-3] = 0.0;
        lhs.StartRow[0] = rhs.StartRow[0] = 0;

        for ( counter=k=0 ; k<cell->ndend ; k++ ) {
            bnow = cell->dendlist[k].root;
            Generate_Dendrite( bnow, &counter);
        }
        printf("Generated numerical representation of dendrite\n");
        getchar( );

        for ( k=0 ; k<cell->ndend ; k++ ) {
            bnow = cell->dendlist[k].root;
            lhs.a[counter] = 0.5*bnow->d[1]*bnow->x[1];
            rhs.a[counter] = -(bnow->d[0])*(bnow->d[1])/bnow->x[1];
            lhs.col[counter] = rhs.col[counter] = bnow->first;
            lhs.a[3*nodes-3] += 1.5*(bnow->d[0])*bnow->x[1];
            rhs.a[3*nodes-3] += (bnow->d[0])*(bnow->d[1])/bnow->x[1];
            counter++;
        }
        lhs.col[3*nodes-3] = rhs.col[3*nodes-3] = nodes-1;
        lhs.StartRow[nodes] = rhs.StartRow[nodes] = 3*nodes-2;

        for( k=0 ; k<cell->ndend ; k++ ) Input_Current(cell->dendlist[k].root);

        fp = fopen("out.res", "w");
        fclose(fp);
        printf("Constructed sparse matricesn");
        getchar( );

        dt = 1.0/((double) NT);
        for ( k=0 ; k<3*nodes-2 ; k++ ) {
            rhs.a[k] *= GA;
            rhs.a[k] += GM*lhs.a[k];
            lhs.a[k] *= CM;
            rhs.a[k] *= 0.5*dt;
            lhs.a[k] += rhs.a[k];
            rhs.a[k] = lhs.a[k]-2.0*rhs.a[k];
        }

/*  Add capacitive term of soma */
        lhs.a[3*nodes-3] += (4.0*AreaOfSoma/pi)*(CS+0.5*GS*dt);
        rhs.a[3*nodes-3] += (4.0*AreaOfSoma/pi)*(CS-0.5*GS*dt);

        if ( start ) {
            v = (double *) malloc( (nodes)*sizeof(double) );
            x = (double *) malloc( (nodes)*sizeof(double) );
        }
        for ( k=0 ; k<nodes ; k++ ) v[k] = x[k] = 0.0;
        begin = 1;
        LU_Factor(&lhs, &begin);

/*  Initialise temporal integration */
        tnow = 0.0;
        tout = DT;
        while ( tnow < TEND ) {
            tnow += dt;
            Matrix_Vector_Multiply(&rhs,v,x);
            x[nodes-1] -= 4.0*dt*SIN/pi;
            tmp = 4.0*dt/pi;
            for ( k=0 ; k<nodes ; k++ ) v[k] = x[k];
            if ( tnow < T ) {
                for ( k=0 ; k<cell->ndend ; k++ ) Assign_Current(cell->dendlist[k].root, x, tmp);
            }
            LU_Solve( &lhs, v, x );

            if ( tnow > tout && tnow <= T ) {
                printf("\rReached time %5.1lf ms\t", tout);
                for ( vs=0.0,k=M ; k>=0 ; k-- ) {
                    arg = tnow*eta[k];
                    if ( arg > 20.0 ) {
                        tmp = 1.0;
                    } else {
                        tmp = (1.0-exp(-arg));
                    }
                    vs -= tmp*chi[k];
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
                for ( vs=0.0,k=M ; k>=0 ; k-- ) {
                    arg = (tnow-T)*eta[k];
                    if ( arg < 20.0 ) {
                        tmp = exp(-arg);
                        arg = tnow*eta[k];
                        if ( arg < 20.0 ) tmp -= exp(-arg);
                        vs -= tmp*chi[k];
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

        free(amp);
        free(loc);
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
    free( b->d );
    free( b->x );
    if ( b->conlist ) {
        prevcon = b->conlist;
        do {
            nextcon = prevcon->next;
            free(prevcon);
            prevcon = nextcon;
        } while ( prevcon );
    }
    free (b);
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
    if ( con ) {
        amp[(*ncon)] = con->amp;
        loc[(*ncon)] = (con->xp)/sqrt(b->diam);
        btmp = b->parent;
        while ( btmp ) {
            loc[(*ncon)] += btmp->elen;
            btmp = btmp->parent;
        }
        (*ncon)++;
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
    *FirstNode += (bnow->nd)-1;
    bnow->first = *FirstNode-1;
    return;
}


 /***************************************************
        Function to constuct sparse matrices
  ***************************************************/
void Generate_Dendrite( branch *bnow, int *counter)
{
    int i, k, n;
    extern SparseMatrix lhs, rhs;
    branch *btmp;
    double SumL, SumR;

/* Step 1 - Recurse to the end of the dendrite */
    if ( bnow->child ) Generate_Dendrite( bnow->child, counter);
    if ( bnow->peer ) Generate_Dendrite( bnow->peer, counter);

/* Step 2 - Build matrix entries for distal node of branch */
    n = (bnow->nd)-1;
    k = (bnow->first)-(bnow->nd)+2;
    if ( bnow->child ) {
        btmp = bnow->child;
        SumR = SumL = 0.0;
        while ( btmp ) {
            lhs.a[*counter] = 0.5*btmp->d[1]*btmp->x[1];
            rhs.a[*counter] = -(btmp->d[0]*btmp->d[1])/btmp->x[1];
            lhs.col[*counter] = rhs.col[*counter] = btmp->first;
            SumL += 1.5*btmp->d[0]*btmp->x[1];
            SumR += (btmp->d[0]*btmp->d[1])/btmp->x[1];
            (*counter)++;
            btmp = btmp->peer;
        }
        lhs.a[*counter] = SumL+1.5*bnow->d[n]*(bnow->x[n]-bnow->x[n-1]);
        rhs.a[*counter] = SumR+(bnow->d[n-1]*bnow->d[n])/(bnow->x[n]-bnow->x[n-1]);
        lhs.col[*counter] = rhs.col[*counter] = k;
        (*counter)++;
        lhs.a[*counter] = 0.5*(bnow->d[n-1])*(bnow->x[n]-bnow->x[n-1]);
        rhs.a[*counter] = -(bnow->d[n-1]*bnow->d[n])/(bnow->x[n]-bnow->x[n-1]);
        if ( k == bnow->first ) {
            lhs.col[*counter] = rhs.col[*counter] = bnow->junct;
        } else {
            lhs.col[*counter] = rhs.col[*counter] = k+1;
        }
        (*counter)++;
        lhs.StartRow[k+1] = rhs.StartRow[k+1] = *counter;
    } else {
        lhs.a[*counter] = 1.5*(bnow->d[n])*(bnow->x[n]-bnow->x[n-1]);
        rhs.a[*counter] = (bnow->d[n-1]*bnow->d[n])/(bnow->x[n]-bnow->x[n-1]);
        lhs.col[*counter] = rhs.col[*counter] = k;
        (*counter)++;
        lhs.a[*counter] = 0.5*(bnow->d[n-1])*(bnow->x[n]-bnow->x[n-1]);
        rhs.a[*counter] = -(bnow->d[n-1]*bnow->d[n])/(bnow->x[n]-bnow->x[n-1]);
        if ( k == bnow->first ) {
            lhs.col[*counter] = rhs.col[*counter] = bnow->junct;
        } else {
            lhs.col[*counter] = rhs.col[*counter] = k+1;
        }
        (*counter)++;
        lhs.StartRow[k+1] = rhs.StartRow[k+1] = *counter;
    }

/* Step 3 - Build matrix entries for internal node of branch */
    for ( i=(bnow->nd)-2 ; i>0 ; i-- ) {
        k = bnow->first+1-i;
        lhs.a[*counter] = 0.5*(bnow->d[i+1])*(bnow->x[i+1]-bnow->x[i]);
        rhs.a[*counter] = -(bnow->d[i]*bnow->d[i+1])/(bnow->x[i+1]-bnow->x[i]);
        lhs.col[*counter] = rhs.col[*counter] = k-1;
        (*counter)++;
        lhs.a[*counter] = 1.5*(bnow->d[i])*(bnow->x[i+1]-bnow->x[i-1]);
        rhs.a[*counter] = (bnow->d[i-1]*bnow->d[i])/(bnow->x[i]-bnow->x[i-1])
                          +(bnow->d[i]*bnow->d[i+1])/(bnow->x[i+1]-bnow->x[i]);
        lhs.col[*counter] = rhs.col[*counter] = k;
        (*counter)++;
        lhs.a[*counter] = 0.5*(bnow->d[i-1])*(bnow->x[i]-bnow->x[i-1]);
        rhs.a[*counter] = -(bnow->d[i-1]*bnow->d[i])/(bnow->x[i]-bnow->x[i-1]);
        lhs.col[*counter] = rhs.col[*counter] = k+1;
        if ( k == bnow->first ) {
            lhs.col[*counter] = rhs.col[*counter] = bnow->junct;
        } else {
            lhs.col[*counter] = rhs.col[*counter] = k+1;
        }
        (*counter)++;
        lhs.StartRow[k+1] = rhs.StartRow[k+1] = *counter;
    }
    return;
}


 /***********************************************
       Function to input current to dendrite
  ***********************************************/
void Input_Current( branch *bnow )
{
    int k;
    double tmp;
    contact *con;

    if ( bnow->child ) Input_Current( bnow->child );
    if ( bnow->peer ) Input_Current( bnow->peer );

    con = bnow->conlist;
    while ( con ) {
        if ( con->xp <= bnow->x[1] ) {
            con->xl = bnow->junct;
            con->xr = bnow->first;
            con->frac = 1.0-(con->xp)/(bnow->x[1]);
        } else {
            k = 1;
            while ( con->xp > bnow->x[k] ) k++;
            con->xl = bnow->first-k+2;
            con->xr = bnow->first-k+1;
            con->frac = (bnow->x[k]-con->xp)/(bnow->x[k]-bnow->x[k-1]);
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
        x[con->xl] -= fac*(con->frac)*(con->amp);
        x[con->xr] -= fac*(1.0-(con->frac))*(con->amp);
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
    int n, k;
    double hused;

    n = ((int) ceil((b->plen)/(*h)))+1;
    hused = (b->plen)/((double) n-1);
    b->d = (double *) malloc( n*sizeof(double));
    b->x = (double *) malloc( n*sizeof(double));
    b->x[0] = 0.0;
    for ( k=1 ; k<(n-1) ; k++ ) b->x[k] = ((double) k)*hused;
    b->x[n-1] = b->plen;
    for ( k=0 ; k<n ; k++ ) b->d[k] = b->diam;
    b->nd = n;

    if ( b->child != NULL) Assign_Branch_Nodes( b->child, h);
    if ( b->peer != NULL) Assign_Branch_Nodes( b->peer, h);
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

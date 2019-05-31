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
        int node;                   /* Node of contact */

        struct contact_t *next;     /* Address of next contact */
} contact;


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
        int nc;                     /* Number of compartments specifying branch */
        int junct;                  /* Junction node of the branch */
        int first;                  /* Internal node connected to junction */

/*  Contact information */
        contact *conlist;           /* Branch contact */
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
        Count_Contacts( branch *, branch *);

double  branch_length( branch *, branch *),
        ran(unsigned int *, unsigned int *, unsigned int *);

void    Destroy_Test_Neuron(neuron *),
        Remove_Branch( branch **, branch *),
        Destroy_Test_Dendrite( branch *),
        Assign_Branch_Nodes( branch *, double *),
        Enumerate_Nodes( branch *, int *),
        Matrix_Vector_Multiply( SparseMatrix *, double *, double *),
        Matrix_Malloc( SparseMatrix *, int, int),
        Matrix_Free( SparseMatrix *),
        Generate_Dendrite(branch *, int *),
        Build_Test_Dendrite( branch **, branch *),
        LU_Factor(SparseMatrix *, int *),
        LU_Solve( SparseMatrix *, double *, double *),
        Input_Current( branch *),
        Assign_Current( branch *, double *, double );

/* Global definitions */
#define             CS    1.0
#define             GS    0.091
#define             GA    14.286
#define             CM    1.0
#define             GM    0.091
#define         OUTPUT    "Old100.dat"
#define           TEND    10.0
#define           NSIM    2000          /* Simulations to be done */
#define             NT    1000
#define             DT    1.0
#define          NODES    100
#define          NSEED    7             /* Seed for random number generator */
#define          FSEED    "Old100.ran"  /* History of random number generator */

/* Parameters for exact solution */
#define           NCON    75            /* Number of contacts */
#define            AMP    2.0e-5
#define            SIN    0.0e-3
#define              T    10.0
#define             RS    0.002

/* Global Variables */
SparseMatrix lhs, rhs;
double pi;
unsigned int ix, iy, iz;

int main( int argc, char **argv )
{
    extern unsigned int ix, iy, iz;
    extern double pi;
    int j, k, start, begin, nodes, n, FirstNode, repeat, maxstep,
        nstep, nb, ncon, nsim, connected, counter, NumberOfInput;
    double *v, *x, max, AreaOfSoma, frac, arg, sum,
           tmp, dt, h, len, CellLength, LocusContact;
    void srand( unsigned int);
    neuron *cell;
    contact *newcon, *oldcon, *con;
    extern SparseMatrix lhs, rhs;
    branch *bnow, *bold, *bnew, *FirstBranch, *CellFirstBranch;
    char word[100];
    FILE *fp;

/*  Initialise simulation counter */
    nsim = 1;
    start = 1;
    pi = 4.0*atan(1.0);
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
    maxstep = ((int) 1000.0*T);
    AreaOfSoma = 4.0*pi*RS*RS;
    if ( argc != 2 ) {
        printf("\n Invoke program with load <input>\n");
        return 1;
    } else {
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
                bnew->conlist = NULL;
                if ( bold != NULL) {
                    bold->child = bnew;
                } else {
                    CellFirstBranch = bnew;
                }
                bnew->parent = bold;
                fscanf(fp,"%lf %lf %lf", &(bnew->xl), &(bnew->yl), &(bnew->zl));
                fscanf(fp,"%lf %lf %lf", &(bnew->xr), &(bnew->yr), &(bnew->zr));
                fscanf(fp,"%lf %lf", &(bnew->plen), &(bnew->diam));
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
    while ( nsim <= NSIM ) {
        printf("\rSimulation %d", nsim);

/*  Step 1. - Generate a copy of the branch list */
        bnow = CellFirstBranch;
        bold = NULL;
        while ( bnow != NULL ) {
            bnew = (branch *) malloc( sizeof(branch) );
            bnew->id = bnow->id;
            bnew->xl = bnow->xl;
            bnew->yl = bnow->yl;
            bnew->zl = bnow->zl;
            bnew->xr = bnow->xr;
            bnew->yr = bnow->yr;
            bnew->zr = bnow->zr;
            bnew->plen = bnow->plen;
            bnew->diam = bnow->diam;
            bnew->peer = NULL;
            bnew->child = NULL;
            bnew->conlist = NULL;
            if ( bold != NULL) {
                bold->child = bnew;
            } else {
                FirstBranch = bnew;
            }
            bnew->parent = bold;
            bold = bnew;
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
            newcon->amp = AMP;
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
        if ( start ) printf("\n\nTree contains %d individual dendrite(s) ...\n", n);
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

/*  STEP 7A. - Count number on inputs on Cell */
        NumberOfInput = 0;
        for ( k=0 ; k<cell->ndend ; k++ ) {
            bnow = cell->dendlist[k].root;
            NumberOfInput += Count_Contacts( cell->dendlist[k].root, bnow);
        }

/*  STEP 8. - Count dendritic segments */
        for ( nb=k=0 ; k<cell->ndend ; k++ ) {
            bnow = cell->dendlist[k].root;
            nb += Count_Branches( bnow, bnow);
        }
        h = CellLength/((double) NODES-nb);
        for ( k=0 ; k<cell->ndend ; k++ ) Assign_Branch_Nodes( cell->dendlist[k].root, &h);

/*  STEP 9. - Enumerate Nodes */
        FirstNode = 0;
        for ( k=0 ; k<cell->ndend ; k++ ) Enumerate_Nodes( cell->dendlist[k].root, &FirstNode );
        for ( k=0 ; k<cell->ndend ; k++ ) cell->dendlist[k].root->junct = FirstNode;
        if ( start ) printf("Number of nodes is %d\n", FirstNode+1);

/*  STEP 10. - Construct Sparse Matrices */
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
        lhs.n = rhs.n = nodes;

/*  STEP 11. - Do somal node */
        lhs.a[3*nodes-3] = rhs.a[3*nodes-3] = 0.0;
        for ( k=0 ; k<cell->ndend ; k++ ) {
            bnow = cell->dendlist[k].root;
            lhs.a[counter] = 0.0;
            rhs.a[counter] = -0.25*pi*pow(bnow->diam,2)/(bnow->hseg);
            rhs.a[3*nodes-3] += 0.25*pi*pow(bnow->diam,2)/(bnow->hseg);
            lhs.a[3*nodes-3] += 0.5*pi*(bnow->diam)*(bnow->hseg);
            lhs.col[counter] = rhs.col[counter] = bnow->first;
            counter++;
        }
        lhs.col[3*nodes-3] = rhs.col[3*nodes-3] = nodes-1;
        lhs.StartRow[nodes] = rhs.StartRow[nodes] = 3*nodes-2;

/*  STEP 12. - Initialise input currents */
        for ( k=0 ; k<cell->ndend ; k++ ) Input_Current(cell->dendlist[k].root);
        dt = 1.0/((double) NT);
        for ( k=0 ; k<3*nodes-2 ; k++ ) {
            rhs.a[k] = 0.5*dt*(GA*rhs.a[k]+GM*lhs.a[k]);
            rhs.a[k] = CM*lhs.a[k]-rhs.a[k];
            lhs.a[k] = 2.0*CM*lhs.a[k]-rhs.a[k];
        }

/*  Add capacitive term of soma */
        lhs.a[3*nodes-3] += AreaOfSoma*(CS+0.5*GS*dt);
        rhs.a[3*nodes-3] += AreaOfSoma*(CS-0.5*GS*dt);
        if ( start ) {
            v = (double *) malloc( (nodes)*sizeof(double) );
            x = (double *) malloc( (nodes)*sizeof(double) );
            start = !start;
        }
        for ( k=0 ; k<nodes ; k++ ) v[k] = 0.0;
        begin = 1;
        LU_Factor(&lhs, &begin);
        if ( nsim == 1 ) {
            fp = fopen(OUTPUT, "w");
            fclose(fp);
        }

/*  Initialise temporal integration */
        nstep = 0;
        while ( nstep < maxstep ) {
            Matrix_Vector_Multiply( &rhs, v, x);
            x[nodes-1] -= dt*SIN;
            if ( nstep < maxstep ) {
                for ( k=0 ; k<cell->ndend ; k++ ) Assign_Current(cell->dendlist[k].root, x, dt);
            }
            LU_Solve( &lhs, v, x );
            nstep++;
            if ( nstep%1000 == 0 ) {
                fp = fopen(OUTPUT,"a");
                fprintf(fp,"%20.15lf",v[nodes-1]);
                fclose(fp);
            }
        }
        fp = fopen(OUTPUT, "a");
        fprintf(fp,"\n");
        fclose(fp);
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
        Function to constuct sparse matrices
  ******************************************************/
void Generate_Dendrite( branch *b, int *counter)
{
    int j, k, nc, CurrentNode;
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
            lhs.col[*counter] = rhs.col[*counter] = btmp->first;
            SumL += 0.5*pi*(btmp->diam)*(btmp->hseg);
            SumR += 0.25*pi*pow(btmp->diam,2)/(btmp->hseg);
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
    CurrentNode++;
    for ( j=(b->nc)-1 ; j>0 ; j-- ) {
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
        CurrentNode++;
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


 /******************************************************
             Function to destroy a NEURON
  ******************************************************/
void Destroy_Test_Neuron(neuron *cell)
{
    int k;

    for ( k=0 ; k<cell->ndend ; k++ ) Destroy_Test_Dendrite( cell->dendlist[k].root );
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
           Function to input current to dendrite
  **********************************************************/
void Input_Current( branch *b )
{
    int k;
    double ratio;
    contact *con;

    if ( b->child != NULL ) Input_Current( b->child );
    if ( b->peer != NULL ) Input_Current( b->peer );

    con = b->conlist;
    while ( con ) {
        ratio = (con->xp)/(b->hseg);
        if ( fmod(ratio,1.0) >= 0.5 ) {
            k = ((int) ceil(ratio));
        } else {
            k = ((int) floor(ratio));
        }
        if ( k == 0 ) {
            con->node = b->junct;
        } else {
            con->node = b->first-k+1;
        }
        con = con->next;
    }
    return;
}


 /*********************************************************
                 Function to assign current
  *********************************************************/
void Assign_Current(branch *b, double *x, double fac )
{
    contact *con;

    if ( b->child ) Assign_Current(b->child, x, fac );
    if ( b->peer ) Assign_Current(b->peer, x, fac );

    con = b->conlist;
    while ( con ) {
        x[con->node] -= fac*(con->amp);
        con = con->next;
    }
    return;
}


 /*********************************************************
            Function to assign branch nodes
  *********************************************************/
void Assign_Branch_Nodes( branch *b, double *h )
{
    b->nc = ((int) ceil((b->plen)/(*h)));
    b->hseg = (b->plen)/((double) b->nc);

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

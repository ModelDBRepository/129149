#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

 /***************************************************************
    Function to generate exact solutions and locations of stimuli
  ***************************************************************/

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
        int nc;                     /* Number of compartments specifying branch */
        int id;                     /* Identifies branch for NEURON calcs */
        double xl;                  /* X-coordinate of lefthand endpoint */
        double yl;                  /* Y-coordinate of lefthand endpoint */
        double zl;                  /* Z-coordinate of lefthand endpoint */
        double xr;                  /* X-coordinate of righthand endpoint */
        double yr;                  /* Y-coordinate of righthand endpoint */
        double zr;                  /* Z-coordinate of righthand endpoint */
        double diam;                /* Branch diameter (cm) */
        double plen;                /* Branch length (cm) */
        double elen;                /* Branch length (eu) */

        double *l;                  /* Length of dendritic segments (cm) */

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
        Output_Current( branch *, FILE *, FILE *),
        Remove_Branch( branch **, branch *),
        Destroy_Test_Neuron( neuron *),
        Destroy_Test_Dendrite( branch *),
        Find_Contacts( branch *, double *, double *, int *),
        Assign_Branch_Nodes( branch *, double *),
        Enumerate_Nodes( branch *, int *),
        Assign_Current( branch *, double *, double );

/* Global definitions */
#define             RS    0.002
#define             GA    14.286
#define             CM    1.0
#define             GM    0.091
#define        OUTPUT1    "InputCurrents.dat"
#define        OUTPUT2    "InputDendrite.dat"
#define      EXACTSOLN    "ExactSolution.dat"
#define           NSIM    2000      /* Simulations to be done */
#define             DT    1.0
#define          NODES    400
#define          NSEED    7         /* Seed for random number generator */
#define          FSEED    "GenRan75.ran"  /* History of random number generator */

/* Parameters for exact solution */
#define           NCON    75        /* Number of contacts */
#define            AMP    2.0e-5
#define            SIN    0.0e-3
#define              T    10.1
#define              M    1000

/* Global Variables */
unsigned int ix, iy, iz;

int main( int argc, char **argv )
{
    extern unsigned int ix, iy, iz;
    int k, j, id, start, n, nc, i, in, new,
        ncon, first, FirstNode, NumberOfInput;
    int counter, nb, nsim, connected;
    double max, *cval, *eta, *eval, AreaOfSoma, gama, *chi,
           xold, xnew, fac, arg, sum, tmp, vs, pi, tnow,
           len, h, CableDiameter, ElectrotonicLength, *amp, *loc,
           input, CableLength, dx, CellLength, LocusContact;
    void srand( unsigned int);
    neuron *cell;
    contact *newcon, *oldcon, *con;
    branch *bnow, *bold, *bnew, *FirstBranch, *CellFirstBranch;
    char word[20];
    FILE *fp, *fp1, *fp2;

/*  Initialise simulation counter */
    nsim = 1;
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
            printf("\n Test Neuron file not found");
            return 1;
        }
        pi = 4.0*atan(1.0);

/*  Get branch data */
        bold = NULL;
        while  ( fscanf(fp,"%s",word) != EOF ) {
            if ( strcmp(word,"Branch") == 0 || strcmp(word,"branch") == 0 ) {
                fscanf(fp,"%d", &k );
                bnew = (branch *) malloc( sizeof(branch) );
                bnew->id = k;
                bnew->peer = NULL;
                bnew->child = NULL;
                bnew->l = NULL;
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

/*  Step 0. - Start simulation procedure */
    start = 1;
    if ( nsim == 1 ) {
        new = first = 1;
    } else {
        new = first = 0;
    }
    while ( nsim <= NSIM ) {
        printf("\r Doing simulation %d", nsim);

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
            bnew->l = NULL;
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
            Find_Contacts(cell->dendlist[k].root, loc, amp, &ncon);
        }
        tmp = sqrt(CableDiameter)/CableLength;
        for ( k=0 ; k<ncon ; k++ ) loc[k] *= tmp;
        if ( ncon == 0 ) {
            printf("\n No contact found - Fatal error!");
            return 1;
        }

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
            chi[k] /= (fac*eta[k]*(1.0+gama*cval[k]*cval[k]));
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
        if ( start ) printf("\nNumber of nodes is %d\n", FirstNode+1);

/*  STEP 11. - Construct current input */
        if ( new ) {
            fp1 = fopen(OUTPUT1,"w");
            fp2 = fopen(OUTPUT2,"w");
            new = 0;
        } else {
            fp1 = fopen(OUTPUT1,"a");
            fp2 = fopen(OUTPUT2,"a");
        }
        for ( k=0 ; k<cell->ndend ; k++ ) {
            Output_Current( cell->dendlist[k].root, fp1, fp2 );
        }
        fprintf(fp1,"\n");
        fprintf(fp2,"\n");
        fclose(fp1);
        fclose(fp2);

/*  STEP 12. - Construct exact solution */
        for ( tnow=1.0 ; tnow<T ; tnow += 1.0 ) {
            for ( vs=0.0,k=M ; k>=0 ; k-- ) {
                arg = tnow*eta[k];
                if ( arg > 20.0 ) {
                    tmp = 1.0;
                } else {
                    tmp = (1.0-exp(-arg));
                }
                vs -= tmp*chi[k];
            }
            if ( first ) {
                fp = fopen(EXACTSOLN, "w");
                first = 0;
            } else {
                fp = fopen(EXACTSOLN, "a");
            }
            fprintf(fp,"%12.6lf",vs);
            fclose(fp);
        }
        fp = fopen(EXACTSOLN, "a");
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


 /***************************************************
    Function to output location of current input
  ***************************************************/
void Output_Current( branch *b, FILE *fp1, FILE *fp2 )
{
    int i;
    double frac;
    contact *con;

    if ( b->child ) Output_Current(b->child, fp1, fp2);
    if ( b->peer ) Output_Current(b->peer, fp1, fp2);

    con = b->conlist;
    while ( con ) {
        frac = (con->xp)/(b->plen);
        fprintf(fp1,"%10.7lf", frac);
        fprintf(fp2,"%4d", b->id);
        con = con->next;
    }
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
    free( b->l );
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
int Count_Contacts( branch *head, branch *bnow)
{
    static int n;
    contact *con;

    if ( head == bnow ) n = 0;
    if ( bnow != NULL ) {
        if ( bnow->child ) Count_Contacts( head, bnow->child);
        if ( bnow->peer ) Count_Contacts( head, bnow->peer);
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



 /**************************************************
            Function to assign branch nodes
  **************************************************/
void Assign_Branch_Nodes( branch *b, double *h )
{
    int nc, k;
    double hseg;

    nc = ((int) ceil((b->plen)/(*h)));
    hseg = (b->plen)/((double) nc);
    b->l = (double *) malloc( nc*sizeof(double));
    for ( k=0 ; k<nc ; k++ ) b->l[k] = hseg;
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

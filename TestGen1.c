#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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
        double pl;                  /* Measurement of physical length (micron) */
        double amp;                 /* Strength of contact */

        int xl;                     /* Left hand node */
        int xr;                     /* Right hand node */
        double frac;                /* Fraction of input to left hand node */

        struct contact_t *next;     /* Address of next contact */
} contact;


typedef struct soma_t
{
/*  Static biophysical properties of soma */
        double cs;                  /* Somal membrane capacitance (mu F/cm^2) */
        double ga;                  /* Intracellular conductance of soma (mS/cm) */
        double gs;                  /* Membrane conductance of soma (mS/cm^2) */

/*  Contact information */
        contact *conlist;           /* List of contacts */
} soma;



typedef struct branch_t
{
/*  Connectivity of branch */
        struct branch_t *parent;  /* Pointer to parent branch */
        struct branch_t *child;   /* Pointer to child branch */
        struct branch_t *peer;    /* Pointer to a peer branch */

/*  Physical properties of branch */
        int    nd;                /* Number of nodes */
        double xs;                /* X-coord of somal endpoint */
        double ys;                /* Y-coord of somal endpoint */
        double zs;                /* Z-coord of somal endpoint */
        double ds;                /* Diameter of somal endpoint */
        double xd;                /* X-coord of distal endpoint */
        double yd;                /* Y-coord of distal endpoint */
        double zd;                /* Z-coord of distal endpoint */
        double dd;                /* Diameter of distal endpoint */
        double *d;                /* Diameter information */
        double *x;                /* Location information */

/*  Biophysical properties of branch */
        double cm;                /* Dendritic membrane capacitance (mu F/cm^2) */
        double ga;                /* Intracellular conductance of dendrite (mS/cm) */
        double gm;                /* Membrane conductance of dendrite (mS/cm^2) */

/*  Node information for spatial representation */
        int nodes;                /* Total number nodes spanning branch */
        int junct;                /* Junction node of the branch */
        int first;                /* Internal node connected to junction */

/*  Contact information */
        contact *conlist;         /* List of contacts */
} branch;



typedef struct dendrite_t
{
        branch *root;             /* Pointer to root branch of dendrite */
        double plen;              /* Total length of dendrite */
} dendrite;


typedef struct neuron_t
{
        int ndend;                  /* Number of dendrites */
        dendrite *dendlist;         /* Pointer to an array of dendrites */
        soma *s;                    /* Soma structure */
} neuron;


/* Function type declarations */
neuron *LoadTestNeuron(char *);

void BuildTestDendrite( branch **, branch *),
     RemoveBranch( branch **, branch *),
     DestroyTestNeuron( neuron * ),
     DestroyTestDendrite( branch *);

void BuildContactInfo(contact *, branch *, branch **);

int CountBranches( branch *, branch *),
    CountContacts( branch *, branch *);

double branch_length( branch *, branch *);

int count_terminal_branches( branch *, branch *);
void output_properties( branch * );
void enumerate_nodes( branch *, int *);

void FreeSparseMatrix( SparseMatrix *),
     MatrixVectorMultiplication( SparseMatrix *, double *, double *),
     SparseMatrixMalloc( SparseMatrix *, int, int),
     LU_Factor( SparseMatrix *, int *),
     LU_Solve( SparseMatrix *, double *, double *);

void Generate_Dendrite(branch *, int *);

void Input_Current( branch *);
void Assign_Current( branch *, double *, double );

void count_dendritic_length ( branch *, double *);
void phys_lengths( branch *, double * );

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
#define           TEND    2.0
#define             NT    100
#define             DT    0.1
#define            SIN    0.0e-3
#define              T    100.0
#define             RS    0.002
#define         NNODES    100
#define         SENTRY    50.0


/* Global Variables */
SparseMatrix lhs, rhs, lhs1;

int main( int argc, char **argv )
{
    int i, k, id, start, nodes, nc, in, FirstNode, counter,
        nb, nspk, first;
    double jo, je, jn, ho, he, hn, mo, me, mn, no, ne, nn,
           vo, ve, vn, xo, xn;
    double pi, as;
    double *v, *x, max, gama, vs, frac, arg,
           sum, tmp, tmp1, tmp2, dt, tnow, tout, left,
           rite, gval, cval, len, h, vmid;
    double g_na = 120.0, g_k = 36.0, g_l = 0.3, v_na = 55.0,
           v_k = -72.0, v_l = -49.416, fact=25.0;
    neuron *n;
    extern SparseMatrix lhs, rhs, lhs1, ptr;
    branch *bnow;
    FILE *fp;

/*  Compute ancillary variables */
    pi = 4.0*atan(1.0);
    as = 4.0*pi*RS*RS;

/*  Compute Equilibrium Potential and equilibrium states */
    vo = -62.0;
    mo = alfa_m(vo)/(alfa_m(vo)+beta_m(vo));
    ho = alfa_h(vo)/(alfa_h(vo)+beta_h(vo));
    no = alfa_n(vo)/(alfa_n(vo)+beta_n(vo));
    jo = g_na*pow(mo,3)*ho*(vo-v_na)+g_k*pow(no,4)*(vo-v_k)+g_l*(vo-v_l);

    vn = -58.0;
    mn = alfa_m(vn)/(alfa_m(vn)+beta_m(vn));
    hn = alfa_h(vn)/(alfa_h(vn)+beta_h(vn));
    nn = alfa_n(vn)/(alfa_n(vn)+beta_n(vn));
    jn = g_na*pow(mn,3)*hn*(vn-v_na)+g_k*pow(nn,4)*(vn-v_k)+g_l*(vn-v_l);

    if ( jo*jn > 0.0 ) {
        printf(" No zero found \n");
        return(0);
    } else {
        while ( fabs(vo-vn) > 5.e-7 ) {
            ve = 0.5*(vo+vn);
            me = alfa_m(ve)/(alfa_m(ve)+beta_m(ve));
            he = alfa_h(ve)/(alfa_h(ve)+beta_h(ve));
            ne = alfa_n(ve)/(alfa_n(ve)+beta_n(ve));
            je = g_na*pow(me,3)*heq*(ve-v_na)+g_k*pow(ne,4)*(ve-v_k)+g_l*(veq-v_l);
            if ( je*jo > 0.0 ) {
                vo = ve;
            } else {
                vn = ve;
            }
        }
    }
    vn = ve;

    mn = alfa_m(vn)/(alfa_m(vn)+beta_m(vn));
    hn = alfa_h(vn)/(alfa_h(vn)+beta_h(vn));
    nn = alfa_n(vn)/(alfa_n(vn)+beta_n(vn));

    v_na -= ve;
    v_k -= ve;
    v_l -= ve;

/*  Load sampled neuron */
    if ( argc != 2 ) {
        printf("\n Invoke program with load <input>\n");
        return(1);
    } else {
        n = LoadTestNeuron( argv[1] );
        if ( !n ) {
            printf("\n Failed to find test neuron\n");
            return(1);
        }
        len = 0.0;
        nb = 0;
        for ( i=0 ; i<n->ndend ; i++) count_branch( n->dendlist[i].root, &nb );
        for ( i=0 ; i<n->ndend ; i++) count_dendritic_length( n->dendlist[i].root, &length );
        h = len/((double) NNODES - nb);
        for( i=0 ; i<n->ndend ; i++) phys_lengths( n->dendlist[i].root, &h );
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

    mat_malloc( &lhs1, nodes, 3*nodes-2);

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

    dt = 1.0/((double) NT);

    for ( i=0 ; i<3*nodes-2 ; i++ ) {
        rhs.a[i] *= GA;
        rhs.a[i] += GM*lhs.a[i];
        lhs.a[i] *= CM;
        rhs.a[i] *= 0.5*dt;
        lhs.a[i] += rhs.a[i];
        rhs.a[i] = lhs.a[i] - 2.0*rhs.a[i];
    }
    left = lhs.a[3*nodes-3]+(4.0*as/pi)*CS;
    rite = rhs.a[3*nodes-3]+(4.0*as/pi)*CS;

    v = (double *) malloc( (nodes)*sizeof(double) );
    x = (double *) malloc( (nodes)*sizeof(double) );
    for ( i=0 ; i<nodes ; i++ ) v[i] = 0.0;

/*  Initialise temporal integration */

    tnow = 0.0;
    tout = DT;
    vold = vnew = vmid = 0.0;
    nspk = 0;
    first = 1;

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

        lhs.a[3*nodes-3] = left + gval;
        rhs.a[3*nodes-3] = rite - gval;

        for ( k=0 ; k<3*nodes-2 ; k++ ) {
            lhs1.a[k] = lhs.a[k];
            lhs1.col[k] = lhs.col[k];
        }
        for ( k=0 ; k<=nodes ; k++ ) lhs1.start_row[k] = lhs.start_row[k];

        OldLU_Factor(&lhs);
        fp = fopen("out-l2.dat","w");
        for ( i=0 ; i<nodes ; i++ ) {
            for ( k=lhs.l->start_row[i] ; k<lhs.l->start_row[i+1] ; k++ ) {
                fprintf(fp,"%d \t (%d,%d) \t %lf\n", k, i, lhs.l->col[k], lhs.l->a[k]);
            }
        }
        fclose(fp);
        fp = fopen("out-u2.dat","w");
        for ( i=0 ; i<nodes ; i++ ) {
            for ( k=lhs.u->start_row[i] ; k<lhs.u->start_row[i+1] ; k++ ) {
                fprintf(fp,"%d \t (%d,%d) \t %lf\n", k, i, lhs.u->col[k], lhs.u->a[k]);
            }
        }
        fclose(fp);
        //return 0;

        LU_Factor(&lhs1, &first);
        fp = fopen("out-l1.dat","w");
        for ( i=0 ; i<nodes ; i++ ) {
            for ( k=lhs1.l->start_row[i] ; k<lhs1.l->start_row[i+1] ; k++ ) {
                fprintf(fp,"%d \t (%d,%d) \t %lf\n", k, i, lhs1.l->col[k], lhs1.l->a[k]);
            }
        }
        fclose(fp);
        fp = fopen("out-u1.dat","w");
        for ( i=0 ; i<nodes ; i++ ) {
            for ( k=lhs1.u->start_row[i] ; k<lhs1.u->start_row[i+1] ; k++ ) {
                fprintf(fp,"%d \t (%d,%d) \t %lf\n", k, i, lhs1.u->col[k], lhs1.u->a[k]);
            }
        }
        fclose(fp);

        tmp = 0.0;
        for ( i=0 ; i<nodes ; i++ ) {
            for ( k=lhs.l->start_row[i] ; k<lhs.l->start_row[i+1] ; k++ ) {
                if ( fabs(lhs1.l->a[k]-lhs.l->a[k]) > tmp ) {
                    tmp = fabs(lhs1.l->a[k]-lhs.l->a[k]);
                }
            }
        }
        if ( tmp > 5.e-12 ) {
            printf("\n trouble %lf %lf",tnow, tmp);
            getchar( );
        }

       // return 0;
        mat_vec_mult(&rhs,v,x);
        x[nodes-1] -= 4.0*dt*SIN/pi;
        x[nodes-1] += cval;
        tmp = 4.0*dt/pi;
        if (tnow<T) for (i=0;i<n->ndend;i++) Assign_Current(n->dendlist[i].root,x,tmp);

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
    }

/*  Count contacts
    for ( n=k=0 ; k<n->ndend ; k++ ) {
        nc += count_contacts( n->dendlist[k].root, n->dendlist[k].root);
    }
    printf("\n Located %d contacts on dendrites", nc);
    printf("\n      Located %d contacts on soma", n->s->ncon);
    printf("\n");  */

    DestroyTestNeuron( n );
    return(0);
}


 /*************************************************
          Function To Load A Test Neuron
  *************************************************/
neuron *LoadTestNeuron(char *filename)
{
    int j, k, ncon, n, id, connected, ignored;
    double tmp, piby2, xl, xr, dl, dr, px, py, pz, min, radius, dx;
    neuron *cell;
    contact *newcon;
    branch *bold, *bnew, *FirstBranch;
    char temp[100];
    FILE *fp;

/*  STEP 1. - Open neuron data file */
    printf("\nOpening file %s\n",filename);
    if ( (fp=fopen(filename,"r"))==NULL ) return NULL;

/*  STEP 2. - Get memory for neuron structure */
    cell = (neuron *) malloc( sizeof(neuron) );

/*  STEP 3. - Get branch and contact data */
    bo = NULL;
    while  ( fscanf(input,"%s", temp)!=EOF ) {
        if ( strcmp(temp,"Branch") == 0 || strcmp(temp,"branch") == 0 ) {
            fscanf(fp, "%d", &n);
            printf("Found a branch\n");
            bnew = (branch *) malloc( sizeof(branch) );
            bnew->nd = n;
            if ( bold ) {
                bold->child = bnew;
            } else {
                FirstBranch = bnew;
            }
            bnew->parent = bold;
            bnew->peer = NULL;
            bnew->child = NULL;

/*  STEP 3b. - Initialise branch */
            bnew->d = (double *) malloc( n*sizeof(double) );
            bnew->x = (double *) malloc( n*sizeof(double) );
            bnew->cm = CM;
            bnew->gm = GM;
            bnew->ga = GA;
            bnew->conlist = NULL;

/*  STEP 3c. - Read branch morphology */
            fscanf(fp,"%lf %lf %lf %lf", &(bnew->xs), &(bnew->ys), &(bnew->zs), &(bnew->ds) );
            fscanf(fp,"%lf %lf %lf %lf", &(bnew->xd), &(bnew->yd), &(bnew->zd), &(bnew->dd) );
            fscanf(fp,"%lf", &len );
            dx = len/((double) n-1 );
            dd = (bnew->dd-bnew->ds)/((double) n-1);
            for ( j=0 ; j<n ; j++ ) {
                bnew->x[j] = dx*((double) j);
                bnew->d[j] = ds+dd*((double) j);
            }
            bold = bnew;
        } else if ( strcmp(temp, "Marker") == 0 || strcmp(temp, "marker") == 0 ) {

/*  STEP 3d. - Initialise marker */
            printf("Found and initialised a branch contact\n");
            newcon = (contact *) malloc( sizeof(contact) );
            newcon->next = NULL;
            fscanf(fp,"%lf %lf", &newcon->xp, &newcon->amp );
            if ( bnew->conlist == NULL ) {
                bnew->conlist = newcon;
            } else {
                oldcon = bnew->conlist;
                while ( oldcon->next ) oldcon = oldcon->next;
                oldcon->next = newcon;
            }
        } else {
            printf("Unrecognised dendritic feature\n");
        }
    }
    fclose(fp);

/*  STEP 4. - Count dendritic branches at soma */
    bold = FirstBranch;
    n = 0;
    while ( bold ) {
        bnew = FirstBranch;
        do {
            tmp = pow(bold->xs-bnew->xd,2)+pow(bold->ys-bnew->yd,2)
                  +pow(bold->zs-bnew->zd,2);
            connected = ( tmp < 0.01 );
            bnew = bnew->child;
        } while ( bnew && !connected );
        if ( !connected ) n++;
        bold = bold->child;
    }
    cell->ndend = n;
    printf("\n\nTree contains %d individual dendrite(s) ...\n", n);

/*  STEP 5. - Identify somal dendrites but extract nothing */
    cell->dendlist = (dendrite *) malloc( (cell->ndend)*sizeof(dendrite) );
    bold = FirstBranch;
    n = 0;
    while ( n < cell->ndend ) {
        bnew = FirstBranch;
        do {
            tmp = pow(bold->xs-bnew->xd,2)+pow(bold->ys-bnew->yd,2)
                  +pow(bold->zs-bnew->zd,2);
            connected = ( tmp < 0.01 );
            bnew = bnew->child;
        } while ( bnew && !connected );
        if ( !connected ) cell->dendlist[n++].root = bold;
        bold = bold->child;
    }

/*  STEP 6. - Extract root of each dendrite from dendrite list */
    for ( k=0 ; k<cell->ndend ; k++ ) {
        bold = cell->dendlist[k].root;
        RemoveBranch( &FirstBranch, bold);
    }

/*  STEP 7. - Build each test dendrite from its root branch */
    for ( k=0 ; k<cell->ndend ; k++ ) {
        BuildTestDendrite( &FirstBranch, cell->dendlist[k].root);
    }
    if ( FirstBranch != NULL ) printf("\nWarning: Unconnected branch segments still exist\n");
    return cell;
}


 /**************************************************
    Function to remove a branch from a branch list
  **************************************************/
void RemoveBranch( branch **head, branch *b)
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


 /********************************************************
      Function to build a test dendrite from its root
  ********************************************************/
void BuildTestDendrite( branch **head, branch *root)
{
    double tmp;
    branch *bnow, *bnext, *btmp;

    bnow = *head;
    while ( bnow != NULL ) {

/*  Store bnow's child in case it's corrupted */
        bnext = bnow->child;

/*  Search if proximal end of bnow is connected to distal end of root */
        tmp = pow(bnow->xs-root->xd,2)+pow(bnow->ys-root->yd,2)+
              pow(bnow->zs-root->zd,2);
        if ( tmp <= 0.01 ) {

/*  Take bnow out of the branch list */
            remove_branch( head, bnow);

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

/*  Iterate through remaining tree */
    if ( root->child ) BuildTestDendrite( head, root->child);
    if ( root->peer ) BuildTestDendrite( head, root->peer);
    return;
}


 /*****************************************************
             Function to destroy a NEURON
  *****************************************************/
void DestroyTestNeuron( neuron *cell)
{
    int k;
    contact *prevcon, *nextcon;

/*  Free Soma */
    if ( cell->s != NULL ) {
        prevcon = cell->s->conlist;
        while ( prevcon ) {
            nextcon = prevcon->next;
            free ( prevcon );
            prevcon = nextcon;
        }
        free ( cell->s );
    }
    for ( k=0 ; k<cell->ndend ; k++ ) DestroyTestDendrite( cell->dendlist[i].root );
    free(cell);
    return;
}


 /************************************************
          Function to destroy Test DENDRITE
  ************************************************/
void DestroyTestDendrite( branch *b )
{
    contact *prevcon, *nextcon;

    if ( b->child ) DestroyTestDendrite(b->child);
    if ( b->peer ) DestroyTestDendrite(b->peer);
    free( b->x );
    free( b->d );
    prevcon = b->conlist;
    while ( prevcon ) {
        nextcon = prevcon->next;
        free ( prevcon );
        prevcon = nextcon;
    }
    free( b );
    return;
}


 /****************************************************
          Function to count number of branches
          from current branch to dendritic tip
  ****************************************************/
int CountBranches( branch *bstart, branch *bnow)
{
    static int n;

    if ( bstart == bnow ) n = 0;
    if ( bnow ) {
        if ( bnow->child ) CountBranches( bstart, bnow->child);
        if ( bnow->peer ) CountBranches( bstart, bnow->peer);
        n++;
    }
    return n;
}


 /*****************************************************
      Function to count number of contacts
      from current branch to the dendritic tip.
  *****************************************************/
int CountContacts( branch *bstart, branch *bnow)
{
    static int n;
    contact *con;

    if ( bstart == bnow ) n = 0;
    if ( bnow ) {
        if ( bnow->child ) CountContacts(bstart, bnow->child);
        if ( bnow->peer ) CountContacts(bstart, bnow->peer);
        con = bnow->conlist;
        while ( con ) {
            n++;
            con = con->next;
        }
    }
    return n;
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



/****************************************************************

               Function to output branch diameters

****************************************************************/

void output_properties( branch *b )
{
    int i, k;
    static int start=1;
    double dold, dnew, len, xold, yold, zold, xnew, ynew, znew, dx, dy, dz, size;
    branch *bran;
    FILE *fp;

    if ( b->child ) output_properties(b->child);
    if ( b->peer ) output_properties(b->peer);
    if ( start ) {
        fp = fopen("output","w");
        start = 0;
    } else {
        fp = fopen("output","a");
        fprintf(fp,"\n");
    }

/*  Outputs branch lengths, diameters, surface areas etc.
    for ( k=0 ; k<b->nobs ; k++ ) {
        fprintf(fp,"%6.2lf \t %6.2lf \t %6.2lf \t %6.2lf \n", b->pl[k], b->d[k], b->sa[k], b->el[k]);
    }*/

/*  Decomposes branches into lengths of uniform diameter
    len = xold = b->pl[1];
    dold = b->d[1];
    for ( k=2 ; k<b->nobs ; k++ ) {
        xnew = b->pl[k];
        dnew = b->d[k];
        if ( dnew != dold ) {
            len += 0.5*(xnew-xold);
            fprintf(fp,"%6.2lf \t %6.2lf \n", len, dold);
            len = 0.5*(xnew-xold);
        } else {
            len += xnew-xold;
        }
        xold = xnew;
        dold = dnew;
    }
    fprintf(fp,"%6.2lf \t %6.2lf \n", len, dold); */

/*  Constructs diameters of a branch and its children/peers */
    if ( b->child ) {
        fprintf(fp,"%6.2lf \t %6.2lf \t", b->d[(b->nobs)-1], b->child->d[1]);
        bran = b->child;
        while ( bran->peer ) {
            bran = bran->peer;
            fprintf(fp,"%6.2lf \t", bran->d[1]);
        }
    }

/*  Prints out branch lengths
    printf("\nBranch length %6.2lf, %6.2lf, %6.2lf", b->p_len, b->d[0], b->d[b->nobs-1] );
    getchar( ); */
    fclose(fp);
    return;
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


 /***************************************************
        Allocate memory to a sparse matrix
  ***************************************************/
void SparseMatrixMalloc( SparseMatrix *a, int n, int w)
{
    a->a = (double *) malloc( w*sizeof(double) );
    a->col = (int *) malloc( w*sizeof(int) );
    a->StartRow = (int *) malloc( (n+1)*sizeof(int) );
    a->n = n;
    a->l = malloc( sizeof(SparseMatrix) );
    a->u = malloc( sizeof(SparseMatrix) );
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


 /********************************************************
     Multiplies sparse matrix a[ ][ ] with vector v[ ]
  ********************************************************/
void MatrixVectorMultiplication( SparseMatrix *a, double *v , double *b)
{
    int i, j, k, n;

    n = a->n;
    for ( j=0 ; j<n ; j++ ) {
        k = a->StartRow[j+1];
        for ( b[j]=0.0,i=(a->StartRow[j]) ; i<k ; i++ ) {
            b[j] += (a->a[i])*v[a->col[i]];
        }
    }
    return;
}


 /**********************************************
      De-allocates memory of a sparse matrix
  **********************************************/
void FreeSparseMatrix( SparseMatrix *a)
{
    free(a->a);
    free(a->col);
    free(a->StartRow);
    free(a);
}


 /***************************************************************
            Function To Factorise A Sparse Matrix
 ***************************************************************/
void LU_Factor( SparseMatrix *m, int *start)
{
    double tmp, sum;
    int i, j, k, r, n, cl, cu, nrow;

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

/*  Step 3. - Fill first row of L and U */
    m->l->a[0] = 1.0;
    for ( k=0 ; k < m->u->StartRow[1] ; k++ ) m->u->a[k] = m->a[k];

/*  Step 4. - Fill remaining entries row by row */
    cl = 1;
    k = cu = m->u->StartRow[1];
    for ( i=1 ; i<n ; i++ ) {
        while ( m->col[k] < i ) { // Fill lower matrix
            sum = m->a[k];
            for ( j=m->l->StartRow[i] ; j<cl ; j++ ) {
                nrow = m->u->StartRow[m->l->col[j]];
                while ( m->u->col[nrow] < m->col[k] ) nrow++;
                sum -= (m->l->a[j])*(m->u->a[nrow]);
            }
            nrow = m->u->StartRow[m->l->col[cl]];
            while ( m->u->col[nrow] < m->col[k] ) nrow++;
            m->l->a[cl++] = sum/(m->u->a[nrow]);
            k++;
        }

        m->l->a[cl++] = 1.0;       // Diagonal entry of lower
        while ( m->col[k] >= i && k < m->start_row[i+1] ) { // Fill upper matrix
            sum = m->a[k];
            for ( j=m->l->StartRow[i] ; j<m->l->StartRow[i+1]-1; j++ ) {
                nrow = m->u->StartRow[m->l->col[j]];
                while ( m->u->col[nrow] < m->col[k] && nrow < m->u->start_row[m->l->col[j]+1] ) nrow++;
                sum -= (m->l->a[j])*(m->u->a[nrow]);
            }
            k++;
            m->u->a[cu++] = sum;
        }
    }
    return;
}


 /****************************************************
         Function to Solve the matrix problem
  ****************************************************/
void LU_Solve( SparseMatrix *m, double *x, double *b )
{
    int i, j;
    double *z;

    z = (double *) malloc( (m->n)*sizeof(double) );

    for ( i=0 ; i<m->n ; i++ ) {
        z[i] = b[i];
        for ( j=m->l->StartRow[i] ; j<m->l->StartRow[i+1]-1 ; j++ ) {
            z[i] -= (m->l->a[j])*(z[(m->l->col[j])]);
        }
        z[i] /= m->l->a[m->l->StartRow[i+1]-1];
    }

    for ( i=(m->n)-1 ; i>=0 ; i-- ) {
        x[i] = z[i];
        for ( j=m->u->StartRow[i]+1 ; j<m->u->StartRow[i+1] ; j++ ) {
            x[i] -= (m->u->a[j])*(x[m->u->col[j]]);
        }
        x[i] /= m->u->a[m->u->StartRow[i]];
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



 /**********************************************************
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


 /**********************************************************
               ALPHA for ACTIVATION OF SODIUM
  **********************************************************/
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


 /**********************************************************
                    BETA for ACTIVATION OF SODIUM
  **********************************************************/
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


 /***********************************************************
              ALPHA for INACTIVATION OF SODIUM
  ***********************************************************/
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


 /**********************************************************
               ALPHA for ACTIVATION OF POTASSIUM
  **********************************************************/
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


 /*********************************************************
              BETA for ACTIVATION OF POTASSIUM
  *********************************************************/
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

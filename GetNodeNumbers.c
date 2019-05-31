#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

 /*****************************************************
    Function to obtain Node Numbers for the numerical
    analysis of the New Compartmental Model.
  *****************************************************/

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

        double *l;                  /* Length of dendritic segments (cm) */

/*  Node information for spatial representation */
        int nodes;                  /* Total number nodes spanning branch */
        int junct;                  /* Junction node of the branch */
        int first;                  /* Internal node connected to junction */

/*  Contact information */
        contact *conlist;           /* Branch contact */
} branch;


/* Global definitions */
#define           INFO    "Info500.dat"
#define          NODES    500

int main( int argc, char **argv )
{
    int k, nb;
    double h, CellLength;
    contact *newcon, *con;
    branch *bold, *bnew, *FirstBranch;
    char word[20];
    FILE *fp;

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

/*  Get branch data */
        bold = NULL;
        while  ( fscanf(fp,"%s",word) != EOF ) {
            if ( strcmp(word,"Branch") == 0 || strcmp(word,"branch") == 0 ) {
                bnew = (branch *) malloc( sizeof(branch) );
                fscanf(fp,"%d", &(bnew->id) );
                bnew->peer = NULL;
                bnew->child = NULL;
                bnew->l = NULL;
                bnew->conlist = NULL;
                if ( bold != NULL) {
                    bold->child = bnew;
                } else {
                    FirstBranch = bnew;
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
    bold = FirstBranch;
    while ( bold ) {
       CellLength += bold->plen;
       bold = bold->child;
    }

/*  STEP 2. - Count branches */
    bold = FirstBranch;
    nb = 0;
    while ( bold ) {
        nb++;
        bold = bold->child;
    }

/*  STEP 3. - Compute number of compartments */
    h = CellLength/((double) NODES-nb);
    bold = FirstBranch;
    while ( bold ) {
        bold->nc = ((int) ceil((bold->plen)/h));
        bold = bold->child;
    }

/*  STEP 4. - Output branch properties */
    fp = fopen(INFO,"w");
    for ( k=0 ; k<nb ; k++ ) {
        bold = FirstBranch;
        while ( bold->id != k ) bold = bold->child;
        fprintf(fp,"%3d \t %3d\n", bold->id, bold->nc);
    }
    fclose(fp);
    return 0;
}

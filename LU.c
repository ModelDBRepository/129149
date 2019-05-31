#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define N 16
#define ALPHA 50

typedef struct sparse_mat_t
{
        double *a;
        int *col;
        int *start_row;
        int n;
        struct sparse_mat_t *l;
        struct sparse_mat_t *u;

} sparse_mat;

void oldLU_Factor(sparse_mat *);

void LU_Factor(sparse_mat *, int *);
double me( sparse_mat *, int , int );
void mat_malloc( sparse_mat *, int , int );
void LU_Solve(sparse_mat *, double *, double *);

void main (void)
{
    int i, k, start=1;
    double b[N-1], x[N-1], h, v[N-1], norm;
    FILE *fp;

    sparse_mat matrix;

    mat_malloc(&matrix, N-1 , 3*N-5 );
    h = 1.0/((double) N);

    for ( k=0 ; k<N-1 ; k++ ) {
        x[k] = ((double) k+1)/((double) N);
        b[k] = h*h*sin(ALPHA*x[k]);
        v[k] = (x[k]*sin(ALPHA) - sin(ALPHA*x[k]))/pow(ALPHA,2);
    }

    matrix.a[0] = -2.0;
    matrix.a[1] = 1.0;
    matrix.col[0] = 0;
    matrix.col[1] = 1;
    matrix.start_row[0] = 0;

    for (k=1; k<N-2; k++) {
        matrix.start_row[k] = 3*k-1;

        matrix.a[3*k-1] = 1.0;
        matrix.col[3*k-1] = k-1;
        matrix.a[3*k] = -2.0;
        matrix.col[3*k] = k;
        matrix.a[3*k+1] = 1.0;
        matrix.col[3*k+1] = k+1;
    }

    matrix.start_row[N-2] = 3*N-7;

    matrix.a[3*N-7] = 1.0;
    matrix.col[3*N-7] = N-3;
    matrix.a[3*N-6] = -2.0;
    matrix.col[3*N-6] = N-2;

    matrix.start_row[N-1] = 3*N-5;

// oldLU_Factor(&matrix);
    LU_Factor(&matrix, &start);
    LU_Solve(&matrix,b,b);

    fp = fopen("out.res", "w");
    norm = 0.0;
    for ( i=0 ; i<N-2 ; i++ ) {
        norm += pow(b[i] - v[i],2);
        fprintf(fp,"%12.6lf \t %d \n",matrix.u->a[i], matrix.u->col[i]);
    }
    norm /= ((double) N-1);
    norm = sqrt(norm);
    printf("%10.7lf",norm);
    return;
}


void LU_Factor(sparse_mat *m, int *start) {

    double tmp, sum;
    int i, cl, cu, k, kk, j, n, nrow, repeat;

/*  Step 1. - Identify matrix dimension */
    n = m->n;

/* Step 2. - Fill column and row vectors for triangular matrices */
    if ( *start ) {
        cl = cu = 0;
        for ( i=k=0 ; i<n ; i++ ) {
            m->l->start_row[i] = cl;
            m->u->start_row[i] = cu;
            while ( m->col[k] < i ) m->l->col[cl++] = m->col[k++];
            m->l->col[cl++] = m->u->col[cu++] = m->col[k++];
            while ( m->col[k] > i && k < 3*n-2  ) m->u->col[cu++] = m->col[k++];
        }
        m->l->start_row[n] = cl;
        m->u->start_row[n] = cu;
        *start = 0;

/* Step 2. - Fill main diagonal of lower triangular matrix */
        for ( i=0 ; i<n ; i++ ) m->l->a[m->l->start_row[i+1]-1] = 1.0;
    }

/* Step 3. - Fill in matrix entries for first row of upper triangualar
             matrix and first column of lower triangular matrix */
    for ( i=0 ; i < m->u->start_row[1] ; i++ ) m->u->a[i] = m->a[i];
    for ( i=1 ; i < n ; i++ ) {
        k =  m->l->start_row[i];
        if ( m->l->col[k] == 0 ) m->l->a[k] = m->a[m->start_row[i]]/m->a[0];
    }

/* Step 4. - Fill in remaining upper and lower triangular matrix entries */
    for ( i=1 ; i<n ; i++ ) {
        for ( j=m->u->start_row[i] ; j < m->u->start_row[i+1] ; j++ ) {
            sum = 0.0;
            for ( k=m->l->start_row[i] ; k < m->l->start_row[i+1]-1 ; k++ ) {
                nrow = m->u->start_row[m->l->col[k]];
                while ( m->u->col[nrow] < m->u->col[j] ) nrow++;
                sum += (m->l->a[k])*(m->u->a[nrow]);
            }
            nrow = j+m->l->start_row[i+1]-i-1;
            m->u->a[j] = m->a[nrow]-sum;
        }
        for ( j=i+1 ; j<n ; j++ ) {
            k = m->l->start_row[j];
            repeat = 1;
            do {
                if ( m->l->col[k] > i ) {
                    repeat = 0;
                } else if ( m->l->col[k] == i ) {
                    sum = 0.0;
                    kk = m->l->start_row[j];
                    while ( m->l->col[kk] < i ) {
                        nrow = m->u->start_row[m->l->col[kk]];
                        while ( m->u->col[nrow] < i ) nrow++;
                        sum += (m->l->a[kk])*(m->u->a[nrow]);
                    }
                    nrow = kk+m->u->start_row[j]-j-1;
                    m->l->a[k] = (m->a[nrow]-sum)/(m->u->a[m->u->start_row[i]]);
                    repeat = 0;
                } else {
                    k++;
                }
            } while ( repeat && k < m->l->start_row[j+1]-1 );
        }
    }
    return;
}


double me( sparse_mat *m, int row, int col ) {

    int i;

    for (i = m->start_row[row]; i < m->start_row[row + 1]; i++) if(m->col[i] == col) return m->a[i];
    return 0.0;
}


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



void oldLU_Factor(sparse_mat *m) {

    double tmp;
    int i,k,r,n;

    m->u->start_row[0] = m->l->start_row[0] = 0;
    m->u->col[0] = m->l->col[0] = 0;

    /* Step.1 - Fill in column and row vectors for lower triangualr matrix */


    for ( i=r=k = 0 ; i<m->n ; i++ ) {
        k = m->start_row[i];
        while ( m->col[k] <= i ) {
            m->l->col[r] = m->col[k];
            k++;
            r++;
        }
        m->l->start_row[i+1] = r;
    }

    printf("Step one complete\n");

    /* Step.2 - Fill in column and row vectors for upper triangular matrix */

    for ( i=r=k = 0 ; i<m->n ; i++ ) {
        k = m->start_row[i];
        while ( m->col[k] < i ) k++;
        while ( k < m->start_row[i+1] ) {
            m->u->col[r] = m->col[k];
            k++;
            r++;
        }
        m->u->start_row[i+1] = r;
    }

    printf("Step two complete\n");

    /* Step.3a - Fill in matrix entries for first row of upper triangualar
                 matrix and first column of lower traingular matrix */

    for ( i=0 ; i < m->u->start_row[1] ; i++ ) m->u->a[i] = m->a[i];
    for ( i=0 ; i < m->n ; i++ ) {
        if ( me(m,i,0) != (0.0) ) m->l->a[m->l->start_row[i]] = (me(m,i,0))/(me(m,0,0));
    }

    printf("Step three A complete\n");
    /* Step.3b - Fill in remaining upper and lower triangular matrix entries */

    for ( i = 1 ; i<m->n ; i++ ) {
        for ( n = i; n<m->n ; n++ ) {
            tmp = 0.0;
            for ( k=0; k<i; k++) tmp += (me(m->l,i,k))*(me(m->u,k,n));
            tmp = me(m,i,n) - tmp;
            for (r=m->u->start_row[i]; r<m->u->start_row[i+1]; r++) if(m->u->col[r] == n) m->u->a[r] = tmp;
            tmp = 0.0;
            for (k=0 ; k<i; k++) tmp += (me(m->l,n,k))*(me(m->u,k,i));
            if (me(m->u,i,i) != (0.0) ) tmp = (me(m,n,i) - tmp)/me(m->u,i,i);
            for ( r = m->l->start_row[n] ; r<m->l->start_row[n+1] ; r++ ) {
                if( m->l->col[r] == i && i != n ) {
                    m->l->a[r] = tmp;
                } else if ( i == n ) {
                    m->l->a[m->l->start_row[i+1] - 1 ] = 1;
                }
            }
        }
    }
    return;
}

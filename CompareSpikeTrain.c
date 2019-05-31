#include <stdio.h>
#include <stdlib.h>
#include <math.h>


 /**********************************************************
        Program to compare spike timings in two trains
  **********************************************************/

/* Global definitions */
#define         INPUT1    "KenTest100.dat"
#define         INPUT2    "JayTest100.dat"


void main( void )
{
    int     j, k, time, nspk1, nspk2, *SpikeTrain1, *SpikeTrain2;
    int     common, isolated, repeat;
    FILE    *fp;

/*  Step 1 - Read first input */
    if ( (fp=fopen(INPUT1,"r")) == NULL ) {
        printf("\n Cannot find file %s", INPUT1);
        return;
    } else {
        nspk1 = 0;
        while ( fscanf(fp,"%d", &k) != EOF ) nspk1++;
        SpikeTrain1 = (int *) malloc( nspk1*sizeof(int) );
        rewind(fp);
        for ( k=0 ; k<nspk1 ; k++ ) fscanf(fp,"%d", &SpikeTrain1[k]);
        fclose(fp);
    }
    printf("\n         No. Spikes in Train 1 %d", nspk1);

/*  Step 2 - Read second input */
    if ( (fp=fopen(INPUT2,"r")) == NULL ) {
        printf("\n Cannot find file %s", INPUT2);
        return;
    } else {
        nspk2 = 0;
        while ( fscanf(fp,"%d", &k) != EOF ) nspk2++;
        SpikeTrain2 = (int *) malloc( nspk2*sizeof(int) );
        rewind(fp);
        for ( k=0 ; k<nspk2 ; k++ ) fscanf(fp,"%d", &SpikeTrain2[k]);
        fclose(fp);
    }
    printf("\n         No. Spikes in Train 2 %d", nspk2);

/*  Step 3 - Compare spike trains */
    common = isolated = 0;
    for ( k=0 ; k<nspk1 ; k++ ) {
        time = SpikeTrain1[k];
        j = 0;
        do {
            repeat = ( abs(time-SpikeTrain2[j]) >= 1 );
            j++;
        } while ( j<nspk2 && repeat );
        if ( !repeat ) {
            common++;
        } else {
            isolated++;
        }
    }

    printf("\n\n");
    printf("\n             No. Common spikes %d", common);
    printf("\nNo. Isolated spikes in Train 1 %d", isolated);
    printf("\n\n");

/*  Step 4 - Clean up */
    free(SpikeTrain1);
    free(SpikeTrain2);

    return;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

 /***************************************************************
    Function to generate exact solutions and locations of stimuli
  ***************************************************************/

/* Global definitions */
#define         OUTPUT1    "Stim.dat"
#define         OUTPUT2    "Func.dat"

/* Parameters for exact solution */
#define           NCON    75    /* Number of contacts */
#define            AMP    2.0e-5

int main( void )
{
    int k, j, ncon;
    double max;
    char word[100], txt[4];
    FILE *fp1, *fp2;

/*  Write out stimulus file */
    fp1 = fopen(OUTPUT1,"w");
    fp2 = fopen(OUTPUT2,"w");
    for ( ncon=1 ; ncon<=NCON ; ncon++ ) {

/*  Generate stimulus input */
        k = ncon;
        txt[0] = 48+k/100;
        k = k%100;
        txt[1] = 48+k/10;
        k = k%10;
        txt[2] = 48+k;
        txt[3] = '\0';
        strcpy(word,"objref stim");
        strcat(word,txt);
        fprintf(fp1,"%s\n",word);

/*  Generate stimulus - line 1 */
        strcpy(word,"dend[$1] stim");
        strcat(word,txt);
        strcat(word," = new IClamp($2)");
        fprintf(fp2,"\t%s\n",word);

/*  Generate stimulus - line 2 */
        strcpy(word,"stim");
        strcat(word,txt);
        strcat(word,".del = 0");
        fprintf(fp2,"\t%s\n",word);

/*  Generate stimulus - line 3 */
        strcpy(word,"stim");
        strcat(word,txt);
        strcat(word,".dur = 10");
        fprintf(fp2,"\t%s\n",word);

/*  Generate stimulus - line 4 */
        strcpy(word,"stim");
        strcat(word,txt);
        strcat(word,".amp = 0.02");
        fprintf(fp2,"\t%s\n\n",word);
    }
    fclose(fp1);
    fclose(fp2);
    return 0;
}

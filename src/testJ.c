/* Fonction pour calculer la trace des valeurs absolues des produits 
   entre variables de X et leur lag vector */

#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include "adesub.h"

/***********************************************************************/
double traceabsXtdLXq (double **X, double **L, double *d, double *q)
/*  Produit matriciel XtDLXQ avec LX comme lag.matrix   */
{
    /* Declarations de variables C locales */
    int j, i, lig, col;
    double **auxi, **A, trace;
    
    
    
    /* Allocation memoire pour les variables C locales */
    lig = X[0][0];
    col = X[1][0];
    taballoc(&auxi, lig, col);
    taballoc(&A, col, col);
    
    
    /* Calcul de LX */
    prodmatABC(L, X, auxi);
    
    /* Calcul de DLX */
    for (i=1;i<=lig;i++) {
        for (j=1;j<=col;j++) {
            auxi[i][j] = auxi[i][j] * d[i];
        }       
    }
    
    /* Calcul de XtDLX */
    prodmatAtBC(X,auxi,A);

    /* Trace en valeurs absolues */ 
    for(i=1;i<=col;i++){
      if(A[i][i] < 0.0) A[i][i] = -A[i][i];
    }
    
    /* Calcul de trace(XtDLXQ) */
    trace=0;
    for (i=1;i<=col;i++) {
        trace = trace + A[i][i] * q[i];
    }
    
    /* Libération des réservations locales */
    freetab (auxi);
    freetab (A);
    return(trace);
}


/* Test sur la trace des valeurs absolues de XtDLXQ 
   Complementaire du test de multispati */
void testJ (int *npermut, int *lig1, int *col1, double *tab, double *mat, double *lw, double *cw, double *inersim){ 
/* Declarations de variables C locales */
    
    int         i, j, k, lig, col, nper, *numero;
    double      **X, **L, **Xperm;
    double      *d, *q, *dperm;
    
/* Allocation memoire pour les variables C locales */
    
    nper = *npermut;
    lig = *lig1;
    col= *col1;
    
    

    taballoc(&X, lig, col);
    taballoc(&L, lig, lig);
    taballoc(&Xperm, lig, col);
    vecintalloc (&numero, lig);
    vecalloc(&dperm, lig);
    vecalloc(&d, lig);
    vecalloc(&q, col);
    
/* On recopie les objets R dans les variables C locales */

    k = 0;
    for (j=1; j<=col; j++) {
        for (i=1; i<=lig; i++) {
            X[i][j] = tab[k];
            k = k + 1;
        }
    }

    k = 0;
    for (i=1; i<=lig; i++) {
        for (j=1; j<=lig; j++) {
            L[j][i] = mat[k];
            k = k + 1;
        }
    }

    k=0;
    for (i=1; i<=lig; i++) {
        d[i]=lw[k];
        k = k + 1;
    }
    
    k=0;
    for (i=1; i<=col; i++) {
        q[i]=cw[k];
        k = k + 1;
    }
    
/* On calcul la valeur observée */
    inersim[0]=traceabsXtdLXq(X, L, d, q);
    
/* On calcul les valeurs pour chaque simulation */

    for (j=1; j<=nper; j++) {
        getpermutation(numero, j);
        matpermut(X, numero ,Xperm);
        vecpermut(d, numero, dperm);
        inersim[j]=traceabsXtdLXq(Xperm, L, dperm, q);
    }
    
    
/* Libération des réservations locales */
    
freetab(X);
freetab(L);
freetab(Xperm);
freeintvec(numero);
freevec(dperm);
freevec(d);
freevec(q);
}

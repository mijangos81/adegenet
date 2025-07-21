/* ================================================================
 GLfunctions.cpp   – adegenet (OpenMP version of GLdotProd)
 ================================================================ */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <R_ext/Print.h>
#include "snpbin.h"

#ifdef _OPENMP
#include <omp.h>
#endif

/* ----------------------------------------------------------------
 1.  GLdotProd  – all pairwise dot products (packed lower-tri ‖ diag)
 ------------------------------------------------------------------- */

/* *** FIX 1: use struct snpbin (no typedef in header) ************/
typedef double (*dotfun_t)(struct snpbin*, struct snpbin*,
                double*,        double*);

extern "C"
void GLdotProd(unsigned char *gen,
               int *nbvecperind, int *byteveclength, int *nbnaperind,
               int *naposi,      int *nind,         int *nloc,
               int *ploidy,
               double *mean,     double *sd,
               short  *freq,     double *res)
{
  /* 1. guard against near-zero variance */
  for (int i = 0; i < *nloc; ++i)
    if (sd[i] < NEARZERO) sd[i] = 1.0;
    
    /* 2. internal C representation */
    genlightC dat = genlightTogenlightC(gen, nbvecperind, byteveclength,
                                        nbnaperind, naposi, nind,
                                        nloc,      ploidy);
    
    /* *** FIX 2: dotfun_t now matches original prototypes ********/
    dotfun_t dotprod = (*freq)
      ? snpbin_dotprod_freq
    : snpbin_dotprod_int;
    
    const int    N           = *nind;
    const size_t diag_offset = (size_t)N * (N - 1) / 2;
    
#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none)      \
    shared(dat, N, mean, sd, res, dotprod)
#endif
      for (int i = 0; i < N - 1; ++i) {
        size_t base = (size_t)i * N - (size_t)i * (i + 1) / 2;
        for (int j = i + 1; j < N; ++j) {
          size_t k = base + (size_t)(j - i - 1);
          res[k]   = dotprod(&dat.x[i], &dat.x[j], mean, sd);
        }
      }
      
#ifdef _OPENMP
#pragma omp parallel for schedule(static) default(none)      \
      shared(dat, N, mean, sd, res, dotprod, diag_offset)
#endif
        for (int i = 0; i < N; ++i) {
          res[diag_offset + i] =
            dotprod(&dat.x[i], &dat.x[i], mean, sd);
        }
}

/* ----------------------------------------------------------------
 2.  GLsumInt  – per-locus integer sums (unchanged)
 ------------------------------------------------------------------- */
extern "C"
void GLsumInt(unsigned char *gen,
              int *nbvecperind, int *byteveclength, int *nbnaperind,
              int *naposi, int *nind, int *nloc, int *ploidy,
              int *res)
{
  genlightC dat = genlightTogenlightC(gen, nbvecperind, byteveclength,
                                      nbnaperind, naposi, nind,
                                      nloc,      ploidy);
  
  int *buf = (int*)calloc(*nloc, sizeof(int));
  for (int i = 0; i < *nind; ++i) {
    snpbin2intvec(&dat.x[i], buf);
    for (int j = 0; j < *nloc; ++j)
      if (!snpbin_isna(&dat.x[i], j)) res[j] += buf[j];
  }
  free(buf);
}

/* ----------------------------------------------------------------
 3.  GLsumFreq – per-locus frequency sums (unchanged)
 ------------------------------------------------------------------- */
extern "C"
void GLsumFreq(unsigned char *gen,
               int *nbvecperind, int *byteveclength, int *nbnaperind,
               int *naposi, int *nind, int *nloc, int *ploidy,
               double *res)
{
  genlightC dat = genlightTogenlightC(gen, nbvecperind, byteveclength,
                                      nbnaperind, naposi, nind,
                                      nloc,      ploidy);
  
  double *buf = (double*)calloc(*nloc, sizeof(double));
  for (int i = 0; i < *nind; ++i) {
    snpbin2freq(&dat.x[i], buf);
    for (int j = 0; j < *nloc; ++j)
      if (!snpbin_isna(&dat.x[i], j)) res[j] += buf[j];
  }
  free(buf);
}

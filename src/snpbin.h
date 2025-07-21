#ifndef SNPBIN_H
#define SNPBIN_H
/* ================================================================
 snpbin.h  –  core bit-coded genotype structures for adegenet
 ---------------------------------------------------------------
 • Provides the snpbin and genlightC data structures
 • Declares all C helpers used in *.c / *.cpp sources
 • Ensures C linkage when included from C++ code
 ================================================================ */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <stdint.h>      /* fixed-width ints if you need them   */

#define NEARZERO 1e-10
#define TRUE       1
#define FALSE      0

/* ----------------------------------------------------------------
 Compact storage for one individual’s SNPs
 ------------------------------------------------------------------- */
struct snpbin {
  unsigned char *bytevec;      /* concatenated bit vectors        */
int *byteveclength;          /* bytes per vector                */
int *bytevecnb;              /* number of byte vectors per ind. */
int *nloc;                   /* number of loci (SNPs)           */
int *nanb;                   /* length of naposi                */
int *naposi;                 /* NA positions                    */
int *ploidy;                 /* 1, 2, …                         */
};

/* ----------------------------------------------------------------
 Array-of-structs wrapper used by genlight ⇒ C conversion
 ------------------------------------------------------------------- */
struct genlightC {
  struct snpbin *x;   /* array of length *nind   */
int           *nind;
};

/* ----------------------------------------------------------------
 Public C API – make visible to C++ with unmangled names
 ------------------------------------------------------------------- */
#ifdef __cplusplus
extern "C" {
#endif
  
  /* ==== low-level bit/byte helpers ================================= */
  void  byteToBinInt   (unsigned char in, int    *out);
  void  byteToBinDouble(unsigned char in, double *out);
  void  bytesToBinInt  (unsigned char *vecbytes, int *vecsize,
                        int *vecres);
  struct snpbin makesnpbin(unsigned char *bytevec,
                           int *byteveclength, int *bytevecnb,
                           int *nloc, int *nanb,
                           int *naposi, int *ploidy);
  
  /* ==== conversion helpers (inflate / deflate) ===================== */
  void  bytesToInt   (unsigned char *vecbytes, int *veclength,
                      int *nbvec,     int *vecres,  int *reslength);
  void  bytesToDouble(unsigned char *vecbytes, int *veclength,
                      int *nbvec,     double *vecres, int *reslength);
  void  binIntToBytes(int *vecsnp, int *vecsize,
                      unsigned char *vecres, int *ressize);
  
  /* ==== snpbin “methods” =========================================== */
  int   nLoc              (struct snpbin *x);
  int   ploidy            (struct snpbin *x);
  void  snpbin2intvec     (struct snpbin *x, int    *out);
  void  snpbin2freq       (struct snpbin *x, double *out);
  void  printsnpbin       (struct snpbin *x);
  short snpbin_isna       (struct snpbin *x, int i);
  
  double snpbin_dotprod_int (struct snpbin *x, struct snpbin *y,
                             double *mean, double *sd);
  double snpbin_dotprod_freq(struct snpbin *x, struct snpbin *y,
                             double *mean, double *sd);
  
  /* ==== glue: R genlight object → C struct ========================= */
  struct genlightC genlightTogenlightC(unsigned char *gen,
                                       int *nbvecperind,   int *byteveclength,
                                       int *nbnaperind,    int *naposi,
                                       int *nind,          int *nloc,
                                       int *ploidy);
  
  /* ==== test helpers =============================================== */
  void testRaw(unsigned char *a, int *n);
  
#ifdef __cplusplus
}   /* extern "C" */
#endif
  
#endif  /* SNPBIN_H */
  
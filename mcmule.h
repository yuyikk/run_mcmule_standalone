extern "C" {
#ifndef MCMULE_H
#define MCMULE_H
#ifdef __GFORTRAN__
#define __OPENLOOPS_BUILD_DIR ""
#define __OPENLOOPS_INSTALL_DIR "/builds/mule-tools/mcmule/usr/lib64"
#else
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef double real;
typedef uint32_t BOOL;

typedef void (*quantfunc)(real**,real*, real*, real*, real*, real*, real*, real*);
typedef void (*usereventfunc)(real*, int*);
typedef void (*inituserfunc)();

void mcmule_set_observable(int number_hist, int number_bins, real *lower_bounds, real *upper_bounds,
      quantfunc measurement_function, usereventfunc user_integration, inituserfunc user_initialisation,
      int user_integration_dimension, int names_length, int filenamesuffix_length);

void mcmule_initflavour(char *flavour, double *scms);

void mcmule_runmcmule(int ncall_ad, int itmx_ad, int ncall, int itmx,
    int initial_ran_seed, real xicut1, real xicut2,
    char *piece, char *tfilename);

typedef void (*proton_ff)(real*q2, real*Ge, real*Gm);
typedef void (*scalar_ff)(real*q2, real*Sff);

extern proton_ff __nucl_protonff_MOD_proton_dipole;
extern proton_ff __nucl_protonff_MOD_proton_monopole;

extern proton_ff __nucl_protonff_MOD_the_proton_ff;
#define mcmule_protonff __nucl_protonff_MOD_the_proton_ff
extern scalar_ff __nucl_protonff_MOD_the_scalar_ff;
#define mcmule_scalarff __nucl_protonff_MOD_the_scalar_ff
extern real __nucl_protonff_MOD_kappa;
#define mcmule_protonff_kappa __nucl_protonff_MOD_kappa
extern real __nucl_protonff_MOD_lambda;
#define mcmule_protonff_lambda __nucl_protonff_MOD_lambda

extern char __global_def_MOD_hvpmodel[11];
#define mcmule_hvpmodel __global_def_MOD_hvpmodel

extern real __global_def_MOD_gf;
#define mcmule_Gf __global_def_MOD_gf
extern real __global_def_MOD_alpha;
#define mcmule_alpha __global_def_MOD_alpha
extern real __global_def_MOD_sw2;
#define mcmule_sw2 __global_def_MOD_sw2
extern real __global_def_MOD_mm;
#define mcmule_Mm __global_def_MOD_mm
extern real __global_def_MOD_me;
#define mcmule_Me __global_def_MOD_me
extern real __global_def_MOD_mt;
#define mcmule_Mtau __global_def_MOD_mt
extern real __global_def_MOD_scms;
#define mcmule_scms __global_def_MOD_scms
extern real __global_def_MOD_clj;
#define mcmule_CLj __global_def_MOD_clj
extern real __global_def_MOD_crj;
#define mcmule_CRj __global_def_MOD_crj
extern real __global_def_MOD_mj;
#define mcmule_Mj __global_def_MOD_mj
extern real *__global_def_MOD_pol1;
#define mcmule_pol1 __global_def_MOD_pol1
extern real *__global_def_MOD_pol2;
#define mcmule_pol2 __global_def_MOD_pol2
extern real __global_def_MOD_musq;
#define mcmule_musq __global_def_MOD_musq
extern int __global_def_MOD_nel;
#define mcmule_Nel __global_def_MOD_nel
extern int __global_def_MOD_nmu;
#define mcmule_Nmu __global_def_MOD_nmu
extern int __global_def_MOD_ntau;
#define mcmule_Ntau __global_def_MOD_ntau
extern int __global_def_MOD_nhad;
#define mcmule_Nhad __global_def_MOD_nhad
extern BOOL *__user_dummy_MOD_pass_cut;
#define mcmule_pass_cut __user_dummy_MOD_pass_cut
extern real __user_dummy_MOD_userweight;
#define mcmule_userweight __user_dummy_MOD_userweight
extern char *__user_dummy_MOD_names;
#define mcmule_names __user_dummy_MOD_names
extern char *__user_dummy_MOD_filenamesuffix;
#define mcmule_filenamesuffix __user_dummy_MOD_filenamesuffix
extern int __user_dummy_MOD_nr_q;
#define mcmule_nr_q __user_dummy_MOD_nr_q
extern size_t mcmule_namelength;

#define MCMULE_SET_NAME(q, name) strcpy(mcmule_names+q*mcmule_namelength, name);

#endif
#endif
}

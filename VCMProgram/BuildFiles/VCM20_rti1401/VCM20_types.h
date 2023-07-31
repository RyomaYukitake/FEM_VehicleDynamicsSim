/*
 * VCM20_types.h
 *
 * Sponsored License - for use in support of a program or activity
 * sponsored by MathWorks.  Not for government, commercial or other
 * non-sponsored organizational use.
 *
 * Code generation for model "VCM20".
 *
 * Model version              : 8.145
 * Simulink Coder version : 9.6 (R2021b) 14-May-2021
 * C source code generated on : Mon Jul 31 14:28:36 2023
 *
 * Target selection: rti1401.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Custom Processor->Custom
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_VCM20_types_h_
#define RTW_HEADER_VCM20_types_h_
#include "rtwtypes.h"
#include "multiword_types.h"

/* Model Code Variants */
#ifndef struct_tag_bLXxAw6j3QavLR7OIJgfL
#define struct_tag_bLXxAw6j3QavLR7OIJgfL

struct tag_bLXxAw6j3QavLR7OIJgfL
{
  int32_T isInitialized;
  boolean_T isSetupComplete;
  real_T pCumSum;
  real_T pCumSumRev[999];
  real_T pCumRevIndex;
  real_T pModValueRev;
};

#endif                                 /* struct_tag_bLXxAw6j3QavLR7OIJgfL */

#ifndef typedef_g_dsp_internal_SlidingWindowA_T
#define typedef_g_dsp_internal_SlidingWindowA_T

typedef struct tag_bLXxAw6j3QavLR7OIJgfL g_dsp_internal_SlidingWindowA_T;

#endif                             /* typedef_g_dsp_internal_SlidingWindowA_T */

#ifndef struct_tag_BlgwLpgj2bjudmbmVKWwDE
#define struct_tag_BlgwLpgj2bjudmbmVKWwDE

struct tag_BlgwLpgj2bjudmbmVKWwDE
{
  uint32_T f1[8];
};

#endif                                 /* struct_tag_BlgwLpgj2bjudmbmVKWwDE */

#ifndef typedef_cell_wrap_VCM20_T
#define typedef_cell_wrap_VCM20_T

typedef struct tag_BlgwLpgj2bjudmbmVKWwDE cell_wrap_VCM20_T;

#endif                                 /* typedef_cell_wrap_VCM20_T */

#ifndef struct_tag_clciNXyUID92GYviiGZ0AE
#define struct_tag_clciNXyUID92GYviiGZ0AE

struct tag_clciNXyUID92GYviiGZ0AE
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  boolean_T TunablePropsChanged;
  cell_wrap_VCM20_T inputVarSize;
  g_dsp_internal_SlidingWindowA_T *pStatistic;
  int32_T NumChannels;
  g_dsp_internal_SlidingWindowA_T _pobj0;
};

#endif                                 /* struct_tag_clciNXyUID92GYviiGZ0AE */

#ifndef typedef_dsp_simulink_MovingAverage_VC_T
#define typedef_dsp_simulink_MovingAverage_VC_T

typedef struct tag_clciNXyUID92GYviiGZ0AE dsp_simulink_MovingAverage_VC_T;

#endif                             /* typedef_dsp_simulink_MovingAverage_VC_T */

#ifndef struct_tag_4CTZbzBcRF4SeaOibUrdpG
#define struct_tag_4CTZbzBcRF4SeaOibUrdpG

struct tag_4CTZbzBcRF4SeaOibUrdpG
{
  int32_T isInitialized;
  boolean_T isSetupComplete;
  real32_T pCumSum;
  real32_T pCumSumRev[199];
  real32_T pCumRevIndex;
  real32_T pModValueRev;
};

#endif                                 /* struct_tag_4CTZbzBcRF4SeaOibUrdpG */

#ifndef typedef_g_dsp_internal_SlidingWindo_b_T
#define typedef_g_dsp_internal_SlidingWindo_b_T

typedef struct tag_4CTZbzBcRF4SeaOibUrdpG g_dsp_internal_SlidingWindo_b_T;

#endif                             /* typedef_g_dsp_internal_SlidingWindo_b_T */

#ifndef struct_tag_RbQaF1XgHKtTPQpya6FBpG
#define struct_tag_RbQaF1XgHKtTPQpya6FBpG

struct tag_RbQaF1XgHKtTPQpya6FBpG
{
  boolean_T matlabCodegenIsDeleted;
  int32_T isInitialized;
  boolean_T isSetupComplete;
  boolean_T TunablePropsChanged;
  cell_wrap_VCM20_T inputVarSize;
  g_dsp_internal_SlidingWindo_b_T *pStatistic;
  int32_T NumChannels;
  g_dsp_internal_SlidingWindo_b_T _pobj0;
};

#endif                                 /* struct_tag_RbQaF1XgHKtTPQpya6FBpG */

#ifndef typedef_dsp_simulink_MovingAverage_b_T
#define typedef_dsp_simulink_MovingAverage_b_T

typedef struct tag_RbQaF1XgHKtTPQpya6FBpG dsp_simulink_MovingAverage_b_T;

#endif                              /* typedef_dsp_simulink_MovingAverage_b_T */

/* Parameters (default storage) */
typedef struct P_VCM20_T_ P_VCM20_T;

/* Forward declaration for rtModel */
typedef struct tag_RTM_VCM20_T RT_MODEL_VCM20_T;

#endif                                 /* RTW_HEADER_VCM20_types_h_ */

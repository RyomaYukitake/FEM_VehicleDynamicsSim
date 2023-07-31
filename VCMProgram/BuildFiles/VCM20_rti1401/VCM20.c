/*
 * VCM20.c
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

#include "VCM20_trc_ptr.h"
#include "VCM20.h"
#include "VCM20_private.h"

/* Named constants for Chart: '<S172>/Chart' */
#define VCM20_IN_Base                  ((uint8_T)1U)
#define VCM20_IN_LaunchGearProtection  ((uint8_T)2U)
#define VCM20_IN_LaunchLEDoff          ((uint8_T)1U)
#define VCM20_IN_LaunchLEDon           ((uint8_T)2U)
#define VCM20_IN_LaunchReady           ((uint8_T)3U)
#define VCM20_IN_LaunchReady_BrakeON   ((uint8_T)4U)
#define VCM20_IN_LaunchReady_BrakeOff  ((uint8_T)5U)
#define VCM20_IN_NO_ACTIVE_CHILD       ((uint8_T)0U)
#define VCM20_IN_OFF                   ((uint8_T)6U)

/* Named constants for Chart: '<S16>/Steer Chart' */
#define VCM20_IN_Display1              ((uint8_T)1U)
#define VCM20_IN_Gain1                 ((uint8_T)2U)
#define VCM20_IN_Gain2                 ((uint8_T)3U)
#define VCM20_IN_Minus_Gain1           ((uint8_T)4U)
#define VCM20_IN_Minus_Gain2           ((uint8_T)5U)
#define VCM20_IN_Plus_Gain1            ((uint8_T)6U)
#define VCM20_IN_Plus_Gain2            ((uint8_T)7U)
#define VCM20_IN_Select1               ((uint8_T)8U)
#define VCM20_IN_Select2               ((uint8_T)9U)
#define VCM20_IN_Select3               ((uint8_T)10U)
#define VCM20_IN_Start                 ((uint8_T)11U)

/* Block signals (default storage) */
B_VCM20_T VCM20_B;

/* Continuous states */
X_VCM20_T VCM20_X;

/* Block states (default storage) */
DW_VCM20_T VCM20_DW;

/* Real-time model */
static RT_MODEL_VCM20_T VCM20_M_;
RT_MODEL_VCM20_T *const VCM20_M = &VCM20_M_;
static void rate_scheduler(void);

/* n-D Spline interpolation function */
real_T look_SplNBinSZcd(uint32_T numDims, const real_T* u, const
  rt_LUTSplineWork * const SWork)
{
  /*
   *   n-D column-major table lookup operating on real_T with:
   *       - Spline interpolation
   *       - Spline extrapolation
   *       - Binary breakpoint search
   *       - Index search starts at the same place each time
   */
  rt_LUTnWork * const TWork_look = SWork->m_TWork;
  real_T* const fraction = (real_T*) TWork_look->m_bpLambda;
  uint32_T* const bpIdx = TWork_look->m_bpIndex;
  const uint32_T* const maxIndex = TWork_look->m_maxIndex;
  uint32_T k;
  for (k = 0U; k < numDims; k++) {
    const real_T* const bpData = ((const real_T * const *)
      TWork_look->m_bpDataSet)[k];
    bpIdx[k] = plook_binx(u[k], bpData, maxIndex[k], &fraction[k]);
  }

  return(intrp_NSplcd(numDims, SWork, 3U));
}

/*
 * Second derivative initialization function for spline
 * for last dimension.
 */
void rt_Spline2Derivd(const real_T *x, const real_T *y, uint32_T n, real_T *u,
                      real_T *y2)
{
  real_T p, qn, sig, un;
  uint32_T n1, i, k;
  n1 = n - 1U;
  y2[0U] = 0.0;
  u[0U] = 0.0;
  for (i = 1U; i < n1; i++) {
    real_T dxm1 = x[i] - x[i - 1U];
    real_T dxp1 = x[i + 1U] - x[i];
    real_T dxpm = dxp1 + dxm1;
    sig = dxm1 / dxpm;
    p = (sig * y2[i - 1U]) + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = ((y[i + 1U] - y[i]) / dxp1) - ((y[i] - y[i - 1U]) / dxm1);
    u[i] = (((6.0 * u[i]) / dxpm) - (sig * u[i - 1U])) / p;
  }

  qn = 0.0;
  un = 0.0;
  y2[n1] = (un - (qn * u[n1 - 1U])) / ((qn * y2[n1 - 1U]) + 1.0);
  for (k = n1; k > 0U; k--) {
    y2[k-1U] = (y2[k-1U] * y2[k]) + u[k-1U];
  }

  return;
}

/* n-D natural spline calculation function */
real_T intrp_NSplcd(uint32_T numDims, const rt_LUTSplineWork * const splWork,
                    uint32_T extrapMethod)
{
  uint32_T il;
  uint32_T iu, k, i;
  real_T h, s, p, smsq, pmsq;

  /* intermediate results work areas "this" and "next" */
  const rt_LUTnWork *TWork_interp = (const rt_LUTnWork *)splWork->m_TWork;
  const real_T *fraction = (real_T *) TWork_interp->m_bpLambda;
  const real_T *yp = (real_T *) TWork_interp->m_tableData;
  real_T *yyA = (real_T *) splWork->m_yyA;
  real_T *yyB = (real_T *) splWork->m_yyB;
  real_T *yy2 = (real_T *) splWork->m_yy2;
  real_T *up = (real_T *) splWork->m_up;
  real_T *y2 = (real_T *) splWork->m_y2;
  uint8_T* reCalc = splWork->m_reCalc;
  real_T *dp = (real_T *) splWork->m_preBp0AndTable;
  const real_T **bpDataSet = (const real_T **) TWork_interp->m_bpDataSet;
  const real_T *xp = bpDataSet[0U];
  real_T *yy = yyA;
  uint32_T bufBank = 0U;
  uint32_T len = TWork_interp->m_maxIndex[0U] + 1U;

  /* compare bp0 and table to see whether they get changed */
  {
    /* compare the bp0 data */
    if (memcmp(dp, xp,
               len * sizeof(real_T)) != 0) {
      *reCalc = 1;
      (void) memcpy(dp, xp,
                    len * sizeof(real_T));
    }

    /* compare the table data */
    dp = &(dp[len]);
    if (memcmp(dp, yp,
               len * splWork->m_numYWorkElts[0U] * sizeof(real_T)) != 0) {
      *reCalc = 1;
      (void) memcpy(dp, yp,
                    len * splWork->m_numYWorkElts[0U] * sizeof(real_T));
    }
  }

  if (*reCalc == 1) {
    /* If table and bps are tunable calculate 1st dim 2nd deriv */
    /* Generate first dimension's second derivatives */
    for (i = 0U; i < splWork->m_numYWorkElts[0U]; i++) {
      rt_Spline2Derivd(xp, yp, len, up, y2);
      yp = &yp[len];
      y2 = &y2[len];
    }

    /* Set pointers back to beginning */
    yp = (const real_T *) TWork_interp->m_tableData;
    y2 = (real_T *) splWork->m_y2;
  }

  *reCalc = 0;

  /* Generate at-point splines in each dimension */
  for (k = 0U; k < numDims; k++ ) {
    /* this dimension's input setup */
    xp = bpDataSet[k];
    len = TWork_interp->m_maxIndex[k] + 1U;
    il = TWork_interp->m_bpIndex[k];
    iu = il + 1U;
    h = xp[iu] - xp[il];
    p = fraction[k];
    s = 1.0 - p;
    pmsq = p * ((p*p) - 1.0);
    smsq = s * ((s*s) - 1.0);

    /*
     * Calculate spline curves for input in this
     * dimension at each value of the higher
     * other dimensions\' points in the table.
     */
    if ((p > 1.0) && (extrapMethod == 2U) ) {
      real_T slope;
      for (i = 0U; i < splWork->m_numYWorkElts[k]; i++) {
        slope = (yp[iu] - yp[il]) + ((y2[il]*h*h)*(1.0/6.0));
        yy[i] = yp[iu] + (slope * (p-1.0));
        yp = &yp[len];
        y2 = &y2[len];
      }
    } else if ((p < 0.0) && (extrapMethod == 2U) ) {
      real_T slope;
      for (i = 0U; i < splWork->m_numYWorkElts[k]; i++) {
        slope = (yp[iu] - yp[il]) - ((y2[iu]*h*h)*(1.0/6.0));
        yy[i] = yp[il] + (slope * p);
        yp = &yp[len];
        y2 = &y2[len];
      }
    } else {
      for (i = 0U; i < splWork->m_numYWorkElts[k]; i++) {
        yy[i] = yp[il] + p * (yp[iu] - yp[il]) +
          ((smsq * y2[il] + pmsq * y2[iu])*h*h)*(1.0/6.0);
        yp = &yp[len];
        y2 = &y2[len];
      }
    }

    /* set pointers to new result and calculate second derivatives */
    yp = yy;
    y2 = yy2;
    if (splWork->m_numYWorkElts[k+1U] > 0U ) {
      uint32_T nextLen = TWork_interp->m_maxIndex[k+1U] + 1U;
      const real_T *nextXp = bpDataSet[k+1U];
      for (i = 0U; i < splWork->m_numYWorkElts[k+1U]; i++) {
        rt_Spline2Derivd(nextXp, yp, nextLen, up, y2);
        yp = &yp[nextLen];
        y2 = &y2[nextLen];
      }
    }

    /*
     * Set work vectors yp, y2 and yy for next iteration;
     * the yy just calculated becomes the yp in the
     * next iteration, y2 was just calculated for these
     * new points and the yy buffer is swapped to the space
     * for storing the next iteration\'s results.
     */
    yp = yy;
    y2 = yy2;

    /*
     * Swap buffers for next dimension and
     * toggle bufBank for next iteration.
     */
    if (bufBank == 0U) {
      yy = yyA;
      bufBank = 1U;
    } else {
      yy = yyB;
      bufBank = 0U;
    }
  }

  return( yp[0U] );
}

real_T look1_binlxpw(real_T u0, const real_T bp0[], const real_T table[],
                     uint32_T maxIndex)
{
  real_T frac;
  real_T yL_0d0;
  uint32_T bpIdx;
  uint32_T iLeft;
  uint32_T iRght;

  /* Column-major Lookup 1-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = (u0 - bp0[0U]) / (bp0[1U] - bp0[0U]);
  } else if (u0 < bp0[maxIndex]) {
    /* Binary Search */
    bpIdx = maxIndex >> 1U;
    iLeft = 0U;
    iRght = maxIndex;
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex - 1U;
    frac = (u0 - bp0[maxIndex - 1U]) / (bp0[maxIndex] - bp0[maxIndex - 1U]);
  }

  /* Column-major Interpolation 1-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  yL_0d0 = table[iLeft];
  return (table[iLeft + 1U] - yL_0d0) * frac + yL_0d0;
}

real_T look2_binlxpw(real_T u0, real_T u1, const real_T bp0[], const real_T bp1[],
                     const real_T table[], const uint32_T maxIndex[], uint32_T
                     stride)
{
  real_T fractions[2];
  real_T frac;
  real_T yL_0d0;
  real_T yL_0d1;
  uint32_T bpIndices[2];
  uint32_T bpIdx;
  uint32_T iLeft;
  uint32_T iRght;

  /* Column-major Lookup 2-D
     Search method: 'binary'
     Use previous index: 'off'
     Interpolation method: 'Linear point-slope'
     Extrapolation method: 'Linear'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u0 <= bp0[0U]) {
    iLeft = 0U;
    frac = (u0 - bp0[0U]) / (bp0[1U] - bp0[0U]);
  } else if (u0 < bp0[maxIndex[0U]]) {
    /* Binary Search */
    bpIdx = maxIndex[0U] >> 1U;
    iLeft = 0U;
    iRght = maxIndex[0U];
    while (iRght - iLeft > 1U) {
      if (u0 < bp0[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u0 - bp0[iLeft]) / (bp0[iLeft + 1U] - bp0[iLeft]);
  } else {
    iLeft = maxIndex[0U] - 1U;
    frac = (u0 - bp0[maxIndex[0U] - 1U]) / (bp0[maxIndex[0U]] - bp0[maxIndex[0U]
      - 1U]);
  }

  fractions[0U] = frac;
  bpIndices[0U] = iLeft;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u1 <= bp1[0U]) {
    iLeft = 0U;
    frac = (u1 - bp1[0U]) / (bp1[1U] - bp1[0U]);
  } else if (u1 < bp1[maxIndex[1U]]) {
    /* Binary Search */
    bpIdx = maxIndex[1U] >> 1U;
    iLeft = 0U;
    iRght = maxIndex[1U];
    while (iRght - iLeft > 1U) {
      if (u1 < bp1[bpIdx]) {
        iRght = bpIdx;
      } else {
        iLeft = bpIdx;
      }

      bpIdx = (iRght + iLeft) >> 1U;
    }

    frac = (u1 - bp1[iLeft]) / (bp1[iLeft + 1U] - bp1[iLeft]);
  } else {
    iLeft = maxIndex[1U] - 1U;
    frac = (u1 - bp1[maxIndex[1U] - 1U]) / (bp1[maxIndex[1U]] - bp1[maxIndex[1U]
      - 1U]);
  }

  /* Column-major Interpolation 2-D
     Interpolation method: 'Linear point-slope'
     Use last breakpoint for index at or above upper limit: 'off'
     Overflow mode: 'portable wrapping'
   */
  bpIdx = iLeft * stride + bpIndices[0U];
  yL_0d0 = table[bpIdx];
  yL_0d0 += (table[bpIdx + 1U] - yL_0d0) * fractions[0U];
  bpIdx += stride;
  yL_0d1 = table[bpIdx];
  return (((table[bpIdx + 1U] - yL_0d1) * fractions[0U] + yL_0d1) - yL_0d0) *
    frac + yL_0d0;
}

uint32_T plook_binx(real_T u, const real_T bp[], uint32_T maxIndex, real_T
                    *fraction)
{
  uint32_T bpIndex;

  /* Prelookup - Index and Fraction
     Index Search method: 'binary'
     Extrapolation method: 'Linear'
     Use previous index: 'off'
     Use last breakpoint for index at or above upper limit: 'off'
     Remove protection against out-of-range input in generated code: 'off'
   */
  if (u <= bp[0U]) {
    bpIndex = 0U;
    *fraction = (u - bp[0U]) / (bp[1U] - bp[0U]);
  } else if (u < bp[maxIndex]) {
    bpIndex = binsearch_u32d(u, bp, maxIndex >> 1U, maxIndex);
    *fraction = (u - bp[bpIndex]) / (bp[bpIndex + 1U] - bp[bpIndex]);
  } else {
    bpIndex = maxIndex - 1U;
    *fraction = (u - bp[maxIndex - 1U]) / (bp[maxIndex] - bp[maxIndex - 1U]);
  }

  return bpIndex;
}

uint32_T binsearch_u32d(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex)
{
  uint32_T bpIdx;
  uint32_T bpIndex;
  uint32_T iRght;

  /* Binary Search */
  bpIdx = startIndex;
  bpIndex = 0U;
  iRght = maxIndex;
  while (iRght - bpIndex > 1U) {
    if (u < bp[bpIdx]) {
      iRght = bpIdx;
    } else {
      bpIndex = bpIdx;
    }

    bpIdx = (iRght + bpIndex) >> 1U;
  }

  return bpIndex;
}

/*
 *         This function updates active task flag for each subrate.
 *         The function is called at model base rate, hence the
 *         generated code self-manages all its subrates.
 */
static void rate_scheduler(void)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (VCM20_M->Timing.TaskCounters.TID[2])++;
  if ((VCM20_M->Timing.TaskCounters.TID[2]) > 9) {/* Sample time: [0.01s, 0.0s] */
    VCM20_M->Timing.TaskCounters.TID[2] = 0;
  }

  (VCM20_M->Timing.TaskCounters.TID[3])++;
  if ((VCM20_M->Timing.TaskCounters.TID[3]) > 99) {/* Sample time: [0.1s, 0.0s] */
    VCM20_M->Timing.TaskCounters.TID[3] = 0;
  }
}

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 8;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  VCM20_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* System initialize for atomic system: */
void VCM20_MovingAverage_Init(DW_MovingAverage_VCM20_T *localDW)
{
  dsp_simulink_MovingAverage_VC_T *obj;
  g_dsp_internal_SlidingWindowA_T *obj_0;
  int32_T i;

  /* InitializeConditions for MATLABSystem: '<S70>/Moving Average' */
  obj = &localDW->obj;
  obj_0 = obj->pStatistic;
  if (obj_0->isInitialized == 1) {
    obj_0->pCumSum = 0.0;
    for (i = 0; i < 999; i++) {
      obj_0->pCumSumRev[i] = 0.0;
    }

    obj_0->pCumRevIndex = 1.0;
    obj_0->pModValueRev = 0.0;
  }

  /* End of InitializeConditions for MATLABSystem: '<S70>/Moving Average' */
}

/* Start for atomic system: */
void VCM20_MovingAverage_Start(DW_MovingAverage_VCM20_T *localDW)
{
  dsp_simulink_MovingAverage_VC_T *b_obj;
  dsp_simulink_MovingAverage_VC_T *obj;
  g_dsp_internal_SlidingWindowA_T *iobj_0;

  /* Start for MATLABSystem: '<S70>/Moving Average' */
  localDW->obj.matlabCodegenIsDeleted = true;
  b_obj = &localDW->obj;
  b_obj->isInitialized = 0;
  b_obj->NumChannels = -1;
  b_obj->matlabCodegenIsDeleted = false;
  localDW->objisempty = true;
  b_obj = &localDW->obj;
  b_obj->isSetupComplete = false;
  b_obj->isInitialized = 1;
  obj = b_obj;
  obj->NumChannels = 1;
  iobj_0 = &obj->_pobj0;
  iobj_0->isInitialized = 0;
  iobj_0->isInitialized = 0;
  obj->pStatistic = iobj_0;
  b_obj->isSetupComplete = true;
  b_obj->TunablePropsChanged = false;
}

/* Output and update for atomic system: */
void VCM20_MovingAverage(real_T rtu_0, B_MovingAverage_VCM20_T *localB,
  DW_MovingAverage_VCM20_T *localDW)
{
  dsp_simulink_MovingAverage_VC_T *obj;
  dsp_simulink_MovingAverage_VC_T *obj_0;
  g_dsp_internal_SlidingWindowA_T *obj_1;
  g_dsp_internal_SlidingWindowA_T *obj_2;
  g_dsp_internal_SlidingWindowA_T *obj_3;
  g_dsp_internal_SlidingWindowA_T *obj_4;
  real_T csum;
  real_T cumRevIndex;
  real_T modValueRev;
  real_T tmp;
  real_T z;
  int32_T i;

  /* MATLABSystem: '<S70>/Moving Average' */
  obj = &localDW->obj;
  obj_0 = obj;
  if (obj_0->TunablePropsChanged) {
    obj_0->TunablePropsChanged = false;
  }

  obj_1 = obj->pStatistic;
  if (obj_1->isInitialized != 1) {
    obj_2 = obj_1;
    obj_3 = obj_2;
    obj_3->isSetupComplete = false;
    obj_3->isInitialized = 1;
    obj_4 = obj_3;
    obj_4->pCumSum = 0.0;
    for (i = 0; i < 999; i++) {
      obj_4->pCumSumRev[i] = 0.0;
    }

    obj_4->pCumRevIndex = 1.0;
    obj_4->pModValueRev = 0.0;
    obj_3->isSetupComplete = true;
    obj_2->pCumSum = 0.0;
    for (i = 0; i < 999; i++) {
      obj_2->pCumSumRev[i] = 0.0;
    }

    obj_2->pCumRevIndex = 1.0;
    obj_2->pModValueRev = 0.0;
  }

  cumRevIndex = obj_1->pCumRevIndex;
  csum = obj_1->pCumSum;
  for (i = 0; i < 999; i++) {
    localB->csumrev[i] = obj_1->pCumSumRev[i];
  }

  modValueRev = obj_1->pModValueRev;
  z = 0.0;
  tmp = 0.0;
  csum += rtu_0;
  if (modValueRev == 0.0) {
    z = localB->csumrev[(int32_T)cumRevIndex - 1] + csum;
  }

  localB->csumrev[(int32_T)cumRevIndex - 1] = rtu_0;
  if (cumRevIndex != 999.0) {
    cumRevIndex++;
  } else {
    cumRevIndex = 1.0;
    csum = 0.0;
    for (i = 997; i >= 0; i--) {
      localB->csumrev[i] += localB->csumrev[i + 1];
    }
  }

  if (modValueRev == 0.0) {
    tmp = z / 1000.0;
  }

  if (modValueRev > 0.0) {
    modValueRev--;
  } else {
    modValueRev = 0.0;
  }

  obj_1->pCumSum = csum;
  for (i = 0; i < 999; i++) {
    obj_1->pCumSumRev[i] = localB->csumrev[i];
  }

  obj_1->pCumRevIndex = cumRevIndex;
  obj_1->pModValueRev = modValueRev;

  /* MATLABSystem: '<S70>/Moving Average' */
  localB->MovingAverage = tmp;
}

/* Termination for atomic system: */
void VCM20_MovingAverage_Term(DW_MovingAverage_VCM20_T *localDW)
{
  dsp_simulink_MovingAverage_VC_T *obj;
  g_dsp_internal_SlidingWindowA_T *obj_0;

  /* Terminate for MATLABSystem: '<S70>/Moving Average' */
  obj = &localDW->obj;
  if (!obj->matlabCodegenIsDeleted) {
    obj->matlabCodegenIsDeleted = true;
    if ((obj->isInitialized == 1) && obj->isSetupComplete) {
      obj_0 = obj->pStatistic;
      if (obj_0->isInitialized == 1) {
        obj_0->isInitialized = 2;
      }

      obj->NumChannels = -1;
    }
  }

  /* End of Terminate for MATLABSystem: '<S70>/Moving Average' */
}

real_T rt_roundd_snf(real_T u)
{
  real_T y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

uint32_T MWDSP_EPH_R_D(real_T evt, uint32_T *sta)
{
  int32_T curState;
  int32_T lastzcevent;
  int32_T newState;
  int32_T newStateR;
  uint32_T previousState;
  uint32_T retVal;

  /* S-Function (sdspcount2): '<S114>/Counter' */
  /* Detect rising edge events */
  previousState = *sta;
  retVal = 0U;
  lastzcevent = 0;
  newState = 5;
  newStateR = 5;
  if (evt > 0.0) {
    curState = 2;
  } else {
    curState = !(evt < 0.0);
  }

  if (previousState == 5U) {
    newStateR = curState;
  } else if ((uint32_T)curState != previousState) {
    if (previousState == 3U) {
      if ((uint32_T)curState == 1U) {
        newStateR = 1;
      } else {
        lastzcevent = 2;
        previousState = 1U;
      }
    }

    if (previousState == 4U) {
      if ((uint32_T)curState == 1U) {
        newStateR = 1;
      } else {
        lastzcevent = 3;
        previousState = 1U;
      }
    }

    if ((previousState == 1U) && ((uint32_T)curState == 2U)) {
      retVal = 2U;
    }

    if (previousState == 0U) {
      retVal = 2U;
    }

    if (retVal == (uint32_T)lastzcevent) {
      retVal = 0U;
    }

    if (((uint32_T)curState == 1U) && (retVal == 2U)) {
      newState = 3;
    } else {
      newState = curState;
    }
  }

  if ((uint32_T)newStateR != 5U) {
    *sta = (uint32_T)newStateR;
    retVal = 0U;
  }

  if ((uint32_T)newState != 5U) {
    *sta = (uint32_T)newState;
  }

  /* End of S-Function (sdspcount2): '<S114>/Counter' */
  return retVal;
}

uint32_T MWDSP_EPH_R_B(boolean_T evt, uint32_T *sta)
{
  int32_T curState;
  int32_T lastzcevent;
  int32_T newState;
  int32_T newStateR;
  uint32_T previousState;
  uint32_T retVal;

  /* S-Function (sdspcount2): '<S114>/Counter' */
  /* Detect rising edge events */
  previousState = *sta;
  retVal = 0U;
  lastzcevent = 0;
  newState = 5;
  newStateR = 5;
  if (evt) {
    curState = 2;
  } else {
    curState = 1;
  }

  if (previousState == 5U) {
    newStateR = curState;
  } else if ((uint32_T)curState != previousState) {
    if (previousState == 3U) {
      if ((uint32_T)curState == 1U) {
        newStateR = 1;
      } else {
        lastzcevent = 2;
        previousState = 1U;
      }
    }

    if (previousState == 4U) {
      if ((uint32_T)curState == 1U) {
        newStateR = 1;
      } else {
        lastzcevent = 3;
        previousState = 1U;
      }
    }

    if ((previousState == 1U) && ((uint32_T)curState == 2U)) {
      retVal = 2U;
    }

    if (previousState == 0U) {
      retVal = 2U;
    }

    if (retVal == (uint32_T)lastzcevent) {
      retVal = 0U;
    }

    if (((uint32_T)curState == 1U) && (retVal == 2U)) {
      newState = 3;
    } else {
      newState = curState;
    }
  }

  if ((uint32_T)newStateR != 5U) {
    *sta = (uint32_T)newStateR;
    retVal = 0U;
  }

  if ((uint32_T)newState != 5U) {
    *sta = (uint32_T)newState;
  }

  /* End of S-Function (sdspcount2): '<S114>/Counter' */
  return retVal;
}

/* Model output function */
void VCM20_output(void)
{
  dsp_simulink_MovingAverage_b_T *obj;
  dsp_simulink_MovingAverage_b_T *obj_0;
  g_dsp_internal_SlidingWindo_b_T *obj_1;
  g_dsp_internal_SlidingWindo_b_T *obj_2;
  g_dsp_internal_SlidingWindo_b_T *obj_3;
  g_dsp_internal_SlidingWindo_b_T *obj_4;
  real_T u0_0;
  real_T u1;
  real_T u1_0;
  real_T u2;
  int32_T i;
  real32_T csum;
  real32_T cumRevIndex;
  real32_T modValueRev;
  real32_T tmp;
  real32_T u0;
  real32_T z;
  boolean_T out;
  static const uint8_T b[8] = { 73U, 77U, 85U, 83U, 116U, 97U, 116U, 32U };

  static const uint8_T c[8] = { 71U, 97U, 105U, 110U, 49U, 32U, 32U, 32U };

  static const uint8_T d[8] = { 71U, 97U, 105U, 110U, 50U, 32U, 32U, 32U };

  if (rtmIsMajorTimeStep(VCM20_M)) {
    /* set solver stop time */
    if (!(VCM20_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&VCM20_M->solverInfo, ((VCM20_M->Timing.clockTickH0
        + 1) * VCM20_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&VCM20_M->solverInfo, ((VCM20_M->Timing.clockTick0 +
        1) * VCM20_M->Timing.stepSize0 + VCM20_M->Timing.clockTickH0 *
        VCM20_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(VCM20_M)) {
    VCM20_M->Timing.t[0] = rtsiGetT(&VCM20_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<Root>/Saturation' incorporates:
     *  Constant: '<Root>/MaxTorque'
     */
    u0_0 = VCM20_P.MaxTorque_Value;
    u1 = VCM20_P.Saturation_LowerSat_f;
    u2 = VCM20_P.Saturation_UpperSat_n;
    if (u0_0 > u2) {
      /* Saturate: '<Root>/Saturation' */
      VCM20_B.Saturation = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<Root>/Saturation' */
      VCM20_B.Saturation = u1;
    } else {
      /* Saturate: '<Root>/Saturation' */
      VCM20_B.Saturation = u0_0;
    }

    /* End of Saturate: '<Root>/Saturation' */

    /* Gain: '<S344>/to m//s' incorporates:
     *  Constant: '<Root>/Speed limitation (kph)'
     */
    VCM20_B.toms = VCM20_P.toms_Gain * VCM20_P.Speedlimitationkph_Value;

    /* Gain: '<S344>/toTire_rpm' */
    u1_0 = 60.0 / (6.2831853071795862 * VCM20_P.Rwf);

    /* Gain: '<S344>/toTire_rpm' */
    VCM20_B.toTire_rpm = u1_0 * VCM20_B.toms;

    /* Gain: '<S344>/Gear' */
    VCM20_B.Gear = VCM20_P.gear_f * VCM20_B.toTire_rpm;

    /* MinMax: '<S19>/Min' incorporates:
     *  Constant: '<Root>/Speed limitation (rpm)'
     */
    u0_0 = VCM20_P.Speedlimitationrpm_Value;
    u1_0 = VCM20_B.Gear;
    if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S19>/Min' */
    VCM20_B.Min = u1_0;

    /* Saturate: '<Root>/Saturation2' */
    u0_0 = VCM20_B.Min;
    u1 = VCM20_P.Saturation2_LowerSat_b;
    u2 = VCM20_P.Saturation2_UpperSat_n;
    if (u0_0 > u2) {
      /* Saturate: '<Root>/Saturation2' */
      VCM20_B.Saturation2 = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<Root>/Saturation2' */
      VCM20_B.Saturation2 = u1;
    } else {
      /* Saturate: '<Root>/Saturation2' */
      VCM20_B.Saturation2 = u0_0;
    }

    /* End of Saturate: '<Root>/Saturation2' */

    /* S-Function (rti_commonblock): '<S108>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:313 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X139]);

      /* Convert timestamp */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X139]->processed) {
        can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X139]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X139]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X139]->timestamp > 0.0) {
        if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X139]->processed) {
          CAN_Msg = can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X139]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "VELOCITY_X" (0|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o1 = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "VELOCITY_Y" (16|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o2 = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "VELOCITY_Z" (32|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o3 = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* Gain: '<S14>/Gain8' */
    VCM20_B.VELOCITY_X = VCM20_P.Gain8_Gain * VCM20_B.SFunction1_o1;

    /* Gain: '<S14>/Gain1' */
    VCM20_B.VELOCITY_X_p = VCM20_P.Gain1_Gain_f * VCM20_B.VELOCITY_X;

    /* Product: '<S290>/Product' */
    VCM20_B.Product = VCM20_B.VELOCITY_X_p * VCM20_B.VELOCITY_X_p;

    /* Sqrt: '<S290>/Sqrt' */
    VCM20_B.Sqrt = sqrt(VCM20_B.Product);

    /* Gain: '<S290>/Time msec' incorporates:
     *  Constant: '<S290>/Constant1'
     */
    u1_0 = VCM20_P.g * 0.001;

    /* Gain: '<S290>/Time msec' */
    VCM20_B.Timemsec = VCM20_P.Timemsec_Gain * u1_0;

    /* Sum: '<S290>/Add' */
    VCM20_B.Add = VCM20_B.Sqrt + VCM20_B.Timemsec;

    /* Gain: '<S290>/toTire_rpm' */
    u1_0 = 60.0 / (6.2831853071795862 * VCM20_P.Rwf);

    /* Gain: '<S290>/toTire_rpm' */
    VCM20_B.toTire_rpm_b = u1_0 * VCM20_B.Add;

    /* Gain: '<S290>/Gear' */
    VCM20_B.Gear_g = VCM20_P.gear_f * VCM20_B.toTire_rpm_b;

    /* RelationalOperator: '<S309>/LowerRelop1' */
    VCM20_B.LowerRelop1 = (VCM20_B.Gear_g > VCM20_B.Saturation2);

    /* RelationalOperator: '<S309>/UpperRelop' incorporates:
     *  Constant: '<S290>/Constant'
     */
    VCM20_B.UpperRelop = (VCM20_B.Gear_g < VCM20_P.Constant_Value_ia);

    /* Switch: '<S309>/Switch' */
    if (VCM20_B.UpperRelop) {
      /* Switch: '<S309>/Switch' incorporates:
       *  Constant: '<S290>/Constant'
       */
      VCM20_B.Switch = VCM20_P.Constant_Value_ia;
    } else {
      /* Switch: '<S309>/Switch' */
      VCM20_B.Switch = VCM20_B.Gear_g;
    }

    /* End of Switch: '<S309>/Switch' */

    /* Switch: '<S309>/Switch2' */
    if (VCM20_B.LowerRelop1) {
      /* Switch: '<S309>/Switch2' */
      VCM20_B.Switch2 = VCM20_B.Saturation2;
    } else {
      /* Switch: '<S309>/Switch2' */
      VCM20_B.Switch2 = VCM20_B.Switch;
    }

    /* End of Switch: '<S309>/Switch2' */

    /* Sum: '<S310>/Add' incorporates:
     *  Constant: '<S290>/Rear_sliprate_ref'
     *  Constant: '<S310>/Constant2'
     */
    VCM20_B.Add_m = VCM20_P.Constant2_Value - VCM20_P.Rear_sliprate_ref_Value;

    /* Product: '<S310>/Divide' */
    VCM20_B.Divide = VCM20_B.Switch2 / VCM20_B.Add_m;
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[3] == 0) {
    /* S-Function (rti_commonblock): '<S78>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "Steer SW" Id:512 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X200]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X200]->processed) {
        can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X200]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X200]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X200]->timestamp > 0.0) {
        if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X200]->processed) {
          CAN_Msg = can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X200]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "SW1" (0|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o1_o = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "SW2" (1|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 1;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o2_k = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "SW3" (2|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 2;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o3_d = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "SW5" (3|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 3;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o4 = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "SW6" (4|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 4;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o5 = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "SW7" (5|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 5;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o6 = ((real_T) CAN_Sgn.UnsignedSgn);
          }
        }
      }
    }

    /* Logic: '<S11>/LaunchSW' */
    VCM20_B.LaunchSW = !(VCM20_B.SFunction1_o5 != 0.0);
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* S-Function (rti_commonblock): '<S69>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "Brakes" Id:514 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X202]);

      /* Convert timestamp */
      if (can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X202]->processed) {
        can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X202]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X202]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X202]->timestamp > 0.0) {
        if (can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X202]->processed) {
          CAN_Msg = can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X202]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "Brake Pressure R" (0|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o1_l = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "Brake Pressure F" (16|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o2_a = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "Brake SW" (32|8, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.UnsignedSgn &= 0x000000FF;
            VCM20_B.SFunction1_o3_l = ((real_T) CAN_Sgn.UnsignedSgn);
          }
        }
      }
    }

    /* Logic: '<S11>/Brake SW' */
    VCM20_B.BrakeSW = !(VCM20_B.SFunction1_o3_l != 0.0);
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S21>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "MC_11" Id:643 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X283]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X283]->processed) {
        can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X283]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X283]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X283]->timestamp > 0.0) {
        if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X283]->processed) {
          CAN_Msg = can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X283]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "MC1_bSystemReady" (8|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o1_m = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_bError" (9|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 1;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o2_c = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_bWarn" (10|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 2;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o3_i = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_bQuitDcOn" (11|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 3;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o4_d = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_bDcOn" (12|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 4;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o5_g = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_bQuitInverterOn" (13|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 5;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o6_f = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_bInverterOn" (14|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 6;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o7 = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_bDerating" (15|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 7;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o8 = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_Actual speed value in 1/rpm" (16|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o9 = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "MC1_DC bus voltage" (32|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o10 = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_Actual torque" (48|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[7];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o11 = ((real_T) CAN_Sgn.UnsignedSgn);
          }
        }
      }
    }

    /* S-Function (rti_commonblock): '<S27>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "MC_31" Id:647 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X287]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X287]->processed) {
        can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X287]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X287]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X287]->timestamp > 0.0) {
        if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X287]->processed) {
          CAN_Msg = can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X287]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "MC3_bSystemReady" (8|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o1_b = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_bError" (9|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 1;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o2_e = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_bWarn" (10|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 2;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o3_o = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_bQuitDcOn" (11|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 3;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o4_g = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_bDcOn" (12|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 4;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o5_p = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_bQuitInverterOn" (13|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 5;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o6_h = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_bInverterOn" (14|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 6;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o7_n = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_bDerating" (15|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 7;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o8_g = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_Actual speed value in 1/rpm" (16|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o9_j = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "MC3_DC bus voltage" (32|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o10_n = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_Actual torque" (48|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[7];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o11_n = ((real_T) CAN_Sgn.UnsignedSgn);
          }
        }
      }
    }

    /* S-Function (rti_commonblock): '<S24>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "MC_21" Id:644 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X284]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X284]->processed) {
        can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X284]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X284]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X284]->timestamp > 0.0) {
        if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X284]->processed) {
          CAN_Msg = can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X284]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "MC2_bSystemReady" (8|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o1_e = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_bError" (9|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 1;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o2_b = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_bWarn" (10|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 2;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o3_k = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_bQuitDcOn" (11|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 3;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o4_j = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_bDcOn" (12|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 4;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o5_e = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_bQuitInverterOn" (13|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 5;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o6_n = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_bInverterOn" (14|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 6;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o7_j = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_bDerating" (15|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 7;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o8_o = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_Actual speed value in 1/rpm" (16|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o9_n = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "MC2_DC bus voltage" (32|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o10_ny = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_Actual torque" (48|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[7];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o11_g = ((real_T) CAN_Sgn.UnsignedSgn);
          }
        }
      }
    }

    /* S-Function (rti_commonblock): '<S30>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "MC_41" Id:648 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X288]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X288]->processed) {
        can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X288]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X288]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X288]->timestamp > 0.0) {
        if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X288]->processed) {
          CAN_Msg = can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X288]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "MC4_bSystemReady" (8|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o1_g = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_bError" (9|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 1;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o2_f = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_bWarn" (10|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 2;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o3_n = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_bQuitDcOn" (11|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 3;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o4_n = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_bDcOn" (12|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 4;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o5_j = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_bQuitInverterOn" (13|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 5;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o6_j = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_bInverterOn" (14|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 6;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o7_i = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_bDerating" (15|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 7;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o8_p = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_Actual speed value in 1/rpm" (16|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o9_g = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "MC4_DC bus voltage" (32|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o10_l = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_Actual torque" (48|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[7];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o11_p = ((real_T) CAN_Sgn.UnsignedSgn);
          }
        }
      }
    }

    /* MinMax: '<S277>/Min' */
    u0_0 = VCM20_B.SFunction1_o9;
    u1_0 = VCM20_B.SFunction1_o9_j;
    if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    u1 = VCM20_B.SFunction1_o9_n;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.SFunction1_o9_g;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    /* MinMax: '<S277>/Min' */
    VCM20_B.maxrpm = u1_0;
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* S-Function (rti_commonblock): '<S68>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "APPS_Steer" Id:513 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X201]);

      /* Convert timestamp */
      if (can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X201]->processed) {
        can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X201]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X201]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X201]->timestamp > 0.0) {
        if (can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X201]->processed) {
          VCM20_B.SFunction1_o5_l = (real_T)
            can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X201]->processed;
          CAN_Msg = can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X201]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "APPS1" (0|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o1_n = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "APPS2" (16|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o2_j = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "Steer_Stroke1" (32|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o3_d2 = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "Steer_Stroke2" (48|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[7];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o4_h = ((real_T) CAN_Sgn.UnsignedSgn);
          }
        }
      }

      if (!can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X201]->processed) {
        /* ... set RX status to 0 because no new message has arrived */
        VCM20_B.SFunction1_o5_l = 0.0;
      }
    }

    /* Sum: '<S70>/Sum1' incorporates:
     *  Constant: '<S70>/Throttle1_MIN'
     */
    VCM20_B.Sum1 = VCM20_B.SFunction1_o1_n - VCM20_P.Throttle1_MIN_Value;

    /* Sum: '<S70>/Sum3' incorporates:
     *  Constant: '<S70>/Throttle1_MAX'
     *  Constant: '<S70>/Throttle1_MIN'
     */
    VCM20_B.Sum3 = VCM20_P.Throttle1_MAX_Value - VCM20_P.Throttle1_MIN_Value;

    /* Product: '<S70>/Divide' */
    VCM20_B.Divide_j = VCM20_B.Sum1 / VCM20_B.Sum3;

    /* Sum: '<S71>/Sum1' incorporates:
     *  Constant: '<S71>/Throttle2_MIN'
     */
    VCM20_B.Sum1_p = VCM20_B.SFunction1_o2_j - VCM20_P.Throttle2_MIN_Value;

    /* Sum: '<S71>/Sum3' incorporates:
     *  Constant: '<S71>/Throttle2_MAX'
     *  Constant: '<S71>/Throttle2_MIN'
     */
    VCM20_B.Sum3_h = VCM20_P.Throttle2_MAX_Value - VCM20_P.Throttle2_MIN_Value;

    /* Product: '<S71>/Divide' */
    VCM20_B.Divide_o = VCM20_B.Sum1_p / VCM20_B.Sum3_h;

    /* MinMax: '<S67>/AppsMni' */
    u0_0 = VCM20_B.Divide_j;
    u1_0 = VCM20_B.Divide_o;
    if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S67>/AppsMni' */
    VCM20_B.AppsMni = u1_0;

    /* Sum: '<S67>/Add' */
    VCM20_B.Add_j = VCM20_B.Divide_j - VCM20_B.Divide_o;

    /* Abs: '<S67>/Abs' */
    VCM20_B.Abs = fabs(VCM20_B.Add_j);

    /* RelationalOperator: '<S83>/Compare' incorporates:
     *  Constant: '<S83>/Constant'
     */
    VCM20_B.Compare = (VCM20_B.Abs >= VCM20_P.CompareToConstant_const);

    /* Delay: '<S82>/Delay' */
    VCM20_B.Delay_b = VCM20_B.Compare;

    /* Delay: '<S82>/Delay1' */
    VCM20_B.Delay1_c = VCM20_DW.Delay1_DSTATE_c;

    /* Delay: '<S82>/Delay2' */
    VCM20_B.Delay2_p = VCM20_DW.Delay2_DSTATE_a[0];

    /* Delay: '<S82>/Delay3' */
    VCM20_B.Delay3_c1 = VCM20_DW.Delay3_DSTATE_p[0];

    /* Delay: '<S82>/Delay4' */
    VCM20_B.Delay4_or = VCM20_DW.Delay4_DSTATE_o[0];

    /* Logic: '<S82>/Logical Operator' */
    VCM20_B.LogicalOperator = (VCM20_B.Delay_b && VCM20_B.Delay1_c &&
      VCM20_B.Delay2_p && VCM20_B.Delay3_c1 && VCM20_B.Delay4_or);

    /* Logic: '<S67>/Logical Operator3' */
    VCM20_B.LogicalOperator3 = VCM20_B.LogicalOperator;

    /* DataTypeConversion: '<S81>/Data Type Conversion' */
    VCM20_B.DataTypeConversion_o = VCM20_B.LogicalOperator3;

    /* MATLABSystem: '<S81>/Moving Average' */
    u0 = VCM20_B.DataTypeConversion_o;
    obj = &VCM20_DW.obj;
    obj_0 = obj;
    if (obj_0->TunablePropsChanged) {
      obj_0->TunablePropsChanged = false;
    }

    obj_1 = obj->pStatistic;
    if (obj_1->isInitialized != 1) {
      obj_2 = obj_1;
      obj_3 = obj_2;
      obj_3->isSetupComplete = false;
      obj_3->isInitialized = 1;
      obj_4 = obj_3;
      obj_4->pCumSum = 0.0F;
      for (i = 0; i < 199; i++) {
        obj_4->pCumSumRev[i] = 0.0F;
      }

      obj_4->pCumRevIndex = 1.0F;
      obj_4->pModValueRev = 0.0F;
      obj_3->isSetupComplete = true;
      obj_2->pCumSum = 0.0F;
      for (i = 0; i < 199; i++) {
        obj_2->pCumSumRev[i] = 0.0F;
      }

      obj_2->pCumRevIndex = 1.0F;
      obj_2->pModValueRev = 0.0F;
    }

    cumRevIndex = obj_1->pCumRevIndex;
    csum = obj_1->pCumSum;
    for (i = 0; i < 199; i++) {
      VCM20_B.csumrev[i] = obj_1->pCumSumRev[i];
    }

    modValueRev = obj_1->pModValueRev;
    z = 0.0F;
    tmp = 0.0F;
    csum += u0;
    if (modValueRev == 0.0F) {
      z = VCM20_B.csumrev[(int32_T)cumRevIndex - 1] + csum;
    }

    VCM20_B.csumrev[(int32_T)cumRevIndex - 1] = u0;
    if (cumRevIndex != 199.0F) {
      cumRevIndex++;
    } else {
      cumRevIndex = 1.0F;
      csum = 0.0F;
      for (i = 197; i >= 0; i--) {
        VCM20_B.csumrev[i] += VCM20_B.csumrev[i + 1];
      }
    }

    if (modValueRev == 0.0F) {
      tmp = z / 200.0F;
    }

    if (modValueRev > 0.0F) {
      modValueRev--;
    } else {
      modValueRev = 0.0F;
    }

    obj_1->pCumSum = csum;
    for (i = 0; i < 199; i++) {
      obj_1->pCumSumRev[i] = VCM20_B.csumrev[i];
    }

    obj_1->pCumRevIndex = cumRevIndex;
    obj_1->pModValueRev = modValueRev;

    /* MATLABSystem: '<S81>/Moving Average' */
    VCM20_B.MovingAverage = tmp;

    /* RelationalOperator: '<S88>/Compare' incorporates:
     *  Constant: '<S88>/Constant'
     */
    VCM20_B.Compare_n = (VCM20_B.MovingAverage >
                         VCM20_P.CompareToConstant_const_k);

    /* Logic: '<S67>/Logical Operator1' */
    VCM20_B.LogicalOperator1 = !VCM20_B.Compare_n;

    /* Logic: '<S67>/AND' */
    VCM20_B.AND = (VCM20_B.LogicalOperator1 && (VCM20_B.SFunction1_o5_l != 0.0));

    /* Product: '<S67>/Divide' */
    VCM20_B.Divide_p = VCM20_B.AppsMni * (real_T)VCM20_B.AND;

    /* Saturate: '<S67>/Saturation' */
    u0_0 = VCM20_B.Divide_p;
    u1 = VCM20_P.Saturation_LowerSat_k;
    u2 = VCM20_P.Saturation_UpperSat_p;
    if (u0_0 > u2) {
      /* Saturate: '<S67>/Saturation' */
      VCM20_B.Saturation_d = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S67>/Saturation' */
      VCM20_B.Saturation_d = u1;
    } else {
      /* Saturate: '<S67>/Saturation' */
      VCM20_B.Saturation_d = u0_0;
    }

    /* End of Saturate: '<S67>/Saturation' */

    /* Saturate: '<S121>/Saturation1' */
    u0_0 = VCM20_B.Saturation_d;
    u1 = VCM20_P.Saturation1_LowerSat_g;
    u2 = VCM20_P.Saturation1_UpperSat_d;
    if (u0_0 > u2) {
      /* Saturate: '<S121>/Saturation1' */
      VCM20_B.APPSsig = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S121>/Saturation1' */
      VCM20_B.APPSsig = u1;
    } else {
      /* Saturate: '<S121>/Saturation1' */
      VCM20_B.APPSsig = u0_0;
    }

    /* End of Saturate: '<S121>/Saturation1' */

    /* Sum: '<S121>/Add' incorporates:
     *  Constant: '<S15>/Trottle_Pedal bias'
     */
    VCM20_B.Add_a = VCM20_B.APPSsig - VCM20_P.Trottle_Pedalbias_Value;
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<Root>/Saturation3' incorporates:
     *  Constant: '<Root>/CoastTorque'
     */
    u0_0 = VCM20_P.CoastTorque_Value;
    u1 = VCM20_P.Saturation3_LowerSat_d;
    u2 = VCM20_P.Saturation3_UpperSat_a;
    if (u0_0 > u2) {
      /* Saturate: '<Root>/Saturation3' */
      VCM20_B.Saturation3 = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<Root>/Saturation3' */
      VCM20_B.Saturation3 = u1;
    } else {
      /* Saturate: '<Root>/Saturation3' */
      VCM20_B.Saturation3 = u0_0;
    }

    /* End of Saturate: '<Root>/Saturation3' */

    /* RelationalOperator: '<S169>/Compare' incorporates:
     *  Constant: '<S169>/Constant'
     */
    VCM20_B.Compare_d = (VCM20_B.Saturation3 == VCM20_P.Constant_Value_m);

    /* Switch: '<S167>/Switch' */
    if (VCM20_B.Compare_d) {
      /* Switch: '<S167>/Switch' incorporates:
       *  Constant: '<S167>/Constant'
       */
      VCM20_B.Switch_c = VCM20_P.Constant_Value_e;
    } else {
      /* Switch: '<S167>/Switch' */
      VCM20_B.Switch_c = VCM20_B.Saturation3;
    }

    /* End of Switch: '<S167>/Switch' */
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* Switch: '<S121>/Switch' */
    if (VCM20_B.Add_a > VCM20_P.Switch_Threshold) {
      /* Fcn: '<S121>/Fcn' incorporates:
       *  Constant: '<S15>/Trottle_Pedal bias'
       */
      VCM20_B.Fcn_h = 1.0 - VCM20_P.Trottle_Pedalbias_Value;

      /* Product: '<S121>/Divide' */
      VCM20_B.Divide_jo = VCM20_B.Add_a / VCM20_B.Fcn_h;

      /* Product: '<S121>/Divide2' */
      VCM20_B.Divide2_n = VCM20_B.Saturation * VCM20_B.Divide_jo;

      /* Switch: '<S121>/Switch' */
      VCM20_B.totaltrqreq = VCM20_B.Divide2_n;
    } else {
      /* Product: '<S121>/Divide3' incorporates:
       *  Constant: '<S15>/Trottle_Pedal bias'
       */
      VCM20_B.Divide3_a5 = VCM20_B.Add_a / VCM20_P.Trottle_Pedalbias_Value;

      /* Product: '<S121>/Divide1' */
      VCM20_B.Divide1_j = VCM20_B.Divide3_a5 * VCM20_B.Switch_c;

      /* Switch: '<S121>/Switch' */
      VCM20_B.totaltrqreq = VCM20_B.Divide1_j;
    }

    /* End of Switch: '<S121>/Switch' */
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Gain: '<S121>/Gain1' */
    VCM20_B.Gain1 = VCM20_P.Gain1_Gain_d * VCM20_B.Switch_c;

    /* Saturate: '<Root>/Saturation1' incorporates:
     *  Constant: '<Root>/MinTorque'
     */
    u0_0 = VCM20_P.MinTorque_Value;
    u1 = VCM20_P.Saturation1_LowerSat_i;
    u2 = VCM20_P.Saturation1_UpperSat_m;
    if (u0_0 > u2) {
      /* Saturate: '<Root>/Saturation1' */
      VCM20_B.Saturation1 = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<Root>/Saturation1' */
      VCM20_B.Saturation1 = u1;
    } else {
      /* Saturate: '<Root>/Saturation1' */
      VCM20_B.Saturation1 = u0_0;
    }

    /* End of Saturate: '<Root>/Saturation1' */

    /* Sum: '<S121>/Add3' */
    VCM20_B.MinTorque = VCM20_B.Saturation1 - VCM20_B.Switch_c;

    /* Gain: '<S121>/Gain2' */
    VCM20_B.Gain2 = VCM20_P.Gain2_Gain * VCM20_B.Saturation1;
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* Switch: '<S121>/Switch2' */
    if (VCM20_B.totaltrqreq > VCM20_P.Switch2_Threshold) {
      /* Switch: '<S121>/Switch2' */
      VCM20_B.Switch2_e = VCM20_B.totaltrqreq;
    } else {
      /* Sum: '<S76>/Sum3' incorporates:
       *  Constant: '<S76>/BrakePressF_MAX'
       *  Constant: '<S76>/BrakePressF_MIN'
       */
      VCM20_B.Sum3_hp = VCM20_P.BrakePressF_MAX_Value -
        VCM20_P.BrakePressF_MIN_Value;

      /* Sum: '<S76>/Sum1' incorporates:
       *  Constant: '<S76>/BrakePressF_MIN'
       */
      VCM20_B.Sum1_fo = VCM20_B.SFunction1_o2_a - VCM20_P.BrakePressF_MIN_Value;

      /* Product: '<S76>/Divide' */
      VCM20_B.Divide_ef = VCM20_B.Sum1_fo / VCM20_B.Sum3_hp;

      /* Saturate: '<S121>/Saturation2' */
      u0_0 = VCM20_B.Divide_ef;
      u1 = VCM20_P.Saturation2_LowerSat;
      u2 = VCM20_P.Saturation2_UpperSat;
      if (u0_0 > u2) {
        /* Saturate: '<S121>/Saturation2' */
        VCM20_B.Saturation2_m = u2;
      } else if (u0_0 < u1) {
        /* Saturate: '<S121>/Saturation2' */
        VCM20_B.Saturation2_m = u1;
      } else {
        /* Saturate: '<S121>/Saturation2' */
        VCM20_B.Saturation2_m = u0_0;
      }

      /* End of Saturate: '<S121>/Saturation2' */

      /* Sum: '<S121>/Add4' incorporates:
       *  Constant: '<S15>/BrakeEnd_pedal bias'
       */
      VCM20_B.Add4_c = VCM20_P.BrakeEnd_pedalbias_Value - VCM20_B.Saturation2_m;

      /* Switch: '<S121>/Switch3' */
      if (VCM20_B.Add4_c > VCM20_P.Switch3_Threshold) {
        /* Sum: '<S121>/Add1' incorporates:
         *  Constant: '<S15>/BrakeStart_Pedal bias'
         */
        VCM20_B.Add1_jj = VCM20_B.Saturation2_m -
          VCM20_P.BrakeStart_Pedalbias_Value;

        /* Switch: '<S121>/Switch1' */
        if (VCM20_B.Add1_jj > VCM20_P.Switch1_Threshold) {
          /* Sum: '<S121>/Add5' incorporates:
           *  Constant: '<S15>/BrakeEnd_pedal bias'
           *  Constant: '<S15>/BrakeStart_Pedal bias'
           */
          VCM20_B.Add5 = VCM20_P.BrakeEnd_pedalbias_Value -
            VCM20_P.BrakeStart_Pedalbias_Value;

          /* Product: '<S121>/Divide4' */
          VCM20_B.Divide4_i = VCM20_B.Add1_jj / VCM20_B.Add5;

          /* Product: '<S121>/Divide6' */
          VCM20_B.Divide6 = VCM20_B.MinTorque * VCM20_B.Divide4_i;

          /* Switch: '<S121>/Switch1' */
          VCM20_B.totaltrqreq_g = VCM20_B.Divide6;
        } else {
          /* Switch: '<S121>/Switch1' incorporates:
           *  Constant: '<S121>/Constant'
           */
          VCM20_B.totaltrqreq_g = VCM20_P.Constant_Value_p;
        }

        /* End of Switch: '<S121>/Switch1' */

        /* Switch: '<S121>/Switch3' */
        VCM20_B.Switch3_ad2 = VCM20_B.totaltrqreq_g;
      } else {
        /* Switch: '<S121>/Switch3' */
        VCM20_B.Switch3_ad2 = VCM20_B.MinTorque;
      }

      /* End of Switch: '<S121>/Switch3' */

      /* Gain: '<S121>/Gain' */
      VCM20_B.Gain_a2 = VCM20_P.Gain_Gain * VCM20_B.Switch3_ad2;

      /* Product: '<S121>/Divide8' */
      VCM20_B.Divide8 = VCM20_B.totaltrqreq / VCM20_B.Gain1;

      /* Saturate: '<S121>/Saturation' */
      u0_0 = VCM20_B.Divide8;
      u1 = VCM20_P.Saturation_LowerSat;
      u2 = VCM20_P.Saturation_UpperSat;
      if (u0_0 > u2) {
        /* Saturate: '<S121>/Saturation' */
        VCM20_B.Saturation_fq = u2;
      } else if (u0_0 < u1) {
        /* Saturate: '<S121>/Saturation' */
        VCM20_B.Saturation_fq = u1;
      } else {
        /* Saturate: '<S121>/Saturation' */
        VCM20_B.Saturation_fq = u0_0;
      }

      /* End of Saturate: '<S121>/Saturation' */

      /* Product: '<S121>/Divide9' */
      VCM20_B.Divide9 = VCM20_B.Saturation_fq * VCM20_B.Gain_a2;

      /* Sum: '<S121>/Add2' */
      VCM20_B.Add2_j = VCM20_B.totaltrqreq + VCM20_B.Divide9;

      /* Switch: '<S121>/Switch2' */
      VCM20_B.Switch2_e = VCM20_B.Add2_j;
    }

    /* End of Switch: '<S121>/Switch2' */

    /* RelationalOperator: '<S168>/LowerRelop1' */
    VCM20_B.LowerRelop1_n = (VCM20_B.Switch2_e > VCM20_B.Saturation);

    /* RelationalOperator: '<S168>/UpperRelop' */
    VCM20_B.UpperRelop_h = (VCM20_B.Switch2_e < VCM20_B.Gain2);

    /* Switch: '<S168>/Switch' */
    if (VCM20_B.UpperRelop_h) {
      /* Switch: '<S168>/Switch' */
      VCM20_B.Switch_i = VCM20_B.Gain2;
    } else {
      /* Switch: '<S168>/Switch' */
      VCM20_B.Switch_i = VCM20_B.Switch2_e;
    }

    /* End of Switch: '<S168>/Switch' */

    /* Switch: '<S168>/Switch2' */
    if (VCM20_B.LowerRelop1_n) {
      /* Switch: '<S168>/Switch2' */
      VCM20_B.Switch2_l = VCM20_B.Saturation;
    } else {
      /* Switch: '<S168>/Switch2' */
      VCM20_B.Switch2_l = VCM20_B.Switch_i;
    }

    /* End of Switch: '<S168>/Switch2' */

    /* Sum: '<S127>/Sum' incorporates:
     *  Constant: '<S15>/Trottle_Pedal bias'
     */
    VCM20_B.Sum = VCM20_B.APPSsig - VCM20_P.Trottle_Pedalbias_Value;

    /* RelationalOperator: '<S128>/Compare' incorporates:
     *  Constant: '<S128>/Constant'
     */
    VCM20_B.Compare_k = (VCM20_B.Sum >= VCM20_P.Constant_Value_o);

    /* Logic: '<S127>/Logical Operator' */
    VCM20_B.LogicalOperator_e = (VCM20_B.Compare_k && VCM20_B.BrakeSW);

    /* Delay: '<S127>/Delay' */
    VCM20_B.Delay_p = VCM20_DW.Delay_DSTATE_d;

    /* Logic: '<S127>/Logical Operator2' */
    VCM20_B.LogicalOperator2 = (VCM20_B.LogicalOperator_e || VCM20_B.Delay_p);

    /* Logic: '<S127>/Logical Operator1' */
    VCM20_B.LogicalOperator1_o = !VCM20_B.LogicalOperator2;

    /* Switch: '<S115>/Brake Over Ride Sw' */
    if (VCM20_B.LogicalOperator1_o) {
      /* Switch: '<S115>/Brake Over Ride Sw' */
      VCM20_B.BrakeOverRideSw = VCM20_B.Switch2_l;
    } else {
      /* Switch: '<S115>/Brake Over Ride Sw' incorporates:
       *  Constant: '<S115>/Constant1'
       */
      VCM20_B.BrakeOverRideSw = VCM20_P.Constant1_Value;
    }

    /* End of Switch: '<S115>/Brake Over Ride Sw' */
  }

  /* TransferFcn: '<S174>/Transfer Fcn' */
  VCM20_B.TransferFcn = 0.0;
  VCM20_B.TransferFcn += VCM20_P.TransferFcn_C * VCM20_X.TransferFcn_CSTATE;

  /* Gain: '<S281>/hg//L' */
  VCM20_B.hgL = VCM20_P.hgL_Gain * VCM20_B.TransferFcn;

  /* Product: '<S281>/Divide' incorporates:
   *  Constant: '<S281>/M'
   */
  VCM20_B.Divide_ob = VCM20_B.hgL * VCM20_P.M_Value;

  /* Gain: '<S281>/Gain' */
  VCM20_B.Gain = VCM20_P.Gain_Gain_n * VCM20_B.Divide_ob;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Switch: '<S126>/Switch1' incorporates:
     *  Constant: '<Root>/MC1_sw'
     */
    if (VCM20_P.MC1_sw_Value > VCM20_P.Switch1_Threshold_p) {
      /* Switch: '<S126>/Switch1' */
      VCM20_B.Switch1 = VCM20_B.SFunction1_o6_f;
    } else {
      /* Switch: '<S126>/Switch1' incorporates:
       *  Constant: '<S126>/Constant'
       */
      VCM20_B.Switch1 = VCM20_P.Constant_Value_pl;
    }

    /* End of Switch: '<S126>/Switch1' */

    /* Switch: '<S126>/Switch2' incorporates:
     *  Constant: '<Root>/MC2_sw'
     */
    if (VCM20_P.MC2_sw_Value > VCM20_P.Switch2_Threshold_n) {
      /* Switch: '<S126>/Switch2' */
      VCM20_B.Switch2_a = VCM20_B.SFunction1_o6_n;
    } else {
      /* Switch: '<S126>/Switch2' incorporates:
       *  Constant: '<S126>/Constant'
       */
      VCM20_B.Switch2_a = VCM20_P.Constant_Value_pl;
    }

    /* End of Switch: '<S126>/Switch2' */

    /* Switch: '<S126>/Switch3' incorporates:
     *  Constant: '<Root>/MC3_sw'
     */
    if (VCM20_P.MC3_sw_Value > VCM20_P.Switch3_Threshold_d) {
      /* Switch: '<S126>/Switch3' */
      VCM20_B.Switch3 = VCM20_B.SFunction1_o6_h;
    } else {
      /* Switch: '<S126>/Switch3' incorporates:
       *  Constant: '<S126>/Constant'
       */
      VCM20_B.Switch3 = VCM20_P.Constant_Value_pl;
    }

    /* End of Switch: '<S126>/Switch3' */

    /* Switch: '<S126>/Switch4' incorporates:
     *  Constant: '<Root>/MC4_sw'
     */
    if (VCM20_P.MC4_sw_Value > VCM20_P.Switch4_Threshold) {
      /* Switch: '<S126>/Switch4' */
      VCM20_B.Switch4 = VCM20_B.SFunction1_o6_j;
    } else {
      /* Switch: '<S126>/Switch4' incorporates:
       *  Constant: '<S126>/Constant'
       */
      VCM20_B.Switch4 = VCM20_P.Constant_Value_pl;
    }

    /* End of Switch: '<S126>/Switch4' */

    /* Logic: '<S126>/All Inverter Enable' */
    VCM20_B.AllInverterEnable = ((VCM20_B.Switch1 != 0.0) && (VCM20_B.Switch2_a
      != 0.0) && (VCM20_B.Switch3 != 0.0) && (VCM20_B.Switch4 != 0.0));

    /* Switch: '<S113>/Switch' */
    if (VCM20_B.AllInverterEnable) {
      /* Switch: '<S113>/Switch' */
      VCM20_B.Switch_j = VCM20_B.BrakeOverRideSw;
    } else {
      /* Switch: '<S113>/Switch' incorporates:
       *  Constant: '<S113>/Constant'
       */
      VCM20_B.Switch_j = VCM20_P.Constant_Value;
    }

    /* End of Switch: '<S113>/Switch' */

    /* Chart: '<S172>/Chart' incorporates:
     *  Constant: '<S172>/GearProtectionTime (s)'
     *  Constant: '<S172>/GearProtectionTorqueLimit'
     */
    if (VCM20_DW.temporalCounter_i1 < MAX_uint32_T) {
      VCM20_DW.temporalCounter_i1++;
    }

    /* Gateway: VCM main/VDC/Launch/Chart */
    /* During: VCM main/VDC/Launch/Chart */
    if (VCM20_DW.is_active_c1_VCM20 == 0U) {
      /* Entry: VCM main/VDC/Launch/Chart */
      VCM20_DW.is_active_c1_VCM20 = 1U;

      /* Entry Internal: VCM main/VDC/Launch/Chart */
      /* Transition: '<S272>:2' */
      VCM20_DW.is_c1_VCM20 = VCM20_IN_Base;

      /* Entry 'Base': '<S272>:1' */
      VCM20_B.LaunchTorqueLimit = VCM20_B.Saturation;
      VCM20_B.LaunchLED = 0.0;
    } else {
      switch (VCM20_DW.is_c1_VCM20) {
       case VCM20_IN_Base:
        VCM20_B.LaunchLED = 0.0;

        /* During 'Base': '<S272>:1' */
        if (VCM20_B.LaunchSW && VCM20_B.BrakeSW && (VCM20_B.maxrpm < 10.0) &&
            (VCM20_B.Switch_j < 0.3)) {
          /* Transition: '<S272>:4' */
          VCM20_DW.is_c1_VCM20 = VCM20_IN_LaunchReady;
          VCM20_DW.temporalCounter_i1 = 0U;

          /* Entry 'LaunchReady': '<S272>:3' */
          VCM20_B.LaunchTorqueLimit = VCM20_B.Saturation;
        }
        break;

       case VCM20_IN_LaunchGearProtection:
        VCM20_B.LaunchLED = 1.0;

        /* During 'LaunchGearProtection': '<S272>:5' */
        if (VCM20_DW.temporalCounter_i1 >= (uint32_T)ceil
            (VCM20_P.GearProtectionTimes_Value * 1000.0)) {
          /* Transition: '<S272>:35' */
          VCM20_DW.is_c1_VCM20 = VCM20_IN_Base;

          /* Entry 'Base': '<S272>:1' */
          VCM20_B.LaunchTorqueLimit = VCM20_B.Saturation;
          VCM20_B.LaunchLED = 0.0;
        }
        break;

       case VCM20_IN_LaunchReady:
        /* During 'LaunchReady': '<S272>:3' */
        out = ((VCM20_DW.temporalCounter_i1 >= 2000U) && VCM20_B.BrakeSW);
        if (out) {
          /* Transition: '<S272>:8' */
          VCM20_DW.is_c1_VCM20 = VCM20_IN_LaunchReady_BrakeON;
          VCM20_DW.temporalCounter_i1 = 0U;

          /* Entry 'LaunchReady_BrakeON': '<S272>:11' */
          VCM20_B.LaunchTorqueLimit = 0.0;
          VCM20_B.LaunchLED = 1.0;
        } else if (!VCM20_B.LaunchSW) {
          /* Transition: '<S272>:23' */
          /* Transition: '<S272>:26' */
          VCM20_DW.is_c1_VCM20 = VCM20_IN_OFF;

          /* Entry 'OFF': '<S272>:31' */
          VCM20_B.LaunchTorqueLimit = 0.0;
        }
        break;

       case VCM20_IN_LaunchReady_BrakeON:
        VCM20_B.LaunchLED = 1.0;

        /* During 'LaunchReady_BrakeON': '<S272>:11' */
        if (!VCM20_B.LaunchSW) {
          /* Transition: '<S272>:15' */
          /* Transition: '<S272>:25' */
          /* Transition: '<S272>:26' */
          VCM20_DW.is_c1_VCM20 = VCM20_IN_OFF;

          /* Entry 'OFF': '<S272>:31' */
          VCM20_B.LaunchTorqueLimit = 0.0;
        } else {
          out = ((VCM20_DW.temporalCounter_i1 >= 1000U) && ((!VCM20_B.BrakeSW) &&
                  (VCM20_B.Switch_j > 5.0)));
          if (out) {
            /* Transition: '<S272>:14' */
            VCM20_DW.is_c1_VCM20 = VCM20_IN_LaunchReady_BrakeOff;

            /* Entry 'LaunchReady_BrakeOff': '<S272>:13' */
            VCM20_B.LaunchTorqueLimit = 0.0;

            /* Entry Internal 'LaunchReady_BrakeOff': '<S272>:13' */
            /* Transition: '<S272>:43' */
            VCM20_DW.is_LaunchReady_BrakeOff = VCM20_IN_LaunchLEDon;
            VCM20_DW.temporalCounter_i1 = 0U;

            /* Entry 'LaunchLEDon': '<S272>:45' */
            VCM20_B.LaunchLED = 1.0;
          }
        }
        break;

       case VCM20_IN_LaunchReady_BrakeOff:
        /* During 'LaunchReady_BrakeOff': '<S272>:13' */
        if ((VCM20_B.Switch_j < 4.0) || VCM20_B.BrakeSW) {
          /* Transition: '<S272>:20' */
          /* Exit Internal 'LaunchReady_BrakeOff': '<S272>:13' */
          VCM20_DW.is_LaunchReady_BrakeOff = VCM20_IN_NO_ACTIVE_CHILD;
          VCM20_DW.is_c1_VCM20 = VCM20_IN_LaunchReady_BrakeON;
          VCM20_DW.temporalCounter_i1 = 0U;

          /* Entry 'LaunchReady_BrakeON': '<S272>:11' */
          VCM20_B.LaunchTorqueLimit = 0.0;
          VCM20_B.LaunchLED = 1.0;
        } else if (!VCM20_B.LaunchSW) {
          /* Transition: '<S272>:27' */
          /* Exit Internal 'LaunchReady_BrakeOff': '<S272>:13' */
          VCM20_DW.is_LaunchReady_BrakeOff = VCM20_IN_NO_ACTIVE_CHILD;
          VCM20_DW.is_c1_VCM20 = VCM20_IN_LaunchGearProtection;
          VCM20_DW.temporalCounter_i1 = 0U;

          /* Entry 'LaunchGearProtection': '<S272>:5' */
          VCM20_B.LaunchTorqueLimit = VCM20_P.GearProtectionTorqueLimit_Value;
          VCM20_B.LaunchLED = 1.0;
        } else if (VCM20_DW.is_LaunchReady_BrakeOff == 1) {
          VCM20_B.LaunchLED = 0.0;

          /* During 'LaunchLEDoff': '<S272>:42' */
          if (VCM20_DW.temporalCounter_i1 >= 200U) {
            /* Transition: '<S272>:44' */
            VCM20_DW.is_LaunchReady_BrakeOff = VCM20_IN_LaunchLEDon;
            VCM20_DW.temporalCounter_i1 = 0U;

            /* Entry 'LaunchLEDon': '<S272>:45' */
            VCM20_B.LaunchLED = 1.0;
          }
        } else {
          VCM20_B.LaunchLED = 1.0;

          /* During 'LaunchLEDon': '<S272>:45' */
          if (VCM20_DW.temporalCounter_i1 >= 200U) {
            /* Transition: '<S272>:46' */
            VCM20_DW.is_LaunchReady_BrakeOff = VCM20_IN_LaunchLEDoff;
            VCM20_DW.temporalCounter_i1 = 0U;

            /* Entry 'LaunchLEDoff': '<S272>:42' */
            VCM20_B.LaunchLED = 0.0;
          }
        }
        break;

       default:
        /* During 'OFF': '<S272>:31' */
        if (VCM20_B.Switch_j < 0.5) {
          /* Transition: '<S272>:32' */
          VCM20_DW.is_c1_VCM20 = VCM20_IN_Base;

          /* Entry 'Base': '<S272>:1' */
          VCM20_B.LaunchTorqueLimit = VCM20_B.Saturation;
          VCM20_B.LaunchLED = 0.0;
        }
        break;
      }
    }

    /* End of Chart: '<S172>/Chart' */

    /* Switch: '<S172>/Switch' incorporates:
     *  Constant: '<Root>/LaunchOn'
     */
    if (VCM20_P.LaunchOn_Value > VCM20_P.Switch_Threshold_b) {
      /* Switch: '<S172>/Switch' */
      VCM20_B.Switch_ci = VCM20_B.LaunchTorqueLimit;
    } else {
      /* Switch: '<S172>/Switch' */
      VCM20_B.Switch_ci = VCM20_B.Saturation;
    }

    /* End of Switch: '<S172>/Switch' */

    /* Gain: '<S281>/Gain12' incorporates:
     *  Constant: '<S281>/M'
     */
    VCM20_B.m_f = VCM20_P.Gain12_Gain * VCM20_P.M_Value;
  }

  /* Sum: '<S281>/Add' */
  VCM20_B.Add_jj = VCM20_B.m_f - VCM20_B.Gain;

  /* Saturate: '<S281>/Saturation1' */
  u0_0 = VCM20_B.Add_jj;
  u1 = VCM20_P.Saturation1_LowerSat_gq;
  u2 = VCM20_P.Saturation1_UpperSat_i;
  if (u0_0 > u2) {
    /* Saturate: '<S281>/Saturation1' */
    VCM20_B.Saturation1_o = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S281>/Saturation1' */
    VCM20_B.Saturation1_o = u1;
  } else {
    /* Saturate: '<S281>/Saturation1' */
    VCM20_B.Saturation1_o = u0_0;
  }

  /* End of Saturate: '<S281>/Saturation1' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Gain: '<S281>/Gain13' incorporates:
     *  Constant: '<S281>/M'
     */
    VCM20_B.m_r = VCM20_P.Gain13_Gain * VCM20_P.M_Value;
  }

  /* Sum: '<S281>/Add1' */
  VCM20_B.Add1 = VCM20_B.Gain + VCM20_B.m_r;

  /* Saturate: '<S281>/Saturation2' */
  u0_0 = VCM20_B.Add1;
  u1 = VCM20_P.Saturation2_LowerSat_a;
  u2 = VCM20_P.Saturation2_UpperSat_f;
  if (u0_0 > u2) {
    /* Saturate: '<S281>/Saturation2' */
    VCM20_B.Saturation2_d = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S281>/Saturation2' */
    VCM20_B.Saturation2_d = u1;
  } else {
    /* Saturate: '<S281>/Saturation2' */
    VCM20_B.Saturation2_d = u0_0;
  }

  /* End of Saturate: '<S281>/Saturation2' */

  /* Product: '<S281>/Divide1' */
  VCM20_B.Divide1 = VCM20_B.Saturation1_o / VCM20_B.Saturation2_d;

  /* Sum: '<S281>/Add2' incorporates:
   *  Constant: '<S281>/Constant9'
   */
  VCM20_B.Add2 = VCM20_B.Divide1 - VCM20_P.Constant9_Value;

  /* Sum: '<S281>/Add3' incorporates:
   *  Constant: '<S281>/Constant10'
   */
  VCM20_B.Add3 = VCM20_B.Divide1 + VCM20_P.Constant10_Value;

  /* Product: '<S281>/Divide2' */
  VCM20_B.Divide2 = VCM20_B.Add2 / VCM20_B.Add3;

  /* Gain: '<S281>/Gain1' */
  VCM20_B.Gain1_c = VCM20_P.Gain1_Gain_k * VCM20_B.Divide2;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* RelationalOperator: '<S282>/LowerRelop1' */
    VCM20_B.LowerRelop1_p = (VCM20_B.Switch_j > VCM20_B.Saturation);

    /* RelationalOperator: '<S282>/UpperRelop' */
    VCM20_B.UpperRelop_f = (VCM20_B.Switch_j < VCM20_B.Saturation1);

    /* Switch: '<S282>/Switch' */
    if (VCM20_B.UpperRelop_f) {
      /* Switch: '<S282>/Switch' */
      VCM20_B.Switch_o = VCM20_B.Saturation1;
    } else {
      /* Switch: '<S282>/Switch' */
      VCM20_B.Switch_o = VCM20_B.Switch_j;
    }

    /* End of Switch: '<S282>/Switch' */

    /* Switch: '<S282>/Switch2' */
    if (VCM20_B.LowerRelop1_p) {
      /* Switch: '<S282>/Switch2' */
      VCM20_B.Switch2_f = VCM20_B.Saturation;
    } else {
      /* Switch: '<S282>/Switch2' */
      VCM20_B.Switch2_f = VCM20_B.Switch_o;
    }

    /* End of Switch: '<S282>/Switch2' */
  }

  /* Product: '<S281>/Divide3' */
  VCM20_B.Divide3 = VCM20_B.Gain1_c * VCM20_B.Switch2_f;

  /* Product: '<S281>/Divide4' incorporates:
   *  Constant: '<S174>/Trq F//R Dist Gain'
   */
  VCM20_B.Divide4 = VCM20_P.TrqFRDistGain_Value * VCM20_B.Divide3;

  /* Sum: '<S281>/Add5' */
  VCM20_B.torque_front = VCM20_B.Switch2_f - VCM20_B.Divide4;

  /* Sum: '<S278>/Add9' */
  VCM20_B.Add9 = VCM20_B.torque_front - VCM20_B.Saturation;

  /* Saturate: '<S278>/Saturation3' */
  u0_0 = VCM20_B.Add9;
  u1 = VCM20_P.Saturation3_LowerSat_o;
  u2 = VCM20_P.Saturation3_UpperSat_m;
  if (u0_0 > u2) {
    /* Saturate: '<S278>/Saturation3' */
    VCM20_B.Saturation3_g = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S278>/Saturation3' */
    VCM20_B.Saturation3_g = u1;
  } else {
    /* Saturate: '<S278>/Saturation3' */
    VCM20_B.Saturation3_g = u0_0;
  }

  /* End of Saturate: '<S278>/Saturation3' */

  /* Sum: '<S281>/Add6' */
  VCM20_B.torque_rear = VCM20_B.Divide4 + VCM20_B.Switch2_f;

  /* Sum: '<S278>/Add4' */
  VCM20_B.Add4 = VCM20_B.torque_rear - VCM20_B.Saturation;

  /* Saturate: '<S278>/Saturation' */
  u0_0 = VCM20_B.Add4;
  u1 = VCM20_P.Saturation_LowerSat_h;
  u2 = VCM20_P.Saturation_UpperSat_c;
  if (u0_0 > u2) {
    /* Saturate: '<S278>/Saturation' */
    VCM20_B.Saturation_b = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S278>/Saturation' */
    VCM20_B.Saturation_b = u1;
  } else {
    /* Saturate: '<S278>/Saturation' */
    VCM20_B.Saturation_b = u0_0;
  }

  /* End of Saturate: '<S278>/Saturation' */

  /* Sum: '<S278>/Add8' */
  VCM20_B.Add8 = VCM20_B.torque_front + VCM20_B.Saturation_b;

  /* Sum: '<S278>/Add10' */
  VCM20_B.Add10 = VCM20_B.Add8 - VCM20_B.Saturation3_g;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Gain: '<S279>/Gain' */
    VCM20_B.Gain_d = VCM20_P.Gain_Gain_l * VCM20_B.Saturation1;
  }

  /* Sum: '<S279>/Add9' */
  VCM20_B.Add9_o = VCM20_B.Add10 - VCM20_B.Gain_d;

  /* Saturate: '<S279>/Saturation3' */
  u0_0 = VCM20_B.Add9_o;
  u1 = VCM20_P.Saturation3_LowerSat_oy;
  u2 = VCM20_P.Saturation3_UpperSat_i;
  if (u0_0 > u2) {
    /* Saturate: '<S279>/Saturation3' */
    VCM20_B.Saturation3_e = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S279>/Saturation3' */
    VCM20_B.Saturation3_e = u1;
  } else {
    /* Saturate: '<S279>/Saturation3' */
    VCM20_B.Saturation3_e = u0_0;
  }

  /* End of Saturate: '<S279>/Saturation3' */

  /* Sum: '<S278>/Add7' */
  VCM20_B.Add7 = VCM20_B.torque_rear - VCM20_B.Saturation_b;

  /* Sum: '<S278>/Add11' */
  VCM20_B.Add11 = VCM20_B.Saturation3_g + VCM20_B.Add7;

  /* Sum: '<S279>/Add4' */
  VCM20_B.Add4_e = VCM20_B.Add11 - VCM20_B.Gain_d;

  /* Saturate: '<S279>/Saturation' */
  u0_0 = VCM20_B.Add4_e;
  u1 = VCM20_P.Saturation_LowerSat_p;
  u2 = VCM20_P.Saturation_UpperSat_pl;
  if (u0_0 > u2) {
    /* Saturate: '<S279>/Saturation' */
    VCM20_B.Saturation_h = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S279>/Saturation' */
    VCM20_B.Saturation_h = u1;
  } else {
    /* Saturate: '<S279>/Saturation' */
    VCM20_B.Saturation_h = u0_0;
  }

  /* End of Saturate: '<S279>/Saturation' */

  /* Sum: '<S279>/Add7' */
  VCM20_B.Add7_i = VCM20_B.Add11 - VCM20_B.Saturation_h;

  /* Sum: '<S279>/Add11' */
  VCM20_B.torque_rear_d = VCM20_B.Saturation3_e + VCM20_B.Add7_i;

  /* Switch: '<S174>/Switch1' incorporates:
   *  Constant: '<Root>/F//RTrqDistOn'
   */
  if (VCM20_P.FRTrqDistOn_Value > VCM20_P.Switch1_Threshold_b) {
    /* Switch: '<S174>/Switch1' */
    VCM20_B.Switch1_g = VCM20_B.torque_rear_d;
  } else {
    /* Switch: '<S174>/Switch1' */
    VCM20_B.Switch1_g = VCM20_B.Switch2_f;
  }

  /* End of Switch: '<S174>/Switch1' */

  /* Gain: '<S174>/Gain4' */
  VCM20_B.RRtorque = VCM20_P.Gain4_Gain * VCM20_B.Switch1_g;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Sum: '<S73>/Add1' incorporates:
     *  Constant: '<S73>/SWAS2_signal_0deg'
     *  Constant: '<S73>/SWAS2_signal_plus'
     */
    VCM20_B.Add1_j = VCM20_P.SWAS2_signal_plus_Value -
      VCM20_P.SWAS2_signal_0deg_Value;
  }

  /* TransferFcn: '<S73>/Transfer Fcn' */
  VCM20_B.TransferFcn_f = 0.0;
  VCM20_B.TransferFcn_f += VCM20_P.TransferFcn_C_p *
    VCM20_X.TransferFcn_CSTATE_k;

  /* Sum: '<S73>/Add' incorporates:
   *  Constant: '<S73>/SWAS2_signal_0deg'
   */
  VCM20_B.Add_i = VCM20_B.TransferFcn_f - VCM20_P.SWAS2_signal_0deg_Value;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Sum: '<S73>/Add2' incorporates:
     *  Constant: '<S73>/SWAS2_signal_0deg'
     *  Constant: '<S73>/SWAS2_signal_minus'
     */
    VCM20_B.Add2_h = VCM20_P.SWAS2_signal_minus_Value -
      VCM20_P.SWAS2_signal_0deg_Value;
  }

  /* Switch: '<S73>/Switch' */
  if (VCM20_B.Add_i > VCM20_P.Switch_Threshold_n) {
    /* Product: '<S73>/Divide' */
    VCM20_B.Divide_kt = 1.0 / VCM20_B.Add1_j * VCM20_B.Add_i;

    /* Gain: '<S73>/Gain1' */
    VCM20_B.Gain1_n = VCM20_P.Gain1_Gain_m * VCM20_B.Divide_kt;

    /* Gain: '<S73>/RightMax_deg' */
    VCM20_B.RightMax_deg = VCM20_P.RightMax_deg_Gain_c * VCM20_B.Gain1_n;

    /* Switch: '<S73>/Switch' */
    VCM20_B.Switch_m = VCM20_B.RightMax_deg;
  } else {
    /* Product: '<S73>/Divide1' */
    VCM20_B.Divide1_g = VCM20_B.Add_i / VCM20_B.Add2_h;

    /* Gain: '<S73>/LeftMax_deg' */
    VCM20_B.LeftMax_deg = VCM20_P.LeftMax_deg_Gain_h * VCM20_B.Divide1_g;

    /* Switch: '<S73>/Switch' */
    VCM20_B.Switch_m = VCM20_B.LeftMax_deg;
  }

  /* End of Switch: '<S73>/Switch' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* RelationalOperator: '<S93>/Compare' incorporates:
     *  Constant: '<S93>/Constant'
     */
    VCM20_B.Compare_c = (VCM20_B.SFunction1_o4_h >=
                         VCM20_P.CompareToConstant1_const);

    /* RelationalOperator: '<S94>/Compare' incorporates:
     *  Constant: '<S94>/Constant'
     */
    VCM20_B.Compare_o = (VCM20_B.SFunction1_o4_h <=
                         VCM20_P.CompareToConstant3_const);

    /* Logic: '<S77>/Logical Operator' */
    VCM20_B.LogicalOperator_g = (VCM20_B.Compare_c && VCM20_B.Compare_o);

    /* Delay: '<S92>/Delay' */
    VCM20_B.Delay_ps = VCM20_B.LogicalOperator_g;

    /* Delay: '<S92>/Delay1' */
    VCM20_B.Delay1_b = VCM20_DW.Delay1_DSTATE_ci;

    /* Delay: '<S92>/Delay2' */
    VCM20_B.Delay2_o = VCM20_DW.Delay2_DSTATE_j[0];

    /* Delay: '<S92>/Delay3' */
    VCM20_B.Delay3_b = VCM20_DW.Delay3_DSTATE_g[0];

    /* Delay: '<S92>/Delay4' */
    VCM20_B.Delay4_m = VCM20_DW.Delay4_DSTATE_g[0];

    /* Delay: '<S92>/Delay5' */
    VCM20_B.Delay5_nu = VCM20_DW.Delay5_DSTATE_g[0];

    /* Delay: '<S92>/Delay6' */
    VCM20_B.Delay6_j = VCM20_DW.Delay6_DSTATE_lu[0];

    /* Delay: '<S92>/Delay7' */
    VCM20_B.Delay7_c = VCM20_DW.Delay7_DSTATE_nq[0];

    /* Delay: '<S92>/Delay8' */
    VCM20_B.Delay8_f = VCM20_DW.Delay8_DSTATE_k[0];

    /* Delay: '<S92>/Delay9' */
    VCM20_B.Delay9_k = VCM20_DW.Delay9_DSTATE_o[0];

    /* Logic: '<S92>/Logical Operator' */
    VCM20_B.LogicalOperator_i = (VCM20_B.Delay_ps && VCM20_B.Delay1_b &&
      VCM20_B.Delay2_o && VCM20_B.Delay3_b && VCM20_B.Delay4_m &&
      VCM20_B.Delay5_nu && VCM20_B.Delay6_j && VCM20_B.Delay7_c &&
      VCM20_B.Delay8_f && VCM20_B.Delay9_k);

    /* Logic: '<S77>/AND' */
    VCM20_B.AND_i = (VCM20_B.LogicalOperator_i && (VCM20_B.SFunction1_o5_l !=
      0.0));
  }

  /* Product: '<S77>/Divide' */
  VCM20_B.Divide_i = VCM20_B.Switch_m * (real_T)VCM20_B.AND_i;

  /* Saturate: '<S77>/Saturation' */
  u0_0 = VCM20_B.Divide_i;
  u1 = VCM20_P.Saturation_LowerSat_c;
  u2 = VCM20_P.Saturation_UpperSat_f;
  if (u0_0 > u2) {
    /* Saturate: '<S77>/Saturation' */
    VCM20_B.Saturation_m = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S77>/Saturation' */
    VCM20_B.Saturation_m = u1;
  } else {
    /* Saturate: '<S77>/Saturation' */
    VCM20_B.Saturation_m = u0_0;
  }

  /* End of Saturate: '<S77>/Saturation' */

  /* Gain: '<S77>/SteerAngle_deg' */
  VCM20_B.SteerAngle_deg = VCM20_P.SteerAngle_deg_Gain * VCM20_B.Saturation_m;

  /* Gain: '<S77>/SteerAngle_rad' */
  VCM20_B.SteerAngle_rad = VCM20_P.SteerAngle_rad_Gain * VCM20_B.SteerAngle_deg;

  /* Gain: '<S216>/Gain' */
  VCM20_B.TargetYawMoment = VCM20_P.Gain_Gain_b * VCM20_B.SteerAngle_rad;

  /* Gain: '<S211>/Gain12' */
  VCM20_B.Gain12 = VCM20_P.Gain12_Gain_d * VCM20_B.TargetYawMoment;

  /* Gain: '<S211>/Gain13' */
  VCM20_B.Gain13 = VCM20_P.Gain13_Gain_f * VCM20_B.Gain12;

  /* Gain: '<S211>/Gain14' */
  VCM20_B.Gain14 = VCM20_P.Gain14_Gain * VCM20_B.Gain13;

  /* Gain: '<S211>/Gain4' */
  VCM20_B.Gain4 = VCM20_P.Gain4_Gain_n * VCM20_B.Gain14;

  /* Saturate: '<S211>/Saturation' */
  u0_0 = VCM20_B.Gain4;
  u1 = VCM20_P.Saturation_LowerSat_n4;
  u2 = VCM20_P.Saturation_UpperSat_e;
  if (u0_0 > u2) {
    /* Saturate: '<S211>/Saturation' */
    VCM20_B.DYC_deltaTorque = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S211>/Saturation' */
    VCM20_B.DYC_deltaTorque = u1;
  } else {
    /* Saturate: '<S211>/Saturation' */
    VCM20_B.DYC_deltaTorque = u0_0;
  }

  /* End of Saturate: '<S211>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Switch: '<S132>/Switch1' incorporates:
     *  Constant: '<Root>/MC1_sw'
     */
    if (VCM20_P.MC1_sw_Value > VCM20_P.Switch1_Threshold_p4) {
      /* Switch: '<S132>/Switch1' */
      VCM20_B.Switch1_o = VCM20_B.SFunction1_o6_f;
    } else {
      /* Switch: '<S132>/Switch1' incorporates:
       *  Constant: '<S132>/Constant'
       */
      VCM20_B.Switch1_o = VCM20_P.Constant_Value_b;
    }

    /* End of Switch: '<S132>/Switch1' */

    /* Switch: '<S132>/Switch2' incorporates:
     *  Constant: '<Root>/MC2_sw'
     */
    if (VCM20_P.MC2_sw_Value > VCM20_P.Switch2_Threshold_d) {
      /* Switch: '<S132>/Switch2' */
      VCM20_B.Switch2_j = VCM20_B.SFunction1_o6_n;
    } else {
      /* Switch: '<S132>/Switch2' incorporates:
       *  Constant: '<S132>/Constant'
       */
      VCM20_B.Switch2_j = VCM20_P.Constant_Value_b;
    }

    /* End of Switch: '<S132>/Switch2' */

    /* Switch: '<S132>/Switch3' incorporates:
     *  Constant: '<Root>/MC3_sw'
     */
    if (VCM20_P.MC3_sw_Value > VCM20_P.Switch3_Threshold_b) {
      /* Switch: '<S132>/Switch3' */
      VCM20_B.Switch3_k = VCM20_B.SFunction1_o6_h;
    } else {
      /* Switch: '<S132>/Switch3' incorporates:
       *  Constant: '<S132>/Constant'
       */
      VCM20_B.Switch3_k = VCM20_P.Constant_Value_b;
    }

    /* End of Switch: '<S132>/Switch3' */

    /* Switch: '<S132>/Switch4' incorporates:
     *  Constant: '<Root>/MC4_sw'
     */
    if (VCM20_P.MC4_sw_Value > VCM20_P.Switch4_Threshold_f) {
      /* Switch: '<S132>/Switch4' */
      VCM20_B.Switch4_d = VCM20_B.SFunction1_o6_j;
    } else {
      /* Switch: '<S132>/Switch4' incorporates:
       *  Constant: '<S132>/Constant'
       */
      VCM20_B.Switch4_d = VCM20_P.Constant_Value_b;
    }

    /* End of Switch: '<S132>/Switch4' */

    /* Logic: '<S116>/All_enable' */
    VCM20_B.All_enable = ((VCM20_B.Switch1_o != 0.0) && (VCM20_B.Switch2_j !=
      0.0) && (VCM20_B.Switch3_k != 0.0) && (VCM20_B.Switch4_d != 0.0));

    /* Logic: '<S212>/Logical Operator' incorporates:
     *  Constant: '<Root>/DYCOn'
     */
    VCM20_B.LogicalOperator_k = ((VCM20_P.DYCOn_Value != 0.0) &&
      VCM20_B.All_enable);
  }

  /* Switch: '<S212>/Switch' */
  if (VCM20_B.LogicalOperator_k) {
    /* Sum: '<S212>/Add' */
    VCM20_B.torque_lp = VCM20_B.RRtorque + VCM20_B.DYC_deltaTorque;

    /* Switch: '<S212>/Switch' */
    VCM20_B.Switch_g = VCM20_B.torque_lp;
  } else {
    /* Switch: '<S212>/Switch' */
    VCM20_B.Switch_g = VCM20_B.RRtorque;
  }

  /* End of Switch: '<S212>/Switch' */

  /* RelationalOperator: '<S273>/LowerRelop1' */
  VCM20_B.LowerRelop1_e = (VCM20_B.Switch_g > VCM20_B.Switch_ci);

  /* RelationalOperator: '<S273>/UpperRelop' */
  VCM20_B.UpperRelop_n = (VCM20_B.Switch_g < VCM20_B.Saturation1);

  /* Switch: '<S273>/Switch' */
  if (VCM20_B.UpperRelop_n) {
    /* Switch: '<S273>/Switch' */
    VCM20_B.Switch_f = VCM20_B.Saturation1;
  } else {
    /* Switch: '<S273>/Switch' */
    VCM20_B.Switch_f = VCM20_B.Switch_g;
  }

  /* End of Switch: '<S273>/Switch' */

  /* Switch: '<S273>/Switch2' */
  if (VCM20_B.LowerRelop1_e) {
    /* Switch: '<S273>/Switch2' */
    VCM20_B.Switch2_h = VCM20_B.Switch_ci;
  } else {
    /* Switch: '<S273>/Switch2' */
    VCM20_B.Switch2_h = VCM20_B.Switch_f;
  }

  /* End of Switch: '<S273>/Switch2' */

  /* RelationalOperator: '<S305>/Compare' incorporates:
   *  Constant: '<S305>/Constant'
   */
  VCM20_B.Compare_cz = (VCM20_B.Switch2_h == VCM20_P.Constant_Value_jm);

  /* Gain: '<S174>/Gain5' */
  VCM20_B.RLtorque = VCM20_P.Gain5_Gain * VCM20_B.Switch1_g;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Logic: '<S213>/Logical Operator' incorporates:
     *  Constant: '<Root>/DYCOn'
     */
    VCM20_B.LogicalOperator_m = ((VCM20_P.DYCOn_Value != 0.0) &&
      VCM20_B.All_enable);
  }

  /* Switch: '<S213>/Switch' */
  if (VCM20_B.LogicalOperator_m) {
    /* Sum: '<S213>/Add' */
    VCM20_B.torque_f1 = VCM20_B.RLtorque + VCM20_B.DYC_deltaTorque;

    /* Switch: '<S213>/Switch' */
    VCM20_B.Switch_cd = VCM20_B.torque_f1;
  } else {
    /* Switch: '<S213>/Switch' */
    VCM20_B.Switch_cd = VCM20_B.RLtorque;
  }

  /* End of Switch: '<S213>/Switch' */

  /* RelationalOperator: '<S274>/LowerRelop1' */
  VCM20_B.LowerRelop1_ej = (VCM20_B.Switch_cd > VCM20_B.Switch_ci);

  /* RelationalOperator: '<S274>/UpperRelop' */
  VCM20_B.UpperRelop_p = (VCM20_B.Switch_cd < VCM20_B.Saturation1);

  /* Switch: '<S274>/Switch' */
  if (VCM20_B.UpperRelop_p) {
    /* Switch: '<S274>/Switch' */
    VCM20_B.Switch_k = VCM20_B.Saturation1;
  } else {
    /* Switch: '<S274>/Switch' */
    VCM20_B.Switch_k = VCM20_B.Switch_cd;
  }

  /* End of Switch: '<S274>/Switch' */

  /* Switch: '<S274>/Switch2' */
  if (VCM20_B.LowerRelop1_ej) {
    /* Switch: '<S274>/Switch2' */
    VCM20_B.Switch2_i = VCM20_B.Switch_ci;
  } else {
    /* Switch: '<S274>/Switch2' */
    VCM20_B.Switch2_i = VCM20_B.Switch_k;
  }

  /* End of Switch: '<S274>/Switch2' */

  /* RelationalOperator: '<S306>/Compare' incorporates:
   *  Constant: '<S306>/Constant'
   */
  VCM20_B.Compare_p = (VCM20_B.Switch2_i == VCM20_P.Constant_Value_ow);

  /* Sum: '<S279>/Add8' */
  VCM20_B.Add8_l = VCM20_B.Add10 + VCM20_B.Saturation_h;

  /* Sum: '<S279>/Add10' */
  VCM20_B.torque_front_f = VCM20_B.Add8_l - VCM20_B.Saturation3_e;

  /* Switch: '<S174>/Switch' incorporates:
   *  Constant: '<Root>/F//RTrqDistOn'
   */
  if (VCM20_P.FRTrqDistOn_Value > VCM20_P.Switch_Threshold_l) {
    /* Switch: '<S174>/Switch' */
    VCM20_B.Switch_m0 = VCM20_B.torque_front_f;
  } else {
    /* Switch: '<S174>/Switch' */
    VCM20_B.Switch_m0 = VCM20_B.Switch2_f;
  }

  /* End of Switch: '<S174>/Switch' */

  /* Gain: '<S174>/Gain10' */
  VCM20_B.FRtorque = VCM20_P.Gain10_Gain * VCM20_B.Switch_m0;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Logic: '<S214>/Logical Operator' incorporates:
     *  Constant: '<Root>/DYCOn'
     */
    VCM20_B.LogicalOperator_d = ((VCM20_P.DYCOn_Value != 0.0) &&
      VCM20_B.All_enable);
  }

  /* Switch: '<S214>/Switch' */
  if (VCM20_B.LogicalOperator_d) {
    /* Sum: '<S214>/Add' */
    VCM20_B.torque_lj = VCM20_B.FRtorque + VCM20_B.DYC_deltaTorque;

    /* Switch: '<S214>/Switch' */
    VCM20_B.Switch_ci2 = VCM20_B.torque_lj;
  } else {
    /* Switch: '<S214>/Switch' */
    VCM20_B.Switch_ci2 = VCM20_B.FRtorque;
  }

  /* End of Switch: '<S214>/Switch' */

  /* RelationalOperator: '<S275>/LowerRelop1' */
  VCM20_B.LowerRelop1_g = (VCM20_B.Switch_ci2 > VCM20_B.Switch_ci);

  /* RelationalOperator: '<S275>/UpperRelop' */
  VCM20_B.UpperRelop_f5 = (VCM20_B.Switch_ci2 < VCM20_B.Saturation1);

  /* Switch: '<S275>/Switch' */
  if (VCM20_B.UpperRelop_f5) {
    /* Switch: '<S275>/Switch' */
    VCM20_B.Switch_h = VCM20_B.Saturation1;
  } else {
    /* Switch: '<S275>/Switch' */
    VCM20_B.Switch_h = VCM20_B.Switch_ci2;
  }

  /* End of Switch: '<S275>/Switch' */

  /* Switch: '<S275>/Switch2' */
  if (VCM20_B.LowerRelop1_g) {
    /* Switch: '<S275>/Switch2' */
    VCM20_B.Switch2_d = VCM20_B.Switch_ci;
  } else {
    /* Switch: '<S275>/Switch2' */
    VCM20_B.Switch2_d = VCM20_B.Switch_h;
  }

  /* End of Switch: '<S275>/Switch2' */

  /* RelationalOperator: '<S307>/Compare' incorporates:
   *  Constant: '<S307>/Constant'
   */
  VCM20_B.Compare_p1 = (VCM20_B.Switch2_d == VCM20_P.Constant_Value_ct);

  /* Gain: '<S174>/Gain11' */
  VCM20_B.FLtorque = VCM20_P.Gain11_Gain * VCM20_B.Switch_m0;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Logic: '<S215>/Logical Operator' incorporates:
     *  Constant: '<Root>/DYCOn'
     */
    VCM20_B.LogicalOperator_ex = ((VCM20_P.DYCOn_Value != 0.0) &&
      VCM20_B.All_enable);
  }

  /* Switch: '<S215>/Switch' */
  if (VCM20_B.LogicalOperator_ex) {
    /* Sum: '<S215>/Add' */
    VCM20_B.torque_ga = VCM20_B.FLtorque + VCM20_B.DYC_deltaTorque;

    /* Switch: '<S215>/Switch' */
    VCM20_B.Switch_e = VCM20_B.torque_ga;
  } else {
    /* Switch: '<S215>/Switch' */
    VCM20_B.Switch_e = VCM20_B.FLtorque;
  }

  /* End of Switch: '<S215>/Switch' */

  /* RelationalOperator: '<S276>/LowerRelop1' */
  VCM20_B.LowerRelop1_b = (VCM20_B.Switch_e > VCM20_B.Switch_ci);

  /* RelationalOperator: '<S276>/UpperRelop' */
  VCM20_B.UpperRelop_ff = (VCM20_B.Switch_e < VCM20_B.Saturation1);

  /* Switch: '<S276>/Switch' */
  if (VCM20_B.UpperRelop_ff) {
    /* Switch: '<S276>/Switch' */
    VCM20_B.Switch_ed = VCM20_B.Saturation1;
  } else {
    /* Switch: '<S276>/Switch' */
    VCM20_B.Switch_ed = VCM20_B.Switch_e;
  }

  /* End of Switch: '<S276>/Switch' */

  /* Switch: '<S276>/Switch2' */
  if (VCM20_B.LowerRelop1_b) {
    /* Switch: '<S276>/Switch2' */
    VCM20_B.Switch2_c = VCM20_B.Switch_ci;
  } else {
    /* Switch: '<S276>/Switch2' */
    VCM20_B.Switch2_c = VCM20_B.Switch_ed;
  }

  /* End of Switch: '<S276>/Switch2' */

  /* RelationalOperator: '<S308>/Compare' incorporates:
   *  Constant: '<S308>/Constant'
   */
  VCM20_B.Compare_j = (VCM20_B.Switch2_c == VCM20_P.Constant_Value_n);

  /* Logic: '<S290>/Logical Operator' */
  VCM20_B.LogicalOperator_l = (VCM20_B.Compare_cz && VCM20_B.Compare_p &&
    VCM20_B.Compare_p1 && VCM20_B.Compare_j);

  /* Logic: '<S290>/Logical Operator1' */
  VCM20_B.LogicalOperator1_f = !VCM20_B.LogicalOperator_l;

  /* Switch: '<S290>/Switch' */
  if (VCM20_B.LogicalOperator1_f) {
    /* Switch: '<S290>/Switch' */
    VCM20_B.Switch_p = VCM20_B.Divide;
  } else {
    /* Switch: '<S290>/Switch' incorporates:
     *  Constant: '<S290>/Constant2'
     */
    VCM20_B.Switch_p = VCM20_P.Constant2_Value_e;
  }

  /* End of Switch: '<S290>/Switch' */

  /* Lookup_n-D: '<S283>/Fd' incorporates:
   *  Switch: '<S290>/Switch'
   */
  VCM20_B.T_loss = look1_binlxpw(VCM20_B.Switch_p, VCM20_P.Fd_bp01Data,
    VCM20_P.Fd_tableData, 8U);
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Sum: '<S295>/Add' incorporates:
     *  Constant: '<S283>/Target SlipRate'
     *  Constant: '<S295>/Constant2'
     */
    VCM20_B.Add_h = VCM20_P.Constant2_Value_o - VCM20_P.TargetSlipRate_Value;

    /* Product: '<S295>/Divide' incorporates:
     *  Constant: '<S283>/Target Acc'
     */
    u1_0 = VCM20_P.g;

    /* Product: '<S295>/Divide' */
    VCM20_B.Divide_h = u1_0 / VCM20_B.Add_h;

    /* Gain: '<S283>/1// tire R1' */
    u1_0 = 1.0 / VCM20_P.Rwf;

    /* Gain: '<S283>/1// tire R1' */
    VCM20_B.omega = u1_0 * VCM20_B.Divide_h;

    /* Gain: '<S283>/wheel Inertia' */
    VCM20_B.wheelInertia = VCM20_P.Iwf * VCM20_B.omega;

    /* Gain: '<S283>/1//Gear' */
    u1_0 = 1.0 / VCM20_P.gear_f;

    /* Gain: '<S283>/1//Gear' */
    VCM20_B.Iomega = u1_0 * VCM20_B.wheelInertia;

    /* Delay: '<S296>/Delay' */
    VCM20_B.Delay = VCM20_DW.Delay_DSTATE;

    /* Gain: '<S296>/FR Load Diff' */
    u1_0 = VCM20_P.M * VCM20_P.hg / VCM20_P.l / 2.0;

    /* Gain: '<S296>/FR Load Diff' */
    VCM20_B.FRLoadDiff = u1_0 * VCM20_B.Delay;

    /* Saturate: '<S296>/Saturation' */
    u1_0 = VCM20_P.W0_R * -9.8;
    u2 = VCM20_P.W0_F * 9.8;
    u0_0 = VCM20_B.FRLoadDiff;
    if (u0_0 > u2) {
      /* Saturate: '<S296>/Saturation' */
      VCM20_B.Saturation_g = u2;
    } else if (u0_0 < u1_0) {
      /* Saturate: '<S296>/Saturation' */
      VCM20_B.Saturation_g = u1_0;
    } else {
      /* Saturate: '<S296>/Saturation' */
      VCM20_B.Saturation_g = u0_0;
    }

    /* End of Saturate: '<S296>/Saturation' */

    /* Sum: '<S296>/Add1' incorporates:
     *  Constant: '<S296>/Front_ini_Load'
     */
    u1_0 = VCM20_P.W0_F * VCM20_P.g;

    /* Sum: '<S296>/Add1' */
    VCM20_B.FrontLoad = u1_0 - VCM20_B.Saturation_g;

    /* Gain: '<S298>/a1*W' */
    u1_0 = VCM20_P.afx_f[1];

    /* Gain: '<S298>/a1*W' */
    VCM20_B.a1W = u1_0 * VCM20_B.FrontLoad;

    /* Sum: '<S298>/Sum' incorporates:
     *  Constant: '<S298>/a2'
     */
    u1_0 = VCM20_P.afx_f[2];

    /* Sum: '<S298>/Sum' */
    VCM20_B.Sum_o = VCM20_B.a1W + u1_0;

    /* Product: '<S298>/D' */
    VCM20_B.D = VCM20_B.Sum_o * VCM20_B.FrontLoad;

    /* Product: '<S298>/Divide' incorporates:
     *  Constant: '<S298>/a4'
     */
    u1_0 = VCM20_P.afx_f[4];

    /* Product: '<S298>/Divide' */
    VCM20_B.Divide_m = VCM20_B.FrontLoad / u1_0;

    /* Trigonometry: '<S298>/atan(W//a4)' */
    VCM20_B.atanWa4 = atan(VCM20_B.Divide_m);

    /* Gain: '<S298>/Gain' */
    VCM20_B.Gain_o = VCM20_P.Gain_Gain_n3 * VCM20_B.atanWa4;

    /* Trigonometry: '<S298>/sin(2*atan(w//a4))' */
    VCM20_B.sin2atanwa4 = sin(VCM20_B.Gain_o);

    /* Gain: '<S298>/BCD' */
    u1_0 = VCM20_P.afx_f[3];

    /* Gain: '<S298>/BCD' */
    VCM20_B.BCD = u1_0 * VCM20_B.sin2atanwa4;

    /* Product: '<S298>/Product2' incorporates:
     *  Constant: '<S298>/C'
     */
    u1_0 = VCM20_P.afx_f[0];

    /* Product: '<S298>/Product2' */
    VCM20_B.Product2 = VCM20_B.D * u1_0;

    /* RelationalOperator: '<S300>/Compare' incorporates:
     *  Constant: '<S300>/Constant'
     */
    VCM20_B.Compare_f = (VCM20_B.Product2 == VCM20_P.Constant_Value_cg);

    /* Switch: '<S298>/Switch' */
    if (VCM20_B.Compare_f) {
      /* Switch: '<S298>/Switch' incorporates:
       *  Constant: '<S298>/Constant1'
       */
      VCM20_B.Switch_hy = VCM20_P.Constant1_Value_i;
    } else {
      /* Switch: '<S298>/Switch' */
      VCM20_B.Switch_hy = VCM20_B.Product2;
    }

    /* End of Switch: '<S298>/Switch' */

    /* Product: '<S298>/B' */
    VCM20_B.B = VCM20_B.BCD / VCM20_B.Switch_hy;

    /* Product: '<S298>/Product' incorporates:
     *  Constant: '<S283>/Target SlipRate'
     */
    VCM20_B.Product_i = VCM20_B.B * VCM20_P.TargetSlipRate_Value;

    /* Trigonometry: '<S298>/atan(B*É¿)' */
    VCM20_B.atanB = atan(VCM20_B.Product_i);

    /* Sum: '<S298>/Sum2' */
    VCM20_B.Sum2 = VCM20_B.Product_i - VCM20_B.atanB;

    /* Gain: '<S298>/a5*W' */
    u1_0 = VCM20_P.afx_f[5];

    /* Gain: '<S298>/a5*W' */
    VCM20_B.a5W = u1_0 * VCM20_B.FrontLoad;

    /* Sum: '<S298>/Sum1' incorporates:
     *  Constant: '<S298>/a6'
     */
    u1_0 = VCM20_P.afx_f[6];

    /* Sum: '<S298>/Sum1' */
    VCM20_B.Sum1_g = VCM20_B.a5W + u1_0;

    /* Product: '<S298>/E*(B*SA-atan(B*SA))' */
    VCM20_B.EBSAatanBSA = VCM20_B.Sum2 * VCM20_B.Sum1_g;

    /* Sum: '<S298>/Sum3' */
    VCM20_B.Sum3_i = VCM20_B.Product_i - VCM20_B.EBSAatanBSA;

    /* Trigonometry: '<S298>/Trigonometric Function' */
    VCM20_B.TrigonometricFunction = atan(VCM20_B.Sum3_i);

    /* Product: '<S298>/Product1' incorporates:
     *  Constant: '<S298>/C'
     */
    u1_0 = VCM20_P.afx_f[0];

    /* Product: '<S298>/Product1' */
    VCM20_B.Product1 = VCM20_B.TrigonometricFunction * u1_0;

    /* Trigonometry: '<S298>/Trigonometric Function2' */
    VCM20_B.TrigonometricFunction2 = sin(VCM20_B.Product1);

    /* Product: '<S298>/+D*cos(delta)' */
    VCM20_B.Dcosdelta = VCM20_B.D * VCM20_B.TrigonometricFunction2;

    /* Gain: '<S298>/Gain1' */
    VCM20_B.Gain1_a = VCM20_P.mu_tire_F * VCM20_B.Dcosdelta;

    /* Gain: '<S283>/Rtire' */
    VCM20_B.Rtire = VCM20_P.Rwf * VCM20_B.Gain1_a;

    /* Gain: '<S283>/1//Gear1' */
    u1_0 = 1.0 / VCM20_P.gear_f;

    /* Gain: '<S283>/1//Gear1' */
    VCM20_B.rFFr = u1_0 * VCM20_B.Rtire;
  }

  /* Sum: '<S283>/Add' */
  VCM20_B.Add_l = VCM20_B.T_loss + VCM20_B.Iomega;

  /* Sum: '<S283>/Add1' */
  VCM20_B.Add1_jx = VCM20_B.Add_l + VCM20_B.rFFr;

  /* Sum: '<S294>/Add3' */
  VCM20_B.Add3_h = VCM20_B.Switch_p - VCM20_B.SFunction1_o9_g;

  /* Gain: '<S294>/P gain' */
  VCM20_B.Pgain = VCM20_P.Pgain_Gain * VCM20_B.Add3_h;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* DiscreteIntegrator: '<S294>/Discrete-Time Integrator' */
    VCM20_B.DiscreteTimeIntegrator = VCM20_DW.DiscreteTimeIntegrator_DSTATE;
  }

  /* Sum: '<S294>/Add2' */
  VCM20_B.FB = VCM20_B.Pgain + VCM20_B.DiscreteTimeIntegrator;

  /* Sum: '<S294>/Add' */
  VCM20_B.torque = VCM20_B.Add1_jx + VCM20_B.FB;

  /* RelationalOperator: '<S317>/LowerRelop1' */
  VCM20_B.LowerRelop1_bx = (VCM20_B.torque > VCM20_B.Saturation);

  /* RelationalOperator: '<S317>/UpperRelop' incorporates:
   *  Constant: '<S294>/Constant'
   */
  VCM20_B.UpperRelop_l = (VCM20_B.torque < VCM20_P.Constant_Value_fi);

  /* Switch: '<S317>/Switch' */
  if (VCM20_B.UpperRelop_l) {
    /* Switch: '<S317>/Switch' incorporates:
     *  Constant: '<S294>/Constant'
     */
    VCM20_B.Switch_b = VCM20_P.Constant_Value_fi;
  } else {
    /* Switch: '<S317>/Switch' */
    VCM20_B.Switch_b = VCM20_B.torque;
  }

  /* End of Switch: '<S317>/Switch' */

  /* Switch: '<S317>/Switch2' */
  if (VCM20_B.LowerRelop1_bx) {
    /* Switch: '<S317>/Switch2' */
    VCM20_B.Switch2_k = VCM20_B.Saturation;
  } else {
    /* Switch: '<S317>/Switch2' */
    VCM20_B.Switch2_k = VCM20_B.Switch_b;
  }

  /* End of Switch: '<S317>/Switch2' */

  /* Switch: '<S294>/Switch' incorporates:
   *  Constant: '<Root>/TCOn'
   */
  if (VCM20_P.TCOn_Value > VCM20_P.Switch_Threshold_g) {
    /* Switch: '<S294>/Switch' */
    VCM20_B.torque_c = VCM20_B.Switch2_k;
  } else {
    /* Switch: '<S294>/Switch' */
    VCM20_B.torque_c = VCM20_B.Saturation;
  }

  /* End of Switch: '<S294>/Switch' */

  /* RelationalOperator: '<S289>/LowerRelop1' */
  VCM20_B.LowerRelop1_pb = (VCM20_B.Switch2_c > VCM20_B.torque_c);
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Gain: '<S284>/toTire_rpm' */
    VCM20_B.toTire_rpm_l = VCM20_P.toTire_rpm_Gain * VCM20_B.VELOCITY_X_p;

    /* Gain: '<S284>/Gear' */
    VCM20_B.Gear_d = VCM20_P.Gear_Gain * VCM20_B.toTire_rpm_l;

    /* Sum: '<S302>/Add' incorporates:
     *  Constant: '<S284>/Front sliprate_inBrake_ref'
     *  Constant: '<S302>/Constant2'
     */
    VCM20_B.Add_m0 = VCM20_P.Constant2_Value_l -
      VCM20_P.Frontsliprate_inBrake_ref_Value;

    /* Product: '<S302>/Divide' */
    VCM20_B.Divide_b = VCM20_B.Gear_d / VCM20_B.Add_m0;

    /* Gain: '<S284>/FrontGain' */
    VCM20_B.FrontGain = VCM20_P.FrontGain_Gain * VCM20_B.Divide_b;

    /* RelationalOperator: '<S301>/LowerRelop1' */
    VCM20_B.LowerRelop1_d = (VCM20_B.FrontGain > VCM20_B.Saturation2);

    /* RelationalOperator: '<S301>/UpperRelop' incorporates:
     *  Constant: '<S284>/Constant'
     */
    VCM20_B.UpperRelop_m = (VCM20_B.FrontGain < VCM20_P.Constant_Value_l);

    /* Switch: '<S301>/Switch' */
    if (VCM20_B.UpperRelop_m) {
      /* Switch: '<S301>/Switch' incorporates:
       *  Constant: '<S284>/Constant'
       */
      VCM20_B.Switch_bt = VCM20_P.Constant_Value_l;
    } else {
      /* Switch: '<S301>/Switch' */
      VCM20_B.Switch_bt = VCM20_B.FrontGain;
    }

    /* End of Switch: '<S301>/Switch' */

    /* Switch: '<S301>/Switch2' */
    if (VCM20_B.LowerRelop1_d) {
      /* Switch: '<S301>/Switch2' */
      VCM20_B.Switch2_o = VCM20_B.Saturation2;
    } else {
      /* Switch: '<S301>/Switch2' */
      VCM20_B.Switch2_o = VCM20_B.Switch_bt;
    }

    /* End of Switch: '<S301>/Switch2' */

    /* Sum: '<S294>/Add1' */
    VCM20_B.Add1_g = VCM20_B.SFunction1_o9_g - VCM20_B.Switch2_o;

    /* RelationalOperator: '<S318>/LowerRelop1' incorporates:
     *  Constant: '<S294>/Constant'
     */
    VCM20_B.LowerRelop1_ez = (VCM20_B.Add1_g > VCM20_P.Constant_Value_fi);

    /* Gain: '<S294>/Gain1' */
    VCM20_B.Gain1_l = VCM20_P.Gain1_Gain_dw * VCM20_B.Saturation1;

    /* RelationalOperator: '<S318>/UpperRelop' */
    VCM20_B.UpperRelop_fr = (VCM20_B.Add1_g < VCM20_B.Gain1_l);

    /* Switch: '<S318>/Switch' */
    if (VCM20_B.UpperRelop_fr) {
      /* Switch: '<S318>/Switch' */
      VCM20_B.Switch_oi = VCM20_B.Gain1_l;
    } else {
      /* Switch: '<S318>/Switch' */
      VCM20_B.Switch_oi = VCM20_B.Add1_g;
    }

    /* End of Switch: '<S318>/Switch' */

    /* Switch: '<S318>/Switch2' */
    if (VCM20_B.LowerRelop1_ez) {
      /* Switch: '<S318>/Switch2' incorporates:
       *  Constant: '<S294>/Constant'
       */
      VCM20_B.Switch2_e1 = VCM20_P.Constant_Value_fi;
    } else {
      /* Switch: '<S318>/Switch2' */
      VCM20_B.Switch2_e1 = VCM20_B.Switch_oi;
    }

    /* End of Switch: '<S318>/Switch2' */

    /* Switch: '<S294>/Switch1' incorporates:
     *  Constant: '<Root>/ABSOn'
     */
    if (VCM20_P.ABSOn_Value > VCM20_P.Switch1_Threshold_f) {
      /* Switch: '<S294>/Switch1' */
      VCM20_B.Switch1_d = VCM20_B.Switch2_e1;
    } else {
      /* Gain: '<S294>/Gain' */
      VCM20_B.Gain_h = VCM20_P.Gain_Gain_i * VCM20_B.Saturation1;

      /* Switch: '<S294>/Switch1' */
      VCM20_B.Switch1_d = VCM20_B.Gain_h;
    }

    /* End of Switch: '<S294>/Switch1' */
  }

  /* RelationalOperator: '<S289>/UpperRelop' */
  VCM20_B.UpperRelop_g = (VCM20_B.Switch2_c < VCM20_B.Switch1_d);

  /* Switch: '<S289>/Switch' */
  if (VCM20_B.UpperRelop_g) {
    /* Switch: '<S289>/Switch' */
    VCM20_B.Switch_fd = VCM20_B.Switch1_d;
  } else {
    /* Switch: '<S289>/Switch' */
    VCM20_B.Switch_fd = VCM20_B.Switch2_c;
  }

  /* End of Switch: '<S289>/Switch' */

  /* Switch: '<S289>/Switch2' */
  if (VCM20_B.LowerRelop1_pb) {
    /* Switch: '<S289>/Switch2' */
    VCM20_B.Switch2_ot = VCM20_B.torque_c;
  } else {
    /* Switch: '<S289>/Switch2' */
    VCM20_B.Switch2_ot = VCM20_B.Switch_fd;
  }

  /* End of Switch: '<S289>/Switch2' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Sum: '<S296>/Add2' incorporates:
     *  Constant: '<S296>/Rear_ini_Load2'
     */
    u1_0 = VCM20_P.W0_R * VCM20_P.g;

    /* Sum: '<S296>/Add2' */
    VCM20_B.RearLoad = VCM20_B.Saturation_g + u1_0;

    /* Gain: '<S297>/a1*W' */
    u1_0 = VCM20_P.afx_f[1];

    /* Gain: '<S297>/a1*W' */
    VCM20_B.a1W_p = u1_0 * VCM20_B.RearLoad;

    /* Sum: '<S297>/Sum' incorporates:
     *  Constant: '<S297>/a2'
     */
    u1_0 = VCM20_P.afx_f[2];

    /* Sum: '<S297>/Sum' */
    VCM20_B.Sum_k = VCM20_B.a1W_p + u1_0;

    /* Product: '<S297>/D' */
    VCM20_B.D_h = VCM20_B.Sum_k * VCM20_B.RearLoad;

    /* Product: '<S297>/Divide' incorporates:
     *  Constant: '<S297>/a4'
     */
    u1_0 = VCM20_P.afx_f[4];

    /* Product: '<S297>/Divide' */
    VCM20_B.Divide_e = VCM20_B.RearLoad / u1_0;

    /* Trigonometry: '<S297>/atan(W//a4)' */
    VCM20_B.atanWa4_j = atan(VCM20_B.Divide_e);

    /* Gain: '<S297>/Gain' */
    VCM20_B.Gain_e = VCM20_P.Gain_Gain_f * VCM20_B.atanWa4_j;

    /* Trigonometry: '<S297>/sin(2*atan(w//a4))' */
    VCM20_B.sin2atanwa4_k = sin(VCM20_B.Gain_e);

    /* Gain: '<S297>/BCD' */
    u1_0 = VCM20_P.afx_f[3];

    /* Gain: '<S297>/BCD' */
    VCM20_B.BCD_a = u1_0 * VCM20_B.sin2atanwa4_k;

    /* Product: '<S297>/Product2' incorporates:
     *  Constant: '<S297>/C'
     */
    u1_0 = VCM20_P.afx_f[0];

    /* Product: '<S297>/Product2' */
    VCM20_B.Product2_p = VCM20_B.D_h * u1_0;

    /* RelationalOperator: '<S299>/Compare' incorporates:
     *  Constant: '<S299>/Constant'
     */
    VCM20_B.Compare_jm = (VCM20_B.Product2_p == VCM20_P.Constant_Value_f);

    /* Switch: '<S297>/Switch' */
    if (VCM20_B.Compare_jm) {
      /* Switch: '<S297>/Switch' incorporates:
       *  Constant: '<S297>/Constant1'
       */
      VCM20_B.Switch_bc = VCM20_P.Constant1_Value_b;
    } else {
      /* Switch: '<S297>/Switch' */
      VCM20_B.Switch_bc = VCM20_B.Product2_p;
    }

    /* End of Switch: '<S297>/Switch' */

    /* Product: '<S297>/B' */
    VCM20_B.B_p = VCM20_B.BCD_a / VCM20_B.Switch_bc;

    /* Product: '<S297>/Product' incorporates:
     *  Constant: '<S283>/Target SlipRate'
     */
    VCM20_B.Product_e = VCM20_B.B_p * VCM20_P.TargetSlipRate_Value;

    /* Trigonometry: '<S297>/atan(B*É¿)' */
    VCM20_B.atanB_i = atan(VCM20_B.Product_e);

    /* Sum: '<S297>/Sum2' */
    VCM20_B.Sum2_i = VCM20_B.Product_e - VCM20_B.atanB_i;

    /* Gain: '<S297>/a5*W' */
    u1_0 = VCM20_P.afx_f[5];

    /* Gain: '<S297>/a5*W' */
    VCM20_B.a5W_d = u1_0 * VCM20_B.RearLoad;

    /* Sum: '<S297>/Sum1' incorporates:
     *  Constant: '<S297>/a6'
     */
    u1_0 = VCM20_P.afx_f[6];

    /* Sum: '<S297>/Sum1' */
    VCM20_B.Sum1_f = VCM20_B.a5W_d + u1_0;

    /* Product: '<S297>/E*(B*SA-atan(B*SA))' */
    VCM20_B.EBSAatanBSA_b = VCM20_B.Sum2_i * VCM20_B.Sum1_f;

    /* Sum: '<S297>/Sum3' */
    VCM20_B.Sum3_p = VCM20_B.Product_e - VCM20_B.EBSAatanBSA_b;

    /* Trigonometry: '<S297>/Trigonometric Function' */
    VCM20_B.TrigonometricFunction_k = atan(VCM20_B.Sum3_p);

    /* Product: '<S297>/Product1' incorporates:
     *  Constant: '<S297>/C'
     */
    u1_0 = VCM20_P.afx_f[0];

    /* Product: '<S297>/Product1' */
    VCM20_B.Product1_n = VCM20_B.TrigonometricFunction_k * u1_0;

    /* Trigonometry: '<S297>/Trigonometric Function2' */
    VCM20_B.TrigonometricFunction2_g = sin(VCM20_B.Product1_n);

    /* Product: '<S297>/+D*cos(delta)' */
    VCM20_B.Dcosdelta_i = VCM20_B.D_h * VCM20_B.TrigonometricFunction2_g;

    /* Gain: '<S297>/Gain1' */
    VCM20_B.Gain1_cs = VCM20_P.mu_tire_F * VCM20_B.Dcosdelta_i;

    /* Gain: '<S283>/Rtire1' */
    u1_0 = 1.0 / VCM20_P.gear_r;

    /* Gain: '<S283>/Rtire1' */
    VCM20_B.Rtire1 = u1_0 * VCM20_B.Gain1_cs;

    /* Gain: '<S283>/1//Gear2' */
    u1_0 = 1.0 / VCM20_P.gear_r;

    /* Gain: '<S283>/1//Gear2' */
    VCM20_B.rFRr = u1_0 * VCM20_B.Rtire1;

    /* DiscreteIntegrator: '<S291>/Discrete-Time Integrator' */
    VCM20_B.DiscreteTimeIntegrator_h = VCM20_DW.DiscreteTimeIntegrator_DSTATE_e;
  }

  /* Sum: '<S283>/Add2' */
  VCM20_B.Add2_g = VCM20_B.Add_l + VCM20_B.rFRr;

  /* Sum: '<S291>/Add3' */
  VCM20_B.Add3_k = VCM20_B.Switch_p - VCM20_B.SFunction1_o9_n;

  /* Gain: '<S291>/P gain' */
  VCM20_B.Pgain_i = VCM20_P.Pgain_Gain_p * VCM20_B.Add3_k;

  /* Sum: '<S291>/Add2' */
  VCM20_B.FB_f = VCM20_B.Pgain_i + VCM20_B.DiscreteTimeIntegrator_h;

  /* Sum: '<S291>/Add' */
  VCM20_B.torque_b = VCM20_B.Add2_g + VCM20_B.FB_f;

  /* RelationalOperator: '<S311>/LowerRelop1' */
  VCM20_B.LowerRelop1_l = (VCM20_B.torque_b > VCM20_B.Saturation);

  /* RelationalOperator: '<S311>/UpperRelop' incorporates:
   *  Constant: '<S291>/Constant'
   */
  VCM20_B.UpperRelop_c = (VCM20_B.torque_b < VCM20_P.Constant_Value_bj);

  /* Switch: '<S311>/Switch' */
  if (VCM20_B.UpperRelop_c) {
    /* Switch: '<S311>/Switch' incorporates:
     *  Constant: '<S291>/Constant'
     */
    VCM20_B.Switch_a = VCM20_P.Constant_Value_bj;
  } else {
    /* Switch: '<S311>/Switch' */
    VCM20_B.Switch_a = VCM20_B.torque_b;
  }

  /* End of Switch: '<S311>/Switch' */

  /* Switch: '<S311>/Switch2' */
  if (VCM20_B.LowerRelop1_l) {
    /* Switch: '<S311>/Switch2' */
    VCM20_B.Switch2_jk = VCM20_B.Saturation;
  } else {
    /* Switch: '<S311>/Switch2' */
    VCM20_B.Switch2_jk = VCM20_B.Switch_a;
  }

  /* End of Switch: '<S311>/Switch2' */

  /* Switch: '<S291>/Switch' incorporates:
   *  Constant: '<Root>/TCOn'
   */
  if (VCM20_P.TCOn_Value > VCM20_P.Switch_Threshold_f) {
    /* Switch: '<S291>/Switch' */
    VCM20_B.torque_f = VCM20_B.Switch2_jk;
  } else {
    /* Switch: '<S291>/Switch' */
    VCM20_B.torque_f = VCM20_B.Saturation;
  }

  /* End of Switch: '<S291>/Switch' */

  /* RelationalOperator: '<S286>/LowerRelop1' */
  VCM20_B.LowerRelop1_m = (VCM20_B.Switch2_h > VCM20_B.torque_f);
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Gain: '<S285>/toTire_rpm' */
    VCM20_B.toTire_rpm_a = VCM20_P.toTire_rpm_Gain_l * VCM20_B.VELOCITY_X_p;

    /* Gain: '<S285>/Gear' */
    VCM20_B.Gear_n = VCM20_P.Gear_Gain_m * VCM20_B.toTire_rpm_a;

    /* Sum: '<S304>/Add' incorporates:
     *  Constant: '<S285>/Rear sliprate_inBrake_ref'
     *  Constant: '<S304>/Constant2'
     */
    VCM20_B.Add_ib = VCM20_P.Constant2_Value_g -
      VCM20_P.Rearsliprate_inBrake_ref_Value;

    /* Product: '<S304>/Divide' */
    VCM20_B.Divide_k = VCM20_B.Gear_n / VCM20_B.Add_ib;

    /* RelationalOperator: '<S303>/LowerRelop1' */
    VCM20_B.LowerRelop1_mw = (VCM20_B.Divide_k > VCM20_B.Saturation2);

    /* RelationalOperator: '<S303>/UpperRelop' incorporates:
     *  Constant: '<S285>/Constant'
     */
    VCM20_B.UpperRelop_pt = (VCM20_B.Divide_k < VCM20_P.Constant_Value_m0);

    /* Switch: '<S303>/Switch' */
    if (VCM20_B.UpperRelop_pt) {
      /* Switch: '<S303>/Switch' incorporates:
       *  Constant: '<S285>/Constant'
       */
      VCM20_B.Switch_pe = VCM20_P.Constant_Value_m0;
    } else {
      /* Switch: '<S303>/Switch' */
      VCM20_B.Switch_pe = VCM20_B.Divide_k;
    }

    /* End of Switch: '<S303>/Switch' */

    /* Switch: '<S303>/Switch2' */
    if (VCM20_B.LowerRelop1_mw) {
      /* Switch: '<S303>/Switch2' */
      VCM20_B.Switch2_m = VCM20_B.Saturation2;
    } else {
      /* Switch: '<S303>/Switch2' */
      VCM20_B.Switch2_m = VCM20_B.Switch_pe;
    }

    /* End of Switch: '<S303>/Switch2' */

    /* Sum: '<S291>/Add1' */
    VCM20_B.Add1_i = VCM20_B.SFunction1_o9_n - VCM20_B.Switch2_m;

    /* RelationalOperator: '<S312>/LowerRelop1' incorporates:
     *  Constant: '<S291>/Constant'
     */
    VCM20_B.LowerRelop1_f = (VCM20_B.Add1_i > VCM20_P.Constant_Value_bj);

    /* Gain: '<S291>/Gain1' */
    VCM20_B.Gain1_i = VCM20_P.Gain1_Gain_j * VCM20_B.Saturation1;

    /* RelationalOperator: '<S312>/UpperRelop' */
    VCM20_B.UpperRelop_d = (VCM20_B.Add1_i < VCM20_B.Gain1_i);

    /* Switch: '<S312>/Switch' */
    if (VCM20_B.UpperRelop_d) {
      /* Switch: '<S312>/Switch' */
      VCM20_B.Switch_cr = VCM20_B.Gain1_i;
    } else {
      /* Switch: '<S312>/Switch' */
      VCM20_B.Switch_cr = VCM20_B.Add1_i;
    }

    /* End of Switch: '<S312>/Switch' */

    /* Switch: '<S312>/Switch2' */
    if (VCM20_B.LowerRelop1_f) {
      /* Switch: '<S312>/Switch2' incorporates:
       *  Constant: '<S291>/Constant'
       */
      VCM20_B.Switch2_jy = VCM20_P.Constant_Value_bj;
    } else {
      /* Switch: '<S312>/Switch2' */
      VCM20_B.Switch2_jy = VCM20_B.Switch_cr;
    }

    /* End of Switch: '<S312>/Switch2' */

    /* Switch: '<S291>/Switch1' incorporates:
     *  Constant: '<Root>/ABSOn'
     */
    if (VCM20_P.ABSOn_Value > VCM20_P.Switch1_Threshold_p2) {
      /* Switch: '<S291>/Switch1' */
      VCM20_B.Switch1_e = VCM20_B.Switch2_jy;
    } else {
      /* Gain: '<S291>/Gain' */
      VCM20_B.Gain_f = VCM20_P.Gain_Gain_h * VCM20_B.Saturation1;

      /* Switch: '<S291>/Switch1' */
      VCM20_B.Switch1_e = VCM20_B.Gain_f;
    }

    /* End of Switch: '<S291>/Switch1' */
  }

  /* RelationalOperator: '<S286>/UpperRelop' */
  VCM20_B.UpperRelop_gu = (VCM20_B.Switch2_h < VCM20_B.Switch1_e);

  /* Switch: '<S286>/Switch' */
  if (VCM20_B.UpperRelop_gu) {
    /* Switch: '<S286>/Switch' */
    VCM20_B.Switch_bb = VCM20_B.Switch1_e;
  } else {
    /* Switch: '<S286>/Switch' */
    VCM20_B.Switch_bb = VCM20_B.Switch2_h;
  }

  /* End of Switch: '<S286>/Switch' */

  /* Switch: '<S286>/Switch2' */
  if (VCM20_B.LowerRelop1_m) {
    /* Switch: '<S286>/Switch2' */
    VCM20_B.Switch2_g = VCM20_B.torque_f;
  } else {
    /* Switch: '<S286>/Switch2' */
    VCM20_B.Switch2_g = VCM20_B.Switch_bb;
  }

  /* End of Switch: '<S286>/Switch2' */

  /* Saturate: '<S186>/Saturation' */
  u0_0 = VCM20_B.Switch2_g;
  u1 = VCM20_P.Saturation_LowerSat_i;
  u2 = VCM20_P.Saturation_UpperSat_l;
  if (u0_0 > u2) {
    /* Saturate: '<S186>/Saturation' */
    VCM20_B.Saturation_a = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S186>/Saturation' */
    VCM20_B.Saturation_a = u1;
  } else {
    /* Saturate: '<S186>/Saturation' */
    VCM20_B.Saturation_a = u0_0;
  }

  /* End of Saturate: '<S186>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S186>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_n;
    u1 = VCM20_P.Saturation1_LowerSat_k;
    u2 = VCM20_P.Saturation1_UpperSat_d5;
    if (u0_0 > u2) {
      /* Saturate: '<S186>/Saturation1' */
      VCM20_B.Saturation1_i = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S186>/Saturation1' */
      VCM20_B.Saturation1_i = u1;
    } else {
      /* Saturate: '<S186>/Saturation1' */
      VCM20_B.Saturation1_i = u0_0;
    }

    /* End of Saturate: '<S186>/Saturation1' */

    /* DiscreteIntegrator: '<S292>/Discrete-Time Integrator' */
    VCM20_B.DiscreteTimeIntegrator_i = VCM20_DW.DiscreteTimeIntegrator_DSTATE_m;
  }

  /* Product: '<S186>/Product' */
  VCM20_B.FRMotorPower = VCM20_B.Saturation_a * VCM20_B.Saturation1_i;

  /* Gain: '<S186>/Gain' */
  VCM20_B.FRMotorPower_W = VCM20_P.Gain_Gain_e * VCM20_B.FRMotorPower;

  /* Lookup_n-D: '<S186>/2-D Lookup Table' incorporates:
   *  Saturate: '<S186>/Saturation'
   *  Saturate: '<S186>/Saturation1'
   */
  VCM20_B.FREfficency = look2_binlxpw(VCM20_B.Saturation_a,
    VCM20_B.Saturation1_i, VCM20_P.uDLookupTable_bp01Data,
    VCM20_P.uDLookupTable_bp02Data, VCM20_P.uDLookupTable_tableData,
    VCM20_P.uDLookupTable_maxIndex, 11U);

  /* Product: '<S186>/Divide' */
  VCM20_B.FRPowerConsumption = VCM20_B.FRMotorPower_W / VCM20_B.FREfficency;

  /* Sum: '<S292>/Add3' */
  VCM20_B.Add3_k0 = VCM20_B.Switch_p - VCM20_B.SFunction1_o9_j;

  /* Gain: '<S292>/P gain' */
  VCM20_B.Pgain_a = VCM20_P.Pgain_Gain_h * VCM20_B.Add3_k0;

  /* Sum: '<S292>/Add2' */
  VCM20_B.FB_d = VCM20_B.Pgain_a + VCM20_B.DiscreteTimeIntegrator_i;

  /* Sum: '<S292>/Add' */
  VCM20_B.torque_g = VCM20_B.Add2_g + VCM20_B.FB_d;

  /* RelationalOperator: '<S313>/LowerRelop1' */
  VCM20_B.LowerRelop1_k = (VCM20_B.torque_g > VCM20_B.Saturation);

  /* RelationalOperator: '<S313>/UpperRelop' incorporates:
   *  Constant: '<S292>/Constant'
   */
  VCM20_B.UpperRelop_b = (VCM20_B.torque_g < VCM20_P.Constant_Value_lq);

  /* Switch: '<S313>/Switch' */
  if (VCM20_B.UpperRelop_b) {
    /* Switch: '<S313>/Switch' incorporates:
     *  Constant: '<S292>/Constant'
     */
    VCM20_B.Switch_m0p = VCM20_P.Constant_Value_lq;
  } else {
    /* Switch: '<S313>/Switch' */
    VCM20_B.Switch_m0p = VCM20_B.torque_g;
  }

  /* End of Switch: '<S313>/Switch' */

  /* Switch: '<S313>/Switch2' */
  if (VCM20_B.LowerRelop1_k) {
    /* Switch: '<S313>/Switch2' */
    VCM20_B.Switch2_az = VCM20_B.Saturation;
  } else {
    /* Switch: '<S313>/Switch2' */
    VCM20_B.Switch2_az = VCM20_B.Switch_m0p;
  }

  /* End of Switch: '<S313>/Switch2' */

  /* Switch: '<S292>/Switch' incorporates:
   *  Constant: '<Root>/TCOn'
   */
  if (VCM20_P.TCOn_Value > VCM20_P.Switch_Threshold_o) {
    /* Switch: '<S292>/Switch' */
    VCM20_B.torque_p = VCM20_B.Switch2_az;
  } else {
    /* Switch: '<S292>/Switch' */
    VCM20_B.torque_p = VCM20_B.Saturation;
  }

  /* End of Switch: '<S292>/Switch' */

  /* RelationalOperator: '<S287>/LowerRelop1' */
  VCM20_B.LowerRelop1_mx = (VCM20_B.Switch2_i > VCM20_B.torque_p);
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Sum: '<S292>/Add1' */
    VCM20_B.Add1_n = VCM20_B.SFunction1_o9_j - VCM20_B.Switch2_m;

    /* RelationalOperator: '<S314>/LowerRelop1' incorporates:
     *  Constant: '<S292>/Constant'
     */
    VCM20_B.LowerRelop1_n1 = (VCM20_B.Add1_n > VCM20_P.Constant_Value_lq);

    /* Gain: '<S292>/Gain1' */
    VCM20_B.Gain1_b = VCM20_P.Gain1_Gain_c * VCM20_B.Saturation1;

    /* RelationalOperator: '<S314>/UpperRelop' */
    VCM20_B.UpperRelop_o = (VCM20_B.Add1_n < VCM20_B.Gain1_b);

    /* Switch: '<S314>/Switch' */
    if (VCM20_B.UpperRelop_o) {
      /* Switch: '<S314>/Switch' */
      VCM20_B.Switch_n = VCM20_B.Gain1_b;
    } else {
      /* Switch: '<S314>/Switch' */
      VCM20_B.Switch_n = VCM20_B.Add1_n;
    }

    /* End of Switch: '<S314>/Switch' */

    /* Switch: '<S314>/Switch2' */
    if (VCM20_B.LowerRelop1_n1) {
      /* Switch: '<S314>/Switch2' incorporates:
       *  Constant: '<S292>/Constant'
       */
      VCM20_B.Switch2_p = VCM20_P.Constant_Value_lq;
    } else {
      /* Switch: '<S314>/Switch2' */
      VCM20_B.Switch2_p = VCM20_B.Switch_n;
    }

    /* End of Switch: '<S314>/Switch2' */

    /* Switch: '<S292>/Switch1' incorporates:
     *  Constant: '<Root>/ABSOn'
     */
    if (VCM20_P.ABSOn_Value > VCM20_P.Switch1_Threshold_d) {
      /* Switch: '<S292>/Switch1' */
      VCM20_B.Switch1_p = VCM20_B.Switch2_p;
    } else {
      /* Gain: '<S292>/Gain' */
      VCM20_B.Gain_b = VCM20_P.Gain_Gain_j * VCM20_B.Saturation1;

      /* Switch: '<S292>/Switch1' */
      VCM20_B.Switch1_p = VCM20_B.Gain_b;
    }

    /* End of Switch: '<S292>/Switch1' */
  }

  /* RelationalOperator: '<S287>/UpperRelop' */
  VCM20_B.UpperRelop_n2 = (VCM20_B.Switch2_i < VCM20_B.Switch1_p);

  /* Switch: '<S287>/Switch' */
  if (VCM20_B.UpperRelop_n2) {
    /* Switch: '<S287>/Switch' */
    VCM20_B.Switch_b4 = VCM20_B.Switch1_p;
  } else {
    /* Switch: '<S287>/Switch' */
    VCM20_B.Switch_b4 = VCM20_B.Switch2_i;
  }

  /* End of Switch: '<S287>/Switch' */

  /* Switch: '<S287>/Switch2' */
  if (VCM20_B.LowerRelop1_mx) {
    /* Switch: '<S287>/Switch2' */
    VCM20_B.Switch2_lc = VCM20_B.torque_p;
  } else {
    /* Switch: '<S287>/Switch2' */
    VCM20_B.Switch2_lc = VCM20_B.Switch_b4;
  }

  /* End of Switch: '<S287>/Switch2' */

  /* Saturate: '<S185>/Saturation' */
  u0_0 = VCM20_B.Switch2_lc;
  u1 = VCM20_P.Saturation_LowerSat_b;
  u2 = VCM20_P.Saturation_UpperSat_k;
  if (u0_0 > u2) {
    /* Saturate: '<S185>/Saturation' */
    VCM20_B.Saturation_aa = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S185>/Saturation' */
    VCM20_B.Saturation_aa = u1;
  } else {
    /* Saturate: '<S185>/Saturation' */
    VCM20_B.Saturation_aa = u0_0;
  }

  /* End of Saturate: '<S185>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S185>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_g;
    u1 = VCM20_P.Saturation1_LowerSat_o;
    u2 = VCM20_P.Saturation1_UpperSat_g;
    if (u0_0 > u2) {
      /* Saturate: '<S185>/Saturation1' */
      VCM20_B.Saturation1_e = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S185>/Saturation1' */
      VCM20_B.Saturation1_e = u1;
    } else {
      /* Saturate: '<S185>/Saturation1' */
      VCM20_B.Saturation1_e = u0_0;
    }

    /* End of Saturate: '<S185>/Saturation1' */

    /* DiscreteIntegrator: '<S293>/Discrete-Time Integrator' */
    VCM20_B.DiscreteTimeIntegrator_a = VCM20_DW.DiscreteTimeIntegrator_DSTATE_o;
  }

  /* Product: '<S185>/Product' */
  VCM20_B.FLMotorPower = VCM20_B.Saturation_aa * VCM20_B.Saturation1_e;

  /* Gain: '<S185>/Gain' */
  VCM20_B.FLMotorPower_W = VCM20_P.Gain_Gain_g * VCM20_B.FLMotorPower;

  /* Lookup_n-D: '<S185>/2-D Lookup Table' incorporates:
   *  Saturate: '<S185>/Saturation'
   *  Saturate: '<S185>/Saturation1'
   */
  VCM20_B.FLEfficency = look2_binlxpw(VCM20_B.Saturation_aa,
    VCM20_B.Saturation1_e, VCM20_P.uDLookupTable_bp01Data_i,
    VCM20_P.uDLookupTable_bp02Data_a, VCM20_P.uDLookupTable_tableData_o,
    VCM20_P.uDLookupTable_maxIndex_k, 11U);

  /* Product: '<S185>/Divide' */
  VCM20_B.FLPowerConsumption = VCM20_B.FLMotorPower_W / VCM20_B.FLEfficency;

  /* Sum: '<S293>/Add3' */
  VCM20_B.Add3_f = VCM20_B.Switch_p - VCM20_B.SFunction1_o9;

  /* Gain: '<S293>/P gain' */
  VCM20_B.Pgain_m = VCM20_P.Pgain_Gain_o * VCM20_B.Add3_f;

  /* Sum: '<S293>/Add2' */
  VCM20_B.FB_h = VCM20_B.Pgain_m + VCM20_B.DiscreteTimeIntegrator_a;

  /* Sum: '<S293>/Add' */
  VCM20_B.torque_i = VCM20_B.Add1_jx + VCM20_B.FB_h;

  /* RelationalOperator: '<S315>/LowerRelop1' */
  VCM20_B.LowerRelop1_h = (VCM20_B.torque_i > VCM20_B.Saturation);

  /* RelationalOperator: '<S315>/UpperRelop' incorporates:
   *  Constant: '<S293>/Constant'
   */
  VCM20_B.UpperRelop_ba = (VCM20_B.torque_i < VCM20_P.Constant_Value_fe);

  /* Switch: '<S315>/Switch' */
  if (VCM20_B.UpperRelop_ba) {
    /* Switch: '<S315>/Switch' incorporates:
     *  Constant: '<S293>/Constant'
     */
    VCM20_B.Switch_eh = VCM20_P.Constant_Value_fe;
  } else {
    /* Switch: '<S315>/Switch' */
    VCM20_B.Switch_eh = VCM20_B.torque_i;
  }

  /* End of Switch: '<S315>/Switch' */

  /* Switch: '<S315>/Switch2' */
  if (VCM20_B.LowerRelop1_h) {
    /* Switch: '<S315>/Switch2' */
    VCM20_B.Switch2_j3 = VCM20_B.Saturation;
  } else {
    /* Switch: '<S315>/Switch2' */
    VCM20_B.Switch2_j3 = VCM20_B.Switch_eh;
  }

  /* End of Switch: '<S315>/Switch2' */

  /* Switch: '<S293>/Switch' incorporates:
   *  Constant: '<Root>/TCOn'
   */
  if (VCM20_P.TCOn_Value > VCM20_P.Switch_Threshold_i) {
    /* Switch: '<S293>/Switch' */
    VCM20_B.torque_fc = VCM20_B.Switch2_j3;
  } else {
    /* Switch: '<S293>/Switch' */
    VCM20_B.torque_fc = VCM20_B.Saturation;
  }

  /* End of Switch: '<S293>/Switch' */

  /* RelationalOperator: '<S288>/LowerRelop1' */
  VCM20_B.LowerRelop1_gs = (VCM20_B.Switch2_d > VCM20_B.torque_fc);
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Sum: '<S293>/Add1' */
    VCM20_B.Add1_k = VCM20_B.SFunction1_o9 - VCM20_B.Switch2_o;

    /* RelationalOperator: '<S316>/LowerRelop1' incorporates:
     *  Constant: '<S293>/Constant'
     */
    VCM20_B.LowerRelop1_hq = (VCM20_B.Add1_k > VCM20_P.Constant_Value_fe);

    /* Gain: '<S293>/Gain1' */
    VCM20_B.Gain1_lu = VCM20_P.Gain1_Gain_d4 * VCM20_B.Saturation1;

    /* RelationalOperator: '<S316>/UpperRelop' */
    VCM20_B.UpperRelop_m4 = (VCM20_B.Add1_k < VCM20_B.Gain1_lu);

    /* Switch: '<S316>/Switch' */
    if (VCM20_B.UpperRelop_m4) {
      /* Switch: '<S316>/Switch' */
      VCM20_B.Switch_ah = VCM20_B.Gain1_lu;
    } else {
      /* Switch: '<S316>/Switch' */
      VCM20_B.Switch_ah = VCM20_B.Add1_k;
    }

    /* End of Switch: '<S316>/Switch' */

    /* Switch: '<S316>/Switch2' */
    if (VCM20_B.LowerRelop1_hq) {
      /* Switch: '<S316>/Switch2' incorporates:
       *  Constant: '<S293>/Constant'
       */
      VCM20_B.Switch2_h3 = VCM20_P.Constant_Value_fe;
    } else {
      /* Switch: '<S316>/Switch2' */
      VCM20_B.Switch2_h3 = VCM20_B.Switch_ah;
    }

    /* End of Switch: '<S316>/Switch2' */

    /* Switch: '<S293>/Switch1' incorporates:
     *  Constant: '<Root>/ABSOn'
     */
    if (VCM20_P.ABSOn_Value > VCM20_P.Switch1_Threshold_fn) {
      /* Switch: '<S293>/Switch1' */
      VCM20_B.Switch1_k = VCM20_B.Switch2_h3;
    } else {
      /* Gain: '<S293>/Gain' */
      VCM20_B.Gain_a = VCM20_P.Gain_Gain_a * VCM20_B.Saturation1;

      /* Switch: '<S293>/Switch1' */
      VCM20_B.Switch1_k = VCM20_B.Gain_a;
    }

    /* End of Switch: '<S293>/Switch1' */
  }

  /* RelationalOperator: '<S288>/UpperRelop' */
  VCM20_B.UpperRelop_gj = (VCM20_B.Switch2_d < VCM20_B.Switch1_k);

  /* Switch: '<S288>/Switch' */
  if (VCM20_B.UpperRelop_gj) {
    /* Switch: '<S288>/Switch' */
    VCM20_B.Switch_ma = VCM20_B.Switch1_k;
  } else {
    /* Switch: '<S288>/Switch' */
    VCM20_B.Switch_ma = VCM20_B.Switch2_d;
  }

  /* End of Switch: '<S288>/Switch' */

  /* Switch: '<S288>/Switch2' */
  if (VCM20_B.LowerRelop1_gs) {
    /* Switch: '<S288>/Switch2' */
    VCM20_B.Switch2_hr = VCM20_B.torque_fc;
  } else {
    /* Switch: '<S288>/Switch2' */
    VCM20_B.Switch2_hr = VCM20_B.Switch_ma;
  }

  /* End of Switch: '<S288>/Switch2' */

  /* Saturate: '<S188>/Saturation' */
  u0_0 = VCM20_B.Switch2_hr;
  u1 = VCM20_P.Saturation_LowerSat_l;
  u2 = VCM20_P.Saturation_UpperSat_j;
  if (u0_0 > u2) {
    /* Saturate: '<S188>/Saturation' */
    VCM20_B.Saturation_p = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S188>/Saturation' */
    VCM20_B.Saturation_p = u1;
  } else {
    /* Saturate: '<S188>/Saturation' */
    VCM20_B.Saturation_p = u0_0;
  }

  /* End of Saturate: '<S188>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S188>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_j;
    u1 = VCM20_P.Saturation1_LowerSat_ov;
    u2 = VCM20_P.Saturation1_UpperSat_c;
    if (u0_0 > u2) {
      /* Saturate: '<S188>/Saturation1' */
      VCM20_B.Saturation1_h = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S188>/Saturation1' */
      VCM20_B.Saturation1_h = u1;
    } else {
      /* Saturate: '<S188>/Saturation1' */
      VCM20_B.Saturation1_h = u0_0;
    }

    /* End of Saturate: '<S188>/Saturation1' */
  }

  /* Product: '<S188>/Product' */
  VCM20_B.RRMotorPower = VCM20_B.Saturation_p * VCM20_B.Saturation1_h;

  /* Gain: '<S188>/Gain' */
  VCM20_B.RRMotorPower_W = VCM20_P.Gain_Gain_p * VCM20_B.RRMotorPower;

  /* Lookup_n-D: '<S188>/2-D Lookup Table' incorporates:
   *  Saturate: '<S188>/Saturation'
   *  Saturate: '<S188>/Saturation1'
   */
  VCM20_B.RREfficency = look2_binlxpw(VCM20_B.Saturation_p,
    VCM20_B.Saturation1_h, VCM20_P.uDLookupTable_bp01Data_l,
    VCM20_P.uDLookupTable_bp02Data_e, VCM20_P.uDLookupTable_tableData_h,
    VCM20_P.uDLookupTable_maxIndex_p, 11U);

  /* Product: '<S188>/Divide' */
  VCM20_B.RRPowerConsumption = VCM20_B.RRMotorPower_W / VCM20_B.RREfficency;

  /* Saturate: '<S187>/Saturation' */
  u0_0 = VCM20_B.Switch2_ot;
  u1 = VCM20_P.Saturation_LowerSat_cm;
  u2 = VCM20_P.Saturation_UpperSat_o;
  if (u0_0 > u2) {
    /* Saturate: '<S187>/Saturation' */
    VCM20_B.Saturation_c = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S187>/Saturation' */
    VCM20_B.Saturation_c = u1;
  } else {
    /* Saturate: '<S187>/Saturation' */
    VCM20_B.Saturation_c = u0_0;
  }

  /* End of Saturate: '<S187>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S187>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9;
    u1 = VCM20_P.Saturation1_LowerSat_e;
    u2 = VCM20_P.Saturation1_UpperSat_m2;
    if (u0_0 > u2) {
      /* Saturate: '<S187>/Saturation1' */
      VCM20_B.Saturation1_a = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S187>/Saturation1' */
      VCM20_B.Saturation1_a = u1;
    } else {
      /* Saturate: '<S187>/Saturation1' */
      VCM20_B.Saturation1_a = u0_0;
    }

    /* End of Saturate: '<S187>/Saturation1' */
  }

  /* Product: '<S187>/Product' */
  VCM20_B.RLMotorPower = VCM20_B.Saturation_c * VCM20_B.Saturation1_a;

  /* Gain: '<S187>/Gain' */
  VCM20_B.RLMotorPower_W = VCM20_P.Gain_Gain_hv * VCM20_B.RLMotorPower;

  /* Lookup_n-D: '<S187>/2-D Lookup Table' incorporates:
   *  Saturate: '<S187>/Saturation'
   *  Saturate: '<S187>/Saturation1'
   */
  VCM20_B.RLEfficency = look2_binlxpw(VCM20_B.Saturation_c,
    VCM20_B.Saturation1_a, VCM20_P.uDLookupTable_bp01Data_iu,
    VCM20_P.uDLookupTable_bp02Data_l, VCM20_P.uDLookupTable_tableData_e,
    VCM20_P.uDLookupTable_maxIndex_h, 11U);

  /* Product: '<S187>/Divide' */
  VCM20_B.RLPowerConsumption = VCM20_B.RLMotorPower_W / VCM20_B.RLEfficency;

  /* Sum: '<S178>/Sum' */
  VCM20_B.Sum_j = ((VCM20_B.FRPowerConsumption + VCM20_B.FLPowerConsumption) +
                   VCM20_B.RRPowerConsumption) + VCM20_B.RLPowerConsumption;

  /* Saturate: '<S178>/Saturation' */
  u0_0 = VCM20_B.Sum_j;
  u1 = VCM20_P.Saturation_LowerSat_e;
  u2 = VCM20_P.Saturation_UpperSat_g;
  if (u0_0 > u2) {
    /* Saturate: '<S178>/Saturation' */
    VCM20_B.Saturation_n = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S178>/Saturation' */
    VCM20_B.Saturation_n = u1;
  } else {
    /* Saturate: '<S178>/Saturation' */
    VCM20_B.Saturation_n = u0_0;
  }

  /* End of Saturate: '<S178>/Saturation' */

  /* Product: '<S178>/Divide3' */
  VCM20_B.Divide3_p = VCM20_B.Switch2_ot / VCM20_B.Saturation_n;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Gain: '<S170>/Eff_INV' incorporates:
     *  Constant: '<S170>/Max Power'
     */
    VCM20_B.Eff_INV = VCM20_P.Eff_INV_Gain * VCM20_P.MaxPower_Value;
  }

  /* Product: '<S178>/Product3' */
  VCM20_B.torque_ba = VCM20_B.Divide3_p * VCM20_B.Eff_INV;

  /* Saturate: '<S178>/Saturation4' */
  u0_0 = VCM20_B.torque_ba;
  u1 = VCM20_P.Saturation4_LowerSat;
  u2 = VCM20_P.Saturation4_UpperSat;
  if (u0_0 > u2) {
    /* Saturate: '<S178>/Saturation4' */
    VCM20_B.torque_gi = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S178>/Saturation4' */
    VCM20_B.torque_gi = u1;
  } else {
    /* Saturate: '<S178>/Saturation4' */
    VCM20_B.torque_gi = u0_0;
  }

  /* End of Saturate: '<S178>/Saturation4' */

  /* Product: '<S178>/Divide' */
  VCM20_B.Divide_iz = VCM20_B.Switch2_g / VCM20_B.Saturation_n;

  /* Product: '<S178>/Product' */
  VCM20_B.torque_j = VCM20_B.Divide_iz * VCM20_B.Eff_INV;

  /* Saturate: '<S178>/Saturation1' */
  u0_0 = VCM20_B.torque_j;
  u1 = VCM20_P.Saturation1_LowerSat_gh;
  u2 = VCM20_P.Saturation1_UpperSat_p;
  if (u0_0 > u2) {
    /* Saturate: '<S178>/Saturation1' */
    VCM20_B.torque_po = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S178>/Saturation1' */
    VCM20_B.torque_po = u1;
  } else {
    /* Saturate: '<S178>/Saturation1' */
    VCM20_B.torque_po = u0_0;
  }

  /* End of Saturate: '<S178>/Saturation1' */

  /* Saturate: '<S190>/Saturation' */
  u0_0 = VCM20_B.torque_po;
  u1 = VCM20_P.Saturation_LowerSat_m;
  u2 = VCM20_P.Saturation_UpperSat_eq;
  if (u0_0 > u2) {
    /* Saturate: '<S190>/Saturation' */
    VCM20_B.Saturation_d2 = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S190>/Saturation' */
    VCM20_B.Saturation_d2 = u1;
  } else {
    /* Saturate: '<S190>/Saturation' */
    VCM20_B.Saturation_d2 = u0_0;
  }

  /* End of Saturate: '<S190>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S190>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_n;
    u1 = VCM20_P.Saturation1_LowerSat_gn;
    u2 = VCM20_P.Saturation1_UpperSat_b;
    if (u0_0 > u2) {
      /* Saturate: '<S190>/Saturation1' */
      VCM20_B.Saturation1_hi = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S190>/Saturation1' */
      VCM20_B.Saturation1_hi = u1;
    } else {
      /* Saturate: '<S190>/Saturation1' */
      VCM20_B.Saturation1_hi = u0_0;
    }

    /* End of Saturate: '<S190>/Saturation1' */
  }

  /* Product: '<S190>/Product' */
  VCM20_B.FRMotorPower_l = VCM20_B.Saturation_d2 * VCM20_B.Saturation1_hi;

  /* Gain: '<S190>/Gain' */
  VCM20_B.FRMotorPower_W_n = VCM20_P.Gain_Gain_bm * VCM20_B.FRMotorPower_l;

  /* Lookup_n-D: '<S190>/2-D Lookup Table' incorporates:
   *  Saturate: '<S190>/Saturation'
   *  Saturate: '<S190>/Saturation1'
   */
  VCM20_B.FREfficency_j = look2_binlxpw(VCM20_B.Saturation_d2,
    VCM20_B.Saturation1_hi, VCM20_P.uDLookupTable_bp01Data_id,
    VCM20_P.uDLookupTable_bp02Data_n, VCM20_P.uDLookupTable_tableData_hp,
    VCM20_P.uDLookupTable_maxIndex_n, 11U);

  /* Product: '<S190>/Divide' */
  VCM20_B.FRPowerConsumption_p = VCM20_B.FRMotorPower_W_n /
    VCM20_B.FREfficency_j;

  /* Product: '<S178>/Divide1' */
  VCM20_B.Divide1_d = VCM20_B.Switch2_lc / VCM20_B.Saturation_n;

  /* Product: '<S178>/Product1' */
  VCM20_B.torque_jv = VCM20_B.Divide1_d * VCM20_B.Eff_INV;

  /* Saturate: '<S178>/Saturation2' */
  u0_0 = VCM20_B.torque_jv;
  u1 = VCM20_P.Saturation2_LowerSat_h;
  u2 = VCM20_P.Saturation2_UpperSat_o;
  if (u0_0 > u2) {
    /* Saturate: '<S178>/Saturation2' */
    VCM20_B.torque_m = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S178>/Saturation2' */
    VCM20_B.torque_m = u1;
  } else {
    /* Saturate: '<S178>/Saturation2' */
    VCM20_B.torque_m = u0_0;
  }

  /* End of Saturate: '<S178>/Saturation2' */

  /* Saturate: '<S189>/Saturation' */
  u0_0 = VCM20_B.torque_m;
  u1 = VCM20_P.Saturation_LowerSat_ir;
  u2 = VCM20_P.Saturation_UpperSat_h;
  if (u0_0 > u2) {
    /* Saturate: '<S189>/Saturation' */
    VCM20_B.Saturation_l = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S189>/Saturation' */
    VCM20_B.Saturation_l = u1;
  } else {
    /* Saturate: '<S189>/Saturation' */
    VCM20_B.Saturation_l = u0_0;
  }

  /* End of Saturate: '<S189>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S189>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_g;
    u1 = VCM20_P.Saturation1_LowerSat_d;
    u2 = VCM20_P.Saturation1_UpperSat_l;
    if (u0_0 > u2) {
      /* Saturate: '<S189>/Saturation1' */
      VCM20_B.Saturation1_ae = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S189>/Saturation1' */
      VCM20_B.Saturation1_ae = u1;
    } else {
      /* Saturate: '<S189>/Saturation1' */
      VCM20_B.Saturation1_ae = u0_0;
    }

    /* End of Saturate: '<S189>/Saturation1' */
  }

  /* Product: '<S189>/Product' */
  VCM20_B.FLMotorPower_g = VCM20_B.Saturation_l * VCM20_B.Saturation1_ae;

  /* Gain: '<S189>/Gain' */
  VCM20_B.FLMotorPower_W_b = VCM20_P.Gain_Gain_gh * VCM20_B.FLMotorPower_g;

  /* Lookup_n-D: '<S189>/2-D Lookup Table' incorporates:
   *  Saturate: '<S189>/Saturation'
   *  Saturate: '<S189>/Saturation1'
   */
  VCM20_B.FLEfficency_p = look2_binlxpw(VCM20_B.Saturation_l,
    VCM20_B.Saturation1_ae, VCM20_P.uDLookupTable_bp01Data_h,
    VCM20_P.uDLookupTable_bp02Data_j, VCM20_P.uDLookupTable_tableData_d,
    VCM20_P.uDLookupTable_maxIndex_pw, 11U);

  /* Product: '<S189>/Divide' */
  VCM20_B.FLPowerConsumption_d = VCM20_B.FLMotorPower_W_b /
    VCM20_B.FLEfficency_p;

  /* Product: '<S178>/Divide2' */
  VCM20_B.Divide2_a = VCM20_B.Switch2_hr / VCM20_B.Saturation_n;

  /* Product: '<S178>/Product2' */
  VCM20_B.torque_iy = VCM20_B.Divide2_a * VCM20_B.Eff_INV;

  /* Saturate: '<S178>/Saturation3' */
  u0_0 = VCM20_B.torque_iy;
  u1 = VCM20_P.Saturation3_LowerSat_c;
  u2 = VCM20_P.Saturation3_UpperSat_k;
  if (u0_0 > u2) {
    /* Saturate: '<S178>/Saturation3' */
    VCM20_B.torque_bu = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S178>/Saturation3' */
    VCM20_B.torque_bu = u1;
  } else {
    /* Saturate: '<S178>/Saturation3' */
    VCM20_B.torque_bu = u0_0;
  }

  /* End of Saturate: '<S178>/Saturation3' */

  /* Saturate: '<S192>/Saturation' */
  u0_0 = VCM20_B.torque_bu;
  u1 = VCM20_P.Saturation_LowerSat_a;
  u2 = VCM20_P.Saturation_UpperSat_e2;
  if (u0_0 > u2) {
    /* Saturate: '<S192>/Saturation' */
    VCM20_B.Saturation_go = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S192>/Saturation' */
    VCM20_B.Saturation_go = u1;
  } else {
    /* Saturate: '<S192>/Saturation' */
    VCM20_B.Saturation_go = u0_0;
  }

  /* End of Saturate: '<S192>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S192>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_j;
    u1 = VCM20_P.Saturation1_LowerSat_h;
    u2 = VCM20_P.Saturation1_UpperSat_bu;
    if (u0_0 > u2) {
      /* Saturate: '<S192>/Saturation1' */
      VCM20_B.Saturation1_ia = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S192>/Saturation1' */
      VCM20_B.Saturation1_ia = u1;
    } else {
      /* Saturate: '<S192>/Saturation1' */
      VCM20_B.Saturation1_ia = u0_0;
    }

    /* End of Saturate: '<S192>/Saturation1' */
  }

  /* Product: '<S192>/Product' */
  VCM20_B.RRMotorPower_b = VCM20_B.Saturation_go * VCM20_B.Saturation1_ia;

  /* Gain: '<S192>/Gain' */
  VCM20_B.RRMotorPower_W_m = VCM20_P.Gain_Gain_m * VCM20_B.RRMotorPower_b;

  /* Lookup_n-D: '<S192>/2-D Lookup Table' incorporates:
   *  Saturate: '<S192>/Saturation'
   *  Saturate: '<S192>/Saturation1'
   */
  VCM20_B.RREfficency_d = look2_binlxpw(VCM20_B.Saturation_go,
    VCM20_B.Saturation1_ia, VCM20_P.uDLookupTable_bp01Data_lj,
    VCM20_P.uDLookupTable_bp02Data_p, VCM20_P.uDLookupTable_tableData_g,
    VCM20_P.uDLookupTable_maxIndex_j, 11U);

  /* Product: '<S192>/Divide' */
  VCM20_B.RRPowerConsumption_k = VCM20_B.RRMotorPower_W_m /
    VCM20_B.RREfficency_d;

  /* Saturate: '<S191>/Saturation' */
  u0_0 = VCM20_B.torque_gi;
  u1 = VCM20_P.Saturation_LowerSat_j;
  u2 = VCM20_P.Saturation_UpperSat_b;
  if (u0_0 > u2) {
    /* Saturate: '<S191>/Saturation' */
    VCM20_B.Saturation_go2 = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S191>/Saturation' */
    VCM20_B.Saturation_go2 = u1;
  } else {
    /* Saturate: '<S191>/Saturation' */
    VCM20_B.Saturation_go2 = u0_0;
  }

  /* End of Saturate: '<S191>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S191>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9;
    u1 = VCM20_P.Saturation1_LowerSat_o0;
    u2 = VCM20_P.Saturation1_UpperSat_o;
    if (u0_0 > u2) {
      /* Saturate: '<S191>/Saturation1' */
      VCM20_B.Saturation1_p = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S191>/Saturation1' */
      VCM20_B.Saturation1_p = u1;
    } else {
      /* Saturate: '<S191>/Saturation1' */
      VCM20_B.Saturation1_p = u0_0;
    }

    /* End of Saturate: '<S191>/Saturation1' */
  }

  /* Product: '<S191>/Product' */
  VCM20_B.RLMotorPower_h = VCM20_B.Saturation_go2 * VCM20_B.Saturation1_p;

  /* Gain: '<S191>/Gain' */
  VCM20_B.RLMotorPower_W_a = VCM20_P.Gain_Gain_a1 * VCM20_B.RLMotorPower_h;

  /* Lookup_n-D: '<S191>/2-D Lookup Table' incorporates:
   *  Saturate: '<S191>/Saturation'
   *  Saturate: '<S191>/Saturation1'
   */
  VCM20_B.RLEfficency_j = look2_binlxpw(VCM20_B.Saturation_go2,
    VCM20_B.Saturation1_p, VCM20_P.uDLookupTable_bp01Data_d,
    VCM20_P.uDLookupTable_bp02Data_b, VCM20_P.uDLookupTable_tableData_j,
    VCM20_P.uDLookupTable_maxIndex_pr, 11U);

  /* Product: '<S191>/Divide' */
  VCM20_B.RLPowerConsumption_e = VCM20_B.RLMotorPower_W_a /
    VCM20_B.RLEfficency_j;

  /* Sum: '<S179>/Sum' */
  VCM20_B.Sum_a = ((VCM20_B.FRPowerConsumption_p + VCM20_B.FLPowerConsumption_d)
                   + VCM20_B.RRPowerConsumption_k) +
    VCM20_B.RLPowerConsumption_e;

  /* Saturate: '<S179>/Saturation' */
  u0_0 = VCM20_B.Sum_a;
  u1 = VCM20_P.Saturation_LowerSat_px;
  u2 = VCM20_P.Saturation_UpperSat_fe;
  if (u0_0 > u2) {
    /* Saturate: '<S179>/Saturation' */
    VCM20_B.Saturation_gs = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S179>/Saturation' */
    VCM20_B.Saturation_gs = u1;
  } else {
    /* Saturate: '<S179>/Saturation' */
    VCM20_B.Saturation_gs = u0_0;
  }

  /* End of Saturate: '<S179>/Saturation' */

  /* Product: '<S179>/Divide3' */
  VCM20_B.Divide3_g = VCM20_B.torque_gi / VCM20_B.Saturation_gs;

  /* Product: '<S179>/Product3' */
  VCM20_B.torque_a = VCM20_B.Divide3_g * VCM20_B.Eff_INV;

  /* Saturate: '<S179>/Saturation4' */
  u0_0 = VCM20_B.torque_a;
  u1 = VCM20_P.Saturation4_LowerSat_m;
  u2 = VCM20_P.Saturation4_UpperSat_j;
  if (u0_0 > u2) {
    /* Saturate: '<S179>/Saturation4' */
    VCM20_B.torque_o = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S179>/Saturation4' */
    VCM20_B.torque_o = u1;
  } else {
    /* Saturate: '<S179>/Saturation4' */
    VCM20_B.torque_o = u0_0;
  }

  /* End of Saturate: '<S179>/Saturation4' */

  /* Product: '<S179>/Divide' */
  VCM20_B.Divide_ir = VCM20_B.torque_po / VCM20_B.Saturation_gs;

  /* Product: '<S179>/Product' */
  VCM20_B.torque_az = VCM20_B.Divide_ir * VCM20_B.Eff_INV;

  /* Saturate: '<S179>/Saturation1' */
  u0_0 = VCM20_B.torque_az;
  u1 = VCM20_P.Saturation1_LowerSat_d1;
  u2 = VCM20_P.Saturation1_UpperSat_h;
  if (u0_0 > u2) {
    /* Saturate: '<S179>/Saturation1' */
    VCM20_B.torque_l = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S179>/Saturation1' */
    VCM20_B.torque_l = u1;
  } else {
    /* Saturate: '<S179>/Saturation1' */
    VCM20_B.torque_l = u0_0;
  }

  /* End of Saturate: '<S179>/Saturation1' */

  /* Saturate: '<S194>/Saturation' */
  u0_0 = VCM20_B.torque_l;
  u1 = VCM20_P.Saturation_LowerSat_d;
  u2 = VCM20_P.Saturation_UpperSat_fd;
  if (u0_0 > u2) {
    /* Saturate: '<S194>/Saturation' */
    VCM20_B.Saturation_e = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S194>/Saturation' */
    VCM20_B.Saturation_e = u1;
  } else {
    /* Saturate: '<S194>/Saturation' */
    VCM20_B.Saturation_e = u0_0;
  }

  /* End of Saturate: '<S194>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S194>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_n;
    u1 = VCM20_P.Saturation1_LowerSat_n;
    u2 = VCM20_P.Saturation1_UpperSat_e;
    if (u0_0 > u2) {
      /* Saturate: '<S194>/Saturation1' */
      VCM20_B.Saturation1_hf = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S194>/Saturation1' */
      VCM20_B.Saturation1_hf = u1;
    } else {
      /* Saturate: '<S194>/Saturation1' */
      VCM20_B.Saturation1_hf = u0_0;
    }

    /* End of Saturate: '<S194>/Saturation1' */
  }

  /* Product: '<S194>/Product' */
  VCM20_B.FRMotorPower_n = VCM20_B.Saturation_e * VCM20_B.Saturation1_hf;

  /* Gain: '<S194>/Gain' */
  VCM20_B.FRMotorPower_W_o = VCM20_P.Gain_Gain_br * VCM20_B.FRMotorPower_n;

  /* Lookup_n-D: '<S194>/2-D Lookup Table' incorporates:
   *  Saturate: '<S194>/Saturation'
   *  Saturate: '<S194>/Saturation1'
   */
  VCM20_B.FREfficency_p = look2_binlxpw(VCM20_B.Saturation_e,
    VCM20_B.Saturation1_hf, VCM20_P.uDLookupTable_bp01Data_i4,
    VCM20_P.uDLookupTable_bp02Data_f, VCM20_P.uDLookupTable_tableData_p,
    VCM20_P.uDLookupTable_maxIndex_a, 11U);

  /* Product: '<S194>/Divide' */
  VCM20_B.FRPowerConsumption_i = VCM20_B.FRMotorPower_W_o /
    VCM20_B.FREfficency_p;

  /* Product: '<S179>/Divide1' */
  VCM20_B.Divide1_i = VCM20_B.torque_m / VCM20_B.Saturation_gs;

  /* Product: '<S179>/Product1' */
  VCM20_B.torque_lg = VCM20_B.Divide1_i * VCM20_B.Eff_INV;

  /* Saturate: '<S179>/Saturation2' */
  u0_0 = VCM20_B.torque_lg;
  u1 = VCM20_P.Saturation2_LowerSat_i;
  u2 = VCM20_P.Saturation2_UpperSat_p;
  if (u0_0 > u2) {
    /* Saturate: '<S179>/Saturation2' */
    VCM20_B.torque_e = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S179>/Saturation2' */
    VCM20_B.torque_e = u1;
  } else {
    /* Saturate: '<S179>/Saturation2' */
    VCM20_B.torque_e = u0_0;
  }

  /* End of Saturate: '<S179>/Saturation2' */

  /* Saturate: '<S193>/Saturation' */
  u0_0 = VCM20_B.torque_e;
  u1 = VCM20_P.Saturation_LowerSat_ns;
  u2 = VCM20_P.Saturation_UpperSat_fg;
  if (u0_0 > u2) {
    /* Saturate: '<S193>/Saturation' */
    VCM20_B.Saturation_hc = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S193>/Saturation' */
    VCM20_B.Saturation_hc = u1;
  } else {
    /* Saturate: '<S193>/Saturation' */
    VCM20_B.Saturation_hc = u0_0;
  }

  /* End of Saturate: '<S193>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S193>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_g;
    u1 = VCM20_P.Saturation1_LowerSat_g4;
    u2 = VCM20_P.Saturation1_UpperSat_de;
    if (u0_0 > u2) {
      /* Saturate: '<S193>/Saturation1' */
      VCM20_B.Saturation1_h5 = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S193>/Saturation1' */
      VCM20_B.Saturation1_h5 = u1;
    } else {
      /* Saturate: '<S193>/Saturation1' */
      VCM20_B.Saturation1_h5 = u0_0;
    }

    /* End of Saturate: '<S193>/Saturation1' */
  }

  /* Product: '<S193>/Product' */
  VCM20_B.FLMotorPower_o = VCM20_B.Saturation_hc * VCM20_B.Saturation1_h5;

  /* Gain: '<S193>/Gain' */
  VCM20_B.FLMotorPower_W_c = VCM20_P.Gain_Gain_fy * VCM20_B.FLMotorPower_o;

  /* Lookup_n-D: '<S193>/2-D Lookup Table' incorporates:
   *  Saturate: '<S193>/Saturation'
   *  Saturate: '<S193>/Saturation1'
   */
  VCM20_B.FLEfficency_o = look2_binlxpw(VCM20_B.Saturation_hc,
    VCM20_B.Saturation1_h5, VCM20_P.uDLookupTable_bp01Data_b,
    VCM20_P.uDLookupTable_bp02Data_be, VCM20_P.uDLookupTable_tableData_l,
    VCM20_P.uDLookupTable_maxIndex_nc, 11U);

  /* Product: '<S193>/Divide' */
  VCM20_B.FLPowerConsumption_i = VCM20_B.FLMotorPower_W_c /
    VCM20_B.FLEfficency_o;

  /* Product: '<S179>/Divide2' */
  VCM20_B.Divide2_e = VCM20_B.torque_bu / VCM20_B.Saturation_gs;

  /* Product: '<S179>/Product2' */
  VCM20_B.torque_f4 = VCM20_B.Divide2_e * VCM20_B.Eff_INV;

  /* Saturate: '<S179>/Saturation3' */
  u0_0 = VCM20_B.torque_f4;
  u1 = VCM20_P.Saturation3_LowerSat_p;
  u2 = VCM20_P.Saturation3_UpperSat_k0;
  if (u0_0 > u2) {
    /* Saturate: '<S179>/Saturation3' */
    VCM20_B.torque_bc = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S179>/Saturation3' */
    VCM20_B.torque_bc = u1;
  } else {
    /* Saturate: '<S179>/Saturation3' */
    VCM20_B.torque_bc = u0_0;
  }

  /* End of Saturate: '<S179>/Saturation3' */

  /* Saturate: '<S196>/Saturation' */
  u0_0 = VCM20_B.torque_bc;
  u1 = VCM20_P.Saturation_LowerSat_o;
  u2 = VCM20_P.Saturation_UpperSat_m;
  if (u0_0 > u2) {
    /* Saturate: '<S196>/Saturation' */
    VCM20_B.Saturation_gy = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S196>/Saturation' */
    VCM20_B.Saturation_gy = u1;
  } else {
    /* Saturate: '<S196>/Saturation' */
    VCM20_B.Saturation_gy = u0_0;
  }

  /* End of Saturate: '<S196>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S196>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_j;
    u1 = VCM20_P.Saturation1_LowerSat_nu;
    u2 = VCM20_P.Saturation1_UpperSat_k;
    if (u0_0 > u2) {
      /* Saturate: '<S196>/Saturation1' */
      VCM20_B.Saturation1_c = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S196>/Saturation1' */
      VCM20_B.Saturation1_c = u1;
    } else {
      /* Saturate: '<S196>/Saturation1' */
      VCM20_B.Saturation1_c = u0_0;
    }

    /* End of Saturate: '<S196>/Saturation1' */
  }

  /* Product: '<S196>/Product' */
  VCM20_B.RRMotorPower_j = VCM20_B.Saturation_gy * VCM20_B.Saturation1_c;

  /* Gain: '<S196>/Gain' */
  VCM20_B.RRMotorPower_W_h = VCM20_P.Gain_Gain_c * VCM20_B.RRMotorPower_j;

  /* Lookup_n-D: '<S196>/2-D Lookup Table' incorporates:
   *  Saturate: '<S196>/Saturation'
   *  Saturate: '<S196>/Saturation1'
   */
  VCM20_B.RREfficency_n = look2_binlxpw(VCM20_B.Saturation_gy,
    VCM20_B.Saturation1_c, VCM20_P.uDLookupTable_bp01Data_m,
    VCM20_P.uDLookupTable_bp02Data_jg, VCM20_P.uDLookupTable_tableData_ji,
    VCM20_P.uDLookupTable_maxIndex_pq, 11U);

  /* Product: '<S196>/Divide' */
  VCM20_B.RRPowerConsumption_b = VCM20_B.RRMotorPower_W_h /
    VCM20_B.RREfficency_n;

  /* Saturate: '<S195>/Saturation' */
  u0_0 = VCM20_B.torque_o;
  u1 = VCM20_P.Saturation_LowerSat_lr;
  u2 = VCM20_P.Saturation_UpperSat_j2;
  if (u0_0 > u2) {
    /* Saturate: '<S195>/Saturation' */
    VCM20_B.Saturation_gc = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S195>/Saturation' */
    VCM20_B.Saturation_gc = u1;
  } else {
    /* Saturate: '<S195>/Saturation' */
    VCM20_B.Saturation_gc = u0_0;
  }

  /* End of Saturate: '<S195>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S195>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9;
    u1 = VCM20_P.Saturation1_LowerSat_j;
    u2 = VCM20_P.Saturation1_UpperSat_hj;
    if (u0_0 > u2) {
      /* Saturate: '<S195>/Saturation1' */
      VCM20_B.Saturation1_f = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S195>/Saturation1' */
      VCM20_B.Saturation1_f = u1;
    } else {
      /* Saturate: '<S195>/Saturation1' */
      VCM20_B.Saturation1_f = u0_0;
    }

    /* End of Saturate: '<S195>/Saturation1' */
  }

  /* Product: '<S195>/Product' */
  VCM20_B.RLMotorPower_m = VCM20_B.Saturation_gc * VCM20_B.Saturation1_f;

  /* Gain: '<S195>/Gain' */
  VCM20_B.RLMotorPower_W_i = VCM20_P.Gain_Gain_lb * VCM20_B.RLMotorPower_m;

  /* Lookup_n-D: '<S195>/2-D Lookup Table' incorporates:
   *  Saturate: '<S195>/Saturation'
   *  Saturate: '<S195>/Saturation1'
   */
  VCM20_B.RLEfficency_f = look2_binlxpw(VCM20_B.Saturation_gc,
    VCM20_B.Saturation1_f, VCM20_P.uDLookupTable_bp01Data_j,
    VCM20_P.uDLookupTable_bp02Data_c, VCM20_P.uDLookupTable_tableData_ek,
    VCM20_P.uDLookupTable_maxIndex_g, 11U);

  /* Product: '<S195>/Divide' */
  VCM20_B.RLPowerConsumption_i = VCM20_B.RLMotorPower_W_i /
    VCM20_B.RLEfficency_f;

  /* Sum: '<S180>/Sum' */
  VCM20_B.Sum_a5 = ((VCM20_B.FRPowerConsumption_i + VCM20_B.FLPowerConsumption_i)
                    + VCM20_B.RRPowerConsumption_b) +
    VCM20_B.RLPowerConsumption_i;

  /* Saturate: '<S180>/Saturation' */
  u0_0 = VCM20_B.Sum_a5;
  u1 = VCM20_P.Saturation_LowerSat_au;
  u2 = VCM20_P.Saturation_UpperSat_ow;
  if (u0_0 > u2) {
    /* Saturate: '<S180>/Saturation' */
    VCM20_B.Saturation_gm = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S180>/Saturation' */
    VCM20_B.Saturation_gm = u1;
  } else {
    /* Saturate: '<S180>/Saturation' */
    VCM20_B.Saturation_gm = u0_0;
  }

  /* End of Saturate: '<S180>/Saturation' */

  /* Product: '<S180>/Divide3' */
  VCM20_B.Divide3_l = VCM20_B.torque_o / VCM20_B.Saturation_gm;

  /* Product: '<S180>/Product3' */
  VCM20_B.torque_i5 = VCM20_B.Divide3_l * VCM20_B.Eff_INV;

  /* Saturate: '<S180>/Saturation4' */
  u0_0 = VCM20_B.torque_i5;
  u1 = VCM20_P.Saturation4_LowerSat_b;
  u2 = VCM20_P.Saturation4_UpperSat_f;
  if (u0_0 > u2) {
    /* Saturate: '<S180>/Saturation4' */
    VCM20_B.torque_or = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S180>/Saturation4' */
    VCM20_B.torque_or = u1;
  } else {
    /* Saturate: '<S180>/Saturation4' */
    VCM20_B.torque_or = u0_0;
  }

  /* End of Saturate: '<S180>/Saturation4' */

  /* Product: '<S180>/Divide' */
  VCM20_B.Divide_mt = VCM20_B.torque_l / VCM20_B.Saturation_gm;

  /* Product: '<S180>/Product' */
  VCM20_B.torque_d = VCM20_B.Divide_mt * VCM20_B.Eff_INV;

  /* Saturate: '<S180>/Saturation1' */
  u0_0 = VCM20_B.torque_d;
  u1 = VCM20_P.Saturation1_LowerSat_ob;
  u2 = VCM20_P.Saturation1_UpperSat_n;
  if (u0_0 > u2) {
    /* Saturate: '<S180>/Saturation1' */
    VCM20_B.torque_d0 = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S180>/Saturation1' */
    VCM20_B.torque_d0 = u1;
  } else {
    /* Saturate: '<S180>/Saturation1' */
    VCM20_B.torque_d0 = u0_0;
  }

  /* End of Saturate: '<S180>/Saturation1' */

  /* Saturate: '<S198>/Saturation' */
  u0_0 = VCM20_B.torque_d0;
  u1 = VCM20_P.Saturation_LowerSat_lg;
  u2 = VCM20_P.Saturation_UpperSat_om;
  if (u0_0 > u2) {
    /* Saturate: '<S198>/Saturation' */
    VCM20_B.Saturation_e2 = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S198>/Saturation' */
    VCM20_B.Saturation_e2 = u1;
  } else {
    /* Saturate: '<S198>/Saturation' */
    VCM20_B.Saturation_e2 = u0_0;
  }

  /* End of Saturate: '<S198>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S198>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_n;
    u1 = VCM20_P.Saturation1_LowerSat_ey;
    u2 = VCM20_P.Saturation1_UpperSat_gq;
    if (u0_0 > u2) {
      /* Saturate: '<S198>/Saturation1' */
      VCM20_B.Saturation1_j = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S198>/Saturation1' */
      VCM20_B.Saturation1_j = u1;
    } else {
      /* Saturate: '<S198>/Saturation1' */
      VCM20_B.Saturation1_j = u0_0;
    }

    /* End of Saturate: '<S198>/Saturation1' */
  }

  /* Product: '<S198>/Product' */
  VCM20_B.FRMotorPower_h = VCM20_B.Saturation_e2 * VCM20_B.Saturation1_j;

  /* Gain: '<S198>/Gain' */
  VCM20_B.FRMotorPower_W_h = VCM20_P.Gain_Gain_f0 * VCM20_B.FRMotorPower_h;

  /* Lookup_n-D: '<S198>/2-D Lookup Table' incorporates:
   *  Saturate: '<S198>/Saturation'
   *  Saturate: '<S198>/Saturation1'
   */
  VCM20_B.FREfficency_n = look2_binlxpw(VCM20_B.Saturation_e2,
    VCM20_B.Saturation1_j, VCM20_P.uDLookupTable_bp01Data_k,
    VCM20_P.uDLookupTable_bp02Data_d, VCM20_P.uDLookupTable_tableData_ld,
    VCM20_P.uDLookupTable_maxIndex_k1, 11U);

  /* Product: '<S198>/Divide' */
  VCM20_B.FRPowerConsumption_m = VCM20_B.FRMotorPower_W_h /
    VCM20_B.FREfficency_n;

  /* Product: '<S180>/Divide1' */
  VCM20_B.Divide1_k = VCM20_B.torque_e / VCM20_B.Saturation_gm;

  /* Product: '<S180>/Product1' */
  VCM20_B.torque_ae = VCM20_B.Divide1_k * VCM20_B.Eff_INV;

  /* Saturate: '<S180>/Saturation2' */
  u0_0 = VCM20_B.torque_ae;
  u1 = VCM20_P.Saturation2_LowerSat_d;
  u2 = VCM20_P.Saturation2_UpperSat_c;
  if (u0_0 > u2) {
    /* Saturate: '<S180>/Saturation2' */
    VCM20_B.torque_mg = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S180>/Saturation2' */
    VCM20_B.torque_mg = u1;
  } else {
    /* Saturate: '<S180>/Saturation2' */
    VCM20_B.torque_mg = u0_0;
  }

  /* End of Saturate: '<S180>/Saturation2' */

  /* Saturate: '<S197>/Saturation' */
  u0_0 = VCM20_B.torque_mg;
  u1 = VCM20_P.Saturation_LowerSat_cl;
  u2 = VCM20_P.Saturation_UpperSat_a2;
  if (u0_0 > u2) {
    /* Saturate: '<S197>/Saturation' */
    VCM20_B.Saturation_o = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S197>/Saturation' */
    VCM20_B.Saturation_o = u1;
  } else {
    /* Saturate: '<S197>/Saturation' */
    VCM20_B.Saturation_o = u0_0;
  }

  /* End of Saturate: '<S197>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S197>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_g;
    u1 = VCM20_P.Saturation1_LowerSat_b;
    u2 = VCM20_P.Saturation1_UpperSat_bup;
    if (u0_0 > u2) {
      /* Saturate: '<S197>/Saturation1' */
      VCM20_B.Saturation1_f1 = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S197>/Saturation1' */
      VCM20_B.Saturation1_f1 = u1;
    } else {
      /* Saturate: '<S197>/Saturation1' */
      VCM20_B.Saturation1_f1 = u0_0;
    }

    /* End of Saturate: '<S197>/Saturation1' */
  }

  /* Product: '<S197>/Product' */
  VCM20_B.FLMotorPower_k = VCM20_B.Saturation_o * VCM20_B.Saturation1_f1;

  /* Gain: '<S197>/Gain' */
  VCM20_B.FLMotorPower_W_bv = VCM20_P.Gain_Gain_ik * VCM20_B.FLMotorPower_k;

  /* Lookup_n-D: '<S197>/2-D Lookup Table' incorporates:
   *  Saturate: '<S197>/Saturation'
   *  Saturate: '<S197>/Saturation1'
   */
  VCM20_B.FLEfficency_n = look2_binlxpw(VCM20_B.Saturation_o,
    VCM20_B.Saturation1_f1, VCM20_P.uDLookupTable_bp01Data_bz,
    VCM20_P.uDLookupTable_bp02Data_fw, VCM20_P.uDLookupTable_tableData_of,
    VCM20_P.uDLookupTable_maxIndex_l, 11U);

  /* Product: '<S197>/Divide' */
  VCM20_B.FLPowerConsumption_k = VCM20_B.FLMotorPower_W_bv /
    VCM20_B.FLEfficency_n;

  /* Product: '<S180>/Divide2' */
  VCM20_B.Divide2_ey = VCM20_B.torque_bc / VCM20_B.Saturation_gm;

  /* Product: '<S180>/Product2' */
  VCM20_B.torque_ms = VCM20_B.Divide2_ey * VCM20_B.Eff_INV;

  /* Saturate: '<S180>/Saturation3' */
  u0_0 = VCM20_B.torque_ms;
  u1 = VCM20_P.Saturation3_LowerSat_l;
  u2 = VCM20_P.Saturation3_UpperSat_e;
  if (u0_0 > u2) {
    /* Saturate: '<S180>/Saturation3' */
    VCM20_B.torque_k = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S180>/Saturation3' */
    VCM20_B.torque_k = u1;
  } else {
    /* Saturate: '<S180>/Saturation3' */
    VCM20_B.torque_k = u0_0;
  }

  /* End of Saturate: '<S180>/Saturation3' */

  /* Saturate: '<S200>/Saturation' */
  u0_0 = VCM20_B.torque_k;
  u1 = VCM20_P.Saturation_LowerSat_hb;
  u2 = VCM20_P.Saturation_UpperSat_nt;
  if (u0_0 > u2) {
    /* Saturate: '<S200>/Saturation' */
    VCM20_B.Saturation_co = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S200>/Saturation' */
    VCM20_B.Saturation_co = u1;
  } else {
    /* Saturate: '<S200>/Saturation' */
    VCM20_B.Saturation_co = u0_0;
  }

  /* End of Saturate: '<S200>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S200>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_j;
    u1 = VCM20_P.Saturation1_LowerSat_nn;
    u2 = VCM20_P.Saturation1_UpperSat_a;
    if (u0_0 > u2) {
      /* Saturate: '<S200>/Saturation1' */
      VCM20_B.Saturation1_ey = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S200>/Saturation1' */
      VCM20_B.Saturation1_ey = u1;
    } else {
      /* Saturate: '<S200>/Saturation1' */
      VCM20_B.Saturation1_ey = u0_0;
    }

    /* End of Saturate: '<S200>/Saturation1' */
  }

  /* Product: '<S200>/Product' */
  VCM20_B.RRMotorPower_p = VCM20_B.Saturation_co * VCM20_B.Saturation1_ey;

  /* Gain: '<S200>/Gain' */
  VCM20_B.RRMotorPower_W_i = VCM20_P.Gain_Gain_lm * VCM20_B.RRMotorPower_p;

  /* Lookup_n-D: '<S200>/2-D Lookup Table' incorporates:
   *  Saturate: '<S200>/Saturation'
   *  Saturate: '<S200>/Saturation1'
   */
  VCM20_B.RREfficency_ne = look2_binlxpw(VCM20_B.Saturation_co,
    VCM20_B.Saturation1_ey, VCM20_P.uDLookupTable_bp01Data_f,
    VCM20_P.uDLookupTable_bp02Data_cb, VCM20_P.uDLookupTable_tableData_e5,
    VCM20_P.uDLookupTable_maxIndex_nr, 11U);

  /* Product: '<S200>/Divide' */
  VCM20_B.RRPowerConsumption_bl = VCM20_B.RRMotorPower_W_i /
    VCM20_B.RREfficency_ne;

  /* Saturate: '<S199>/Saturation' */
  u0_0 = VCM20_B.torque_or;
  u1 = VCM20_P.Saturation_LowerSat_jy;
  u2 = VCM20_P.Saturation_UpperSat_i;
  if (u0_0 > u2) {
    /* Saturate: '<S199>/Saturation' */
    VCM20_B.Saturation_c2 = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S199>/Saturation' */
    VCM20_B.Saturation_c2 = u1;
  } else {
    /* Saturate: '<S199>/Saturation' */
    VCM20_B.Saturation_c2 = u0_0;
  }

  /* End of Saturate: '<S199>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S199>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9;
    u1 = VCM20_P.Saturation1_LowerSat_m;
    u2 = VCM20_P.Saturation1_UpperSat_ex;
    if (u0_0 > u2) {
      /* Saturate: '<S199>/Saturation1' */
      VCM20_B.Saturation1_b = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S199>/Saturation1' */
      VCM20_B.Saturation1_b = u1;
    } else {
      /* Saturate: '<S199>/Saturation1' */
      VCM20_B.Saturation1_b = u0_0;
    }

    /* End of Saturate: '<S199>/Saturation1' */
  }

  /* Product: '<S199>/Product' */
  VCM20_B.RLMotorPower_j = VCM20_B.Saturation_c2 * VCM20_B.Saturation1_b;

  /* Gain: '<S199>/Gain' */
  VCM20_B.RLMotorPower_W_m = VCM20_P.Gain_Gain_fl * VCM20_B.RLMotorPower_j;

  /* Lookup_n-D: '<S199>/2-D Lookup Table' incorporates:
   *  Saturate: '<S199>/Saturation'
   *  Saturate: '<S199>/Saturation1'
   */
  VCM20_B.RLEfficency_b = look2_binlxpw(VCM20_B.Saturation_c2,
    VCM20_B.Saturation1_b, VCM20_P.uDLookupTable_bp01Data_n,
    VCM20_P.uDLookupTable_bp02Data_n1, VCM20_P.uDLookupTable_tableData_ot,
    VCM20_P.uDLookupTable_maxIndex_b, 11U);

  /* Product: '<S199>/Divide' */
  VCM20_B.RLPowerConsumption_d = VCM20_B.RLMotorPower_W_m /
    VCM20_B.RLEfficency_b;

  /* Sum: '<S181>/Sum' */
  VCM20_B.Sum_e = ((VCM20_B.FRPowerConsumption_m + VCM20_B.FLPowerConsumption_k)
                   + VCM20_B.RRPowerConsumption_bl) +
    VCM20_B.RLPowerConsumption_d;

  /* Saturate: '<S181>/Saturation' */
  u0_0 = VCM20_B.Sum_e;
  u1 = VCM20_P.Saturation_LowerSat_pj;
  u2 = VCM20_P.Saturation_UpperSat_d;
  if (u0_0 > u2) {
    /* Saturate: '<S181>/Saturation' */
    VCM20_B.Saturation_ak = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S181>/Saturation' */
    VCM20_B.Saturation_ak = u1;
  } else {
    /* Saturate: '<S181>/Saturation' */
    VCM20_B.Saturation_ak = u0_0;
  }

  /* End of Saturate: '<S181>/Saturation' */

  /* Product: '<S181>/Divide3' */
  VCM20_B.Divide3_n = VCM20_B.torque_or / VCM20_B.Saturation_ak;

  /* Product: '<S181>/Product3' */
  VCM20_B.torque_m1 = VCM20_B.Divide3_n * VCM20_B.Eff_INV;

  /* Saturate: '<S181>/Saturation4' */
  u0_0 = VCM20_B.torque_m1;
  u1 = VCM20_P.Saturation4_LowerSat_l;
  u2 = VCM20_P.Saturation4_UpperSat_b;
  if (u0_0 > u2) {
    /* Saturate: '<S181>/Saturation4' */
    VCM20_B.torque_a4 = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S181>/Saturation4' */
    VCM20_B.torque_a4 = u1;
  } else {
    /* Saturate: '<S181>/Saturation4' */
    VCM20_B.torque_a4 = u0_0;
  }

  /* End of Saturate: '<S181>/Saturation4' */

  /* Product: '<S181>/Divide' */
  VCM20_B.Divide_d = VCM20_B.torque_d0 / VCM20_B.Saturation_ak;

  /* Product: '<S181>/Product' */
  VCM20_B.torque_m3 = VCM20_B.Divide_d * VCM20_B.Eff_INV;

  /* Saturate: '<S181>/Saturation1' */
  u0_0 = VCM20_B.torque_m3;
  u1 = VCM20_P.Saturation1_LowerSat_ng;
  u2 = VCM20_P.Saturation1_UpperSat_j;
  if (u0_0 > u2) {
    /* Saturate: '<S181>/Saturation1' */
    VCM20_B.torque_cl = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S181>/Saturation1' */
    VCM20_B.torque_cl = u1;
  } else {
    /* Saturate: '<S181>/Saturation1' */
    VCM20_B.torque_cl = u0_0;
  }

  /* End of Saturate: '<S181>/Saturation1' */

  /* Saturate: '<S202>/Saturation' */
  u0_0 = VCM20_B.torque_cl;
  u1 = VCM20_P.Saturation_LowerSat_jd;
  u2 = VCM20_P.Saturation_UpperSat_du;
  if (u0_0 > u2) {
    /* Saturate: '<S202>/Saturation' */
    VCM20_B.Saturation_p5 = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S202>/Saturation' */
    VCM20_B.Saturation_p5 = u1;
  } else {
    /* Saturate: '<S202>/Saturation' */
    VCM20_B.Saturation_p5 = u0_0;
  }

  /* End of Saturate: '<S202>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S202>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_n;
    u1 = VCM20_P.Saturation1_LowerSat_a;
    u2 = VCM20_P.Saturation1_UpperSat_ij;
    if (u0_0 > u2) {
      /* Saturate: '<S202>/Saturation1' */
      VCM20_B.Saturation1_i3 = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S202>/Saturation1' */
      VCM20_B.Saturation1_i3 = u1;
    } else {
      /* Saturate: '<S202>/Saturation1' */
      VCM20_B.Saturation1_i3 = u0_0;
    }

    /* End of Saturate: '<S202>/Saturation1' */
  }

  /* Product: '<S202>/Product' */
  VCM20_B.FRMotorPower_h3 = VCM20_B.Saturation_p5 * VCM20_B.Saturation1_i3;

  /* Gain: '<S202>/Gain' */
  VCM20_B.FRMotorPower_W_f = VCM20_P.Gain_Gain_d * VCM20_B.FRMotorPower_h3;

  /* Lookup_n-D: '<S202>/2-D Lookup Table' incorporates:
   *  Saturate: '<S202>/Saturation'
   *  Saturate: '<S202>/Saturation1'
   */
  VCM20_B.FREfficency_f = look2_binlxpw(VCM20_B.Saturation_p5,
    VCM20_B.Saturation1_i3, VCM20_P.uDLookupTable_bp01Data_o,
    VCM20_P.uDLookupTable_bp02Data_i, VCM20_P.uDLookupTable_tableData_n,
    VCM20_P.uDLookupTable_maxIndex_ga, 11U);

  /* Product: '<S202>/Divide' */
  VCM20_B.FRPowerConsumption_c = VCM20_B.FRMotorPower_W_f /
    VCM20_B.FREfficency_f;

  /* Product: '<S181>/Divide1' */
  VCM20_B.Divide1_iv = VCM20_B.torque_mg / VCM20_B.Saturation_ak;

  /* Product: '<S181>/Product1' */
  VCM20_B.torque_ju = VCM20_B.Divide1_iv * VCM20_B.Eff_INV;

  /* Saturate: '<S181>/Saturation2' */
  u0_0 = VCM20_B.torque_ju;
  u1 = VCM20_P.Saturation2_LowerSat_f;
  u2 = VCM20_P.Saturation2_UpperSat_j;
  if (u0_0 > u2) {
    /* Saturate: '<S181>/Saturation2' */
    VCM20_B.torque_fy = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S181>/Saturation2' */
    VCM20_B.torque_fy = u1;
  } else {
    /* Saturate: '<S181>/Saturation2' */
    VCM20_B.torque_fy = u0_0;
  }

  /* End of Saturate: '<S181>/Saturation2' */

  /* Saturate: '<S201>/Saturation' */
  u0_0 = VCM20_B.torque_fy;
  u1 = VCM20_P.Saturation_LowerSat_mb;
  u2 = VCM20_P.Saturation_UpperSat_aq;
  if (u0_0 > u2) {
    /* Saturate: '<S201>/Saturation' */
    VCM20_B.Saturation_i = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S201>/Saturation' */
    VCM20_B.Saturation_i = u1;
  } else {
    /* Saturate: '<S201>/Saturation' */
    VCM20_B.Saturation_i = u0_0;
  }

  /* End of Saturate: '<S201>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S201>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_g;
    u1 = VCM20_P.Saturation1_LowerSat_l;
    u2 = VCM20_P.Saturation1_UpperSat_ee;
    if (u0_0 > u2) {
      /* Saturate: '<S201>/Saturation1' */
      VCM20_B.Saturation1_d = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S201>/Saturation1' */
      VCM20_B.Saturation1_d = u1;
    } else {
      /* Saturate: '<S201>/Saturation1' */
      VCM20_B.Saturation1_d = u0_0;
    }

    /* End of Saturate: '<S201>/Saturation1' */
  }

  /* Product: '<S201>/Product' */
  VCM20_B.FLMotorPower_d = VCM20_B.Saturation_i * VCM20_B.Saturation1_d;

  /* Gain: '<S201>/Gain' */
  VCM20_B.FLMotorPower_W_o = VCM20_P.Gain_Gain_lg * VCM20_B.FLMotorPower_d;

  /* Lookup_n-D: '<S201>/2-D Lookup Table' incorporates:
   *  Saturate: '<S201>/Saturation'
   *  Saturate: '<S201>/Saturation1'
   */
  VCM20_B.FLEfficency_i = look2_binlxpw(VCM20_B.Saturation_i,
    VCM20_B.Saturation1_d, VCM20_P.uDLookupTable_bp01Data_c,
    VCM20_P.uDLookupTable_bp02Data_h, VCM20_P.uDLookupTable_tableData_gk,
    VCM20_P.uDLookupTable_maxIndex_o, 11U);

  /* Product: '<S201>/Divide' */
  VCM20_B.FLPowerConsumption_l = VCM20_B.FLMotorPower_W_o /
    VCM20_B.FLEfficency_i;

  /* Product: '<S181>/Divide2' */
  VCM20_B.Divide2_j = VCM20_B.torque_k / VCM20_B.Saturation_ak;

  /* Product: '<S181>/Product2' */
  VCM20_B.torque_n = VCM20_B.Divide2_j * VCM20_B.Eff_INV;

  /* Saturate: '<S181>/Saturation3' */
  u0_0 = VCM20_B.torque_n;
  u1 = VCM20_P.Saturation3_LowerSat_ds;
  u2 = VCM20_P.Saturation3_UpperSat_ma;
  if (u0_0 > u2) {
    /* Saturate: '<S181>/Saturation3' */
    VCM20_B.torque_ag = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S181>/Saturation3' */
    VCM20_B.torque_ag = u1;
  } else {
    /* Saturate: '<S181>/Saturation3' */
    VCM20_B.torque_ag = u0_0;
  }

  /* End of Saturate: '<S181>/Saturation3' */

  /* Saturate: '<S204>/Saturation' */
  u0_0 = VCM20_B.torque_ag;
  u1 = VCM20_P.Saturation_LowerSat_ib;
  u2 = VCM20_P.Saturation_UpperSat_ntr;
  if (u0_0 > u2) {
    /* Saturate: '<S204>/Saturation' */
    VCM20_B.Saturation_ha = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S204>/Saturation' */
    VCM20_B.Saturation_ha = u1;
  } else {
    /* Saturate: '<S204>/Saturation' */
    VCM20_B.Saturation_ha = u0_0;
  }

  /* End of Saturate: '<S204>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S204>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9_j;
    u1 = VCM20_P.Saturation1_LowerSat_kj;
    u2 = VCM20_P.Saturation1_UpperSat_ae;
    if (u0_0 > u2) {
      /* Saturate: '<S204>/Saturation1' */
      VCM20_B.Saturation1_fg = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S204>/Saturation1' */
      VCM20_B.Saturation1_fg = u1;
    } else {
      /* Saturate: '<S204>/Saturation1' */
      VCM20_B.Saturation1_fg = u0_0;
    }

    /* End of Saturate: '<S204>/Saturation1' */
  }

  /* Product: '<S204>/Product' */
  VCM20_B.RRMotorPower_m = VCM20_B.Saturation_ha * VCM20_B.Saturation1_fg;

  /* Gain: '<S204>/Gain' */
  VCM20_B.RRMotorPower_W_e = VCM20_P.Gain_Gain_o * VCM20_B.RRMotorPower_m;

  /* Lookup_n-D: '<S204>/2-D Lookup Table' incorporates:
   *  Saturate: '<S204>/Saturation'
   *  Saturate: '<S204>/Saturation1'
   */
  VCM20_B.RREfficency_d1 = look2_binlxpw(VCM20_B.Saturation_ha,
    VCM20_B.Saturation1_fg, VCM20_P.uDLookupTable_bp01Data_mg,
    VCM20_P.uDLookupTable_bp02Data_cu, VCM20_P.uDLookupTable_tableData_dj,
    VCM20_P.uDLookupTable_maxIndex_ln, 11U);

  /* Product: '<S204>/Divide' */
  VCM20_B.RRPowerConsumption_kk = VCM20_B.RRMotorPower_W_e /
    VCM20_B.RREfficency_d1;

  /* Saturate: '<S203>/Saturation' */
  u0_0 = VCM20_B.torque_a4;
  u1 = VCM20_P.Saturation_LowerSat_kd;
  u2 = VCM20_P.Saturation_UpperSat_aw;
  if (u0_0 > u2) {
    /* Saturate: '<S203>/Saturation' */
    VCM20_B.Saturation_f = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S203>/Saturation' */
    VCM20_B.Saturation_f = u1;
  } else {
    /* Saturate: '<S203>/Saturation' */
    VCM20_B.Saturation_f = u0_0;
  }

  /* End of Saturate: '<S203>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Saturate: '<S203>/Saturation1' */
    u0_0 = VCM20_B.SFunction1_o9;
    u1 = VCM20_P.Saturation1_LowerSat_ii;
    u2 = VCM20_P.Saturation1_UpperSat_e1;
    if (u0_0 > u2) {
      /* Saturate: '<S203>/Saturation1' */
      VCM20_B.Saturation1_fj = u2;
    } else if (u0_0 < u1) {
      /* Saturate: '<S203>/Saturation1' */
      VCM20_B.Saturation1_fj = u1;
    } else {
      /* Saturate: '<S203>/Saturation1' */
      VCM20_B.Saturation1_fj = u0_0;
    }

    /* End of Saturate: '<S203>/Saturation1' */
  }

  /* Product: '<S203>/Product' */
  VCM20_B.RLMotorPower_p = VCM20_B.Saturation_f * VCM20_B.Saturation1_fj;

  /* Gain: '<S203>/Gain' */
  VCM20_B.RLMotorPower_W_n = VCM20_P.Gain_Gain_bj * VCM20_B.RLMotorPower_p;

  /* Lookup_n-D: '<S203>/2-D Lookup Table' incorporates:
   *  Saturate: '<S203>/Saturation'
   *  Saturate: '<S203>/Saturation1'
   */
  VCM20_B.RLEfficency_a = look2_binlxpw(VCM20_B.Saturation_f,
    VCM20_B.Saturation1_fj, VCM20_P.uDLookupTable_bp01Data_d0,
    VCM20_P.uDLookupTable_bp02Data_m, VCM20_P.uDLookupTable_tableData_hn,
    VCM20_P.uDLookupTable_maxIndex_of, 11U);

  /* Product: '<S203>/Divide' */
  VCM20_B.RLPowerConsumption_b = VCM20_B.RLMotorPower_W_n /
    VCM20_B.RLEfficency_a;

  /* Sum: '<S182>/Sum' */
  VCM20_B.Sum_b = ((VCM20_B.FRPowerConsumption_c + VCM20_B.FLPowerConsumption_l)
                   + VCM20_B.RRPowerConsumption_kk) +
    VCM20_B.RLPowerConsumption_b;

  /* Saturate: '<S182>/Saturation' */
  u0_0 = VCM20_B.Sum_b;
  u1 = VCM20_P.Saturation_LowerSat_o1;
  u2 = VCM20_P.Saturation_UpperSat_lo;
  if (u0_0 > u2) {
    /* Saturate: '<S182>/Saturation' */
    VCM20_B.Saturation_l5 = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S182>/Saturation' */
    VCM20_B.Saturation_l5 = u1;
  } else {
    /* Saturate: '<S182>/Saturation' */
    VCM20_B.Saturation_l5 = u0_0;
  }

  /* End of Saturate: '<S182>/Saturation' */

  /* Product: '<S182>/Divide3' */
  VCM20_B.Divide3_d = VCM20_B.torque_a4 / VCM20_B.Saturation_l5;

  /* Product: '<S182>/Product3' */
  VCM20_B.torque_dc = VCM20_B.Divide3_d * VCM20_B.Eff_INV;

  /* Saturate: '<S182>/Saturation4' */
  u0_0 = VCM20_B.torque_dc;
  u1 = VCM20_P.Saturation4_LowerSat_p;
  u2 = VCM20_P.Saturation4_UpperSat_o;
  if (u0_0 > u2) {
    /* Saturate: '<S182>/Saturation4' */
    VCM20_B.torque_ej = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S182>/Saturation4' */
    VCM20_B.torque_ej = u1;
  } else {
    /* Saturate: '<S182>/Saturation4' */
    VCM20_B.torque_ej = u0_0;
  }

  /* End of Saturate: '<S182>/Saturation4' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S112>/S-Function1' */
    /* This comment workarounds a code generation problem */
    {
      /* dSPACE I/O Board DS1401STDADCT4 #1 Unit:ADC Group:ADC */
      adc_tp4_single_new_read(ADC_TP4_1_MODULE_ADDR,
        ADC_TP4_CH4,
        (dsfloat *)&VCM20_B.SFunction1);
    }

    /* Fcn: '<S14>/Fcn1' */
    VCM20_B.Current = (VCM20_B.SFunction1 - 2.5) / 7.9;

    /* Product: '<S173>/Product' */
    VCM20_B.Product_a = VCM20_B.SFunction1_o10 * VCM20_B.Current;

    /* Sum: '<S183>/Sum' incorporates:
     *  Constant: '<S183>/Constant'
     */
    VCM20_B.Sum_aa = VCM20_P.Constant_Value_lo - VCM20_B.Product_a;

    /* RelationalOperator: '<S206>/Compare' incorporates:
     *  Constant: '<S206>/Constant'
     */
    VCM20_B.Compare_b = (VCM20_B.Sum_aa < VCM20_P.Constant_Value_h);

    /* Gain: '<S183>/Gain' */
    VCM20_B.Gain_p = VCM20_P.Gain_Gain_fw * VCM20_B.Sum_aa;

    /* Product: '<S183>/Divide' */
    VCM20_B.Divide_hj = (real_T)VCM20_B.Compare_b * VCM20_B.Gain_p;
  }

  /* Sum: '<S183>/Sum4' */
  VCM20_B.Sum4 = VCM20_B.torque_ej + VCM20_B.Divide_hj;

  /* Switch: '<S183>/Switch2' */
  if (VCM20_B.Switch2_ot > VCM20_P.Switch2_Threshold_m) {
    /* MinMax: '<S183>/Min2' */
    u0_0 = VCM20_B.Switch2_ot;
    u1_0 = VCM20_B.Sum4;
    if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S183>/Min2' */
    VCM20_B.Min2_a = u1_0;

    /* Switch: '<S183>/Switch2' */
    VCM20_B.torque_lx = VCM20_B.Min2_a;
  } else {
    /* MinMax: '<S183>/Max2' */
    u0_0 = VCM20_B.Switch2_ot;
    u1_0 = VCM20_B.Sum4;
    if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S183>/Max2' */
    VCM20_B.Max2 = u1_0;

    /* Switch: '<S183>/Switch2' */
    VCM20_B.torque_lx = VCM20_B.Max2;
  }

  /* End of Switch: '<S183>/Switch2' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* RelationalOperator: '<S210>/Compare' incorporates:
     *  Constant: '<S210>/Constant'
     */
    VCM20_B.Compare_pg = (VCM20_B.SFunction1_o9 <=
                          VCM20_P.CompareToConstant3_const_h);
  }

  /* Switch: '<S184>/Switch3' incorporates:
   *  Constant: '<Root>/80kWOn'
   */
  if (VCM20_P.u0kWOn_Value > VCM20_P.Switch3_Threshold_l) {
    /* Switch: '<S205>/Switch7' */
    if (VCM20_B.Compare_pg) {
      /* Switch: '<S205>/Switch7' */
      VCM20_B.FLtrq_k = VCM20_B.torque_lx;
    } else {
      /* Saturate: '<S205>/Saturation3' */
      u0_0 = VCM20_B.torque_lx;
      u1 = VCM20_P.Saturation3_LowerSat;
      u2 = VCM20_P.Saturation3_UpperSat;
      if (u0_0 > u2) {
        /* Saturate: '<S205>/Saturation3' */
        VCM20_B.Saturation3_b = u2;
      } else if (u0_0 < u1) {
        /* Saturate: '<S205>/Saturation3' */
        VCM20_B.Saturation3_b = u1;
      } else {
        /* Saturate: '<S205>/Saturation3' */
        VCM20_B.Saturation3_b = u0_0;
      }

      /* End of Saturate: '<S205>/Saturation3' */

      /* Switch: '<S205>/Switch7' */
      VCM20_B.FLtrq_k = VCM20_B.Saturation3_b;
    }

    /* End of Switch: '<S205>/Switch7' */

    /* Switch: '<S184>/Switch3' */
    VCM20_B.FLtorque_c = VCM20_B.FLtrq_k;
  } else {
    /* Switch: '<S184>/Switch3' */
    VCM20_B.FLtorque_c = VCM20_B.Switch2_ot;
  }

  /* End of Switch: '<S184>/Switch3' */

  /* Gain: '<S122>/Gain3' */
  VCM20_B.FLtrq = VCM20_P.Gain3_Gain * VCM20_B.FLtorque_c;

  /* Gain: '<S37>/Gain30' */
  VCM20_B.Gain30 = VCM20_P.Gain30_Gain * VCM20_B.FLtrq;

  /* Product: '<S182>/Divide2' */
  VCM20_B.Divide2_k = VCM20_B.torque_ag / VCM20_B.Saturation_l5;

  /* Product: '<S182>/Product2' */
  VCM20_B.torque_gr = VCM20_B.Divide2_k * VCM20_B.Eff_INV;

  /* Saturate: '<S182>/Saturation3' */
  u0_0 = VCM20_B.torque_gr;
  u1 = VCM20_P.Saturation3_LowerSat_b;
  u2 = VCM20_P.Saturation3_UpperSat_kj;
  if (u0_0 > u2) {
    /* Saturate: '<S182>/Saturation3' */
    VCM20_B.torque_fj = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S182>/Saturation3' */
    VCM20_B.torque_fj = u1;
  } else {
    /* Saturate: '<S182>/Saturation3' */
    VCM20_B.torque_fj = u0_0;
  }

  /* End of Saturate: '<S182>/Saturation3' */

  /* Sum: '<S183>/Sum3' */
  VCM20_B.Sum3_n = VCM20_B.torque_fj + VCM20_B.Divide_hj;

  /* Switch: '<S183>/Switch3' */
  if (VCM20_B.Switch2_hr > VCM20_P.Switch3_Threshold_le) {
    /* MinMax: '<S183>/Min3' */
    u0_0 = VCM20_B.Switch2_hr;
    u1_0 = VCM20_B.Sum3_n;
    if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S183>/Min3' */
    VCM20_B.Min3_h = u1_0;

    /* Switch: '<S183>/Switch3' */
    VCM20_B.torque_ok = VCM20_B.Min3_h;
  } else {
    /* MinMax: '<S183>/Max3' */
    u0_0 = VCM20_B.Switch2_hr;
    u1_0 = VCM20_B.Sum3_n;
    if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S183>/Max3' */
    VCM20_B.Max3 = u1_0;

    /* Switch: '<S183>/Switch3' */
    VCM20_B.torque_ok = VCM20_B.Max3;
  }

  /* End of Switch: '<S183>/Switch3' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* RelationalOperator: '<S209>/Compare' incorporates:
     *  Constant: '<S209>/Constant'
     */
    VCM20_B.Compare_a = (VCM20_B.SFunction1_o9_j <=
                         VCM20_P.CompareToConstant2_const);
  }

  /* Switch: '<S184>/Switch2' incorporates:
   *  Constant: '<Root>/80kWOn'
   */
  if (VCM20_P.u0kWOn_Value > VCM20_P.Switch2_Threshold_g) {
    /* Switch: '<S205>/Switch6' */
    if (VCM20_B.Compare_a) {
      /* Switch: '<S205>/Switch6' */
      VCM20_B.FRtrq_i = VCM20_B.torque_ok;
    } else {
      /* Saturate: '<S205>/Saturation2' */
      u0_0 = VCM20_B.torque_ok;
      u1 = VCM20_P.Saturation2_LowerSat_g;
      u2 = VCM20_P.Saturation2_UpperSat_b;
      if (u0_0 > u2) {
        /* Saturate: '<S205>/Saturation2' */
        VCM20_B.Saturation2_f = u2;
      } else if (u0_0 < u1) {
        /* Saturate: '<S205>/Saturation2' */
        VCM20_B.Saturation2_f = u1;
      } else {
        /* Saturate: '<S205>/Saturation2' */
        VCM20_B.Saturation2_f = u0_0;
      }

      /* End of Saturate: '<S205>/Saturation2' */

      /* Switch: '<S205>/Switch6' */
      VCM20_B.FRtrq_i = VCM20_B.Saturation2_f;
    }

    /* End of Switch: '<S205>/Switch6' */

    /* Switch: '<S184>/Switch2' */
    VCM20_B.FRtorque_a = VCM20_B.FRtrq_i;
  } else {
    /* Switch: '<S184>/Switch2' */
    VCM20_B.FRtorque_a = VCM20_B.Switch2_hr;
  }

  /* End of Switch: '<S184>/Switch2' */

  /* Gain: '<S122>/Gain2' */
  VCM20_B.FRtrq = VCM20_P.Gain2_Gain_o * VCM20_B.FRtorque_a;

  /* Gain: '<S37>/Gain31' */
  VCM20_B.Gain31 = VCM20_P.Gain31_Gain * VCM20_B.FRtrq;

  /* Product: '<S182>/Divide' */
  VCM20_B.Divide_d5 = VCM20_B.torque_cl / VCM20_B.Saturation_l5;

  /* Product: '<S182>/Product' */
  VCM20_B.torque_dz = VCM20_B.Divide_d5 * VCM20_B.Eff_INV;

  /* Saturate: '<S182>/Saturation1' */
  u0_0 = VCM20_B.torque_dz;
  u1 = VCM20_P.Saturation1_LowerSat_c;
  u2 = VCM20_P.Saturation1_UpperSat_f;
  if (u0_0 > u2) {
    /* Saturate: '<S182>/Saturation1' */
    VCM20_B.torque_h = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S182>/Saturation1' */
    VCM20_B.torque_h = u1;
  } else {
    /* Saturate: '<S182>/Saturation1' */
    VCM20_B.torque_h = u0_0;
  }

  /* End of Saturate: '<S182>/Saturation1' */

  /* Sum: '<S183>/Sum1' */
  VCM20_B.Sum1_l = VCM20_B.torque_h + VCM20_B.Divide_hj;

  /* Switch: '<S183>/Switch' */
  if (VCM20_B.Switch2_g > VCM20_P.Switch_Threshold_in) {
    /* MinMax: '<S183>/Min' */
    u0_0 = VCM20_B.Switch2_g;
    u1_0 = VCM20_B.Sum1_l;
    if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S183>/Min' */
    VCM20_B.Min_b = u1_0;

    /* Switch: '<S183>/Switch' */
    VCM20_B.torque_co = VCM20_B.Min_b;
  } else {
    /* MinMax: '<S183>/Max' */
    u0_0 = VCM20_B.Switch2_g;
    u1_0 = VCM20_B.Sum1_l;
    if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S183>/Max' */
    VCM20_B.Max = u1_0;

    /* Switch: '<S183>/Switch' */
    VCM20_B.torque_co = VCM20_B.Max;
  }

  /* End of Switch: '<S183>/Switch' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* RelationalOperator: '<S207>/Compare' incorporates:
     *  Constant: '<S207>/Constant'
     */
    VCM20_B.Compare_h = (VCM20_B.SFunction1_o9_n <=
                         VCM20_P.CompareToConstant_const_d);
  }

  /* Switch: '<S184>/Switch' incorporates:
   *  Constant: '<Root>/80kWOn'
   */
  if (VCM20_P.u0kWOn_Value > VCM20_P.Switch_Threshold_h) {
    /* Switch: '<S205>/Switch4' */
    if (VCM20_B.Compare_h) {
      /* Switch: '<S205>/Switch4' */
      VCM20_B.RRtrq_c = VCM20_B.torque_co;
    } else {
      /* Saturate: '<S205>/Saturation' */
      u0_0 = VCM20_B.torque_co;
      u1 = VCM20_P.Saturation_LowerSat_n;
      u2 = VCM20_P.Saturation_UpperSat_a;
      if (u0_0 > u2) {
        /* Saturate: '<S205>/Saturation' */
        VCM20_B.Saturation_ch = u2;
      } else if (u0_0 < u1) {
        /* Saturate: '<S205>/Saturation' */
        VCM20_B.Saturation_ch = u1;
      } else {
        /* Saturate: '<S205>/Saturation' */
        VCM20_B.Saturation_ch = u0_0;
      }

      /* End of Saturate: '<S205>/Saturation' */

      /* Switch: '<S205>/Switch4' */
      VCM20_B.RRtrq_c = VCM20_B.Saturation_ch;
    }

    /* End of Switch: '<S205>/Switch4' */

    /* Switch: '<S184>/Switch' */
    VCM20_B.RRtorque_e = VCM20_B.RRtrq_c;
  } else {
    /* Switch: '<S184>/Switch' */
    VCM20_B.RRtorque_e = VCM20_B.Switch2_g;
  }

  /* End of Switch: '<S184>/Switch' */

  /* Gain: '<S122>/Gain' */
  VCM20_B.RRtrq = VCM20_P.Gain_Gain_oh * VCM20_B.RRtorque_e;

  /* Gain: '<S37>/Gain32' */
  VCM20_B.Gain32 = VCM20_P.Gain32_Gain * VCM20_B.RRtrq;

  /* Product: '<S182>/Divide1' */
  VCM20_B.Divide1_h = VCM20_B.torque_fy / VCM20_B.Saturation_l5;

  /* Product: '<S182>/Product1' */
  VCM20_B.torque_ka = VCM20_B.Divide1_h * VCM20_B.Eff_INV;

  /* Saturate: '<S182>/Saturation2' */
  u0_0 = VCM20_B.torque_ka;
  u1 = VCM20_P.Saturation2_LowerSat_ir;
  u2 = VCM20_P.Saturation2_UpperSat_jl;
  if (u0_0 > u2) {
    /* Saturate: '<S182>/Saturation2' */
    VCM20_B.torque_ku = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S182>/Saturation2' */
    VCM20_B.torque_ku = u1;
  } else {
    /* Saturate: '<S182>/Saturation2' */
    VCM20_B.torque_ku = u0_0;
  }

  /* End of Saturate: '<S182>/Saturation2' */

  /* Sum: '<S183>/Sum2' */
  VCM20_B.Sum2_g = VCM20_B.torque_ku + VCM20_B.Divide_hj;

  /* Switch: '<S183>/Switch1' */
  if (VCM20_B.Switch2_lc > VCM20_P.Switch1_Threshold_d1) {
    /* MinMax: '<S183>/Min1' */
    u0_0 = VCM20_B.Switch2_lc;
    u1_0 = VCM20_B.Sum2_g;
    if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S183>/Min1' */
    VCM20_B.Min1_l = u1_0;

    /* Switch: '<S183>/Switch1' */
    VCM20_B.torque_nn = VCM20_B.Min1_l;
  } else {
    /* MinMax: '<S183>/Max1' */
    u0_0 = VCM20_B.Switch2_lc;
    u1_0 = VCM20_B.Sum2_g;
    if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S183>/Max1' */
    VCM20_B.Max1 = u1_0;

    /* Switch: '<S183>/Switch1' */
    VCM20_B.torque_nn = VCM20_B.Max1;
  }

  /* End of Switch: '<S183>/Switch1' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* RelationalOperator: '<S208>/Compare' incorporates:
     *  Constant: '<S208>/Constant'
     */
    VCM20_B.Compare_cc = (VCM20_B.SFunction1_o9_g <=
                          VCM20_P.CompareToConstant1_const_l);
  }

  /* Switch: '<S184>/Switch1' incorporates:
   *  Constant: '<Root>/80kWOn'
   */
  if (VCM20_P.u0kWOn_Value > VCM20_P.Switch1_Threshold_c) {
    /* Switch: '<S205>/Switch5' */
    if (VCM20_B.Compare_cc) {
      /* Switch: '<S205>/Switch5' */
      VCM20_B.RLtrq_c = VCM20_B.torque_nn;
    } else {
      /* Saturate: '<S205>/Saturation1' */
      u0_0 = VCM20_B.torque_nn;
      u1 = VCM20_P.Saturation1_LowerSat;
      u2 = VCM20_P.Saturation1_UpperSat;
      if (u0_0 > u2) {
        /* Saturate: '<S205>/Saturation1' */
        VCM20_B.Saturation1_pl = u2;
      } else if (u0_0 < u1) {
        /* Saturate: '<S205>/Saturation1' */
        VCM20_B.Saturation1_pl = u1;
      } else {
        /* Saturate: '<S205>/Saturation1' */
        VCM20_B.Saturation1_pl = u0_0;
      }

      /* End of Saturate: '<S205>/Saturation1' */

      /* Switch: '<S205>/Switch5' */
      VCM20_B.RLtrq_c = VCM20_B.Saturation1_pl;
    }

    /* End of Switch: '<S205>/Switch5' */

    /* Switch: '<S184>/Switch1' */
    VCM20_B.RLtorque_d = VCM20_B.RLtrq_c;
  } else {
    /* Switch: '<S184>/Switch1' */
    VCM20_B.RLtorque_d = VCM20_B.Switch2_lc;
  }

  /* End of Switch: '<S184>/Switch1' */

  /* Gain: '<S122>/Gain1' */
  VCM20_B.RLtrq = VCM20_P.Gain1_Gain_e * VCM20_B.RLtorque_d;

  /* Gain: '<S37>/Gain33' */
  VCM20_B.Gain33 = VCM20_P.Gain33_Gain * VCM20_B.RLtrq;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S57>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s004" Id:1284 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X504]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X504]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X504]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X504]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "TargetTorqueFL" (0|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain30 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain30 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain30 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "TargetTorqueFR" (16|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain31 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain31 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain31 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "TargetTorqueRR" (32|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain32 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain32 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain32 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "TargetTorqueRL" (48|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain33 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain33 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain33 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X504], 8,
                       &(CAN_Msg[0]), delayTime);
    }
  }

  /* Gain: '<S37>/Gain1' */
  VCM20_B.Gain1_e = VCM20_P.Gain1_Gain_b * VCM20_B.DYC_deltaTorque;

  /* TransferFcn: '<S220>/Transfer Fcn' */
  VCM20_B.vms = 0.0;
  VCM20_B.vms += VCM20_P.TransferFcn_C_l * VCM20_X.TransferFcn_CSTATE_j;

  /* Product: '<S222>/Product5' */
  VCM20_B.Product5 = VCM20_B.vms * VCM20_B.vms;

  /* Product: '<S222>/Divide3' incorporates:
   *  Constant: '<S216>/StabilityFactor_ref'
   */
  VCM20_B.Divide3_j = VCM20_B.Product5 * VCM20_P.StabilityFactor_ref_Value;

  /* Sum: '<S222>/Add4' incorporates:
   *  Constant: '<S222>/Constant3'
   */
  VCM20_B.Add4_j = VCM20_P.Constant3_Value + VCM20_B.Divide3_j;

  /* Product: '<S222>/Divide2' */
  VCM20_B.Divide2_o = VCM20_B.vms / VCM20_B.Add4_j;

  /* Gain: '<S222>/Gain6' */
  VCM20_B.Gain6 = VCM20_P.Gain6_Gain * VCM20_B.Divide2_o;

  /* Product: '<S222>/Divide' */
  VCM20_B.Divide_a = VCM20_B.Gain6 * VCM20_B.SteerAngle_rad;

  /* Gain: '<S37>/Gain2' */
  VCM20_B.Gain2_b = VCM20_P.Gain2_Gain_c * VCM20_B.Divide_a;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Gain: '<S176>/Gain' */
    VCM20_B.powerconsumptionkW = VCM20_P.Gain_Gain_k * VCM20_B.Product_a;

    /* Gain: '<S37>/Gain3' */
    VCM20_B.Gain3 = VCM20_P.Gain3_Gain_c * VCM20_B.powerconsumptionkW;

    /* S-Function (rti_commonblock): '<S58>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s004" Id:1313 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X521]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X521]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X521]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X521]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "DYC_deltaTorque" (0|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain1_e ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain1_e ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain1_e ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "r_ref" (16|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain2_b ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain2_b ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain2_b ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "PowerConsumption" (32|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain3 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain3 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain3 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X521], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* Gain: '<S322>/toTire_rpm' */
    VCM20_B.toTire_rpm_c = VCM20_P.toTire_rpm_Gain_g * VCM20_B.VELOCITY_X_p;

    /* Gain: '<S322>/Gear' */
    VCM20_B.Gear_nh = VCM20_P.Gear_Gain_e * VCM20_B.toTire_rpm_c;

    /* Sum: '<S322>/Add' */
    VCM20_B.Add_i3 = VCM20_B.SFunction1_o9_g - VCM20_B.Gear_nh;

    /* RelationalOperator: '<S327>/Compare' incorporates:
     *  Constant: '<S327>/Constant'
     */
    VCM20_B.Compare_f3 = (VCM20_B.SFunction1_o9_g == VCM20_P.Constant_Value_h3);

    /* Switch: '<S326>/Switch' */
    if (VCM20_B.Compare_f3) {
      /* Switch: '<S326>/Switch' incorporates:
       *  Constant: '<S326>/Constant'
       */
      VCM20_B.Switch_i0 = VCM20_P.Constant_Value_g;
    } else {
      /* Switch: '<S326>/Switch' */
      VCM20_B.Switch_i0 = VCM20_B.SFunction1_o9_g;
    }

    /* End of Switch: '<S326>/Switch' */

    /* Product: '<S322>/Divide' */
    VCM20_B.FL_Sliprate = VCM20_B.Add_i3 / VCM20_B.Switch_i0;

    /* Gain: '<S37>/Gain4' */
    VCM20_B.Gain4_j = VCM20_P.Gain4_Gain_k * VCM20_B.FL_Sliprate;

    /* Gain: '<S323>/toTire_rpm' */
    VCM20_B.toTire_rpm_lr = VCM20_P.toTire_rpm_Gain_b * VCM20_B.VELOCITY_X_p;

    /* Gain: '<S323>/Gear' */
    VCM20_B.Gear_dd = VCM20_P.Gear_Gain_p * VCM20_B.toTire_rpm_lr;

    /* Sum: '<S323>/Add' */
    VCM20_B.Add_k = VCM20_B.SFunction1_o9 - VCM20_B.Gear_dd;

    /* RelationalOperator: '<S329>/Compare' incorporates:
     *  Constant: '<S329>/Constant'
     */
    VCM20_B.Compare_m = (VCM20_B.SFunction1_o9 == VCM20_P.Constant_Value_dl);

    /* Switch: '<S328>/Switch' */
    if (VCM20_B.Compare_m) {
      /* Switch: '<S328>/Switch' incorporates:
       *  Constant: '<S328>/Constant'
       */
      VCM20_B.Switch_ahz = VCM20_P.Constant_Value_c;
    } else {
      /* Switch: '<S328>/Switch' */
      VCM20_B.Switch_ahz = VCM20_B.SFunction1_o9;
    }

    /* End of Switch: '<S328>/Switch' */

    /* Product: '<S323>/Divide' */
    VCM20_B.FR_Sliprate = VCM20_B.Add_k / VCM20_B.Switch_ahz;

    /* Gain: '<S37>/Gain5' */
    VCM20_B.Gain5 = VCM20_P.Gain5_Gain_l * VCM20_B.FR_Sliprate;

    /* Gain: '<S324>/toTire_rpm' */
    VCM20_B.toTire_rpm_f = VCM20_P.toTire_rpm_Gain_m * VCM20_B.VELOCITY_X_p;

    /* Gain: '<S324>/Gear' */
    VCM20_B.Gear_m = VCM20_P.Gear_Gain_j * VCM20_B.toTire_rpm_f;

    /* Sum: '<S324>/Add' */
    VCM20_B.Add_n = VCM20_B.SFunction1_o9_j - VCM20_B.Gear_m;

    /* RelationalOperator: '<S331>/Compare' incorporates:
     *  Constant: '<S331>/Constant'
     */
    VCM20_B.Compare_pp = (VCM20_B.SFunction1_o9_j == VCM20_P.Constant_Value_i);

    /* Switch: '<S330>/Switch' */
    if (VCM20_B.Compare_pp) {
      /* Switch: '<S330>/Switch' incorporates:
       *  Constant: '<S330>/Constant'
       */
      VCM20_B.Switch_il = VCM20_P.Constant_Value_d;
    } else {
      /* Switch: '<S330>/Switch' */
      VCM20_B.Switch_il = VCM20_B.SFunction1_o9_j;
    }

    /* End of Switch: '<S330>/Switch' */

    /* Product: '<S324>/Divide' */
    VCM20_B.RL_Sliprate = VCM20_B.Add_n / VCM20_B.Switch_il;

    /* Gain: '<S37>/Gain6' */
    VCM20_B.Gain6_l = VCM20_P.Gain6_Gain_b * VCM20_B.RL_Sliprate;

    /* Gain: '<S325>/toTire_rpm' */
    VCM20_B.toTire_rpm_a5 = VCM20_P.toTire_rpm_Gain_h * VCM20_B.VELOCITY_X_p;

    /* Gain: '<S325>/Gear' */
    VCM20_B.Gear_dp = VCM20_P.Gear_Gain_f * VCM20_B.toTire_rpm_a5;

    /* Sum: '<S325>/Add' */
    VCM20_B.Add_ji = VCM20_B.SFunction1_o9_n - VCM20_B.Gear_dp;

    /* RelationalOperator: '<S333>/Compare' incorporates:
     *  Constant: '<S333>/Constant'
     */
    VCM20_B.Compare_fx = (VCM20_B.SFunction1_o9_n == VCM20_P.Constant_Value_in);

    /* Switch: '<S332>/Switch' */
    if (VCM20_B.Compare_fx) {
      /* Switch: '<S332>/Switch' incorporates:
       *  Constant: '<S332>/Constant'
       */
      VCM20_B.Switch_np = VCM20_P.Constant_Value_j;
    } else {
      /* Switch: '<S332>/Switch' */
      VCM20_B.Switch_np = VCM20_B.SFunction1_o9_n;
    }

    /* End of Switch: '<S332>/Switch' */

    /* Product: '<S325>/Divide' */
    VCM20_B.RR_Sliprate = VCM20_B.Add_ji / VCM20_B.Switch_np;

    /* Gain: '<S37>/Gain7' */
    VCM20_B.Gain7 = VCM20_P.Gain7_Gain * VCM20_B.RR_Sliprate;

    /* S-Function (rti_commonblock): '<S59>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s004" Id:1314 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X522]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X522]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X522]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X522]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "FL_Sliprate" (0|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain4_j ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain4_j ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain4_j ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "FR_Sliprate" (16|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain5 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain5 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain5 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "RL_Sliprate" (32|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain6_l ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain6_l ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain6_l ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "RR_Sliprate" (48|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain7 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain7 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain7 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X522], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* DataTypeConversion: '<S37>/Cast3' incorporates:
     *  Constant: '<Root>/TCOn'
     */
    VCM20_B.Cast3 = VCM20_P.TCOn_Value;

    /* Gain: '<S37>/Gain22' */
    VCM20_B.Gain22 = VCM20_P.Gain22_Gain * VCM20_B.Cast3;

    /* DataTypeConversion: '<S37>/Cast4' incorporates:
     *  Constant: '<Root>/ABSOn'
     */
    VCM20_B.Cast4 = VCM20_P.ABSOn_Value;

    /* Gain: '<S37>/Gain23' */
    VCM20_B.Gain23 = VCM20_P.Gain23_Gain * VCM20_B.Cast4;

    /* DataTypeConversion: '<S37>/Cast5' incorporates:
     *  Constant: '<Root>/DYCOn'
     */
    VCM20_B.Cast5 = VCM20_P.DYCOn_Value;

    /* Gain: '<S37>/Gain24' */
    VCM20_B.Gain24 = VCM20_P.Gain24_Gain * VCM20_B.Cast5;

    /* DataTypeConversion: '<S37>/Cast6' incorporates:
     *  Constant: '<Root>/80kWOn'
     */
    VCM20_B.Cast6 = VCM20_P.u0kWOn_Value;

    /* Gain: '<S37>/Gain25' */
    VCM20_B.Gain25 = VCM20_P.Gain25_Gain * VCM20_B.Cast6;
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* Logic: '<S11>/VDCSW' */
    VCM20_B.VDCSW = true;

    /* Logic: '<S320>/Logical Operator1' incorporates:
     *  Constant: '<Root>/VDC_On'
     */
    VCM20_B.VDCStatus = (VCM20_B.VDCSW && (VCM20_P.VDC_On_Value != 0.0));

    /* DataTypeConversion: '<S37>/Cast' */
    VCM20_B.Cast = VCM20_B.VDCStatus;

    /* Gain: '<S37>/Gain18' */
    VCM20_B.Gain18 = VCM20_P.Gain18_Gain * VCM20_B.Cast;

    /* DataTypeConversion: '<S37>/Cast1' */
    VCM20_B.Cast1 = VCM20_B.LogicalOperator2;

    /* Gain: '<S37>/Gain20' */
    VCM20_B.Gain20 = VCM20_P.Gain20_Gain * VCM20_B.Cast1;

    /* Logic: '<S67>/NOT' */
    VCM20_B.NOT = !VCM20_B.AND;

    /* DataTypeConversion: '<S37>/Cast2' */
    VCM20_B.Cast2 = VCM20_B.NOT;

    /* Gain: '<S37>/Gain21' */
    VCM20_B.Gain21 = VCM20_P.Gain21_Gain * VCM20_B.Cast2;

    /* S-Function (rti_commonblock): '<S60>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s004" Id:1793 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X701]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X701]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X701]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X701]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "VDCstatus" (0|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.Gain18 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "BrakeOverRide" (1|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.Gain20 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 1;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "APPSError" (2|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.Gain21 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 2;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "TC" (3|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.Gain22 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "ABS" (4|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.Gain23 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 4;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "DYC" (5|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.Gain24 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 5;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "80kWlimit" (6|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.Gain25 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 6;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X701], 8,
                       &(CAN_Msg[0]), delayTime);
    }
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S49>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "TX Message" Id:1281 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X501]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X501]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X501]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X501]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "Actual Speed RL" (0|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.SFunction1_o9_j ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o9_j ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o9_j ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "Actual Speed FR" (16|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.SFunction1_o9 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o9 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o9 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "Actual Speed FL" (32|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.SFunction1_o9_g ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o9_g ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o9_g ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "Actual Speed RR" (48|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.SFunction1_o9_n ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o9_n ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o9_n ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X501], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* Gain: '<S36>/Gain19' */
    VCM20_B.Gain19 = VCM20_P.Gain19_Gain * VCM20_B.SFunction1_o11_n;

    /* Gain: '<S36>/Gain26' */
    VCM20_B.Gain26 = VCM20_P.Gain26_Gain * VCM20_B.SFunction1_o11;

    /* Gain: '<S36>/Gain27' */
    VCM20_B.Gain27 = VCM20_P.Gain27_Gain * VCM20_B.SFunction1_o11_p;

    /* Gain: '<S36>/Gain28' */
    VCM20_B.Gain28 = VCM20_P.Gain28_Gain * VCM20_B.SFunction1_o11_g;

    /* S-Function (rti_commonblock): '<S50>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s002" Id:1282 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X502]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X502]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X502]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X502]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "Actual torque RL" (0|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain19 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain19 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain19 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "Actual torque FR" (16|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain26 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain26 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain26 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "Actual torque FL" (32|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain27 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain27 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain27 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "Actual torque RR" (48|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain28 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain28 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain28 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X502], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* S-Function (rti_commonblock): '<S51>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s005" Id:1285 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X505]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X505]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X505]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X505]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "bSystemReady RL" (8|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o1_b ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bError RL" (9|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o2_e ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 1;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bWarn RL" (10|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o3_o ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 2;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bQuitDcOn RL" (11|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o4_g ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 3;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bDcOn RL" (12|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o5_p ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 4;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bQuitInverterOn RL" (13|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o6_h ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 5;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bInverterOn RL" (14|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o7_n ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 6;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bDerating RL" (15|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o8_g ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 7;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bSystemReady FR" (24|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o1_m ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bError FR" (25|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o2_c ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 1;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bWarn FR" (26|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o3_i ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 2;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bQuitDcOn FR" (27|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o4_d ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 3;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bDcOn FR" (28|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o5_g ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 4;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bQuitInverterOn FR" (29|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o6_f ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 5;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bInverterOn FR" (30|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o7 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 6;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bDerating FR" (31|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o8 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 7;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bSystemReady FL" (40|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o1_g ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bError FL" (41|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o2_f ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 1;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bWarn FL" (42|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o3_n ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 2;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bQuitDcOn FL" (43|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o4_n ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 3;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bDcOn FL" (44|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o5_j ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 4;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bQuitInverterOn FL" (45|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o6_j ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 5;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bInverterOn FL" (46|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o7_i ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 6;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bDerating FL" (47|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o8_p ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 7;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bSystemReady RR" (56|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o1_e ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bError RR" (57|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o2_b ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 1;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bWarn RR" (58|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o3_k ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 2;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bQuitDcOn RR" (59|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o4_j ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 3;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bDcOn RR" (60|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o5_e ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 4;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bQuitInverterOn RR" (61|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o6_n ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 5;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bInverterOn RR" (62|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o7_j ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 6;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "bDerating RR" (63|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o8_o ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 7;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte0;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X505], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* S-Function (rti_commonblock): '<S28>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "MC_32" Id:649 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X289]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X289]->processed) {
        can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X289]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X289]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X289]->timestamp > 0.0) {
        if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X289]->processed) {
          CAN_Msg = can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X289]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "MC3_Motor temperature" (0|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o1_h = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_Inverter cold plate temperature" (16|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o2_d = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_IGBT temperature" (32|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o3_o3 = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_Diagnostic number" (48|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[7];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o4_dd = ((real_T) CAN_Sgn.UnsignedSgn);
          }
        }
      }
    }

    /* S-Function (rti_commonblock): '<S22>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "MC_12" Id:645 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X285]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X285]->processed) {
        can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X285]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X285]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X285]->timestamp > 0.0) {
        if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X285]->processed) {
          CAN_Msg = can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X285]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "MC1_Motor temperature" (0|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o1_ml = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_Inverter cold plate temperature" (16|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o2_o = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_IGBT temperature" (32|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o3_e = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_Diagnostic number" (48|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[7];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o4_c = ((real_T) CAN_Sgn.UnsignedSgn);
          }
        }
      }
    }

    /* S-Function (rti_commonblock): '<S31>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "MC_42" Id:650 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X28A]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X28A]->processed) {
        can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X28A]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X28A]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X28A]->timestamp > 0.0) {
        if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X28A]->processed) {
          CAN_Msg = can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X28A]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "MC4_Motor temperature" (0|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o1_i = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_Inverter cold plate temperature" (16|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o2_n = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_IGBT temperature" (32|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o3_f = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_Diagnostic number" (48|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[7];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o4_e = ((real_T) CAN_Sgn.UnsignedSgn);
          }
        }
      }
    }

    /* S-Function (rti_commonblock): '<S25>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "MC_22" Id:646 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X286]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X286]->processed) {
        can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X286]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X286]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X286]->timestamp > 0.0) {
        if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X286]->processed) {
          CAN_Msg = can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X286]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "MC2_Motor temperature" (0|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o1_eq = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_Inverter cold plate temperature" (16|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o2_oc = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_IGBT temperature" (32|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o3_d2j = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_Diagnostic number" (48|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[7];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o4_f = ((real_T) CAN_Sgn.UnsignedSgn);
          }
        }
      }
    }

    /* S-Function (rti_commonblock): '<S52>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s006" Id:1286 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X506]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X506]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X506]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X506]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "Motor Temp RL" (0|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o1_h ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "Motor Temp FR" (16|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o1_ml ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "Motor Temp FL" (32|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o1_i ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "Motor Temp RR" (48|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o1_eq ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X506], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* S-Function (rti_commonblock): '<S53>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s007" Id:1287 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X507]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X507]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X507]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X507]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "IGBT Temp RL" (0|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o3_o3 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "IGBT Temp FR" (16|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o3_e ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "IGBT Temp FL" (32|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o3_f ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "IGBT Temp RR" (48|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o3_d2j ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X507], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* S-Function (rti_commonblock): '<S54>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s008" Id:1288 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X508]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X508]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X508]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X508]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "Diag Number RL" (0|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o4_dd ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "Diag Number FR" (16|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o4_c ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "Diag Number FL" (32|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o4_e ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "Diag Number RR" (48|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o4_f ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X508], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* S-Function (rti_commonblock): '<S29>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "MC_33" Id:773 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X305]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X305]->processed) {
        can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X305]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X305]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X305]->timestamp > 0.0) {
        if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X305]->processed) {
          CAN_Msg = can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X305]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "MC3_Error info 1" (0|8, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.UnsignedSgn &= 0x000000FF;
            VCM20_B.SFunction1_o1_j = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_Error info 2" (8|8, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x000000FF;
            VCM20_B.SFunction1_o2_p = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_Error info 3" (16|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o3_h = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC3_Actual speed value in 0;0001/rpm" (32|32, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[7];
            VCM20_B.SFunction1_o4_hn = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* S-Function (rti_commonblock): '<S23>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "MC_13" Id:769 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X301]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X301]->processed) {
        can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X301]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X301]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X301]->timestamp > 0.0) {
        if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X301]->processed) {
          CAN_Msg = can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X301]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "MC1_Error info 1" (0|8, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.UnsignedSgn &= 0x000000FF;
            VCM20_B.SFunction1_o1_c = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_Error info 2" (8|8, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x000000FF;
            VCM20_B.SFunction1_o2_an = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_Error info 3" (16|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o3_g = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC1_Actual speed value in 0;0001/rpm" (32|32, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[7];
            VCM20_B.SFunction1_o4_jb = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* MinMax: '<S36>/Min' */
    u0_0 = VCM20_B.SFunction1_o2_o;
    u1_0 = VCM20_B.SFunction1_o2_oc;
    if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    u1 = VCM20_B.SFunction1_o2_d;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.SFunction1_o2_n;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    /* MinMax: '<S36>/Min' */
    VCM20_B.Min_a = u1_0;

    /* S-Function (rti_commonblock): '<S55>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s009" Id:1289 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X509]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X509]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X509]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X509]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "Info1 RL" (0|8, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o1_j ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Info2 RL" (8|8, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o2_p ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Info3 RL" (16|8, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o3_h ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Info1 FR" (24|8, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o1_c ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Info2 FR" (32|8, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o2_an ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Info3 FR" (40|8, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o3_g ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "ColdPlateTemp" (48|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.Min_a ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X509], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* S-Function (rti_commonblock): '<S32>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "MC_43" Id:774 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X306]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X306]->processed) {
        can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X306]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X306]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X306]->timestamp > 0.0) {
        if (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X306]->processed) {
          CAN_Msg = can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X306]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "MC4_Error info 1" (0|8, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.UnsignedSgn &= 0x000000FF;
            VCM20_B.SFunction1_o1_f = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_Error info 2" (8|8, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x000000FF;
            VCM20_B.SFunction1_o2_o4 = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_Error info 3" (16|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o3_io = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC4_Actual speed value in 0;0001/rpm" (32|32, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[7];
            VCM20_B.SFunction1_o4_ea = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* S-Function (rti_commonblock): '<S26>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "MC_23" Id:770 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X302]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X302]->processed) {
        can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X302]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X302]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X302]->timestamp > 0.0) {
        if (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X302]->processed) {
          CAN_Msg = can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X302]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "MC2_Error info 1" (0|8, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.UnsignedSgn &= 0x000000FF;
            VCM20_B.SFunction1_o1_om = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_Error info 2" (8|8, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x000000FF;
            VCM20_B.SFunction1_o2_ao = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_Error info 3" (16|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o3_iq = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "MC2_Actual speed value in 0;0001/rpm" (32|32, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[7];
            VCM20_B.SFunction1_o4_p = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* MinMax: '<S36>/Min1' */
    u0_0 = VCM20_B.SFunction1_o10;
    u1_0 = VCM20_B.SFunction1_o10_ny;
    if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    u1 = VCM20_B.SFunction1_o10_n;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.SFunction1_o10_l;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    /* MinMax: '<S36>/Min1' */
    VCM20_B.Min1 = u1_0;

    /* S-Function (rti_commonblock): '<S56>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s00A" Id:1290 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50A]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50A]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50A]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50A]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "Info1 FL" (0|8, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o1_f ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Info2 FL" (8|8, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o2_o4 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Info3 FL" (16|8, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o3_io ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Info1 RR" (24|8, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o1_om ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Info2 RR" (32|8, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o2_ao ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Info3 RR" (40|8, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SFunction1_o3_iq ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "DCBusVoltage" (48|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.Min1 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50A], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* S-Function (rti_commonblock): '<S105>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:289 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X121]);

      /* Convert timestamp */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X121]->processed) {
        can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X121]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X121]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X121]->timestamp > 0.0) {
        if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X121]->processed) {
          CAN_Msg = can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X121]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "ACCEL_X" (0|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o1_m3 = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "ACCEL_Y" (16|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o2_i = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "ACCEL_Z" (32|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o3_nw = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* Gain: '<S14>/ACCEL_X_Gain' */
    VCM20_B.ACCEL_X = VCM20_P.ACCEL_X_Gain_Gain * VCM20_B.SFunction1_o1_m3;

    /* Gain: '<S14>/Gain' */
    VCM20_B.ACCEL_X_m = VCM20_P.Gain_Gain_dv * VCM20_B.ACCEL_X;

    /* Gain: '<S35>/Gain3' */
    VCM20_B.ACCEL_X_mq = VCM20_P.Gain3_Gain_l * VCM20_B.ACCEL_X_m;

    /* Gain: '<S14>/ACCEL_Y_Gain' */
    VCM20_B.ACCEL_Y = VCM20_P.ACCEL_Y_Gain_Gain * VCM20_B.SFunction1_o2_i;

    /* Gain: '<S35>/Gain5' */
    VCM20_B.ACCEL_Y_b = VCM20_P.Gain5_Gain_o * VCM20_B.ACCEL_Y;

    /* Gain: '<S14>/ACCEL_Z_Gain' */
    VCM20_B.ACCEL_Z = VCM20_P.ACCEL_Z_Gain_Gain * VCM20_B.SFunction1_o3_nw;

    /* Gain: '<S35>/Gain4' */
    VCM20_B.ACCEL_Y_e = VCM20_P.Gain4_Gain_p * VCM20_B.ACCEL_Z;

    /* S-Function (rti_commonblock): '<S43>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "TX Message" Id:1291 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50B]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50B]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50B]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50B]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "ACCEL X" (0|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.ACCEL_X_mq ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.ACCEL_X_mq ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.ACCEL_X_mq ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "ACCEL Y" (16|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.ACCEL_Y_b ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.ACCEL_Y_b ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.ACCEL_Y_b ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "ACCEL Z" (32|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.ACCEL_Y_e ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.ACCEL_Y_e ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.ACCEL_Y_e ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50B], 6,
                       &(CAN_Msg[0]), delayTime);
    }

    /* S-Function (rti_commonblock): '<S106>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:290 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X122]);

      /* Convert timestamp */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X122]->processed) {
        can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X122]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X122]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X122]->timestamp > 0.0) {
        if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X122]->processed) {
          CAN_Msg = can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X122]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "GYRO_X" (0|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o1_bx = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "GYRO_Y" (16|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o2_bg = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "GYRO_Z" (32|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o3_b = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* Gain: '<S14>/Gain3' */
    VCM20_B.GYRO_X = VCM20_P.Gain3_Gain_p * VCM20_B.SFunction1_o1_bx;

    /* Gain: '<S35>/Gain6' */
    VCM20_B.GYRO_X_o = VCM20_P.Gain6_Gain_p * VCM20_B.GYRO_X;

    /* Gain: '<S14>/Gain5' */
    VCM20_B.GYRO_Y = VCM20_P.Gain5_Gain_p * VCM20_B.SFunction1_o2_bg;

    /* Gain: '<S35>/Gain8' */
    VCM20_B.GYRO_Y_o = VCM20_P.Gain8_Gain_f * VCM20_B.GYRO_Y;

    /* Gain: '<S14>/Gain4' */
    VCM20_B.GYRO_Z = VCM20_P.Gain4_Gain_ns * VCM20_B.SFunction1_o3_b;

    /* Gain: '<S35>/Gain7' */
    VCM20_B.GYRO_Z_m = VCM20_P.Gain7_Gain_m * VCM20_B.GYRO_Z;

    /* S-Function (rti_commonblock): '<S44>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "TX Message" Id:1292 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50C]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50C]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50C]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50C]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "GYRO X" (0|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.GYRO_X_o ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.GYRO_X_o ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.GYRO_X_o ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "GYRO Y" (16|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.GYRO_Y_o ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.GYRO_Y_o ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.GYRO_Y_o ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "GYRO Z" (32|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.GYRO_Z_m ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.GYRO_Z_m ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.GYRO_Z_m ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50C], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* Gain: '<S35>/Gain9' */
    VCM20_B.VELOCITY_X_p5 = VCM20_P.Gain9_Gain * VCM20_B.VELOCITY_X_p;

    /* Gain: '<S14>/Gain10' */
    VCM20_B.VELOCITY_Y = VCM20_P.Gain10_Gain_o * VCM20_B.SFunction1_o2;

    /* Gain: '<S35>/Gain10' */
    VCM20_B.VELOCITY_Y_j = VCM20_P.Gain10_Gain_m * VCM20_B.VELOCITY_Y;

    /* Gain: '<S14>/Gain9' */
    VCM20_B.VELOCITY_Z = VCM20_P.Gain9_Gain_o * VCM20_B.SFunction1_o3;

    /* Gain: '<S35>/Gain14' */
    VCM20_B.VELOCITY_Z_a = VCM20_P.Gain14_Gain_e * VCM20_B.VELOCITY_Z;

    /* S-Function (rti_commonblock): '<S45>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "TX Message" Id:1293 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50D]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50D]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50D]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50D]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "VELOCITY X" (0|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.VELOCITY_X_p5 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.VELOCITY_X_p5 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.VELOCITY_X_p5 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "VELOCITY Y" (16|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.VELOCITY_Y_j ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.VELOCITY_Y_j ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.VELOCITY_Y_j ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "VELOCITY Z" (32|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.VELOCITY_Z_a ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.VELOCITY_Z_a ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.VELOCITY_Z_a ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50D], 6,
                       &(CAN_Msg[0]), delayTime);
    }

    /* S-Function (rti_commonblock): '<S111>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:544 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X220]);

      /* Convert timestamp */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X220]->processed) {
        can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X220]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X220]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X220]->timestamp > 0.0) {
        if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X220]->processed) {
          CAN_Msg = can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X220]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "ANGLE_TRACK" (0|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o1_hc = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "ANGLE_SLIP" (16|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o2_nj = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "CURVATURE_RADIUS" (32|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o3_dd = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* Gain: '<S14>/Gain24' */
    VCM20_B.ANGLE_TRACK = VCM20_P.Gain24_Gain_g * VCM20_B.SFunction1_o1_hc;

    /* Gain: '<S35>/Gain11' */
    VCM20_B.ANGLE_TRACK_g = VCM20_P.Gain11_Gain_e * VCM20_B.ANGLE_TRACK;

    /* Gain: '<S14>/Gain23' */
    VCM20_B.ANGLE_SLIP = VCM20_P.Gain23_Gain_p * VCM20_B.SFunction1_o2_nj;

    /* Gain: '<S35>/Gain13' */
    VCM20_B.ANGLE_SLIP_d = VCM20_P.Gain13_Gain_a * VCM20_B.ANGLE_SLIP;

    /* Gain: '<S14>/Gain25' */
    VCM20_B.CURVATURE_RADIUS = VCM20_P.Gain25_Gain_a * VCM20_B.SFunction1_o3_dd;

    /* Gain: '<S35>/Gain12' */
    VCM20_B.CURVATURE_RADIUS_a = VCM20_P.Gain12_Gain_m *
      VCM20_B.CURVATURE_RADIUS;

    /* S-Function (rti_commonblock): '<S46>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "TX Message" Id:1294 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50E]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50E]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50E]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50E]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "ANGLE_TRACK" (0|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.ANGLE_TRACK_g ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.ANGLE_TRACK_g ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.ANGLE_TRACK_g ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "ANGLE_SLIP" (16|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.ANGLE_SLIP_d ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.ANGLE_SLIP_d ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.ANGLE_SLIP_d ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "CURVATURE_RADIUS" (32|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.CURVATURE_RADIUS_a ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.CURVATURE_RADIUS_a ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.CURVATURE_RADIUS_a ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50E], 6,
                       &(CAN_Msg[0]), delayTime);
    }

    /* S-Function (rti_commonblock): '<S107>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:308 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X134]);

      /* Convert timestamp */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X134]->processed) {
        can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X134]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X134]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X134]->timestamp > 0.0) {
        if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X134]->processed) {
          CAN_Msg = can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X134]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "LATITUDE" (0|32, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[3];
            VCM20_B.SFunction1_o1_nu = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "LONGITUDE" (32|32, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[7];
            VCM20_B.SFunction1_o2_ej = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* Gain: '<S14>/Gain21' */
    VCM20_B.LATITUDE = VCM20_P.Gain21_Gain_n * VCM20_B.SFunction1_o1_nu;

    /* Gain: '<S35>/Gain16' */
    VCM20_B.LATITUDE_c = VCM20_P.Gain16_Gain * VCM20_B.LATITUDE;

    /* Gain: '<S14>/Gain22' */
    VCM20_B.LONGITUDE = VCM20_P.Gain22_Gain_d * VCM20_B.SFunction1_o2_ej;

    /* Gain: '<S35>/Gain15' */
    VCM20_B.LONGITUDE_a = VCM20_P.Gain15_Gain * VCM20_B.LONGITUDE;

    /* S-Function (rti_commonblock): '<S47>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "TX Message" Id:1295 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50F]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50F]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50F]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50F]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "LATITUDE" (0|32, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.LATITUDE_c ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.LATITUDE_c ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.LATITUDE_c ) + 0.5);

        {
          UInt32 tmp;
          tmp = CAN_Sgn.SgnBytes.Byte0;
          CAN_Sgn.SgnBytes.Byte0 = CAN_Sgn.SgnBytes.Byte3;
          CAN_Sgn.SgnBytes.Byte3 = tmp;
          tmp = CAN_Sgn.SgnBytes.Byte1;
          CAN_Sgn.SgnBytes.Byte1 = CAN_Sgn.SgnBytes.Byte2;
          CAN_Sgn.SgnBytes.Byte2 = tmp;
        }

        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte1;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LONGITUDE" (32|32, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.LONGITUDE_a ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.LONGITUDE_a ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.LONGITUDE_a ) + 0.5);

        {
          UInt32 tmp;
          tmp = CAN_Sgn.SgnBytes.Byte0;
          CAN_Sgn.SgnBytes.Byte0 = CAN_Sgn.SgnBytes.Byte3;
          CAN_Sgn.SgnBytes.Byte3 = tmp;
          tmp = CAN_Sgn.SgnBytes.Byte1;
          CAN_Sgn.SgnBytes.Byte1 = CAN_Sgn.SgnBytes.Byte2;
          CAN_Sgn.SgnBytes.Byte2 = tmp;
        }

        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte2;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte1;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte0;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50F], 8,
                       &(CAN_Msg[0]), delayTime);
    }
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[3] == 0) {
    /* S-Function (rti_commonblock): '<S104>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:258 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X102]);

      /* Convert timestamp */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X102]->processed) {
        can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X102]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X102]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X102]->timestamp > 0.0) {
        if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X102]->processed) {
          CAN_Msg = can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X102]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "SOLUTION_MODE" (0|4, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.UnsignedSgn &= 0x0000000F;
            VCM20_B.SFunction1_o1_or = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "VELOCITY_VALID" (6|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 6;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o2_l = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "POSITION_VALID" (7|1, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) >> 7;
            CAN_Sgn.UnsignedSgn &= 0x00000001;
            VCM20_B.SFunction1_o3_nj = ((real_T) CAN_Sgn.UnsignedSgn);
          }
        }
      }
    }

    /* Gain: '<S14>/Gain18' */
    VCM20_B.SOLUTION_MODE = VCM20_P.Gain18_Gain_n * VCM20_B.SFunction1_o1_or;

    /* Gain: '<S35>/Gain2' */
    VCM20_B.SOLUTION_MODE_g = VCM20_P.Gain2_Gain_h * VCM20_B.SOLUTION_MODE;

    /* Gain: '<S14>/Gain19' */
    VCM20_B.VELOCITY_VALID = VCM20_P.Gain19_Gain_f * VCM20_B.SFunction1_o2_l;

    /* Gain: '<S35>/Gain1' */
    VCM20_B.VELOCITY_VALID_b = VCM20_P.Gain1_Gain_el * VCM20_B.VELOCITY_VALID;

    /* Gain: '<S14>/Gain20' */
    VCM20_B.POSITION_VALID = VCM20_P.Gain20_Gain_m * VCM20_B.SFunction1_o3_nj;

    /* Gain: '<S35>/Gain18' */
    VCM20_B.POSITION_VALID_a = VCM20_P.Gain18_Gain_l * VCM20_B.POSITION_VALID;

    /* S-Function (rti_commonblock): '<S48>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "TX Message" Id:1296 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X510]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X510]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X510]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X510]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "SOLUTION_MODE" (0|4, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.SOLUTION_MODE_g ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x0000000F;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "VELOCITY_VALID" (6|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.VELOCITY_VALID_b ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 6;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "POSITION_VALID" (7|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.POSITION_VALID_a ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 7;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X510], 6,
                       &(CAN_Msg[0]), delayTime);
    }
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* Gain: '<S34>/Gain' */
    VCM20_B.Gain_n = VCM20_P.Gain_Gain_ic * VCM20_B.Saturation_d;

    /* Gain: '<S34>/Gain2' */
    VCM20_B.Gain2_l = VCM20_P.Gain2_Gain_m * VCM20_B.SFunction1_o2_a;

    /* Gain: '<S34>/Gain1' */
    VCM20_B.Gain1_j = VCM20_P.Gain1_Gain_g * VCM20_B.SFunction1_o1_l;
  }

  /* Gain: '<S34>/Gain29' */
  VCM20_B.Gain29 = VCM20_P.Gain29_Gain * VCM20_B.SteerAngle_rad;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S41>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s003" Id:1283 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X503]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X503]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X503]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X503]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "throttle" (0|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.Gain_n ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "Steering angle" (16|16, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Gain29 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain29 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Gain29 ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.SignedSgn &= 0xFFFF0000;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "BrakePressF" (32|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.Gain2_l ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte2;

        /* ...... "BrakePressR" (48|16, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.Gain1_j ) + 0.5);
        CAN_Sgn.SgnBytes.Byte3 = CAN_Sgn.SgnBytes.Byte0;
        CAN_Sgn.SgnBytes.Byte2 = CAN_Sgn.SgnBytes.Byte1;
        CAN_Sgn.UnsignedSgn &= 0xFFFF0000;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte3;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte2;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X503], 8,
                       &(CAN_Msg[0]), delayTime);
    }
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* DataTypeConversion: '<S34>/Data Type Conversion8' */
    VCM20_B.DataTypeConversion8 = VCM20_B.BrakeSW;

    /* S-Function (rti_commonblock): '<S42>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s003" Id:1797 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X705]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X705]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X705]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X705]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "BrakeSW" (0|1, standard signal, unsigned int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion8 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X705], 8,
                       &(CAN_Msg[0]), delayTime);
    }
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S20>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:1025 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X401]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X401]->processed) {
        can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X401]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X401]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X401]->timestamp > 0.0) {
        if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X401]->processed) {
          CAN_Msg = can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X401]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "TS_Voltage" (0|8, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SignedSgn &= 0x000000FF;
            if (CAN_Sgn.SignedSgn >> 7) {
              CAN_Sgn.SignedSgn |= 0xFFFFFF00;
            }

            VCM20_B.SFunction1_o1_lm = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "MaxCellVoltage" (8|8, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[1];
            CAN_Sgn.SignedSgn &= 0x000000FF;
            if (CAN_Sgn.SignedSgn >> 7) {
              CAN_Sgn.SignedSgn |= 0xFFFFFF00;
            }

            VCM20_B.SFunction1_o2_lh = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "MinCellVoltage" (16|8, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SignedSgn &= 0x000000FF;
            if (CAN_Sgn.SignedSgn >> 7) {
              CAN_Sgn.SignedSgn |= 0xFFFFFF00;
            }

            VCM20_B.SFunction1_o3_ol = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "MaxCellTemperture" (24|8, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x000000FF;
            if (CAN_Sgn.SignedSgn >> 7) {
              CAN_Sgn.SignedSgn |= 0xFFFFFF00;
            }

            VCM20_B.SFunction1_o4_i = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "MinCellTemperture" (32|8, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SignedSgn &= 0x000000FF;
            if (CAN_Sgn.SignedSgn >> 7) {
              CAN_Sgn.SignedSgn |= 0xFFFFFF00;
            }

            VCM20_B.SFunction1_o5_h = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "Error_info" (40|8, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[5];
            CAN_Sgn.SignedSgn &= 0x000000FF;
            if (CAN_Sgn.SignedSgn >> 7) {
              CAN_Sgn.SignedSgn |= 0xFFFFFF00;
            }

            VCM20_B.SFunction1_o6_d = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "Error_Count" (48|8, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[6];
            CAN_Sgn.SignedSgn &= 0x000000FF;
            if (CAN_Sgn.SignedSgn >> 7) {
              CAN_Sgn.SignedSgn |= 0xFFFFFF00;
            }

            VCM20_B.SFunction1_o7_k = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "AMS_Status" (56|8, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[7];
            CAN_Sgn.SignedSgn &= 0x000000FF;
            if (CAN_Sgn.SignedSgn >> 7) {
              CAN_Sgn.SignedSgn |= 0xFFFFFF00;
            }

            VCM20_B.SFunction1_o8_go = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* S-Function (rti_commonblock): '<S40>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s005" Id:1794 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X702]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X702]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X702]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X702]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "TS_Voltage" (0|8, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.SFunction1_o1_lm ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o1_lm ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o1_lm ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x000000FF;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "MaxCellVoltage" (8|8, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.SFunction1_o2_lh ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o2_lh ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o2_lh ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x000000FF;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "MinCellVoltage" (16|8, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.SFunction1_o3_ol ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o3_ol ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o3_ol ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x000000FF;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "MaxCellTemperture" (24|8, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.SFunction1_o4_i ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o4_i ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o4_i ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x000000FF;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "MinCellTemperture" (32|8, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.SFunction1_o5_h ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o5_h ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o5_h ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x000000FF;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Erroe_info" (40|8, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.SFunction1_o6_d ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o6_d ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o6_d ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x000000FF;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Error_Count" (48|8, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.SFunction1_o7_k ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o7_k ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o7_k ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x000000FF;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "AMS_Status" (56|8, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.SFunction1_o8_go ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o8_go ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.SFunction1_o8_go ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x000000FF;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte0;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X702], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* S-Function (rti_commonblock): '<S342>/S-Function1' */
    /* This comment workarounds a code generation problem */
    {
      /* dSPACE I/O Board DS1401STDADCT4 #1 Unit:ADC Group:ADC */
      adc_tp4_single_new_read(ADC_TP4_1_MODULE_ADDR,
        ADC_TP4_CH9,
        (dsfloat *)&VCM20_B.SFunction1_m);
    }

    /* Product: '<S339>/Product' incorporates:
     *  Constant: '<S339>/Constant'
     */
    VCM20_B.Product_g = VCM20_P.Constant_Value_b2 * VCM20_B.SFunction1_m;

    /* Sum: '<S339>/Add' incorporates:
     *  Constant: '<S339>/Constant1'
     */
    VCM20_B.reikyaku_temp2 = VCM20_B.Product_g - VCM20_P.Constant1_Value_f;

    /* S-Function (rti_commonblock): '<S343>/S-Function1' */
    /* This comment workarounds a code generation problem */
    {
      /* dSPACE I/O Board DS1401STDADCT4 #1 Unit:ADC Group:ADC */
      adc_tp4_single_new_read(ADC_TP4_1_MODULE_ADDR,
        ADC_TP4_CH10,
        (dsfloat *)&VCM20_B.SFunction1_c);
    }

    /* Product: '<S340>/Product' incorporates:
     *  Constant: '<S340>/Constant'
     */
    VCM20_B.Product_b = VCM20_P.Constant_Value_mx * VCM20_B.SFunction1_c;

    /* Sum: '<S340>/Add' incorporates:
     *  Constant: '<S340>/Constant1'
     */
    VCM20_B.reikyaku_temp3 = VCM20_B.Product_b - VCM20_P.Constant1_Value_l;

    /* S-Function (rti_commonblock): '<S61>/S-Function1' incorporates:
     *  Constant: '<S38>/Constant'
     */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s005" Id:1795 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X703]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X703]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X703]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X703]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "coolor_temp1" (0|8, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_P.Constant_Value_it ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_P.Constant_Value_it ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_P.Constant_Value_it ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x000000FF;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "coolor_temp2" (8|8, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.reikyaku_temp2 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.reikyaku_temp2 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.reikyaku_temp2 ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x000000FF;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "coolor_temp3" (16|8, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.reikyaku_temp3 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.reikyaku_temp3 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.reikyaku_temp3 ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x000000FF;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte0;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X703], 8,
                       &(CAN_Msg[0]), delayTime);
    }
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[3] == 0) {
    /* S-Function (rti_commonblock): '<S95>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* DataTypeConversion: '<S39>/Data Type Conversion10' */
    VCM20_B.DataTypeConversion10 = VCM20_B.SFunction1_b;

    /* S-Function (rti_commonblock): '<S96>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* DataTypeConversion: '<S39>/Data Type Conversion1' */
    VCM20_B.DataTypeConversion1 = VCM20_B.SFunction1_o;

    /* S-Function (rti_commonblock): '<S97>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* DataTypeConversion: '<S39>/Data Type Conversion2' */
    VCM20_B.DataTypeConversion2 = VCM20_B.SFunction1_j;

    /* S-Function (rti_commonblock): '<S62>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "EVO4s005" Id:1796 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X704]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X704]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X704]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X704]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "AMS_shutdownsig" (0|1, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.DataTypeConversion10 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.DataTypeConversion10 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.DataTypeConversion10 ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x00000001;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "BSPD_shutdownsig" (1|1, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.DataTypeConversion1 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.DataTypeConversion1 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.DataTypeConversion1 ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x00000001;
        CAN_Sgn.SignedSgn = ((UInt32)CAN_Sgn.SignedSgn) << 1;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "IMD_shutdownsig" (2|1, standard signal, signed int, big endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.DataTypeConversion2 ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.DataTypeConversion2 ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.DataTypeConversion2 ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x00000001;
        CAN_Sgn.SignedSgn = ((UInt32)CAN_Sgn.SignedSgn) << 2;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X704], 8,
                       &(CAN_Msg[0]), delayTime);
    }
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Switch: '<S133>/Switch1' incorporates:
     *  Constant: '<Root>/MC1_sw'
     */
    if (VCM20_P.MC1_sw_Value > VCM20_P.Switch1_Threshold_n) {
      /* Switch: '<S133>/Switch1' */
      VCM20_B.Switch1_l = VCM20_B.SFunction1_o4_d;
    } else {
      /* Switch: '<S133>/Switch1' incorporates:
       *  Constant: '<S133>/Constant'
       */
      VCM20_B.Switch1_l = VCM20_P.Constant_Value_n4;
    }

    /* End of Switch: '<S133>/Switch1' */

    /* Switch: '<S133>/Switch2' incorporates:
     *  Constant: '<Root>/MC2_sw'
     */
    if (VCM20_P.MC2_sw_Value > VCM20_P.Switch2_Threshold_h) {
      /* Switch: '<S133>/Switch2' */
      VCM20_B.Switch2_eg = VCM20_B.SFunction1_o4_j;
    } else {
      /* Switch: '<S133>/Switch2' incorporates:
       *  Constant: '<S133>/Constant'
       */
      VCM20_B.Switch2_eg = VCM20_P.Constant_Value_n4;
    }

    /* End of Switch: '<S133>/Switch2' */

    /* Switch: '<S133>/Switch3' incorporates:
     *  Constant: '<Root>/MC3_sw'
     */
    if (VCM20_P.MC3_sw_Value > VCM20_P.Switch3_Threshold_dy) {
      /* Switch: '<S133>/Switch3' */
      VCM20_B.Switch3_a = VCM20_B.SFunction1_o4_g;
    } else {
      /* Switch: '<S133>/Switch3' incorporates:
       *  Constant: '<S133>/Constant'
       */
      VCM20_B.Switch3_a = VCM20_P.Constant_Value_n4;
    }

    /* End of Switch: '<S133>/Switch3' */

    /* Switch: '<S133>/Switch4' incorporates:
     *  Constant: '<Root>/MC4_sw'
     */
    if (VCM20_P.MC4_sw_Value > VCM20_P.Switch4_Threshold_n) {
      /* Switch: '<S133>/Switch4' */
      VCM20_B.Switch4_i = VCM20_B.SFunction1_o4_n;
    } else {
      /* Switch: '<S133>/Switch4' incorporates:
       *  Constant: '<S133>/Constant'
       */
      VCM20_B.Switch4_i = VCM20_P.Constant_Value_n4;
    }

    /* End of Switch: '<S133>/Switch4' */

    /* Logic: '<S116>/All_QuitDcOn' */
    VCM20_B.All_QuitDcOn = ((VCM20_B.Switch1_l != 0.0) && (VCM20_B.Switch2_eg !=
      0.0) && (VCM20_B.Switch3_a != 0.0) && (VCM20_B.Switch4_i != 0.0));
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[3] == 0) {
    /* Logic: '<S11>/RtoDSW' */
    VCM20_B.RtoDSW = !(VCM20_B.SFunction1_o4 != 0.0);

    /* Logic: '<S11>/ErrorRessetSW' */
    VCM20_B.ErrorResetSW = !(VCM20_B.SFunction1_o6 != 0.0);

    /* Logic: '<S11>/ErrorResetSW1' incorporates:
     *  Constant: '<S11>/ErrorReset_ControllDesk'
     */
    VCM20_B.ErrorResetSW_f = (VCM20_B.ErrorResetSW ||
      (VCM20_P.ErrorReset_ControllDesk_Value != 0.0));

    /* DataTypeConversion: '<S114>/Data Type Conversion' */
    VCM20_B.AMK_bErrorReset = VCM20_B.ErrorResetSW_f;
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Logic: '<S116>/ReadyToDrive' */
    VCM20_B.ReadyToDrive = (VCM20_B.All_QuitDcOn && VCM20_B.BrakeSW &&
      VCM20_B.RtoDSW);

    /* Logic: '<S116>/ReadyToDrive_or' */
    VCM20_B.ReadyToDrive_or = (VCM20_B.ReadyToDrive || VCM20_B.All_enable);

    /* Switch: '<S131>/Switch1' incorporates:
     *  Constant: '<Root>/MC1_sw'
     */
    if (VCM20_P.MC1_sw_Value > VCM20_P.Switch1_Threshold_dd) {
      /* Switch: '<S131>/Switch1' */
      VCM20_B.Switch1_a = VCM20_B.SFunction1_o2_c;
    } else {
      /* Switch: '<S131>/Switch1' incorporates:
       *  Constant: '<S131>/Constant'
       */
      VCM20_B.Switch1_a = VCM20_P.Constant_Value_i4;
    }

    /* End of Switch: '<S131>/Switch1' */

    /* Switch: '<S131>/Switch2' incorporates:
     *  Constant: '<Root>/MC2_sw'
     */
    if (VCM20_P.MC2_sw_Value > VCM20_P.Switch2_Threshold_m1) {
      /* Switch: '<S131>/Switch2' */
      VCM20_B.Switch2_fm = VCM20_B.SFunction1_o2_b;
    } else {
      /* Switch: '<S131>/Switch2' incorporates:
       *  Constant: '<S131>/Constant'
       */
      VCM20_B.Switch2_fm = VCM20_P.Constant_Value_i4;
    }

    /* End of Switch: '<S131>/Switch2' */

    /* Switch: '<S131>/Switch3' incorporates:
     *  Constant: '<Root>/MC3_sw'
     */
    if (VCM20_P.MC3_sw_Value > VCM20_P.Switch3_Threshold_o) {
      /* Switch: '<S131>/Switch3' */
      VCM20_B.Switch3_d = VCM20_B.SFunction1_o2_e;
    } else {
      /* Switch: '<S131>/Switch3' incorporates:
       *  Constant: '<S131>/Constant'
       */
      VCM20_B.Switch3_d = VCM20_P.Constant_Value_i4;
    }

    /* End of Switch: '<S131>/Switch3' */

    /* Switch: '<S131>/Switch4' incorporates:
     *  Constant: '<Root>/MC4_sw'
     */
    if (VCM20_P.MC4_sw_Value > VCM20_P.Switch4_Threshold_g) {
      /* Switch: '<S131>/Switch4' */
      VCM20_B.Switch4_g = VCM20_B.SFunction1_o2_f;
    } else {
      /* Switch: '<S131>/Switch4' incorporates:
       *  Constant: '<S131>/Constant'
       */
      VCM20_B.Switch4_g = VCM20_P.Constant_Value_i4;
    }

    /* End of Switch: '<S131>/Switch4' */

    /* Logic: '<S116>/All_Error' */
    VCM20_B.All_Error = ((!(VCM20_B.Switch1_a != 0.0)) && (!(VCM20_B.Switch2_fm
      != 0.0)) && (!(VCM20_B.Switch3_d != 0.0)) && (!(VCM20_B.Switch4_g != 0.0)));

    /* Logic: '<S116>/ReadyToDrive_and' */
    VCM20_B.ReadyToDrive_and = (VCM20_B.All_QuitDcOn && VCM20_B.ReadyToDrive_or &&
      VCM20_B.All_Error);

    /* Switch: '<S124>/Switch1' incorporates:
     *  Constant: '<Root>/MC1_sw'
     */
    if (VCM20_P.MC1_sw_Value > VCM20_P.Switch1_Threshold_pk) {
      /* Switch: '<S124>/Switch1' */
      VCM20_B.Switch1_h = VCM20_B.ReadyToDrive_and;
    } else {
      /* Switch: '<S124>/Switch1' incorporates:
       *  Constant: '<S124>/Constant'
       */
      VCM20_B.Switch1_h = VCM20_P.Constant_Value_ou;
    }

    /* End of Switch: '<S124>/Switch1' */

    /* DataTypeConversion: '<S5>/Data Type Conversion' */
    VCM20_B.DataTypeConversion = VCM20_B.Switch1_h;

    /* Switch: '<S134>/Switch1' incorporates:
     *  Constant: '<Root>/MC1_sw'
     */
    if (VCM20_P.MC1_sw_Value > VCM20_P.Switch1_Threshold_nh) {
      /* Switch: '<S134>/Switch1' */
      VCM20_B.Switch1_dh = VCM20_B.SFunction1_o1_m;
    } else {
      /* Switch: '<S134>/Switch1' incorporates:
       *  Constant: '<S134>/Constant'
       */
      VCM20_B.Switch1_dh = VCM20_P.Constant_Value_cc;
    }

    /* End of Switch: '<S134>/Switch1' */

    /* Switch: '<S134>/Switch2' incorporates:
     *  Constant: '<Root>/MC2_sw'
     */
    if (VCM20_P.MC2_sw_Value > VCM20_P.Switch2_Threshold_gb) {
      /* Switch: '<S134>/Switch2' */
      VCM20_B.Switch2_fz = VCM20_B.SFunction1_o1_e;
    } else {
      /* Switch: '<S134>/Switch2' incorporates:
       *  Constant: '<S134>/Constant'
       */
      VCM20_B.Switch2_fz = VCM20_P.Constant_Value_cc;
    }

    /* End of Switch: '<S134>/Switch2' */

    /* Switch: '<S134>/Switch3' incorporates:
     *  Constant: '<Root>/MC3_sw'
     */
    if (VCM20_P.MC3_sw_Value > VCM20_P.Switch3_Threshold_j) {
      /* Switch: '<S134>/Switch3' */
      VCM20_B.Switch3_ad = VCM20_B.SFunction1_o1_b;
    } else {
      /* Switch: '<S134>/Switch3' incorporates:
       *  Constant: '<S134>/Constant'
       */
      VCM20_B.Switch3_ad = VCM20_P.Constant_Value_cc;
    }

    /* End of Switch: '<S134>/Switch3' */

    /* Switch: '<S134>/Switch4' incorporates:
     *  Constant: '<Root>/MC4_sw'
     */
    if (VCM20_P.MC4_sw_Value > VCM20_P.Switch4_Threshold_b) {
      /* Switch: '<S134>/Switch4' */
      VCM20_B.Switch4_ie = VCM20_B.SFunction1_o1_g;
    } else {
      /* Switch: '<S134>/Switch4' incorporates:
       *  Constant: '<S134>/Constant'
       */
      VCM20_B.Switch4_ie = VCM20_P.Constant_Value_cc;
    }

    /* End of Switch: '<S134>/Switch4' */

    /* Logic: '<S116>/All_sbm' */
    VCM20_B.All_sbm = ((VCM20_B.Switch1_dh != 0.0) && (VCM20_B.Switch2_fz != 0.0)
                       && (VCM20_B.Switch3_ad != 0.0) && (VCM20_B.Switch4_ie !=
      0.0));

    /* MinMax: '<S116>/Min DC bus Voltage' */
    u0_0 = VCM20_B.SFunction1_o10;
    u1_0 = VCM20_B.SFunction1_o10_ny;
    if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    u1 = VCM20_B.SFunction1_o10_n;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.SFunction1_o10_l;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    /* MinMax: '<S116>/Min DC bus Voltage' */
    VCM20_B.MinDCbusVoltage = u1_0;

    /* RelationalOperator: '<S130>/Compare' incorporates:
     *  Constant: '<S130>/Constant'
     */
    VCM20_B.Compare_jr = (VCM20_B.MinDCbusVoltage >=
                          VCM20_P.CompareToConstant_const_l);

    /* Logic: '<S116>/DcOn' */
    VCM20_B.DcOn = (VCM20_B.All_sbm && VCM20_B.Compare_jr);

    /* DataTypeConversion: '<S5>/Data Type Conversion2' */
    VCM20_B.DataTypeConversion2_g = VCM20_B.DcOn;

    /* Switch: '<S123>/Switch1' incorporates:
     *  Constant: '<Root>/MC1_sw'
     */
    if (VCM20_P.MC1_sw_Value > VCM20_P.Switch1_Threshold_g) {
      /* Switch: '<S123>/Switch1' */
      VCM20_B.Switch1_m = VCM20_B.ReadyToDrive_and;
    } else {
      /* Switch: '<S123>/Switch1' incorporates:
       *  Constant: '<S123>/Constant'
       */
      VCM20_B.Switch1_m = VCM20_P.Constant_Value_bn;
    }

    /* End of Switch: '<S123>/Switch1' */

    /* DataTypeConversion: '<S5>/Data Type Conversion1' */
    VCM20_B.DataTypeConversion1_i = VCM20_B.Switch1_m;

    /* Gain: '<S122>/Gain8' */
    VCM20_B.FRspeed = VCM20_P.Gain8_Gain_g * VCM20_B.Saturation2;

    /* MinMax: '<S125>/Min' */
    u0_0 = VCM20_B.FRspeed;
    u1_0 = VCM20_B.Saturation2;
    if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S125>/Min' */
    VCM20_B.Min_h = u1_0;

    /* Switch: '<S136>/Switch' incorporates:
     *  Concatenate: '<S136>/Vector Concatenate'
     *  Constant: '<Root>/80kWOn'
     *  Constant: '<S136>/all power limit'
     *  Constant: '<S136>/each power limit '
     */
    if (VCM20_P.u0kWOn_Value > VCM20_P.Switch_Threshold_a) {
      VCM20_B.VectorConcatenate[0] = VCM20_P.allpowerlimit_Value;
    } else {
      VCM20_B.VectorConcatenate[0] = VCM20_P.eachpowerlimit_Value;
    }

    /* End of Switch: '<S136>/Switch' */

    /* Saturate: '<S136>/Saturation' incorporates:
     *  Concatenate: '<S136>/Vector Concatenate'
     */
    u0_0 = VCM20_B.SFunction1_o9;
    u1 = VCM20_P.Saturation_LowerSat_eu;
    u2 = VCM20_P.Saturation_UpperSat_a4;
    if (u0_0 > u2) {
      u0_0 = u2;
    } else if (u0_0 < u1) {
      u0_0 = u1;
    }

    VCM20_B.VectorConcatenate[1] = u0_0;

    /* End of Saturate: '<S136>/Saturation' */

    /* Fcn: '<S136>/Fcn' */
    VCM20_B.Fcn = VCM20_B.VectorConcatenate[0] / (6.2831853071795862 *
      VCM20_B.VectorConcatenate[1] / 60.0);
  }

  /* Switch: '<S125>/MC1' incorporates:
   *  Constant: '<Root>/MC1_sw'
   */
  if (VCM20_P.MC1_sw_Value > VCM20_P.MC1_Threshold) {
    /* Switch: '<S125>/MC1' */
    VCM20_B.MC1 = VCM20_B.FRtrq;
  } else {
    /* Switch: '<S125>/MC1' incorporates:
     *  Constant: '<S125>/Constant'
     */
    VCM20_B.MC1 = VCM20_P.Constant_Value_cj;
  }

  /* End of Switch: '<S125>/MC1' */

  /* RelationalOperator: '<S142>/LowerRelop1' */
  VCM20_B.LowerRelop1_i = (VCM20_B.MC1 > VCM20_B.Fcn);

  /* RelationalOperator: '<S142>/UpperRelop' incorporates:
   *  Constant: '<S136>/Constant1'
   */
  VCM20_B.UpperRelop_pg = (VCM20_B.MC1 < VCM20_P.Constant1_Value_h);

  /* Switch: '<S142>/Switch' */
  if (VCM20_B.UpperRelop_pg) {
    /* Switch: '<S142>/Switch' incorporates:
     *  Constant: '<S136>/Constant1'
     */
    VCM20_B.Switch_nf = VCM20_P.Constant1_Value_h;
  } else {
    /* Switch: '<S142>/Switch' */
    VCM20_B.Switch_nf = VCM20_B.MC1;
  }

  /* End of Switch: '<S142>/Switch' */

  /* Switch: '<S142>/Switch2' */
  if (VCM20_B.LowerRelop1_i) {
    /* Switch: '<S142>/Switch2' */
    VCM20_B.Switch2_n = VCM20_B.Fcn;
  } else {
    /* Switch: '<S142>/Switch2' */
    VCM20_B.Switch2_n = VCM20_B.Switch_nf;
  }

  /* End of Switch: '<S142>/Switch2' */

  /* Gain: '<S117>/Gain2' */
  VCM20_B.Gain2_n = VCM20_P.Gain2_Gain_oe * VCM20_B.Switch2_n;

  /* Quantizer: '<S117>/Quantizer1' */
  u1_0 = VCM20_B.Gain2_n;

  /* Quantizer: '<S117>/Quantizer1' */
  VCM20_B.Quantizer1 = rt_roundd_snf(u1_0 / VCM20_P.Quantizer1_Interval) *
    VCM20_P.Quantizer1_Interval;

  /* Switch: '<S138>/Switch' */
  if (VCM20_B.Quantizer1 > VCM20_P.Switch_Threshold_j) {
    /* Switch: '<S138>/Switch' */
    VCM20_B.Switch_as = VCM20_B.Min_h;
  } else {
    /* Switch: '<S138>/Switch' incorporates:
     *  Constant: '<S138>/Constant1'
     */
    VCM20_B.Switch_as = VCM20_P.Constant1_Value_a;
  }

  /* End of Switch: '<S138>/Switch' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Delay: '<S138>/Delay' */
    VCM20_B.Delay_g = VCM20_DW.Delay_DSTATE_h;
  }

  /* MinMax: '<S138>/MinMax' */
  u0_0 = VCM20_B.Switch_as;
  u1_0 = VCM20_B.Delay_g;
  if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S138>/MinMax' */
  VCM20_B.MinMax = u1_0;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Delay: '<S140>/Delay' */
    VCM20_B.Delay_o = VCM20_DW.Delay_DSTATE_k;
  }

  /* MinMax: '<S140>/MinMax' */
  u0_0 = VCM20_B.Quantizer1;
  u1_0 = VCM20_B.Delay_o;
  if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S140>/MinMax' */
  VCM20_B.MinMax_c = u1_0;

  /* MinMax: '<S137>/MinMax' incorporates:
   *  Constant: '<S137>/M set point MAX 0.1%Mn'
   */
  u0_0 = VCM20_B.MinMax_c;
  u1_0 = VCM20_P.MsetpointMAX01Mn_Value;
  if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S137>/MinMax' */
  VCM20_B.Positivetorquelimit = u1_0;

  /* RelationalOperator: '<S141>/LowerRelop1' incorporates:
   *  Constant: '<S135>/power limit'
   */
  VCM20_B.LowerRelop1_nq = (VCM20_B.MC1 > VCM20_P.powerlimit_Value);
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S135>/Constant' incorporates:
     *  Concatenate: '<S135>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_a[0] = VCM20_P.Constant_Value_gc;

    /* Saturate: '<S135>/Saturation' incorporates:
     *  Concatenate: '<S135>/Vector Concatenate'
     */
    u0_0 = VCM20_B.SFunction1_o9;
    u1 = VCM20_P.Saturation_LowerSat_pt;
    u2 = VCM20_P.Saturation_UpperSat_lu;
    if (u0_0 > u2) {
      u0_0 = u2;
    } else if (u0_0 < u1) {
      u0_0 = u1;
    }

    VCM20_B.VectorConcatenate_a[1] = u0_0;

    /* End of Saturate: '<S135>/Saturation' */

    /* Fcn: '<S135>/Fcn' */
    VCM20_B.Fcn_p = VCM20_B.VectorConcatenate_a[0] / (6.2831853071795862 *
      VCM20_B.VectorConcatenate_a[1] / 60.0);
  }

  /* RelationalOperator: '<S141>/UpperRelop' */
  VCM20_B.UpperRelop_db = (VCM20_B.MC1 < VCM20_B.Fcn_p);

  /* Switch: '<S141>/Switch' */
  if (VCM20_B.UpperRelop_db) {
    /* Switch: '<S141>/Switch' */
    VCM20_B.Switch_e3 = VCM20_B.Fcn_p;
  } else {
    /* Switch: '<S141>/Switch' */
    VCM20_B.Switch_e3 = VCM20_B.MC1;
  }

  /* End of Switch: '<S141>/Switch' */

  /* Switch: '<S141>/Switch2' */
  if (VCM20_B.LowerRelop1_nq) {
    /* Switch: '<S141>/Switch2' incorporates:
     *  Constant: '<S135>/power limit'
     */
    VCM20_B.Switch2_di = VCM20_P.powerlimit_Value;
  } else {
    /* Switch: '<S141>/Switch2' */
    VCM20_B.Switch2_di = VCM20_B.Switch_e3;
  }

  /* End of Switch: '<S141>/Switch2' */

  /* Gain: '<S117>/Gain6' */
  VCM20_B.Gain6_lv = VCM20_P.Gain6_Gain_k * VCM20_B.Switch2_di;

  /* Quantizer: '<S117>/Quantizer2' */
  u1_0 = VCM20_B.Gain6_lv;

  /* Quantizer: '<S117>/Quantizer2' */
  VCM20_B.Quantizer2 = rt_roundd_snf(u1_0 / VCM20_P.Quantizer2_Interval) *
    VCM20_P.Quantizer2_Interval;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Delay: '<S139>/Delay' */
    VCM20_B.Delay_l = VCM20_DW.Delay_DSTATE_b;
  }

  /* MinMax: '<S139>/MinMax' */
  u0_0 = VCM20_B.Quantizer2;
  u1_0 = VCM20_B.Delay_l;
  if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S139>/MinMax' */
  VCM20_B.MinMax_co = u1_0;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S63>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "VCM_11" Id:388 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X184]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X184]->processed) {
        can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X184]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X184]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "AMK_bInverterOn" (8|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "AMK_bDcOn" (9|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion2_g ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 1;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "AMK_bEnable" (10|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion1_i ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 2;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "AMK_bErrorReset" (11|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.AMK_bErrorReset ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 3;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Speed Setpoint in 1/rpm" (16|16, standard signal, signed int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.MinMax ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x0000FFFF;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte1;

        /* ...... "Positive toruque limit in 0,1 % Mn" (32|16, standard signal, signed int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Positivetorquelimit ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Positivetorquelimit ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Positivetorquelimit ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x0000FFFF;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte1;

        /* ...... "Negative torque limit in 0,1% Mn" (48|16, standard signal, signed int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.MinMax_co ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_co ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_co ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x0000FFFF;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte1;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X184], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* Switch: '<S124>/Switch2' incorporates:
     *  Constant: '<Root>/MC2_sw'
     */
    if (VCM20_P.MC2_sw_Value > VCM20_P.Switch2_Threshold_f) {
      /* Switch: '<S124>/Switch2' */
      VCM20_B.Switch2_aw = VCM20_B.ReadyToDrive_and;
    } else {
      /* Switch: '<S124>/Switch2' incorporates:
       *  Constant: '<S124>/Constant'
       */
      VCM20_B.Switch2_aw = VCM20_P.Constant_Value_ou;
    }

    /* End of Switch: '<S124>/Switch2' */

    /* DataTypeConversion: '<S5>/Data Type Conversion3' */
    VCM20_B.DataTypeConversion3 = VCM20_B.Switch2_aw;

    /* DataTypeConversion: '<S5>/Data Type Conversion5' */
    VCM20_B.DataTypeConversion5 = VCM20_B.DcOn;

    /* Switch: '<S123>/Switch2' incorporates:
     *  Constant: '<Root>/MC2_sw'
     */
    if (VCM20_P.MC2_sw_Value > VCM20_P.Switch2_Threshold_e) {
      /* Switch: '<S123>/Switch2' */
      VCM20_B.Switch2_g1 = VCM20_B.ReadyToDrive_and;
    } else {
      /* Switch: '<S123>/Switch2' incorporates:
       *  Constant: '<S123>/Constant'
       */
      VCM20_B.Switch2_g1 = VCM20_P.Constant_Value_bn;
    }

    /* End of Switch: '<S123>/Switch2' */

    /* DataTypeConversion: '<S5>/Data Type Conversion4' */
    VCM20_B.DataTypeConversion4 = VCM20_B.Switch2_g1;

    /* Gain: '<S122>/Gain6' */
    VCM20_B.RRspeed = VCM20_P.Gain6_Gain_n * VCM20_B.Saturation2;

    /* MinMax: '<S125>/Min1' */
    u0_0 = VCM20_B.RRspeed;
    u1_0 = VCM20_B.Saturation2;
    if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S125>/Min1' */
    VCM20_B.Min1_i = u1_0;

    /* Switch: '<S144>/Switch' incorporates:
     *  Concatenate: '<S144>/Vector Concatenate'
     *  Constant: '<Root>/80kWOn'
     *  Constant: '<S144>/all power limit'
     *  Constant: '<S144>/each power limit '
     */
    if (VCM20_P.u0kWOn_Value > VCM20_P.Switch_Threshold_lh) {
      VCM20_B.VectorConcatenate_b[0] = VCM20_P.allpowerlimit_Value_g;
    } else {
      VCM20_B.VectorConcatenate_b[0] = VCM20_P.eachpowerlimit_Value_b;
    }

    /* End of Switch: '<S144>/Switch' */

    /* Saturate: '<S144>/Saturation' incorporates:
     *  Concatenate: '<S144>/Vector Concatenate'
     */
    u0_0 = VCM20_B.SFunction1_o9_n;
    u1 = VCM20_P.Saturation_LowerSat_h3;
    u2 = VCM20_P.Saturation_UpperSat_lm;
    if (u0_0 > u2) {
      u0_0 = u2;
    } else if (u0_0 < u1) {
      u0_0 = u1;
    }

    VCM20_B.VectorConcatenate_b[1] = u0_0;

    /* End of Saturate: '<S144>/Saturation' */

    /* Fcn: '<S144>/Fcn' */
    VCM20_B.Fcn_l = VCM20_B.VectorConcatenate_b[0] / (6.2831853071795862 *
      VCM20_B.VectorConcatenate_b[1] / 60.0);
  }

  /* Switch: '<S125>/MC2' incorporates:
   *  Constant: '<Root>/MC2_sw'
   */
  if (VCM20_P.MC2_sw_Value > VCM20_P.MC2_Threshold) {
    /* Switch: '<S125>/MC2' */
    VCM20_B.MC2 = VCM20_B.RRtrq;
  } else {
    /* Switch: '<S125>/MC2' incorporates:
     *  Constant: '<S125>/Constant'
     */
    VCM20_B.MC2 = VCM20_P.Constant_Value_cj;
  }

  /* End of Switch: '<S125>/MC2' */

  /* RelationalOperator: '<S150>/LowerRelop1' */
  VCM20_B.LowerRelop1_mv = (VCM20_B.MC2 > VCM20_B.Fcn_l);

  /* RelationalOperator: '<S150>/UpperRelop' incorporates:
   *  Constant: '<S144>/Constant1'
   */
  VCM20_B.UpperRelop_pv = (VCM20_B.MC2 < VCM20_P.Constant1_Value_ie);

  /* Switch: '<S150>/Switch' */
  if (VCM20_B.UpperRelop_pv) {
    /* Switch: '<S150>/Switch' incorporates:
     *  Constant: '<S144>/Constant1'
     */
    VCM20_B.Switch_l = VCM20_P.Constant1_Value_ie;
  } else {
    /* Switch: '<S150>/Switch' */
    VCM20_B.Switch_l = VCM20_B.MC2;
  }

  /* End of Switch: '<S150>/Switch' */

  /* Switch: '<S150>/Switch2' */
  if (VCM20_B.LowerRelop1_mv) {
    /* Switch: '<S150>/Switch2' */
    VCM20_B.Switch2_pb = VCM20_B.Fcn_l;
  } else {
    /* Switch: '<S150>/Switch2' */
    VCM20_B.Switch2_pb = VCM20_B.Switch_l;
  }

  /* End of Switch: '<S150>/Switch2' */

  /* Gain: '<S118>/Gain2' */
  VCM20_B.Gain2_bd = VCM20_P.Gain2_Gain_l * VCM20_B.Switch2_pb;

  /* Quantizer: '<S118>/Quantizer1' */
  u1_0 = VCM20_B.Gain2_bd;

  /* Quantizer: '<S118>/Quantizer1' */
  VCM20_B.Quantizer1_p = rt_roundd_snf(u1_0 / VCM20_P.Quantizer1_Interval_d) *
    VCM20_P.Quantizer1_Interval_d;

  /* Switch: '<S146>/Switch' */
  if (VCM20_B.Quantizer1_p > VCM20_P.Switch_Threshold_c) {
    /* Switch: '<S146>/Switch' */
    VCM20_B.Switch_kd = VCM20_B.Min1_i;
  } else {
    /* Switch: '<S146>/Switch' incorporates:
     *  Constant: '<S146>/Constant1'
     */
    VCM20_B.Switch_kd = VCM20_P.Constant1_Value_is;
  }

  /* End of Switch: '<S146>/Switch' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Delay: '<S146>/Delay' */
    VCM20_B.Delay_ox = VCM20_DW.Delay_DSTATE_e;
  }

  /* MinMax: '<S146>/MinMax' */
  u0_0 = VCM20_B.Switch_kd;
  u1_0 = VCM20_B.Delay_ox;
  if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S146>/MinMax' */
  VCM20_B.MinMax_g = u1_0;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Delay: '<S148>/Delay' */
    VCM20_B.Delay_gx = VCM20_DW.Delay_DSTATE_p;
  }

  /* MinMax: '<S148>/MinMax' */
  u0_0 = VCM20_B.Quantizer1_p;
  u1_0 = VCM20_B.Delay_gx;
  if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S148>/MinMax' */
  VCM20_B.MinMax_i = u1_0;

  /* MinMax: '<S145>/MinMax' incorporates:
   *  Constant: '<S145>/M set point MAX 0.1%Mn'
   */
  u0_0 = VCM20_B.MinMax_i;
  u1_0 = VCM20_P.MsetpointMAX01Mn_Value_p;
  if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S145>/MinMax' */
  VCM20_B.Positivetorquelimit_j = u1_0;

  /* RelationalOperator: '<S149>/LowerRelop1' incorporates:
   *  Constant: '<S143>/power limit'
   */
  VCM20_B.LowerRelop1_g3 = (VCM20_B.MC2 > VCM20_P.powerlimit_Value_a);
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S143>/Constant' incorporates:
     *  Concatenate: '<S143>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_p[0] = VCM20_P.Constant_Value_i0;

    /* Saturate: '<S143>/Saturation' incorporates:
     *  Concatenate: '<S143>/Vector Concatenate'
     */
    u0_0 = VCM20_B.SFunction1_o9_n;
    u1 = VCM20_P.Saturation_LowerSat_np;
    u2 = VCM20_P.Saturation_UpperSat_eqy;
    if (u0_0 > u2) {
      u0_0 = u2;
    } else if (u0_0 < u1) {
      u0_0 = u1;
    }

    VCM20_B.VectorConcatenate_p[1] = u0_0;

    /* End of Saturate: '<S143>/Saturation' */

    /* Fcn: '<S143>/Fcn' */
    VCM20_B.Fcn_k = VCM20_B.VectorConcatenate_p[0] / (6.2831853071795862 *
      VCM20_B.VectorConcatenate_p[1] / 60.0);
  }

  /* RelationalOperator: '<S149>/UpperRelop' */
  VCM20_B.UpperRelop_cy = (VCM20_B.MC2 < VCM20_B.Fcn_k);

  /* Switch: '<S149>/Switch' */
  if (VCM20_B.UpperRelop_cy) {
    /* Switch: '<S149>/Switch' */
    VCM20_B.Switch_bm = VCM20_B.Fcn_k;
  } else {
    /* Switch: '<S149>/Switch' */
    VCM20_B.Switch_bm = VCM20_B.MC2;
  }

  /* End of Switch: '<S149>/Switch' */

  /* Switch: '<S149>/Switch2' */
  if (VCM20_B.LowerRelop1_g3) {
    /* Switch: '<S149>/Switch2' incorporates:
     *  Constant: '<S143>/power limit'
     */
    VCM20_B.Switch2_pv = VCM20_P.powerlimit_Value_a;
  } else {
    /* Switch: '<S149>/Switch2' */
    VCM20_B.Switch2_pv = VCM20_B.Switch_bm;
  }

  /* End of Switch: '<S149>/Switch2' */

  /* Gain: '<S118>/Gain6' */
  VCM20_B.Gain6_e = VCM20_P.Gain6_Gain_e * VCM20_B.Switch2_pv;

  /* Quantizer: '<S118>/Quantizer2' */
  u1_0 = VCM20_B.Gain6_e;

  /* Quantizer: '<S118>/Quantizer2' */
  VCM20_B.Quantizer2_p = rt_roundd_snf(u1_0 / VCM20_P.Quantizer2_Interval_p) *
    VCM20_P.Quantizer2_Interval_p;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Delay: '<S147>/Delay' */
    VCM20_B.Delay_j = VCM20_DW.Delay_DSTATE_i;
  }

  /* MinMax: '<S147>/MinMax' */
  u0_0 = VCM20_B.Quantizer2_p;
  u1_0 = VCM20_B.Delay_j;
  if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S147>/MinMax' */
  VCM20_B.MinMax_m = u1_0;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S64>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "VCM_21" Id:389 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X185]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X185]->processed) {
        can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X185]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X185]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "AMK_bInverterOn" (8|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion3 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "AMK_bDcOn" (9|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion5 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 1;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "AMK_bEnable" (10|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion4 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 2;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "AMK_bErrorReset" (11|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.AMK_bErrorReset ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 3;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Speed Setpoint in 1/rpm" (16|16, standard signal, signed int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.MinMax_g ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_g ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_g ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x0000FFFF;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte1;

        /* ...... "Positive toruque limit in 0,1 % Mn" (32|16, standard signal, signed int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Positivetorquelimit_j ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Positivetorquelimit_j ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Positivetorquelimit_j ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x0000FFFF;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte1;

        /* ...... "Negative torque limit in 0,1% Mn" (48|16, standard signal, signed int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.MinMax_m ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_m ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_m ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x0000FFFF;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte1;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X185], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* Switch: '<S124>/Switch3' incorporates:
     *  Constant: '<Root>/MC3_sw'
     */
    if (VCM20_P.MC3_sw_Value > VCM20_P.Switch3_Threshold_p) {
      /* Switch: '<S124>/Switch3' */
      VCM20_B.Switch3_j = VCM20_B.ReadyToDrive_and;
    } else {
      /* Switch: '<S124>/Switch3' incorporates:
       *  Constant: '<S124>/Constant'
       */
      VCM20_B.Switch3_j = VCM20_P.Constant_Value_ou;
    }

    /* End of Switch: '<S124>/Switch3' */

    /* DataTypeConversion: '<S5>/Data Type Conversion6' */
    VCM20_B.DataTypeConversion6 = VCM20_B.Switch3_j;

    /* DataTypeConversion: '<S5>/Data Type Conversion7' */
    VCM20_B.DataTypeConversion7 = VCM20_B.DcOn;

    /* Switch: '<S123>/Switch3' incorporates:
     *  Constant: '<Root>/MC3_sw'
     */
    if (VCM20_P.MC3_sw_Value > VCM20_P.Switch3_Threshold_f) {
      /* Switch: '<S123>/Switch3' */
      VCM20_B.Switch3_n = VCM20_B.ReadyToDrive_and;
    } else {
      /* Switch: '<S123>/Switch3' incorporates:
       *  Constant: '<S123>/Constant'
       */
      VCM20_B.Switch3_n = VCM20_P.Constant_Value_bn;
    }

    /* End of Switch: '<S123>/Switch3' */

    /* DataTypeConversion: '<S5>/Data Type Conversion8' */
    VCM20_B.DataTypeConversion8_l = VCM20_B.Switch3_n;

    /* Gain: '<S122>/Gain7' */
    VCM20_B.RLspeed = VCM20_P.Gain7_Gain_h * VCM20_B.Saturation2;

    /* MinMax: '<S125>/Min2' */
    u0_0 = VCM20_B.RLspeed;
    u1_0 = VCM20_B.Saturation2;
    if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S125>/Min2' */
    VCM20_B.Min2 = u1_0;

    /* Switch: '<S152>/Switch' incorporates:
     *  Concatenate: '<S152>/Vector Concatenate'
     *  Constant: '<Root>/80kWOn'
     *  Constant: '<S152>/all power limit'
     *  Constant: '<S152>/each power limit '
     */
    if (VCM20_P.u0kWOn_Value > VCM20_P.Switch_Threshold_m) {
      VCM20_B.VectorConcatenate_j[0] = VCM20_P.allpowerlimit_Value_i;
    } else {
      VCM20_B.VectorConcatenate_j[0] = VCM20_P.eachpowerlimit_Value_k;
    }

    /* End of Switch: '<S152>/Switch' */

    /* Saturate: '<S152>/Saturation' incorporates:
     *  Concatenate: '<S152>/Vector Concatenate'
     */
    u0_0 = VCM20_B.SFunction1_o9_j;
    u1 = VCM20_P.Saturation_LowerSat_bh;
    u2 = VCM20_P.Saturation_UpperSat_lx;
    if (u0_0 > u2) {
      u0_0 = u2;
    } else if (u0_0 < u1) {
      u0_0 = u1;
    }

    VCM20_B.VectorConcatenate_j[1] = u0_0;

    /* End of Saturate: '<S152>/Saturation' */

    /* Fcn: '<S152>/Fcn' */
    VCM20_B.Fcn_b = VCM20_B.VectorConcatenate_j[0] / (6.2831853071795862 *
      VCM20_B.VectorConcatenate_j[1] / 60.0);
  }

  /* Switch: '<S125>/MC3' incorporates:
   *  Constant: '<Root>/MC3_sw'
   */
  if (VCM20_P.MC3_sw_Value > VCM20_P.MC3_Threshold) {
    /* Switch: '<S125>/MC3' */
    VCM20_B.MC3 = VCM20_B.RLtrq;
  } else {
    /* Switch: '<S125>/MC3' incorporates:
     *  Constant: '<S125>/Constant'
     */
    VCM20_B.MC3 = VCM20_P.Constant_Value_cj;
  }

  /* End of Switch: '<S125>/MC3' */

  /* RelationalOperator: '<S158>/LowerRelop1' */
  VCM20_B.LowerRelop1_j = (VCM20_B.MC3 > VCM20_B.Fcn_b);

  /* RelationalOperator: '<S158>/UpperRelop' incorporates:
   *  Constant: '<S152>/Constant1'
   */
  VCM20_B.UpperRelop_lf = (VCM20_B.MC3 < VCM20_P.Constant1_Value_c);

  /* Switch: '<S158>/Switch' */
  if (VCM20_B.UpperRelop_lf) {
    /* Switch: '<S158>/Switch' incorporates:
     *  Constant: '<S152>/Constant1'
     */
    VCM20_B.Switch_j2 = VCM20_P.Constant1_Value_c;
  } else {
    /* Switch: '<S158>/Switch' */
    VCM20_B.Switch_j2 = VCM20_B.MC3;
  }

  /* End of Switch: '<S158>/Switch' */

  /* Switch: '<S158>/Switch2' */
  if (VCM20_B.LowerRelop1_j) {
    /* Switch: '<S158>/Switch2' */
    VCM20_B.Switch2_oq = VCM20_B.Fcn_b;
  } else {
    /* Switch: '<S158>/Switch2' */
    VCM20_B.Switch2_oq = VCM20_B.Switch_j2;
  }

  /* End of Switch: '<S158>/Switch2' */

  /* Gain: '<S119>/Gain2' */
  VCM20_B.Gain2_h = VCM20_P.Gain2_Gain_e * VCM20_B.Switch2_oq;

  /* Quantizer: '<S119>/Quantizer1' */
  u1_0 = VCM20_B.Gain2_h;

  /* Quantizer: '<S119>/Quantizer1' */
  VCM20_B.Quantizer1_m = rt_roundd_snf(u1_0 / VCM20_P.Quantizer1_Interval_n) *
    VCM20_P.Quantizer1_Interval_n;

  /* Switch: '<S154>/Switch' */
  if (VCM20_B.Quantizer1_m > VCM20_P.Switch_Threshold_fg) {
    /* Switch: '<S154>/Switch' */
    VCM20_B.Switch_k1 = VCM20_B.Min2;
  } else {
    /* Switch: '<S154>/Switch' incorporates:
     *  Constant: '<S154>/Constant1'
     */
    VCM20_B.Switch_k1 = VCM20_P.Constant1_Value_n;
  }

  /* End of Switch: '<S154>/Switch' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Delay: '<S154>/Delay' */
    VCM20_B.Delay_m = VCM20_DW.Delay_DSTATE_bc;
  }

  /* MinMax: '<S154>/MinMax' */
  u0_0 = VCM20_B.Switch_k1;
  u1_0 = VCM20_B.Delay_m;
  if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S154>/MinMax' */
  VCM20_B.MinMax_h = u1_0;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Delay: '<S156>/Delay' */
    VCM20_B.Delay_i = VCM20_DW.Delay_DSTATE_kj;
  }

  /* MinMax: '<S156>/MinMax' */
  u0_0 = VCM20_B.Quantizer1_m;
  u1_0 = VCM20_B.Delay_i;
  if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S156>/MinMax' */
  VCM20_B.MinMax_k = u1_0;

  /* MinMax: '<S153>/MinMax' incorporates:
   *  Constant: '<S153>/M set point MAX 0.1%Mn'
   */
  u0_0 = VCM20_B.MinMax_k;
  u1_0 = VCM20_P.MsetpointMAX01Mn_Value_d;
  if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S153>/MinMax' */
  VCM20_B.Positivetorquelimit_m = u1_0;

  /* RelationalOperator: '<S157>/LowerRelop1' incorporates:
   *  Constant: '<S151>/power limit'
   */
  VCM20_B.LowerRelop1_mh = (VCM20_B.MC3 > VCM20_P.powerlimit_Value_d);
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S151>/Constant' incorporates:
     *  Concatenate: '<S151>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_ju[0] = VCM20_P.Constant_Value_nt;

    /* Saturate: '<S151>/Saturation' incorporates:
     *  Concatenate: '<S151>/Vector Concatenate'
     */
    u0_0 = VCM20_B.SFunction1_o9_j;
    u1 = VCM20_P.Saturation_LowerSat_ms;
    u2 = VCM20_P.Saturation_UpperSat_ep;
    if (u0_0 > u2) {
      u0_0 = u2;
    } else if (u0_0 < u1) {
      u0_0 = u1;
    }

    VCM20_B.VectorConcatenate_ju[1] = u0_0;

    /* End of Saturate: '<S151>/Saturation' */

    /* Fcn: '<S151>/Fcn' */
    VCM20_B.Fcn_c = VCM20_B.VectorConcatenate_ju[0] / (6.2831853071795862 *
      VCM20_B.VectorConcatenate_ju[1] / 60.0);
  }

  /* RelationalOperator: '<S157>/UpperRelop' */
  VCM20_B.UpperRelop_lq = (VCM20_B.MC3 < VCM20_B.Fcn_c);

  /* Switch: '<S157>/Switch' */
  if (VCM20_B.UpperRelop_lq) {
    /* Switch: '<S157>/Switch' */
    VCM20_B.Switch_c1 = VCM20_B.Fcn_c;
  } else {
    /* Switch: '<S157>/Switch' */
    VCM20_B.Switch_c1 = VCM20_B.MC3;
  }

  /* End of Switch: '<S157>/Switch' */

  /* Switch: '<S157>/Switch2' */
  if (VCM20_B.LowerRelop1_mh) {
    /* Switch: '<S157>/Switch2' incorporates:
     *  Constant: '<S151>/power limit'
     */
    VCM20_B.Switch2_cp = VCM20_P.powerlimit_Value_d;
  } else {
    /* Switch: '<S157>/Switch2' */
    VCM20_B.Switch2_cp = VCM20_B.Switch_c1;
  }

  /* End of Switch: '<S157>/Switch2' */

  /* Gain: '<S119>/Gain6' */
  VCM20_B.Gain6_p = VCM20_P.Gain6_Gain_pn * VCM20_B.Switch2_cp;

  /* Quantizer: '<S119>/Quantizer2' */
  u1_0 = VCM20_B.Gain6_p;

  /* Quantizer: '<S119>/Quantizer2' */
  VCM20_B.Quantizer2_n = rt_roundd_snf(u1_0 / VCM20_P.Quantizer2_Interval_k) *
    VCM20_P.Quantizer2_Interval_k;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Delay: '<S155>/Delay' */
    VCM20_B.Delay_m5 = VCM20_DW.Delay_DSTATE_hs;
  }

  /* MinMax: '<S155>/MinMax' */
  u0_0 = VCM20_B.Quantizer2_n;
  u1_0 = VCM20_B.Delay_m5;
  if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S155>/MinMax' */
  VCM20_B.MinMax_d = u1_0;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S65>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "VCM01" Id:392 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X188]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X188]->processed) {
        can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X188]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X188]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "AMK_bInverterOn" (8|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion6 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "AMK_bDcOn" (9|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion7 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 1;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "AMK_bEnable" (10|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion8_l ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 2;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "AMK_bErrorReset" (11|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.AMK_bErrorReset ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 3;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Speed Setpoint in 1/rpm" (16|16, standard signal, signed int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.MinMax_h ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_h ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_h ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x0000FFFF;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte1;

        /* ...... "Positive toruque limit in 0,1 % Mn" (32|16, standard signal, signed int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Positivetorquelimit_m ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Positivetorquelimit_m ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Positivetorquelimit_m ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x0000FFFF;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte1;

        /* ...... "Negative torque limit in 0,1% Mn" (48|16, standard signal, signed int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.MinMax_d ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_d ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_d ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x0000FFFF;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte1;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X188], 8,
                       &(CAN_Msg[0]), delayTime);
    }

    /* Switch: '<S124>/Switch4' incorporates:
     *  Constant: '<Root>/MC4_sw'
     */
    if (VCM20_P.MC4_sw_Value > VCM20_P.Switch4_Threshold_a) {
      /* Switch: '<S124>/Switch4' */
      VCM20_B.Switch4_b = VCM20_B.ReadyToDrive_and;
    } else {
      /* Switch: '<S124>/Switch4' incorporates:
       *  Constant: '<S124>/Constant'
       */
      VCM20_B.Switch4_b = VCM20_P.Constant_Value_ou;
    }

    /* End of Switch: '<S124>/Switch4' */

    /* DataTypeConversion: '<S5>/Data Type Conversion9' */
    VCM20_B.DataTypeConversion9 = VCM20_B.Switch4_b;

    /* DataTypeConversion: '<S5>/Data Type Conversion11' */
    VCM20_B.DataTypeConversion11 = VCM20_B.DcOn;

    /* Switch: '<S123>/Switch4' incorporates:
     *  Constant: '<Root>/MC4_sw'
     */
    if (VCM20_P.MC4_sw_Value > VCM20_P.Switch4_Threshold_l) {
      /* Switch: '<S123>/Switch4' */
      VCM20_B.Switch4_f = VCM20_B.ReadyToDrive_and;
    } else {
      /* Switch: '<S123>/Switch4' incorporates:
       *  Constant: '<S123>/Constant'
       */
      VCM20_B.Switch4_f = VCM20_P.Constant_Value_bn;
    }

    /* End of Switch: '<S123>/Switch4' */

    /* DataTypeConversion: '<S5>/Data Type Conversion10' */
    VCM20_B.DataTypeConversion10_l = VCM20_B.Switch4_f;

    /* Gain: '<S122>/Gain9' */
    VCM20_B.FLspeed = VCM20_P.Gain9_Gain_e * VCM20_B.Saturation2;

    /* MinMax: '<S125>/Min3' */
    u0_0 = VCM20_B.FLspeed;
    u1_0 = VCM20_B.Saturation2;
    if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    /* MinMax: '<S125>/Min3' */
    VCM20_B.Min3 = u1_0;

    /* Switch: '<S160>/Switch' incorporates:
     *  Concatenate: '<S160>/Vector Concatenate'
     *  Constant: '<Root>/80kWOn'
     *  Constant: '<S160>/all power limit'
     *  Constant: '<S160>/each power limit '
     */
    if (VCM20_P.u0kWOn_Value > VCM20_P.Switch_Threshold_gk) {
      VCM20_B.VectorConcatenate_m[0] = VCM20_P.allpowerlimit_Value_h;
    } else {
      VCM20_B.VectorConcatenate_m[0] = VCM20_P.eachpowerlimit_Value_d;
    }

    /* End of Switch: '<S160>/Switch' */

    /* Saturate: '<S160>/Saturation' incorporates:
     *  Concatenate: '<S160>/Vector Concatenate'
     */
    u0_0 = VCM20_B.SFunction1_o9_g;
    u1 = VCM20_P.Saturation_LowerSat_fq;
    u2 = VCM20_P.Saturation_UpperSat_n3;
    if (u0_0 > u2) {
      u0_0 = u2;
    } else if (u0_0 < u1) {
      u0_0 = u1;
    }

    VCM20_B.VectorConcatenate_m[1] = u0_0;

    /* End of Saturate: '<S160>/Saturation' */

    /* Fcn: '<S160>/Fcn' */
    VCM20_B.Fcn_o = VCM20_B.VectorConcatenate_m[0] / (6.2831853071795862 *
      VCM20_B.VectorConcatenate_m[1] / 60.0);
  }

  /* Switch: '<S125>/MC4' incorporates:
   *  Constant: '<Root>/MC4_sw'
   */
  if (VCM20_P.MC4_sw_Value > VCM20_P.MC4_Threshold) {
    /* Switch: '<S125>/MC4' */
    VCM20_B.MC4 = VCM20_B.FLtrq;
  } else {
    /* Switch: '<S125>/MC4' incorporates:
     *  Constant: '<S125>/Constant'
     */
    VCM20_B.MC4 = VCM20_P.Constant_Value_cj;
  }

  /* End of Switch: '<S125>/MC4' */

  /* RelationalOperator: '<S166>/LowerRelop1' */
  VCM20_B.LowerRelop1_ig = (VCM20_B.MC4 > VCM20_B.Fcn_o);

  /* RelationalOperator: '<S166>/UpperRelop' incorporates:
   *  Constant: '<S160>/Constant1'
   */
  VCM20_B.UpperRelop_md = (VCM20_B.MC4 < VCM20_P.Constant1_Value_d);

  /* Switch: '<S166>/Switch' */
  if (VCM20_B.UpperRelop_md) {
    /* Switch: '<S166>/Switch' incorporates:
     *  Constant: '<S160>/Constant1'
     */
    VCM20_B.Switch_ge = VCM20_P.Constant1_Value_d;
  } else {
    /* Switch: '<S166>/Switch' */
    VCM20_B.Switch_ge = VCM20_B.MC4;
  }

  /* End of Switch: '<S166>/Switch' */

  /* Switch: '<S166>/Switch2' */
  if (VCM20_B.LowerRelop1_ig) {
    /* Switch: '<S166>/Switch2' */
    VCM20_B.Switch2_ff = VCM20_B.Fcn_o;
  } else {
    /* Switch: '<S166>/Switch2' */
    VCM20_B.Switch2_ff = VCM20_B.Switch_ge;
  }

  /* End of Switch: '<S166>/Switch2' */

  /* Gain: '<S120>/Gain2' */
  VCM20_B.Gain2_b1 = VCM20_P.Gain2_Gain_d * VCM20_B.Switch2_ff;

  /* Quantizer: '<S120>/Quantizer1' */
  u1_0 = VCM20_B.Gain2_b1;

  /* Quantizer: '<S120>/Quantizer1' */
  VCM20_B.Quantizer1_n = rt_roundd_snf(u1_0 / VCM20_P.Quantizer1_Interval_f) *
    VCM20_P.Quantizer1_Interval_f;

  /* Switch: '<S162>/Switch' */
  if (VCM20_B.Quantizer1_n > VCM20_P.Switch_Threshold_ik) {
    /* Switch: '<S162>/Switch' */
    VCM20_B.Switch_f4 = VCM20_B.Min3;
  } else {
    /* Switch: '<S162>/Switch' incorporates:
     *  Constant: '<S162>/Constant1'
     */
    VCM20_B.Switch_f4 = VCM20_P.Constant1_Value_g;
  }

  /* End of Switch: '<S162>/Switch' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Delay: '<S162>/Delay' */
    VCM20_B.Delay_gq = VCM20_DW.Delay_DSTATE_pw;
  }

  /* MinMax: '<S162>/MinMax' */
  u0_0 = VCM20_B.Switch_f4;
  u1_0 = VCM20_B.Delay_gq;
  if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S162>/MinMax' */
  VCM20_B.MinMax_l = u1_0;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Delay: '<S164>/Delay' */
    VCM20_B.Delay_n = VCM20_DW.Delay_DSTATE_hw;
  }

  /* MinMax: '<S164>/MinMax' */
  u0_0 = VCM20_B.Quantizer1_n;
  u1_0 = VCM20_B.Delay_n;
  if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S164>/MinMax' */
  VCM20_B.MinMax_lq = u1_0;

  /* MinMax: '<S161>/MinMax' incorporates:
   *  Constant: '<S161>/M set point MAX 0.1%Mn'
   */
  u0_0 = VCM20_B.MinMax_lq;
  u1_0 = VCM20_P.MsetpointMAX01Mn_Value_a;
  if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S161>/MinMax' */
  VCM20_B.Positivetorquelimit_p = u1_0;

  /* RelationalOperator: '<S165>/LowerRelop1' incorporates:
   *  Constant: '<S159>/power limit'
   */
  VCM20_B.LowerRelop1_pi = (VCM20_B.MC4 > VCM20_P.powerlimit_Value_j);
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S159>/Constant' incorporates:
     *  Concatenate: '<S159>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_g[0] = VCM20_P.Constant_Value_mr;

    /* Saturate: '<S159>/Saturation' incorporates:
     *  Concatenate: '<S159>/Vector Concatenate'
     */
    u0_0 = VCM20_B.SFunction1_o9_g;
    u1 = VCM20_P.Saturation_LowerSat_hg;
    u2 = VCM20_P.Saturation_UpperSat_hu;
    if (u0_0 > u2) {
      u0_0 = u2;
    } else if (u0_0 < u1) {
      u0_0 = u1;
    }

    VCM20_B.VectorConcatenate_g[1] = u0_0;

    /* End of Saturate: '<S159>/Saturation' */

    /* Fcn: '<S159>/Fcn' */
    VCM20_B.Fcn_oe = VCM20_B.VectorConcatenate_g[0] / (6.2831853071795862 *
      VCM20_B.VectorConcatenate_g[1] / 60.0);
  }

  /* RelationalOperator: '<S165>/UpperRelop' */
  VCM20_B.UpperRelop_lc = (VCM20_B.MC4 < VCM20_B.Fcn_oe);

  /* Switch: '<S165>/Switch' */
  if (VCM20_B.UpperRelop_lc) {
    /* Switch: '<S165>/Switch' */
    VCM20_B.Switch_oi0 = VCM20_B.Fcn_oe;
  } else {
    /* Switch: '<S165>/Switch' */
    VCM20_B.Switch_oi0 = VCM20_B.MC4;
  }

  /* End of Switch: '<S165>/Switch' */

  /* Switch: '<S165>/Switch2' */
  if (VCM20_B.LowerRelop1_pi) {
    /* Switch: '<S165>/Switch2' incorporates:
     *  Constant: '<S159>/power limit'
     */
    VCM20_B.Switch2_lu = VCM20_P.powerlimit_Value_j;
  } else {
    /* Switch: '<S165>/Switch2' */
    VCM20_B.Switch2_lu = VCM20_B.Switch_oi0;
  }

  /* End of Switch: '<S165>/Switch2' */

  /* Gain: '<S120>/Gain6' */
  VCM20_B.Gain6_f = VCM20_P.Gain6_Gain_l * VCM20_B.Switch2_lu;

  /* Quantizer: '<S120>/Quantizer2' */
  u1_0 = VCM20_B.Gain6_f;

  /* Quantizer: '<S120>/Quantizer2' */
  VCM20_B.Quantizer2_e = rt_roundd_snf(u1_0 / VCM20_P.Quantizer2_Interval_po) *
    VCM20_P.Quantizer2_Interval_po;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Delay: '<S163>/Delay' */
    VCM20_B.Delay_a = VCM20_DW.Delay_DSTATE_m;
  }

  /* MinMax: '<S163>/MinMax' */
  u0_0 = VCM20_B.Quantizer2_e;
  u1_0 = VCM20_B.Delay_a;
  if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
    u1_0 = u0_0;
  }

  /* MinMax: '<S163>/MinMax' */
  VCM20_B.MinMax_iq = u1_0;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S66>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "VCM_41" Id:393 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X189]);

      /* Convert timestamp */
      if (can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X189]->processed) {
        can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X189]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X189]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "AMK_bInverterOn" (8|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion9 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "AMK_bDcOn" (9|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion11 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 1;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "AMK_bEnable" (10|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion10_l ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 2;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "AMK_bErrorReset" (11|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.AMK_bErrorReset ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 3;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "Speed Setpoint in 1/rpm" (16|16, standard signal, signed int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.MinMax_l ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_l ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_l ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x0000FFFF;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte1;

        /* ...... "Positive toruque limit in 0,1 % Mn" (32|16, standard signal, signed int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.Positivetorquelimit_p ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Positivetorquelimit_p ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.Positivetorquelimit_p ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x0000FFFF;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte1;

        /* ...... "Negative torque limit in 0,1% Mn" (48|16, standard signal, signed int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        if (( VCM20_B.MinMax_iq ) < -0.5)
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_iq ) - 0.5);
        else
          CAN_Sgn.SignedSgn = (Int32) (( VCM20_B.MinMax_iq ) + 0.5);
        CAN_Sgn.SignedSgn &= 0x0000FFFF;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte1;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X189], 8,
                       &(CAN_Msg[0]), delayTime);
    }
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* S-Function (rti_commonblock): '<S334>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* write output state value to digital output channel 4 on port 1 */
    dio_tp4_digout_write(DIO_TP4_1_MODULE_ADDR, 1, DIO_TP4_MASK_CH4, (UInt16)
                         (VCM20_B.BrakeSW << 3));
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[3] == 0) {
    /* Logic: '<S11>/MinusSW' */
    VCM20_B.MinusSW = !(VCM20_B.SFunction1_o3_d != 0.0);

    /* Logic: '<S11>/PlusSW' */
    VCM20_B.PlusSW = !(VCM20_B.SFunction1_o2_k != 0.0);

    /* Logic: '<S11>/SelectSW' */
    VCM20_B.SelectSW = !(VCM20_B.SFunction1_o1_o != 0.0);

    /* Chart: '<S16>/Steer Chart' */
    /* Gateway: VCMOutput/Steer Chart */
    /* During: VCMOutput/Steer Chart */
    if (VCM20_DW.is_active_c3_VCM20 == 0U) {
      /* Entry: VCMOutput/Steer Chart */
      VCM20_DW.is_active_c3_VCM20 = 1U;

      /* Entry Internal: VCMOutput/Steer Chart */
      /* Transition: '<S335>:2' */
      VCM20_DW.is_c3_VCM20 = VCM20_IN_Start;

      /* Entry 'Start': '<S335>:1' */
      VCM20_B.gain1 = 5.0;
      VCM20_DW.gain1max = 20.0;
      VCM20_DW.gain1min = 0.0;
      VCM20_B.gain2 = 10.0;
      VCM20_DW.gain2max = 15.0;
      VCM20_DW.gain2min = 5.0;
    } else {
      switch (VCM20_DW.is_c3_VCM20) {
       case VCM20_IN_Display1:
        /* During 'Display1': '<S335>:39' */
        if (VCM20_B.SelectSW) {
          /* Transition: '<S335>:40' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Select1;
        }
        break;

       case VCM20_IN_Gain1:
        /* During 'Gain1': '<S335>:3' */
        if (VCM20_B.PlusSW && (VCM20_B.gain1 < VCM20_DW.gain1max)) {
          /* Transition: '<S335>:6' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Plus_Gain1;

          /* Entry 'Plus_Gain1': '<S335>:5' */
          VCM20_B.gain1++;
          VCM20_B.LCDnumber = VCM20_B.gain1;
        } else if (VCM20_B.MinusSW && (VCM20_B.gain1 > VCM20_DW.gain1min)) {
          /* Transition: '<S335>:8' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Minus_Gain1;

          /* Entry 'Minus_Gain1': '<S335>:9' */
          VCM20_B.gain1--;
          VCM20_B.LCDnumber = VCM20_B.gain1;
        } else if (VCM20_B.SelectSW) {
          /* Transition: '<S335>:12' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Select2;
        }
        break;

       case VCM20_IN_Gain2:
        /* During 'Gain2': '<S335>:14' */
        if (VCM20_B.PlusSW && (VCM20_B.gain2 < VCM20_DW.gain2max)) {
          /* Transition: '<S335>:18' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Plus_Gain2;

          /* Entry 'Plus_Gain2': '<S335>:19' */
          VCM20_B.gain2++;
          VCM20_B.LCDnumber = VCM20_B.gain2;
        } else if (VCM20_B.MinusSW && (VCM20_B.gain2 > VCM20_DW.gain2min)) {
          /* Transition: '<S335>:16' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Minus_Gain2;

          /* Entry 'Minus_Gain2': '<S335>:17' */
          VCM20_B.gain2--;
          VCM20_B.LCDnumber = VCM20_B.gain2;
        } else if (VCM20_B.SelectSW) {
          /* Transition: '<S335>:22' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Select3;
        }
        break;

       case VCM20_IN_Minus_Gain1:
        /* During 'Minus_Gain1': '<S335>:9' */
        if (!VCM20_B.MinusSW) {
          /* Transition: '<S335>:10' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Gain1;

          /* Entry 'Gain1': '<S335>:3' */
          for (i = 0; i < 8; i++) {
            VCM20_B.LCDtext[i] = c[i];
          }

          VCM20_B.LCDnumber = VCM20_B.gain1;
        }
        break;

       case VCM20_IN_Minus_Gain2:
        /* During 'Minus_Gain2': '<S335>:17' */
        if (!VCM20_B.MinusSW) {
          /* Transition: '<S335>:15' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Gain2;

          /* Entry 'Gain2': '<S335>:14' */
          for (i = 0; i < 8; i++) {
            VCM20_B.LCDtext[i] = d[i];
          }

          VCM20_B.LCDnumber = VCM20_B.gain2;
        }
        break;

       case VCM20_IN_Plus_Gain1:
        /* During 'Plus_Gain1': '<S335>:5' */
        if (!VCM20_B.PlusSW) {
          /* Transition: '<S335>:7' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Gain1;

          /* Entry 'Gain1': '<S335>:3' */
          for (i = 0; i < 8; i++) {
            VCM20_B.LCDtext[i] = c[i];
          }

          VCM20_B.LCDnumber = VCM20_B.gain1;
        }
        break;

       case VCM20_IN_Plus_Gain2:
        /* During 'Plus_Gain2': '<S335>:19' */
        if (!VCM20_B.PlusSW) {
          /* Transition: '<S335>:21' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Gain2;

          /* Entry 'Gain2': '<S335>:14' */
          for (i = 0; i < 8; i++) {
            VCM20_B.LCDtext[i] = d[i];
          }

          VCM20_B.LCDnumber = VCM20_B.gain2;
        }
        break;

       case VCM20_IN_Select1:
        /* During 'Select1': '<S335>:38' */
        if (!VCM20_B.SelectSW) {
          /* Transition: '<S335>:41' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Gain1;

          /* Entry 'Gain1': '<S335>:3' */
          for (i = 0; i < 8; i++) {
            VCM20_B.LCDtext[i] = c[i];
          }

          VCM20_B.LCDnumber = VCM20_B.gain1;
        }
        break;

       case VCM20_IN_Select2:
        /* During 'Select2': '<S335>:31' */
        if (!VCM20_B.SelectSW) {
          /* Transition: '<S335>:20' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Gain2;

          /* Entry 'Gain2': '<S335>:14' */
          for (i = 0; i < 8; i++) {
            VCM20_B.LCDtext[i] = d[i];
          }

          VCM20_B.LCDnumber = VCM20_B.gain2;
        }
        break;

       case VCM20_IN_Select3:
        /* During 'Select3': '<S335>:32' */
        if (!VCM20_B.SelectSW) {
          /* Transition: '<S335>:23' */
          VCM20_DW.is_c3_VCM20 = VCM20_IN_Display1;

          /* Entry 'Display1': '<S335>:39' */
          for (i = 0; i < 8; i++) {
            VCM20_B.LCDtext[i] = b[i];
          }

          VCM20_B.LCDnumber = VCM20_B.SOLUTION_MODE;
        }
        break;

       default:
        /* During 'Start': '<S335>:1' */
        /* Transition: '<S335>:4' */
        VCM20_DW.is_c3_VCM20 = VCM20_IN_Display1;

        /* Entry 'Display1': '<S335>:39' */
        for (i = 0; i < 8; i++) {
          VCM20_B.LCDtext[i] = b[i];
        }

        VCM20_B.LCDnumber = VCM20_B.SOLUTION_MODE;
        break;
      }
    }

    /* End of Chart: '<S16>/Steer Chart' */
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* DataTypeConversion: '<S16>/Data Type Conversion14' */
    VCM20_B.DataTypeConversion14 = VCM20_B.LaunchLED;

    /* DataTypeConversion: '<S16>/Data Type Conversion13' */
    VCM20_B.DataTypeConversion13 = VCM20_B.All_enable;

    /* DataTypeConversion: '<S16>/Data Type Conversion12' */
    VCM20_B.DataTypeConversion12 = VCM20_B.All_QuitDcOn;
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* DataTypeConversion: '<S16>/Data Type Conversion11' */
    VCM20_B.DataTypeConversion11_f = VCM20_B.BrakeSW;
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[3] == 0) {
    /* Logic: '<S16>/NOT' */
    VCM20_B.NOT_b = !VCM20_B.SFunction1_o;

    /* DataTypeConversion: '<S16>/Data Type Conversion8' */
    VCM20_B.DataTypeConversion8_f = VCM20_B.NOT_b;

    /* Logic: '<S16>/NOT1' */
    VCM20_B.NOT1 = !VCM20_B.SFunction1_b;

    /* DataTypeConversion: '<S16>/Data Type Conversion9' */
    VCM20_B.DataTypeConversion9_o = VCM20_B.NOT1;

    /* Logic: '<S16>/NOT2' */
    VCM20_B.NOT2 = !VCM20_B.SFunction1_j;

    /* DataTypeConversion: '<S16>/Data Type Conversion10' */
    VCM20_B.DataTypeConversion10_a = VCM20_B.NOT2;

    /* S-Function (rti_commonblock): '<S336>/S-Function1' incorporates:
     *  Constant: '<S16>/LED1'
     *  Constant: '<S16>/LED2'
     *  Constant: '<S16>/LED3'
     *  Constant: '<S16>/LED5'
     *  Constant: '<S16>/LED_Displaymode'
     *  Constant: '<S16>/LED_brightness'
     *  Constant: '<S16>/LED_number2'
     *  Constant: '<S16>/buzzer'
     */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "Steer LED" Id:528 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X210]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X210]->processed) {
        can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X210]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X210]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "LED1_LCD" (0|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_P.LED1_Value ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LED2_blue" (1|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_P.LED2_Value ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 1;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LED3_green" (2|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_P.LED3_Value ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 2;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LED4_yellow" (3|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion14 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 3;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LED5_blue" (4|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_P.LED5_Value ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 4;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LED6_green" (5|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion13 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 5;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LED7_yellow" (6|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion12 ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 6;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LED8_red" (7|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion11_f ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 7;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LED9_red" (8|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion8_f ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LED10_red" (9|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion9_o ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 1;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LED11_red" (10|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion10_a ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 2;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "buzzer" (11|1, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_P.buzzer_Value ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x00000001;
        CAN_Sgn.UnsignedSgn = ((UInt32)CAN_Sgn.UnsignedSgn) << 3;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LED_brightness" (16|8, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_P.LED_brightness_Value ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LCD_number1" (24|16, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.LCDnumber ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte1;

        /* ...... "LCD_Displaymode" (40|8, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_P.LED_Displaymode_Value ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LCD_number2" (48|16, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_P.LED_number2_Value ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte1;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X210], 8,
                       &(CAN_Msg[0]), delayTime);
    }
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[3] == 0) {
    /* DataTypeConversion: '<S16>/Data Type Conversion' */
    VCM20_B.DataTypeConversion_l = VCM20_B.LCDtext[0];

    /* DataTypeConversion: '<S16>/Data Type Conversion1' */
    VCM20_B.DataTypeConversion1_o = VCM20_B.LCDtext[1];

    /* DataTypeConversion: '<S16>/Data Type Conversion2' */
    VCM20_B.DataTypeConversion2_h = VCM20_B.LCDtext[2];

    /* DataTypeConversion: '<S16>/Data Type Conversion3' */
    VCM20_B.DataTypeConversion3_i = VCM20_B.LCDtext[3];

    /* DataTypeConversion: '<S16>/Data Type Conversion4' */
    VCM20_B.DataTypeConversion4_b = VCM20_B.LCDtext[4];

    /* DataTypeConversion: '<S16>/Data Type Conversion5' */
    VCM20_B.DataTypeConversion5_h = VCM20_B.LCDtext[5];

    /* DataTypeConversion: '<S16>/Data Type Conversion6' */
    VCM20_B.DataTypeConversion6_d = VCM20_B.LCDtext[6];

    /* DataTypeConversion: '<S16>/Data Type Conversion7' */
    VCM20_B.DataTypeConversion7_j = VCM20_B.LCDtext[7];

    /* S-Function (rti_commonblock): '<S337>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "Steer LCD" Id:529 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X211]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X211]->processed) {
        can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X211]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X211]->timestamp);
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "LCD_word1" (0|8, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion_l ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LCD_word2" (8|8, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion1_o ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LCD_word3" (16|8, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion2_h ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LCD_word4" (24|8, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion3_i ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LCD_word5" (32|8, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion4_b ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[4] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LCD_word6" (40|8, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion5_h ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[5] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LCD_word7" (48|8, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion6_d ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[6] |= CAN_Sgn.SgnBytes.Byte0;

        /* ...... "LCD_word8" (56|8, standard signal, unsigned int, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.UnsignedSgn = (UInt32) (( VCM20_B.DataTypeConversion7_j ) + 0.5);
        CAN_Sgn.UnsignedSgn &= 0x000000FF;
        CAN_Msg[7] |= CAN_Sgn.SgnBytes.Byte0;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X211], 8,
                       &(CAN_Msg[0]), delayTime);
    }
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Sum: '<S72>/Add1' incorporates:
     *  Constant: '<S72>/SWAS1_signal_0deg'
     *  Constant: '<S72>/SWAS1_signal_plus'
     */
    VCM20_B.Add1_jn = VCM20_P.SWAS1_signal_plus_Value -
      VCM20_P.SWAS1_signal_0deg_Value;
  }

  /* TransferFcn: '<S72>/Transfer Fcn' */
  VCM20_B.TransferFcn_a = 0.0;
  VCM20_B.TransferFcn_a += VCM20_P.TransferFcn_C_a *
    VCM20_X.TransferFcn_CSTATE_e;

  /* Sum: '<S72>/Add' incorporates:
   *  Constant: '<S72>/SWAS1_signal_0deg'
   */
  VCM20_B.Add_b = VCM20_B.TransferFcn_a - VCM20_P.SWAS1_signal_0deg_Value;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Sum: '<S72>/Add2' incorporates:
     *  Constant: '<S72>/SWAS1_signal_0deg'
     *  Constant: '<S72>/SWAS1_signal_minus'
     */
    VCM20_B.Add2_o = VCM20_P.SWAS1_signal_minus_Value -
      VCM20_P.SWAS1_signal_0deg_Value;
  }

  /* Switch: '<S72>/Switch' */
  if (VCM20_B.Add_b > VCM20_P.Switch_Threshold_ob) {
    /* Product: '<S72>/Divide' */
    VCM20_B.Divide_jq = 1.0 / VCM20_B.Add1_jn * VCM20_B.Add_b;

    /* Gain: '<S72>/LeftMax_deg' */
    VCM20_B.LeftMax_deg_a = VCM20_P.LeftMax_deg_Gain * VCM20_B.Divide_jq;

    /* Switch: '<S72>/Switch' */
    VCM20_B.Switch_co = VCM20_B.LeftMax_deg_a;
  } else {
    /* Product: '<S72>/Divide1' */
    VCM20_B.Divide1_hi = VCM20_B.Add_b / VCM20_B.Add2_o;

    /* Gain: '<S72>/Gain1' */
    VCM20_B.Gain1_f = VCM20_P.Gain1_Gain * VCM20_B.Divide1_hi;

    /* Gain: '<S72>/RightMax_deg' */
    VCM20_B.RightMax_deg_d = VCM20_P.RightMax_deg_Gain * VCM20_B.Gain1_f;

    /* Switch: '<S72>/Switch' */
    VCM20_B.Switch_co = VCM20_B.RightMax_deg_d;
  }

  /* End of Switch: '<S72>/Switch' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* Gain: '<S11>/Gain' */
    VCM20_B.APPS1_noncalibrated = VCM20_P.Gain_Gain_ap * VCM20_B.SFunction1_o1_n;

    /* RelationalOperator: '<S84>/Compare' incorporates:
     *  Constant: '<S84>/Constant'
     */
    VCM20_B.Compare_np = (VCM20_B.APPS1_noncalibrated >=
                          VCM20_P.CompareToConstant1_const_f);

    /* RelationalOperator: '<S86>/Compare' incorporates:
     *  Constant: '<S86>/Constant'
     */
    VCM20_B.Compare_ky = (VCM20_B.APPS1_noncalibrated <=
                          VCM20_P.CompareToConstant3_const_f);

    /* Logic: '<S67>/Logical Operator' */
    VCM20_B.LogicalOperator_b = (VCM20_B.Compare_np && VCM20_B.Compare_ky);

    /* Delay: '<S79>/Delay' */
    VCM20_B.Delay_k = VCM20_B.LogicalOperator_b;

    /* Delay: '<S79>/Delay1' */
    VCM20_B.Delay1_l = VCM20_DW.Delay1_DSTATE_p;

    /* Delay: '<S79>/Delay2' */
    VCM20_B.Delay2_l = VCM20_DW.Delay2_DSTATE_e[0];

    /* Delay: '<S79>/Delay3' */
    VCM20_B.Delay3_o = VCM20_DW.Delay3_DSTATE_n[0];

    /* Delay: '<S79>/Delay4' */
    VCM20_B.Delay4_g = VCM20_DW.Delay4_DSTATE_j[0];

    /* Delay: '<S79>/Delay5' */
    VCM20_B.Delay5_e = VCM20_DW.Delay5_DSTATE_d[0];

    /* Delay: '<S79>/Delay6' */
    VCM20_B.Delay6_b = VCM20_DW.Delay6_DSTATE_i[0];

    /* Delay: '<S79>/Delay7' */
    VCM20_B.Delay7_cb = VCM20_DW.Delay7_DSTATE_g[0];

    /* Delay: '<S79>/Delay8' */
    VCM20_B.Delay8_c = VCM20_DW.Delay8_DSTATE_c[0];

    /* Delay: '<S79>/Delay9' */
    VCM20_B.Delay9_c = VCM20_DW.Delay9_DSTATE_h[0];

    /* Logic: '<S79>/Logical Operator' */
    VCM20_B.LogicalOperator_n = (VCM20_B.Delay_k && VCM20_B.Delay1_l &&
      VCM20_B.Delay2_l && VCM20_B.Delay3_o && VCM20_B.Delay4_g &&
      VCM20_B.Delay5_e && VCM20_B.Delay6_b && VCM20_B.Delay7_cb &&
      VCM20_B.Delay8_c && VCM20_B.Delay9_c);

    /* Gain: '<S11>/Gain1' */
    VCM20_B.APPS2_noncalibrated = VCM20_P.Gain1_Gain_l * VCM20_B.SFunction1_o2_j;

    /* RelationalOperator: '<S85>/Compare' incorporates:
     *  Constant: '<S85>/Constant'
     */
    VCM20_B.Compare_fi = (VCM20_B.APPS2_noncalibrated >=
                          VCM20_P.CompareToConstant2_const_p);

    /* RelationalOperator: '<S87>/Compare' incorporates:
     *  Constant: '<S87>/Constant'
     */
    VCM20_B.Compare_g = (VCM20_B.APPS2_noncalibrated <=
                         VCM20_P.CompareToConstant4_const);

    /* Logic: '<S67>/Logical Operator2' */
    VCM20_B.LogicalOperator2_n = (VCM20_B.Compare_fi && VCM20_B.Compare_g);

    /* Delay: '<S80>/Delay' */
    VCM20_B.Delay_d = VCM20_B.LogicalOperator2_n;

    /* Delay: '<S80>/Delay1' */
    VCM20_B.Delay1_lz = VCM20_DW.Delay1_DSTATE_f;

    /* Delay: '<S80>/Delay2' */
    VCM20_B.Delay2_a = VCM20_DW.Delay2_DSTATE_k[0];

    /* Delay: '<S80>/Delay3' */
    VCM20_B.Delay3_a = VCM20_DW.Delay3_DSTATE_gb[0];

    /* Delay: '<S80>/Delay4' */
    VCM20_B.Delay4_gj = VCM20_DW.Delay4_DSTATE_a[0];

    /* Delay: '<S80>/Delay5' */
    VCM20_B.Delay5_h = VCM20_DW.Delay5_DSTATE_k[0];

    /* Delay: '<S80>/Delay6' */
    VCM20_B.Delay6_jm = VCM20_DW.Delay6_DSTATE_n[0];

    /* Delay: '<S80>/Delay7' */
    VCM20_B.Delay7_f = VCM20_DW.Delay7_DSTATE_m[0];

    /* Delay: '<S80>/Delay8' */
    VCM20_B.Delay8_l = VCM20_DW.Delay8_DSTATE_d[0];

    /* Delay: '<S80>/Delay9' */
    VCM20_B.Delay9_d = VCM20_DW.Delay9_DSTATE_p[0];

    /* Logic: '<S80>/Logical Operator' */
    VCM20_B.LogicalOperator_j = (VCM20_B.Delay_d && VCM20_B.Delay1_lz &&
      VCM20_B.Delay2_a && VCM20_B.Delay3_a && VCM20_B.Delay4_gj &&
      VCM20_B.Delay5_h && VCM20_B.Delay6_jm && VCM20_B.Delay7_f &&
      VCM20_B.Delay8_l && VCM20_B.Delay9_d);
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* Delay: '<S89>/Delay10' */
    VCM20_B.Delay10 = VCM20_B.SFunction1_o1_n;

    /* Delay: '<S89>/Delay11' */
    VCM20_B.Delay11 = VCM20_DW.Delay11_DSTATE;

    /* Delay: '<S89>/Delay12' */
    VCM20_B.Delay12 = VCM20_DW.Delay12_DSTATE[0];

    /* Delay: '<S89>/Delay13' */
    VCM20_B.Delay13 = VCM20_DW.Delay13_DSTATE[0];

    /* Delay: '<S89>/Delay14' */
    VCM20_B.Delay14 = VCM20_DW.Delay14_DSTATE[0];

    /* Delay: '<S89>/Delay15' */
    VCM20_B.Delay15 = VCM20_DW.Delay15_DSTATE[0];

    /* Delay: '<S89>/Delay16' */
    VCM20_B.Delay16 = VCM20_DW.Delay16_DSTATE[0];

    /* Delay: '<S89>/Delay17' */
    VCM20_B.Delay17 = VCM20_DW.Delay17_DSTATE[0];

    /* Delay: '<S89>/Delay18' */
    VCM20_B.Delay18 = VCM20_DW.Delay18_DSTATE[0];

    /* Delay: '<S89>/Delay19' */
    VCM20_B.Delay19 = VCM20_DW.Delay19_DSTATE[0];

    /* MinMax: '<S89>/Min1' */
    u0_0 = VCM20_B.Delay10;
    u1_0 = VCM20_B.Delay11;
    if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    u1 = VCM20_B.Delay12;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay13;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay14;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay15;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay16;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay17;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay18;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay19;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    /* MinMax: '<S89>/Min1' */
    VCM20_B.Min1_g = u1_0;

    /* Delay: '<S89>/Delay' */
    VCM20_B.Delay_ih = VCM20_B.SFunction1_o1_n;

    /* Delay: '<S89>/Delay1' */
    VCM20_B.Delay1 = VCM20_DW.Delay1_DSTATE;

    /* Delay: '<S89>/Delay2' */
    VCM20_B.Delay2 = VCM20_DW.Delay2_DSTATE[0];

    /* Delay: '<S89>/Delay3' */
    VCM20_B.Delay3 = VCM20_DW.Delay3_DSTATE[0];

    /* Delay: '<S89>/Delay4' */
    VCM20_B.Delay4 = VCM20_DW.Delay4_DSTATE[0];

    /* Delay: '<S89>/Delay5' */
    VCM20_B.Delay5 = VCM20_DW.Delay5_DSTATE[0];

    /* Delay: '<S89>/Delay6' */
    VCM20_B.Delay6 = VCM20_DW.Delay6_DSTATE[0];

    /* Delay: '<S89>/Delay7' */
    VCM20_B.Delay7 = VCM20_DW.Delay7_DSTATE[0];

    /* Delay: '<S89>/Delay8' */
    VCM20_B.Delay8 = VCM20_DW.Delay8_DSTATE[0];

    /* Delay: '<S89>/Delay9' */
    VCM20_B.Delay9 = VCM20_DW.Delay9_DSTATE[0];

    /* MinMax: '<S89>/Min' */
    u0_0 = VCM20_B.Delay_ih;
    u1_0 = VCM20_B.Delay1;
    if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    u1 = VCM20_B.Delay2;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay3;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay4;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay5;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay6;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay7;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay8;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay9;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    /* MinMax: '<S89>/Min' */
    VCM20_B.Min_k = u1_0;
    VCM20_MovingAverage(VCM20_B.SFunction1_o1_n, &VCM20_B.MovingAverage_p,
                        &VCM20_DW.MovingAverage_p);

    /* Delay: '<S90>/Delay10' */
    VCM20_B.Delay10_h = VCM20_B.SFunction1_o2_j;

    /* Delay: '<S90>/Delay11' */
    VCM20_B.Delay11_f = VCM20_DW.Delay11_DSTATE_o;

    /* Delay: '<S90>/Delay12' */
    VCM20_B.Delay12_l = VCM20_DW.Delay12_DSTATE_i[0];

    /* Delay: '<S90>/Delay13' */
    VCM20_B.Delay13_h = VCM20_DW.Delay13_DSTATE_p[0];

    /* Delay: '<S90>/Delay14' */
    VCM20_B.Delay14_n = VCM20_DW.Delay14_DSTATE_k[0];

    /* Delay: '<S90>/Delay15' */
    VCM20_B.Delay15_p = VCM20_DW.Delay15_DSTATE_n[0];

    /* Delay: '<S90>/Delay16' */
    VCM20_B.Delay16_d = VCM20_DW.Delay16_DSTATE_c[0];

    /* Delay: '<S90>/Delay17' */
    VCM20_B.Delay17_i = VCM20_DW.Delay17_DSTATE_e[0];

    /* Delay: '<S90>/Delay18' */
    VCM20_B.Delay18_n = VCM20_DW.Delay18_DSTATE_a[0];

    /* Delay: '<S90>/Delay19' */
    VCM20_B.Delay19_j = VCM20_DW.Delay19_DSTATE_h[0];

    /* MinMax: '<S90>/Min1' */
    u0_0 = VCM20_B.Delay10_h;
    u1_0 = VCM20_B.Delay11_f;
    if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    u1 = VCM20_B.Delay12_l;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay13_h;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay14_n;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay15_p;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay16_d;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay17_i;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay18_n;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay19_j;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    /* MinMax: '<S90>/Min1' */
    VCM20_B.Min1_e = u1_0;

    /* Delay: '<S90>/Delay' */
    VCM20_B.Delay_nn = VCM20_B.SFunction1_o2_j;

    /* Delay: '<S90>/Delay1' */
    VCM20_B.Delay1_i = VCM20_DW.Delay1_DSTATE_i;

    /* Delay: '<S90>/Delay2' */
    VCM20_B.Delay2_g = VCM20_DW.Delay2_DSTATE_p[0];

    /* Delay: '<S90>/Delay3' */
    VCM20_B.Delay3_c = VCM20_DW.Delay3_DSTATE_o[0];

    /* Delay: '<S90>/Delay4' */
    VCM20_B.Delay4_o = VCM20_DW.Delay4_DSTATE_i[0];

    /* Delay: '<S90>/Delay5' */
    VCM20_B.Delay5_n = VCM20_DW.Delay5_DSTATE_m[0];

    /* Delay: '<S90>/Delay6' */
    VCM20_B.Delay6_a = VCM20_DW.Delay6_DSTATE_p[0];

    /* Delay: '<S90>/Delay7' */
    VCM20_B.Delay7_n = VCM20_DW.Delay7_DSTATE_k[0];

    /* Delay: '<S90>/Delay8' */
    VCM20_B.Delay8_j = VCM20_DW.Delay8_DSTATE_j[0];

    /* Delay: '<S90>/Delay9' */
    VCM20_B.Delay9_g = VCM20_DW.Delay9_DSTATE_a[0];

    /* MinMax: '<S90>/Min' */
    u0_0 = VCM20_B.Delay_nn;
    u1_0 = VCM20_B.Delay1_i;
    if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    u1 = VCM20_B.Delay2_g;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay3_c;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay4_o;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay5_n;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay6_a;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay7_n;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay8_j;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay9_g;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    /* MinMax: '<S90>/Min' */
    VCM20_B.Min_j = u1_0;
    VCM20_MovingAverage(VCM20_B.SFunction1_o2_j, &VCM20_B.MovingAverage_pn,
                        &VCM20_DW.MovingAverage_pn);
  }

  VCM20_MovingAverage(VCM20_B.TransferFcn_a, &VCM20_B.MovingAverage_pna,
                      &VCM20_DW.MovingAverage_pna);
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  VCM20_MovingAverage(VCM20_B.TransferFcn_f, &VCM20_B.MovingAverage_pnae,
                      &VCM20_DW.MovingAverage_pnae);
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* Delay: '<S91>/Delay10' */
    VCM20_B.Delay10_m = VCM20_B.SFunction1_o2_a;

    /* Delay: '<S91>/Delay11' */
    VCM20_B.Delay11_m = VCM20_DW.Delay11_DSTATE_k;

    /* Delay: '<S91>/Delay12' */
    VCM20_B.Delay12_e = VCM20_DW.Delay12_DSTATE_h[0];

    /* Delay: '<S91>/Delay13' */
    VCM20_B.Delay13_f = VCM20_DW.Delay13_DSTATE_j[0];

    /* Delay: '<S91>/Delay14' */
    VCM20_B.Delay14_m = VCM20_DW.Delay14_DSTATE_l[0];

    /* Delay: '<S91>/Delay15' */
    VCM20_B.Delay15_c = VCM20_DW.Delay15_DSTATE_i[0];

    /* Delay: '<S91>/Delay16' */
    VCM20_B.Delay16_f = VCM20_DW.Delay16_DSTATE_b[0];

    /* Delay: '<S91>/Delay17' */
    VCM20_B.Delay17_c = VCM20_DW.Delay17_DSTATE_n[0];

    /* Delay: '<S91>/Delay18' */
    VCM20_B.Delay18_k = VCM20_DW.Delay18_DSTATE_h[0];

    /* Delay: '<S91>/Delay19' */
    VCM20_B.Delay19_d = VCM20_DW.Delay19_DSTATE_n[0];

    /* MinMax: '<S91>/Min1' */
    u0_0 = VCM20_B.Delay10_m;
    u1_0 = VCM20_B.Delay11_m;
    if ((u0_0 >= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    u1 = VCM20_B.Delay12_e;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay13_f;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay14_m;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay15_c;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay16_f;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay17_c;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay18_k;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay19_d;
    if ((!(u1_0 >= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    /* MinMax: '<S91>/Min1' */
    VCM20_B.Min1_i2 = u1_0;

    /* Delay: '<S91>/Delay' */
    VCM20_B.Delay_oz = VCM20_B.SFunction1_o2_a;

    /* Delay: '<S91>/Delay1' */
    VCM20_B.Delay1_io = VCM20_DW.Delay1_DSTATE_h;

    /* Delay: '<S91>/Delay2' */
    VCM20_B.Delay2_c = VCM20_DW.Delay2_DSTATE_l[0];

    /* Delay: '<S91>/Delay3' */
    VCM20_B.Delay3_m = VCM20_DW.Delay3_DSTATE_o3[0];

    /* Delay: '<S91>/Delay4' */
    VCM20_B.Delay4_p = VCM20_DW.Delay4_DSTATE_il[0];

    /* Delay: '<S91>/Delay5' */
    VCM20_B.Delay5_j = VCM20_DW.Delay5_DSTATE_mb[0];

    /* Delay: '<S91>/Delay6' */
    VCM20_B.Delay6_a2 = VCM20_DW.Delay6_DSTATE_l[0];

    /* Delay: '<S91>/Delay7' */
    VCM20_B.Delay7_k = VCM20_DW.Delay7_DSTATE_n[0];

    /* Delay: '<S91>/Delay8' */
    VCM20_B.Delay8_b = VCM20_DW.Delay8_DSTATE_jc[0];

    /* Delay: '<S91>/Delay9' */
    VCM20_B.Delay9_a = VCM20_DW.Delay9_DSTATE_e[0];

    /* MinMax: '<S91>/Min' */
    u0_0 = VCM20_B.Delay_oz;
    u1_0 = VCM20_B.Delay1_io;
    if ((u0_0 <= u1_0) || rtIsNaN(u1_0)) {
      u1_0 = u0_0;
    }

    u1 = VCM20_B.Delay2_c;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay3_m;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay4_p;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay5_j;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay6_a2;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay7_k;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay8_b;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    u1 = VCM20_B.Delay9_a;
    if ((!(u1_0 <= u1)) && (!rtIsNaN(u1))) {
      u1_0 = u1;
    }

    /* MinMax: '<S91>/Min' */
    VCM20_B.Min_c = u1_0;
    VCM20_MovingAverage(VCM20_B.SFunction1_o2_a, &VCM20_B.MovingAverage_pnaev,
                        &VCM20_DW.MovingAverage_pnaev);

    /* Gain: '<S11>/Gain2' */
    VCM20_B.Gain2_m = VCM20_P.Gain2_Gain_e2 * VCM20_B.Saturation_d;

    /* Logic: '<S77>/Logical Operator1' */
    VCM20_B.LogicalOperator1_b = !VCM20_B.AND_i;
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S101>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:514 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X202]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X202]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X202]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X202]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X202]->timestamp > 0.0) {
        if (can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X202]->processed) {
          CAN_Msg = can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X202]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "LapTime" (0|16, standard signal, signed int, big endian) */
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[1];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte0 = CAN_Sgn.SgnBytes.Byte3;
            CAN_Sgn.SgnBytes.Byte1 = CAN_Sgn.SgnBytes.Byte2;
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o1_cw = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "BestTime" (16|16, standard signal, signed int, big endian) */
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[3];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte0 = CAN_Sgn.SgnBytes.Byte3;
            CAN_Sgn.SgnBytes.Byte1 = CAN_Sgn.SgnBytes.Byte2;
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o2_jm = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "+-BestTime" (32|16, standard signal, signed int, big endian) */
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[5];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte0 = CAN_Sgn.SgnBytes.Byte3;
            CAN_Sgn.SgnBytes.Byte1 = CAN_Sgn.SgnBytes.Byte2;
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o3_dc = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "LapNumber" (48|8, standard signal, signed int, big endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[6];
            CAN_Sgn.SignedSgn &= 0x000000FF;
            if (CAN_Sgn.SignedSgn >> 7) {
              CAN_Sgn.SignedSgn |= 0xFFFFFF00;
            }

            VCM20_B.SFunction1_o4_b = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* S-Function (rti_commonblock): '<S102>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:513 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X201]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X201]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X201]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X201]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X201]->timestamp > 0.0) {
        if (can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X201]->processed) {
          CAN_Msg = can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X201]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "Current_mV" (0|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o1_k = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "yaw rate" (16|16, standard signal, signed int, big endian) */
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[3];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte0 = CAN_Sgn.SgnBytes.Byte3;
            CAN_Sgn.SgnBytes.Byte1 = CAN_Sgn.SgnBytes.Byte2;
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o2_ka = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "Latelal acc" (32|16, standard signal, signed int, big endian) */
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[5];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte0 = CAN_Sgn.SgnBytes.Byte3;
            CAN_Sgn.SgnBytes.Byte1 = CAN_Sgn.SgnBytes.Byte2;
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o3_fn = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "Inline acc" (48|16, standard signal, signed int, big endian) */
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[7];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte0 = CAN_Sgn.SgnBytes.Byte3;
            CAN_Sgn.SgnBytes.Byte1 = CAN_Sgn.SgnBytes.Byte2;
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o4_bh = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* S-Function (rti_commonblock): '<S109>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:337 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X151]);

      /* Convert timestamp */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X151]->processed) {
        can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X151]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X151]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X151]->timestamp > 0.0) {
        if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X151]->processed) {
          CAN_Msg = can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X151]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "MAG_X" (0|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o1_hy = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "MAG_Y" (16|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o2_g = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "MAG_Z" (32|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o3_ht = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* S-Function (rti_commonblock): '<S110>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:369 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X171]);

      /* Convert timestamp */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X171]->processed) {
        can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X171]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X171]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X171]->timestamp > 0.0) {
        if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X171]->processed) {
          CAN_Msg = can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X171]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "VEL_N" (0|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o1_bu = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "VEL_E" (16|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o2_id = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "VEL_D" (32|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            VCM20_B.SFunction1_o3_bh = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }
    }

    /* Gain: '<S14>/Gain11' */
    VCM20_B.MAG_X = VCM20_P.Gain11_Gain_p * VCM20_B.SFunction1_o1_hy;

    /* Gain: '<S14>/Gain12' */
    VCM20_B.MAG_Z = VCM20_P.Gain12_Gain_h * VCM20_B.SFunction1_o3_ht;

    /* Gain: '<S14>/Gain13' */
    VCM20_B.MAG_Y = VCM20_P.Gain13_Gain_m * VCM20_B.SFunction1_o2_g;

    /* Gain: '<S14>/Gain14' */
    VCM20_B.VEL_E = VCM20_P.Gain14_Gain_j * VCM20_B.SFunction1_o2_id;

    /* Gain: '<S14>/Gain15' */
    VCM20_B.VEL_N = VCM20_P.Gain15_Gain_c * VCM20_B.SFunction1_o1_bu;

    /* Gain: '<S14>/Gain16' */
    VCM20_B.VEL_D = VCM20_P.Gain16_Gain_j * VCM20_B.SFunction1_o3_bh;

    /* Gain: '<S14>/Gain17' */
    VCM20_B.yawrate = VCM20_P.Gain17_Gain * VCM20_B.SFunction1_o2_ka;

    /* Gain: '<S14>/Gain6' */
    VCM20_B.Inlineacc = VCM20_P.Gain6_Gain_b1 * VCM20_B.SFunction1_o4_bh;

    /* Gain: '<S14>/Gain7' */
    VCM20_B.Latelalacc = VCM20_P.Gain7_Gain_e * VCM20_B.SFunction1_o3_fn;

    /* S-Function (rti_commonblock): '<S103>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:256 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X100]);

      /* Convert timestamp */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X100]->processed) {
        can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X100]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X100]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X100]->timestamp > 0.0) {
        if (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X100]->processed) {
          CAN_Msg = can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X100]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "TIME_STAMP" (0|32, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[3];
            VCM20_B.SFunction1_o1_n1 = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "GENERAL_STATUS" (32|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o2_h = ((real_T) CAN_Sgn.UnsignedSgn);

            /* ...... "CLOCK_STATUS" (48|16, standard signal, unsigned int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[7];
            CAN_Sgn.UnsignedSgn &= 0x0000FFFF;
            VCM20_B.SFunction1_o3_ll = ((real_T) CAN_Sgn.UnsignedSgn);
          }
        }
      }
    }

    /* Logic: '<S114>/Logical Operator1' */
    VCM20_B.LogicalOperator1_j = ((VCM20_B.SFunction1_o2_c != 0.0) ||
      (VCM20_B.SFunction1_o2_b != 0.0) || (VCM20_B.SFunction1_o2_e != 0.0) ||
      (VCM20_B.SFunction1_o2_f != 0.0));

    /* DiscretePulseGenerator: '<S114>/Pulse Generator' */
    VCM20_B.PulseGenerator = (VCM20_DW.clockTickCounter <
      VCM20_P.PulseGenerator_Duty) && (VCM20_DW.clockTickCounter >= 0) ?
      VCM20_P.PulseGenerator_Amp : 0.0;

    /* DiscretePulseGenerator: '<S114>/Pulse Generator' */
    if (VCM20_DW.clockTickCounter >= VCM20_P.PulseGenerator_Period - 1.0) {
      VCM20_DW.clockTickCounter = 0;
    } else {
      VCM20_DW.clockTickCounter++;
    }

    /* S-Function (sdspcount2): '<S114>/Counter' */
    VCM20_B.Counter = false;

    /* S-Function (sdspcount2): '<S114>/Counter' */
    if (MWDSP_EPH_R_D(VCM20_B.PulseGenerator, &VCM20_DW.Counter_RstEphState) !=
        0U) {
      VCM20_DW.Counter_Count = VCM20_P.Counter_InitialCount;
    }

    if (MWDSP_EPH_R_B(VCM20_B.LogicalOperator1_j, &VCM20_DW.Counter_ClkEphState)
        != 0U) {
      if (VCM20_DW.Counter_Count < ((uint8_T)5U)) {
        VCM20_DW.Counter_Count++;
      } else {
        VCM20_DW.Counter_Count = 0U;
      }
    }

    /* S-Function (sdspcount2): '<S114>/Counter' */
    VCM20_B.Counter = ((VCM20_DW.Counter_Count == VCM20_P.Counter_HitValue) ||
                       VCM20_B.Counter);

    /* Delay: '<S114>/Delay' */
    VCM20_B.Delay_mz = VCM20_DW.Delay_DSTATE_f;

    /* Logic: '<S114>/Logical Operator4' */
    VCM20_B.LogicalOperator4 = (VCM20_B.Counter || VCM20_B.Delay_mz);

    /* Logic: '<S114>/Logical Operator3' */
    VCM20_B.LogicalOperator3_j = !VCM20_B.LogicalOperator4;

    /* Logic: '<S114>/Logical Operator2' */
    VCM20_B.LogicalOperator2_m = (VCM20_B.LogicalOperator1_j &&
      VCM20_B.LogicalOperator3_j);

    /* Logic: '<S114>/Logical Operator' */
    VCM20_B.LogicalOperator_ij = (VCM20_B.LogicalOperator2_m ||
      VCM20_B.ErrorResetSW_f);
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* Gain: '<S127>/Gain' incorporates:
     *  Constant: '<S15>/Trottle_Pedal bias'
     */
    VCM20_B.Gain_nx = VCM20_P.Gain_Gain_il * VCM20_P.Trottle_Pedalbias_Value;

    /* Sum: '<S127>/Sum1' */
    VCM20_B.Sum1_b = VCM20_B.APPSsig - VCM20_B.Gain_nx;

    /* RelationalOperator: '<S129>/Compare' incorporates:
     *  Constant: '<S129>/Constant'
     */
    VCM20_B.Compare_bz = (VCM20_B.Sum1_b >= VCM20_P.Constant_Value_a);

    /* Logic: '<S127>/Logical Operator3' */
    VCM20_B.LogicalOperator3_d = (VCM20_B.LogicalOperator2 && VCM20_B.Compare_bz);
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* Product: '<S221>/Product5' */
  VCM20_B.Product5_n = VCM20_B.vms * VCM20_B.vms;

  /* Product: '<S221>/Divide3' incorporates:
   *  Constant: '<S216>/stabilityfactor_noDYC'
   */
  VCM20_B.Divide3_a = VCM20_B.Product5_n * VCM20_P.stabilityfactor_noDYC_Value;

  /* Sum: '<S221>/Add4' incorporates:
   *  Constant: '<S221>/Constant3'
   */
  VCM20_B.Add4_g = VCM20_P.Constant3_Value_k + VCM20_B.Divide3_a;

  /* Product: '<S221>/Divide2' */
  VCM20_B.Divide2_m = VCM20_B.vms / VCM20_B.Add4_g;

  /* Gain: '<S221>/Gain6' */
  VCM20_B.Gain6_k = VCM20_P.Gain6_Gain_nl * VCM20_B.Divide2_m;

  /* Product: '<S221>/Divide' */
  VCM20_B.Divide_a2 = VCM20_B.Gain6_k * VCM20_B.SteerAngle_rad;

  /* Sum: '<S216>/Sum2' */
  VCM20_B.Sum2_j = VCM20_B.Divide_a - VCM20_B.Divide_a2;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* TransferFcn: '<S217>/Transfer Fcn1' */
  VCM20_B.TransferFcn1 = 0.0;
  VCM20_B.TransferFcn1 += VCM20_P.TransferFcn1_C * VCM20_X.TransferFcn1_CSTATE;

  /* Sum: '<S217>/Add1' */
  VCM20_B.Add1_f = VCM20_B.Divide_a - VCM20_B.TransferFcn1;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* Gain: '<S249>/Derivative Gain' */
  VCM20_B.DerivativeGain = VCM20_P.PIDController_D * VCM20_B.Add1_f;

  /* Integrator: '<S250>/Filter' */
  VCM20_B.Filter = VCM20_X.Filter_CSTATE;

  /* Sum: '<S250>/SumD' */
  VCM20_B.SumD = VCM20_B.DerivativeGain - VCM20_B.Filter;

  /* Gain: '<S252>/Integral Gain' */
  VCM20_B.IntegralGain = VCM20_P.PIDController_I * VCM20_B.Add1_f;

  /* Integrator: '<S255>/Integrator' */
  VCM20_B.Integrator = VCM20_X.Integrator_CSTATE;

  /* Gain: '<S258>/Filter Coefficient' */
  VCM20_B.FilterCoefficient = VCM20_P.PIDController_N * VCM20_B.SumD;

  /* Gain: '<S260>/Proportional Gain' */
  VCM20_B.ProportionalGain = VCM20_P.PIDController_P * VCM20_B.Add1_f;

  /* Sum: '<S264>/Sum' */
  VCM20_B.Sum_g = (VCM20_B.ProportionalGain + VCM20_B.Integrator) +
    VCM20_B.FilterCoefficient;

  /* Gain: '<S218>/m//s to km//h' */
  VCM20_B.v = VCM20_P.mstokmh_Gain * VCM20_B.vms;

  /* Saturate: '<S218>/Saturation' */
  u0_0 = VCM20_B.v;
  u1 = VCM20_P.Saturation_LowerSat_dn;
  u2 = VCM20_P.Saturation_UpperSat_d4;
  if (u0_0 > u2) {
    /* Saturate: '<S218>/Saturation' */
    VCM20_B.Saturation_oa = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S218>/Saturation' */
    VCM20_B.Saturation_oa = u1;
  } else {
    /* Saturate: '<S218>/Saturation' */
    VCM20_B.Saturation_oa = u0_0;
  }

  /* End of Saturate: '<S218>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S218>/FFGain_0km//h' incorporates:
     *  Concatenate: '<S218>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_k[0] = VCM20_P.FFGain_0kmh_Value;

    /* Constant: '<S218>/FFGain_20km//h' incorporates:
     *  Concatenate: '<S218>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_k[1] = VCM20_P.FFGain_20kmh_Value;

    /* Constant: '<S218>/FFGain_40km//h' incorporates:
     *  Concatenate: '<S218>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_k[2] = VCM20_P.FFGain_40kmh_Value;

    /* Constant: '<S218>/FFGain_60km//h' incorporates:
     *  Concatenate: '<S218>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_k[3] = VCM20_P.FFGain_60kmh_Value;

    /* Constant: '<S218>/FFGain_80km//h' incorporates:
     *  Concatenate: '<S218>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_k[4] = VCM20_P.FFGain_80kmh_Value;
  }

  /* Lookup_n-D: '<S218>/1-D Lookup Table' */
  /*
   * About '<S218>/1-D Lookup Table':
   *       Table size:  5
   *    Interpolation:  Spline
   *    Extrapolation:  Spline
   *   Breakpt Search:  Binary
   *    Breakpt Cache:  OFF
   */
  VCM20_B.uDLookupTable = look_SplNBinSZcd(1U, &VCM20_B.Saturation_oa,
    (rt_LUTSplineWork*)&VCM20_DW.SWork[0]);

  /* Gain: '<S218>/Gain_DYC_FF' */
  VCM20_B.Gain_DYC_FF = VCM20_P.Gain_DYC_FF_Gain * VCM20_B.Sum2_j;

  /* Product: '<S218>/Divide' */
  VCM20_B.Divide_b5 = VCM20_B.Gain_DYC_FF * VCM20_B.uDLookupTable;

  /* Gain: '<S219>/m//s to km//h' */
  VCM20_B.v_b = VCM20_P.mstokmh_Gain_i * VCM20_B.vms;

  /* Saturate: '<S219>/Saturation' */
  u0_0 = VCM20_B.v_b;
  u1 = VCM20_P.Saturation_LowerSat_mf;
  u2 = VCM20_P.Saturation_UpperSat_pi;
  if (u0_0 > u2) {
    /* Saturate: '<S219>/Saturation' */
    VCM20_B.Saturation_k = u2;
  } else if (u0_0 < u1) {
    /* Saturate: '<S219>/Saturation' */
    VCM20_B.Saturation_k = u1;
  } else {
    /* Saturate: '<S219>/Saturation' */
    VCM20_B.Saturation_k = u0_0;
  }

  /* End of Saturate: '<S219>/Saturation' */
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Constant: '<S219>/FFdiffGain_0km//h' incorporates:
     *  Concatenate: '<S219>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_gm[0] = VCM20_P.FFdiffGain_0kmh_Value;

    /* Constant: '<S219>/FFdiffGain_20km//h' incorporates:
     *  Concatenate: '<S219>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_gm[1] = VCM20_P.FFdiffGain_20kmh_Value;

    /* Constant: '<S219>/FFdiffGain_40km//h' incorporates:
     *  Concatenate: '<S219>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_gm[2] = VCM20_P.FFdiffGain_40kmh_Value;

    /* Constant: '<S219>/FFdiffGain_60km//h' incorporates:
     *  Concatenate: '<S219>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_gm[3] = VCM20_P.FFdiffGain_60kmh_Value;

    /* Constant: '<S219>/FFdiffGain_80km//h' incorporates:
     *  Concatenate: '<S219>/Vector Concatenate'
     */
    VCM20_B.VectorConcatenate_gm[4] = VCM20_P.FFdiffGain_80kmh_Value;
  }

  /* Lookup_n-D: '<S219>/1-D Lookup Table' */
  /*
   * About '<S219>/1-D Lookup Table':
   *       Table size:  5
   *    Interpolation:  Spline
   *    Extrapolation:  Spline
   *   Breakpt Search:  Binary
   *    Breakpt Cache:  OFF
   */
  VCM20_B.uDLookupTable_f = look_SplNBinSZcd(1U, &VCM20_B.Saturation_k,
    (rt_LUTSplineWork*)&VCM20_DW.SWork_m[0]);

  /* TransferFcn: '<S219>/Transfer Fcn2' */
  VCM20_B.TransferFcn2 = 0.0;
  VCM20_B.TransferFcn2 += VCM20_P.TransferFcn2_C * VCM20_X.TransferFcn2_CSTATE;
  VCM20_B.TransferFcn2 += VCM20_P.TransferFcn2_D * VCM20_B.Sum2_j;

  /* Gain: '<S219>/I' */
  VCM20_B.I = VCM20_P.I_Gain * VCM20_B.TransferFcn2;

  /* Gain: '<S219>/FF_diff_Gain' */
  VCM20_B.FF_diff_Gain = VCM20_P.FF_diff_Gain_Gain * VCM20_B.I;

  /* Product: '<S219>/Divide' */
  VCM20_B.Divide_m3 = VCM20_B.FF_diff_Gain * VCM20_B.uDLookupTable_f;

  /* Sum: '<S216>/Sum' */
  VCM20_B.Sum_f = VCM20_B.Divide_b5;

  /* Sum: '<S216>/Sum1' */
  VCM20_B.Sum1_d = VCM20_B.Sum_f + VCM20_B.Sum_g;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Math: '<S220>/Square1' */
    VCM20_B.Square1 = VCM20_B.VELOCITY_X_p * VCM20_B.VELOCITY_X_p;

    /* Math: '<S220>/Square2' */
    VCM20_B.Square2 = VCM20_B.VELOCITY_Y * VCM20_B.VELOCITY_Y;

    /* Sum: '<S220>/Add' */
    VCM20_B.Add_f = VCM20_B.Square1 + VCM20_B.Square2;

    /* Sqrt: '<S220>/Sqrt' */
    VCM20_B.vms_f = sqrt(VCM20_B.Add_f);
  }

  /* Gain: '<S175>/Gain' */
  VCM20_B.RR = VCM20_P.Gain_Gain_jv * VCM20_B.Switch2_g;

  /* Gain: '<S175>/Gain1' */
  VCM20_B.RL = VCM20_P.Gain1_Gain_m4 * VCM20_B.Switch2_lc;

  /* Gain: '<S175>/Gain2' */
  VCM20_B.FR = VCM20_P.Gain2_Gain_i * VCM20_B.Switch2_hr;

  /* Gain: '<S175>/Gain3' */
  VCM20_B.FL = VCM20_P.Gain3_Gain_i * VCM20_B.Switch2_ot;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* Gain: '<S291>/I gain' */
  VCM20_B.Igain = VCM20_P.Igain_Gain * VCM20_B.Pgain_i;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* Gain: '<S292>/I gain' */
  VCM20_B.Igain_b = VCM20_P.Igain_Gain_b * VCM20_B.Pgain_a;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* Gain: '<S293>/I gain' */
  VCM20_B.Igain_a = VCM20_P.Igain_Gain_k * VCM20_B.Pgain_m;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
  }

  /* Gain: '<S294>/I gain' */
  VCM20_B.Igain_j = VCM20_P.Igain_Gain_l * VCM20_B.Pgain;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S341>/S-Function1' */
    /* This comment workarounds a code generation problem */
    {
      /* dSPACE I/O Board DS1401STDADCT4 #1 Unit:ADC Group:ADC */
      adc_tp4_single_new_read(ADC_TP4_1_MODULE_ADDR,
        ADC_TP4_CH8,
        (dsfloat *)&VCM20_B.SFunction1_n);
    }

    /* Product: '<S338>/Product' incorporates:
     *  Constant: '<S338>/Constant'
     */
    VCM20_B.Product_n = VCM20_P.Constant_Value_o3 * VCM20_B.SFunction1_n;

    /* Sum: '<S338>/Add' incorporates:
     *  Constant: '<S338>/Constant1'
     */
    VCM20_B.reikyaku_temp1 = VCM20_B.Product_n - VCM20_P.Constant1_Value_m;

    /* S-Function (rti_commonblock): '<S6>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* S-Function (rti_commonblock): '<S7>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* S-Function (rti_commonblock): '<S8>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* S-Function (rti_commonblock): '<S9>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* S-Function (rti_commonblock): '<S10>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* S-Function (rti_commonblock): '<S17>/S-Function1' */
    /* This comment workarounds a code generation problem */
  }
}

/* Model update function */
void VCM20_update(void)
{
  int_T idxDelay;
  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[2] == 0) {
    /* Update for Delay: '<S82>/Delay1' */
    VCM20_DW.Delay1_DSTATE_c = VCM20_B.Compare;

    /* Update for Delay: '<S82>/Delay2' */
    VCM20_DW.Delay2_DSTATE_a[0] = VCM20_DW.Delay2_DSTATE_a[1];
    VCM20_DW.Delay2_DSTATE_a[1] = VCM20_B.Compare;

    /* Update for Delay: '<S82>/Delay3' */
    VCM20_DW.Delay3_DSTATE_p[0] = VCM20_DW.Delay3_DSTATE_p[1];
    VCM20_DW.Delay3_DSTATE_p[1] = VCM20_DW.Delay3_DSTATE_p[2];
    VCM20_DW.Delay3_DSTATE_p[2] = VCM20_B.Compare;

    /* Update for Delay: '<S82>/Delay4' */
    VCM20_DW.Delay4_DSTATE_o[0] = VCM20_DW.Delay4_DSTATE_o[1];
    VCM20_DW.Delay4_DSTATE_o[1] = VCM20_DW.Delay4_DSTATE_o[2];
    VCM20_DW.Delay4_DSTATE_o[2] = VCM20_DW.Delay4_DSTATE_o[3];
    VCM20_DW.Delay4_DSTATE_o[3] = VCM20_B.Compare;

    /* Update for Delay: '<S127>/Delay' */
    VCM20_DW.Delay_DSTATE_d = VCM20_B.LogicalOperator3_d;

    /* Update for Delay: '<S92>/Delay1' */
    VCM20_DW.Delay1_DSTATE_ci = VCM20_B.LogicalOperator_g;

    /* Update for Delay: '<S92>/Delay2' */
    VCM20_DW.Delay2_DSTATE_j[0] = VCM20_DW.Delay2_DSTATE_j[1];
    VCM20_DW.Delay2_DSTATE_j[1] = VCM20_B.LogicalOperator_g;

    /* Update for Delay: '<S92>/Delay3' */
    VCM20_DW.Delay3_DSTATE_g[0] = VCM20_DW.Delay3_DSTATE_g[1];
    VCM20_DW.Delay3_DSTATE_g[1] = VCM20_DW.Delay3_DSTATE_g[2];
    VCM20_DW.Delay3_DSTATE_g[2] = VCM20_B.LogicalOperator_g;

    /* Update for Delay: '<S92>/Delay4' */
    VCM20_DW.Delay4_DSTATE_g[0] = VCM20_DW.Delay4_DSTATE_g[1];
    VCM20_DW.Delay4_DSTATE_g[1] = VCM20_DW.Delay4_DSTATE_g[2];
    VCM20_DW.Delay4_DSTATE_g[2] = VCM20_DW.Delay4_DSTATE_g[3];
    VCM20_DW.Delay4_DSTATE_g[3] = VCM20_B.LogicalOperator_g;

    /* Update for Delay: '<S92>/Delay5' */
    VCM20_DW.Delay5_DSTATE_g[0] = VCM20_DW.Delay5_DSTATE_g[1];
    VCM20_DW.Delay5_DSTATE_g[1] = VCM20_DW.Delay5_DSTATE_g[2];
    VCM20_DW.Delay5_DSTATE_g[2] = VCM20_DW.Delay5_DSTATE_g[3];
    VCM20_DW.Delay5_DSTATE_g[3] = VCM20_DW.Delay5_DSTATE_g[4];
    VCM20_DW.Delay5_DSTATE_g[4] = VCM20_B.LogicalOperator_g;

    /* Update for Delay: '<S92>/Delay6' */
    for (idxDelay = 0; idxDelay < 5; idxDelay++) {
      VCM20_DW.Delay6_DSTATE_lu[idxDelay] = VCM20_DW.Delay6_DSTATE_lu[idxDelay +
        1];
    }

    VCM20_DW.Delay6_DSTATE_lu[5] = VCM20_B.LogicalOperator_g;

    /* End of Update for Delay: '<S92>/Delay6' */

    /* Update for Delay: '<S92>/Delay7' */
    for (idxDelay = 0; idxDelay < 6; idxDelay++) {
      VCM20_DW.Delay7_DSTATE_nq[idxDelay] = VCM20_DW.Delay7_DSTATE_nq[idxDelay +
        1];
    }

    VCM20_DW.Delay7_DSTATE_nq[6] = VCM20_B.LogicalOperator_g;

    /* End of Update for Delay: '<S92>/Delay7' */

    /* Update for Delay: '<S92>/Delay8' */
    for (idxDelay = 0; idxDelay < 7; idxDelay++) {
      VCM20_DW.Delay8_DSTATE_k[idxDelay] = VCM20_DW.Delay8_DSTATE_k[idxDelay + 1];
    }

    VCM20_DW.Delay8_DSTATE_k[7] = VCM20_B.LogicalOperator_g;

    /* End of Update for Delay: '<S92>/Delay8' */

    /* Update for Delay: '<S92>/Delay9' */
    for (idxDelay = 0; idxDelay < 8; idxDelay++) {
      VCM20_DW.Delay9_DSTATE_o[idxDelay] = VCM20_DW.Delay9_DSTATE_o[idxDelay + 1];
    }

    VCM20_DW.Delay9_DSTATE_o[8] = VCM20_B.LogicalOperator_g;

    /* End of Update for Delay: '<S92>/Delay9' */

    /* Update for Delay: '<S79>/Delay1' */
    VCM20_DW.Delay1_DSTATE_p = VCM20_B.LogicalOperator_b;

    /* Update for Delay: '<S79>/Delay2' */
    VCM20_DW.Delay2_DSTATE_e[0] = VCM20_DW.Delay2_DSTATE_e[1];
    VCM20_DW.Delay2_DSTATE_e[1] = VCM20_B.LogicalOperator_b;

    /* Update for Delay: '<S79>/Delay3' */
    VCM20_DW.Delay3_DSTATE_n[0] = VCM20_DW.Delay3_DSTATE_n[1];
    VCM20_DW.Delay3_DSTATE_n[1] = VCM20_DW.Delay3_DSTATE_n[2];
    VCM20_DW.Delay3_DSTATE_n[2] = VCM20_B.LogicalOperator_b;

    /* Update for Delay: '<S79>/Delay4' */
    VCM20_DW.Delay4_DSTATE_j[0] = VCM20_DW.Delay4_DSTATE_j[1];
    VCM20_DW.Delay4_DSTATE_j[1] = VCM20_DW.Delay4_DSTATE_j[2];
    VCM20_DW.Delay4_DSTATE_j[2] = VCM20_DW.Delay4_DSTATE_j[3];
    VCM20_DW.Delay4_DSTATE_j[3] = VCM20_B.LogicalOperator_b;

    /* Update for Delay: '<S79>/Delay5' */
    VCM20_DW.Delay5_DSTATE_d[0] = VCM20_DW.Delay5_DSTATE_d[1];
    VCM20_DW.Delay5_DSTATE_d[1] = VCM20_DW.Delay5_DSTATE_d[2];
    VCM20_DW.Delay5_DSTATE_d[2] = VCM20_DW.Delay5_DSTATE_d[3];
    VCM20_DW.Delay5_DSTATE_d[3] = VCM20_DW.Delay5_DSTATE_d[4];
    VCM20_DW.Delay5_DSTATE_d[4] = VCM20_B.LogicalOperator_b;

    /* Update for Delay: '<S79>/Delay6' */
    for (idxDelay = 0; idxDelay < 5; idxDelay++) {
      VCM20_DW.Delay6_DSTATE_i[idxDelay] = VCM20_DW.Delay6_DSTATE_i[idxDelay + 1];
    }

    VCM20_DW.Delay6_DSTATE_i[5] = VCM20_B.LogicalOperator_b;

    /* End of Update for Delay: '<S79>/Delay6' */

    /* Update for Delay: '<S79>/Delay7' */
    for (idxDelay = 0; idxDelay < 6; idxDelay++) {
      VCM20_DW.Delay7_DSTATE_g[idxDelay] = VCM20_DW.Delay7_DSTATE_g[idxDelay + 1];
    }

    VCM20_DW.Delay7_DSTATE_g[6] = VCM20_B.LogicalOperator_b;

    /* End of Update for Delay: '<S79>/Delay7' */

    /* Update for Delay: '<S79>/Delay8' */
    for (idxDelay = 0; idxDelay < 7; idxDelay++) {
      VCM20_DW.Delay8_DSTATE_c[idxDelay] = VCM20_DW.Delay8_DSTATE_c[idxDelay + 1];
    }

    VCM20_DW.Delay8_DSTATE_c[7] = VCM20_B.LogicalOperator_b;

    /* End of Update for Delay: '<S79>/Delay8' */

    /* Update for Delay: '<S79>/Delay9' */
    for (idxDelay = 0; idxDelay < 8; idxDelay++) {
      VCM20_DW.Delay9_DSTATE_h[idxDelay] = VCM20_DW.Delay9_DSTATE_h[idxDelay + 1];
    }

    VCM20_DW.Delay9_DSTATE_h[8] = VCM20_B.LogicalOperator_b;

    /* End of Update for Delay: '<S79>/Delay9' */

    /* Update for Delay: '<S80>/Delay1' */
    VCM20_DW.Delay1_DSTATE_f = VCM20_B.LogicalOperator2_n;

    /* Update for Delay: '<S80>/Delay2' */
    VCM20_DW.Delay2_DSTATE_k[0] = VCM20_DW.Delay2_DSTATE_k[1];
    VCM20_DW.Delay2_DSTATE_k[1] = VCM20_B.LogicalOperator2_n;

    /* Update for Delay: '<S80>/Delay3' */
    VCM20_DW.Delay3_DSTATE_gb[0] = VCM20_DW.Delay3_DSTATE_gb[1];
    VCM20_DW.Delay3_DSTATE_gb[1] = VCM20_DW.Delay3_DSTATE_gb[2];
    VCM20_DW.Delay3_DSTATE_gb[2] = VCM20_B.LogicalOperator2_n;

    /* Update for Delay: '<S80>/Delay4' */
    VCM20_DW.Delay4_DSTATE_a[0] = VCM20_DW.Delay4_DSTATE_a[1];
    VCM20_DW.Delay4_DSTATE_a[1] = VCM20_DW.Delay4_DSTATE_a[2];
    VCM20_DW.Delay4_DSTATE_a[2] = VCM20_DW.Delay4_DSTATE_a[3];
    VCM20_DW.Delay4_DSTATE_a[3] = VCM20_B.LogicalOperator2_n;

    /* Update for Delay: '<S80>/Delay5' */
    VCM20_DW.Delay5_DSTATE_k[0] = VCM20_DW.Delay5_DSTATE_k[1];
    VCM20_DW.Delay5_DSTATE_k[1] = VCM20_DW.Delay5_DSTATE_k[2];
    VCM20_DW.Delay5_DSTATE_k[2] = VCM20_DW.Delay5_DSTATE_k[3];
    VCM20_DW.Delay5_DSTATE_k[3] = VCM20_DW.Delay5_DSTATE_k[4];
    VCM20_DW.Delay5_DSTATE_k[4] = VCM20_B.LogicalOperator2_n;

    /* Update for Delay: '<S80>/Delay6' */
    for (idxDelay = 0; idxDelay < 5; idxDelay++) {
      VCM20_DW.Delay6_DSTATE_n[idxDelay] = VCM20_DW.Delay6_DSTATE_n[idxDelay + 1];
    }

    VCM20_DW.Delay6_DSTATE_n[5] = VCM20_B.LogicalOperator2_n;

    /* End of Update for Delay: '<S80>/Delay6' */

    /* Update for Delay: '<S80>/Delay7' */
    for (idxDelay = 0; idxDelay < 6; idxDelay++) {
      VCM20_DW.Delay7_DSTATE_m[idxDelay] = VCM20_DW.Delay7_DSTATE_m[idxDelay + 1];
    }

    VCM20_DW.Delay7_DSTATE_m[6] = VCM20_B.LogicalOperator2_n;

    /* End of Update for Delay: '<S80>/Delay7' */

    /* Update for Delay: '<S80>/Delay8' */
    for (idxDelay = 0; idxDelay < 7; idxDelay++) {
      VCM20_DW.Delay8_DSTATE_d[idxDelay] = VCM20_DW.Delay8_DSTATE_d[idxDelay + 1];
    }

    VCM20_DW.Delay8_DSTATE_d[7] = VCM20_B.LogicalOperator2_n;

    /* End of Update for Delay: '<S80>/Delay8' */

    /* Update for Delay: '<S80>/Delay9' */
    for (idxDelay = 0; idxDelay < 8; idxDelay++) {
      VCM20_DW.Delay9_DSTATE_p[idxDelay] = VCM20_DW.Delay9_DSTATE_p[idxDelay + 1];
    }

    VCM20_DW.Delay9_DSTATE_p[8] = VCM20_B.LogicalOperator2_n;

    /* End of Update for Delay: '<S80>/Delay9' */

    /* Update for Delay: '<S89>/Delay11' */
    VCM20_DW.Delay11_DSTATE = VCM20_B.SFunction1_o1_n;

    /* Update for Delay: '<S89>/Delay12' */
    VCM20_DW.Delay12_DSTATE[0] = VCM20_DW.Delay12_DSTATE[1];
    VCM20_DW.Delay12_DSTATE[1] = VCM20_B.SFunction1_o1_n;

    /* Update for Delay: '<S89>/Delay13' */
    VCM20_DW.Delay13_DSTATE[0] = VCM20_DW.Delay13_DSTATE[1];
    VCM20_DW.Delay13_DSTATE[1] = VCM20_DW.Delay13_DSTATE[2];
    VCM20_DW.Delay13_DSTATE[2] = VCM20_B.SFunction1_o1_n;

    /* Update for Delay: '<S89>/Delay14' */
    VCM20_DW.Delay14_DSTATE[0] = VCM20_DW.Delay14_DSTATE[1];
    VCM20_DW.Delay14_DSTATE[1] = VCM20_DW.Delay14_DSTATE[2];
    VCM20_DW.Delay14_DSTATE[2] = VCM20_DW.Delay14_DSTATE[3];
    VCM20_DW.Delay14_DSTATE[3] = VCM20_B.SFunction1_o1_n;

    /* Update for Delay: '<S89>/Delay15' */
    VCM20_DW.Delay15_DSTATE[0] = VCM20_DW.Delay15_DSTATE[1];
    VCM20_DW.Delay15_DSTATE[1] = VCM20_DW.Delay15_DSTATE[2];
    VCM20_DW.Delay15_DSTATE[2] = VCM20_DW.Delay15_DSTATE[3];
    VCM20_DW.Delay15_DSTATE[3] = VCM20_DW.Delay15_DSTATE[4];
    VCM20_DW.Delay15_DSTATE[4] = VCM20_B.SFunction1_o1_n;

    /* Update for Delay: '<S89>/Delay16' */
    for (idxDelay = 0; idxDelay < 5; idxDelay++) {
      VCM20_DW.Delay16_DSTATE[idxDelay] = VCM20_DW.Delay16_DSTATE[idxDelay + 1];
    }

    VCM20_DW.Delay16_DSTATE[5] = VCM20_B.SFunction1_o1_n;

    /* End of Update for Delay: '<S89>/Delay16' */

    /* Update for Delay: '<S89>/Delay17' */
    for (idxDelay = 0; idxDelay < 6; idxDelay++) {
      VCM20_DW.Delay17_DSTATE[idxDelay] = VCM20_DW.Delay17_DSTATE[idxDelay + 1];
    }

    VCM20_DW.Delay17_DSTATE[6] = VCM20_B.SFunction1_o1_n;

    /* End of Update for Delay: '<S89>/Delay17' */

    /* Update for Delay: '<S89>/Delay18' */
    for (idxDelay = 0; idxDelay < 7; idxDelay++) {
      VCM20_DW.Delay18_DSTATE[idxDelay] = VCM20_DW.Delay18_DSTATE[idxDelay + 1];
    }

    VCM20_DW.Delay18_DSTATE[7] = VCM20_B.SFunction1_o1_n;

    /* End of Update for Delay: '<S89>/Delay18' */

    /* Update for Delay: '<S89>/Delay19' */
    for (idxDelay = 0; idxDelay < 8; idxDelay++) {
      VCM20_DW.Delay19_DSTATE[idxDelay] = VCM20_DW.Delay19_DSTATE[idxDelay + 1];
    }

    VCM20_DW.Delay19_DSTATE[8] = VCM20_B.SFunction1_o1_n;

    /* End of Update for Delay: '<S89>/Delay19' */

    /* Update for Delay: '<S89>/Delay1' */
    VCM20_DW.Delay1_DSTATE = VCM20_B.SFunction1_o1_n;

    /* Update for Delay: '<S89>/Delay2' */
    VCM20_DW.Delay2_DSTATE[0] = VCM20_DW.Delay2_DSTATE[1];
    VCM20_DW.Delay2_DSTATE[1] = VCM20_B.SFunction1_o1_n;

    /* Update for Delay: '<S89>/Delay3' */
    VCM20_DW.Delay3_DSTATE[0] = VCM20_DW.Delay3_DSTATE[1];
    VCM20_DW.Delay3_DSTATE[1] = VCM20_DW.Delay3_DSTATE[2];
    VCM20_DW.Delay3_DSTATE[2] = VCM20_B.SFunction1_o1_n;

    /* Update for Delay: '<S89>/Delay4' */
    VCM20_DW.Delay4_DSTATE[0] = VCM20_DW.Delay4_DSTATE[1];
    VCM20_DW.Delay4_DSTATE[1] = VCM20_DW.Delay4_DSTATE[2];
    VCM20_DW.Delay4_DSTATE[2] = VCM20_DW.Delay4_DSTATE[3];
    VCM20_DW.Delay4_DSTATE[3] = VCM20_B.SFunction1_o1_n;

    /* Update for Delay: '<S89>/Delay5' */
    VCM20_DW.Delay5_DSTATE[0] = VCM20_DW.Delay5_DSTATE[1];
    VCM20_DW.Delay5_DSTATE[1] = VCM20_DW.Delay5_DSTATE[2];
    VCM20_DW.Delay5_DSTATE[2] = VCM20_DW.Delay5_DSTATE[3];
    VCM20_DW.Delay5_DSTATE[3] = VCM20_DW.Delay5_DSTATE[4];
    VCM20_DW.Delay5_DSTATE[4] = VCM20_B.SFunction1_o1_n;

    /* Update for Delay: '<S89>/Delay6' */
    for (idxDelay = 0; idxDelay < 5; idxDelay++) {
      VCM20_DW.Delay6_DSTATE[idxDelay] = VCM20_DW.Delay6_DSTATE[idxDelay + 1];
    }

    VCM20_DW.Delay6_DSTATE[5] = VCM20_B.SFunction1_o1_n;

    /* End of Update for Delay: '<S89>/Delay6' */

    /* Update for Delay: '<S89>/Delay7' */
    for (idxDelay = 0; idxDelay < 6; idxDelay++) {
      VCM20_DW.Delay7_DSTATE[idxDelay] = VCM20_DW.Delay7_DSTATE[idxDelay + 1];
    }

    VCM20_DW.Delay7_DSTATE[6] = VCM20_B.SFunction1_o1_n;

    /* End of Update for Delay: '<S89>/Delay7' */

    /* Update for Delay: '<S89>/Delay8' */
    for (idxDelay = 0; idxDelay < 7; idxDelay++) {
      VCM20_DW.Delay8_DSTATE[idxDelay] = VCM20_DW.Delay8_DSTATE[idxDelay + 1];
    }

    VCM20_DW.Delay8_DSTATE[7] = VCM20_B.SFunction1_o1_n;

    /* End of Update for Delay: '<S89>/Delay8' */

    /* Update for Delay: '<S89>/Delay9' */
    for (idxDelay = 0; idxDelay < 8; idxDelay++) {
      VCM20_DW.Delay9_DSTATE[idxDelay] = VCM20_DW.Delay9_DSTATE[idxDelay + 1];
    }

    VCM20_DW.Delay9_DSTATE[8] = VCM20_B.SFunction1_o1_n;

    /* End of Update for Delay: '<S89>/Delay9' */

    /* Update for Delay: '<S90>/Delay11' */
    VCM20_DW.Delay11_DSTATE_o = VCM20_B.SFunction1_o2_j;

    /* Update for Delay: '<S90>/Delay12' */
    VCM20_DW.Delay12_DSTATE_i[0] = VCM20_DW.Delay12_DSTATE_i[1];
    VCM20_DW.Delay12_DSTATE_i[1] = VCM20_B.SFunction1_o2_j;

    /* Update for Delay: '<S90>/Delay13' */
    VCM20_DW.Delay13_DSTATE_p[0] = VCM20_DW.Delay13_DSTATE_p[1];
    VCM20_DW.Delay13_DSTATE_p[1] = VCM20_DW.Delay13_DSTATE_p[2];
    VCM20_DW.Delay13_DSTATE_p[2] = VCM20_B.SFunction1_o2_j;

    /* Update for Delay: '<S90>/Delay14' */
    VCM20_DW.Delay14_DSTATE_k[0] = VCM20_DW.Delay14_DSTATE_k[1];
    VCM20_DW.Delay14_DSTATE_k[1] = VCM20_DW.Delay14_DSTATE_k[2];
    VCM20_DW.Delay14_DSTATE_k[2] = VCM20_DW.Delay14_DSTATE_k[3];
    VCM20_DW.Delay14_DSTATE_k[3] = VCM20_B.SFunction1_o2_j;

    /* Update for Delay: '<S90>/Delay15' */
    VCM20_DW.Delay15_DSTATE_n[0] = VCM20_DW.Delay15_DSTATE_n[1];
    VCM20_DW.Delay15_DSTATE_n[1] = VCM20_DW.Delay15_DSTATE_n[2];
    VCM20_DW.Delay15_DSTATE_n[2] = VCM20_DW.Delay15_DSTATE_n[3];
    VCM20_DW.Delay15_DSTATE_n[3] = VCM20_DW.Delay15_DSTATE_n[4];
    VCM20_DW.Delay15_DSTATE_n[4] = VCM20_B.SFunction1_o2_j;

    /* Update for Delay: '<S90>/Delay16' */
    for (idxDelay = 0; idxDelay < 5; idxDelay++) {
      VCM20_DW.Delay16_DSTATE_c[idxDelay] = VCM20_DW.Delay16_DSTATE_c[idxDelay +
        1];
    }

    VCM20_DW.Delay16_DSTATE_c[5] = VCM20_B.SFunction1_o2_j;

    /* End of Update for Delay: '<S90>/Delay16' */

    /* Update for Delay: '<S90>/Delay17' */
    for (idxDelay = 0; idxDelay < 6; idxDelay++) {
      VCM20_DW.Delay17_DSTATE_e[idxDelay] = VCM20_DW.Delay17_DSTATE_e[idxDelay +
        1];
    }

    VCM20_DW.Delay17_DSTATE_e[6] = VCM20_B.SFunction1_o2_j;

    /* End of Update for Delay: '<S90>/Delay17' */

    /* Update for Delay: '<S90>/Delay18' */
    for (idxDelay = 0; idxDelay < 7; idxDelay++) {
      VCM20_DW.Delay18_DSTATE_a[idxDelay] = VCM20_DW.Delay18_DSTATE_a[idxDelay +
        1];
    }

    VCM20_DW.Delay18_DSTATE_a[7] = VCM20_B.SFunction1_o2_j;

    /* End of Update for Delay: '<S90>/Delay18' */

    /* Update for Delay: '<S90>/Delay19' */
    for (idxDelay = 0; idxDelay < 8; idxDelay++) {
      VCM20_DW.Delay19_DSTATE_h[idxDelay] = VCM20_DW.Delay19_DSTATE_h[idxDelay +
        1];
    }

    VCM20_DW.Delay19_DSTATE_h[8] = VCM20_B.SFunction1_o2_j;

    /* End of Update for Delay: '<S90>/Delay19' */

    /* Update for Delay: '<S90>/Delay1' */
    VCM20_DW.Delay1_DSTATE_i = VCM20_B.SFunction1_o2_j;

    /* Update for Delay: '<S90>/Delay2' */
    VCM20_DW.Delay2_DSTATE_p[0] = VCM20_DW.Delay2_DSTATE_p[1];
    VCM20_DW.Delay2_DSTATE_p[1] = VCM20_B.SFunction1_o2_j;

    /* Update for Delay: '<S90>/Delay3' */
    VCM20_DW.Delay3_DSTATE_o[0] = VCM20_DW.Delay3_DSTATE_o[1];
    VCM20_DW.Delay3_DSTATE_o[1] = VCM20_DW.Delay3_DSTATE_o[2];
    VCM20_DW.Delay3_DSTATE_o[2] = VCM20_B.SFunction1_o2_j;

    /* Update for Delay: '<S90>/Delay4' */
    VCM20_DW.Delay4_DSTATE_i[0] = VCM20_DW.Delay4_DSTATE_i[1];
    VCM20_DW.Delay4_DSTATE_i[1] = VCM20_DW.Delay4_DSTATE_i[2];
    VCM20_DW.Delay4_DSTATE_i[2] = VCM20_DW.Delay4_DSTATE_i[3];
    VCM20_DW.Delay4_DSTATE_i[3] = VCM20_B.SFunction1_o2_j;

    /* Update for Delay: '<S90>/Delay5' */
    VCM20_DW.Delay5_DSTATE_m[0] = VCM20_DW.Delay5_DSTATE_m[1];
    VCM20_DW.Delay5_DSTATE_m[1] = VCM20_DW.Delay5_DSTATE_m[2];
    VCM20_DW.Delay5_DSTATE_m[2] = VCM20_DW.Delay5_DSTATE_m[3];
    VCM20_DW.Delay5_DSTATE_m[3] = VCM20_DW.Delay5_DSTATE_m[4];
    VCM20_DW.Delay5_DSTATE_m[4] = VCM20_B.SFunction1_o2_j;

    /* Update for Delay: '<S90>/Delay6' */
    for (idxDelay = 0; idxDelay < 5; idxDelay++) {
      VCM20_DW.Delay6_DSTATE_p[idxDelay] = VCM20_DW.Delay6_DSTATE_p[idxDelay + 1];
    }

    VCM20_DW.Delay6_DSTATE_p[5] = VCM20_B.SFunction1_o2_j;

    /* End of Update for Delay: '<S90>/Delay6' */

    /* Update for Delay: '<S90>/Delay7' */
    for (idxDelay = 0; idxDelay < 6; idxDelay++) {
      VCM20_DW.Delay7_DSTATE_k[idxDelay] = VCM20_DW.Delay7_DSTATE_k[idxDelay + 1];
    }

    VCM20_DW.Delay7_DSTATE_k[6] = VCM20_B.SFunction1_o2_j;

    /* End of Update for Delay: '<S90>/Delay7' */

    /* Update for Delay: '<S90>/Delay8' */
    for (idxDelay = 0; idxDelay < 7; idxDelay++) {
      VCM20_DW.Delay8_DSTATE_j[idxDelay] = VCM20_DW.Delay8_DSTATE_j[idxDelay + 1];
    }

    VCM20_DW.Delay8_DSTATE_j[7] = VCM20_B.SFunction1_o2_j;

    /* End of Update for Delay: '<S90>/Delay8' */

    /* Update for Delay: '<S90>/Delay9' */
    for (idxDelay = 0; idxDelay < 8; idxDelay++) {
      VCM20_DW.Delay9_DSTATE_a[idxDelay] = VCM20_DW.Delay9_DSTATE_a[idxDelay + 1];
    }

    VCM20_DW.Delay9_DSTATE_a[8] = VCM20_B.SFunction1_o2_j;

    /* End of Update for Delay: '<S90>/Delay9' */

    /* Update for Delay: '<S91>/Delay11' */
    VCM20_DW.Delay11_DSTATE_k = VCM20_B.SFunction1_o2_a;

    /* Update for Delay: '<S91>/Delay12' */
    VCM20_DW.Delay12_DSTATE_h[0] = VCM20_DW.Delay12_DSTATE_h[1];
    VCM20_DW.Delay12_DSTATE_h[1] = VCM20_B.SFunction1_o2_a;

    /* Update for Delay: '<S91>/Delay13' */
    VCM20_DW.Delay13_DSTATE_j[0] = VCM20_DW.Delay13_DSTATE_j[1];
    VCM20_DW.Delay13_DSTATE_j[1] = VCM20_DW.Delay13_DSTATE_j[2];
    VCM20_DW.Delay13_DSTATE_j[2] = VCM20_B.SFunction1_o2_a;

    /* Update for Delay: '<S91>/Delay14' */
    VCM20_DW.Delay14_DSTATE_l[0] = VCM20_DW.Delay14_DSTATE_l[1];
    VCM20_DW.Delay14_DSTATE_l[1] = VCM20_DW.Delay14_DSTATE_l[2];
    VCM20_DW.Delay14_DSTATE_l[2] = VCM20_DW.Delay14_DSTATE_l[3];
    VCM20_DW.Delay14_DSTATE_l[3] = VCM20_B.SFunction1_o2_a;

    /* Update for Delay: '<S91>/Delay15' */
    VCM20_DW.Delay15_DSTATE_i[0] = VCM20_DW.Delay15_DSTATE_i[1];
    VCM20_DW.Delay15_DSTATE_i[1] = VCM20_DW.Delay15_DSTATE_i[2];
    VCM20_DW.Delay15_DSTATE_i[2] = VCM20_DW.Delay15_DSTATE_i[3];
    VCM20_DW.Delay15_DSTATE_i[3] = VCM20_DW.Delay15_DSTATE_i[4];
    VCM20_DW.Delay15_DSTATE_i[4] = VCM20_B.SFunction1_o2_a;

    /* Update for Delay: '<S91>/Delay16' */
    for (idxDelay = 0; idxDelay < 5; idxDelay++) {
      VCM20_DW.Delay16_DSTATE_b[idxDelay] = VCM20_DW.Delay16_DSTATE_b[idxDelay +
        1];
    }

    VCM20_DW.Delay16_DSTATE_b[5] = VCM20_B.SFunction1_o2_a;

    /* End of Update for Delay: '<S91>/Delay16' */

    /* Update for Delay: '<S91>/Delay17' */
    for (idxDelay = 0; idxDelay < 6; idxDelay++) {
      VCM20_DW.Delay17_DSTATE_n[idxDelay] = VCM20_DW.Delay17_DSTATE_n[idxDelay +
        1];
    }

    VCM20_DW.Delay17_DSTATE_n[6] = VCM20_B.SFunction1_o2_a;

    /* End of Update for Delay: '<S91>/Delay17' */

    /* Update for Delay: '<S91>/Delay18' */
    for (idxDelay = 0; idxDelay < 7; idxDelay++) {
      VCM20_DW.Delay18_DSTATE_h[idxDelay] = VCM20_DW.Delay18_DSTATE_h[idxDelay +
        1];
    }

    VCM20_DW.Delay18_DSTATE_h[7] = VCM20_B.SFunction1_o2_a;

    /* End of Update for Delay: '<S91>/Delay18' */

    /* Update for Delay: '<S91>/Delay19' */
    for (idxDelay = 0; idxDelay < 8; idxDelay++) {
      VCM20_DW.Delay19_DSTATE_n[idxDelay] = VCM20_DW.Delay19_DSTATE_n[idxDelay +
        1];
    }

    VCM20_DW.Delay19_DSTATE_n[8] = VCM20_B.SFunction1_o2_a;

    /* End of Update for Delay: '<S91>/Delay19' */

    /* Update for Delay: '<S91>/Delay1' */
    VCM20_DW.Delay1_DSTATE_h = VCM20_B.SFunction1_o2_a;

    /* Update for Delay: '<S91>/Delay2' */
    VCM20_DW.Delay2_DSTATE_l[0] = VCM20_DW.Delay2_DSTATE_l[1];
    VCM20_DW.Delay2_DSTATE_l[1] = VCM20_B.SFunction1_o2_a;

    /* Update for Delay: '<S91>/Delay3' */
    VCM20_DW.Delay3_DSTATE_o3[0] = VCM20_DW.Delay3_DSTATE_o3[1];
    VCM20_DW.Delay3_DSTATE_o3[1] = VCM20_DW.Delay3_DSTATE_o3[2];
    VCM20_DW.Delay3_DSTATE_o3[2] = VCM20_B.SFunction1_o2_a;

    /* Update for Delay: '<S91>/Delay4' */
    VCM20_DW.Delay4_DSTATE_il[0] = VCM20_DW.Delay4_DSTATE_il[1];
    VCM20_DW.Delay4_DSTATE_il[1] = VCM20_DW.Delay4_DSTATE_il[2];
    VCM20_DW.Delay4_DSTATE_il[2] = VCM20_DW.Delay4_DSTATE_il[3];
    VCM20_DW.Delay4_DSTATE_il[3] = VCM20_B.SFunction1_o2_a;

    /* Update for Delay: '<S91>/Delay5' */
    VCM20_DW.Delay5_DSTATE_mb[0] = VCM20_DW.Delay5_DSTATE_mb[1];
    VCM20_DW.Delay5_DSTATE_mb[1] = VCM20_DW.Delay5_DSTATE_mb[2];
    VCM20_DW.Delay5_DSTATE_mb[2] = VCM20_DW.Delay5_DSTATE_mb[3];
    VCM20_DW.Delay5_DSTATE_mb[3] = VCM20_DW.Delay5_DSTATE_mb[4];
    VCM20_DW.Delay5_DSTATE_mb[4] = VCM20_B.SFunction1_o2_a;

    /* Update for Delay: '<S91>/Delay6' */
    for (idxDelay = 0; idxDelay < 5; idxDelay++) {
      VCM20_DW.Delay6_DSTATE_l[idxDelay] = VCM20_DW.Delay6_DSTATE_l[idxDelay + 1];
    }

    VCM20_DW.Delay6_DSTATE_l[5] = VCM20_B.SFunction1_o2_a;

    /* End of Update for Delay: '<S91>/Delay6' */

    /* Update for Delay: '<S91>/Delay7' */
    for (idxDelay = 0; idxDelay < 6; idxDelay++) {
      VCM20_DW.Delay7_DSTATE_n[idxDelay] = VCM20_DW.Delay7_DSTATE_n[idxDelay + 1];
    }

    VCM20_DW.Delay7_DSTATE_n[6] = VCM20_B.SFunction1_o2_a;

    /* End of Update for Delay: '<S91>/Delay7' */

    /* Update for Delay: '<S91>/Delay8' */
    for (idxDelay = 0; idxDelay < 7; idxDelay++) {
      VCM20_DW.Delay8_DSTATE_jc[idxDelay] = VCM20_DW.Delay8_DSTATE_jc[idxDelay +
        1];
    }

    VCM20_DW.Delay8_DSTATE_jc[7] = VCM20_B.SFunction1_o2_a;

    /* End of Update for Delay: '<S91>/Delay8' */

    /* Update for Delay: '<S91>/Delay9' */
    for (idxDelay = 0; idxDelay < 8; idxDelay++) {
      VCM20_DW.Delay9_DSTATE_e[idxDelay] = VCM20_DW.Delay9_DSTATE_e[idxDelay + 1];
    }

    VCM20_DW.Delay9_DSTATE_e[8] = VCM20_B.SFunction1_o2_a;

    /* End of Update for Delay: '<S91>/Delay9' */
  }

  if (rtmIsMajorTimeStep(VCM20_M) &&
      VCM20_M->Timing.TaskCounters.TID[1] == 0) {
    /* Update for Delay: '<S296>/Delay' incorporates:
     *  Constant: '<S283>/Target Acc'
     */
    VCM20_DW.Delay_DSTATE = VCM20_P.g;

    /* Update for DiscreteIntegrator: '<S294>/Discrete-Time Integrator' */
    VCM20_DW.DiscreteTimeIntegrator_DSTATE +=
      VCM20_P.DiscreteTimeIntegrator_gainval * VCM20_B.Igain_j;

    /* Update for DiscreteIntegrator: '<S291>/Discrete-Time Integrator' */
    VCM20_DW.DiscreteTimeIntegrator_DSTATE_e +=
      VCM20_P.DiscreteTimeIntegrator_gainva_j * VCM20_B.Igain;

    /* Update for DiscreteIntegrator: '<S292>/Discrete-Time Integrator' */
    VCM20_DW.DiscreteTimeIntegrator_DSTATE_m +=
      VCM20_P.DiscreteTimeIntegrator_gainva_m * VCM20_B.Igain_b;

    /* Update for DiscreteIntegrator: '<S293>/Discrete-Time Integrator' */
    VCM20_DW.DiscreteTimeIntegrator_DSTATE_o +=
      VCM20_P.DiscreteTimeIntegrator_gainv_ji * VCM20_B.Igain_a;

    /* Update for Delay: '<S138>/Delay' */
    VCM20_DW.Delay_DSTATE_h = VCM20_B.Switch_as;

    /* Update for Delay: '<S140>/Delay' */
    VCM20_DW.Delay_DSTATE_k = VCM20_B.Quantizer1;

    /* Update for Delay: '<S139>/Delay' */
    VCM20_DW.Delay_DSTATE_b = VCM20_B.Quantizer2;

    /* Update for Delay: '<S146>/Delay' */
    VCM20_DW.Delay_DSTATE_e = VCM20_B.Switch_kd;

    /* Update for Delay: '<S148>/Delay' */
    VCM20_DW.Delay_DSTATE_p = VCM20_B.Quantizer1_p;

    /* Update for Delay: '<S147>/Delay' */
    VCM20_DW.Delay_DSTATE_i = VCM20_B.Quantizer2_p;

    /* Update for Delay: '<S154>/Delay' */
    VCM20_DW.Delay_DSTATE_bc = VCM20_B.Switch_k1;

    /* Update for Delay: '<S156>/Delay' */
    VCM20_DW.Delay_DSTATE_kj = VCM20_B.Quantizer1_m;

    /* Update for Delay: '<S155>/Delay' */
    VCM20_DW.Delay_DSTATE_hs = VCM20_B.Quantizer2_n;

    /* Update for Delay: '<S162>/Delay' */
    VCM20_DW.Delay_DSTATE_pw = VCM20_B.Switch_f4;

    /* Update for Delay: '<S164>/Delay' */
    VCM20_DW.Delay_DSTATE_hw = VCM20_B.Quantizer1_n;

    /* Update for Delay: '<S163>/Delay' */
    VCM20_DW.Delay_DSTATE_m = VCM20_B.Quantizer2_e;

    /* Update for Delay: '<S114>/Delay' */
    VCM20_DW.Delay_DSTATE_f = VCM20_B.LogicalOperator4;
  }

  if (rtmIsMajorTimeStep(VCM20_M)) {
    rt_ertODEUpdateContinuousStates(&VCM20_M->solverInfo);
  }

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++VCM20_M->Timing.clockTick0)) {
    ++VCM20_M->Timing.clockTickH0;
  }

  VCM20_M->Timing.t[0] = rtsiGetSolverStopTime(&VCM20_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.001s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The resolution of this integer timer is 0.001, which is the step size
     * of the task. Size of "clockTick1" ensures timer will not overflow during the
     * application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    VCM20_M->Timing.clockTick1++;
    if (!VCM20_M->Timing.clockTick1) {
      VCM20_M->Timing.clockTickH1++;
    }
  }

  rate_scheduler();
}

/* Derivatives for root system: '<Root>' */
void VCM20_derivatives(void)
{
  XDot_VCM20_T *_rtXdot;
  _rtXdot = ((XDot_VCM20_T *) VCM20_M->derivs);

  /* Derivatives for TransferFcn: '<S174>/Transfer Fcn' */
  _rtXdot->TransferFcn_CSTATE = 0.0;
  _rtXdot->TransferFcn_CSTATE += VCM20_P.TransferFcn_A *
    VCM20_X.TransferFcn_CSTATE;
  _rtXdot->TransferFcn_CSTATE += VCM20_B.ACCEL_X_m;

  /* Derivatives for TransferFcn: '<S73>/Transfer Fcn' */
  _rtXdot->TransferFcn_CSTATE_k = 0.0;
  _rtXdot->TransferFcn_CSTATE_k += VCM20_P.TransferFcn_A_c *
    VCM20_X.TransferFcn_CSTATE_k;
  _rtXdot->TransferFcn_CSTATE_k += VCM20_B.SFunction1_o4_h;

  /* Derivatives for TransferFcn: '<S220>/Transfer Fcn' */
  _rtXdot->TransferFcn_CSTATE_j = 0.0;
  _rtXdot->TransferFcn_CSTATE_j += VCM20_P.TransferFcn_A_e *
    VCM20_X.TransferFcn_CSTATE_j;
  _rtXdot->TransferFcn_CSTATE_j += VCM20_B.vms_f;

  /* Derivatives for TransferFcn: '<S72>/Transfer Fcn' */
  _rtXdot->TransferFcn_CSTATE_e = 0.0;
  _rtXdot->TransferFcn_CSTATE_e += VCM20_P.TransferFcn_A_a *
    VCM20_X.TransferFcn_CSTATE_e;
  _rtXdot->TransferFcn_CSTATE_e += VCM20_B.SFunction1_o3_d2;

  /* Derivatives for TransferFcn: '<S217>/Transfer Fcn1' */
  _rtXdot->TransferFcn1_CSTATE = 0.0;
  _rtXdot->TransferFcn1_CSTATE += VCM20_P.TransferFcn1_A *
    VCM20_X.TransferFcn1_CSTATE;
  _rtXdot->TransferFcn1_CSTATE += VCM20_B.GYRO_Z;

  /* Derivatives for Integrator: '<S250>/Filter' */
  _rtXdot->Filter_CSTATE = VCM20_B.FilterCoefficient;

  /* Derivatives for Integrator: '<S255>/Integrator' */
  _rtXdot->Integrator_CSTATE = VCM20_B.IntegralGain;

  /* Derivatives for TransferFcn: '<S219>/Transfer Fcn2' */
  _rtXdot->TransferFcn2_CSTATE = 0.0;
  _rtXdot->TransferFcn2_CSTATE += VCM20_P.TransferFcn2_A *
    VCM20_X.TransferFcn2_CSTATE;
  _rtXdot->TransferFcn2_CSTATE += VCM20_B.Sum2_j;
}

/* Model initialize function */
void VCM20_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* non-finite (run-time) assignments */
  VCM20_P.Saturation1_LowerSat = rtMinusInf;
  VCM20_P.Saturation2_LowerSat_g = rtMinusInf;
  VCM20_P.Saturation3_LowerSat = rtMinusInf;
  VCM20_P.Saturation_LowerSat_n = rtMinusInf;
  VCM20_P.Saturation3_UpperSat_m = rtInf;
  VCM20_P.Saturation_UpperSat_c = rtInf;
  VCM20_P.Saturation3_LowerSat_oy = rtMinusInf;
  VCM20_P.Saturation_LowerSat_p = rtMinusInf;
  VCM20_P.Saturation_UpperSat_g = rtInf;
  VCM20_P.Saturation_UpperSat_fe = rtInf;
  VCM20_P.Saturation_UpperSat_ow = rtInf;
  VCM20_P.Saturation_UpperSat_d = rtInf;
  VCM20_P.Saturation_UpperSat_lo = rtInf;

  /* initialize real-time model */
  (void) memset((void *)VCM20_M, 0,
                sizeof(RT_MODEL_VCM20_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&VCM20_M->solverInfo, &VCM20_M->Timing.simTimeStep);
    rtsiSetTPtr(&VCM20_M->solverInfo, &rtmGetTPtr(VCM20_M));
    rtsiSetStepSizePtr(&VCM20_M->solverInfo, &VCM20_M->Timing.stepSize0);
    rtsiSetdXPtr(&VCM20_M->solverInfo, &VCM20_M->derivs);
    rtsiSetContStatesPtr(&VCM20_M->solverInfo, (real_T **) &VCM20_M->contStates);
    rtsiSetNumContStatesPtr(&VCM20_M->solverInfo, &VCM20_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&VCM20_M->solverInfo,
      &VCM20_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&VCM20_M->solverInfo,
      &VCM20_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&VCM20_M->solverInfo,
      &VCM20_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&VCM20_M->solverInfo, (&rtmGetErrorStatus(VCM20_M)));
    rtsiSetRTModelPtr(&VCM20_M->solverInfo, VCM20_M);
  }

  rtsiSetSimTimeStep(&VCM20_M->solverInfo, MAJOR_TIME_STEP);
  VCM20_M->intgData.f[0] = VCM20_M->odeF[0];
  VCM20_M->contStates = ((X_VCM20_T *) &VCM20_X);
  rtsiSetSolverData(&VCM20_M->solverInfo, (void *)&VCM20_M->intgData);
  rtsiSetSolverName(&VCM20_M->solverInfo,"ode1");
  rtmSetTPtr(VCM20_M, &VCM20_M->Timing.tArray[0]);
  VCM20_M->Timing.stepSize0 = 0.001;

  /* block I/O */
  (void) memset(((void *) &VCM20_B), 0,
                sizeof(B_VCM20_T));

  /* states (continuous) */
  {
    (void) memset((void *)&VCM20_X, 0,
                  sizeof(X_VCM20_T));
  }

  /* states (dwork) */
  (void) memset((void *)&VCM20_DW, 0,
                sizeof(DW_VCM20_T));

  {
    /* user code (registration function declaration) */
    /*Initialize global TRC pointers. */
    VCM20_rti_init_trc_pointers();
  }

  {
    dsp_simulink_MovingAverage_b_T *b_obj;
    dsp_simulink_MovingAverage_b_T *obj;
    g_dsp_internal_SlidingWindo_b_T *iobj_0;

    /* Start for MATLABSystem: '<S81>/Moving Average' */
    VCM20_DW.obj.matlabCodegenIsDeleted = true;
    b_obj = &VCM20_DW.obj;
    b_obj->isInitialized = 0;
    b_obj->NumChannels = -1;
    b_obj->matlabCodegenIsDeleted = false;
    VCM20_DW.objisempty = true;
    b_obj = &VCM20_DW.obj;
    b_obj->isSetupComplete = false;
    b_obj->isInitialized = 1;
    obj = b_obj;
    obj->NumChannels = 1;
    iobj_0 = &obj->_pobj0;
    iobj_0->isInitialized = 0;
    iobj_0->isInitialized = 0;
    obj->pStatistic = iobj_0;
    b_obj->isSetupComplete = true;
    b_obj->TunablePropsChanged = false;
    VCM20_MovingAverage_Start(&VCM20_DW.MovingAverage_p);
    VCM20_MovingAverage_Start(&VCM20_DW.MovingAverage_pn);
    VCM20_MovingAverage_Start(&VCM20_DW.MovingAverage_pna);
    VCM20_MovingAverage_Start(&VCM20_DW.MovingAverage_pnae);
    VCM20_MovingAverage_Start(&VCM20_DW.MovingAverage_pnaev);

    /* Start for DiscretePulseGenerator: '<S114>/Pulse Generator' */
    VCM20_DW.clockTickCounter = 0;

    /* Start for Lookup_n-D: '<S218>/1-D Lookup Table' */
    {
      rt_LUTnWork *TWork_start = (rt_LUTnWork *) &VCM20_DW.TWork[0];
      void **bpDataSet = (void **) &VCM20_DW.m_bpDataSet;
      TWork_start->m_dimSizes = (const uint32_T *)
        &VCM20_P.uDLookupTable_dimSizes;
      TWork_start->m_tableData = (void *) &VCM20_B.VectorConcatenate_k[0];
      TWork_start->m_bpDataSet = bpDataSet;
      TWork_start->m_bpIndex = &VCM20_DW.m_bpIndex;
      TWork_start->m_bpLambda = &VCM20_DW.m_bpLambda;
      TWork_start->m_maxIndex = (const uint32_T *)
        &VCM20_P.uDLookupTable_maxIndex_e;
      bpDataSet[0] = (void *) VCM20_P.uDLookupTable_bp01Data_h3;
    }

    {
      const real_T **bpDataSet;
      const real_T *xp, *yp;
      real_T *dp;
      uint32_T len;
      const rt_LUTnWork *TWork_interp;
      rt_LUTSplineWork *rt_SplWk = (rt_LUTSplineWork*)&VCM20_DW.SWork[0];
      rt_SplWk->m_TWork = (rt_LUTnWork*)&VCM20_DW.TWork[0];
      rt_SplWk->m_yyA = &VCM20_DW.m_yyA;
      rt_SplWk->m_yyB = &VCM20_DW.m_yyB;
      rt_SplWk->m_yy2 = &VCM20_DW.m_yy2;
      rt_SplWk->m_up = &VCM20_DW.m_up[0];
      rt_SplWk->m_y2 = &VCM20_DW.m_y2[0];
      rt_SplWk->m_numYWorkElts = VCM20_P.uDLookupTable_numYWorkElts;
      rt_SplWk->m_reCalc = &VCM20_DW.reCalcSecDerivFirstDimCoeffs;
      rt_SplWk->m_preBp0AndTable = &VCM20_DW.prevBp0AndTableData[0];
      *rt_SplWk->m_reCalc = 1;

      /* cache table data and first breakpoint data */
      TWork_interp = (const rt_LUTnWork *)rt_SplWk->m_TWork;
      bpDataSet = (const real_T **) TWork_interp->m_bpDataSet;
      xp = bpDataSet[0U];
      len = TWork_interp->m_maxIndex[0U] + 1U;
      dp = (real_T *) rt_SplWk->m_preBp0AndTable;
      yp = (real_T *) TWork_interp->m_tableData;
      (void) memcpy(dp, xp,
                    len * sizeof(real_T));
      dp = &(dp[len]);

      /* save the table data */
      (void) memcpy(dp, yp,
                    len * rt_SplWk->m_numYWorkElts[0U] * sizeof(real_T));
    }

    /* Start for Lookup_n-D: '<S219>/1-D Lookup Table' */
    {
      rt_LUTnWork *TWork_start = (rt_LUTnWork *) &VCM20_DW.TWork_c[0];
      void **bpDataSet = (void **) &VCM20_DW.m_bpDataSet_o;
      TWork_start->m_dimSizes = (const uint32_T *)
        &VCM20_P.uDLookupTable_dimSizes_e;
      TWork_start->m_tableData = (void *) &VCM20_B.VectorConcatenate_gm[0];
      TWork_start->m_bpDataSet = bpDataSet;
      TWork_start->m_bpIndex = &VCM20_DW.m_bpIndex_h;
      TWork_start->m_bpLambda = &VCM20_DW.m_bpLambda_h;
      TWork_start->m_maxIndex = (const uint32_T *)
        &VCM20_P.uDLookupTable_maxIndex_oo;
      bpDataSet[0] = (void *) VCM20_P.uDLookupTable_bp01Data_lf;
    }

    {
      const real_T **bpDataSet;
      const real_T *xp, *yp;
      real_T *dp;
      uint32_T len;
      const rt_LUTnWork *TWork_interp;
      rt_LUTSplineWork *rt_SplWk = (rt_LUTSplineWork*)&VCM20_DW.SWork_m[0];
      rt_SplWk->m_TWork = (rt_LUTnWork*)&VCM20_DW.TWork_c[0];
      rt_SplWk->m_yyA = &VCM20_DW.m_yyA_j;
      rt_SplWk->m_yyB = &VCM20_DW.m_yyB_m;
      rt_SplWk->m_yy2 = &VCM20_DW.m_yy2_p;
      rt_SplWk->m_up = &VCM20_DW.m_up_m[0];
      rt_SplWk->m_y2 = &VCM20_DW.m_y2_k[0];
      rt_SplWk->m_numYWorkElts = VCM20_P.uDLookupTable_numYWorkElts_i;
      rt_SplWk->m_reCalc = &VCM20_DW.reCalcSecDerivFirstDimCoeffs_i;
      rt_SplWk->m_preBp0AndTable = &VCM20_DW.prevBp0AndTableData_m[0];
      *rt_SplWk->m_reCalc = 1;

      /* cache table data and first breakpoint data */
      TWork_interp = (const rt_LUTnWork *)rt_SplWk->m_TWork;
      bpDataSet = (const real_T **) TWork_interp->m_bpDataSet;
      xp = bpDataSet[0U];
      len = TWork_interp->m_maxIndex[0U] + 1U;
      dp = (real_T *) rt_SplWk->m_preBp0AndTable;
      yp = (real_T *) TWork_interp->m_tableData;
      (void) memcpy(dp, xp,
                    len * sizeof(real_T));
      dp = &(dp[len]);

      /* save the table data */
      (void) memcpy(dp, yp,
                    len * rt_SplWk->m_numYWorkElts[0U] * sizeof(real_T));
    }
  }

  {
    dsp_simulink_MovingAverage_b_T *obj;
    g_dsp_internal_SlidingWindo_b_T *obj_0;
    int32_T i;

    /* InitializeConditions for Delay: '<S82>/Delay1' */
    VCM20_DW.Delay1_DSTATE_c = VCM20_P.Delay1_InitialCondition_e;

    /* InitializeConditions for Delay: '<S82>/Delay2' */
    VCM20_DW.Delay2_DSTATE_a[0] = VCM20_P.Delay2_InitialCondition_hk;
    VCM20_DW.Delay2_DSTATE_a[1] = VCM20_P.Delay2_InitialCondition_hk;

    /* InitializeConditions for Delay: '<S82>/Delay3' */
    VCM20_DW.Delay3_DSTATE_p[0] = VCM20_P.Delay3_InitialCondition_c;
    VCM20_DW.Delay3_DSTATE_p[1] = VCM20_P.Delay3_InitialCondition_c;
    VCM20_DW.Delay3_DSTATE_p[2] = VCM20_P.Delay3_InitialCondition_c;

    /* InitializeConditions for Delay: '<S82>/Delay4' */
    VCM20_DW.Delay4_DSTATE_o[0] = VCM20_P.Delay4_InitialCondition_l;
    VCM20_DW.Delay4_DSTATE_o[1] = VCM20_P.Delay4_InitialCondition_l;
    VCM20_DW.Delay4_DSTATE_o[2] = VCM20_P.Delay4_InitialCondition_l;
    VCM20_DW.Delay4_DSTATE_o[3] = VCM20_P.Delay4_InitialCondition_l;

    /* InitializeConditions for Delay: '<S127>/Delay' */
    VCM20_DW.Delay_DSTATE_d = VCM20_P.Delay_InitialCondition_lq;

    /* InitializeConditions for TransferFcn: '<S174>/Transfer Fcn' */
    VCM20_X.TransferFcn_CSTATE = 0.0;

    /* InitializeConditions for TransferFcn: '<S73>/Transfer Fcn' */
    VCM20_X.TransferFcn_CSTATE_k = 0.0;

    /* InitializeConditions for Delay: '<S92>/Delay1' */
    VCM20_DW.Delay1_DSTATE_ci = VCM20_P.Delay1_InitialCondition_j;

    /* InitializeConditions for Delay: '<S92>/Delay2' */
    VCM20_DW.Delay2_DSTATE_j[0] = VCM20_P.Delay2_InitialCondition_n;
    VCM20_DW.Delay2_DSTATE_j[1] = VCM20_P.Delay2_InitialCondition_n;

    /* InitializeConditions for Delay: '<S92>/Delay3' */
    VCM20_DW.Delay3_DSTATE_g[0] = VCM20_P.Delay3_InitialCondition_b;
    VCM20_DW.Delay3_DSTATE_g[1] = VCM20_P.Delay3_InitialCondition_b;
    VCM20_DW.Delay3_DSTATE_g[2] = VCM20_P.Delay3_InitialCondition_b;

    /* InitializeConditions for Delay: '<S92>/Delay4' */
    VCM20_DW.Delay4_DSTATE_g[0] = VCM20_P.Delay4_InitialCondition_a;
    VCM20_DW.Delay4_DSTATE_g[1] = VCM20_P.Delay4_InitialCondition_a;
    VCM20_DW.Delay4_DSTATE_g[2] = VCM20_P.Delay4_InitialCondition_a;
    VCM20_DW.Delay4_DSTATE_g[3] = VCM20_P.Delay4_InitialCondition_a;

    /* InitializeConditions for Delay: '<S92>/Delay5' */
    for (i = 0; i < 5; i++) {
      VCM20_DW.Delay5_DSTATE_g[i] = VCM20_P.Delay5_InitialCondition_b;
    }

    /* End of InitializeConditions for Delay: '<S92>/Delay5' */

    /* InitializeConditions for Delay: '<S92>/Delay6' */
    for (i = 0; i < 6; i++) {
      VCM20_DW.Delay6_DSTATE_lu[i] = VCM20_P.Delay6_InitialCondition_i;
    }

    /* End of InitializeConditions for Delay: '<S92>/Delay6' */

    /* InitializeConditions for Delay: '<S92>/Delay7' */
    for (i = 0; i < 7; i++) {
      VCM20_DW.Delay7_DSTATE_nq[i] = VCM20_P.Delay7_InitialCondition_h;
    }

    /* End of InitializeConditions for Delay: '<S92>/Delay7' */

    /* InitializeConditions for Delay: '<S92>/Delay8' */
    for (i = 0; i < 8; i++) {
      VCM20_DW.Delay8_DSTATE_k[i] = VCM20_P.Delay8_InitialCondition_o;
    }

    /* End of InitializeConditions for Delay: '<S92>/Delay8' */

    /* InitializeConditions for Delay: '<S92>/Delay9' */
    for (i = 0; i < 9; i++) {
      VCM20_DW.Delay9_DSTATE_o[i] = VCM20_P.Delay9_InitialCondition_f;
    }

    /* End of InitializeConditions for Delay: '<S92>/Delay9' */

    /* InitializeConditions for Delay: '<S296>/Delay' */
    VCM20_DW.Delay_DSTATE = VCM20_P.Delay_InitialCondition;

    /* InitializeConditions for DiscreteIntegrator: '<S294>/Discrete-Time Integrator' */
    VCM20_DW.DiscreteTimeIntegrator_DSTATE = VCM20_P.DiscreteTimeIntegrator_IC;

    /* InitializeConditions for DiscreteIntegrator: '<S291>/Discrete-Time Integrator' */
    VCM20_DW.DiscreteTimeIntegrator_DSTATE_e =
      VCM20_P.DiscreteTimeIntegrator_IC_d;

    /* InitializeConditions for DiscreteIntegrator: '<S292>/Discrete-Time Integrator' */
    VCM20_DW.DiscreteTimeIntegrator_DSTATE_m =
      VCM20_P.DiscreteTimeIntegrator_IC_a;

    /* InitializeConditions for DiscreteIntegrator: '<S293>/Discrete-Time Integrator' */
    VCM20_DW.DiscreteTimeIntegrator_DSTATE_o =
      VCM20_P.DiscreteTimeIntegrator_IC_o;

    /* InitializeConditions for TransferFcn: '<S220>/Transfer Fcn' */
    VCM20_X.TransferFcn_CSTATE_j = 0.0;

    /* InitializeConditions for Delay: '<S138>/Delay' */
    VCM20_DW.Delay_DSTATE_h = VCM20_P.Delay_InitialCondition_i;

    /* InitializeConditions for Delay: '<S140>/Delay' */
    VCM20_DW.Delay_DSTATE_k = VCM20_P.Delay_InitialCondition_p;

    /* InitializeConditions for Delay: '<S139>/Delay' */
    VCM20_DW.Delay_DSTATE_b = VCM20_P.Delay_InitialCondition_l;

    /* InitializeConditions for Delay: '<S146>/Delay' */
    VCM20_DW.Delay_DSTATE_e = VCM20_P.Delay_InitialCondition_n;

    /* InitializeConditions for Delay: '<S148>/Delay' */
    VCM20_DW.Delay_DSTATE_p = VCM20_P.Delay_InitialCondition_g;

    /* InitializeConditions for Delay: '<S147>/Delay' */
    VCM20_DW.Delay_DSTATE_i = VCM20_P.Delay_InitialCondition_j;

    /* InitializeConditions for Delay: '<S154>/Delay' */
    VCM20_DW.Delay_DSTATE_bc = VCM20_P.Delay_InitialCondition_m;

    /* InitializeConditions for Delay: '<S156>/Delay' */
    VCM20_DW.Delay_DSTATE_kj = VCM20_P.Delay_InitialCondition_m0;

    /* InitializeConditions for Delay: '<S155>/Delay' */
    VCM20_DW.Delay_DSTATE_hs = VCM20_P.Delay_InitialCondition_i3;

    /* InitializeConditions for Delay: '<S162>/Delay' */
    VCM20_DW.Delay_DSTATE_pw = VCM20_P.Delay_InitialCondition_mz;

    /* InitializeConditions for Delay: '<S164>/Delay' */
    VCM20_DW.Delay_DSTATE_hw = VCM20_P.Delay_InitialCondition_p2;

    /* InitializeConditions for Delay: '<S163>/Delay' */
    VCM20_DW.Delay_DSTATE_m = VCM20_P.Delay_InitialCondition_a;

    /* InitializeConditions for TransferFcn: '<S72>/Transfer Fcn' */
    VCM20_X.TransferFcn_CSTATE_e = 0.0;

    /* InitializeConditions for Delay: '<S79>/Delay1' */
    VCM20_DW.Delay1_DSTATE_p = VCM20_P.Delay1_InitialCondition_n;

    /* InitializeConditions for Delay: '<S79>/Delay2' */
    VCM20_DW.Delay2_DSTATE_e[0] = VCM20_P.Delay2_InitialCondition_nz;
    VCM20_DW.Delay2_DSTATE_e[1] = VCM20_P.Delay2_InitialCondition_nz;

    /* InitializeConditions for Delay: '<S79>/Delay3' */
    VCM20_DW.Delay3_DSTATE_n[0] = VCM20_P.Delay3_InitialCondition_j;
    VCM20_DW.Delay3_DSTATE_n[1] = VCM20_P.Delay3_InitialCondition_j;
    VCM20_DW.Delay3_DSTATE_n[2] = VCM20_P.Delay3_InitialCondition_j;

    /* InitializeConditions for Delay: '<S79>/Delay4' */
    VCM20_DW.Delay4_DSTATE_j[0] = VCM20_P.Delay4_InitialCondition_al;
    VCM20_DW.Delay4_DSTATE_j[1] = VCM20_P.Delay4_InitialCondition_al;
    VCM20_DW.Delay4_DSTATE_j[2] = VCM20_P.Delay4_InitialCondition_al;
    VCM20_DW.Delay4_DSTATE_j[3] = VCM20_P.Delay4_InitialCondition_al;

    /* InitializeConditions for Delay: '<S79>/Delay5' */
    for (i = 0; i < 5; i++) {
      VCM20_DW.Delay5_DSTATE_d[i] = VCM20_P.Delay5_InitialCondition_m;
    }

    /* End of InitializeConditions for Delay: '<S79>/Delay5' */

    /* InitializeConditions for Delay: '<S79>/Delay6' */
    for (i = 0; i < 6; i++) {
      VCM20_DW.Delay6_DSTATE_i[i] = VCM20_P.Delay6_InitialCondition_d;
    }

    /* End of InitializeConditions for Delay: '<S79>/Delay6' */

    /* InitializeConditions for Delay: '<S79>/Delay7' */
    for (i = 0; i < 7; i++) {
      VCM20_DW.Delay7_DSTATE_g[i] = VCM20_P.Delay7_InitialCondition_f;
    }

    /* End of InitializeConditions for Delay: '<S79>/Delay7' */

    /* InitializeConditions for Delay: '<S79>/Delay8' */
    for (i = 0; i < 8; i++) {
      VCM20_DW.Delay8_DSTATE_c[i] = VCM20_P.Delay8_InitialCondition_m;
    }

    /* End of InitializeConditions for Delay: '<S79>/Delay8' */

    /* InitializeConditions for Delay: '<S79>/Delay9' */
    for (i = 0; i < 9; i++) {
      VCM20_DW.Delay9_DSTATE_h[i] = VCM20_P.Delay9_InitialCondition_fa;
    }

    /* End of InitializeConditions for Delay: '<S79>/Delay9' */

    /* InitializeConditions for Delay: '<S80>/Delay1' */
    VCM20_DW.Delay1_DSTATE_f = VCM20_P.Delay1_InitialCondition_l;

    /* InitializeConditions for Delay: '<S80>/Delay2' */
    VCM20_DW.Delay2_DSTATE_k[0] = VCM20_P.Delay2_InitialCondition_j;
    VCM20_DW.Delay2_DSTATE_k[1] = VCM20_P.Delay2_InitialCondition_j;

    /* InitializeConditions for Delay: '<S80>/Delay3' */
    VCM20_DW.Delay3_DSTATE_gb[0] = VCM20_P.Delay3_InitialCondition_i;
    VCM20_DW.Delay3_DSTATE_gb[1] = VCM20_P.Delay3_InitialCondition_i;
    VCM20_DW.Delay3_DSTATE_gb[2] = VCM20_P.Delay3_InitialCondition_i;

    /* InitializeConditions for Delay: '<S80>/Delay4' */
    VCM20_DW.Delay4_DSTATE_a[0] = VCM20_P.Delay4_InitialCondition_i;
    VCM20_DW.Delay4_DSTATE_a[1] = VCM20_P.Delay4_InitialCondition_i;
    VCM20_DW.Delay4_DSTATE_a[2] = VCM20_P.Delay4_InitialCondition_i;
    VCM20_DW.Delay4_DSTATE_a[3] = VCM20_P.Delay4_InitialCondition_i;

    /* InitializeConditions for Delay: '<S80>/Delay5' */
    for (i = 0; i < 5; i++) {
      VCM20_DW.Delay5_DSTATE_k[i] = VCM20_P.Delay5_InitialCondition_e;
    }

    /* End of InitializeConditions for Delay: '<S80>/Delay5' */

    /* InitializeConditions for Delay: '<S80>/Delay6' */
    for (i = 0; i < 6; i++) {
      VCM20_DW.Delay6_DSTATE_n[i] = VCM20_P.Delay6_InitialCondition_k;
    }

    /* End of InitializeConditions for Delay: '<S80>/Delay6' */

    /* InitializeConditions for Delay: '<S80>/Delay7' */
    for (i = 0; i < 7; i++) {
      VCM20_DW.Delay7_DSTATE_m[i] = VCM20_P.Delay7_InitialCondition_mc;
    }

    /* End of InitializeConditions for Delay: '<S80>/Delay7' */

    /* InitializeConditions for Delay: '<S80>/Delay8' */
    for (i = 0; i < 8; i++) {
      VCM20_DW.Delay8_DSTATE_d[i] = VCM20_P.Delay8_InitialCondition_i;
    }

    /* End of InitializeConditions for Delay: '<S80>/Delay8' */

    /* InitializeConditions for Delay: '<S80>/Delay9' */
    for (i = 0; i < 9; i++) {
      VCM20_DW.Delay9_DSTATE_p[i] = VCM20_P.Delay9_InitialCondition_hf;
    }

    /* End of InitializeConditions for Delay: '<S80>/Delay9' */

    /* InitializeConditions for Delay: '<S89>/Delay11' */
    VCM20_DW.Delay11_DSTATE = VCM20_P.Delay11_InitialCondition;

    /* InitializeConditions for Delay: '<S89>/Delay12' */
    VCM20_DW.Delay12_DSTATE[0] = VCM20_P.Delay12_InitialCondition;
    VCM20_DW.Delay12_DSTATE[1] = VCM20_P.Delay12_InitialCondition;

    /* InitializeConditions for Delay: '<S89>/Delay13' */
    VCM20_DW.Delay13_DSTATE[0] = VCM20_P.Delay13_InitialCondition;
    VCM20_DW.Delay13_DSTATE[1] = VCM20_P.Delay13_InitialCondition;
    VCM20_DW.Delay13_DSTATE[2] = VCM20_P.Delay13_InitialCondition;

    /* InitializeConditions for Delay: '<S89>/Delay14' */
    VCM20_DW.Delay14_DSTATE[0] = VCM20_P.Delay14_InitialCondition;
    VCM20_DW.Delay14_DSTATE[1] = VCM20_P.Delay14_InitialCondition;
    VCM20_DW.Delay14_DSTATE[2] = VCM20_P.Delay14_InitialCondition;
    VCM20_DW.Delay14_DSTATE[3] = VCM20_P.Delay14_InitialCondition;

    /* InitializeConditions for Delay: '<S89>/Delay15' */
    for (i = 0; i < 5; i++) {
      VCM20_DW.Delay15_DSTATE[i] = VCM20_P.Delay15_InitialCondition;
    }

    /* End of InitializeConditions for Delay: '<S89>/Delay15' */

    /* InitializeConditions for Delay: '<S89>/Delay16' */
    for (i = 0; i < 6; i++) {
      VCM20_DW.Delay16_DSTATE[i] = VCM20_P.Delay16_InitialCondition;
    }

    /* End of InitializeConditions for Delay: '<S89>/Delay16' */

    /* InitializeConditions for Delay: '<S89>/Delay17' */
    for (i = 0; i < 7; i++) {
      VCM20_DW.Delay17_DSTATE[i] = VCM20_P.Delay17_InitialCondition;
    }

    /* End of InitializeConditions for Delay: '<S89>/Delay17' */

    /* InitializeConditions for Delay: '<S89>/Delay18' */
    for (i = 0; i < 8; i++) {
      VCM20_DW.Delay18_DSTATE[i] = VCM20_P.Delay18_InitialCondition;
    }

    /* End of InitializeConditions for Delay: '<S89>/Delay18' */

    /* InitializeConditions for Delay: '<S89>/Delay19' */
    for (i = 0; i < 9; i++) {
      VCM20_DW.Delay19_DSTATE[i] = VCM20_P.Delay19_InitialCondition;
    }

    /* End of InitializeConditions for Delay: '<S89>/Delay19' */

    /* InitializeConditions for Delay: '<S89>/Delay1' */
    VCM20_DW.Delay1_DSTATE = VCM20_P.Delay1_InitialCondition;

    /* InitializeConditions for Delay: '<S89>/Delay2' */
    VCM20_DW.Delay2_DSTATE[0] = VCM20_P.Delay2_InitialCondition;
    VCM20_DW.Delay2_DSTATE[1] = VCM20_P.Delay2_InitialCondition;

    /* InitializeConditions for Delay: '<S89>/Delay3' */
    VCM20_DW.Delay3_DSTATE[0] = VCM20_P.Delay3_InitialCondition;
    VCM20_DW.Delay3_DSTATE[1] = VCM20_P.Delay3_InitialCondition;
    VCM20_DW.Delay3_DSTATE[2] = VCM20_P.Delay3_InitialCondition;

    /* InitializeConditions for Delay: '<S89>/Delay4' */
    VCM20_DW.Delay4_DSTATE[0] = VCM20_P.Delay4_InitialCondition;
    VCM20_DW.Delay4_DSTATE[1] = VCM20_P.Delay4_InitialCondition;
    VCM20_DW.Delay4_DSTATE[2] = VCM20_P.Delay4_InitialCondition;
    VCM20_DW.Delay4_DSTATE[3] = VCM20_P.Delay4_InitialCondition;

    /* InitializeConditions for Delay: '<S89>/Delay5' */
    for (i = 0; i < 5; i++) {
      VCM20_DW.Delay5_DSTATE[i] = VCM20_P.Delay5_InitialCondition;
    }

    /* End of InitializeConditions for Delay: '<S89>/Delay5' */

    /* InitializeConditions for Delay: '<S89>/Delay6' */
    for (i = 0; i < 6; i++) {
      VCM20_DW.Delay6_DSTATE[i] = VCM20_P.Delay6_InitialCondition;
    }

    /* End of InitializeConditions for Delay: '<S89>/Delay6' */

    /* InitializeConditions for Delay: '<S89>/Delay7' */
    for (i = 0; i < 7; i++) {
      VCM20_DW.Delay7_DSTATE[i] = VCM20_P.Delay7_InitialCondition;
    }

    /* End of InitializeConditions for Delay: '<S89>/Delay7' */

    /* InitializeConditions for Delay: '<S89>/Delay8' */
    for (i = 0; i < 8; i++) {
      VCM20_DW.Delay8_DSTATE[i] = VCM20_P.Delay8_InitialCondition;
    }

    /* End of InitializeConditions for Delay: '<S89>/Delay8' */

    /* InitializeConditions for Delay: '<S89>/Delay9' */
    for (i = 0; i < 9; i++) {
      VCM20_DW.Delay9_DSTATE[i] = VCM20_P.Delay9_InitialCondition;
    }

    /* End of InitializeConditions for Delay: '<S89>/Delay9' */

    /* InitializeConditions for Delay: '<S90>/Delay11' */
    VCM20_DW.Delay11_DSTATE_o = VCM20_P.Delay11_InitialCondition_n;

    /* InitializeConditions for Delay: '<S90>/Delay12' */
    VCM20_DW.Delay12_DSTATE_i[0] = VCM20_P.Delay12_InitialCondition_d;
    VCM20_DW.Delay12_DSTATE_i[1] = VCM20_P.Delay12_InitialCondition_d;

    /* InitializeConditions for Delay: '<S90>/Delay13' */
    VCM20_DW.Delay13_DSTATE_p[0] = VCM20_P.Delay13_InitialCondition_f;
    VCM20_DW.Delay13_DSTATE_p[1] = VCM20_P.Delay13_InitialCondition_f;
    VCM20_DW.Delay13_DSTATE_p[2] = VCM20_P.Delay13_InitialCondition_f;

    /* InitializeConditions for Delay: '<S90>/Delay14' */
    VCM20_DW.Delay14_DSTATE_k[0] = VCM20_P.Delay14_InitialCondition_e;
    VCM20_DW.Delay14_DSTATE_k[1] = VCM20_P.Delay14_InitialCondition_e;
    VCM20_DW.Delay14_DSTATE_k[2] = VCM20_P.Delay14_InitialCondition_e;
    VCM20_DW.Delay14_DSTATE_k[3] = VCM20_P.Delay14_InitialCondition_e;

    /* InitializeConditions for Delay: '<S90>/Delay15' */
    for (i = 0; i < 5; i++) {
      VCM20_DW.Delay15_DSTATE_n[i] = VCM20_P.Delay15_InitialCondition_l;
    }

    /* End of InitializeConditions for Delay: '<S90>/Delay15' */

    /* InitializeConditions for Delay: '<S90>/Delay16' */
    for (i = 0; i < 6; i++) {
      VCM20_DW.Delay16_DSTATE_c[i] = VCM20_P.Delay16_InitialCondition_g;
    }

    /* End of InitializeConditions for Delay: '<S90>/Delay16' */

    /* InitializeConditions for Delay: '<S90>/Delay17' */
    for (i = 0; i < 7; i++) {
      VCM20_DW.Delay17_DSTATE_e[i] = VCM20_P.Delay17_InitialCondition_a;
    }

    /* End of InitializeConditions for Delay: '<S90>/Delay17' */

    /* InitializeConditions for Delay: '<S90>/Delay18' */
    for (i = 0; i < 8; i++) {
      VCM20_DW.Delay18_DSTATE_a[i] = VCM20_P.Delay18_InitialCondition_j;
    }

    /* End of InitializeConditions for Delay: '<S90>/Delay18' */

    /* InitializeConditions for Delay: '<S90>/Delay19' */
    for (i = 0; i < 9; i++) {
      VCM20_DW.Delay19_DSTATE_h[i] = VCM20_P.Delay19_InitialCondition_c;
    }

    /* End of InitializeConditions for Delay: '<S90>/Delay19' */

    /* InitializeConditions for Delay: '<S90>/Delay1' */
    VCM20_DW.Delay1_DSTATE_i = VCM20_P.Delay1_InitialCondition_g;

    /* InitializeConditions for Delay: '<S90>/Delay2' */
    VCM20_DW.Delay2_DSTATE_p[0] = VCM20_P.Delay2_InitialCondition_m;
    VCM20_DW.Delay2_DSTATE_p[1] = VCM20_P.Delay2_InitialCondition_m;

    /* InitializeConditions for Delay: '<S90>/Delay3' */
    VCM20_DW.Delay3_DSTATE_o[0] = VCM20_P.Delay3_InitialCondition_p;
    VCM20_DW.Delay3_DSTATE_o[1] = VCM20_P.Delay3_InitialCondition_p;
    VCM20_DW.Delay3_DSTATE_o[2] = VCM20_P.Delay3_InitialCondition_p;

    /* InitializeConditions for Delay: '<S90>/Delay4' */
    VCM20_DW.Delay4_DSTATE_i[0] = VCM20_P.Delay4_InitialCondition_j;
    VCM20_DW.Delay4_DSTATE_i[1] = VCM20_P.Delay4_InitialCondition_j;
    VCM20_DW.Delay4_DSTATE_i[2] = VCM20_P.Delay4_InitialCondition_j;
    VCM20_DW.Delay4_DSTATE_i[3] = VCM20_P.Delay4_InitialCondition_j;

    /* InitializeConditions for Delay: '<S90>/Delay5' */
    for (i = 0; i < 5; i++) {
      VCM20_DW.Delay5_DSTATE_m[i] = VCM20_P.Delay5_InitialCondition_c;
    }

    /* End of InitializeConditions for Delay: '<S90>/Delay5' */

    /* InitializeConditions for Delay: '<S90>/Delay6' */
    for (i = 0; i < 6; i++) {
      VCM20_DW.Delay6_DSTATE_p[i] = VCM20_P.Delay6_InitialCondition_e;
    }

    /* End of InitializeConditions for Delay: '<S90>/Delay6' */

    /* InitializeConditions for Delay: '<S90>/Delay7' */
    for (i = 0; i < 7; i++) {
      VCM20_DW.Delay7_DSTATE_k[i] = VCM20_P.Delay7_InitialCondition_e;
    }

    /* End of InitializeConditions for Delay: '<S90>/Delay7' */

    /* InitializeConditions for Delay: '<S90>/Delay8' */
    for (i = 0; i < 8; i++) {
      VCM20_DW.Delay8_DSTATE_j[i] = VCM20_P.Delay8_InitialCondition_f;
    }

    /* End of InitializeConditions for Delay: '<S90>/Delay8' */

    /* InitializeConditions for Delay: '<S90>/Delay9' */
    for (i = 0; i < 9; i++) {
      VCM20_DW.Delay9_DSTATE_a[i] = VCM20_P.Delay9_InitialCondition_h;
    }

    /* End of InitializeConditions for Delay: '<S90>/Delay9' */

    /* InitializeConditions for Delay: '<S91>/Delay11' */
    VCM20_DW.Delay11_DSTATE_k = VCM20_P.Delay11_InitialCondition_e;

    /* InitializeConditions for Delay: '<S91>/Delay12' */
    VCM20_DW.Delay12_DSTATE_h[0] = VCM20_P.Delay12_InitialCondition_o;
    VCM20_DW.Delay12_DSTATE_h[1] = VCM20_P.Delay12_InitialCondition_o;

    /* InitializeConditions for Delay: '<S91>/Delay13' */
    VCM20_DW.Delay13_DSTATE_j[0] = VCM20_P.Delay13_InitialCondition_g;
    VCM20_DW.Delay13_DSTATE_j[1] = VCM20_P.Delay13_InitialCondition_g;
    VCM20_DW.Delay13_DSTATE_j[2] = VCM20_P.Delay13_InitialCondition_g;

    /* InitializeConditions for Delay: '<S91>/Delay14' */
    VCM20_DW.Delay14_DSTATE_l[0] = VCM20_P.Delay14_InitialCondition_g;
    VCM20_DW.Delay14_DSTATE_l[1] = VCM20_P.Delay14_InitialCondition_g;
    VCM20_DW.Delay14_DSTATE_l[2] = VCM20_P.Delay14_InitialCondition_g;
    VCM20_DW.Delay14_DSTATE_l[3] = VCM20_P.Delay14_InitialCondition_g;

    /* InitializeConditions for Delay: '<S91>/Delay15' */
    for (i = 0; i < 5; i++) {
      VCM20_DW.Delay15_DSTATE_i[i] = VCM20_P.Delay15_InitialCondition_o;
    }

    /* End of InitializeConditions for Delay: '<S91>/Delay15' */

    /* InitializeConditions for Delay: '<S91>/Delay16' */
    for (i = 0; i < 6; i++) {
      VCM20_DW.Delay16_DSTATE_b[i] = VCM20_P.Delay16_InitialCondition_o;
    }

    /* End of InitializeConditions for Delay: '<S91>/Delay16' */

    /* InitializeConditions for Delay: '<S91>/Delay17' */
    for (i = 0; i < 7; i++) {
      VCM20_DW.Delay17_DSTATE_n[i] = VCM20_P.Delay17_InitialCondition_o;
    }

    /* End of InitializeConditions for Delay: '<S91>/Delay17' */

    /* InitializeConditions for Delay: '<S91>/Delay18' */
    for (i = 0; i < 8; i++) {
      VCM20_DW.Delay18_DSTATE_h[i] = VCM20_P.Delay18_InitialCondition_d;
    }

    /* End of InitializeConditions for Delay: '<S91>/Delay18' */

    /* InitializeConditions for Delay: '<S91>/Delay19' */
    for (i = 0; i < 9; i++) {
      VCM20_DW.Delay19_DSTATE_n[i] = VCM20_P.Delay19_InitialCondition_a;
    }

    /* End of InitializeConditions for Delay: '<S91>/Delay19' */

    /* InitializeConditions for Delay: '<S91>/Delay1' */
    VCM20_DW.Delay1_DSTATE_h = VCM20_P.Delay1_InitialCondition_k;

    /* InitializeConditions for Delay: '<S91>/Delay2' */
    VCM20_DW.Delay2_DSTATE_l[0] = VCM20_P.Delay2_InitialCondition_h;
    VCM20_DW.Delay2_DSTATE_l[1] = VCM20_P.Delay2_InitialCondition_h;

    /* InitializeConditions for Delay: '<S91>/Delay3' */
    VCM20_DW.Delay3_DSTATE_o3[0] = VCM20_P.Delay3_InitialCondition_o;
    VCM20_DW.Delay3_DSTATE_o3[1] = VCM20_P.Delay3_InitialCondition_o;
    VCM20_DW.Delay3_DSTATE_o3[2] = VCM20_P.Delay3_InitialCondition_o;

    /* InitializeConditions for Delay: '<S91>/Delay4' */
    VCM20_DW.Delay4_DSTATE_il[0] = VCM20_P.Delay4_InitialCondition_m;
    VCM20_DW.Delay4_DSTATE_il[1] = VCM20_P.Delay4_InitialCondition_m;
    VCM20_DW.Delay4_DSTATE_il[2] = VCM20_P.Delay4_InitialCondition_m;
    VCM20_DW.Delay4_DSTATE_il[3] = VCM20_P.Delay4_InitialCondition_m;

    /* InitializeConditions for Delay: '<S91>/Delay5' */
    for (i = 0; i < 5; i++) {
      VCM20_DW.Delay5_DSTATE_mb[i] = VCM20_P.Delay5_InitialCondition_k;
    }

    /* End of InitializeConditions for Delay: '<S91>/Delay5' */

    /* InitializeConditions for Delay: '<S91>/Delay6' */
    for (i = 0; i < 6; i++) {
      VCM20_DW.Delay6_DSTATE_l[i] = VCM20_P.Delay6_InitialCondition_g;
    }

    /* End of InitializeConditions for Delay: '<S91>/Delay6' */

    /* InitializeConditions for Delay: '<S91>/Delay7' */
    for (i = 0; i < 7; i++) {
      VCM20_DW.Delay7_DSTATE_n[i] = VCM20_P.Delay7_InitialCondition_m;
    }

    /* End of InitializeConditions for Delay: '<S91>/Delay7' */

    /* InitializeConditions for Delay: '<S91>/Delay8' */
    for (i = 0; i < 8; i++) {
      VCM20_DW.Delay8_DSTATE_jc[i] = VCM20_P.Delay8_InitialCondition_k;
    }

    /* End of InitializeConditions for Delay: '<S91>/Delay8' */

    /* InitializeConditions for Delay: '<S91>/Delay9' */
    for (i = 0; i < 9; i++) {
      VCM20_DW.Delay9_DSTATE_e[i] = VCM20_P.Delay9_InitialCondition_e;
    }

    /* End of InitializeConditions for Delay: '<S91>/Delay9' */

    /* InitializeConditions for S-Function (sdspcount2): '<S114>/Counter' */
    VCM20_DW.Counter_ClkEphState = 5U;
    VCM20_DW.Counter_RstEphState = 5U;
    VCM20_DW.Counter_Count = VCM20_P.Counter_InitialCount;

    /* InitializeConditions for Delay: '<S114>/Delay' */
    VCM20_DW.Delay_DSTATE_f = VCM20_P.Delay_InitialCondition_e;

    /* InitializeConditions for TransferFcn: '<S217>/Transfer Fcn1' */
    VCM20_X.TransferFcn1_CSTATE = 0.0;

    /* InitializeConditions for Integrator: '<S250>/Filter' */
    VCM20_X.Filter_CSTATE = VCM20_P.PIDController_InitialConditionF;

    /* InitializeConditions for Integrator: '<S255>/Integrator' */
    VCM20_X.Integrator_CSTATE = VCM20_P.PIDController_InitialConditio_l;

    /* InitializeConditions for TransferFcn: '<S219>/Transfer Fcn2' */
    VCM20_X.TransferFcn2_CSTATE = 0.0;

    /* SystemInitialize for Chart: '<S172>/Chart' */
    VCM20_DW.is_LaunchReady_BrakeOff = VCM20_IN_NO_ACTIVE_CHILD;
    VCM20_DW.temporalCounter_i1 = 0U;
    VCM20_DW.is_active_c1_VCM20 = 0U;
    VCM20_DW.is_c1_VCM20 = VCM20_IN_NO_ACTIVE_CHILD;
    VCM20_B.LaunchLED = 0.0;
    VCM20_B.LaunchTorqueLimit = 0.0;

    /* SystemInitialize for Chart: '<S16>/Steer Chart' */
    VCM20_DW.is_active_c3_VCM20 = 0U;
    VCM20_DW.is_c3_VCM20 = VCM20_IN_NO_ACTIVE_CHILD;
    VCM20_DW.gain1max = 0.0;
    VCM20_DW.gain1min = 0.0;
    VCM20_DW.gain2max = 0.0;
    VCM20_DW.gain2min = 0.0;
    VCM20_B.LCDnumber = 0.0;
    for (i = 0; i < 8; i++) {
      VCM20_B.LCDtext[i] = 0U;
    }

    VCM20_B.gain1 = 0.0;
    VCM20_B.gain2 = 0.0;

    /* End of SystemInitialize for Chart: '<S16>/Steer Chart' */

    /* InitializeConditions for MATLABSystem: '<S81>/Moving Average' */
    obj = &VCM20_DW.obj;
    obj_0 = obj->pStatistic;
    if (obj_0->isInitialized == 1) {
      obj_0->pCumSum = 0.0F;
      for (i = 0; i < 199; i++) {
        obj_0->pCumSumRev[i] = 0.0F;
      }

      obj_0->pCumRevIndex = 1.0F;
      obj_0->pModValueRev = 0.0F;
    }

    /* End of InitializeConditions for MATLABSystem: '<S81>/Moving Average' */
    VCM20_MovingAverage_Init(&VCM20_DW.MovingAverage_p);
    VCM20_MovingAverage_Init(&VCM20_DW.MovingAverage_pn);
    VCM20_MovingAverage_Init(&VCM20_DW.MovingAverage_pna);
    VCM20_MovingAverage_Init(&VCM20_DW.MovingAverage_pnae);
    VCM20_MovingAverage_Init(&VCM20_DW.MovingAverage_pnaev);
  }
}

/* Model terminate function */
void VCM20_terminate(void)
{
  dsp_simulink_MovingAverage_b_T *obj;
  g_dsp_internal_SlidingWindo_b_T *obj_0;

  /* Terminate for S-Function (rti_commonblock): '<S108>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:313 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[2][2] = can_tp1_msg_sleep
            (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X139])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S78>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "Steer SW" Id:512 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X200])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S69>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "Brakes" Id:514 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[2][4] = can_tp1_msg_sleep
            (can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X202])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S21>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "MC_11" Id:643 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][1] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X283])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S27>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "MC_31" Id:647 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][3] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X287])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S24>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "MC_21" Id:644 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][1] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X284])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S30>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "MC_41" Id:648 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][3] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X288])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S68>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "APPS_Steer" Id:513 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[2][4] = can_tp1_msg_sleep
            (can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X201])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for MATLABSystem: '<S81>/Moving Average' */
  obj = &VCM20_DW.obj;
  if (!obj->matlabCodegenIsDeleted) {
    obj->matlabCodegenIsDeleted = true;
    if ((obj->isInitialized == 1) && obj->isSetupComplete) {
      obj_0 = obj->pStatistic;
      if (obj_0->isInitialized == 1) {
        obj_0->isInitialized = 2;
      }

      obj->NumChannels = -1;
    }
  }

  /* End of Terminate for MATLABSystem: '<S81>/Moving Average' */
  /* Terminate for S-Function (rti_commonblock): '<S57>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s004" Id:1284 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][1] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X504])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S58>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s004" Id:1313 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][4] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X521])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S59>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s004" Id:1314 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][4] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X522])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S60>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s004" Id:1793 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][4] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X701])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S49>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1281 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][1] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X501])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S50>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s002" Id:1282 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][1] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X502])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S51>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s005" Id:1285 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][2] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X505])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S28>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "MC_32" Id:649 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][3] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X289])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S22>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "MC_12" Id:645 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][2] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X285])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S31>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "MC_42" Id:650 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][4] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X28A])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S25>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "MC_22" Id:646 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][2] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X286])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S52>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s006" Id:1286 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][2] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X506])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S53>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s007" Id:1287 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][2] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X507])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S54>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s008" Id:1288 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][2] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X508])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S29>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "MC_33" Id:773 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][5] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X305])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S23>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "MC_13" Id:769 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][4] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X301])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S55>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s009" Id:1289 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][2] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X509])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S32>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "MC_43" Id:774 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][5] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X306])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S26>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "MC_23" Id:770 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][4] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X302])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S56>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s00A" Id:1290 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][3] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50A])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S105>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:289 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[2][1] = can_tp1_msg_sleep
            (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X121])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S43>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1291 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][3] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50B])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S106>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:290 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[2][1] = can_tp1_msg_sleep
            (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X122])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S44>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1292 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][3] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50C])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S45>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1293 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][3] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50D])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S111>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:544 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[2][5] = can_tp1_msg_sleep
            (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X220])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S46>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1294 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][3] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50E])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S107>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:308 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[2][2] = can_tp1_msg_sleep
            (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X134])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S47>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1295 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][4] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50F])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S104>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:258 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[2][0] = can_tp1_msg_sleep
            (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X102])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S48>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1296 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][4] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X510])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S41>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s003" Id:1283 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][1] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X503])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S42>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s003" Id:1797 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][5] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X705])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S20>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:1025 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][1] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X401])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S40>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s005" Id:1794 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][5] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X702])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S61>/S-Function1' incorporates:
   *  Constant: '<S38>/Constant'
   */

  /* dSPACE RTICAN TX Message Block: "EVO4s005" Id:1795 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][5] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X703])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S62>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "EVO4s005" Id:1796 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][5] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X704])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S63>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "VCM_11" Id:388 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][0] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X184])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S64>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "VCM_21" Id:389 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][0] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X185])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S65>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "VCM01" Id:392 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][0] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X188])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S66>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "VCM_41" Id:393 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[1][1] = can_tp1_msg_sleep
            (can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X189])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S334>/S-Function1' */

  /* disable digital output channel 4 on port 1 *
   * (set to high-impedance), when the simulation terminates       */
  dio_tp4_digout_init(DIO_TP4_1_MODULE_ADDR, 1, DIO_TP4_MASK_CH4,
                      DIO_TP4_LS_DISABLE, DIO_TP4_HS_DISABLE);

  /* Terminate for S-Function (rti_commonblock): '<S336>/S-Function1' incorporates:
   *  Constant: '<S16>/LED1'
   *  Constant: '<S16>/LED2'
   *  Constant: '<S16>/LED3'
   *  Constant: '<S16>/LED5'
   *  Constant: '<S16>/LED_Displaymode'
   *  Constant: '<S16>/LED_brightness'
   *  Constant: '<S16>/LED_number2'
   *  Constant: '<S16>/buzzer'
   */

  /* dSPACE RTICAN TX Message Block: "Steer LED" Id:528 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X210])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S337>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "Steer LCD" Id:529 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X211])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  VCM20_MovingAverage_Term(&VCM20_DW.MovingAverage_p);
  VCM20_MovingAverage_Term(&VCM20_DW.MovingAverage_pn);
  VCM20_MovingAverage_Term(&VCM20_DW.MovingAverage_pna);
  VCM20_MovingAverage_Term(&VCM20_DW.MovingAverage_pnae);
  VCM20_MovingAverage_Term(&VCM20_DW.MovingAverage_pnaev);

  /* Terminate for S-Function (rti_commonblock): '<S101>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:514 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X202])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S102>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:513 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X201])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S109>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:337 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[2][3] = can_tp1_msg_sleep
            (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X151])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S110>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:369 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[2][3] = can_tp1_msg_sleep
            (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X171])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S103>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:256 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[2][0] = can_tp1_msg_sleep
            (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X100])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }
}

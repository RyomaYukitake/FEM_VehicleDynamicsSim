/*
 * VCM20.h
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

#ifndef RTW_HEADER_VCM20_h_
#define RTW_HEADER_VCM20_h_
#include <string.h>
#include <math.h>
#ifndef VCM20_COMMON_INCLUDES_
#define VCM20_COMMON_INCLUDES_
#include <brtenv.h>
#include <rtkernel.h>
#include <rti_assert.h>
#include <rtidefineddatatypes.h>
#include <rtican_ds1401.h>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* VCM20_COMMON_INCLUDES_ */

#include "VCM20_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rtsplntypes.h"
#include "rt_nonfinite.h"
#include "rtGetInf.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

/* Block signals for system '<S70>/Moving Average' */
typedef struct {
  real_T csumrev[999];
  real_T MovingAverage;                /* '<S70>/Moving Average' */
} B_MovingAverage_VCM20_T;

/* Block states (default storage) for system '<S70>/Moving Average' */
typedef struct {
  dsp_simulink_MovingAverage_VC_T obj; /* '<S70>/Moving Average' */
  boolean_T objisempty;                /* '<S70>/Moving Average' */
} DW_MovingAverage_VCM20_T;

/* Block signals (default storage) */
typedef struct {
  real32_T csumrev[199];
  real_T Saturation;                   /* '<Root>/Saturation' */
  real_T toms;                         /* '<S344>/to m//s' */
  real_T toTire_rpm;                   /* '<S344>/toTire_rpm' */
  real_T Gear;                         /* '<S344>/Gear' */
  real_T Min;                          /* '<S19>/Min' */
  real_T Saturation2;                  /* '<Root>/Saturation2' */
  real_T SFunction1_o1;                /* '<S108>/S-Function1' */
  real_T SFunction1_o2;                /* '<S108>/S-Function1' */
  real_T SFunction1_o3;                /* '<S108>/S-Function1' */
  real_T VELOCITY_X;                   /* '<S14>/Gain8' */
  real_T VELOCITY_X_p;                 /* '<S14>/Gain1' */
  real_T Product;                      /* '<S290>/Product' */
  real_T Sqrt;                         /* '<S290>/Sqrt' */
  real_T Timemsec;                     /* '<S290>/Time msec' */
  real_T Add;                          /* '<S290>/Add' */
  real_T toTire_rpm_b;                 /* '<S290>/toTire_rpm' */
  real_T Gear_g;                       /* '<S290>/Gear' */
  real_T Switch;                       /* '<S309>/Switch' */
  real_T Switch2;                      /* '<S309>/Switch2' */
  real_T Add_m;                        /* '<S310>/Add' */
  real_T Divide;                       /* '<S310>/Divide' */
  real_T SFunction1_o1_o;              /* '<S78>/S-Function1' */
  real_T SFunction1_o2_k;              /* '<S78>/S-Function1' */
  real_T SFunction1_o3_d;              /* '<S78>/S-Function1' */
  real_T SFunction1_o4;                /* '<S78>/S-Function1' */
  real_T SFunction1_o5;                /* '<S78>/S-Function1' */
  real_T SFunction1_o6;                /* '<S78>/S-Function1' */
  real_T SFunction1_o1_l;              /* '<S69>/S-Function1' */
  real_T SFunction1_o2_a;              /* '<S69>/S-Function1' */
  real_T SFunction1_o3_l;              /* '<S69>/S-Function1' */
  real_T SFunction1_o1_m;              /* '<S21>/S-Function1' */
  real_T SFunction1_o2_c;              /* '<S21>/S-Function1' */
  real_T SFunction1_o3_i;              /* '<S21>/S-Function1' */
  real_T SFunction1_o4_d;              /* '<S21>/S-Function1' */
  real_T SFunction1_o5_g;              /* '<S21>/S-Function1' */
  real_T SFunction1_o6_f;              /* '<S21>/S-Function1' */
  real_T SFunction1_o7;                /* '<S21>/S-Function1' */
  real_T SFunction1_o8;                /* '<S21>/S-Function1' */
  real_T SFunction1_o9;                /* '<S21>/S-Function1' */
  real_T SFunction1_o10;               /* '<S21>/S-Function1' */
  real_T SFunction1_o11;               /* '<S21>/S-Function1' */
  real_T SFunction1_o1_b;              /* '<S27>/S-Function1' */
  real_T SFunction1_o2_e;              /* '<S27>/S-Function1' */
  real_T SFunction1_o3_o;              /* '<S27>/S-Function1' */
  real_T SFunction1_o4_g;              /* '<S27>/S-Function1' */
  real_T SFunction1_o5_p;              /* '<S27>/S-Function1' */
  real_T SFunction1_o6_h;              /* '<S27>/S-Function1' */
  real_T SFunction1_o7_n;              /* '<S27>/S-Function1' */
  real_T SFunction1_o8_g;              /* '<S27>/S-Function1' */
  real_T SFunction1_o9_j;              /* '<S27>/S-Function1' */
  real_T SFunction1_o10_n;             /* '<S27>/S-Function1' */
  real_T SFunction1_o11_n;             /* '<S27>/S-Function1' */
  real_T SFunction1_o1_e;              /* '<S24>/S-Function1' */
  real_T SFunction1_o2_b;              /* '<S24>/S-Function1' */
  real_T SFunction1_o3_k;              /* '<S24>/S-Function1' */
  real_T SFunction1_o4_j;              /* '<S24>/S-Function1' */
  real_T SFunction1_o5_e;              /* '<S24>/S-Function1' */
  real_T SFunction1_o6_n;              /* '<S24>/S-Function1' */
  real_T SFunction1_o7_j;              /* '<S24>/S-Function1' */
  real_T SFunction1_o8_o;              /* '<S24>/S-Function1' */
  real_T SFunction1_o9_n;              /* '<S24>/S-Function1' */
  real_T SFunction1_o10_ny;            /* '<S24>/S-Function1' */
  real_T SFunction1_o11_g;             /* '<S24>/S-Function1' */
  real_T SFunction1_o1_g;              /* '<S30>/S-Function1' */
  real_T SFunction1_o2_f;              /* '<S30>/S-Function1' */
  real_T SFunction1_o3_n;              /* '<S30>/S-Function1' */
  real_T SFunction1_o4_n;              /* '<S30>/S-Function1' */
  real_T SFunction1_o5_j;              /* '<S30>/S-Function1' */
  real_T SFunction1_o6_j;              /* '<S30>/S-Function1' */
  real_T SFunction1_o7_i;              /* '<S30>/S-Function1' */
  real_T SFunction1_o8_p;              /* '<S30>/S-Function1' */
  real_T SFunction1_o9_g;              /* '<S30>/S-Function1' */
  real_T SFunction1_o10_l;             /* '<S30>/S-Function1' */
  real_T SFunction1_o11_p;             /* '<S30>/S-Function1' */
  real_T maxrpm;                       /* '<S277>/Min' */
  real_T SFunction1_o1_n;              /* '<S68>/S-Function1' */
  real_T SFunction1_o2_j;              /* '<S68>/S-Function1' */
  real_T SFunction1_o3_d2;             /* '<S68>/S-Function1' */
  real_T SFunction1_o4_h;              /* '<S68>/S-Function1' */
  real_T SFunction1_o5_l;              /* '<S68>/S-Function1' */
  real_T Sum1;                         /* '<S70>/Sum1' */
  real_T Sum3;                         /* '<S70>/Sum3' */
  real_T Divide_j;                     /* '<S70>/Divide' */
  real_T Sum1_p;                       /* '<S71>/Sum1' */
  real_T Sum3_h;                       /* '<S71>/Sum3' */
  real_T Divide_o;                     /* '<S71>/Divide' */
  real_T AppsMni;                      /* '<S67>/AppsMni' */
  real_T Add_j;                        /* '<S67>/Add' */
  real_T Abs;                          /* '<S67>/Abs' */
  real_T Divide_p;                     /* '<S67>/Divide' */
  real_T Saturation_d;                 /* '<S67>/Saturation' */
  real_T APPSsig;                      /* '<S121>/Saturation1' */
  real_T Add_a;                        /* '<S121>/Add' */
  real_T Saturation3;                  /* '<Root>/Saturation3' */
  real_T Switch_c;                     /* '<S167>/Switch' */
  real_T totaltrqreq;                  /* '<S121>/Switch' */
  real_T Gain1;                        /* '<S121>/Gain1' */
  real_T Saturation1;                  /* '<Root>/Saturation1' */
  real_T MinTorque;                    /* '<S121>/Add3' */
  real_T Switch2_e;                    /* '<S121>/Switch2' */
  real_T Gain2;                        /* '<S121>/Gain2' */
  real_T Switch_i;                     /* '<S168>/Switch' */
  real_T Switch2_l;                    /* '<S168>/Switch2' */
  real_T Sum;                          /* '<S127>/Sum' */
  real_T BrakeOverRideSw;              /* '<S115>/Brake Over Ride Sw' */
  real_T Switch1;                      /* '<S126>/Switch1' */
  real_T Switch2_a;                    /* '<S126>/Switch2' */
  real_T Switch3;                      /* '<S126>/Switch3' */
  real_T Switch4;                      /* '<S126>/Switch4' */
  real_T Switch_j;                     /* '<S113>/Switch' */
  real_T Switch_ci;                    /* '<S172>/Switch' */
  real_T TransferFcn;                  /* '<S174>/Transfer Fcn' */
  real_T hgL;                          /* '<S281>/hg//L' */
  real_T Divide_ob;                    /* '<S281>/Divide' */
  real_T Gain;                         /* '<S281>/Gain' */
  real_T m_f;                          /* '<S281>/Gain12' */
  real_T Add_jj;                       /* '<S281>/Add' */
  real_T Saturation1_o;                /* '<S281>/Saturation1' */
  real_T m_r;                          /* '<S281>/Gain13' */
  real_T Add1;                         /* '<S281>/Add1' */
  real_T Saturation2_d;                /* '<S281>/Saturation2' */
  real_T Divide1;                      /* '<S281>/Divide1' */
  real_T Add2;                         /* '<S281>/Add2' */
  real_T Add3;                         /* '<S281>/Add3' */
  real_T Divide2;                      /* '<S281>/Divide2' */
  real_T Gain1_c;                      /* '<S281>/Gain1' */
  real_T Switch_o;                     /* '<S282>/Switch' */
  real_T Switch2_f;                    /* '<S282>/Switch2' */
  real_T Divide3;                      /* '<S281>/Divide3' */
  real_T Divide4;                      /* '<S281>/Divide4' */
  real_T torque_front;                 /* '<S281>/Add5' */
  real_T Add9;                         /* '<S278>/Add9' */
  real_T Saturation3_g;                /* '<S278>/Saturation3' */
  real_T torque_rear;                  /* '<S281>/Add6' */
  real_T Add4;                         /* '<S278>/Add4' */
  real_T Saturation_b;                 /* '<S278>/Saturation' */
  real_T Add8;                         /* '<S278>/Add8' */
  real_T Add10;                        /* '<S278>/Add10' */
  real_T Gain_d;                       /* '<S279>/Gain' */
  real_T Add9_o;                       /* '<S279>/Add9' */
  real_T Saturation3_e;                /* '<S279>/Saturation3' */
  real_T Add7;                         /* '<S278>/Add7' */
  real_T Add11;                        /* '<S278>/Add11' */
  real_T Add4_e;                       /* '<S279>/Add4' */
  real_T Saturation_h;                 /* '<S279>/Saturation' */
  real_T Add7_i;                       /* '<S279>/Add7' */
  real_T torque_rear_d;                /* '<S279>/Add11' */
  real_T Switch1_g;                    /* '<S174>/Switch1' */
  real_T RRtorque;                     /* '<S174>/Gain4' */
  real_T Add1_j;                       /* '<S73>/Add1' */
  real_T TransferFcn_f;                /* '<S73>/Transfer Fcn' */
  real_T Add_i;                        /* '<S73>/Add' */
  real_T Add2_h;                       /* '<S73>/Add2' */
  real_T Switch_m;                     /* '<S73>/Switch' */
  real_T Divide_i;                     /* '<S77>/Divide' */
  real_T Saturation_m;                 /* '<S77>/Saturation' */
  real_T SteerAngle_deg;               /* '<S77>/SteerAngle_deg' */
  real_T SteerAngle_rad;               /* '<S77>/SteerAngle_rad' */
  real_T TargetYawMoment;              /* '<S216>/Gain' */
  real_T Gain12;                       /* '<S211>/Gain12' */
  real_T Gain13;                       /* '<S211>/Gain13' */
  real_T Gain14;                       /* '<S211>/Gain14' */
  real_T Gain4;                        /* '<S211>/Gain4' */
  real_T DYC_deltaTorque;              /* '<S211>/Saturation' */
  real_T Switch1_o;                    /* '<S132>/Switch1' */
  real_T Switch2_j;                    /* '<S132>/Switch2' */
  real_T Switch3_k;                    /* '<S132>/Switch3' */
  real_T Switch4_d;                    /* '<S132>/Switch4' */
  real_T Switch_g;                     /* '<S212>/Switch' */
  real_T Switch_f;                     /* '<S273>/Switch' */
  real_T Switch2_h;                    /* '<S273>/Switch2' */
  real_T RLtorque;                     /* '<S174>/Gain5' */
  real_T Switch_cd;                    /* '<S213>/Switch' */
  real_T Switch_k;                     /* '<S274>/Switch' */
  real_T Switch2_i;                    /* '<S274>/Switch2' */
  real_T Add8_l;                       /* '<S279>/Add8' */
  real_T torque_front_f;               /* '<S279>/Add10' */
  real_T Switch_m0;                    /* '<S174>/Switch' */
  real_T FRtorque;                     /* '<S174>/Gain10' */
  real_T Switch_ci2;                   /* '<S214>/Switch' */
  real_T Switch_h;                     /* '<S275>/Switch' */
  real_T Switch2_d;                    /* '<S275>/Switch2' */
  real_T FLtorque;                     /* '<S174>/Gain11' */
  real_T Switch_e;                     /* '<S215>/Switch' */
  real_T Switch_ed;                    /* '<S276>/Switch' */
  real_T Switch2_c;                    /* '<S276>/Switch2' */
  real_T Switch_p;                     /* '<S290>/Switch' */
  real_T T_loss;                       /* '<S283>/Fd' */
  real_T Add_h;                        /* '<S295>/Add' */
  real_T Divide_h;                     /* '<S295>/Divide' */
  real_T omega;                        /* '<S283>/1// tire R1' */
  real_T wheelInertia;                 /* '<S283>/wheel Inertia' */
  real_T Iomega;                       /* '<S283>/1//Gear' */
  real_T Add_l;                        /* '<S283>/Add' */
  real_T Delay;                        /* '<S296>/Delay' */
  real_T FRLoadDiff;                   /* '<S296>/FR Load Diff' */
  real_T Saturation_g;                 /* '<S296>/Saturation' */
  real_T FrontLoad;                    /* '<S296>/Add1' */
  real_T a1W;                          /* '<S298>/a1*W' */
  real_T Sum_o;                        /* '<S298>/Sum' */
  real_T D;                            /* '<S298>/D' */
  real_T Divide_m;                     /* '<S298>/Divide' */
  real_T atanWa4;                      /* '<S298>/atan(W//a4)' */
  real_T Gain_o;                       /* '<S298>/Gain' */
  real_T sin2atanwa4;                  /* '<S298>/sin(2*atan(w//a4))' */
  real_T BCD;                          /* '<S298>/BCD' */
  real_T Product2;                     /* '<S298>/Product2' */
  real_T Switch_hy;                    /* '<S298>/Switch' */
  real_T B;                            /* '<S298>/B' */
  real_T Product_i;                    /* '<S298>/Product' */
  real_T atanB;                        /* '<S298>/atan(B*É¿)' */
  real_T Sum2;                         /* '<S298>/Sum2' */
  real_T a5W;                          /* '<S298>/a5*W' */
  real_T Sum1_g;                       /* '<S298>/Sum1' */
  real_T EBSAatanBSA;                  /* '<S298>/E*(B*SA-atan(B*SA))' */
  real_T Sum3_i;                       /* '<S298>/Sum3' */
  real_T TrigonometricFunction;        /* '<S298>/Trigonometric Function' */
  real_T Product1;                     /* '<S298>/Product1' */
  real_T TrigonometricFunction2;       /* '<S298>/Trigonometric Function2' */
  real_T Dcosdelta;                    /* '<S298>/+D*cos(delta)' */
  real_T Gain1_a;                      /* '<S298>/Gain1' */
  real_T Rtire;                        /* '<S283>/Rtire' */
  real_T rFFr;                         /* '<S283>/1//Gear1' */
  real_T Add1_jx;                      /* '<S283>/Add1' */
  real_T Add3_h;                       /* '<S294>/Add3' */
  real_T Pgain;                        /* '<S294>/P gain' */
  real_T DiscreteTimeIntegrator;       /* '<S294>/Discrete-Time Integrator' */
  real_T FB;                           /* '<S294>/Add2' */
  real_T torque;                       /* '<S294>/Add' */
  real_T Switch_b;                     /* '<S317>/Switch' */
  real_T Switch2_k;                    /* '<S317>/Switch2' */
  real_T torque_c;                     /* '<S294>/Switch' */
  real_T toTire_rpm_l;                 /* '<S284>/toTire_rpm' */
  real_T Gear_d;                       /* '<S284>/Gear' */
  real_T Add_m0;                       /* '<S302>/Add' */
  real_T Divide_b;                     /* '<S302>/Divide' */
  real_T FrontGain;                    /* '<S284>/FrontGain' */
  real_T Switch_bt;                    /* '<S301>/Switch' */
  real_T Switch2_o;                    /* '<S301>/Switch2' */
  real_T Add1_g;                       /* '<S294>/Add1' */
  real_T Gain1_l;                      /* '<S294>/Gain1' */
  real_T Switch_oi;                    /* '<S318>/Switch' */
  real_T Switch2_e1;                   /* '<S318>/Switch2' */
  real_T Switch1_d;                    /* '<S294>/Switch1' */
  real_T Switch_fd;                    /* '<S289>/Switch' */
  real_T Switch2_ot;                   /* '<S289>/Switch2' */
  real_T RearLoad;                     /* '<S296>/Add2' */
  real_T a1W_p;                        /* '<S297>/a1*W' */
  real_T Sum_k;                        /* '<S297>/Sum' */
  real_T D_h;                          /* '<S297>/D' */
  real_T Divide_e;                     /* '<S297>/Divide' */
  real_T atanWa4_j;                    /* '<S297>/atan(W//a4)' */
  real_T Gain_e;                       /* '<S297>/Gain' */
  real_T sin2atanwa4_k;                /* '<S297>/sin(2*atan(w//a4))' */
  real_T BCD_a;                        /* '<S297>/BCD' */
  real_T Product2_p;                   /* '<S297>/Product2' */
  real_T Switch_bc;                    /* '<S297>/Switch' */
  real_T B_p;                          /* '<S297>/B' */
  real_T Product_e;                    /* '<S297>/Product' */
  real_T atanB_i;                      /* '<S297>/atan(B*É¿)' */
  real_T Sum2_i;                       /* '<S297>/Sum2' */
  real_T a5W_d;                        /* '<S297>/a5*W' */
  real_T Sum1_f;                       /* '<S297>/Sum1' */
  real_T EBSAatanBSA_b;                /* '<S297>/E*(B*SA-atan(B*SA))' */
  real_T Sum3_p;                       /* '<S297>/Sum3' */
  real_T TrigonometricFunction_k;      /* '<S297>/Trigonometric Function' */
  real_T Product1_n;                   /* '<S297>/Product1' */
  real_T TrigonometricFunction2_g;     /* '<S297>/Trigonometric Function2' */
  real_T Dcosdelta_i;                  /* '<S297>/+D*cos(delta)' */
  real_T Gain1_cs;                     /* '<S297>/Gain1' */
  real_T Rtire1;                       /* '<S283>/Rtire1' */
  real_T rFRr;                         /* '<S283>/1//Gear2' */
  real_T Add2_g;                       /* '<S283>/Add2' */
  real_T Add3_k;                       /* '<S291>/Add3' */
  real_T Pgain_i;                      /* '<S291>/P gain' */
  real_T DiscreteTimeIntegrator_h;     /* '<S291>/Discrete-Time Integrator' */
  real_T FB_f;                         /* '<S291>/Add2' */
  real_T torque_b;                     /* '<S291>/Add' */
  real_T Switch_a;                     /* '<S311>/Switch' */
  real_T Switch2_jk;                   /* '<S311>/Switch2' */
  real_T torque_f;                     /* '<S291>/Switch' */
  real_T toTire_rpm_a;                 /* '<S285>/toTire_rpm' */
  real_T Gear_n;                       /* '<S285>/Gear' */
  real_T Add_ib;                       /* '<S304>/Add' */
  real_T Divide_k;                     /* '<S304>/Divide' */
  real_T Switch_pe;                    /* '<S303>/Switch' */
  real_T Switch2_m;                    /* '<S303>/Switch2' */
  real_T Add1_i;                       /* '<S291>/Add1' */
  real_T Gain1_i;                      /* '<S291>/Gain1' */
  real_T Switch_cr;                    /* '<S312>/Switch' */
  real_T Switch2_jy;                   /* '<S312>/Switch2' */
  real_T Switch1_e;                    /* '<S291>/Switch1' */
  real_T Switch_bb;                    /* '<S286>/Switch' */
  real_T Switch2_g;                    /* '<S286>/Switch2' */
  real_T Saturation_a;                 /* '<S186>/Saturation' */
  real_T Saturation1_i;                /* '<S186>/Saturation1' */
  real_T FRMotorPower;                 /* '<S186>/Product' */
  real_T FRMotorPower_W;               /* '<S186>/Gain' */
  real_T FREfficency;                  /* '<S186>/2-D Lookup Table' */
  real_T FRPowerConsumption;           /* '<S186>/Divide' */
  real_T Add3_k0;                      /* '<S292>/Add3' */
  real_T Pgain_a;                      /* '<S292>/P gain' */
  real_T DiscreteTimeIntegrator_i;     /* '<S292>/Discrete-Time Integrator' */
  real_T FB_d;                         /* '<S292>/Add2' */
  real_T torque_g;                     /* '<S292>/Add' */
  real_T Switch_m0p;                   /* '<S313>/Switch' */
  real_T Switch2_az;                   /* '<S313>/Switch2' */
  real_T torque_p;                     /* '<S292>/Switch' */
  real_T Add1_n;                       /* '<S292>/Add1' */
  real_T Gain1_b;                      /* '<S292>/Gain1' */
  real_T Switch_n;                     /* '<S314>/Switch' */
  real_T Switch2_p;                    /* '<S314>/Switch2' */
  real_T Switch1_p;                    /* '<S292>/Switch1' */
  real_T Switch_b4;                    /* '<S287>/Switch' */
  real_T Switch2_lc;                   /* '<S287>/Switch2' */
  real_T Saturation_aa;                /* '<S185>/Saturation' */
  real_T Saturation1_e;                /* '<S185>/Saturation1' */
  real_T FLMotorPower;                 /* '<S185>/Product' */
  real_T FLMotorPower_W;               /* '<S185>/Gain' */
  real_T FLEfficency;                  /* '<S185>/2-D Lookup Table' */
  real_T FLPowerConsumption;           /* '<S185>/Divide' */
  real_T Add3_f;                       /* '<S293>/Add3' */
  real_T Pgain_m;                      /* '<S293>/P gain' */
  real_T DiscreteTimeIntegrator_a;     /* '<S293>/Discrete-Time Integrator' */
  real_T FB_h;                         /* '<S293>/Add2' */
  real_T torque_i;                     /* '<S293>/Add' */
  real_T Switch_eh;                    /* '<S315>/Switch' */
  real_T Switch2_j3;                   /* '<S315>/Switch2' */
  real_T torque_fc;                    /* '<S293>/Switch' */
  real_T Add1_k;                       /* '<S293>/Add1' */
  real_T Gain1_lu;                     /* '<S293>/Gain1' */
  real_T Switch_ah;                    /* '<S316>/Switch' */
  real_T Switch2_h3;                   /* '<S316>/Switch2' */
  real_T Switch1_k;                    /* '<S293>/Switch1' */
  real_T Switch_ma;                    /* '<S288>/Switch' */
  real_T Switch2_hr;                   /* '<S288>/Switch2' */
  real_T Saturation_p;                 /* '<S188>/Saturation' */
  real_T Saturation1_h;                /* '<S188>/Saturation1' */
  real_T RRMotorPower;                 /* '<S188>/Product' */
  real_T RRMotorPower_W;               /* '<S188>/Gain' */
  real_T RREfficency;                  /* '<S188>/2-D Lookup Table' */
  real_T RRPowerConsumption;           /* '<S188>/Divide' */
  real_T Saturation_c;                 /* '<S187>/Saturation' */
  real_T Saturation1_a;                /* '<S187>/Saturation1' */
  real_T RLMotorPower;                 /* '<S187>/Product' */
  real_T RLMotorPower_W;               /* '<S187>/Gain' */
  real_T RLEfficency;                  /* '<S187>/2-D Lookup Table' */
  real_T RLPowerConsumption;           /* '<S187>/Divide' */
  real_T Sum_j;                        /* '<S178>/Sum' */
  real_T Saturation_n;                 /* '<S178>/Saturation' */
  real_T Divide3_p;                    /* '<S178>/Divide3' */
  real_T Eff_INV;                      /* '<S170>/Eff_INV' */
  real_T torque_ba;                    /* '<S178>/Product3' */
  real_T torque_gi;                    /* '<S178>/Saturation4' */
  real_T Divide_iz;                    /* '<S178>/Divide' */
  real_T torque_j;                     /* '<S178>/Product' */
  real_T torque_po;                    /* '<S178>/Saturation1' */
  real_T Saturation_d2;                /* '<S190>/Saturation' */
  real_T Saturation1_hi;               /* '<S190>/Saturation1' */
  real_T FRMotorPower_l;               /* '<S190>/Product' */
  real_T FRMotorPower_W_n;             /* '<S190>/Gain' */
  real_T FREfficency_j;                /* '<S190>/2-D Lookup Table' */
  real_T FRPowerConsumption_p;         /* '<S190>/Divide' */
  real_T Divide1_d;                    /* '<S178>/Divide1' */
  real_T torque_jv;                    /* '<S178>/Product1' */
  real_T torque_m;                     /* '<S178>/Saturation2' */
  real_T Saturation_l;                 /* '<S189>/Saturation' */
  real_T Saturation1_ae;               /* '<S189>/Saturation1' */
  real_T FLMotorPower_g;               /* '<S189>/Product' */
  real_T FLMotorPower_W_b;             /* '<S189>/Gain' */
  real_T FLEfficency_p;                /* '<S189>/2-D Lookup Table' */
  real_T FLPowerConsumption_d;         /* '<S189>/Divide' */
  real_T Divide2_a;                    /* '<S178>/Divide2' */
  real_T torque_iy;                    /* '<S178>/Product2' */
  real_T torque_bu;                    /* '<S178>/Saturation3' */
  real_T Saturation_go;                /* '<S192>/Saturation' */
  real_T Saturation1_ia;               /* '<S192>/Saturation1' */
  real_T RRMotorPower_b;               /* '<S192>/Product' */
  real_T RRMotorPower_W_m;             /* '<S192>/Gain' */
  real_T RREfficency_d;                /* '<S192>/2-D Lookup Table' */
  real_T RRPowerConsumption_k;         /* '<S192>/Divide' */
  real_T Saturation_go2;               /* '<S191>/Saturation' */
  real_T Saturation1_p;                /* '<S191>/Saturation1' */
  real_T RLMotorPower_h;               /* '<S191>/Product' */
  real_T RLMotorPower_W_a;             /* '<S191>/Gain' */
  real_T RLEfficency_j;                /* '<S191>/2-D Lookup Table' */
  real_T RLPowerConsumption_e;         /* '<S191>/Divide' */
  real_T Sum_a;                        /* '<S179>/Sum' */
  real_T Saturation_gs;                /* '<S179>/Saturation' */
  real_T Divide3_g;                    /* '<S179>/Divide3' */
  real_T torque_a;                     /* '<S179>/Product3' */
  real_T torque_o;                     /* '<S179>/Saturation4' */
  real_T Divide_ir;                    /* '<S179>/Divide' */
  real_T torque_az;                    /* '<S179>/Product' */
  real_T torque_l;                     /* '<S179>/Saturation1' */
  real_T Saturation_e;                 /* '<S194>/Saturation' */
  real_T Saturation1_hf;               /* '<S194>/Saturation1' */
  real_T FRMotorPower_n;               /* '<S194>/Product' */
  real_T FRMotorPower_W_o;             /* '<S194>/Gain' */
  real_T FREfficency_p;                /* '<S194>/2-D Lookup Table' */
  real_T FRPowerConsumption_i;         /* '<S194>/Divide' */
  real_T Divide1_i;                    /* '<S179>/Divide1' */
  real_T torque_lg;                    /* '<S179>/Product1' */
  real_T torque_e;                     /* '<S179>/Saturation2' */
  real_T Saturation_hc;                /* '<S193>/Saturation' */
  real_T Saturation1_h5;               /* '<S193>/Saturation1' */
  real_T FLMotorPower_o;               /* '<S193>/Product' */
  real_T FLMotorPower_W_c;             /* '<S193>/Gain' */
  real_T FLEfficency_o;                /* '<S193>/2-D Lookup Table' */
  real_T FLPowerConsumption_i;         /* '<S193>/Divide' */
  real_T Divide2_e;                    /* '<S179>/Divide2' */
  real_T torque_f4;                    /* '<S179>/Product2' */
  real_T torque_bc;                    /* '<S179>/Saturation3' */
  real_T Saturation_gy;                /* '<S196>/Saturation' */
  real_T Saturation1_c;                /* '<S196>/Saturation1' */
  real_T RRMotorPower_j;               /* '<S196>/Product' */
  real_T RRMotorPower_W_h;             /* '<S196>/Gain' */
  real_T RREfficency_n;                /* '<S196>/2-D Lookup Table' */
  real_T RRPowerConsumption_b;         /* '<S196>/Divide' */
  real_T Saturation_gc;                /* '<S195>/Saturation' */
  real_T Saturation1_f;                /* '<S195>/Saturation1' */
  real_T RLMotorPower_m;               /* '<S195>/Product' */
  real_T RLMotorPower_W_i;             /* '<S195>/Gain' */
  real_T RLEfficency_f;                /* '<S195>/2-D Lookup Table' */
  real_T RLPowerConsumption_i;         /* '<S195>/Divide' */
  real_T Sum_a5;                       /* '<S180>/Sum' */
  real_T Saturation_gm;                /* '<S180>/Saturation' */
  real_T Divide3_l;                    /* '<S180>/Divide3' */
  real_T torque_i5;                    /* '<S180>/Product3' */
  real_T torque_or;                    /* '<S180>/Saturation4' */
  real_T Divide_mt;                    /* '<S180>/Divide' */
  real_T torque_d;                     /* '<S180>/Product' */
  real_T torque_d0;                    /* '<S180>/Saturation1' */
  real_T Saturation_e2;                /* '<S198>/Saturation' */
  real_T Saturation1_j;                /* '<S198>/Saturation1' */
  real_T FRMotorPower_h;               /* '<S198>/Product' */
  real_T FRMotorPower_W_h;             /* '<S198>/Gain' */
  real_T FREfficency_n;                /* '<S198>/2-D Lookup Table' */
  real_T FRPowerConsumption_m;         /* '<S198>/Divide' */
  real_T Divide1_k;                    /* '<S180>/Divide1' */
  real_T torque_ae;                    /* '<S180>/Product1' */
  real_T torque_mg;                    /* '<S180>/Saturation2' */
  real_T Saturation_o;                 /* '<S197>/Saturation' */
  real_T Saturation1_f1;               /* '<S197>/Saturation1' */
  real_T FLMotorPower_k;               /* '<S197>/Product' */
  real_T FLMotorPower_W_bv;            /* '<S197>/Gain' */
  real_T FLEfficency_n;                /* '<S197>/2-D Lookup Table' */
  real_T FLPowerConsumption_k;         /* '<S197>/Divide' */
  real_T Divide2_ey;                   /* '<S180>/Divide2' */
  real_T torque_ms;                    /* '<S180>/Product2' */
  real_T torque_k;                     /* '<S180>/Saturation3' */
  real_T Saturation_co;                /* '<S200>/Saturation' */
  real_T Saturation1_ey;               /* '<S200>/Saturation1' */
  real_T RRMotorPower_p;               /* '<S200>/Product' */
  real_T RRMotorPower_W_i;             /* '<S200>/Gain' */
  real_T RREfficency_ne;               /* '<S200>/2-D Lookup Table' */
  real_T RRPowerConsumption_bl;        /* '<S200>/Divide' */
  real_T Saturation_c2;                /* '<S199>/Saturation' */
  real_T Saturation1_b;                /* '<S199>/Saturation1' */
  real_T RLMotorPower_j;               /* '<S199>/Product' */
  real_T RLMotorPower_W_m;             /* '<S199>/Gain' */
  real_T RLEfficency_b;                /* '<S199>/2-D Lookup Table' */
  real_T RLPowerConsumption_d;         /* '<S199>/Divide' */
  real_T Sum_e;                        /* '<S181>/Sum' */
  real_T Saturation_ak;                /* '<S181>/Saturation' */
  real_T Divide3_n;                    /* '<S181>/Divide3' */
  real_T torque_m1;                    /* '<S181>/Product3' */
  real_T torque_a4;                    /* '<S181>/Saturation4' */
  real_T Divide_d;                     /* '<S181>/Divide' */
  real_T torque_m3;                    /* '<S181>/Product' */
  real_T torque_cl;                    /* '<S181>/Saturation1' */
  real_T Saturation_p5;                /* '<S202>/Saturation' */
  real_T Saturation1_i3;               /* '<S202>/Saturation1' */
  real_T FRMotorPower_h3;              /* '<S202>/Product' */
  real_T FRMotorPower_W_f;             /* '<S202>/Gain' */
  real_T FREfficency_f;                /* '<S202>/2-D Lookup Table' */
  real_T FRPowerConsumption_c;         /* '<S202>/Divide' */
  real_T Divide1_iv;                   /* '<S181>/Divide1' */
  real_T torque_ju;                    /* '<S181>/Product1' */
  real_T torque_fy;                    /* '<S181>/Saturation2' */
  real_T Saturation_i;                 /* '<S201>/Saturation' */
  real_T Saturation1_d;                /* '<S201>/Saturation1' */
  real_T FLMotorPower_d;               /* '<S201>/Product' */
  real_T FLMotorPower_W_o;             /* '<S201>/Gain' */
  real_T FLEfficency_i;                /* '<S201>/2-D Lookup Table' */
  real_T FLPowerConsumption_l;         /* '<S201>/Divide' */
  real_T Divide2_j;                    /* '<S181>/Divide2' */
  real_T torque_n;                     /* '<S181>/Product2' */
  real_T torque_ag;                    /* '<S181>/Saturation3' */
  real_T Saturation_ha;                /* '<S204>/Saturation' */
  real_T Saturation1_fg;               /* '<S204>/Saturation1' */
  real_T RRMotorPower_m;               /* '<S204>/Product' */
  real_T RRMotorPower_W_e;             /* '<S204>/Gain' */
  real_T RREfficency_d1;               /* '<S204>/2-D Lookup Table' */
  real_T RRPowerConsumption_kk;        /* '<S204>/Divide' */
  real_T Saturation_f;                 /* '<S203>/Saturation' */
  real_T Saturation1_fj;               /* '<S203>/Saturation1' */
  real_T RLMotorPower_p;               /* '<S203>/Product' */
  real_T RLMotorPower_W_n;             /* '<S203>/Gain' */
  real_T RLEfficency_a;                /* '<S203>/2-D Lookup Table' */
  real_T RLPowerConsumption_b;         /* '<S203>/Divide' */
  real_T Sum_b;                        /* '<S182>/Sum' */
  real_T Saturation_l5;                /* '<S182>/Saturation' */
  real_T Divide3_d;                    /* '<S182>/Divide3' */
  real_T torque_dc;                    /* '<S182>/Product3' */
  real_T torque_ej;                    /* '<S182>/Saturation4' */
  real_T SFunction1;                   /* '<S112>/S-Function1' */
  real_T Current;                      /* '<S14>/Fcn1' */
  real_T Product_a;                    /* '<S173>/Product' */
  real_T Sum_aa;                       /* '<S183>/Sum' */
  real_T Gain_p;                       /* '<S183>/Gain' */
  real_T Divide_hj;                    /* '<S183>/Divide' */
  real_T Sum4;                         /* '<S183>/Sum4' */
  real_T torque_lx;                    /* '<S183>/Switch2' */
  real_T FLtorque_c;                   /* '<S184>/Switch3' */
  real_T FLtrq;                        /* '<S122>/Gain3' */
  real_T Gain30;                       /* '<S37>/Gain30' */
  real_T Divide2_k;                    /* '<S182>/Divide2' */
  real_T torque_gr;                    /* '<S182>/Product2' */
  real_T torque_fj;                    /* '<S182>/Saturation3' */
  real_T Sum3_n;                       /* '<S183>/Sum3' */
  real_T torque_ok;                    /* '<S183>/Switch3' */
  real_T FRtorque_a;                   /* '<S184>/Switch2' */
  real_T FRtrq;                        /* '<S122>/Gain2' */
  real_T Gain31;                       /* '<S37>/Gain31' */
  real_T Divide_d5;                    /* '<S182>/Divide' */
  real_T torque_dz;                    /* '<S182>/Product' */
  real_T torque_h;                     /* '<S182>/Saturation1' */
  real_T Sum1_l;                       /* '<S183>/Sum1' */
  real_T torque_co;                    /* '<S183>/Switch' */
  real_T RRtorque_e;                   /* '<S184>/Switch' */
  real_T RRtrq;                        /* '<S122>/Gain' */
  real_T Gain32;                       /* '<S37>/Gain32' */
  real_T Divide1_h;                    /* '<S182>/Divide1' */
  real_T torque_ka;                    /* '<S182>/Product1' */
  real_T torque_ku;                    /* '<S182>/Saturation2' */
  real_T Sum2_g;                       /* '<S183>/Sum2' */
  real_T torque_nn;                    /* '<S183>/Switch1' */
  real_T RLtorque_d;                   /* '<S184>/Switch1' */
  real_T RLtrq;                        /* '<S122>/Gain1' */
  real_T Gain33;                       /* '<S37>/Gain33' */
  real_T Gain1_e;                      /* '<S37>/Gain1' */
  real_T vms;                          /* '<S220>/Transfer Fcn' */
  real_T Product5;                     /* '<S222>/Product5' */
  real_T Divide3_j;                    /* '<S222>/Divide3' */
  real_T Add4_j;                       /* '<S222>/Add4' */
  real_T Divide2_o;                    /* '<S222>/Divide2' */
  real_T Gain6;                        /* '<S222>/Gain6' */
  real_T Divide_a;                     /* '<S222>/Divide' */
  real_T Gain2_b;                      /* '<S37>/Gain2' */
  real_T powerconsumptionkW;           /* '<S176>/Gain' */
  real_T Gain3;                        /* '<S37>/Gain3' */
  real_T toTire_rpm_c;                 /* '<S322>/toTire_rpm' */
  real_T Gear_nh;                      /* '<S322>/Gear' */
  real_T Add_i3;                       /* '<S322>/Add' */
  real_T Switch_i0;                    /* '<S326>/Switch' */
  real_T FL_Sliprate;                  /* '<S322>/Divide' */
  real_T Gain4_j;                      /* '<S37>/Gain4' */
  real_T toTire_rpm_lr;                /* '<S323>/toTire_rpm' */
  real_T Gear_dd;                      /* '<S323>/Gear' */
  real_T Add_k;                        /* '<S323>/Add' */
  real_T Switch_ahz;                   /* '<S328>/Switch' */
  real_T FR_Sliprate;                  /* '<S323>/Divide' */
  real_T Gain5;                        /* '<S37>/Gain5' */
  real_T toTire_rpm_f;                 /* '<S324>/toTire_rpm' */
  real_T Gear_m;                       /* '<S324>/Gear' */
  real_T Add_n;                        /* '<S324>/Add' */
  real_T Switch_il;                    /* '<S330>/Switch' */
  real_T RL_Sliprate;                  /* '<S324>/Divide' */
  real_T Gain6_l;                      /* '<S37>/Gain6' */
  real_T toTire_rpm_a5;                /* '<S325>/toTire_rpm' */
  real_T Gear_dp;                      /* '<S325>/Gear' */
  real_T Add_ji;                       /* '<S325>/Add' */
  real_T Switch_np;                    /* '<S332>/Switch' */
  real_T RR_Sliprate;                  /* '<S325>/Divide' */
  real_T Gain7;                        /* '<S37>/Gain7' */
  real_T Cast;                         /* '<S37>/Cast' */
  real_T Gain18;                       /* '<S37>/Gain18' */
  real_T Cast1;                        /* '<S37>/Cast1' */
  real_T Gain20;                       /* '<S37>/Gain20' */
  real_T Cast2;                        /* '<S37>/Cast2' */
  real_T Gain21;                       /* '<S37>/Gain21' */
  real_T Cast3;                        /* '<S37>/Cast3' */
  real_T Gain22;                       /* '<S37>/Gain22' */
  real_T Cast4;                        /* '<S37>/Cast4' */
  real_T Gain23;                       /* '<S37>/Gain23' */
  real_T Cast5;                        /* '<S37>/Cast5' */
  real_T Gain24;                       /* '<S37>/Gain24' */
  real_T Cast6;                        /* '<S37>/Cast6' */
  real_T Gain25;                       /* '<S37>/Gain25' */
  real_T Gain19;                       /* '<S36>/Gain19' */
  real_T Gain26;                       /* '<S36>/Gain26' */
  real_T Gain27;                       /* '<S36>/Gain27' */
  real_T Gain28;                       /* '<S36>/Gain28' */
  real_T SFunction1_o1_h;              /* '<S28>/S-Function1' */
  real_T SFunction1_o2_d;              /* '<S28>/S-Function1' */
  real_T SFunction1_o3_o3;             /* '<S28>/S-Function1' */
  real_T SFunction1_o4_dd;             /* '<S28>/S-Function1' */
  real_T SFunction1_o1_ml;             /* '<S22>/S-Function1' */
  real_T SFunction1_o2_o;              /* '<S22>/S-Function1' */
  real_T SFunction1_o3_e;              /* '<S22>/S-Function1' */
  real_T SFunction1_o4_c;              /* '<S22>/S-Function1' */
  real_T SFunction1_o1_i;              /* '<S31>/S-Function1' */
  real_T SFunction1_o2_n;              /* '<S31>/S-Function1' */
  real_T SFunction1_o3_f;              /* '<S31>/S-Function1' */
  real_T SFunction1_o4_e;              /* '<S31>/S-Function1' */
  real_T SFunction1_o1_eq;             /* '<S25>/S-Function1' */
  real_T SFunction1_o2_oc;             /* '<S25>/S-Function1' */
  real_T SFunction1_o3_d2j;            /* '<S25>/S-Function1' */
  real_T SFunction1_o4_f;              /* '<S25>/S-Function1' */
  real_T SFunction1_o1_j;              /* '<S29>/S-Function1' */
  real_T SFunction1_o2_p;              /* '<S29>/S-Function1' */
  real_T SFunction1_o3_h;              /* '<S29>/S-Function1' */
  real_T SFunction1_o4_hn;             /* '<S29>/S-Function1' */
  real_T SFunction1_o1_c;              /* '<S23>/S-Function1' */
  real_T SFunction1_o2_an;             /* '<S23>/S-Function1' */
  real_T SFunction1_o3_g;              /* '<S23>/S-Function1' */
  real_T SFunction1_o4_jb;             /* '<S23>/S-Function1' */
  real_T Min_a;                        /* '<S36>/Min' */
  real_T SFunction1_o1_f;              /* '<S32>/S-Function1' */
  real_T SFunction1_o2_o4;             /* '<S32>/S-Function1' */
  real_T SFunction1_o3_io;             /* '<S32>/S-Function1' */
  real_T SFunction1_o4_ea;             /* '<S32>/S-Function1' */
  real_T SFunction1_o1_om;             /* '<S26>/S-Function1' */
  real_T SFunction1_o2_ao;             /* '<S26>/S-Function1' */
  real_T SFunction1_o3_iq;             /* '<S26>/S-Function1' */
  real_T SFunction1_o4_p;              /* '<S26>/S-Function1' */
  real_T Min1;                         /* '<S36>/Min1' */
  real_T SFunction1_o1_m3;             /* '<S105>/S-Function1' */
  real_T SFunction1_o2_i;              /* '<S105>/S-Function1' */
  real_T SFunction1_o3_nw;             /* '<S105>/S-Function1' */
  real_T ACCEL_X;                      /* '<S14>/ACCEL_X_Gain' */
  real_T ACCEL_X_m;                    /* '<S14>/Gain' */
  real_T ACCEL_X_mq;                   /* '<S35>/Gain3' */
  real_T ACCEL_Y;                      /* '<S14>/ACCEL_Y_Gain' */
  real_T ACCEL_Y_b;                    /* '<S35>/Gain5' */
  real_T ACCEL_Z;                      /* '<S14>/ACCEL_Z_Gain' */
  real_T ACCEL_Y_e;                    /* '<S35>/Gain4' */
  real_T SFunction1_o1_bx;             /* '<S106>/S-Function1' */
  real_T SFunction1_o2_bg;             /* '<S106>/S-Function1' */
  real_T SFunction1_o3_b;              /* '<S106>/S-Function1' */
  real_T GYRO_X;                       /* '<S14>/Gain3' */
  real_T GYRO_X_o;                     /* '<S35>/Gain6' */
  real_T GYRO_Y;                       /* '<S14>/Gain5' */
  real_T GYRO_Y_o;                     /* '<S35>/Gain8' */
  real_T GYRO_Z;                       /* '<S14>/Gain4' */
  real_T GYRO_Z_m;                     /* '<S35>/Gain7' */
  real_T VELOCITY_X_p5;                /* '<S35>/Gain9' */
  real_T VELOCITY_Y;                   /* '<S14>/Gain10' */
  real_T VELOCITY_Y_j;                 /* '<S35>/Gain10' */
  real_T VELOCITY_Z;                   /* '<S14>/Gain9' */
  real_T VELOCITY_Z_a;                 /* '<S35>/Gain14' */
  real_T SFunction1_o1_hc;             /* '<S111>/S-Function1' */
  real_T SFunction1_o2_nj;             /* '<S111>/S-Function1' */
  real_T SFunction1_o3_dd;             /* '<S111>/S-Function1' */
  real_T ANGLE_TRACK;                  /* '<S14>/Gain24' */
  real_T ANGLE_TRACK_g;                /* '<S35>/Gain11' */
  real_T ANGLE_SLIP;                   /* '<S14>/Gain23' */
  real_T ANGLE_SLIP_d;                 /* '<S35>/Gain13' */
  real_T CURVATURE_RADIUS;             /* '<S14>/Gain25' */
  real_T CURVATURE_RADIUS_a;           /* '<S35>/Gain12' */
  real_T SFunction1_o1_nu;             /* '<S107>/S-Function1' */
  real_T SFunction1_o2_ej;             /* '<S107>/S-Function1' */
  real_T LATITUDE;                     /* '<S14>/Gain21' */
  real_T LATITUDE_c;                   /* '<S35>/Gain16' */
  real_T LONGITUDE;                    /* '<S14>/Gain22' */
  real_T LONGITUDE_a;                  /* '<S35>/Gain15' */
  real_T SFunction1_o1_or;             /* '<S104>/S-Function1' */
  real_T SFunction1_o2_l;              /* '<S104>/S-Function1' */
  real_T SFunction1_o3_nj;             /* '<S104>/S-Function1' */
  real_T SOLUTION_MODE;                /* '<S14>/Gain18' */
  real_T SOLUTION_MODE_g;              /* '<S35>/Gain2' */
  real_T VELOCITY_VALID;               /* '<S14>/Gain19' */
  real_T VELOCITY_VALID_b;             /* '<S35>/Gain1' */
  real_T POSITION_VALID;               /* '<S14>/Gain20' */
  real_T POSITION_VALID_a;             /* '<S35>/Gain18' */
  real_T Gain_n;                       /* '<S34>/Gain' */
  real_T Gain29;                       /* '<S34>/Gain29' */
  real_T Gain2_l;                      /* '<S34>/Gain2' */
  real_T Gain1_j;                      /* '<S34>/Gain1' */
  real_T DataTypeConversion8;          /* '<S34>/Data Type Conversion8' */
  real_T SFunction1_o1_lm;             /* '<S20>/S-Function1' */
  real_T SFunction1_o2_lh;             /* '<S20>/S-Function1' */
  real_T SFunction1_o3_ol;             /* '<S20>/S-Function1' */
  real_T SFunction1_o4_i;              /* '<S20>/S-Function1' */
  real_T SFunction1_o5_h;              /* '<S20>/S-Function1' */
  real_T SFunction1_o6_d;              /* '<S20>/S-Function1' */
  real_T SFunction1_o7_k;              /* '<S20>/S-Function1' */
  real_T SFunction1_o8_go;             /* '<S20>/S-Function1' */
  real_T SFunction1_m;                 /* '<S342>/S-Function1' */
  real_T Product_g;                    /* '<S339>/Product' */
  real_T reikyaku_temp2;               /* '<S339>/Add' */
  real_T SFunction1_c;                 /* '<S343>/S-Function1' */
  real_T Product_b;                    /* '<S340>/Product' */
  real_T reikyaku_temp3;               /* '<S340>/Add' */
  real_T DataTypeConversion10;         /* '<S39>/Data Type Conversion10' */
  real_T DataTypeConversion1;          /* '<S39>/Data Type Conversion1' */
  real_T DataTypeConversion2;          /* '<S39>/Data Type Conversion2' */
  real_T Switch1_l;                    /* '<S133>/Switch1' */
  real_T Switch2_eg;                   /* '<S133>/Switch2' */
  real_T Switch3_a;                    /* '<S133>/Switch3' */
  real_T Switch4_i;                    /* '<S133>/Switch4' */
  real_T Switch1_a;                    /* '<S131>/Switch1' */
  real_T Switch2_fm;                   /* '<S131>/Switch2' */
  real_T Switch3_d;                    /* '<S131>/Switch3' */
  real_T Switch4_g;                    /* '<S131>/Switch4' */
  real_T Switch1_h;                    /* '<S124>/Switch1' */
  real_T DataTypeConversion;           /* '<S5>/Data Type Conversion' */
  real_T Switch1_dh;                   /* '<S134>/Switch1' */
  real_T Switch2_fz;                   /* '<S134>/Switch2' */
  real_T Switch3_ad;                   /* '<S134>/Switch3' */
  real_T Switch4_ie;                   /* '<S134>/Switch4' */
  real_T MinDCbusVoltage;              /* '<S116>/Min DC bus Voltage' */
  real_T DataTypeConversion2_g;        /* '<S5>/Data Type Conversion2' */
  real_T Switch1_m;                    /* '<S123>/Switch1' */
  real_T DataTypeConversion1_i;        /* '<S5>/Data Type Conversion1' */
  real_T AMK_bErrorReset;              /* '<S114>/Data Type Conversion' */
  real_T FRspeed;                      /* '<S122>/Gain8' */
  real_T Min_h;                        /* '<S125>/Min' */
  real_T VectorConcatenate[2];         /* '<S136>/Vector Concatenate' */
  real_T Fcn;                          /* '<S136>/Fcn' */
  real_T MC1;                          /* '<S125>/MC1' */
  real_T Switch_nf;                    /* '<S142>/Switch' */
  real_T Switch2_n;                    /* '<S142>/Switch2' */
  real_T Gain2_n;                      /* '<S117>/Gain2' */
  real_T Quantizer1;                   /* '<S117>/Quantizer1' */
  real_T Switch_as;                    /* '<S138>/Switch' */
  real_T Delay_g;                      /* '<S138>/Delay' */
  real_T MinMax;                       /* '<S138>/MinMax' */
  real_T Delay_o;                      /* '<S140>/Delay' */
  real_T MinMax_c;                     /* '<S140>/MinMax' */
  real_T Positivetorquelimit;          /* '<S137>/MinMax' */
  real_T VectorConcatenate_a[2];       /* '<S135>/Vector Concatenate' */
  real_T Fcn_p;                        /* '<S135>/Fcn' */
  real_T Switch_e3;                    /* '<S141>/Switch' */
  real_T Switch2_di;                   /* '<S141>/Switch2' */
  real_T Gain6_lv;                     /* '<S117>/Gain6' */
  real_T Quantizer2;                   /* '<S117>/Quantizer2' */
  real_T Delay_l;                      /* '<S139>/Delay' */
  real_T MinMax_co;                    /* '<S139>/MinMax' */
  real_T Switch2_aw;                   /* '<S124>/Switch2' */
  real_T DataTypeConversion3;          /* '<S5>/Data Type Conversion3' */
  real_T DataTypeConversion5;          /* '<S5>/Data Type Conversion5' */
  real_T Switch2_g1;                   /* '<S123>/Switch2' */
  real_T DataTypeConversion4;          /* '<S5>/Data Type Conversion4' */
  real_T RRspeed;                      /* '<S122>/Gain6' */
  real_T Min1_i;                       /* '<S125>/Min1' */
  real_T VectorConcatenate_b[2];       /* '<S144>/Vector Concatenate' */
  real_T Fcn_l;                        /* '<S144>/Fcn' */
  real_T MC2;                          /* '<S125>/MC2' */
  real_T Switch_l;                     /* '<S150>/Switch' */
  real_T Switch2_pb;                   /* '<S150>/Switch2' */
  real_T Gain2_bd;                     /* '<S118>/Gain2' */
  real_T Quantizer1_p;                 /* '<S118>/Quantizer1' */
  real_T Switch_kd;                    /* '<S146>/Switch' */
  real_T Delay_ox;                     /* '<S146>/Delay' */
  real_T MinMax_g;                     /* '<S146>/MinMax' */
  real_T Delay_gx;                     /* '<S148>/Delay' */
  real_T MinMax_i;                     /* '<S148>/MinMax' */
  real_T Positivetorquelimit_j;        /* '<S145>/MinMax' */
  real_T VectorConcatenate_p[2];       /* '<S143>/Vector Concatenate' */
  real_T Fcn_k;                        /* '<S143>/Fcn' */
  real_T Switch_bm;                    /* '<S149>/Switch' */
  real_T Switch2_pv;                   /* '<S149>/Switch2' */
  real_T Gain6_e;                      /* '<S118>/Gain6' */
  real_T Quantizer2_p;                 /* '<S118>/Quantizer2' */
  real_T Delay_j;                      /* '<S147>/Delay' */
  real_T MinMax_m;                     /* '<S147>/MinMax' */
  real_T Switch3_j;                    /* '<S124>/Switch3' */
  real_T DataTypeConversion6;          /* '<S5>/Data Type Conversion6' */
  real_T DataTypeConversion7;          /* '<S5>/Data Type Conversion7' */
  real_T Switch3_n;                    /* '<S123>/Switch3' */
  real_T DataTypeConversion8_l;        /* '<S5>/Data Type Conversion8' */
  real_T RLspeed;                      /* '<S122>/Gain7' */
  real_T Min2;                         /* '<S125>/Min2' */
  real_T VectorConcatenate_j[2];       /* '<S152>/Vector Concatenate' */
  real_T Fcn_b;                        /* '<S152>/Fcn' */
  real_T MC3;                          /* '<S125>/MC3' */
  real_T Switch_j2;                    /* '<S158>/Switch' */
  real_T Switch2_oq;                   /* '<S158>/Switch2' */
  real_T Gain2_h;                      /* '<S119>/Gain2' */
  real_T Quantizer1_m;                 /* '<S119>/Quantizer1' */
  real_T Switch_k1;                    /* '<S154>/Switch' */
  real_T Delay_m;                      /* '<S154>/Delay' */
  real_T MinMax_h;                     /* '<S154>/MinMax' */
  real_T Delay_i;                      /* '<S156>/Delay' */
  real_T MinMax_k;                     /* '<S156>/MinMax' */
  real_T Positivetorquelimit_m;        /* '<S153>/MinMax' */
  real_T VectorConcatenate_ju[2];      /* '<S151>/Vector Concatenate' */
  real_T Fcn_c;                        /* '<S151>/Fcn' */
  real_T Switch_c1;                    /* '<S157>/Switch' */
  real_T Switch2_cp;                   /* '<S157>/Switch2' */
  real_T Gain6_p;                      /* '<S119>/Gain6' */
  real_T Quantizer2_n;                 /* '<S119>/Quantizer2' */
  real_T Delay_m5;                     /* '<S155>/Delay' */
  real_T MinMax_d;                     /* '<S155>/MinMax' */
  real_T Switch4_b;                    /* '<S124>/Switch4' */
  real_T DataTypeConversion9;          /* '<S5>/Data Type Conversion9' */
  real_T DataTypeConversion11;         /* '<S5>/Data Type Conversion11' */
  real_T Switch4_f;                    /* '<S123>/Switch4' */
  real_T DataTypeConversion10_l;       /* '<S5>/Data Type Conversion10' */
  real_T FLspeed;                      /* '<S122>/Gain9' */
  real_T Min3;                         /* '<S125>/Min3' */
  real_T VectorConcatenate_m[2];       /* '<S160>/Vector Concatenate' */
  real_T Fcn_o;                        /* '<S160>/Fcn' */
  real_T MC4;                          /* '<S125>/MC4' */
  real_T Switch_ge;                    /* '<S166>/Switch' */
  real_T Switch2_ff;                   /* '<S166>/Switch2' */
  real_T Gain2_b1;                     /* '<S120>/Gain2' */
  real_T Quantizer1_n;                 /* '<S120>/Quantizer1' */
  real_T Switch_f4;                    /* '<S162>/Switch' */
  real_T Delay_gq;                     /* '<S162>/Delay' */
  real_T MinMax_l;                     /* '<S162>/MinMax' */
  real_T Delay_n;                      /* '<S164>/Delay' */
  real_T MinMax_lq;                    /* '<S164>/MinMax' */
  real_T Positivetorquelimit_p;        /* '<S161>/MinMax' */
  real_T VectorConcatenate_g[2];       /* '<S159>/Vector Concatenate' */
  real_T Fcn_oe;                       /* '<S159>/Fcn' */
  real_T Switch_oi0;                   /* '<S165>/Switch' */
  real_T Switch2_lu;                   /* '<S165>/Switch2' */
  real_T Gain6_f;                      /* '<S120>/Gain6' */
  real_T Quantizer2_e;                 /* '<S120>/Quantizer2' */
  real_T Delay_a;                      /* '<S163>/Delay' */
  real_T MinMax_iq;                    /* '<S163>/MinMax' */
  real_T DataTypeConversion14;         /* '<S16>/Data Type Conversion14' */
  real_T DataTypeConversion13;         /* '<S16>/Data Type Conversion13' */
  real_T DataTypeConversion12;         /* '<S16>/Data Type Conversion12' */
  real_T DataTypeConversion11_f;       /* '<S16>/Data Type Conversion11' */
  real_T DataTypeConversion8_f;        /* '<S16>/Data Type Conversion8' */
  real_T DataTypeConversion9_o;        /* '<S16>/Data Type Conversion9' */
  real_T DataTypeConversion10_a;       /* '<S16>/Data Type Conversion10' */
  real_T DataTypeConversion_l;         /* '<S16>/Data Type Conversion' */
  real_T DataTypeConversion1_o;        /* '<S16>/Data Type Conversion1' */
  real_T DataTypeConversion2_h;        /* '<S16>/Data Type Conversion2' */
  real_T DataTypeConversion3_i;        /* '<S16>/Data Type Conversion3' */
  real_T DataTypeConversion4_b;        /* '<S16>/Data Type Conversion4' */
  real_T DataTypeConversion5_h;        /* '<S16>/Data Type Conversion5' */
  real_T DataTypeConversion6_d;        /* '<S16>/Data Type Conversion6' */
  real_T DataTypeConversion7_j;        /* '<S16>/Data Type Conversion7' */
  real_T Add1_jn;                      /* '<S72>/Add1' */
  real_T TransferFcn_a;                /* '<S72>/Transfer Fcn' */
  real_T Add_b;                        /* '<S72>/Add' */
  real_T Add2_o;                       /* '<S72>/Add2' */
  real_T Switch_co;                    /* '<S72>/Switch' */
  real_T APPS1_noncalibrated;          /* '<S11>/Gain' */
  real_T APPS2_noncalibrated;          /* '<S11>/Gain1' */
  real_T Delay10;                      /* '<S89>/Delay10' */
  real_T Delay11;                      /* '<S89>/Delay11' */
  real_T Delay12;                      /* '<S89>/Delay12' */
  real_T Delay13;                      /* '<S89>/Delay13' */
  real_T Delay14;                      /* '<S89>/Delay14' */
  real_T Delay15;                      /* '<S89>/Delay15' */
  real_T Delay16;                      /* '<S89>/Delay16' */
  real_T Delay17;                      /* '<S89>/Delay17' */
  real_T Delay18;                      /* '<S89>/Delay18' */
  real_T Delay19;                      /* '<S89>/Delay19' */
  real_T Min1_g;                       /* '<S89>/Min1' */
  real_T Delay_ih;                     /* '<S89>/Delay' */
  real_T Delay1;                       /* '<S89>/Delay1' */
  real_T Delay2;                       /* '<S89>/Delay2' */
  real_T Delay3;                       /* '<S89>/Delay3' */
  real_T Delay4;                       /* '<S89>/Delay4' */
  real_T Delay5;                       /* '<S89>/Delay5' */
  real_T Delay6;                       /* '<S89>/Delay6' */
  real_T Delay7;                       /* '<S89>/Delay7' */
  real_T Delay8;                       /* '<S89>/Delay8' */
  real_T Delay9;                       /* '<S89>/Delay9' */
  real_T Min_k;                        /* '<S89>/Min' */
  real_T Delay10_h;                    /* '<S90>/Delay10' */
  real_T Delay11_f;                    /* '<S90>/Delay11' */
  real_T Delay12_l;                    /* '<S90>/Delay12' */
  real_T Delay13_h;                    /* '<S90>/Delay13' */
  real_T Delay14_n;                    /* '<S90>/Delay14' */
  real_T Delay15_p;                    /* '<S90>/Delay15' */
  real_T Delay16_d;                    /* '<S90>/Delay16' */
  real_T Delay17_i;                    /* '<S90>/Delay17' */
  real_T Delay18_n;                    /* '<S90>/Delay18' */
  real_T Delay19_j;                    /* '<S90>/Delay19' */
  real_T Min1_e;                       /* '<S90>/Min1' */
  real_T Delay_nn;                     /* '<S90>/Delay' */
  real_T Delay1_i;                     /* '<S90>/Delay1' */
  real_T Delay2_g;                     /* '<S90>/Delay2' */
  real_T Delay3_c;                     /* '<S90>/Delay3' */
  real_T Delay4_o;                     /* '<S90>/Delay4' */
  real_T Delay5_n;                     /* '<S90>/Delay5' */
  real_T Delay6_a;                     /* '<S90>/Delay6' */
  real_T Delay7_n;                     /* '<S90>/Delay7' */
  real_T Delay8_j;                     /* '<S90>/Delay8' */
  real_T Delay9_g;                     /* '<S90>/Delay9' */
  real_T Min_j;                        /* '<S90>/Min' */
  real_T Delay10_m;                    /* '<S91>/Delay10' */
  real_T Delay11_m;                    /* '<S91>/Delay11' */
  real_T Delay12_e;                    /* '<S91>/Delay12' */
  real_T Delay13_f;                    /* '<S91>/Delay13' */
  real_T Delay14_m;                    /* '<S91>/Delay14' */
  real_T Delay15_c;                    /* '<S91>/Delay15' */
  real_T Delay16_f;                    /* '<S91>/Delay16' */
  real_T Delay17_c;                    /* '<S91>/Delay17' */
  real_T Delay18_k;                    /* '<S91>/Delay18' */
  real_T Delay19_d;                    /* '<S91>/Delay19' */
  real_T Min1_i2;                      /* '<S91>/Min1' */
  real_T Delay_oz;                     /* '<S91>/Delay' */
  real_T Delay1_io;                    /* '<S91>/Delay1' */
  real_T Delay2_c;                     /* '<S91>/Delay2' */
  real_T Delay3_m;                     /* '<S91>/Delay3' */
  real_T Delay4_p;                     /* '<S91>/Delay4' */
  real_T Delay5_j;                     /* '<S91>/Delay5' */
  real_T Delay6_a2;                    /* '<S91>/Delay6' */
  real_T Delay7_k;                     /* '<S91>/Delay7' */
  real_T Delay8_b;                     /* '<S91>/Delay8' */
  real_T Delay9_a;                     /* '<S91>/Delay9' */
  real_T Min_c;                        /* '<S91>/Min' */
  real_T Gain2_m;                      /* '<S11>/Gain2' */
  real_T SFunction1_o1_cw;             /* '<S101>/S-Function1' */
  real_T SFunction1_o2_jm;             /* '<S101>/S-Function1' */
  real_T SFunction1_o3_dc;             /* '<S101>/S-Function1' */
  real_T SFunction1_o4_b;              /* '<S101>/S-Function1' */
  real_T SFunction1_o1_k;              /* '<S102>/S-Function1' */
  real_T SFunction1_o2_ka;             /* '<S102>/S-Function1' */
  real_T SFunction1_o3_fn;             /* '<S102>/S-Function1' */
  real_T SFunction1_o4_bh;             /* '<S102>/S-Function1' */
  real_T SFunction1_o1_hy;             /* '<S109>/S-Function1' */
  real_T SFunction1_o2_g;              /* '<S109>/S-Function1' */
  real_T SFunction1_o3_ht;             /* '<S109>/S-Function1' */
  real_T SFunction1_o1_bu;             /* '<S110>/S-Function1' */
  real_T SFunction1_o2_id;             /* '<S110>/S-Function1' */
  real_T SFunction1_o3_bh;             /* '<S110>/S-Function1' */
  real_T MAG_X;                        /* '<S14>/Gain11' */
  real_T MAG_Z;                        /* '<S14>/Gain12' */
  real_T MAG_Y;                        /* '<S14>/Gain13' */
  real_T VEL_E;                        /* '<S14>/Gain14' */
  real_T VEL_N;                        /* '<S14>/Gain15' */
  real_T VEL_D;                        /* '<S14>/Gain16' */
  real_T yawrate;                      /* '<S14>/Gain17' */
  real_T Inlineacc;                    /* '<S14>/Gain6' */
  real_T Latelalacc;                   /* '<S14>/Gain7' */
  real_T SFunction1_o1_n1;             /* '<S103>/S-Function1' */
  real_T SFunction1_o2_h;              /* '<S103>/S-Function1' */
  real_T SFunction1_o3_ll;             /* '<S103>/S-Function1' */
  real_T PulseGenerator;               /* '<S114>/Pulse Generator' */
  real_T Gain_nx;                      /* '<S127>/Gain' */
  real_T Sum1_b;                       /* '<S127>/Sum1' */
  real_T Product5_n;                   /* '<S221>/Product5' */
  real_T Divide3_a;                    /* '<S221>/Divide3' */
  real_T Add4_g;                       /* '<S221>/Add4' */
  real_T Divide2_m;                    /* '<S221>/Divide2' */
  real_T Gain6_k;                      /* '<S221>/Gain6' */
  real_T Divide_a2;                    /* '<S221>/Divide' */
  real_T Sum2_j;                       /* '<S216>/Sum2' */
  real_T TransferFcn1;                 /* '<S217>/Transfer Fcn1' */
  real_T Add1_f;                       /* '<S217>/Add1' */
  real_T DerivativeGain;               /* '<S249>/Derivative Gain' */
  real_T Filter;                       /* '<S250>/Filter' */
  real_T SumD;                         /* '<S250>/SumD' */
  real_T IntegralGain;                 /* '<S252>/Integral Gain' */
  real_T Integrator;                   /* '<S255>/Integrator' */
  real_T FilterCoefficient;            /* '<S258>/Filter Coefficient' */
  real_T ProportionalGain;             /* '<S260>/Proportional Gain' */
  real_T Sum_g;                        /* '<S264>/Sum' */
  real_T v;                            /* '<S218>/m//s to km//h' */
  real_T Saturation_oa;                /* '<S218>/Saturation' */
  real_T VectorConcatenate_k[5];       /* '<S218>/Vector Concatenate' */
  real_T uDLookupTable;                /* '<S218>/1-D Lookup Table' */
  real_T Gain_DYC_FF;                  /* '<S218>/Gain_DYC_FF' */
  real_T Divide_b5;                    /* '<S218>/Divide' */
  real_T v_b;                          /* '<S219>/m//s to km//h' */
  real_T Saturation_k;                 /* '<S219>/Saturation' */
  real_T VectorConcatenate_gm[5];      /* '<S219>/Vector Concatenate' */
  real_T uDLookupTable_f;              /* '<S219>/1-D Lookup Table' */
  real_T TransferFcn2;                 /* '<S219>/Transfer Fcn2' */
  real_T I;                            /* '<S219>/I' */
  real_T FF_diff_Gain;                 /* '<S219>/FF_diff_Gain' */
  real_T Divide_m3;                    /* '<S219>/Divide' */
  real_T Sum_f;                        /* '<S216>/Sum' */
  real_T Sum1_d;                       /* '<S216>/Sum1' */
  real_T Square1;                      /* '<S220>/Square1' */
  real_T Square2;                      /* '<S220>/Square2' */
  real_T Add_f;                        /* '<S220>/Add' */
  real_T vms_f;                        /* '<S220>/Sqrt' */
  real_T RR;                           /* '<S175>/Gain' */
  real_T RL;                           /* '<S175>/Gain1' */
  real_T FR;                           /* '<S175>/Gain2' */
  real_T FL;                           /* '<S175>/Gain3' */
  real_T Igain;                        /* '<S291>/I gain' */
  real_T Igain_b;                      /* '<S292>/I gain' */
  real_T Igain_a;                      /* '<S293>/I gain' */
  real_T Igain_j;                      /* '<S294>/I gain' */
  real_T SFunction1_n;                 /* '<S341>/S-Function1' */
  real_T Product_n;                    /* '<S338>/Product' */
  real_T reikyaku_temp1;               /* '<S338>/Add' */
  real_T LCDnumber;                    /* '<S16>/Steer Chart' */
  real_T gain1;                        /* '<S16>/Steer Chart' */
  real_T gain2;                        /* '<S16>/Steer Chart' */
  real_T Gain_h;                       /* '<S294>/Gain' */
  real_T Gain_a;                       /* '<S293>/Gain' */
  real_T Gain_b;                       /* '<S292>/Gain' */
  real_T Gain_f;                       /* '<S291>/Gain' */
  real_T LaunchLED;                    /* '<S172>/Chart' */
  real_T LaunchTorqueLimit;            /* '<S172>/Chart' */
  real_T torque_ga;                    /* '<S215>/Add' */
  real_T torque_lj;                    /* '<S214>/Add' */
  real_T torque_f1;                    /* '<S213>/Add' */
  real_T torque_lp;                    /* '<S212>/Add' */
  real_T RRtrq_c;                      /* '<S205>/Switch4' */
  real_T Saturation_ch;                /* '<S205>/Saturation' */
  real_T FLtrq_k;                      /* '<S205>/Switch7' */
  real_T Saturation3_b;                /* '<S205>/Saturation3' */
  real_T FRtrq_i;                      /* '<S205>/Switch6' */
  real_T Saturation2_f;                /* '<S205>/Saturation2' */
  real_T RLtrq_c;                      /* '<S205>/Switch5' */
  real_T Saturation1_pl;               /* '<S205>/Saturation1' */
  real_T Max;                          /* '<S183>/Max' */
  real_T Min_b;                        /* '<S183>/Min' */
  real_T Max3;                         /* '<S183>/Max3' */
  real_T Min3_h;                       /* '<S183>/Min3' */
  real_T Max2;                         /* '<S183>/Max2' */
  real_T Min2_a;                       /* '<S183>/Min2' */
  real_T Max1;                         /* '<S183>/Max1' */
  real_T Min1_l;                       /* '<S183>/Min1' */
  real_T Divide3_a5;                   /* '<S121>/Divide3' */
  real_T Divide1_j;                    /* '<S121>/Divide1' */
  real_T Fcn_h;                        /* '<S121>/Fcn' */
  real_T Divide_jo;                    /* '<S121>/Divide' */
  real_T Divide2_n;                    /* '<S121>/Divide2' */
  real_T Sum3_hp;                      /* '<S76>/Sum3' */
  real_T Sum1_fo;                      /* '<S76>/Sum1' */
  real_T Divide_ef;                    /* '<S76>/Divide' */
  real_T Saturation2_m;                /* '<S121>/Saturation2' */
  real_T Add4_c;                       /* '<S121>/Add4' */
  real_T Switch3_ad2;                  /* '<S121>/Switch3' */
  real_T Gain_a2;                      /* '<S121>/Gain' */
  real_T Divide8;                      /* '<S121>/Divide8' */
  real_T Saturation_fq;                /* '<S121>/Saturation' */
  real_T Divide9;                      /* '<S121>/Divide9' */
  real_T Add2_j;                       /* '<S121>/Add2' */
  real_T Add1_jj;                      /* '<S121>/Add1' */
  real_T totaltrqreq_g;                /* '<S121>/Switch1' */
  real_T Add5;                         /* '<S121>/Add5' */
  real_T Divide4_i;                    /* '<S121>/Divide4' */
  real_T Divide6;                      /* '<S121>/Divide6' */
  real_T Divide1_g;                    /* '<S73>/Divide1' */
  real_T LeftMax_deg;                  /* '<S73>/LeftMax_deg' */
  real_T Divide_kt;                    /* '<S73>/Divide' */
  real_T Gain1_n;                      /* '<S73>/Gain1' */
  real_T RightMax_deg;                 /* '<S73>/RightMax_deg' */
  real_T Divide1_hi;                   /* '<S72>/Divide1' */
  real_T Gain1_f;                      /* '<S72>/Gain1' */
  real_T RightMax_deg_d;               /* '<S72>/RightMax_deg' */
  real_T Divide_jq;                    /* '<S72>/Divide' */
  real_T LeftMax_deg_a;                /* '<S72>/LeftMax_deg' */
  real32_T DataTypeConversion_o;       /* '<S81>/Data Type Conversion' */
  real32_T MovingAverage;              /* '<S81>/Moving Average' */
  uint8_T LCDtext[8];                  /* '<S16>/Steer Chart' */
  boolean_T LowerRelop1;               /* '<S309>/LowerRelop1' */
  boolean_T UpperRelop;                /* '<S309>/UpperRelop' */
  boolean_T LaunchSW;                  /* '<S11>/LaunchSW' */
  boolean_T BrakeSW;                   /* '<S11>/Brake SW' */
  boolean_T Compare;                   /* '<S83>/Compare' */
  boolean_T Delay_b;                   /* '<S82>/Delay' */
  boolean_T Delay1_c;                  /* '<S82>/Delay1' */
  boolean_T Delay2_p;                  /* '<S82>/Delay2' */
  boolean_T Delay3_c1;                 /* '<S82>/Delay3' */
  boolean_T Delay4_or;                 /* '<S82>/Delay4' */
  boolean_T LogicalOperator;           /* '<S82>/Logical Operator' */
  boolean_T LogicalOperator3;          /* '<S67>/Logical Operator3' */
  boolean_T Compare_n;                 /* '<S88>/Compare' */
  boolean_T LogicalOperator1;          /* '<S67>/Logical Operator1' */
  boolean_T AND;                       /* '<S67>/AND' */
  boolean_T Compare_d;                 /* '<S169>/Compare' */
  boolean_T LowerRelop1_n;             /* '<S168>/LowerRelop1' */
  boolean_T UpperRelop_h;              /* '<S168>/UpperRelop' */
  boolean_T Compare_k;                 /* '<S128>/Compare' */
  boolean_T LogicalOperator_e;         /* '<S127>/Logical Operator' */
  boolean_T Delay_p;                   /* '<S127>/Delay' */
  boolean_T LogicalOperator2;          /* '<S127>/Logical Operator2' */
  boolean_T LogicalOperator1_o;        /* '<S127>/Logical Operator1' */
  boolean_T AllInverterEnable;         /* '<S126>/All Inverter Enable' */
  boolean_T LowerRelop1_p;             /* '<S282>/LowerRelop1' */
  boolean_T UpperRelop_f;              /* '<S282>/UpperRelop' */
  boolean_T Compare_c;                 /* '<S93>/Compare' */
  boolean_T Compare_o;                 /* '<S94>/Compare' */
  boolean_T LogicalOperator_g;         /* '<S77>/Logical Operator' */
  boolean_T Delay_ps;                  /* '<S92>/Delay' */
  boolean_T Delay1_b;                  /* '<S92>/Delay1' */
  boolean_T Delay2_o;                  /* '<S92>/Delay2' */
  boolean_T Delay3_b;                  /* '<S92>/Delay3' */
  boolean_T Delay4_m;                  /* '<S92>/Delay4' */
  boolean_T Delay5_nu;                 /* '<S92>/Delay5' */
  boolean_T Delay6_j;                  /* '<S92>/Delay6' */
  boolean_T Delay7_c;                  /* '<S92>/Delay7' */
  boolean_T Delay8_f;                  /* '<S92>/Delay8' */
  boolean_T Delay9_k;                  /* '<S92>/Delay9' */
  boolean_T LogicalOperator_i;         /* '<S92>/Logical Operator' */
  boolean_T AND_i;                     /* '<S77>/AND' */
  boolean_T All_enable;                /* '<S116>/All_enable' */
  boolean_T LogicalOperator_k;         /* '<S212>/Logical Operator' */
  boolean_T LowerRelop1_e;             /* '<S273>/LowerRelop1' */
  boolean_T UpperRelop_n;              /* '<S273>/UpperRelop' */
  boolean_T Compare_cz;                /* '<S305>/Compare' */
  boolean_T LogicalOperator_m;         /* '<S213>/Logical Operator' */
  boolean_T LowerRelop1_ej;            /* '<S274>/LowerRelop1' */
  boolean_T UpperRelop_p;              /* '<S274>/UpperRelop' */
  boolean_T Compare_p;                 /* '<S306>/Compare' */
  boolean_T LogicalOperator_d;         /* '<S214>/Logical Operator' */
  boolean_T LowerRelop1_g;             /* '<S275>/LowerRelop1' */
  boolean_T UpperRelop_f5;             /* '<S275>/UpperRelop' */
  boolean_T Compare_p1;                /* '<S307>/Compare' */
  boolean_T LogicalOperator_ex;        /* '<S215>/Logical Operator' */
  boolean_T LowerRelop1_b;             /* '<S276>/LowerRelop1' */
  boolean_T UpperRelop_ff;             /* '<S276>/UpperRelop' */
  boolean_T Compare_j;                 /* '<S308>/Compare' */
  boolean_T LogicalOperator_l;         /* '<S290>/Logical Operator' */
  boolean_T LogicalOperator1_f;        /* '<S290>/Logical Operator1' */
  boolean_T Compare_f;                 /* '<S300>/Compare' */
  boolean_T LowerRelop1_bx;            /* '<S317>/LowerRelop1' */
  boolean_T UpperRelop_l;              /* '<S317>/UpperRelop' */
  boolean_T LowerRelop1_pb;            /* '<S289>/LowerRelop1' */
  boolean_T LowerRelop1_d;             /* '<S301>/LowerRelop1' */
  boolean_T UpperRelop_m;              /* '<S301>/UpperRelop' */
  boolean_T LowerRelop1_ez;            /* '<S318>/LowerRelop1' */
  boolean_T UpperRelop_fr;             /* '<S318>/UpperRelop' */
  boolean_T UpperRelop_g;              /* '<S289>/UpperRelop' */
  boolean_T Compare_jm;                /* '<S299>/Compare' */
  boolean_T LowerRelop1_l;             /* '<S311>/LowerRelop1' */
  boolean_T UpperRelop_c;              /* '<S311>/UpperRelop' */
  boolean_T LowerRelop1_m;             /* '<S286>/LowerRelop1' */
  boolean_T LowerRelop1_mw;            /* '<S303>/LowerRelop1' */
  boolean_T UpperRelop_pt;             /* '<S303>/UpperRelop' */
  boolean_T LowerRelop1_f;             /* '<S312>/LowerRelop1' */
  boolean_T UpperRelop_d;              /* '<S312>/UpperRelop' */
  boolean_T UpperRelop_gu;             /* '<S286>/UpperRelop' */
  boolean_T LowerRelop1_k;             /* '<S313>/LowerRelop1' */
  boolean_T UpperRelop_b;              /* '<S313>/UpperRelop' */
  boolean_T LowerRelop1_mx;            /* '<S287>/LowerRelop1' */
  boolean_T LowerRelop1_n1;            /* '<S314>/LowerRelop1' */
  boolean_T UpperRelop_o;              /* '<S314>/UpperRelop' */
  boolean_T UpperRelop_n2;             /* '<S287>/UpperRelop' */
  boolean_T LowerRelop1_h;             /* '<S315>/LowerRelop1' */
  boolean_T UpperRelop_ba;             /* '<S315>/UpperRelop' */
  boolean_T LowerRelop1_gs;            /* '<S288>/LowerRelop1' */
  boolean_T LowerRelop1_hq;            /* '<S316>/LowerRelop1' */
  boolean_T UpperRelop_m4;             /* '<S316>/UpperRelop' */
  boolean_T UpperRelop_gj;             /* '<S288>/UpperRelop' */
  boolean_T Compare_b;                 /* '<S206>/Compare' */
  boolean_T Compare_pg;                /* '<S210>/Compare' */
  boolean_T Compare_a;                 /* '<S209>/Compare' */
  boolean_T Compare_h;                 /* '<S207>/Compare' */
  boolean_T Compare_cc;                /* '<S208>/Compare' */
  boolean_T Compare_f3;                /* '<S327>/Compare' */
  boolean_T Compare_m;                 /* '<S329>/Compare' */
  boolean_T Compare_pp;                /* '<S331>/Compare' */
  boolean_T Compare_fx;                /* '<S333>/Compare' */
  boolean_T VDCSW;                     /* '<S11>/VDCSW' */
  boolean_T VDCStatus;                 /* '<S320>/Logical Operator1' */
  boolean_T NOT;                       /* '<S67>/NOT' */
  boolean_T SFunction1_b;              /* '<S95>/S-Function1' */
  boolean_T SFunction1_o;              /* '<S96>/S-Function1' */
  boolean_T SFunction1_j;              /* '<S97>/S-Function1' */
  boolean_T All_QuitDcOn;              /* '<S116>/All_QuitDcOn' */
  boolean_T RtoDSW;                    /* '<S11>/RtoDSW' */
  boolean_T ReadyToDrive;              /* '<S116>/ReadyToDrive' */
  boolean_T ReadyToDrive_or;           /* '<S116>/ReadyToDrive_or' */
  boolean_T All_Error;                 /* '<S116>/All_Error' */
  boolean_T ReadyToDrive_and;          /* '<S116>/ReadyToDrive_and' */
  boolean_T All_sbm;                   /* '<S116>/All_sbm' */
  boolean_T Compare_jr;                /* '<S130>/Compare' */
  boolean_T DcOn;                      /* '<S116>/DcOn' */
  boolean_T ErrorResetSW;              /* '<S11>/ErrorRessetSW' */
  boolean_T ErrorResetSW_f;            /* '<S11>/ErrorResetSW1' */
  boolean_T LowerRelop1_i;             /* '<S142>/LowerRelop1' */
  boolean_T UpperRelop_pg;             /* '<S142>/UpperRelop' */
  boolean_T LowerRelop1_nq;            /* '<S141>/LowerRelop1' */
  boolean_T UpperRelop_db;             /* '<S141>/UpperRelop' */
  boolean_T LowerRelop1_mv;            /* '<S150>/LowerRelop1' */
  boolean_T UpperRelop_pv;             /* '<S150>/UpperRelop' */
  boolean_T LowerRelop1_g3;            /* '<S149>/LowerRelop1' */
  boolean_T UpperRelop_cy;             /* '<S149>/UpperRelop' */
  boolean_T LowerRelop1_j;             /* '<S158>/LowerRelop1' */
  boolean_T UpperRelop_lf;             /* '<S158>/UpperRelop' */
  boolean_T LowerRelop1_mh;            /* '<S157>/LowerRelop1' */
  boolean_T UpperRelop_lq;             /* '<S157>/UpperRelop' */
  boolean_T LowerRelop1_ig;            /* '<S166>/LowerRelop1' */
  boolean_T UpperRelop_md;             /* '<S166>/UpperRelop' */
  boolean_T LowerRelop1_pi;            /* '<S165>/LowerRelop1' */
  boolean_T UpperRelop_lc;             /* '<S165>/UpperRelop' */
  boolean_T MinusSW;                   /* '<S11>/MinusSW' */
  boolean_T PlusSW;                    /* '<S11>/PlusSW' */
  boolean_T SelectSW;                  /* '<S11>/SelectSW' */
  boolean_T NOT_b;                     /* '<S16>/NOT' */
  boolean_T NOT1;                      /* '<S16>/NOT1' */
  boolean_T NOT2;                      /* '<S16>/NOT2' */
  boolean_T Compare_np;                /* '<S84>/Compare' */
  boolean_T Compare_ky;                /* '<S86>/Compare' */
  boolean_T LogicalOperator_b;         /* '<S67>/Logical Operator' */
  boolean_T Delay_k;                   /* '<S79>/Delay' */
  boolean_T Delay1_l;                  /* '<S79>/Delay1' */
  boolean_T Delay2_l;                  /* '<S79>/Delay2' */
  boolean_T Delay3_o;                  /* '<S79>/Delay3' */
  boolean_T Delay4_g;                  /* '<S79>/Delay4' */
  boolean_T Delay5_e;                  /* '<S79>/Delay5' */
  boolean_T Delay6_b;                  /* '<S79>/Delay6' */
  boolean_T Delay7_cb;                 /* '<S79>/Delay7' */
  boolean_T Delay8_c;                  /* '<S79>/Delay8' */
  boolean_T Delay9_c;                  /* '<S79>/Delay9' */
  boolean_T LogicalOperator_n;         /* '<S79>/Logical Operator' */
  boolean_T Compare_fi;                /* '<S85>/Compare' */
  boolean_T Compare_g;                 /* '<S87>/Compare' */
  boolean_T LogicalOperator2_n;        /* '<S67>/Logical Operator2' */
  boolean_T Delay_d;                   /* '<S80>/Delay' */
  boolean_T Delay1_lz;                 /* '<S80>/Delay1' */
  boolean_T Delay2_a;                  /* '<S80>/Delay2' */
  boolean_T Delay3_a;                  /* '<S80>/Delay3' */
  boolean_T Delay4_gj;                 /* '<S80>/Delay4' */
  boolean_T Delay5_h;                  /* '<S80>/Delay5' */
  boolean_T Delay6_jm;                 /* '<S80>/Delay6' */
  boolean_T Delay7_f;                  /* '<S80>/Delay7' */
  boolean_T Delay8_l;                  /* '<S80>/Delay8' */
  boolean_T Delay9_d;                  /* '<S80>/Delay9' */
  boolean_T LogicalOperator_j;         /* '<S80>/Logical Operator' */
  boolean_T LogicalOperator1_b;        /* '<S77>/Logical Operator1' */
  boolean_T LogicalOperator1_j;        /* '<S114>/Logical Operator1' */
  boolean_T Counter;                   /* '<S114>/Counter' */
  boolean_T Delay_mz;                  /* '<S114>/Delay' */
  boolean_T LogicalOperator4;          /* '<S114>/Logical Operator4' */
  boolean_T LogicalOperator3_j;        /* '<S114>/Logical Operator3' */
  boolean_T LogicalOperator2_m;        /* '<S114>/Logical Operator2' */
  boolean_T LogicalOperator_ij;        /* '<S114>/Logical Operator' */
  boolean_T Compare_bz;                /* '<S129>/Compare' */
  boolean_T LogicalOperator3_d;        /* '<S127>/Logical Operator3' */
  B_MovingAverage_VCM20_T MovingAverage_pnaev;/* '<S70>/Moving Average' */
  B_MovingAverage_VCM20_T MovingAverage_pnae;/* '<S70>/Moving Average' */
  B_MovingAverage_VCM20_T MovingAverage_pna;/* '<S70>/Moving Average' */
  B_MovingAverage_VCM20_T MovingAverage_pn;/* '<S70>/Moving Average' */
  B_MovingAverage_VCM20_T MovingAverage_p;/* '<S70>/Moving Average' */
} B_VCM20_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  dsp_simulink_MovingAverage_b_T obj;  /* '<S81>/Moving Average' */
  real_T Delay_DSTATE;                 /* '<S296>/Delay' */
  real_T DiscreteTimeIntegrator_DSTATE;/* '<S294>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_e;/* '<S291>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_m;/* '<S292>/Discrete-Time Integrator' */
  real_T DiscreteTimeIntegrator_DSTATE_o;/* '<S293>/Discrete-Time Integrator' */
  real_T Delay_DSTATE_h;               /* '<S138>/Delay' */
  real_T Delay_DSTATE_k;               /* '<S140>/Delay' */
  real_T Delay_DSTATE_b;               /* '<S139>/Delay' */
  real_T Delay_DSTATE_e;               /* '<S146>/Delay' */
  real_T Delay_DSTATE_p;               /* '<S148>/Delay' */
  real_T Delay_DSTATE_i;               /* '<S147>/Delay' */
  real_T Delay_DSTATE_bc;              /* '<S154>/Delay' */
  real_T Delay_DSTATE_kj;              /* '<S156>/Delay' */
  real_T Delay_DSTATE_hs;              /* '<S155>/Delay' */
  real_T Delay_DSTATE_pw;              /* '<S162>/Delay' */
  real_T Delay_DSTATE_hw;              /* '<S164>/Delay' */
  real_T Delay_DSTATE_m;               /* '<S163>/Delay' */
  real_T Delay11_DSTATE;               /* '<S89>/Delay11' */
  real_T Delay12_DSTATE[2];            /* '<S89>/Delay12' */
  real_T Delay13_DSTATE[3];            /* '<S89>/Delay13' */
  real_T Delay14_DSTATE[4];            /* '<S89>/Delay14' */
  real_T Delay15_DSTATE[5];            /* '<S89>/Delay15' */
  real_T Delay16_DSTATE[6];            /* '<S89>/Delay16' */
  real_T Delay17_DSTATE[7];            /* '<S89>/Delay17' */
  real_T Delay18_DSTATE[8];            /* '<S89>/Delay18' */
  real_T Delay19_DSTATE[9];            /* '<S89>/Delay19' */
  real_T Delay1_DSTATE;                /* '<S89>/Delay1' */
  real_T Delay2_DSTATE[2];             /* '<S89>/Delay2' */
  real_T Delay3_DSTATE[3];             /* '<S89>/Delay3' */
  real_T Delay4_DSTATE[4];             /* '<S89>/Delay4' */
  real_T Delay5_DSTATE[5];             /* '<S89>/Delay5' */
  real_T Delay6_DSTATE[6];             /* '<S89>/Delay6' */
  real_T Delay7_DSTATE[7];             /* '<S89>/Delay7' */
  real_T Delay8_DSTATE[8];             /* '<S89>/Delay8' */
  real_T Delay9_DSTATE[9];             /* '<S89>/Delay9' */
  real_T Delay11_DSTATE_o;             /* '<S90>/Delay11' */
  real_T Delay12_DSTATE_i[2];          /* '<S90>/Delay12' */
  real_T Delay13_DSTATE_p[3];          /* '<S90>/Delay13' */
  real_T Delay14_DSTATE_k[4];          /* '<S90>/Delay14' */
  real_T Delay15_DSTATE_n[5];          /* '<S90>/Delay15' */
  real_T Delay16_DSTATE_c[6];          /* '<S90>/Delay16' */
  real_T Delay17_DSTATE_e[7];          /* '<S90>/Delay17' */
  real_T Delay18_DSTATE_a[8];          /* '<S90>/Delay18' */
  real_T Delay19_DSTATE_h[9];          /* '<S90>/Delay19' */
  real_T Delay1_DSTATE_i;              /* '<S90>/Delay1' */
  real_T Delay2_DSTATE_p[2];           /* '<S90>/Delay2' */
  real_T Delay3_DSTATE_o[3];           /* '<S90>/Delay3' */
  real_T Delay4_DSTATE_i[4];           /* '<S90>/Delay4' */
  real_T Delay5_DSTATE_m[5];           /* '<S90>/Delay5' */
  real_T Delay6_DSTATE_p[6];           /* '<S90>/Delay6' */
  real_T Delay7_DSTATE_k[7];           /* '<S90>/Delay7' */
  real_T Delay8_DSTATE_j[8];           /* '<S90>/Delay8' */
  real_T Delay9_DSTATE_a[9];           /* '<S90>/Delay9' */
  real_T Delay11_DSTATE_k;             /* '<S91>/Delay11' */
  real_T Delay12_DSTATE_h[2];          /* '<S91>/Delay12' */
  real_T Delay13_DSTATE_j[3];          /* '<S91>/Delay13' */
  real_T Delay14_DSTATE_l[4];          /* '<S91>/Delay14' */
  real_T Delay15_DSTATE_i[5];          /* '<S91>/Delay15' */
  real_T Delay16_DSTATE_b[6];          /* '<S91>/Delay16' */
  real_T Delay17_DSTATE_n[7];          /* '<S91>/Delay17' */
  real_T Delay18_DSTATE_h[8];          /* '<S91>/Delay18' */
  real_T Delay19_DSTATE_n[9];          /* '<S91>/Delay19' */
  real_T Delay1_DSTATE_h;              /* '<S91>/Delay1' */
  real_T Delay2_DSTATE_l[2];           /* '<S91>/Delay2' */
  real_T Delay3_DSTATE_o3[3];          /* '<S91>/Delay3' */
  real_T Delay4_DSTATE_il[4];          /* '<S91>/Delay4' */
  real_T Delay5_DSTATE_mb[5];          /* '<S91>/Delay5' */
  real_T Delay6_DSTATE_l[6];           /* '<S91>/Delay6' */
  real_T Delay7_DSTATE_n[7];           /* '<S91>/Delay7' */
  real_T Delay8_DSTATE_jc[8];          /* '<S91>/Delay8' */
  real_T Delay9_DSTATE_e[9];           /* '<S91>/Delay9' */
  real_T m_bpLambda;                   /* '<S218>/1-D Lookup Table' */
  real_T m_yyA;                        /* '<S218>/1-D Lookup Table' */
  real_T m_yyB;                        /* '<S218>/1-D Lookup Table' */
  real_T m_yy2;                        /* '<S218>/1-D Lookup Table' */
  real_T m_up[5];                      /* '<S218>/1-D Lookup Table' */
  real_T m_y2[5];                      /* '<S218>/1-D Lookup Table' */
  real_T prevBp0AndTableData[10];      /* '<S218>/1-D Lookup Table' */
  real_T m_bpLambda_h;                 /* '<S219>/1-D Lookup Table' */
  real_T m_yyA_j;                      /* '<S219>/1-D Lookup Table' */
  real_T m_yyB_m;                      /* '<S219>/1-D Lookup Table' */
  real_T m_yy2_p;                      /* '<S219>/1-D Lookup Table' */
  real_T m_up_m[5];                    /* '<S219>/1-D Lookup Table' */
  real_T m_y2_k[5];                    /* '<S219>/1-D Lookup Table' */
  real_T prevBp0AndTableData_m[10];    /* '<S219>/1-D Lookup Table' */
  real_T gain1max;                     /* '<S16>/Steer Chart' */
  real_T gain1min;                     /* '<S16>/Steer Chart' */
  real_T gain2max;                     /* '<S16>/Steer Chart' */
  real_T gain2min;                     /* '<S16>/Steer Chart' */
  void* m_bpDataSet;                   /* '<S218>/1-D Lookup Table' */
  void* TWork[6];                      /* '<S218>/1-D Lookup Table' */
  void* SWork[9];                      /* '<S218>/1-D Lookup Table' */
  void* m_bpDataSet_o;                 /* '<S219>/1-D Lookup Table' */
  void* TWork_c[6];                    /* '<S219>/1-D Lookup Table' */
  void* SWork_m[9];                    /* '<S219>/1-D Lookup Table' */
  int32_T clockTickCounter;            /* '<S114>/Pulse Generator' */
  uint32_T Counter_ClkEphState;        /* '<S114>/Counter' */
  uint32_T Counter_RstEphState;        /* '<S114>/Counter' */
  uint32_T m_bpIndex;                  /* '<S218>/1-D Lookup Table' */
  uint32_T m_bpIndex_h;                /* '<S219>/1-D Lookup Table' */
  uint32_T temporalCounter_i1;         /* '<S172>/Chart' */
  int_T SFunction1_IWORK[2];           /* '<S6>/S-Function1' */
  int_T SFunction1_IWORK_n[2];         /* '<S7>/S-Function1' */
  int_T SFunction1_IWORK_e[2];         /* '<S8>/S-Function1' */
  int_T SFunction1_IWORK_o[2];         /* '<S9>/S-Function1' */
  int_T SFunction1_IWORK_a[2];         /* '<S10>/S-Function1' */
  int_T SFunction1_IWORK_j[2];         /* '<S17>/S-Function1' */
  boolean_T Delay1_DSTATE_c;           /* '<S82>/Delay1' */
  boolean_T Delay2_DSTATE_a[2];        /* '<S82>/Delay2' */
  boolean_T Delay3_DSTATE_p[3];        /* '<S82>/Delay3' */
  boolean_T Delay4_DSTATE_o[4];        /* '<S82>/Delay4' */
  boolean_T Delay_DSTATE_d;            /* '<S127>/Delay' */
  boolean_T Delay1_DSTATE_ci;          /* '<S92>/Delay1' */
  boolean_T Delay2_DSTATE_j[2];        /* '<S92>/Delay2' */
  boolean_T Delay3_DSTATE_g[3];        /* '<S92>/Delay3' */
  boolean_T Delay4_DSTATE_g[4];        /* '<S92>/Delay4' */
  boolean_T Delay5_DSTATE_g[5];        /* '<S92>/Delay5' */
  boolean_T Delay6_DSTATE_lu[6];       /* '<S92>/Delay6' */
  boolean_T Delay7_DSTATE_nq[7];       /* '<S92>/Delay7' */
  boolean_T Delay8_DSTATE_k[8];        /* '<S92>/Delay8' */
  boolean_T Delay9_DSTATE_o[9];        /* '<S92>/Delay9' */
  boolean_T Delay1_DSTATE_p;           /* '<S79>/Delay1' */
  boolean_T Delay2_DSTATE_e[2];        /* '<S79>/Delay2' */
  boolean_T Delay3_DSTATE_n[3];        /* '<S79>/Delay3' */
  boolean_T Delay4_DSTATE_j[4];        /* '<S79>/Delay4' */
  boolean_T Delay5_DSTATE_d[5];        /* '<S79>/Delay5' */
  boolean_T Delay6_DSTATE_i[6];        /* '<S79>/Delay6' */
  boolean_T Delay7_DSTATE_g[7];        /* '<S79>/Delay7' */
  boolean_T Delay8_DSTATE_c[8];        /* '<S79>/Delay8' */
  boolean_T Delay9_DSTATE_h[9];        /* '<S79>/Delay9' */
  boolean_T Delay1_DSTATE_f;           /* '<S80>/Delay1' */
  boolean_T Delay2_DSTATE_k[2];        /* '<S80>/Delay2' */
  boolean_T Delay3_DSTATE_gb[3];       /* '<S80>/Delay3' */
  boolean_T Delay4_DSTATE_a[4];        /* '<S80>/Delay4' */
  boolean_T Delay5_DSTATE_k[5];        /* '<S80>/Delay5' */
  boolean_T Delay6_DSTATE_n[6];        /* '<S80>/Delay6' */
  boolean_T Delay7_DSTATE_m[7];        /* '<S80>/Delay7' */
  boolean_T Delay8_DSTATE_d[8];        /* '<S80>/Delay8' */
  boolean_T Delay9_DSTATE_p[9];        /* '<S80>/Delay9' */
  boolean_T Delay_DSTATE_f;            /* '<S114>/Delay' */
  uint8_T Counter_Count;               /* '<S114>/Counter' */
  uint8_T reCalcSecDerivFirstDimCoeffs;/* '<S218>/1-D Lookup Table' */
  uint8_T reCalcSecDerivFirstDimCoeffs_i;/* '<S219>/1-D Lookup Table' */
  uint8_T is_active_c3_VCM20;          /* '<S16>/Steer Chart' */
  uint8_T is_c3_VCM20;                 /* '<S16>/Steer Chart' */
  uint8_T is_active_c1_VCM20;          /* '<S172>/Chart' */
  uint8_T is_c1_VCM20;                 /* '<S172>/Chart' */
  uint8_T is_LaunchReady_BrakeOff;     /* '<S172>/Chart' */
  boolean_T objisempty;                /* '<S81>/Moving Average' */
  DW_MovingAverage_VCM20_T MovingAverage_pnaev;/* '<S70>/Moving Average' */
  DW_MovingAverage_VCM20_T MovingAverage_pnae;/* '<S70>/Moving Average' */
  DW_MovingAverage_VCM20_T MovingAverage_pna;/* '<S70>/Moving Average' */
  DW_MovingAverage_VCM20_T MovingAverage_pn;/* '<S70>/Moving Average' */
  DW_MovingAverage_VCM20_T MovingAverage_p;/* '<S70>/Moving Average' */
} DW_VCM20_T;

/* Continuous states (default storage) */
typedef struct {
  real_T TransferFcn_CSTATE;           /* '<S174>/Transfer Fcn' */
  real_T TransferFcn_CSTATE_k;         /* '<S73>/Transfer Fcn' */
  real_T TransferFcn_CSTATE_j;         /* '<S220>/Transfer Fcn' */
  real_T TransferFcn_CSTATE_e;         /* '<S72>/Transfer Fcn' */
  real_T TransferFcn1_CSTATE;          /* '<S217>/Transfer Fcn1' */
  real_T Filter_CSTATE;                /* '<S250>/Filter' */
  real_T Integrator_CSTATE;            /* '<S255>/Integrator' */
  real_T TransferFcn2_CSTATE;          /* '<S219>/Transfer Fcn2' */
} X_VCM20_T;

/* State derivatives (default storage) */
typedef struct {
  real_T TransferFcn_CSTATE;           /* '<S174>/Transfer Fcn' */
  real_T TransferFcn_CSTATE_k;         /* '<S73>/Transfer Fcn' */
  real_T TransferFcn_CSTATE_j;         /* '<S220>/Transfer Fcn' */
  real_T TransferFcn_CSTATE_e;         /* '<S72>/Transfer Fcn' */
  real_T TransferFcn1_CSTATE;          /* '<S217>/Transfer Fcn1' */
  real_T Filter_CSTATE;                /* '<S250>/Filter' */
  real_T Integrator_CSTATE;            /* '<S255>/Integrator' */
  real_T TransferFcn2_CSTATE;          /* '<S219>/Transfer Fcn2' */
} XDot_VCM20_T;

/* State disabled  */
typedef struct {
  boolean_T TransferFcn_CSTATE;        /* '<S174>/Transfer Fcn' */
  boolean_T TransferFcn_CSTATE_k;      /* '<S73>/Transfer Fcn' */
  boolean_T TransferFcn_CSTATE_j;      /* '<S220>/Transfer Fcn' */
  boolean_T TransferFcn_CSTATE_e;      /* '<S72>/Transfer Fcn' */
  boolean_T TransferFcn1_CSTATE;       /* '<S217>/Transfer Fcn1' */
  boolean_T Filter_CSTATE;             /* '<S250>/Filter' */
  boolean_T Integrator_CSTATE;         /* '<S255>/Integrator' */
  boolean_T TransferFcn2_CSTATE;       /* '<S219>/Transfer Fcn2' */
} XDis_VCM20_T;

#ifndef ODE1_INTG
#define ODE1_INTG

/* ODE1 Integration Data */
typedef struct {
  real_T *f[1];                        /* derivatives */
} ODE1_IntgData;

#endif

/* Parameters (default storage) */
struct P_VCM20_T_ {
  real_T Iwf;                          /* Variable: Iwf
                                        * Referenced by: '<S283>/wheel Inertia'
                                        */
  real_T M;                            /* Variable: M
                                        * Referenced by: '<S296>/FR Load Diff'
                                        */
  real_T Rwf;                          /* Variable: Rwf
                                        * Referenced by:
                                        *   '<S344>/toTire_rpm'
                                        *   '<S283>/1// tire R1'
                                        *   '<S283>/Rtire'
                                        *   '<S290>/toTire_rpm'
                                        */
  real_T W0_F;                         /* Variable: W0_F
                                        * Referenced by:
                                        *   '<S296>/Front_ini_Load'
                                        *   '<S296>/Saturation'
                                        */
  real_T W0_R;                         /* Variable: W0_R
                                        * Referenced by:
                                        *   '<S296>/Rear_ini_Load2'
                                        *   '<S296>/Saturation'
                                        */
  real_T afx_f[7];                     /* Variable: afx_f
                                        * Referenced by:
                                        *   '<S297>/C'
                                        *   '<S297>/a2'
                                        *   '<S297>/a4'
                                        *   '<S297>/a6'
                                        *   '<S297>/BCD'
                                        *   '<S297>/a1*W'
                                        *   '<S297>/a5*W'
                                        *   '<S298>/C'
                                        *   '<S298>/a2'
                                        *   '<S298>/a4'
                                        *   '<S298>/a6'
                                        *   '<S298>/BCD'
                                        *   '<S298>/a1*W'
                                        *   '<S298>/a5*W'
                                        */
  real_T g;                            /* Variable: g
                                        * Referenced by:
                                        *   '<S283>/Target Acc'
                                        *   '<S290>/Constant1'
                                        *   '<S296>/Front_ini_Load'
                                        *   '<S296>/Rear_ini_Load2'
                                        */
  real_T gear_f;                       /* Variable: gear_f
                                        * Referenced by:
                                        *   '<S344>/Gear'
                                        *   '<S283>/1//Gear'
                                        *   '<S283>/1//Gear1'
                                        *   '<S290>/Gear'
                                        */
  real_T gear_r;                       /* Variable: gear_r
                                        * Referenced by:
                                        *   '<S283>/1//Gear2'
                                        *   '<S283>/Rtire1'
                                        */
  real_T hg;                           /* Variable: hg
                                        * Referenced by: '<S296>/FR Load Diff'
                                        */
  real_T l;                            /* Variable: l
                                        * Referenced by: '<S296>/FR Load Diff'
                                        */
  real_T mu_tire_F;                    /* Variable: mu_tire_F
                                        * Referenced by:
                                        *   '<S297>/Gain1'
                                        *   '<S298>/Gain1'
                                        */
  real_T PIDController_D;              /* Mask Parameter: PIDController_D
                                        * Referenced by: '<S249>/Derivative Gain'
                                        */
  real_T PIDController_I;              /* Mask Parameter: PIDController_I
                                        * Referenced by: '<S252>/Integral Gain'
                                        */
  real_T PIDController_InitialConditionF;
                              /* Mask Parameter: PIDController_InitialConditionF
                               * Referenced by: '<S250>/Filter'
                               */
  real_T PIDController_InitialConditio_l;
                              /* Mask Parameter: PIDController_InitialConditio_l
                               * Referenced by: '<S255>/Integrator'
                               */
  real_T PIDController_N;              /* Mask Parameter: PIDController_N
                                        * Referenced by: '<S258>/Filter Coefficient'
                                        */
  real_T PIDController_P;              /* Mask Parameter: PIDController_P
                                        * Referenced by: '<S260>/Proportional Gain'
                                        */
  real_T CompareToConstant_const;     /* Mask Parameter: CompareToConstant_const
                                       * Referenced by: '<S83>/Constant'
                                       */
  real_T CompareToConstant1_const;   /* Mask Parameter: CompareToConstant1_const
                                      * Referenced by: '<S93>/Constant'
                                      */
  real_T CompareToConstant3_const;   /* Mask Parameter: CompareToConstant3_const
                                      * Referenced by: '<S94>/Constant'
                                      */
  real_T CompareToConstant3_const_h;
                                   /* Mask Parameter: CompareToConstant3_const_h
                                    * Referenced by: '<S210>/Constant'
                                    */
  real_T CompareToConstant2_const;   /* Mask Parameter: CompareToConstant2_const
                                      * Referenced by: '<S209>/Constant'
                                      */
  real_T CompareToConstant_const_d; /* Mask Parameter: CompareToConstant_const_d
                                     * Referenced by: '<S207>/Constant'
                                     */
  real_T CompareToConstant1_const_l;
                                   /* Mask Parameter: CompareToConstant1_const_l
                                    * Referenced by: '<S208>/Constant'
                                    */
  real_T CompareToConstant_const_l; /* Mask Parameter: CompareToConstant_const_l
                                     * Referenced by: '<S130>/Constant'
                                     */
  real_T CompareToConstant1_const_f;
                                   /* Mask Parameter: CompareToConstant1_const_f
                                    * Referenced by: '<S84>/Constant'
                                    */
  real_T CompareToConstant3_const_f;
                                   /* Mask Parameter: CompareToConstant3_const_f
                                    * Referenced by: '<S86>/Constant'
                                    */
  real_T CompareToConstant2_const_p;
                                   /* Mask Parameter: CompareToConstant2_const_p
                                    * Referenced by: '<S85>/Constant'
                                    */
  real_T CompareToConstant4_const;   /* Mask Parameter: CompareToConstant4_const
                                      * Referenced by: '<S87>/Constant'
                                      */
  real32_T CompareToConstant_const_k;
                                    /* Mask Parameter: CompareToConstant_const_k
                                     * Referenced by: '<S88>/Constant'
                                     */
  uint8_T Counter_HitValue;            /* Mask Parameter: Counter_HitValue
                                        * Referenced by: '<S114>/Counter'
                                        */
  uint8_T Counter_InitialCount;        /* Mask Parameter: Counter_InitialCount
                                        * Referenced by: '<S114>/Counter'
                                        */
  real_T LeftMax_deg_Gain;             /* Expression: 95
                                        * Referenced by: '<S72>/LeftMax_deg'
                                        */
  real_T Gain1_Gain;                   /* Expression: -1
                                        * Referenced by: '<S72>/Gain1'
                                        */
  real_T RightMax_deg_Gain;            /* Expression: 85
                                        * Referenced by: '<S72>/RightMax_deg'
                                        */
  real_T Gain1_Gain_m;                 /* Expression: -1
                                        * Referenced by: '<S73>/Gain1'
                                        */
  real_T RightMax_deg_Gain_c;          /* Expression: 95
                                        * Referenced by: '<S73>/RightMax_deg'
                                        */
  real_T LeftMax_deg_Gain_h;           /* Expression: 85
                                        * Referenced by: '<S73>/LeftMax_deg'
                                        */
  real_T Constant_Value;               /* Expression: 0
                                        * Referenced by: '<S113>/Constant'
                                        */
  real_T Constant1_Value;              /* Expression: 0
                                        * Referenced by: '<S115>/Constant1'
                                        */
  real_T allpowerlimit_Value;          /* Expression: 35000*0.6
                                        * Referenced by: '<S136>/all power limit'
                                        */
  real_T eachpowerlimit_Value;         /* Expression: 80000/4*0.7
                                        * Referenced by: '<S136>/each power limit '
                                        */
  real_T allpowerlimit_Value_g;        /* Expression: 35000*0.6
                                        * Referenced by: '<S144>/all power limit'
                                        */
  real_T eachpowerlimit_Value_b;       /* Expression: 80000/4*0.7
                                        * Referenced by: '<S144>/each power limit '
                                        */
  real_T allpowerlimit_Value_i;        /* Expression: 35000*0.6
                                        * Referenced by: '<S152>/all power limit'
                                        */
  real_T eachpowerlimit_Value_k;       /* Expression: 80000/4*0.7
                                        * Referenced by: '<S152>/each power limit '
                                        */
  real_T allpowerlimit_Value_h;        /* Expression: 35000*0.6
                                        * Referenced by: '<S160>/all power limit'
                                        */
  real_T eachpowerlimit_Value_d;       /* Expression: 80000/4*0.7
                                        * Referenced by: '<S160>/each power limit '
                                        */
  real_T Constant_Value_e;             /* Expression: eps
                                        * Referenced by: '<S167>/Constant'
                                        */
  real_T Constant_Value_p;             /* Expression: 0
                                        * Referenced by: '<S121>/Constant'
                                        */
  real_T BrakeStart_Pedalbias_Value;   /* Expression: 0.2
                                        * Referenced by: '<S15>/BrakeStart_Pedal bias'
                                        */
  real_T Switch1_Threshold;            /* Expression: 0
                                        * Referenced by: '<S121>/Switch1'
                                        */
  real_T BrakePressF_MAX_Value;        /* Expression: 515
                                        * Referenced by: '<S76>/BrakePressF_MAX'
                                        */
  real_T BrakePressF_MIN_Value;        /* Expression: 826
                                        * Referenced by: '<S76>/BrakePressF_MIN'
                                        */
  real_T Saturation2_UpperSat;         /* Expression: 0.5
                                        * Referenced by: '<S121>/Saturation2'
                                        */
  real_T Saturation2_LowerSat;         /* Expression: -0.5
                                        * Referenced by: '<S121>/Saturation2'
                                        */
  real_T BrakeEnd_pedalbias_Value;     /* Expression: 0.8
                                        * Referenced by: '<S15>/BrakeEnd_pedal bias'
                                        */
  real_T Switch3_Threshold;            /* Expression: 0
                                        * Referenced by: '<S121>/Switch3'
                                        */
  real_T Gain_Gain;                    /* Expression: -1
                                        * Referenced by: '<S121>/Gain'
                                        */
  real_T Saturation_UpperSat;          /* Expression: 1
                                        * Referenced by: '<S121>/Saturation'
                                        */
  real_T Saturation_LowerSat;          /* Expression: 0
                                        * Referenced by: '<S121>/Saturation'
                                        */
  real_T Saturation1_UpperSat;         /* Expression: 8.72
                                        * Referenced by: '<S205>/Saturation1'
                                        */
  real_T Saturation1_LowerSat;         /* Expression: -inf
                                        * Referenced by: '<S205>/Saturation1'
                                        */
  real_T Saturation2_UpperSat_b;       /* Expression: 8.72
                                        * Referenced by: '<S205>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_g;       /* Expression: -inf
                                        * Referenced by: '<S205>/Saturation2'
                                        */
  real_T Saturation3_UpperSat;         /* Expression: 8.72
                                        * Referenced by: '<S205>/Saturation3'
                                        */
  real_T Saturation3_LowerSat;         /* Expression: -inf
                                        * Referenced by: '<S205>/Saturation3'
                                        */
  real_T Saturation_UpperSat_a;        /* Expression: 8.72
                                        * Referenced by: '<S205>/Saturation'
                                        */
  real_T Saturation_LowerSat_n;        /* Expression: -inf
                                        * Referenced by: '<S205>/Saturation'
                                        */
  real_T Constant1_Value_b;            /* Expression: eps
                                        * Referenced by: '<S297>/Constant1'
                                        */
  real_T Constant1_Value_i;            /* Expression: eps
                                        * Referenced by: '<S298>/Constant1'
                                        */
  real_T Gain_Gain_h;                  /* Expression: -1
                                        * Referenced by: '<S291>/Gain'
                                        */
  real_T Gain_Gain_j;                  /* Expression: -1
                                        * Referenced by: '<S292>/Gain'
                                        */
  real_T Gain_Gain_a;                  /* Expression: -1
                                        * Referenced by: '<S293>/Gain'
                                        */
  real_T Gain_Gain_i;                  /* Expression: -1
                                        * Referenced by: '<S294>/Gain'
                                        */
  real_T Constant_Value_g;             /* Expression: eps
                                        * Referenced by: '<S326>/Constant'
                                        */
  real_T Constant_Value_c;             /* Expression: eps
                                        * Referenced by: '<S328>/Constant'
                                        */
  real_T Constant_Value_d;             /* Expression: eps
                                        * Referenced by: '<S330>/Constant'
                                        */
  real_T Constant_Value_j;             /* Expression: eps
                                        * Referenced by: '<S332>/Constant'
                                        */
  real_T Constant_Value_o;             /* Expression: 0
                                        * Referenced by: '<S128>/Constant'
                                        */
  real_T Constant_Value_a;             /* Expression: 0
                                        * Referenced by: '<S129>/Constant'
                                        */
  real_T Constant_Value_m;             /* Expression: 0
                                        * Referenced by: '<S169>/Constant'
                                        */
  real_T Constant_Value_h;             /* Expression: 0
                                        * Referenced by: '<S206>/Constant'
                                        */
  real_T Constant_Value_f;             /* Expression: 0
                                        * Referenced by: '<S299>/Constant'
                                        */
  real_T Constant_Value_cg;            /* Expression: 0
                                        * Referenced by: '<S300>/Constant'
                                        */
  real_T Constant_Value_jm;            /* Expression: 0
                                        * Referenced by: '<S305>/Constant'
                                        */
  real_T Constant_Value_ow;            /* Expression: 0
                                        * Referenced by: '<S306>/Constant'
                                        */
  real_T Constant_Value_ct;            /* Expression: 0
                                        * Referenced by: '<S307>/Constant'
                                        */
  real_T Constant_Value_n;             /* Expression: 0
                                        * Referenced by: '<S308>/Constant'
                                        */
  real_T Constant_Value_h3;            /* Expression: 0
                                        * Referenced by: '<S327>/Constant'
                                        */
  real_T Constant_Value_dl;            /* Expression: 0
                                        * Referenced by: '<S329>/Constant'
                                        */
  real_T Constant_Value_i;             /* Expression: 0
                                        * Referenced by: '<S331>/Constant'
                                        */
  real_T Constant_Value_in;            /* Expression: 0
                                        * Referenced by: '<S333>/Constant'
                                        */
  real_T u0kWOn_Value;                 /* Expression: 1
                                        * Referenced by: '<Root>/80kWOn'
                                        */
  real_T ABSOn_Value;                  /* Expression: 0
                                        * Referenced by: '<Root>/ABSOn'
                                        */
  real_T MaxTorque_Value;              /* Expression: 10
                                        * Referenced by: '<Root>/MaxTorque'
                                        */
  real_T Saturation_UpperSat_n;        /* Expression: 21
                                        * Referenced by: '<Root>/Saturation'
                                        */
  real_T Saturation_LowerSat_f;        /* Expression: 0
                                        * Referenced by: '<Root>/Saturation'
                                        */
  real_T Speedlimitationrpm_Value;     /* Expression: 20000
                                        * Referenced by: '<Root>/Speed limitation (rpm)'
                                        */
  real_T Speedlimitationkph_Value;     /* Expression: 120
                                        * Referenced by: '<Root>/Speed limitation (kph)'
                                        */
  real_T toms_Gain;                    /* Expression: 1000/3600
                                        * Referenced by: '<S344>/to m//s'
                                        */
  real_T Saturation2_UpperSat_n;       /* Expression: 19500
                                        * Referenced by: '<Root>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_b;       /* Expression: 0
                                        * Referenced by: '<Root>/Saturation2'
                                        */
  real_T Gain8_Gain;                   /* Expression: 10^-2
                                        * Referenced by: '<S14>/Gain8'
                                        */
  real_T Gain1_Gain_f;                 /* Expression: -1
                                        * Referenced by: '<S14>/Gain1'
                                        */
  real_T Timemsec_Gain;                /* Expression: 80
                                        * Referenced by: '<S290>/Time msec'
                                        */
  real_T Constant_Value_ia;            /* Expression: 100
                                        * Referenced by: '<S290>/Constant'
                                        */
  real_T Constant2_Value;              /* Expression: 1
                                        * Referenced by: '<S310>/Constant2'
                                        */
  real_T Rear_sliprate_ref_Value;      /* Expression: 0.18
                                        * Referenced by: '<S290>/Rear_sliprate_ref'
                                        */
  real_T Throttle1_MIN_Value;          /* Expression: 180
                                        * Referenced by: '<S70>/Throttle1_MIN'
                                        */
  real_T Throttle1_MAX_Value;          /* Expression: 517
                                        * Referenced by: '<S70>/Throttle1_MAX'
                                        */
  real_T Throttle2_MIN_Value;          /* Expression: 829
                                        * Referenced by: '<S71>/Throttle2_MIN'
                                        */
  real_T Throttle2_MAX_Value;          /* Expression: 484
                                        * Referenced by: '<S71>/Throttle2_MAX'
                                        */
  real_T Saturation_UpperSat_p;        /* Expression: 1
                                        * Referenced by: '<S67>/Saturation'
                                        */
  real_T Saturation_LowerSat_k;        /* Expression: 0
                                        * Referenced by: '<S67>/Saturation'
                                        */
  real_T Saturation1_UpperSat_d;       /* Expression: 1
                                        * Referenced by: '<S121>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_g;       /* Expression: 0
                                        * Referenced by: '<S121>/Saturation1'
                                        */
  real_T Trottle_Pedalbias_Value;      /* Expression: 0.05
                                        * Referenced by: '<S15>/Trottle_Pedal bias'
                                        */
  real_T CoastTorque_Value;            /* Expression: 0
                                        * Referenced by: '<Root>/CoastTorque'
                                        */
  real_T Saturation3_UpperSat_a;       /* Expression: 21
                                        * Referenced by: '<Root>/Saturation3'
                                        */
  real_T Saturation3_LowerSat_d;       /* Expression: 0
                                        * Referenced by: '<Root>/Saturation3'
                                        */
  real_T Switch_Threshold;             /* Expression: 0
                                        * Referenced by: '<S121>/Switch'
                                        */
  real_T Gain1_Gain_d;                 /* Expression: -1
                                        * Referenced by: '<S121>/Gain1'
                                        */
  real_T MinTorque_Value;              /* Expression: 0
                                        * Referenced by: '<Root>/MinTorque'
                                        */
  real_T Saturation1_UpperSat_m;       /* Expression: 21
                                        * Referenced by: '<Root>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_i;       /* Expression: 0
                                        * Referenced by: '<Root>/Saturation1'
                                        */
  real_T Switch2_Threshold;            /* Expression: 0
                                        * Referenced by: '<S121>/Switch2'
                                        */
  real_T Gain2_Gain;                   /* Expression: -1
                                        * Referenced by: '<S121>/Gain2'
                                        */
  real_T MC1_sw_Value;                 /* Expression: 1
                                        * Referenced by: '<Root>/MC1_sw'
                                        */
  real_T Constant_Value_pl;            /* Expression: 1
                                        * Referenced by: '<S126>/Constant'
                                        */
  real_T Switch1_Threshold_p;          /* Expression: 0
                                        * Referenced by: '<S126>/Switch1'
                                        */
  real_T MC2_sw_Value;                 /* Expression: 1
                                        * Referenced by: '<Root>/MC2_sw'
                                        */
  real_T Switch2_Threshold_n;          /* Expression: 0
                                        * Referenced by: '<S126>/Switch2'
                                        */
  real_T MC3_sw_Value;                 /* Expression: 1
                                        * Referenced by: '<Root>/MC3_sw'
                                        */
  real_T Switch3_Threshold_d;          /* Expression: 0
                                        * Referenced by: '<S126>/Switch3'
                                        */
  real_T MC4_sw_Value;                 /* Expression: 1
                                        * Referenced by: '<Root>/MC4_sw'
                                        */
  real_T Switch4_Threshold;            /* Expression: 0
                                        * Referenced by: '<S126>/Switch4'
                                        */
  real_T GearProtectionTorqueLimit_Value;/* Expression: 1
                                          * Referenced by: '<S172>/GearProtectionTorqueLimit'
                                          */
  real_T GearProtectionTimes_Value;    /* Expression: 0.5
                                        * Referenced by: '<S172>/GearProtectionTime (s)'
                                        */
  real_T LaunchOn_Value;               /* Expression: 1
                                        * Referenced by: '<Root>/LaunchOn'
                                        */
  real_T Switch_Threshold_b;           /* Expression: 0.5
                                        * Referenced by: '<S172>/Switch'
                                        */
  real_T TrqFRDistGain_Value;          /* Expression: 0.8
                                        * Referenced by: '<S174>/Trq F//R Dist Gain'
                                        */
  real_T Constant9_Value;              /* Expression: 1
                                        * Referenced by: '<S281>/Constant9'
                                        */
  real_T TransferFcn_A;                /* Computed Parameter: TransferFcn_A
                                        * Referenced by: '<S174>/Transfer Fcn'
                                        */
  real_T TransferFcn_C;                /* Computed Parameter: TransferFcn_C
                                        * Referenced by: '<S174>/Transfer Fcn'
                                        */
  real_T hgL_Gain;                     /* Expression: 0.273/1.530
                                        * Referenced by: '<S281>/hg//L'
                                        */
  real_T M_Value;                      /* Expression: 270
                                        * Referenced by: '<S281>/M'
                                        */
  real_T Gain_Gain_n;                  /* Expression: 1/9.8
                                        * Referenced by: '<S281>/Gain'
                                        */
  real_T Gain12_Gain;                  /* Expression: 0.48
                                        * Referenced by: '<S281>/Gain12'
                                        */
  real_T Saturation1_UpperSat_i;       /* Expression: 270
                                        * Referenced by: '<S281>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_gq;      /* Expression: 10
                                        * Referenced by: '<S281>/Saturation1'
                                        */
  real_T Gain13_Gain;                  /* Expression: 0.52
                                        * Referenced by: '<S281>/Gain13'
                                        */
  real_T Saturation2_UpperSat_f;       /* Expression: 270
                                        * Referenced by: '<S281>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_a;       /* Expression: 10
                                        * Referenced by: '<S281>/Saturation2'
                                        */
  real_T Constant10_Value;             /* Expression: 1
                                        * Referenced by: '<S281>/Constant10'
                                        */
  real_T Gain1_Gain_k;                 /* Expression: -1
                                        * Referenced by: '<S281>/Gain1'
                                        */
  real_T Saturation3_UpperSat_m;       /* Expression: inf
                                        * Referenced by: '<S278>/Saturation3'
                                        */
  real_T Saturation3_LowerSat_o;       /* Expression: 0
                                        * Referenced by: '<S278>/Saturation3'
                                        */
  real_T Saturation_UpperSat_c;        /* Expression: inf
                                        * Referenced by: '<S278>/Saturation'
                                        */
  real_T Saturation_LowerSat_h;        /* Expression: 0
                                        * Referenced by: '<S278>/Saturation'
                                        */
  real_T Gain_Gain_l;                  /* Expression: -1
                                        * Referenced by: '<S279>/Gain'
                                        */
  real_T Saturation3_UpperSat_i;       /* Expression: 0
                                        * Referenced by: '<S279>/Saturation3'
                                        */
  real_T Saturation3_LowerSat_oy;      /* Expression: -inf
                                        * Referenced by: '<S279>/Saturation3'
                                        */
  real_T Saturation_UpperSat_pl;       /* Expression: 0
                                        * Referenced by: '<S279>/Saturation'
                                        */
  real_T Saturation_LowerSat_p;        /* Expression: -inf
                                        * Referenced by: '<S279>/Saturation'
                                        */
  real_T FRTrqDistOn_Value;            /* Expression: 0
                                        * Referenced by: '<Root>/F//RTrqDistOn'
                                        */
  real_T Switch1_Threshold_b;          /* Expression: 0.5
                                        * Referenced by: '<S174>/Switch1'
                                        */
  real_T Gain4_Gain;                   /* Expression: 1
                                        * Referenced by: '<S174>/Gain4'
                                        */
  real_T SWAS2_signal_plus_Value;      /* Expression: 789
                                        * Referenced by: '<S73>/SWAS2_signal_plus'
                                        */
  real_T SWAS2_signal_0deg_Value;      /* Expression: 520
                                        * Referenced by: '<S73>/SWAS2_signal_0deg'
                                        */
  real_T TransferFcn_A_c;              /* Computed Parameter: TransferFcn_A_c
                                        * Referenced by: '<S73>/Transfer Fcn'
                                        */
  real_T TransferFcn_C_p;              /* Computed Parameter: TransferFcn_C_p
                                        * Referenced by: '<S73>/Transfer Fcn'
                                        */
  real_T SWAS2_signal_minus_Value;     /* Expression: 252
                                        * Referenced by: '<S73>/SWAS2_signal_minus'
                                        */
  real_T Switch_Threshold_n;           /* Expression: 0
                                        * Referenced by: '<S73>/Switch'
                                        */
  real_T Saturation_UpperSat_f;        /* Expression: 180
                                        * Referenced by: '<S77>/Saturation'
                                        */
  real_T Saturation_LowerSat_c;        /* Expression: -180
                                        * Referenced by: '<S77>/Saturation'
                                        */
  real_T SteerAngle_deg_Gain;          /* Expression: 1/5.65
                                        * Referenced by: '<S77>/SteerAngle_deg'
                                        */
  real_T SteerAngle_rad_Gain;          /* Expression: pi/180
                                        * Referenced by: '<S77>/SteerAngle_rad'
                                        */
  real_T Gain_Gain_b;                  /* Expression: 300
                                        * Referenced by: '<S216>/Gain'
                                        */
  real_T Gain12_Gain_d;                /* Expression: 1/4
                                        * Referenced by: '<S211>/Gain12'
                                        */
  real_T Gain13_Gain_f;                /* Expression: 1/0.625
                                        * Referenced by: '<S211>/Gain13'
                                        */
  real_T Gain14_Gain;                  /* Expression: 1/13.75
                                        * Referenced by: '<S211>/Gain14'
                                        */
  real_T Gain4_Gain_n;                 /* Expression: 18*(25.4/1000)/2
                                        * Referenced by: '<S211>/Gain4'
                                        */
  real_T Saturation_UpperSat_e;        /* Expression: 5
                                        * Referenced by: '<S211>/Saturation'
                                        */
  real_T Saturation_LowerSat_n4;       /* Expression: -5
                                        * Referenced by: '<S211>/Saturation'
                                        */
  real_T DYCOn_Value;                  /* Expression: 0
                                        * Referenced by: '<Root>/DYCOn'
                                        */
  real_T Constant_Value_b;             /* Expression: 1
                                        * Referenced by: '<S132>/Constant'
                                        */
  real_T Switch1_Threshold_p4;         /* Expression: 0
                                        * Referenced by: '<S132>/Switch1'
                                        */
  real_T Switch2_Threshold_d;          /* Expression: 0
                                        * Referenced by: '<S132>/Switch2'
                                        */
  real_T Switch3_Threshold_b;          /* Expression: 0
                                        * Referenced by: '<S132>/Switch3'
                                        */
  real_T Switch4_Threshold_f;          /* Expression: 0
                                        * Referenced by: '<S132>/Switch4'
                                        */
  real_T Gain5_Gain;                   /* Expression: 1
                                        * Referenced by: '<S174>/Gain5'
                                        */
  real_T Switch_Threshold_l;           /* Expression: 0.5
                                        * Referenced by: '<S174>/Switch'
                                        */
  real_T Gain10_Gain;                  /* Expression: 1
                                        * Referenced by: '<S174>/Gain10'
                                        */
  real_T Gain11_Gain;                  /* Expression: 1
                                        * Referenced by: '<S174>/Gain11'
                                        */
  real_T Constant2_Value_e;            /* Expression: 0
                                        * Referenced by: '<S290>/Constant2'
                                        */
  real_T Fd_tableData[9];              /* Expression: [0.55
                                          0.6
                                          0.7
                                          0.8
                                          1
                                          1.5
                                          2
                                          2.5
                                          3]
                                        * Referenced by: '<S283>/Fd'
                                        */
  real_T Fd_bp01Data[9];               /* Expression: [2430
                                          3385
                                          4780
                                          5950
                                          7650
                                          10810
                                          13260
                                          15350
                                          17190]
                                        * Referenced by: '<S283>/Fd'
                                        */
  real_T Constant2_Value_o;            /* Expression: 1
                                        * Referenced by: '<S295>/Constant2'
                                        */
  real_T TargetSlipRate_Value;         /* Expression: 0.18
                                        * Referenced by: '<S283>/Target SlipRate'
                                        */
  real_T Delay_InitialCondition;       /* Expression: 0.0
                                        * Referenced by: '<S296>/Delay'
                                        */
  real_T Gain_Gain_n3;                 /* Expression: 2
                                        * Referenced by: '<S298>/Gain'
                                        */
  real_T Pgain_Gain;                   /* Expression: 0.01
                                        * Referenced by: '<S294>/P gain'
                                        */
  real_T DiscreteTimeIntegrator_gainval;
                           /* Computed Parameter: DiscreteTimeIntegrator_gainval
                            * Referenced by: '<S294>/Discrete-Time Integrator'
                            */
  real_T DiscreteTimeIntegrator_IC;    /* Expression: 0
                                        * Referenced by: '<S294>/Discrete-Time Integrator'
                                        */
  real_T Constant_Value_fi;            /* Expression: 0
                                        * Referenced by: '<S294>/Constant'
                                        */
  real_T TCOn_Value;                   /* Expression: 0
                                        * Referenced by: '<Root>/TCOn'
                                        */
  real_T Switch_Threshold_g;           /* Expression: 0.5
                                        * Referenced by: '<S294>/Switch'
                                        */
  real_T toTire_rpm_Gain;              /* Expression: 60/(2*pi*0.228)
                                        * Referenced by: '<S284>/toTire_rpm'
                                        */
  real_T Gear_Gain;                    /* Expression: 13.75
                                        * Referenced by: '<S284>/Gear'
                                        */
  real_T Constant2_Value_l;            /* Expression: 1
                                        * Referenced by: '<S302>/Constant2'
                                        */
  real_T Frontsliprate_inBrake_ref_Value;/* Expression: 0.18
                                          * Referenced by: '<S284>/Front sliprate_inBrake_ref'
                                          */
  real_T FrontGain_Gain;               /* Expression: 1.2
                                        * Referenced by: '<S284>/FrontGain'
                                        */
  real_T Constant_Value_l;             /* Expression: 0
                                        * Referenced by: '<S284>/Constant'
                                        */
  real_T Gain1_Gain_dw;                /* Expression: -1
                                        * Referenced by: '<S294>/Gain1'
                                        */
  real_T Switch1_Threshold_f;          /* Expression: 0.5
                                        * Referenced by: '<S294>/Switch1'
                                        */
  real_T Gain_Gain_f;                  /* Expression: 2
                                        * Referenced by: '<S297>/Gain'
                                        */
  real_T Pgain_Gain_p;                 /* Expression: 0.01
                                        * Referenced by: '<S291>/P gain'
                                        */
  real_T DiscreteTimeIntegrator_gainva_j;
                          /* Computed Parameter: DiscreteTimeIntegrator_gainva_j
                           * Referenced by: '<S291>/Discrete-Time Integrator'
                           */
  real_T DiscreteTimeIntegrator_IC_d;  /* Expression: 0
                                        * Referenced by: '<S291>/Discrete-Time Integrator'
                                        */
  real_T Constant_Value_bj;            /* Expression: 0
                                        * Referenced by: '<S291>/Constant'
                                        */
  real_T Switch_Threshold_f;           /* Expression: 0.5
                                        * Referenced by: '<S291>/Switch'
                                        */
  real_T toTire_rpm_Gain_l;            /* Expression: 60/(2*pi*0.228)
                                        * Referenced by: '<S285>/toTire_rpm'
                                        */
  real_T Gear_Gain_m;                  /* Expression: 13.75
                                        * Referenced by: '<S285>/Gear'
                                        */
  real_T Constant2_Value_g;            /* Expression: 1
                                        * Referenced by: '<S304>/Constant2'
                                        */
  real_T Rearsliprate_inBrake_ref_Value;/* Expression: 0.18
                                         * Referenced by: '<S285>/Rear sliprate_inBrake_ref'
                                         */
  real_T Constant_Value_m0;            /* Expression: 0
                                        * Referenced by: '<S285>/Constant'
                                        */
  real_T Gain1_Gain_j;                 /* Expression: -1
                                        * Referenced by: '<S291>/Gain1'
                                        */
  real_T Switch1_Threshold_p2;         /* Expression: 0.5
                                        * Referenced by: '<S291>/Switch1'
                                        */
  real_T Saturation_UpperSat_l;        /* Expression: 21
                                        * Referenced by: '<S186>/Saturation'
                                        */
  real_T Saturation_LowerSat_i;        /* Expression: 0
                                        * Referenced by: '<S186>/Saturation'
                                        */
  real_T Saturation1_UpperSat_d5;      /* Expression: 19000
                                        * Referenced by: '<S186>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_k;       /* Expression: 0
                                        * Referenced by: '<S186>/Saturation1'
                                        */
  real_T Gain_Gain_e;                  /* Expression: 2*pi/60
                                        * Referenced by: '<S186>/Gain'
                                        */
  real_T uDLookupTable_tableData[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S186>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S186>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S186>/2-D Lookup Table'
   */
  real_T Pgain_Gain_h;                 /* Expression: 0.01
                                        * Referenced by: '<S292>/P gain'
                                        */
  real_T DiscreteTimeIntegrator_gainva_m;
                          /* Computed Parameter: DiscreteTimeIntegrator_gainva_m
                           * Referenced by: '<S292>/Discrete-Time Integrator'
                           */
  real_T DiscreteTimeIntegrator_IC_a;  /* Expression: 0
                                        * Referenced by: '<S292>/Discrete-Time Integrator'
                                        */
  real_T Constant_Value_lq;            /* Expression: 0
                                        * Referenced by: '<S292>/Constant'
                                        */
  real_T Switch_Threshold_o;           /* Expression: 0.5
                                        * Referenced by: '<S292>/Switch'
                                        */
  real_T Gain1_Gain_c;                 /* Expression: -1
                                        * Referenced by: '<S292>/Gain1'
                                        */
  real_T Switch1_Threshold_d;          /* Expression: 0.5
                                        * Referenced by: '<S292>/Switch1'
                                        */
  real_T Saturation_UpperSat_k;        /* Expression: 21
                                        * Referenced by: '<S185>/Saturation'
                                        */
  real_T Saturation_LowerSat_b;        /* Expression: 0
                                        * Referenced by: '<S185>/Saturation'
                                        */
  real_T Saturation1_UpperSat_g;       /* Expression: 19000
                                        * Referenced by: '<S185>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_o;       /* Expression: 0
                                        * Referenced by: '<S185>/Saturation1'
                                        */
  real_T Gain_Gain_g;                  /* Expression: 2*pi/60
                                        * Referenced by: '<S185>/Gain'
                                        */
  real_T uDLookupTable_tableData_o[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S185>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_i[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S185>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_a[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S185>/2-D Lookup Table'
   */
  real_T Pgain_Gain_o;                 /* Expression: 0.01
                                        * Referenced by: '<S293>/P gain'
                                        */
  real_T DiscreteTimeIntegrator_gainv_ji;
                          /* Computed Parameter: DiscreteTimeIntegrator_gainv_ji
                           * Referenced by: '<S293>/Discrete-Time Integrator'
                           */
  real_T DiscreteTimeIntegrator_IC_o;  /* Expression: 0
                                        * Referenced by: '<S293>/Discrete-Time Integrator'
                                        */
  real_T Constant_Value_fe;            /* Expression: 0
                                        * Referenced by: '<S293>/Constant'
                                        */
  real_T Switch_Threshold_i;           /* Expression: 0.5
                                        * Referenced by: '<S293>/Switch'
                                        */
  real_T Gain1_Gain_d4;                /* Expression: -1
                                        * Referenced by: '<S293>/Gain1'
                                        */
  real_T Switch1_Threshold_fn;         /* Expression: 0.5
                                        * Referenced by: '<S293>/Switch1'
                                        */
  real_T Saturation_UpperSat_j;        /* Expression: 21
                                        * Referenced by: '<S188>/Saturation'
                                        */
  real_T Saturation_LowerSat_l;        /* Expression: 0
                                        * Referenced by: '<S188>/Saturation'
                                        */
  real_T Saturation1_UpperSat_c;       /* Expression: 19000
                                        * Referenced by: '<S188>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_ov;      /* Expression: 0
                                        * Referenced by: '<S188>/Saturation1'
                                        */
  real_T Gain_Gain_p;                  /* Expression: 2*pi/60
                                        * Referenced by: '<S188>/Gain'
                                        */
  real_T uDLookupTable_tableData_h[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S188>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_l[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S188>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_e[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S188>/2-D Lookup Table'
   */
  real_T Saturation_UpperSat_o;        /* Expression: 21
                                        * Referenced by: '<S187>/Saturation'
                                        */
  real_T Saturation_LowerSat_cm;       /* Expression: 0
                                        * Referenced by: '<S187>/Saturation'
                                        */
  real_T Saturation1_UpperSat_m2;      /* Expression: 19000
                                        * Referenced by: '<S187>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_e;       /* Expression: 0
                                        * Referenced by: '<S187>/Saturation1'
                                        */
  real_T Gain_Gain_hv;                 /* Expression: 2*pi/60
                                        * Referenced by: '<S187>/Gain'
                                        */
  real_T uDLookupTable_tableData_e[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S187>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_iu[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S187>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_l[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S187>/2-D Lookup Table'
   */
  real_T Saturation_UpperSat_g;        /* Expression: inf
                                        * Referenced by: '<S178>/Saturation'
                                        */
  real_T Saturation_LowerSat_e;        /* Expression: 1
                                        * Referenced by: '<S178>/Saturation'
                                        */
  real_T MaxPower_Value;               /* Expression: 78000
                                        * Referenced by: '<S170>/Max Power'
                                        */
  real_T Eff_INV_Gain;                 /* Expression: 98/100
                                        * Referenced by: '<S170>/Eff_INV'
                                        */
  real_T Saturation4_UpperSat;         /* Expression: 21
                                        * Referenced by: '<S178>/Saturation4'
                                        */
  real_T Saturation4_LowerSat;         /* Expression: -21
                                        * Referenced by: '<S178>/Saturation4'
                                        */
  real_T Saturation1_UpperSat_p;       /* Expression: 21
                                        * Referenced by: '<S178>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_gh;      /* Expression: -21
                                        * Referenced by: '<S178>/Saturation1'
                                        */
  real_T Saturation_UpperSat_eq;       /* Expression: 21
                                        * Referenced by: '<S190>/Saturation'
                                        */
  real_T Saturation_LowerSat_m;        /* Expression: 0
                                        * Referenced by: '<S190>/Saturation'
                                        */
  real_T Saturation1_UpperSat_b;       /* Expression: 19000
                                        * Referenced by: '<S190>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_gn;      /* Expression: 0
                                        * Referenced by: '<S190>/Saturation1'
                                        */
  real_T Gain_Gain_bm;                 /* Expression: 2*pi/60
                                        * Referenced by: '<S190>/Gain'
                                        */
  real_T uDLookupTable_tableData_hp[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S190>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_id[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S190>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_n[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S190>/2-D Lookup Table'
   */
  real_T Saturation2_UpperSat_o;       /* Expression: 21
                                        * Referenced by: '<S178>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_h;       /* Expression: -21
                                        * Referenced by: '<S178>/Saturation2'
                                        */
  real_T Saturation_UpperSat_h;        /* Expression: 21
                                        * Referenced by: '<S189>/Saturation'
                                        */
  real_T Saturation_LowerSat_ir;       /* Expression: 0
                                        * Referenced by: '<S189>/Saturation'
                                        */
  real_T Saturation1_UpperSat_l;       /* Expression: 19000
                                        * Referenced by: '<S189>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_d;       /* Expression: 0
                                        * Referenced by: '<S189>/Saturation1'
                                        */
  real_T Gain_Gain_gh;                 /* Expression: 2*pi/60
                                        * Referenced by: '<S189>/Gain'
                                        */
  real_T uDLookupTable_tableData_d[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S189>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_h[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S189>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_j[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S189>/2-D Lookup Table'
   */
  real_T Saturation3_UpperSat_k;       /* Expression: 21
                                        * Referenced by: '<S178>/Saturation3'
                                        */
  real_T Saturation3_LowerSat_c;       /* Expression: -21
                                        * Referenced by: '<S178>/Saturation3'
                                        */
  real_T Saturation_UpperSat_e2;       /* Expression: 21
                                        * Referenced by: '<S192>/Saturation'
                                        */
  real_T Saturation_LowerSat_a;        /* Expression: 0
                                        * Referenced by: '<S192>/Saturation'
                                        */
  real_T Saturation1_UpperSat_bu;      /* Expression: 19000
                                        * Referenced by: '<S192>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_h;       /* Expression: 0
                                        * Referenced by: '<S192>/Saturation1'
                                        */
  real_T Gain_Gain_m;                  /* Expression: 2*pi/60
                                        * Referenced by: '<S192>/Gain'
                                        */
  real_T uDLookupTable_tableData_g[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S192>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_lj[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S192>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_p[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S192>/2-D Lookup Table'
   */
  real_T Saturation_UpperSat_b;        /* Expression: 21
                                        * Referenced by: '<S191>/Saturation'
                                        */
  real_T Saturation_LowerSat_j;        /* Expression: 0
                                        * Referenced by: '<S191>/Saturation'
                                        */
  real_T Saturation1_UpperSat_o;       /* Expression: 19000
                                        * Referenced by: '<S191>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_o0;      /* Expression: 0
                                        * Referenced by: '<S191>/Saturation1'
                                        */
  real_T Gain_Gain_a1;                 /* Expression: 2*pi/60
                                        * Referenced by: '<S191>/Gain'
                                        */
  real_T uDLookupTable_tableData_j[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S191>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_d[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S191>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_b[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S191>/2-D Lookup Table'
   */
  real_T Saturation_UpperSat_fe;       /* Expression: inf
                                        * Referenced by: '<S179>/Saturation'
                                        */
  real_T Saturation_LowerSat_px;       /* Expression: 1
                                        * Referenced by: '<S179>/Saturation'
                                        */
  real_T Saturation4_UpperSat_j;       /* Expression: 21
                                        * Referenced by: '<S179>/Saturation4'
                                        */
  real_T Saturation4_LowerSat_m;       /* Expression: -21
                                        * Referenced by: '<S179>/Saturation4'
                                        */
  real_T Saturation1_UpperSat_h;       /* Expression: 21
                                        * Referenced by: '<S179>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_d1;      /* Expression: -21
                                        * Referenced by: '<S179>/Saturation1'
                                        */
  real_T Saturation_UpperSat_fd;       /* Expression: 21
                                        * Referenced by: '<S194>/Saturation'
                                        */
  real_T Saturation_LowerSat_d;        /* Expression: 0
                                        * Referenced by: '<S194>/Saturation'
                                        */
  real_T Saturation1_UpperSat_e;       /* Expression: 19000
                                        * Referenced by: '<S194>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_n;       /* Expression: 0
                                        * Referenced by: '<S194>/Saturation1'
                                        */
  real_T Gain_Gain_br;                 /* Expression: 2*pi/60
                                        * Referenced by: '<S194>/Gain'
                                        */
  real_T uDLookupTable_tableData_p[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S194>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_i4[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S194>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_f[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S194>/2-D Lookup Table'
   */
  real_T Saturation2_UpperSat_p;       /* Expression: 21
                                        * Referenced by: '<S179>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_i;       /* Expression: -21
                                        * Referenced by: '<S179>/Saturation2'
                                        */
  real_T Saturation_UpperSat_fg;       /* Expression: 21
                                        * Referenced by: '<S193>/Saturation'
                                        */
  real_T Saturation_LowerSat_ns;       /* Expression: 0
                                        * Referenced by: '<S193>/Saturation'
                                        */
  real_T Saturation1_UpperSat_de;      /* Expression: 19000
                                        * Referenced by: '<S193>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_g4;      /* Expression: 0
                                        * Referenced by: '<S193>/Saturation1'
                                        */
  real_T Gain_Gain_fy;                 /* Expression: 2*pi/60
                                        * Referenced by: '<S193>/Gain'
                                        */
  real_T uDLookupTable_tableData_l[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S193>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_b[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S193>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_be[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S193>/2-D Lookup Table'
   */
  real_T Saturation3_UpperSat_k0;      /* Expression: 21
                                        * Referenced by: '<S179>/Saturation3'
                                        */
  real_T Saturation3_LowerSat_p;       /* Expression: -21
                                        * Referenced by: '<S179>/Saturation3'
                                        */
  real_T Saturation_UpperSat_m;        /* Expression: 21
                                        * Referenced by: '<S196>/Saturation'
                                        */
  real_T Saturation_LowerSat_o;        /* Expression: 0
                                        * Referenced by: '<S196>/Saturation'
                                        */
  real_T Saturation1_UpperSat_k;       /* Expression: 19000
                                        * Referenced by: '<S196>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_nu;      /* Expression: 0
                                        * Referenced by: '<S196>/Saturation1'
                                        */
  real_T Gain_Gain_c;                  /* Expression: 2*pi/60
                                        * Referenced by: '<S196>/Gain'
                                        */
  real_T uDLookupTable_tableData_ji[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S196>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_m[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S196>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_jg[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S196>/2-D Lookup Table'
   */
  real_T Saturation_UpperSat_j2;       /* Expression: 21
                                        * Referenced by: '<S195>/Saturation'
                                        */
  real_T Saturation_LowerSat_lr;       /* Expression: 0
                                        * Referenced by: '<S195>/Saturation'
                                        */
  real_T Saturation1_UpperSat_hj;      /* Expression: 19000
                                        * Referenced by: '<S195>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_j;       /* Expression: 0
                                        * Referenced by: '<S195>/Saturation1'
                                        */
  real_T Gain_Gain_lb;                 /* Expression: 2*pi/60
                                        * Referenced by: '<S195>/Gain'
                                        */
  real_T uDLookupTable_tableData_ek[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S195>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_j[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S195>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_c[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S195>/2-D Lookup Table'
   */
  real_T Saturation_UpperSat_ow;       /* Expression: inf
                                        * Referenced by: '<S180>/Saturation'
                                        */
  real_T Saturation_LowerSat_au;       /* Expression: 1
                                        * Referenced by: '<S180>/Saturation'
                                        */
  real_T Saturation4_UpperSat_f;       /* Expression: 21
                                        * Referenced by: '<S180>/Saturation4'
                                        */
  real_T Saturation4_LowerSat_b;       /* Expression: -21
                                        * Referenced by: '<S180>/Saturation4'
                                        */
  real_T Saturation1_UpperSat_n;       /* Expression: 21
                                        * Referenced by: '<S180>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_ob;      /* Expression: -21
                                        * Referenced by: '<S180>/Saturation1'
                                        */
  real_T Saturation_UpperSat_om;       /* Expression: 21
                                        * Referenced by: '<S198>/Saturation'
                                        */
  real_T Saturation_LowerSat_lg;       /* Expression: 0
                                        * Referenced by: '<S198>/Saturation'
                                        */
  real_T Saturation1_UpperSat_gq;      /* Expression: 19000
                                        * Referenced by: '<S198>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_ey;      /* Expression: 0
                                        * Referenced by: '<S198>/Saturation1'
                                        */
  real_T Gain_Gain_f0;                 /* Expression: 2*pi/60
                                        * Referenced by: '<S198>/Gain'
                                        */
  real_T uDLookupTable_tableData_ld[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S198>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_k[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S198>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_d[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S198>/2-D Lookup Table'
   */
  real_T Saturation2_UpperSat_c;       /* Expression: 21
                                        * Referenced by: '<S180>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_d;       /* Expression: -21
                                        * Referenced by: '<S180>/Saturation2'
                                        */
  real_T Saturation_UpperSat_a2;       /* Expression: 21
                                        * Referenced by: '<S197>/Saturation'
                                        */
  real_T Saturation_LowerSat_cl;       /* Expression: 0
                                        * Referenced by: '<S197>/Saturation'
                                        */
  real_T Saturation1_UpperSat_bup;     /* Expression: 19000
                                        * Referenced by: '<S197>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_b;       /* Expression: 0
                                        * Referenced by: '<S197>/Saturation1'
                                        */
  real_T Gain_Gain_ik;                 /* Expression: 2*pi/60
                                        * Referenced by: '<S197>/Gain'
                                        */
  real_T uDLookupTable_tableData_of[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S197>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_bz[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S197>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_fw[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S197>/2-D Lookup Table'
   */
  real_T Saturation3_UpperSat_e;       /* Expression: 21
                                        * Referenced by: '<S180>/Saturation3'
                                        */
  real_T Saturation3_LowerSat_l;       /* Expression: -21
                                        * Referenced by: '<S180>/Saturation3'
                                        */
  real_T Saturation_UpperSat_nt;       /* Expression: 21
                                        * Referenced by: '<S200>/Saturation'
                                        */
  real_T Saturation_LowerSat_hb;       /* Expression: 0
                                        * Referenced by: '<S200>/Saturation'
                                        */
  real_T Saturation1_UpperSat_a;       /* Expression: 19000
                                        * Referenced by: '<S200>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_nn;      /* Expression: 0
                                        * Referenced by: '<S200>/Saturation1'
                                        */
  real_T Gain_Gain_lm;                 /* Expression: 2*pi/60
                                        * Referenced by: '<S200>/Gain'
                                        */
  real_T uDLookupTable_tableData_e5[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S200>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_f[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S200>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_cb[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S200>/2-D Lookup Table'
   */
  real_T Saturation_UpperSat_i;        /* Expression: 21
                                        * Referenced by: '<S199>/Saturation'
                                        */
  real_T Saturation_LowerSat_jy;       /* Expression: 0
                                        * Referenced by: '<S199>/Saturation'
                                        */
  real_T Saturation1_UpperSat_ex;      /* Expression: 19000
                                        * Referenced by: '<S199>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_m;       /* Expression: 0
                                        * Referenced by: '<S199>/Saturation1'
                                        */
  real_T Gain_Gain_fl;                 /* Expression: 2*pi/60
                                        * Referenced by: '<S199>/Gain'
                                        */
  real_T uDLookupTable_tableData_ot[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S199>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_n[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S199>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_n1[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S199>/2-D Lookup Table'
   */
  real_T Saturation_UpperSat_d;        /* Expression: inf
                                        * Referenced by: '<S181>/Saturation'
                                        */
  real_T Saturation_LowerSat_pj;       /* Expression: 1
                                        * Referenced by: '<S181>/Saturation'
                                        */
  real_T Saturation4_UpperSat_b;       /* Expression: 21
                                        * Referenced by: '<S181>/Saturation4'
                                        */
  real_T Saturation4_LowerSat_l;       /* Expression: -21
                                        * Referenced by: '<S181>/Saturation4'
                                        */
  real_T Saturation1_UpperSat_j;       /* Expression: 21
                                        * Referenced by: '<S181>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_ng;      /* Expression: -21
                                        * Referenced by: '<S181>/Saturation1'
                                        */
  real_T Saturation_UpperSat_du;       /* Expression: 21
                                        * Referenced by: '<S202>/Saturation'
                                        */
  real_T Saturation_LowerSat_jd;       /* Expression: 0
                                        * Referenced by: '<S202>/Saturation'
                                        */
  real_T Saturation1_UpperSat_ij;      /* Expression: 19000
                                        * Referenced by: '<S202>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_a;       /* Expression: 0
                                        * Referenced by: '<S202>/Saturation1'
                                        */
  real_T Gain_Gain_d;                  /* Expression: 2*pi/60
                                        * Referenced by: '<S202>/Gain'
                                        */
  real_T uDLookupTable_tableData_n[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S202>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_o[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S202>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_i[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S202>/2-D Lookup Table'
   */
  real_T Saturation2_UpperSat_j;       /* Expression: 21
                                        * Referenced by: '<S181>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_f;       /* Expression: -21
                                        * Referenced by: '<S181>/Saturation2'
                                        */
  real_T Saturation_UpperSat_aq;       /* Expression: 21
                                        * Referenced by: '<S201>/Saturation'
                                        */
  real_T Saturation_LowerSat_mb;       /* Expression: 0
                                        * Referenced by: '<S201>/Saturation'
                                        */
  real_T Saturation1_UpperSat_ee;      /* Expression: 19000
                                        * Referenced by: '<S201>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_l;       /* Expression: 0
                                        * Referenced by: '<S201>/Saturation1'
                                        */
  real_T Gain_Gain_lg;                 /* Expression: 2*pi/60
                                        * Referenced by: '<S201>/Gain'
                                        */
  real_T uDLookupTable_tableData_gk[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S201>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_c[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S201>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_h[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S201>/2-D Lookup Table'
   */
  real_T Saturation3_UpperSat_ma;      /* Expression: 21
                                        * Referenced by: '<S181>/Saturation3'
                                        */
  real_T Saturation3_LowerSat_ds;      /* Expression: -21
                                        * Referenced by: '<S181>/Saturation3'
                                        */
  real_T Saturation_UpperSat_ntr;      /* Expression: 21
                                        * Referenced by: '<S204>/Saturation'
                                        */
  real_T Saturation_LowerSat_ib;       /* Expression: 0
                                        * Referenced by: '<S204>/Saturation'
                                        */
  real_T Saturation1_UpperSat_ae;      /* Expression: 19000
                                        * Referenced by: '<S204>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_kj;      /* Expression: 0
                                        * Referenced by: '<S204>/Saturation1'
                                        */
  real_T Gain_Gain_o;                  /* Expression: 2*pi/60
                                        * Referenced by: '<S204>/Gain'
                                        */
  real_T uDLookupTable_tableData_dj[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S204>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_mg[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S204>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_cu[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S204>/2-D Lookup Table'
   */
  real_T Saturation_UpperSat_aw;       /* Expression: 21
                                        * Referenced by: '<S203>/Saturation'
                                        */
  real_T Saturation_LowerSat_kd;       /* Expression: 0
                                        * Referenced by: '<S203>/Saturation'
                                        */
  real_T Saturation1_UpperSat_e1;      /* Expression: 19000
                                        * Referenced by: '<S203>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_ii;      /* Expression: 0
                                        * Referenced by: '<S203>/Saturation1'
                                        */
  real_T Gain_Gain_bj;                 /* Expression: 2*pi/60
                                        * Referenced by: '<S203>/Gain'
                                        */
  real_T uDLookupTable_tableData_hn[110];
  /* Expression: [ 0.6437    0.7133    0.7364    0.7470    0.7543    0.7657    0.7700    0.7708    0.7756    0.7814
     0.5842    0.7048    0.7757    0.8040    0.8201    0.8392    0.8516    0.8544    0.8597    0.8650
     0.4494    0.6081    0.7335    0.7882    0.8194    0.8543    0.8820    0.8888    0.8971    0.9044
     0.3559    0.5190    0.6702    0.7426    0.7854    0.8342    0.8758    0.8865    0.8984    0.9086
     0.2914    0.4478    0.6101    0.6941    0.7457    0.8062    0.8593    0.8758    0.8886    0.9016
     0.2417    0.3871    0.5522    0.6439    0.7024    0.7730    0.8373    0.8548    0.8737    0.8898
     0.2041    0.3371    0.5004    0.5965    0.6599    0.7388    0.8133    0.8342    0.8566    0.8759
     0.1731    0.2940    0.4510    0.5487    0.6155    0.7010    0.7856    0.8097    0.8356    0.8581
     0.1482    0.2575    0.4067    0.5041    0.5728    0.6634    0.7570    0.7840    0.8134    0.8271
     0.1281    0.2267    0.3672    0.4630    0.5325    0.6267    0.7277    0.7575    0.7902    0.7696
     0.1117    0.2005    0.3321    0.4251    0.4944    0.5909    0.6982    0.7306    0.6766    0.6928]
   * Referenced by: '<S203>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp01Data_d0[11];
  /* Expression: [1.3000    2.7000    5.4000    7.9000   10.4000   12.5000   14.4000   16.0000   17.4000   18.5000   19.6000]
   * Referenced by: '<S203>/2-D Lookup Table'
   */
  real_T uDLookupTable_bp02Data_m[10];
  /* Expression: [ 500        1000        2000        3000        4000        6000       10000       12000       15000       19000]
   * Referenced by: '<S203>/2-D Lookup Table'
   */
  real_T Saturation_UpperSat_lo;       /* Expression: inf
                                        * Referenced by: '<S182>/Saturation'
                                        */
  real_T Saturation_LowerSat_o1;       /* Expression: 1
                                        * Referenced by: '<S182>/Saturation'
                                        */
  real_T Saturation4_UpperSat_o;       /* Expression: 21
                                        * Referenced by: '<S182>/Saturation4'
                                        */
  real_T Saturation4_LowerSat_p;       /* Expression: -21
                                        * Referenced by: '<S182>/Saturation4'
                                        */
  real_T Constant_Value_lo;            /* Expression: 80000
                                        * Referenced by: '<S183>/Constant'
                                        */
  real_T Gain_Gain_fw;                 /* Expression: 1/500
                                        * Referenced by: '<S183>/Gain'
                                        */
  real_T Switch2_Threshold_m;          /* Expression: 0
                                        * Referenced by: '<S183>/Switch2'
                                        */
  real_T Switch3_Threshold_l;          /* Expression: 0.5
                                        * Referenced by: '<S184>/Switch3'
                                        */
  real_T Gain3_Gain;                   /* Expression: 1
                                        * Referenced by: '<S122>/Gain3'
                                        */
  real_T Gain30_Gain;                  /* Expression: 1000
                                        * Referenced by: '<S37>/Gain30'
                                        */
  real_T Saturation3_UpperSat_kj;      /* Expression: 21
                                        * Referenced by: '<S182>/Saturation3'
                                        */
  real_T Saturation3_LowerSat_b;       /* Expression: -21
                                        * Referenced by: '<S182>/Saturation3'
                                        */
  real_T Switch3_Threshold_le;         /* Expression: 0
                                        * Referenced by: '<S183>/Switch3'
                                        */
  real_T Switch2_Threshold_g;          /* Expression: 0.5
                                        * Referenced by: '<S184>/Switch2'
                                        */
  real_T Gain2_Gain_o;                 /* Expression: 1
                                        * Referenced by: '<S122>/Gain2'
                                        */
  real_T Gain31_Gain;                  /* Expression: 1000
                                        * Referenced by: '<S37>/Gain31'
                                        */
  real_T Saturation1_UpperSat_f;       /* Expression: 21
                                        * Referenced by: '<S182>/Saturation1'
                                        */
  real_T Saturation1_LowerSat_c;       /* Expression: -21
                                        * Referenced by: '<S182>/Saturation1'
                                        */
  real_T Switch_Threshold_in;          /* Expression: 0
                                        * Referenced by: '<S183>/Switch'
                                        */
  real_T Switch_Threshold_h;           /* Expression: 0.5
                                        * Referenced by: '<S184>/Switch'
                                        */
  real_T Gain_Gain_oh;                 /* Expression: 1
                                        * Referenced by: '<S122>/Gain'
                                        */
  real_T Gain32_Gain;                  /* Expression: 1000
                                        * Referenced by: '<S37>/Gain32'
                                        */
  real_T Saturation2_UpperSat_jl;      /* Expression: 21
                                        * Referenced by: '<S182>/Saturation2'
                                        */
  real_T Saturation2_LowerSat_ir;      /* Expression: -21
                                        * Referenced by: '<S182>/Saturation2'
                                        */
  real_T Switch1_Threshold_d1;         /* Expression: 0
                                        * Referenced by: '<S183>/Switch1'
                                        */
  real_T Switch1_Threshold_c;          /* Expression: 0.5
                                        * Referenced by: '<S184>/Switch1'
                                        */
  real_T Gain1_Gain_e;                 /* Expression: 1
                                        * Referenced by: '<S122>/Gain1'
                                        */
  real_T Gain33_Gain;                  /* Expression: 1000
                                        * Referenced by: '<S37>/Gain33'
                                        */
  real_T Gain1_Gain_b;                 /* Expression: 10^3
                                        * Referenced by: '<S37>/Gain1'
                                        */
  real_T TransferFcn_A_e;              /* Computed Parameter: TransferFcn_A_e
                                        * Referenced by: '<S220>/Transfer Fcn'
                                        */
  real_T TransferFcn_C_l;              /* Computed Parameter: TransferFcn_C_l
                                        * Referenced by: '<S220>/Transfer Fcn'
                                        */
  real_T Constant3_Value;              /* Expression: 1
                                        * Referenced by: '<S222>/Constant3'
                                        */
  real_T StabilityFactor_ref_Value;    /* Expression: 0.00125
                                        * Referenced by: '<S216>/StabilityFactor_ref'
                                        */
  real_T Gain6_Gain;                   /* Expression: 1/1.53
                                        * Referenced by: '<S222>/Gain6'
                                        */
  real_T Gain2_Gain_c;                 /* Expression: 10^3
                                        * Referenced by: '<S37>/Gain2'
                                        */
  real_T Gain_Gain_k;                  /* Expression: 0.001
                                        * Referenced by: '<S176>/Gain'
                                        */
  real_T Gain3_Gain_c;                 /* Expression: 10^3
                                        * Referenced by: '<S37>/Gain3'
                                        */
  real_T toTire_rpm_Gain_g;            /* Expression: 60/(2*pi*0.228)
                                        * Referenced by: '<S322>/toTire_rpm'
                                        */
  real_T Gear_Gain_e;                  /* Expression: 13.75
                                        * Referenced by: '<S322>/Gear'
                                        */
  real_T Gain4_Gain_k;                 /* Expression: 1000
                                        * Referenced by: '<S37>/Gain4'
                                        */
  real_T toTire_rpm_Gain_b;            /* Expression: 60/(2*pi*0.228)
                                        * Referenced by: '<S323>/toTire_rpm'
                                        */
  real_T Gear_Gain_p;                  /* Expression: 13.75
                                        * Referenced by: '<S323>/Gear'
                                        */
  real_T Gain5_Gain_l;                 /* Expression: 1000
                                        * Referenced by: '<S37>/Gain5'
                                        */
  real_T toTire_rpm_Gain_m;            /* Expression: 60/(2*pi*0.228)
                                        * Referenced by: '<S324>/toTire_rpm'
                                        */
  real_T Gear_Gain_j;                  /* Expression: 13.75
                                        * Referenced by: '<S324>/Gear'
                                        */
  real_T Gain6_Gain_b;                 /* Expression: 1000
                                        * Referenced by: '<S37>/Gain6'
                                        */
  real_T toTire_rpm_Gain_h;            /* Expression: 60/(2*pi*0.228)
                                        * Referenced by: '<S325>/toTire_rpm'
                                        */
  real_T Gear_Gain_f;                  /* Expression: 13.75
                                        * Referenced by: '<S325>/Gear'
                                        */
  real_T Gain7_Gain;                   /* Expression: 1000
                                        * Referenced by: '<S37>/Gain7'
                                        */
  real_T VDC_On_Value;                 /* Expression: 1
                                        * Referenced by: '<Root>/VDC_On'
                                        */
  real_T Gain18_Gain;                  /* Expression: 1
                                        * Referenced by: '<S37>/Gain18'
                                        */
  real_T Gain20_Gain;                  /* Expression: 1
                                        * Referenced by: '<S37>/Gain20'
                                        */
  real_T Gain21_Gain;                  /* Expression: 1
                                        * Referenced by: '<S37>/Gain21'
                                        */
  real_T Gain22_Gain;                  /* Expression: 1
                                        * Referenced by: '<S37>/Gain22'
                                        */
  real_T Gain23_Gain;                  /* Expression: 1
                                        * Referenced by: '<S37>/Gain23'
                                        */
  real_T Gain24_Gain;                  /* Expression: 1
                                        * Referenced by: '<S37>/Gain24'
                                        */
  real_T Gain25_Gain;                  /* Expression: 1
                                        * Referenced by: '<S37>/Gain25'
                                        */
  real_T Gain19_Gain;                  /* Expression: 100
                                        * Referenced by: '<S36>/Gain19'
                                        */
  real_T Gain26_Gain;                  /* Expression: 100
                                        * Referenced by: '<S36>/Gain26'
                                        */
  real_T Gain27_Gain;                  /* Expression: 100
                                        * Referenced by: '<S36>/Gain27'
                                        */
  real_T Gain28_Gain;                  /* Expression: 100
                                        * Referenced by: '<S36>/Gain28'
                                        */
  real_T ACCEL_X_Gain_Gain;            /* Expression: 10^-2
                                        * Referenced by: '<S14>/ACCEL_X_Gain'
                                        */
  real_T Gain_Gain_dv;                 /* Expression: -1
                                        * Referenced by: '<S14>/Gain'
                                        */
  real_T Gain3_Gain_l;                 /* Expression: 10^2
                                        * Referenced by: '<S35>/Gain3'
                                        */
  real_T ACCEL_Y_Gain_Gain;            /* Expression: 10^-2
                                        * Referenced by: '<S14>/ACCEL_Y_Gain'
                                        */
  real_T Gain5_Gain_o;                 /* Expression: 10^2
                                        * Referenced by: '<S35>/Gain5'
                                        */
  real_T ACCEL_Z_Gain_Gain;            /* Expression: 10^-2
                                        * Referenced by: '<S14>/ACCEL_Z_Gain'
                                        */
  real_T Gain4_Gain_p;                 /* Expression: 10^2
                                        * Referenced by: '<S35>/Gain4'
                                        */
  real_T Gain3_Gain_p;                 /* Expression: 10^-3
                                        * Referenced by: '<S14>/Gain3'
                                        */
  real_T Gain6_Gain_p;                 /* Expression: 10^3
                                        * Referenced by: '<S35>/Gain6'
                                        */
  real_T Gain5_Gain_p;                 /* Expression: 10^-3
                                        * Referenced by: '<S14>/Gain5'
                                        */
  real_T Gain8_Gain_f;                 /* Expression: 10^3
                                        * Referenced by: '<S35>/Gain8'
                                        */
  real_T Gain4_Gain_ns;                /* Expression: 10^-3
                                        * Referenced by: '<S14>/Gain4'
                                        */
  real_T Gain7_Gain_m;                 /* Expression: 10^3
                                        * Referenced by: '<S35>/Gain7'
                                        */
  real_T Gain9_Gain;                   /* Expression: 10^2
                                        * Referenced by: '<S35>/Gain9'
                                        */
  real_T Gain10_Gain_o;                /* Expression: 10^-2
                                        * Referenced by: '<S14>/Gain10'
                                        */
  real_T Gain10_Gain_m;                /* Expression: 10^2
                                        * Referenced by: '<S35>/Gain10'
                                        */
  real_T Gain9_Gain_o;                 /* Expression: 10^-2
                                        * Referenced by: '<S14>/Gain9'
                                        */
  real_T Gain14_Gain_e;                /* Expression: 10^2
                                        * Referenced by: '<S35>/Gain14'
                                        */
  real_T Gain24_Gain_g;                /* Expression: 10^-4
                                        * Referenced by: '<S14>/Gain24'
                                        */
  real_T Gain11_Gain_e;                /* Expression: 10^4
                                        * Referenced by: '<S35>/Gain11'
                                        */
  real_T Gain23_Gain_p;                /* Expression: 10^-4
                                        * Referenced by: '<S14>/Gain23'
                                        */
  real_T Gain13_Gain_a;                /* Expression: 10^4
                                        * Referenced by: '<S35>/Gain13'
                                        */
  real_T Gain25_Gain_a;                /* Expression: 10^-2
                                        * Referenced by: '<S14>/Gain25'
                                        */
  real_T Gain12_Gain_m;                /* Expression: 10^2
                                        * Referenced by: '<S35>/Gain12'
                                        */
  real_T Gain21_Gain_n;                /* Expression: 10^-7
                                        * Referenced by: '<S14>/Gain21'
                                        */
  real_T Gain16_Gain;                  /* Expression: 10^7
                                        * Referenced by: '<S35>/Gain16'
                                        */
  real_T Gain22_Gain_d;                /* Expression: 10^-7
                                        * Referenced by: '<S14>/Gain22'
                                        */
  real_T Gain15_Gain;                  /* Expression: 10^7
                                        * Referenced by: '<S35>/Gain15'
                                        */
  real_T Gain18_Gain_n;                /* Expression: 1
                                        * Referenced by: '<S14>/Gain18'
                                        */
  real_T Gain2_Gain_h;                 /* Expression: 1
                                        * Referenced by: '<S35>/Gain2'
                                        */
  real_T Gain19_Gain_f;                /* Expression: 1
                                        * Referenced by: '<S14>/Gain19'
                                        */
  real_T Gain1_Gain_el;                /* Expression: 1
                                        * Referenced by: '<S35>/Gain1'
                                        */
  real_T Gain20_Gain_m;                /* Expression: 1
                                        * Referenced by: '<S14>/Gain20'
                                        */
  real_T Gain18_Gain_l;                /* Expression: 1
                                        * Referenced by: '<S35>/Gain18'
                                        */
  real_T Gain_Gain_ic;                 /* Expression: 1000
                                        * Referenced by: '<S34>/Gain'
                                        */
  real_T Gain29_Gain;                  /* Expression: 100
                                        * Referenced by: '<S34>/Gain29'
                                        */
  real_T Gain2_Gain_m;                 /* Expression: 1000
                                        * Referenced by: '<S34>/Gain2'
                                        */
  real_T Gain1_Gain_g;                 /* Expression: 1000
                                        * Referenced by: '<S34>/Gain1'
                                        */
  real_T Constant_Value_it;            /* Expression: 123
                                        * Referenced by: '<S38>/Constant'
                                        */
  real_T Constant_Value_b2;            /* Expression: -21.188
                                        * Referenced by: '<S339>/Constant'
                                        */
  real_T Constant1_Value_f;            /* Expression: 105.95
                                        * Referenced by: '<S339>/Constant1'
                                        */
  real_T Constant_Value_mx;            /* Expression: -21.188
                                        * Referenced by: '<S340>/Constant'
                                        */
  real_T Constant1_Value_l;            /* Expression: 105.95
                                        * Referenced by: '<S340>/Constant1'
                                        */
  real_T Constant_Value_n4;            /* Expression: 1
                                        * Referenced by: '<S133>/Constant'
                                        */
  real_T Switch1_Threshold_n;          /* Expression: 0
                                        * Referenced by: '<S133>/Switch1'
                                        */
  real_T Switch2_Threshold_h;          /* Expression: 0
                                        * Referenced by: '<S133>/Switch2'
                                        */
  real_T Switch3_Threshold_dy;         /* Expression: 0
                                        * Referenced by: '<S133>/Switch3'
                                        */
  real_T Switch4_Threshold_n;          /* Expression: 0
                                        * Referenced by: '<S133>/Switch4'
                                        */
  real_T Constant_Value_i4;            /* Expression: 0
                                        * Referenced by: '<S131>/Constant'
                                        */
  real_T Switch1_Threshold_dd;         /* Expression: 0
                                        * Referenced by: '<S131>/Switch1'
                                        */
  real_T Switch2_Threshold_m1;         /* Expression: 0
                                        * Referenced by: '<S131>/Switch2'
                                        */
  real_T Switch3_Threshold_o;          /* Expression: 0
                                        * Referenced by: '<S131>/Switch3'
                                        */
  real_T Switch4_Threshold_g;          /* Expression: 0
                                        * Referenced by: '<S131>/Switch4'
                                        */
  real_T Constant_Value_ou;            /* Expression: 0
                                        * Referenced by: '<S124>/Constant'
                                        */
  real_T Switch1_Threshold_pk;         /* Expression: 0
                                        * Referenced by: '<S124>/Switch1'
                                        */
  real_T Constant_Value_cc;            /* Expression: 1
                                        * Referenced by: '<S134>/Constant'
                                        */
  real_T Switch1_Threshold_nh;         /* Expression: 0
                                        * Referenced by: '<S134>/Switch1'
                                        */
  real_T Switch2_Threshold_gb;         /* Expression: 0
                                        * Referenced by: '<S134>/Switch2'
                                        */
  real_T Switch3_Threshold_j;          /* Expression: 0
                                        * Referenced by: '<S134>/Switch3'
                                        */
  real_T Switch4_Threshold_b;          /* Expression: 0
                                        * Referenced by: '<S134>/Switch4'
                                        */
  real_T Constant_Value_bn;            /* Expression: 0
                                        * Referenced by: '<S123>/Constant'
                                        */
  real_T Switch1_Threshold_g;          /* Expression: 0
                                        * Referenced by: '<S123>/Switch1'
                                        */
  real_T ErrorReset_ControllDesk_Value;/* Expression: 0
                                        * Referenced by: '<S11>/ErrorReset_ControllDesk'
                                        */
  real_T Gain8_Gain_g;                 /* Expression: 1
                                        * Referenced by: '<S122>/Gain8'
                                        */
  real_T Switch_Threshold_a;           /* Expression: 0.5
                                        * Referenced by: '<S136>/Switch'
                                        */
  real_T Saturation_UpperSat_a4;       /* Expression: 30000
                                        * Referenced by: '<S136>/Saturation'
                                        */
  real_T Saturation_LowerSat_eu;       /* Expression: 1
                                        * Referenced by: '<S136>/Saturation'
                                        */
  real_T Constant_Value_cj;            /* Expression: 0
                                        * Referenced by: '<S125>/Constant'
                                        */
  real_T MC1_Threshold;                /* Expression: 0
                                        * Referenced by: '<S125>/MC1'
                                        */
  real_T Constant1_Value_h;            /* Expression: 0
                                        * Referenced by: '<S136>/Constant1'
                                        */
  real_T Gain2_Gain_oe;                /* Expression: 1000/9.8
                                        * Referenced by: '<S117>/Gain2'
                                        */
  real_T Quantizer1_Interval;          /* Expression: 1
                                        * Referenced by: '<S117>/Quantizer1'
                                        */
  real_T Constant1_Value_a;            /* Expression: 0
                                        * Referenced by: '<S138>/Constant1'
                                        */
  real_T Switch_Threshold_j;           /* Expression: 0
                                        * Referenced by: '<S138>/Switch'
                                        */
  real_T Delay_InitialCondition_i;     /* Expression: 0
                                        * Referenced by: '<S138>/Delay'
                                        */
  real_T Delay_InitialCondition_p;     /* Expression: 0
                                        * Referenced by: '<S140>/Delay'
                                        */
  real_T MsetpointMAX01Mn_Value;       /* Expression: 2142
                                        * Referenced by: '<S137>/M set point MAX 0.1%Mn'
                                        */
  real_T powerlimit_Value;             /* Expression: 0
                                        * Referenced by: '<S135>/power limit'
                                        */
  real_T Constant_Value_gc;            /* Expression: -10000
                                        * Referenced by: '<S135>/Constant'
                                        */
  real_T Saturation_UpperSat_lu;       /* Expression: 30000
                                        * Referenced by: '<S135>/Saturation'
                                        */
  real_T Saturation_LowerSat_pt;       /* Expression: 1
                                        * Referenced by: '<S135>/Saturation'
                                        */
  real_T Gain6_Gain_k;                 /* Expression: 1000/9.8
                                        * Referenced by: '<S117>/Gain6'
                                        */
  real_T Quantizer2_Interval;          /* Expression: 1
                                        * Referenced by: '<S117>/Quantizer2'
                                        */
  real_T Delay_InitialCondition_l;     /* Expression: 0
                                        * Referenced by: '<S139>/Delay'
                                        */
  real_T Switch2_Threshold_f;          /* Expression: 0
                                        * Referenced by: '<S124>/Switch2'
                                        */
  real_T Switch2_Threshold_e;          /* Expression: 0
                                        * Referenced by: '<S123>/Switch2'
                                        */
  real_T Gain6_Gain_n;                 /* Expression: 1
                                        * Referenced by: '<S122>/Gain6'
                                        */
  real_T Switch_Threshold_lh;          /* Expression: 0.5
                                        * Referenced by: '<S144>/Switch'
                                        */
  real_T Saturation_UpperSat_lm;       /* Expression: 30000
                                        * Referenced by: '<S144>/Saturation'
                                        */
  real_T Saturation_LowerSat_h3;       /* Expression: 1
                                        * Referenced by: '<S144>/Saturation'
                                        */
  real_T MC2_Threshold;                /* Expression: 0
                                        * Referenced by: '<S125>/MC2'
                                        */
  real_T Constant1_Value_ie;           /* Expression: 0
                                        * Referenced by: '<S144>/Constant1'
                                        */
  real_T Gain2_Gain_l;                 /* Expression: 1000/9.8
                                        * Referenced by: '<S118>/Gain2'
                                        */
  real_T Quantizer1_Interval_d;        /* Expression: 1
                                        * Referenced by: '<S118>/Quantizer1'
                                        */
  real_T Constant1_Value_is;           /* Expression: 0
                                        * Referenced by: '<S146>/Constant1'
                                        */
  real_T Switch_Threshold_c;           /* Expression: 0
                                        * Referenced by: '<S146>/Switch'
                                        */
  real_T Delay_InitialCondition_n;     /* Expression: 0
                                        * Referenced by: '<S146>/Delay'
                                        */
  real_T Delay_InitialCondition_g;     /* Expression: 0
                                        * Referenced by: '<S148>/Delay'
                                        */
  real_T MsetpointMAX01Mn_Value_p;     /* Expression: 2142
                                        * Referenced by: '<S145>/M set point MAX 0.1%Mn'
                                        */
  real_T powerlimit_Value_a;           /* Expression: 0
                                        * Referenced by: '<S143>/power limit'
                                        */
  real_T Constant_Value_i0;            /* Expression: -10000
                                        * Referenced by: '<S143>/Constant'
                                        */
  real_T Saturation_UpperSat_eqy;      /* Expression: 30000
                                        * Referenced by: '<S143>/Saturation'
                                        */
  real_T Saturation_LowerSat_np;       /* Expression: 1
                                        * Referenced by: '<S143>/Saturation'
                                        */
  real_T Gain6_Gain_e;                 /* Expression: 1000/9.8
                                        * Referenced by: '<S118>/Gain6'
                                        */
  real_T Quantizer2_Interval_p;        /* Expression: 1
                                        * Referenced by: '<S118>/Quantizer2'
                                        */
  real_T Delay_InitialCondition_j;     /* Expression: 0
                                        * Referenced by: '<S147>/Delay'
                                        */
  real_T Switch3_Threshold_p;          /* Expression: 0
                                        * Referenced by: '<S124>/Switch3'
                                        */
  real_T Switch3_Threshold_f;          /* Expression: 0
                                        * Referenced by: '<S123>/Switch3'
                                        */
  real_T Gain7_Gain_h;                 /* Expression: 1
                                        * Referenced by: '<S122>/Gain7'
                                        */
  real_T Switch_Threshold_m;           /* Expression: 0.5
                                        * Referenced by: '<S152>/Switch'
                                        */
  real_T Saturation_UpperSat_lx;       /* Expression: 30000
                                        * Referenced by: '<S152>/Saturation'
                                        */
  real_T Saturation_LowerSat_bh;       /* Expression: 1
                                        * Referenced by: '<S152>/Saturation'
                                        */
  real_T MC3_Threshold;                /* Expression: 0
                                        * Referenced by: '<S125>/MC3'
                                        */
  real_T Constant1_Value_c;            /* Expression: 0
                                        * Referenced by: '<S152>/Constant1'
                                        */
  real_T Gain2_Gain_e;                 /* Expression: 1000/9.8
                                        * Referenced by: '<S119>/Gain2'
                                        */
  real_T Quantizer1_Interval_n;        /* Expression: 1
                                        * Referenced by: '<S119>/Quantizer1'
                                        */
  real_T Constant1_Value_n;            /* Expression: 0
                                        * Referenced by: '<S154>/Constant1'
                                        */
  real_T Switch_Threshold_fg;          /* Expression: 0
                                        * Referenced by: '<S154>/Switch'
                                        */
  real_T Delay_InitialCondition_m;     /* Expression: 0
                                        * Referenced by: '<S154>/Delay'
                                        */
  real_T Delay_InitialCondition_m0;    /* Expression: 0
                                        * Referenced by: '<S156>/Delay'
                                        */
  real_T MsetpointMAX01Mn_Value_d;     /* Expression: 2142
                                        * Referenced by: '<S153>/M set point MAX 0.1%Mn'
                                        */
  real_T powerlimit_Value_d;           /* Expression: 0
                                        * Referenced by: '<S151>/power limit'
                                        */
  real_T Constant_Value_nt;            /* Expression: -10000
                                        * Referenced by: '<S151>/Constant'
                                        */
  real_T Saturation_UpperSat_ep;       /* Expression: 30000
                                        * Referenced by: '<S151>/Saturation'
                                        */
  real_T Saturation_LowerSat_ms;       /* Expression: 1
                                        * Referenced by: '<S151>/Saturation'
                                        */
  real_T Gain6_Gain_pn;                /* Expression: 1000/9.8
                                        * Referenced by: '<S119>/Gain6'
                                        */
  real_T Quantizer2_Interval_k;        /* Expression: 1
                                        * Referenced by: '<S119>/Quantizer2'
                                        */
  real_T Delay_InitialCondition_i3;    /* Expression: 0
                                        * Referenced by: '<S155>/Delay'
                                        */
  real_T Switch4_Threshold_a;          /* Expression: 0
                                        * Referenced by: '<S124>/Switch4'
                                        */
  real_T Switch4_Threshold_l;          /* Expression: 0
                                        * Referenced by: '<S123>/Switch4'
                                        */
  real_T Gain9_Gain_e;                 /* Expression: 1
                                        * Referenced by: '<S122>/Gain9'
                                        */
  real_T Switch_Threshold_gk;          /* Expression: 0.5
                                        * Referenced by: '<S160>/Switch'
                                        */
  real_T Saturation_UpperSat_n3;       /* Expression: 30000
                                        * Referenced by: '<S160>/Saturation'
                                        */
  real_T Saturation_LowerSat_fq;       /* Expression: 1
                                        * Referenced by: '<S160>/Saturation'
                                        */
  real_T MC4_Threshold;                /* Expression: 0
                                        * Referenced by: '<S125>/MC4'
                                        */
  real_T Constant1_Value_d;            /* Expression: 0
                                        * Referenced by: '<S160>/Constant1'
                                        */
  real_T Gain2_Gain_d;                 /* Expression: 1000/9.8
                                        * Referenced by: '<S120>/Gain2'
                                        */
  real_T Quantizer1_Interval_f;        /* Expression: 1
                                        * Referenced by: '<S120>/Quantizer1'
                                        */
  real_T Constant1_Value_g;            /* Expression: 0
                                        * Referenced by: '<S162>/Constant1'
                                        */
  real_T Switch_Threshold_ik;          /* Expression: 0
                                        * Referenced by: '<S162>/Switch'
                                        */
  real_T Delay_InitialCondition_mz;    /* Expression: 0
                                        * Referenced by: '<S162>/Delay'
                                        */
  real_T Delay_InitialCondition_p2;    /* Expression: 0
                                        * Referenced by: '<S164>/Delay'
                                        */
  real_T MsetpointMAX01Mn_Value_a;     /* Expression: 2142
                                        * Referenced by: '<S161>/M set point MAX 0.1%Mn'
                                        */
  real_T powerlimit_Value_j;           /* Expression: 0
                                        * Referenced by: '<S159>/power limit'
                                        */
  real_T Constant_Value_mr;            /* Expression: -10000
                                        * Referenced by: '<S159>/Constant'
                                        */
  real_T Saturation_UpperSat_hu;       /* Expression: 30000
                                        * Referenced by: '<S159>/Saturation'
                                        */
  real_T Saturation_LowerSat_hg;       /* Expression: 1
                                        * Referenced by: '<S159>/Saturation'
                                        */
  real_T Gain6_Gain_l;                 /* Expression: 1000/9.8
                                        * Referenced by: '<S120>/Gain6'
                                        */
  real_T Quantizer2_Interval_po;       /* Expression: 1
                                        * Referenced by: '<S120>/Quantizer2'
                                        */
  real_T Delay_InitialCondition_a;     /* Expression: 0
                                        * Referenced by: '<S163>/Delay'
                                        */
  real_T LED1_Value;                   /* Expression: 1
                                        * Referenced by: '<S16>/LED1'
                                        */
  real_T LED2_Value;                   /* Expression: 0
                                        * Referenced by: '<S16>/LED2'
                                        */
  real_T LED3_Value;                   /* Expression: 0
                                        * Referenced by: '<S16>/LED3'
                                        */
  real_T LED5_Value;                   /* Expression: 0
                                        * Referenced by: '<S16>/LED5'
                                        */
  real_T buzzer_Value;                 /* Expression: 0
                                        * Referenced by: '<S16>/buzzer'
                                        */
  real_T LED_brightness_Value;         /* Expression: 84
                                        * Referenced by: '<S16>/LED_brightness'
                                        */
  real_T LED_Displaymode_Value;        /* Expression: 1
                                        * Referenced by: '<S16>/LED_Displaymode'
                                        */
  real_T LED_number2_Value;            /* Expression: 1
                                        * Referenced by: '<S16>/LED_number2'
                                        */
  real_T SWAS1_signal_plus_Value;      /* Expression: 0.148
                                        * Referenced by: '<S72>/SWAS1_signal_plus'
                                        */
  real_T SWAS1_signal_0deg_Value;      /* Expression: 0.279
                                        * Referenced by: '<S72>/SWAS1_signal_0deg'
                                        */
  real_T TransferFcn_A_a;              /* Computed Parameter: TransferFcn_A_a
                                        * Referenced by: '<S72>/Transfer Fcn'
                                        */
  real_T TransferFcn_C_a;              /* Computed Parameter: TransferFcn_C_a
                                        * Referenced by: '<S72>/Transfer Fcn'
                                        */
  real_T SWAS1_signal_minus_Value;     /* Expression: 0.403
                                        * Referenced by: '<S72>/SWAS1_signal_minus'
                                        */
  real_T Switch_Threshold_ob;          /* Expression: 0
                                        * Referenced by: '<S72>/Switch'
                                        */
  real_T Gain_Gain_ap;                 /* Expression: 1
                                        * Referenced by: '<S11>/Gain'
                                        */
  real_T Gain1_Gain_l;                 /* Expression: 1
                                        * Referenced by: '<S11>/Gain1'
                                        */
  real_T Delay10_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S89>/Delay10'
                                        */
  real_T Delay11_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S89>/Delay11'
                                        */
  real_T Delay12_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S89>/Delay12'
                                        */
  real_T Delay13_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S89>/Delay13'
                                        */
  real_T Delay14_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S89>/Delay14'
                                        */
  real_T Delay15_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S89>/Delay15'
                                        */
  real_T Delay16_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S89>/Delay16'
                                        */
  real_T Delay17_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S89>/Delay17'
                                        */
  real_T Delay18_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S89>/Delay18'
                                        */
  real_T Delay19_InitialCondition;     /* Expression: 0
                                        * Referenced by: '<S89>/Delay19'
                                        */
  real_T Delay_InitialCondition_gx;    /* Expression: 0
                                        * Referenced by: '<S89>/Delay'
                                        */
  real_T Delay1_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S89>/Delay1'
                                        */
  real_T Delay2_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S89>/Delay2'
                                        */
  real_T Delay3_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S89>/Delay3'
                                        */
  real_T Delay4_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S89>/Delay4'
                                        */
  real_T Delay5_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S89>/Delay5'
                                        */
  real_T Delay6_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S89>/Delay6'
                                        */
  real_T Delay7_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S89>/Delay7'
                                        */
  real_T Delay8_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S89>/Delay8'
                                        */
  real_T Delay9_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S89>/Delay9'
                                        */
  real_T Delay10_InitialCondition_l;   /* Expression: 0
                                        * Referenced by: '<S90>/Delay10'
                                        */
  real_T Delay11_InitialCondition_n;   /* Expression: 0
                                        * Referenced by: '<S90>/Delay11'
                                        */
  real_T Delay12_InitialCondition_d;   /* Expression: 0
                                        * Referenced by: '<S90>/Delay12'
                                        */
  real_T Delay13_InitialCondition_f;   /* Expression: 0
                                        * Referenced by: '<S90>/Delay13'
                                        */
  real_T Delay14_InitialCondition_e;   /* Expression: 0
                                        * Referenced by: '<S90>/Delay14'
                                        */
  real_T Delay15_InitialCondition_l;   /* Expression: 0
                                        * Referenced by: '<S90>/Delay15'
                                        */
  real_T Delay16_InitialCondition_g;   /* Expression: 0
                                        * Referenced by: '<S90>/Delay16'
                                        */
  real_T Delay17_InitialCondition_a;   /* Expression: 0
                                        * Referenced by: '<S90>/Delay17'
                                        */
  real_T Delay18_InitialCondition_j;   /* Expression: 0
                                        * Referenced by: '<S90>/Delay18'
                                        */
  real_T Delay19_InitialCondition_c;   /* Expression: 0
                                        * Referenced by: '<S90>/Delay19'
                                        */
  real_T Delay_InitialCondition_k;     /* Expression: 0
                                        * Referenced by: '<S90>/Delay'
                                        */
  real_T Delay1_InitialCondition_g;    /* Expression: 0
                                        * Referenced by: '<S90>/Delay1'
                                        */
  real_T Delay2_InitialCondition_m;    /* Expression: 0
                                        * Referenced by: '<S90>/Delay2'
                                        */
  real_T Delay3_InitialCondition_p;    /* Expression: 0
                                        * Referenced by: '<S90>/Delay3'
                                        */
  real_T Delay4_InitialCondition_j;    /* Expression: 0
                                        * Referenced by: '<S90>/Delay4'
                                        */
  real_T Delay5_InitialCondition_c;    /* Expression: 0
                                        * Referenced by: '<S90>/Delay5'
                                        */
  real_T Delay6_InitialCondition_e;    /* Expression: 0
                                        * Referenced by: '<S90>/Delay6'
                                        */
  real_T Delay7_InitialCondition_e;    /* Expression: 0
                                        * Referenced by: '<S90>/Delay7'
                                        */
  real_T Delay8_InitialCondition_f;    /* Expression: 0
                                        * Referenced by: '<S90>/Delay8'
                                        */
  real_T Delay9_InitialCondition_h;    /* Expression: 0
                                        * Referenced by: '<S90>/Delay9'
                                        */
  real_T Delay10_InitialCondition_c;   /* Expression: 0
                                        * Referenced by: '<S91>/Delay10'
                                        */
  real_T Delay11_InitialCondition_e;   /* Expression: 0
                                        * Referenced by: '<S91>/Delay11'
                                        */
  real_T Delay12_InitialCondition_o;   /* Expression: 0
                                        * Referenced by: '<S91>/Delay12'
                                        */
  real_T Delay13_InitialCondition_g;   /* Expression: 0
                                        * Referenced by: '<S91>/Delay13'
                                        */
  real_T Delay14_InitialCondition_g;   /* Expression: 0
                                        * Referenced by: '<S91>/Delay14'
                                        */
  real_T Delay15_InitialCondition_o;   /* Expression: 0
                                        * Referenced by: '<S91>/Delay15'
                                        */
  real_T Delay16_InitialCondition_o;   /* Expression: 0
                                        * Referenced by: '<S91>/Delay16'
                                        */
  real_T Delay17_InitialCondition_o;   /* Expression: 0
                                        * Referenced by: '<S91>/Delay17'
                                        */
  real_T Delay18_InitialCondition_d;   /* Expression: 0
                                        * Referenced by: '<S91>/Delay18'
                                        */
  real_T Delay19_InitialCondition_a;   /* Expression: 0
                                        * Referenced by: '<S91>/Delay19'
                                        */
  real_T Delay_InitialCondition_nl;    /* Expression: 0
                                        * Referenced by: '<S91>/Delay'
                                        */
  real_T Delay1_InitialCondition_k;    /* Expression: 0
                                        * Referenced by: '<S91>/Delay1'
                                        */
  real_T Delay2_InitialCondition_h;    /* Expression: 0
                                        * Referenced by: '<S91>/Delay2'
                                        */
  real_T Delay3_InitialCondition_o;    /* Expression: 0
                                        * Referenced by: '<S91>/Delay3'
                                        */
  real_T Delay4_InitialCondition_m;    /* Expression: 0
                                        * Referenced by: '<S91>/Delay4'
                                        */
  real_T Delay5_InitialCondition_k;    /* Expression: 0
                                        * Referenced by: '<S91>/Delay5'
                                        */
  real_T Delay6_InitialCondition_g;    /* Expression: 0
                                        * Referenced by: '<S91>/Delay6'
                                        */
  real_T Delay7_InitialCondition_m;    /* Expression: 0
                                        * Referenced by: '<S91>/Delay7'
                                        */
  real_T Delay8_InitialCondition_k;    /* Expression: 0
                                        * Referenced by: '<S91>/Delay8'
                                        */
  real_T Delay9_InitialCondition_e;    /* Expression: 0
                                        * Referenced by: '<S91>/Delay9'
                                        */
  real_T Gain2_Gain_e2;                /* Expression: 1000
                                        * Referenced by: '<S11>/Gain2'
                                        */
  real_T Gain11_Gain_p;                /* Expression: 10^-3
                                        * Referenced by: '<S14>/Gain11'
                                        */
  real_T Gain12_Gain_h;                /* Expression: 10^-3
                                        * Referenced by: '<S14>/Gain12'
                                        */
  real_T Gain13_Gain_m;                /* Expression: 10^-3
                                        * Referenced by: '<S14>/Gain13'
                                        */
  real_T Gain14_Gain_j;                /* Expression: 10^-2
                                        * Referenced by: '<S14>/Gain14'
                                        */
  real_T Gain15_Gain_c;                /* Expression: 10^-2
                                        * Referenced by: '<S14>/Gain15'
                                        */
  real_T Gain16_Gain_j;                /* Expression: 10^-2
                                        * Referenced by: '<S14>/Gain16'
                                        */
  real_T Gain17_Gain;                  /* Expression: 10^-1
                                        * Referenced by: '<S14>/Gain17'
                                        */
  real_T Gain6_Gain_b1;                /* Expression: 10^-2
                                        * Referenced by: '<S14>/Gain6'
                                        */
  real_T Gain7_Gain_e;                 /* Expression: 10^-2
                                        * Referenced by: '<S14>/Gain7'
                                        */
  real_T PulseGenerator_Amp;           /* Expression: 1
                                        * Referenced by: '<S114>/Pulse Generator'
                                        */
  real_T PulseGenerator_Period;     /* Computed Parameter: PulseGenerator_Period
                                     * Referenced by: '<S114>/Pulse Generator'
                                     */
  real_T PulseGenerator_Duty;         /* Computed Parameter: PulseGenerator_Duty
                                       * Referenced by: '<S114>/Pulse Generator'
                                       */
  real_T PulseGenerator_PhaseDelay;    /* Expression: 0
                                        * Referenced by: '<S114>/Pulse Generator'
                                        */
  real_T Gain_Gain_il;                 /* Expression: 1/2
                                        * Referenced by: '<S127>/Gain'
                                        */
  real_T Constant3_Value_k;            /* Expression: 1
                                        * Referenced by: '<S221>/Constant3'
                                        */
  real_T stabilityfactor_noDYC_Value;  /* Expression: 0.00604
                                        * Referenced by: '<S216>/stabilityfactor_noDYC'
                                        */
  real_T Gain6_Gain_nl;                /* Expression: 1/1.53
                                        * Referenced by: '<S221>/Gain6'
                                        */
  real_T TransferFcn1_A;               /* Computed Parameter: TransferFcn1_A
                                        * Referenced by: '<S217>/Transfer Fcn1'
                                        */
  real_T TransferFcn1_C;               /* Computed Parameter: TransferFcn1_C
                                        * Referenced by: '<S217>/Transfer Fcn1'
                                        */
  real_T mstokmh_Gain;                 /* Expression: 3.6
                                        * Referenced by: '<S218>/m//s to km//h'
                                        */
  real_T Saturation_UpperSat_d4;       /* Expression: 120
                                        * Referenced by: '<S218>/Saturation'
                                        */
  real_T Saturation_LowerSat_dn;       /* Expression: 0
                                        * Referenced by: '<S218>/Saturation'
                                        */
  real_T FFGain_0kmh_Value;            /* Expression: 1
                                        * Referenced by: '<S218>/FFGain_0km//h'
                                        */
  real_T FFGain_20kmh_Value;           /* Expression: 1
                                        * Referenced by: '<S218>/FFGain_20km//h'
                                        */
  real_T FFGain_40kmh_Value;           /* Expression: 1
                                        * Referenced by: '<S218>/FFGain_40km//h'
                                        */
  real_T FFGain_60kmh_Value;           /* Expression: 1
                                        * Referenced by: '<S218>/FFGain_60km//h'
                                        */
  real_T FFGain_80kmh_Value;           /* Expression: 1
                                        * Referenced by: '<S218>/FFGain_80km//h'
                                        */
  real_T uDLookupTable_bp01Data_h3[5]; /* Expression: [0 20 40 60 80]
                                        * Referenced by: '<S218>/1-D Lookup Table'
                                        */
  real_T Gain_DYC_FF_Gain;             /* Expression: 265
                                        * Referenced by: '<S218>/Gain_DYC_FF'
                                        */
  real_T mstokmh_Gain_i;               /* Expression: 3.6
                                        * Referenced by: '<S219>/m//s to km//h'
                                        */
  real_T Saturation_UpperSat_pi;       /* Expression: 120
                                        * Referenced by: '<S219>/Saturation'
                                        */
  real_T Saturation_LowerSat_mf;       /* Expression: 0
                                        * Referenced by: '<S219>/Saturation'
                                        */
  real_T FFdiffGain_0kmh_Value;        /* Expression: 1
                                        * Referenced by: '<S219>/FFdiffGain_0km//h'
                                        */
  real_T FFdiffGain_20kmh_Value;       /* Expression: 1
                                        * Referenced by: '<S219>/FFdiffGain_20km//h'
                                        */
  real_T FFdiffGain_40kmh_Value;       /* Expression: 1
                                        * Referenced by: '<S219>/FFdiffGain_40km//h'
                                        */
  real_T FFdiffGain_60kmh_Value;       /* Expression: 1
                                        * Referenced by: '<S219>/FFdiffGain_60km//h'
                                        */
  real_T FFdiffGain_80kmh_Value;       /* Expression: 1
                                        * Referenced by: '<S219>/FFdiffGain_80km//h'
                                        */
  real_T uDLookupTable_bp01Data_lf[5]; /* Expression: [0 20 40 60 80]
                                        * Referenced by: '<S219>/1-D Lookup Table'
                                        */
  real_T TransferFcn2_A;               /* Computed Parameter: TransferFcn2_A
                                        * Referenced by: '<S219>/Transfer Fcn2'
                                        */
  real_T TransferFcn2_C;               /* Computed Parameter: TransferFcn2_C
                                        * Referenced by: '<S219>/Transfer Fcn2'
                                        */
  real_T TransferFcn2_D;               /* Computed Parameter: TransferFcn2_D
                                        * Referenced by: '<S219>/Transfer Fcn2'
                                        */
  real_T I_Gain;                       /* Expression: 112
                                        * Referenced by: '<S219>/I'
                                        */
  real_T FF_diff_Gain_Gain;            /* Expression: 1/3
                                        * Referenced by: '<S219>/FF_diff_Gain'
                                        */
  real_T Gain_Gain_jv;                 /* Expression: 1
                                        * Referenced by: '<S175>/Gain'
                                        */
  real_T Gain1_Gain_m4;                /* Expression: 1
                                        * Referenced by: '<S175>/Gain1'
                                        */
  real_T Gain2_Gain_i;                 /* Expression: 1
                                        * Referenced by: '<S175>/Gain2'
                                        */
  real_T Gain3_Gain_i;                 /* Expression: 1
                                        * Referenced by: '<S175>/Gain3'
                                        */
  real_T Igain_Gain;                   /* Expression: 0
                                        * Referenced by: '<S291>/I gain'
                                        */
  real_T Igain_Gain_b;                 /* Expression: 0
                                        * Referenced by: '<S292>/I gain'
                                        */
  real_T Igain_Gain_k;                 /* Expression: 0
                                        * Referenced by: '<S293>/I gain'
                                        */
  real_T Igain_Gain_l;                 /* Expression: 0
                                        * Referenced by: '<S294>/I gain'
                                        */
  real_T Constant_Value_o3;            /* Expression: -21.188
                                        * Referenced by: '<S338>/Constant'
                                        */
  real_T Constant1_Value_m;            /* Expression: 105.95
                                        * Referenced by: '<S338>/Constant1'
                                        */
  uint32_T uDLookupTable_maxIndex[2];
                                   /* Computed Parameter: uDLookupTable_maxIndex
                                    * Referenced by: '<S186>/2-D Lookup Table'
                                    */
  uint32_T uDLookupTable_maxIndex_k[2];
                                 /* Computed Parameter: uDLookupTable_maxIndex_k
                                  * Referenced by: '<S185>/2-D Lookup Table'
                                  */
  uint32_T uDLookupTable_maxIndex_p[2];
                                 /* Computed Parameter: uDLookupTable_maxIndex_p
                                  * Referenced by: '<S188>/2-D Lookup Table'
                                  */
  uint32_T uDLookupTable_maxIndex_h[2];
                                 /* Computed Parameter: uDLookupTable_maxIndex_h
                                  * Referenced by: '<S187>/2-D Lookup Table'
                                  */
  uint32_T uDLookupTable_maxIndex_n[2];
                                 /* Computed Parameter: uDLookupTable_maxIndex_n
                                  * Referenced by: '<S190>/2-D Lookup Table'
                                  */
  uint32_T uDLookupTable_maxIndex_pw[2];
                                /* Computed Parameter: uDLookupTable_maxIndex_pw
                                 * Referenced by: '<S189>/2-D Lookup Table'
                                 */
  uint32_T uDLookupTable_maxIndex_j[2];
                                 /* Computed Parameter: uDLookupTable_maxIndex_j
                                  * Referenced by: '<S192>/2-D Lookup Table'
                                  */
  uint32_T uDLookupTable_maxIndex_pr[2];
                                /* Computed Parameter: uDLookupTable_maxIndex_pr
                                 * Referenced by: '<S191>/2-D Lookup Table'
                                 */
  uint32_T uDLookupTable_maxIndex_a[2];
                                 /* Computed Parameter: uDLookupTable_maxIndex_a
                                  * Referenced by: '<S194>/2-D Lookup Table'
                                  */
  uint32_T uDLookupTable_maxIndex_nc[2];
                                /* Computed Parameter: uDLookupTable_maxIndex_nc
                                 * Referenced by: '<S193>/2-D Lookup Table'
                                 */
  uint32_T uDLookupTable_maxIndex_pq[2];
                                /* Computed Parameter: uDLookupTable_maxIndex_pq
                                 * Referenced by: '<S196>/2-D Lookup Table'
                                 */
  uint32_T uDLookupTable_maxIndex_g[2];
                                 /* Computed Parameter: uDLookupTable_maxIndex_g
                                  * Referenced by: '<S195>/2-D Lookup Table'
                                  */
  uint32_T uDLookupTable_maxIndex_k1[2];
                                /* Computed Parameter: uDLookupTable_maxIndex_k1
                                 * Referenced by: '<S198>/2-D Lookup Table'
                                 */
  uint32_T uDLookupTable_maxIndex_l[2];
                                 /* Computed Parameter: uDLookupTable_maxIndex_l
                                  * Referenced by: '<S197>/2-D Lookup Table'
                                  */
  uint32_T uDLookupTable_maxIndex_nr[2];
                                /* Computed Parameter: uDLookupTable_maxIndex_nr
                                 * Referenced by: '<S200>/2-D Lookup Table'
                                 */
  uint32_T uDLookupTable_maxIndex_b[2];
                                 /* Computed Parameter: uDLookupTable_maxIndex_b
                                  * Referenced by: '<S199>/2-D Lookup Table'
                                  */
  uint32_T uDLookupTable_maxIndex_ga[2];
                                /* Computed Parameter: uDLookupTable_maxIndex_ga
                                 * Referenced by: '<S202>/2-D Lookup Table'
                                 */
  uint32_T uDLookupTable_maxIndex_o[2];
                                 /* Computed Parameter: uDLookupTable_maxIndex_o
                                  * Referenced by: '<S201>/2-D Lookup Table'
                                  */
  uint32_T uDLookupTable_maxIndex_ln[2];
                                /* Computed Parameter: uDLookupTable_maxIndex_ln
                                 * Referenced by: '<S204>/2-D Lookup Table'
                                 */
  uint32_T uDLookupTable_maxIndex_of[2];
                                /* Computed Parameter: uDLookupTable_maxIndex_of
                                 * Referenced by: '<S203>/2-D Lookup Table'
                                 */
  uint32_T uDLookupTable_maxIndex_e;
                                 /* Computed Parameter: uDLookupTable_maxIndex_e
                                  * Referenced by: '<S218>/1-D Lookup Table'
                                  */
  uint32_T uDLookupTable_dimSizes; /* Computed Parameter: uDLookupTable_dimSizes
                                    * Referenced by: '<S218>/1-D Lookup Table'
                                    */
  uint32_T uDLookupTable_numYWorkElts[2];
                               /* Computed Parameter: uDLookupTable_numYWorkElts
                                * Referenced by: '<S218>/1-D Lookup Table'
                                */
  uint32_T uDLookupTable_maxIndex_oo;
                                /* Computed Parameter: uDLookupTable_maxIndex_oo
                                 * Referenced by: '<S219>/1-D Lookup Table'
                                 */
  uint32_T uDLookupTable_dimSizes_e;
                                 /* Computed Parameter: uDLookupTable_dimSizes_e
                                  * Referenced by: '<S219>/1-D Lookup Table'
                                  */
  uint32_T uDLookupTable_numYWorkElts_i[2];
                             /* Computed Parameter: uDLookupTable_numYWorkElts_i
                              * Referenced by: '<S219>/1-D Lookup Table'
                              */
  boolean_T Delay_InitialCondition_ge;
                                /* Computed Parameter: Delay_InitialCondition_ge
                                 * Referenced by: '<S82>/Delay'
                                 */
  boolean_T Delay1_InitialCondition_e;
                                /* Computed Parameter: Delay1_InitialCondition_e
                                 * Referenced by: '<S82>/Delay1'
                                 */
  boolean_T Delay2_InitialCondition_hk;
                               /* Computed Parameter: Delay2_InitialCondition_hk
                                * Referenced by: '<S82>/Delay2'
                                */
  boolean_T Delay3_InitialCondition_c;
                                /* Computed Parameter: Delay3_InitialCondition_c
                                 * Referenced by: '<S82>/Delay3'
                                 */
  boolean_T Delay4_InitialCondition_l;
                                /* Computed Parameter: Delay4_InitialCondition_l
                                 * Referenced by: '<S82>/Delay4'
                                 */
  boolean_T Delay_InitialCondition_lq;
                                /* Computed Parameter: Delay_InitialCondition_lq
                                 * Referenced by: '<S127>/Delay'
                                 */
  boolean_T Delay_InitialCondition_f;
                                 /* Computed Parameter: Delay_InitialCondition_f
                                  * Referenced by: '<S92>/Delay'
                                  */
  boolean_T Delay1_InitialCondition_j;
                                /* Computed Parameter: Delay1_InitialCondition_j
                                 * Referenced by: '<S92>/Delay1'
                                 */
  boolean_T Delay2_InitialCondition_n;
                                /* Computed Parameter: Delay2_InitialCondition_n
                                 * Referenced by: '<S92>/Delay2'
                                 */
  boolean_T Delay3_InitialCondition_b;
                                /* Computed Parameter: Delay3_InitialCondition_b
                                 * Referenced by: '<S92>/Delay3'
                                 */
  boolean_T Delay4_InitialCondition_a;
                                /* Computed Parameter: Delay4_InitialCondition_a
                                 * Referenced by: '<S92>/Delay4'
                                 */
  boolean_T Delay5_InitialCondition_b;
                                /* Computed Parameter: Delay5_InitialCondition_b
                                 * Referenced by: '<S92>/Delay5'
                                 */
  boolean_T Delay6_InitialCondition_i;
                                /* Computed Parameter: Delay6_InitialCondition_i
                                 * Referenced by: '<S92>/Delay6'
                                 */
  boolean_T Delay7_InitialCondition_h;
                                /* Computed Parameter: Delay7_InitialCondition_h
                                 * Referenced by: '<S92>/Delay7'
                                 */
  boolean_T Delay8_InitialCondition_o;
                                /* Computed Parameter: Delay8_InitialCondition_o
                                 * Referenced by: '<S92>/Delay8'
                                 */
  boolean_T Delay9_InitialCondition_f;
                                /* Computed Parameter: Delay9_InitialCondition_f
                                 * Referenced by: '<S92>/Delay9'
                                 */
  boolean_T Delay_InitialCondition_m3;
                                /* Computed Parameter: Delay_InitialCondition_m3
                                 * Referenced by: '<S79>/Delay'
                                 */
  boolean_T Delay1_InitialCondition_n;
                                /* Computed Parameter: Delay1_InitialCondition_n
                                 * Referenced by: '<S79>/Delay1'
                                 */
  boolean_T Delay2_InitialCondition_nz;
                               /* Computed Parameter: Delay2_InitialCondition_nz
                                * Referenced by: '<S79>/Delay2'
                                */
  boolean_T Delay3_InitialCondition_j;
                                /* Computed Parameter: Delay3_InitialCondition_j
                                 * Referenced by: '<S79>/Delay3'
                                 */
  boolean_T Delay4_InitialCondition_al;
                               /* Computed Parameter: Delay4_InitialCondition_al
                                * Referenced by: '<S79>/Delay4'
                                */
  boolean_T Delay5_InitialCondition_m;
                                /* Computed Parameter: Delay5_InitialCondition_m
                                 * Referenced by: '<S79>/Delay5'
                                 */
  boolean_T Delay6_InitialCondition_d;
                                /* Computed Parameter: Delay6_InitialCondition_d
                                 * Referenced by: '<S79>/Delay6'
                                 */
  boolean_T Delay7_InitialCondition_f;
                                /* Computed Parameter: Delay7_InitialCondition_f
                                 * Referenced by: '<S79>/Delay7'
                                 */
  boolean_T Delay8_InitialCondition_m;
                                /* Computed Parameter: Delay8_InitialCondition_m
                                 * Referenced by: '<S79>/Delay8'
                                 */
  boolean_T Delay9_InitialCondition_fa;
                               /* Computed Parameter: Delay9_InitialCondition_fa
                                * Referenced by: '<S79>/Delay9'
                                */
  boolean_T Delay_InitialCondition_kq;
                                /* Computed Parameter: Delay_InitialCondition_kq
                                 * Referenced by: '<S80>/Delay'
                                 */
  boolean_T Delay1_InitialCondition_l;
                                /* Computed Parameter: Delay1_InitialCondition_l
                                 * Referenced by: '<S80>/Delay1'
                                 */
  boolean_T Delay2_InitialCondition_j;
                                /* Computed Parameter: Delay2_InitialCondition_j
                                 * Referenced by: '<S80>/Delay2'
                                 */
  boolean_T Delay3_InitialCondition_i;
                                /* Computed Parameter: Delay3_InitialCondition_i
                                 * Referenced by: '<S80>/Delay3'
                                 */
  boolean_T Delay4_InitialCondition_i;
                                /* Computed Parameter: Delay4_InitialCondition_i
                                 * Referenced by: '<S80>/Delay4'
                                 */
  boolean_T Delay5_InitialCondition_e;
                                /* Computed Parameter: Delay5_InitialCondition_e
                                 * Referenced by: '<S80>/Delay5'
                                 */
  boolean_T Delay6_InitialCondition_k;
                                /* Computed Parameter: Delay6_InitialCondition_k
                                 * Referenced by: '<S80>/Delay6'
                                 */
  boolean_T Delay7_InitialCondition_mc;
                               /* Computed Parameter: Delay7_InitialCondition_mc
                                * Referenced by: '<S80>/Delay7'
                                */
  boolean_T Delay8_InitialCondition_i;
                                /* Computed Parameter: Delay8_InitialCondition_i
                                 * Referenced by: '<S80>/Delay8'
                                 */
  boolean_T Delay9_InitialCondition_hf;
                               /* Computed Parameter: Delay9_InitialCondition_hf
                                * Referenced by: '<S80>/Delay9'
                                */
  boolean_T Delay_InitialCondition_e;
                                 /* Computed Parameter: Delay_InitialCondition_e
                                  * Referenced by: '<S114>/Delay'
                                  */
};

/* Real-time Model Data Structure */
struct tag_RTM_VCM20_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_VCM20_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeF[1][8];
  ODE1_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    struct {
      uint8_T TID[4];
    } TaskCounters;

    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[4];
  } Timing;
};

/* Block parameters (default storage) */
extern P_VCM20_T VCM20_P;

/* Block signals (default storage) */
extern B_VCM20_T VCM20_B;

/* Continuous states (default storage) */
extern X_VCM20_T VCM20_X;

/* Block states (default storage) */
extern DW_VCM20_T VCM20_DW;

/* Model entry point functions */
extern void VCM20_initialize(void);
extern void VCM20_output(void);
extern void VCM20_update(void);
extern void VCM20_terminate(void);

/* Real-time Model object */
extern RT_MODEL_VCM20_T *const VCM20_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'VCM20'
 * '<S1>'   : 'VCM20/CAN Bundle'
 * '<S2>'   : 'VCM20/CAN from AMS'
 * '<S3>'   : 'VCM20/CAN from Inv'
 * '<S4>'   : 'VCM20/CAN to EVO4S'
 * '<S5>'   : 'VCM20/CAN to Inv'
 * '<S6>'   : 'VCM20/CAN_TYPE1_SETUP_M1_C1'
 * '<S7>'   : 'VCM20/CAN_TYPE1_SETUP_M1_C2'
 * '<S8>'   : 'VCM20/CAN_TYPE1_SETUP_M2_C2'
 * '<S9>'   : 'VCM20/CAN_TYPE1_SETUP_M3_C1'
 * '<S10>'  : 'VCM20/CAN_TYPE1_SETUP_M3_C2'
 * '<S11>'  : 'VCM20/Driver Input'
 * '<S12>'  : 'VCM20/From relay board'
 * '<S13>'  : 'VCM20/RTI Data'
 * '<S14>'  : 'VCM20/Sensor'
 * '<S15>'  : 'VCM20/VCM main'
 * '<S16>'  : 'VCM20/VCMOutput'
 * '<S17>'  : 'VCM20/__CAN_TYPE1_SETUP_M2_C1'
 * '<S18>'  : 'VCM20/sensors from cooling'
 * '<S19>'  : 'VCM20/speedlimit'
 * '<S20>'  : 'VCM20/CAN from AMS/CAN from AMS_401'
 * '<S21>'  : 'VCM20/CAN from Inv/MC_11'
 * '<S22>'  : 'VCM20/CAN from Inv/MC_12'
 * '<S23>'  : 'VCM20/CAN from Inv/MC_13'
 * '<S24>'  : 'VCM20/CAN from Inv/MC_21'
 * '<S25>'  : 'VCM20/CAN from Inv/MC_22'
 * '<S26>'  : 'VCM20/CAN from Inv/MC_23'
 * '<S27>'  : 'VCM20/CAN from Inv/MC_31'
 * '<S28>'  : 'VCM20/CAN from Inv/MC_32'
 * '<S29>'  : 'VCM20/CAN from Inv/MC_33'
 * '<S30>'  : 'VCM20/CAN from Inv/MC_41'
 * '<S31>'  : 'VCM20/CAN from Inv/MC_42'
 * '<S32>'  : 'VCM20/CAN from Inv/MC_43'
 * '<S33>'  : 'VCM20/CAN to EVO4S/AMS to EVO4s'
 * '<S34>'  : 'VCM20/CAN to EVO4S/Driver to EVO4s'
 * '<S35>'  : 'VCM20/CAN to EVO4S/IMU to EVO4s'
 * '<S36>'  : 'VCM20/CAN to EVO4S/Inv to EVO4s'
 * '<S37>'  : 'VCM20/CAN to EVO4S/VCM to EVO4s'
 * '<S38>'  : 'VCM20/CAN to EVO4S/cooler to EVO4s'
 * '<S39>'  : 'VCM20/CAN to EVO4S/relay to EVO4s'
 * '<S40>'  : 'VCM20/CAN to EVO4S/AMS to EVO4s/to EVO4s 702'
 * '<S41>'  : 'VCM20/CAN to EVO4S/Driver to EVO4s/to EVO4s 503'
 * '<S42>'  : 'VCM20/CAN to EVO4S/Driver to EVO4s/to EVO4s 705'
 * '<S43>'  : 'VCM20/CAN to EVO4S/IMU to EVO4s/to EVO4s 50B'
 * '<S44>'  : 'VCM20/CAN to EVO4S/IMU to EVO4s/to EVO4s 50C'
 * '<S45>'  : 'VCM20/CAN to EVO4S/IMU to EVO4s/to EVO4s 50D'
 * '<S46>'  : 'VCM20/CAN to EVO4S/IMU to EVO4s/to EVO4s 50E'
 * '<S47>'  : 'VCM20/CAN to EVO4S/IMU to EVO4s/to EVO4s 50F'
 * '<S48>'  : 'VCM20/CAN to EVO4S/IMU to EVO4s/to EVO4s 510'
 * '<S49>'  : 'VCM20/CAN to EVO4S/Inv to EVO4s/to EVO4s 501'
 * '<S50>'  : 'VCM20/CAN to EVO4S/Inv to EVO4s/to EVO4s 502'
 * '<S51>'  : 'VCM20/CAN to EVO4S/Inv to EVO4s/to EVO4s 505'
 * '<S52>'  : 'VCM20/CAN to EVO4S/Inv to EVO4s/to EVO4s 506'
 * '<S53>'  : 'VCM20/CAN to EVO4S/Inv to EVO4s/to EVO4s 507'
 * '<S54>'  : 'VCM20/CAN to EVO4S/Inv to EVO4s/to EVO4s 508'
 * '<S55>'  : 'VCM20/CAN to EVO4S/Inv to EVO4s/to EVO4s 509'
 * '<S56>'  : 'VCM20/CAN to EVO4S/Inv to EVO4s/to EVO4s 50A'
 * '<S57>'  : 'VCM20/CAN to EVO4S/VCM to EVO4s/to EVO4s 504'
 * '<S58>'  : 'VCM20/CAN to EVO4S/VCM to EVO4s/to EVO4s 521'
 * '<S59>'  : 'VCM20/CAN to EVO4S/VCM to EVO4s/to EVO4s 522'
 * '<S60>'  : 'VCM20/CAN to EVO4S/VCM to EVO4s/to EVO4s 701'
 * '<S61>'  : 'VCM20/CAN to EVO4S/cooler to EVO4s/to EVO4s 703'
 * '<S62>'  : 'VCM20/CAN to EVO4S/relay to EVO4s/to EVO4s 704'
 * '<S63>'  : 'VCM20/CAN to Inv/VCM_11'
 * '<S64>'  : 'VCM20/CAN to Inv/VCM_21'
 * '<S65>'  : 'VCM20/CAN to Inv/VCM_31'
 * '<S66>'  : 'VCM20/CAN to Inv/VCM_41'
 * '<S67>'  : 'VCM20/Driver Input/APPS_SafetyFilter'
 * '<S68>'  : 'VCM20/Driver Input/APPS_Steer'
 * '<S69>'  : 'VCM20/Driver Input/Brakes'
 * '<S70>'  : 'VCM20/Driver Input/Calibration1'
 * '<S71>'  : 'VCM20/Driver Input/Calibration2'
 * '<S72>'  : 'VCM20/Driver Input/Calibration3'
 * '<S73>'  : 'VCM20/Driver Input/Calibration4'
 * '<S74>'  : 'VCM20/Driver Input/Calibration5'
 * '<S75>'  : 'VCM20/Driver Input/Calibration6'
 * '<S76>'  : 'VCM20/Driver Input/Calibration7'
 * '<S77>'  : 'VCM20/Driver Input/SWAS_SafetyFilter1'
 * '<S78>'  : 'VCM20/Driver Input/Steer SW'
 * '<S79>'  : 'VCM20/Driver Input/APPS_SafetyFilter/10ms_and1'
 * '<S80>'  : 'VCM20/Driver Input/APPS_SafetyFilter/10ms_and2'
 * '<S81>'  : 'VCM20/Driver Input/APPS_SafetyFilter/10ms_or'
 * '<S82>'  : 'VCM20/Driver Input/APPS_SafetyFilter/5ms_and'
 * '<S83>'  : 'VCM20/Driver Input/APPS_SafetyFilter/Compare To Constant'
 * '<S84>'  : 'VCM20/Driver Input/APPS_SafetyFilter/Compare To Constant1'
 * '<S85>'  : 'VCM20/Driver Input/APPS_SafetyFilter/Compare To Constant2'
 * '<S86>'  : 'VCM20/Driver Input/APPS_SafetyFilter/Compare To Constant3'
 * '<S87>'  : 'VCM20/Driver Input/APPS_SafetyFilter/Compare To Constant4'
 * '<S88>'  : 'VCM20/Driver Input/APPS_SafetyFilter/10ms_or/Compare To Constant'
 * '<S89>'  : 'VCM20/Driver Input/Calibration1/APPS1_Scope'
 * '<S90>'  : 'VCM20/Driver Input/Calibration2/APPS2_Scope'
 * '<S91>'  : 'VCM20/Driver Input/Calibration7/BrakePressF_Scope'
 * '<S92>'  : 'VCM20/Driver Input/SWAS_SafetyFilter1/10ms_and1'
 * '<S93>'  : 'VCM20/Driver Input/SWAS_SafetyFilter1/Compare To Constant1'
 * '<S94>'  : 'VCM20/Driver Input/SWAS_SafetyFilter1/Compare To Constant3'
 * '<S95>'  : 'VCM20/From relay board/AMS_shutdownsig'
 * '<S96>'  : 'VCM20/From relay board/BSPD_shutdownsig'
 * '<S97>'  : 'VCM20/From relay board/IMD_shutdownsig'
 * '<S98>'  : 'VCM20/RTI Data/RTI Data Store'
 * '<S99>'  : 'VCM20/RTI Data/RTI Data Store/RTI Data Store'
 * '<S100>' : 'VCM20/RTI Data/RTI Data Store/RTI Data Store/RTI Data Store'
 * '<S101>' : 'VCM20/Sensor/CAN from EVO4_1'
 * '<S102>' : 'VCM20/Sensor/CAN from EVO4_201'
 * '<S103>' : 'VCM20/Sensor/CAN_from_IMU_100'
 * '<S104>' : 'VCM20/Sensor/CAN_from_IMU_102'
 * '<S105>' : 'VCM20/Sensor/CAN_from_IMU_121'
 * '<S106>' : 'VCM20/Sensor/CAN_from_IMU_122'
 * '<S107>' : 'VCM20/Sensor/CAN_from_IMU_134'
 * '<S108>' : 'VCM20/Sensor/CAN_from_IMU_139'
 * '<S109>' : 'VCM20/Sensor/CAN_from_IMU_151'
 * '<S110>' : 'VCM20/Sensor/CAN_from_IMU_171'
 * '<S111>' : 'VCM20/Sensor/CAN_from_IMU_220'
 * '<S112>' : 'VCM20/Sensor/Current_mV'
 * '<S113>' : 'VCM20/VCM main/All Inv enable'
 * '<S114>' : 'VCM20/VCM main/Auto Error Reset'
 * '<S115>' : 'VCM20/VCM main/Brake Over Ride'
 * '<S116>' : 'VCM20/VCM main/Sequence'
 * '<S117>' : 'VCM20/VCM main/Speed Torque MC1'
 * '<S118>' : 'VCM20/VCM main/Speed Torque MC2'
 * '<S119>' : 'VCM20/VCM main/Speed Torque MC3'
 * '<S120>' : 'VCM20/VCM main/Speed Torque MC4'
 * '<S121>' : 'VCM20/VCM main/Trq Clc '
 * '<S122>' : 'VCM20/VCM main/VDC'
 * '<S123>' : 'VCM20/VCM main/bEnableOn_SwBox'
 * '<S124>' : 'VCM20/VCM main/bInberterOn_SwBox'
 * '<S125>' : 'VCM20/VCM main/torque limit_SwBox'
 * '<S126>' : 'VCM20/VCM main/All Inv enable/All Inv enable'
 * '<S127>' : 'VCM20/VCM main/Brake Over Ride/Brake Over Ride judge'
 * '<S128>' : 'VCM20/VCM main/Brake Over Ride/Brake Over Ride judge/Compare To Zero'
 * '<S129>' : 'VCM20/VCM main/Brake Over Ride/Brake Over Ride judge/Compare To Zero1'
 * '<S130>' : 'VCM20/VCM main/Sequence/Compare To Constant'
 * '<S131>' : 'VCM20/VCM main/Sequence/bError_SwBox'
 * '<S132>' : 'VCM20/VCM main/Sequence/bInberterOn_SwBox1'
 * '<S133>' : 'VCM20/VCM main/Sequence/bQuitDcOn_SwBox'
 * '<S134>' : 'VCM20/VCM main/Sequence/bSystemReady_SwBox'
 * '<S135>' : 'VCM20/VCM main/Speed Torque MC1/Power limit Negative'
 * '<S136>' : 'VCM20/VCM main/Speed Torque MC1/Power limit Positive'
 * '<S137>' : 'VCM20/VCM main/Speed Torque MC1/Subsystem'
 * '<S138>' : 'VCM20/VCM main/Speed Torque MC1/TargetVelocity'
 * '<S139>' : 'VCM20/VCM main/Speed Torque MC1/Torque delay N'
 * '<S140>' : 'VCM20/VCM main/Speed Torque MC1/Torque delay P'
 * '<S141>' : 'VCM20/VCM main/Speed Torque MC1/Power limit Negative/Saturation Dynamic'
 * '<S142>' : 'VCM20/VCM main/Speed Torque MC1/Power limit Positive/Saturation Dynamic'
 * '<S143>' : 'VCM20/VCM main/Speed Torque MC2/Power limit Negative'
 * '<S144>' : 'VCM20/VCM main/Speed Torque MC2/Power limit Positive'
 * '<S145>' : 'VCM20/VCM main/Speed Torque MC2/Subsystem'
 * '<S146>' : 'VCM20/VCM main/Speed Torque MC2/TargetVelocity'
 * '<S147>' : 'VCM20/VCM main/Speed Torque MC2/Torque delay N'
 * '<S148>' : 'VCM20/VCM main/Speed Torque MC2/Torque delay P'
 * '<S149>' : 'VCM20/VCM main/Speed Torque MC2/Power limit Negative/Saturation Dynamic'
 * '<S150>' : 'VCM20/VCM main/Speed Torque MC2/Power limit Positive/Saturation Dynamic'
 * '<S151>' : 'VCM20/VCM main/Speed Torque MC3/Power limit Negative'
 * '<S152>' : 'VCM20/VCM main/Speed Torque MC3/Power limit Positive'
 * '<S153>' : 'VCM20/VCM main/Speed Torque MC3/Subsystem'
 * '<S154>' : 'VCM20/VCM main/Speed Torque MC3/TargetVelocity'
 * '<S155>' : 'VCM20/VCM main/Speed Torque MC3/Torque delay N'
 * '<S156>' : 'VCM20/VCM main/Speed Torque MC3/Torque delay P'
 * '<S157>' : 'VCM20/VCM main/Speed Torque MC3/Power limit Negative/Saturation Dynamic'
 * '<S158>' : 'VCM20/VCM main/Speed Torque MC3/Power limit Positive/Saturation Dynamic'
 * '<S159>' : 'VCM20/VCM main/Speed Torque MC4/Power limit Negative'
 * '<S160>' : 'VCM20/VCM main/Speed Torque MC4/Power limit Positive'
 * '<S161>' : 'VCM20/VCM main/Speed Torque MC4/Subsystem'
 * '<S162>' : 'VCM20/VCM main/Speed Torque MC4/TargetVelocity'
 * '<S163>' : 'VCM20/VCM main/Speed Torque MC4/Torque delay N'
 * '<S164>' : 'VCM20/VCM main/Speed Torque MC4/Torque delay P'
 * '<S165>' : 'VCM20/VCM main/Speed Torque MC4/Power limit Negative/Saturation Dynamic'
 * '<S166>' : 'VCM20/VCM main/Speed Torque MC4/Power limit Positive/Saturation Dynamic'
 * '<S167>' : 'VCM20/VCM main/Trq Clc /Avoid division by zero1'
 * '<S168>' : 'VCM20/VCM main/Trq Clc /Saturation Dynamic'
 * '<S169>' : 'VCM20/VCM main/Trq Clc /Avoid division by zero1/Compare To Zero'
 * '<S170>' : 'VCM20/VCM main/VDC/80kw_limitation'
 * '<S171>' : 'VCM20/VCM main/VDC/DYC'
 * '<S172>' : 'VCM20/VCM main/VDC/Launch'
 * '<S173>' : 'VCM20/VCM main/VDC/PowerConsumption'
 * '<S174>' : 'VCM20/VCM main/VDC/Torque Front//rear distributor Main'
 * '<S175>' : 'VCM20/VCM main/VDC/Traction Control'
 * '<S176>' : 'VCM20/VCM main/VDC/VCM_Status'
 * '<S177>' : 'VCM20/VCM main/VDC/80kw_limitation/In signal select'
 * '<S178>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_1'
 * '<S179>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_2'
 * '<S180>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_3'
 * '<S181>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_4'
 * '<S182>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_5'
 * '<S183>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_last'
 * '<S184>' : 'VCM20/VCM main/VDC/80kw_limitation/Subsystem'
 * '<S185>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_1/FL Power Culc'
 * '<S186>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_1/FR Power Culc'
 * '<S187>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_1/RL Power Culc'
 * '<S188>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_1/RR Power Culc'
 * '<S189>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_2/FL Power Culc'
 * '<S190>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_2/FR Power Culc'
 * '<S191>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_2/RL Power Culc'
 * '<S192>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_2/RR Power Culc'
 * '<S193>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_3/FL Power Culc'
 * '<S194>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_3/FR Power Culc'
 * '<S195>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_3/RL Power Culc'
 * '<S196>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_3/RR Power Culc'
 * '<S197>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_4/FL Power Culc'
 * '<S198>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_4/FR Power Culc'
 * '<S199>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_4/RL Power Culc'
 * '<S200>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_4/RR Power Culc'
 * '<S201>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_5/FL Power Culc'
 * '<S202>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_5/FR Power Culc'
 * '<S203>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_5/RL Power Culc'
 * '<S204>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_5/RR Power Culc'
 * '<S205>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_last/80kmLimitation_@HighSpeed'
 * '<S206>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_last/Compare To Zero'
 * '<S207>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_last/80kmLimitation_@HighSpeed/Compare To Constant'
 * '<S208>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_last/80kmLimitation_@HighSpeed/Compare To Constant1'
 * '<S209>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_last/80kmLimitation_@HighSpeed/Compare To Constant2'
 * '<S210>' : 'VCM20/VCM main/VDC/80kw_limitation/Power Limit_last/80kmLimitation_@HighSpeed/Compare To Constant3'
 * '<S211>' : 'VCM20/VCM main/VDC/DYC/Moment to É¢T'
 * '<S212>' : 'VCM20/VCM main/VDC/DYC/Torque distributor1'
 * '<S213>' : 'VCM20/VCM main/VDC/DYC/Torque distributor2'
 * '<S214>' : 'VCM20/VCM main/VDC/DYC/Torque distributor3'
 * '<S215>' : 'VCM20/VCM main/VDC/DYC/Torque distributor4'
 * '<S216>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal'
 * '<S217>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB'
 * '<S218>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FF_Gain'
 * '<S219>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FF_diff_Gain'
 * '<S220>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/Velocity Clc'
 * '<S221>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/É¬ to r_noDYC'
 * '<S222>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/É¬ to r_ref'
 * '<S223>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller'
 * '<S224>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Anti-windup'
 * '<S225>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/D Gain'
 * '<S226>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Filter'
 * '<S227>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Filter ICs'
 * '<S228>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/I Gain'
 * '<S229>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Ideal P Gain'
 * '<S230>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Ideal P Gain Fdbk'
 * '<S231>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Integrator'
 * '<S232>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Integrator ICs'
 * '<S233>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/N Copy'
 * '<S234>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/N Gain'
 * '<S235>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/P Copy'
 * '<S236>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Parallel P Gain'
 * '<S237>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Reset Signal'
 * '<S238>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Saturation'
 * '<S239>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Saturation Fdbk'
 * '<S240>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Sum'
 * '<S241>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Sum Fdbk'
 * '<S242>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Tracking Mode'
 * '<S243>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Tracking Mode Sum'
 * '<S244>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Tsamp - Integral'
 * '<S245>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Tsamp - Ngain'
 * '<S246>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/postSat Signal'
 * '<S247>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/preSat Signal'
 * '<S248>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Anti-windup/Passthrough'
 * '<S249>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/D Gain/Internal Parameters'
 * '<S250>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Filter/Cont. Filter'
 * '<S251>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Filter ICs/Internal IC - Filter'
 * '<S252>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/I Gain/Internal Parameters'
 * '<S253>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Ideal P Gain/Passthrough'
 * '<S254>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S255>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Integrator/Continuous'
 * '<S256>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Integrator ICs/Internal IC'
 * '<S257>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/N Copy/Disabled'
 * '<S258>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/N Gain/Internal Parameters'
 * '<S259>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/P Copy/Disabled'
 * '<S260>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Parallel P Gain/Internal Parameters'
 * '<S261>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Reset Signal/Disabled'
 * '<S262>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Saturation/Passthrough'
 * '<S263>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Saturation Fdbk/Disabled'
 * '<S264>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Sum/Sum_PID'
 * '<S265>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Sum Fdbk/Disabled'
 * '<S266>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Tracking Mode/Disabled'
 * '<S267>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Tracking Mode Sum/Passthrough'
 * '<S268>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Tsamp - Integral/Passthrough'
 * '<S269>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/Tsamp - Ngain/Passthrough'
 * '<S270>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/postSat Signal/Forward_Path'
 * '<S271>' : 'VCM20/VCM main/VDC/DYC/YawmomentCal/FB/PID Controller/preSat Signal/Forward_Path'
 * '<S272>' : 'VCM20/VCM main/VDC/Launch/Chart'
 * '<S273>' : 'VCM20/VCM main/VDC/Launch/Saturation Dynamic1'
 * '<S274>' : 'VCM20/VCM main/VDC/Launch/Saturation Dynamic2'
 * '<S275>' : 'VCM20/VCM main/VDC/Launch/Saturation Dynamic3'
 * '<S276>' : 'VCM20/VCM main/VDC/Launch/Saturation Dynamic4'
 * '<S277>' : 'VCM20/VCM main/VDC/Launch/max rpm'
 * '<S278>' : 'VCM20/VCM main/VDC/Torque Front//rear distributor Main/Max Torque Limit'
 * '<S279>' : 'VCM20/VCM main/VDC/Torque Front//rear distributor Main/Min Torque Limit'
 * '<S280>' : 'VCM20/VCM main/VDC/Torque Front//rear distributor Main/Saturation'
 * '<S281>' : 'VCM20/VCM main/VDC/Torque Front//rear distributor Main/Torque Front//rear distributor '
 * '<S282>' : 'VCM20/VCM main/VDC/Torque Front//rear distributor Main/Saturation/Saturation Dynamic1'
 * '<S283>' : 'VCM20/VCM main/VDC/Traction Control/FF'
 * '<S284>' : 'VCM20/VCM main/VDC/Traction Control/Front TargetSpeedClc_inBrake'
 * '<S285>' : 'VCM20/VCM main/VDC/Traction Control/Rear TargetSpeedClc_inBrake'
 * '<S286>' : 'VCM20/VCM main/VDC/Traction Control/Saturation Dynamic1'
 * '<S287>' : 'VCM20/VCM main/VDC/Traction Control/Saturation Dynamic2'
 * '<S288>' : 'VCM20/VCM main/VDC/Traction Control/Saturation Dynamic3'
 * '<S289>' : 'VCM20/VCM main/VDC/Traction Control/Saturation Dynamic4'
 * '<S290>' : 'VCM20/VCM main/VDC/Traction Control/TargetSpeed_inAcc'
 * '<S291>' : 'VCM20/VCM main/VDC/Traction Control/TractionControll1'
 * '<S292>' : 'VCM20/VCM main/VDC/Traction Control/TractionControll2'
 * '<S293>' : 'VCM20/VCM main/VDC/Traction Control/TractionControll3'
 * '<S294>' : 'VCM20/VCM main/VDC/Traction Control/TractionControll4'
 * '<S295>' : 'VCM20/VCM main/VDC/Traction Control/FF/1//(1-SlipRate)'
 * '<S296>' : 'VCM20/VCM main/VDC/Traction Control/FF/FR Load Diff'
 * '<S297>' : 'VCM20/VCM main/VDC/Traction Control/FF/Tire Force F'
 * '<S298>' : 'VCM20/VCM main/VDC/Traction Control/FF/Tire Force F1'
 * '<S299>' : 'VCM20/VCM main/VDC/Traction Control/FF/Tire Force F/Compare To Zero'
 * '<S300>' : 'VCM20/VCM main/VDC/Traction Control/FF/Tire Force F1/Compare To Zero'
 * '<S301>' : 'VCM20/VCM main/VDC/Traction Control/Front TargetSpeedClc_inBrake/Saturation Dynamic'
 * '<S302>' : 'VCM20/VCM main/VDC/Traction Control/Front TargetSpeedClc_inBrake/toTargetspeed'
 * '<S303>' : 'VCM20/VCM main/VDC/Traction Control/Rear TargetSpeedClc_inBrake/Saturation Dynamic'
 * '<S304>' : 'VCM20/VCM main/VDC/Traction Control/Rear TargetSpeedClc_inBrake/toTargetspeed'
 * '<S305>' : 'VCM20/VCM main/VDC/Traction Control/TargetSpeed_inAcc/Compare To Zero'
 * '<S306>' : 'VCM20/VCM main/VDC/Traction Control/TargetSpeed_inAcc/Compare To Zero1'
 * '<S307>' : 'VCM20/VCM main/VDC/Traction Control/TargetSpeed_inAcc/Compare To Zero2'
 * '<S308>' : 'VCM20/VCM main/VDC/Traction Control/TargetSpeed_inAcc/Compare To Zero3'
 * '<S309>' : 'VCM20/VCM main/VDC/Traction Control/TargetSpeed_inAcc/Saturation Dynamic'
 * '<S310>' : 'VCM20/VCM main/VDC/Traction Control/TargetSpeed_inAcc/toTargetspeed'
 * '<S311>' : 'VCM20/VCM main/VDC/Traction Control/TractionControll1/Saturation Dynamic'
 * '<S312>' : 'VCM20/VCM main/VDC/Traction Control/TractionControll1/Saturation Dynamic1'
 * '<S313>' : 'VCM20/VCM main/VDC/Traction Control/TractionControll2/Saturation Dynamic'
 * '<S314>' : 'VCM20/VCM main/VDC/Traction Control/TractionControll2/Saturation Dynamic1'
 * '<S315>' : 'VCM20/VCM main/VDC/Traction Control/TractionControll3/Saturation Dynamic'
 * '<S316>' : 'VCM20/VCM main/VDC/Traction Control/TractionControll3/Saturation Dynamic1'
 * '<S317>' : 'VCM20/VCM main/VDC/Traction Control/TractionControll4/Saturation Dynamic'
 * '<S318>' : 'VCM20/VCM main/VDC/Traction Control/TractionControll4/Saturation Dynamic1'
 * '<S319>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate'
 * '<S320>' : 'VCM20/VCM main/VDC/VCM_Status/VCMStatus_bus'
 * '<S321>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate/Sliprate Calculator'
 * '<S322>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate/Sliprate Calculator/FL_Sliprate Caluclator'
 * '<S323>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate/Sliprate Calculator/FR_Sliprate Caluclator'
 * '<S324>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate/Sliprate Calculator/RL_Sliprate Caluclator'
 * '<S325>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate/Sliprate Calculator/RR_Sliprate Caluclator'
 * '<S326>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate/Sliprate Calculator/FL_Sliprate Caluclator/Avoid division by zero'
 * '<S327>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate/Sliprate Calculator/FL_Sliprate Caluclator/Avoid division by zero/Compare To Zero'
 * '<S328>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate/Sliprate Calculator/FR_Sliprate Caluclator/Avoid division by zero'
 * '<S329>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate/Sliprate Calculator/FR_Sliprate Caluclator/Avoid division by zero/Compare To Zero'
 * '<S330>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate/Sliprate Calculator/RL_Sliprate Caluclator/Avoid division by zero'
 * '<S331>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate/Sliprate Calculator/RL_Sliprate Caluclator/Avoid division by zero/Compare To Zero'
 * '<S332>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate/Sliprate Calculator/RR_Sliprate Caluclator/Avoid division by zero'
 * '<S333>' : 'VCM20/VCM main/VDC/VCM_Status/Sliprate/Sliprate Calculator/RR_Sliprate Caluclator/Avoid division by zero/Compare To Zero'
 * '<S334>' : 'VCM20/VCMOutput/BrakeLight'
 * '<S335>' : 'VCM20/VCMOutput/Steer Chart'
 * '<S336>' : 'VCM20/VCMOutput/to steer'
 * '<S337>' : 'VCM20/VCMOutput/to steer LCD'
 * '<S338>' : 'VCM20/sensors from cooling/Subsystem'
 * '<S339>' : 'VCM20/sensors from cooling/Subsystem1'
 * '<S340>' : 'VCM20/sensors from cooling/Subsystem2'
 * '<S341>' : 'VCM20/sensors from cooling/reikyaku_volt1 '
 * '<S342>' : 'VCM20/sensors from cooling/reikyaku_volt2 '
 * '<S343>' : 'VCM20/sensors from cooling/reikyaku_volt3 '
 * '<S344>' : 'VCM20/speedlimit/kph to motor rpm'
 */
#endif                                 /* RTW_HEADER_VCM20_h_ */

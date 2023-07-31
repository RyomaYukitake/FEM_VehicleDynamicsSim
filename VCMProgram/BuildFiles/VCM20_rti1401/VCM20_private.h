/*
 * VCM20_private.h
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

#ifndef RTW_HEADER_VCM20_private_h_
#define RTW_HEADER_VCM20_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"
#include "VCM20.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmIsMajorTimeStep
#define rtmIsMajorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
#define rtmIsMinorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
#define rtmSetTPtr(rtm, val)           ((rtm)->Timing.t = (val))
#endif

/* ...  variable for information on a CAN channel */
extern can_tp1_canChannel* can_type1_channel_M3_C1;

/* ...  variable for information on a CAN channel */
extern can_tp1_canChannel* can_type1_channel_M3_C2;

/* ... definition of message variable for the RTICAN blocks */
#define CANTP1_M3_NUMMSG               11

extern can_tp1_canMsg* can_type1_msg_M3[CANTP1_M3_NUMMSG];

/* ... variable for taskqueue error checking                  */
extern Int32 rtican_type1_tq_error[CAN_TYPE1_NUM_MODULES]
  [CAN_TYPE1_NUM_TASKQUEUES];

/* ...  variable for information on a CAN channel */
extern can_tp1_canChannel* can_type1_channel_M1_C1;

/* ...  variable for information on a CAN channel */
extern can_tp1_canChannel* can_type1_channel_M1_C2;

/* ... definition of message variable for the RTICAN blocks */
#define CANTP1_M1_NUMMSG               29

extern can_tp1_canMsg* can_type1_msg_M1[CANTP1_M1_NUMMSG];

/* ... variable for taskqueue error checking                  */
extern Int32 rtican_type1_tq_error[CAN_TYPE1_NUM_MODULES]
  [CAN_TYPE1_NUM_TASKQUEUES];

/* ...  variable for information on a CAN channel */
extern can_tp1_canChannel* can_type1_channel_M2_C1;

/* ...  variable for information on a CAN channel */
extern can_tp1_canChannel* can_type1_channel_M2_C2;

/* ... definition of message variable for the RTICAN blocks */
#define CANTP1_M2_NUMMSG               16

extern can_tp1_canMsg* can_type1_msg_M2[CANTP1_M2_NUMMSG];

/* ... variable for taskqueue error checking                  */
extern Int32 rtican_type1_tq_error[CAN_TYPE1_NUM_MODULES]
  [CAN_TYPE1_NUM_TASKQUEUES];

/* Declaration of user indices (CAN_Type1_M3) */
#define CANTP1_M3_C1_RX_STD_0X139      0
#define RX_C1_STD_0X139                0
#undef RX_C1_STD_0X139
#define CANTP1_M3_C1_RX_STD_0X121      1
#define RX_C1_STD_0X121                1
#undef RX_C1_STD_0X121
#define CANTP1_M3_C1_RX_STD_0X122      2
#define RX_C1_STD_0X122                2
#undef RX_C1_STD_0X122
#define CANTP1_M3_C1_RX_STD_0X220      3
#define RX_C1_STD_0X220                3
#undef RX_C1_STD_0X220
#define CANTP1_M3_C1_RX_STD_0X134      4
#define RX_C1_STD_0X134                4
#undef RX_C1_STD_0X134
#define CANTP1_M3_C1_RX_STD_0X102      5
#define RX_C1_STD_0X102                5
#undef RX_C1_STD_0X102
#define CANTP1_M3_C1_RX_STD_0X151      6
#define RX_C1_STD_0X151                6
#undef RX_C1_STD_0X151
#define CANTP1_M3_C1_RX_STD_0X171      7
#define RX_C1_STD_0X171                7
#undef RX_C1_STD_0X171
#define CANTP1_M3_C1_RX_STD_0X100      8
#define RX_C1_STD_0X100                8
#undef RX_C1_STD_0X100
#define CANTP1_M3_C2_RX_STD_0X202      9
#define RX_C2_STD_0X202                9
#undef RX_C2_STD_0X202
#define CANTP1_M3_C2_RX_STD_0X201      10
#define RX_C2_STD_0X201                10
#undef RX_C2_STD_0X201

/* predefine needed TX-definition code to support TX-Custom code */
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M3_C1_STD;
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M3_C1_XTD;

/* predefine needed TX-definition code to support TX-Custom code */
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M3_C2_STD;
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M3_C2_XTD;

/* Declaration of user indices (CAN_Type1_M1) */
#define CANTP1_M1_C1_TX_STD_0X210      0
#define TX_C1_STD_0X210                0
#undef TX_C1_STD_0X210
#define CANTP1_M1_C1_TX_STD_0X211      1
#define TX_C1_STD_0X211                1
#undef TX_C1_STD_0X211
#define CANTP1_M1_C1_RX_STD_0X200      2
#define RX_C1_STD_0X200                2
#undef RX_C1_STD_0X200
#define CANTP1_M1_C1_RX_STD_0X401      3
#define RX_C1_STD_0X401                3
#undef RX_C1_STD_0X401
#define CANTP1_M1_C2_TX_STD_0X504      4
#define TX_C2_STD_0X504                4
#undef TX_C2_STD_0X504
#define CANTP1_M1_C2_TX_STD_0X521      5
#define TX_C2_STD_0X521                5
#undef TX_C2_STD_0X521
#define CANTP1_M1_C2_TX_STD_0X522      6
#define TX_C2_STD_0X522                6
#undef TX_C2_STD_0X522
#define CANTP1_M1_C2_TX_STD_0X701      7
#define TX_C2_STD_0X701                7
#undef TX_C2_STD_0X701
#define CANTP1_M1_C2_TX_STD_0X501      8
#define TX_C2_STD_0X501                8
#undef TX_C2_STD_0X501
#define CANTP1_M1_C2_TX_STD_0X502      9
#define TX_C2_STD_0X502                9
#undef TX_C2_STD_0X502
#define CANTP1_M1_C2_TX_STD_0X505      10
#define TX_C2_STD_0X505                10
#undef TX_C2_STD_0X505
#define CANTP1_M1_C2_TX_STD_0X506      11
#define TX_C2_STD_0X506                11
#undef TX_C2_STD_0X506
#define CANTP1_M1_C2_TX_STD_0X507      12
#define TX_C2_STD_0X507                12
#undef TX_C2_STD_0X507
#define CANTP1_M1_C2_TX_STD_0X508      13
#define TX_C2_STD_0X508                13
#undef TX_C2_STD_0X508
#define CANTP1_M1_C2_TX_STD_0X509      14
#define TX_C2_STD_0X509                14
#undef TX_C2_STD_0X509
#define CANTP1_M1_C2_TX_STD_0X50A      15
#define TX_C2_STD_0X50A                15
#undef TX_C2_STD_0X50A
#define CANTP1_M1_C2_TX_STD_0X50B      16
#define TX_C2_STD_0X50B                16
#undef TX_C2_STD_0X50B
#define CANTP1_M1_C2_TX_STD_0X50C      17
#define TX_C2_STD_0X50C                17
#undef TX_C2_STD_0X50C
#define CANTP1_M1_C2_TX_STD_0X50D      18
#define TX_C2_STD_0X50D                18
#undef TX_C2_STD_0X50D
#define CANTP1_M1_C2_TX_STD_0X50E      19
#define TX_C2_STD_0X50E                19
#undef TX_C2_STD_0X50E
#define CANTP1_M1_C2_TX_STD_0X50F      20
#define TX_C2_STD_0X50F                20
#undef TX_C2_STD_0X50F
#define CANTP1_M1_C2_TX_STD_0X510      21
#define TX_C2_STD_0X510                21
#undef TX_C2_STD_0X510
#define CANTP1_M1_C2_TX_STD_0X503      22
#define TX_C2_STD_0X503                22
#undef TX_C2_STD_0X503
#define CANTP1_M1_C2_TX_STD_0X705      23
#define TX_C2_STD_0X705                23
#undef TX_C2_STD_0X705
#define CANTP1_M1_C2_TX_STD_0X702      24
#define TX_C2_STD_0X702                24
#undef TX_C2_STD_0X702
#define CANTP1_M1_C2_TX_STD_0X703      25
#define TX_C2_STD_0X703                25
#undef TX_C2_STD_0X703
#define CANTP1_M1_C2_TX_STD_0X704      26
#define TX_C2_STD_0X704                26
#undef TX_C2_STD_0X704
#define CANTP1_M1_C2_RX_STD_0X202      27
#define RX_C2_STD_0X202                27
#undef RX_C2_STD_0X202
#define CANTP1_M1_C2_RX_STD_0X201      28
#define RX_C2_STD_0X201                28
#undef RX_C2_STD_0X201

/* predefine needed TX-definition code to support TX-Custom code */
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M1_C1_STD;
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M1_C1_XTD;

/* predefine needed TX-definition code to support TX-Custom code */
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M1_C2_STD;
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M1_C2_XTD;

/* Declaration of user indices (CAN_Type1_M2) */
#define CANTP1_M2_C1_TX_STD_0X184      0
#define TX_C1_STD_0X184                0
#undef TX_C1_STD_0X184
#define CANTP1_M2_C1_TX_STD_0X189      1
#define TX_C1_STD_0X189                1
#undef TX_C1_STD_0X189
#define CANTP1_M2_C1_RX_STD_0X283      2
#define RX_C1_STD_0X283                2
#undef RX_C1_STD_0X283
#define CANTP1_M2_C1_RX_STD_0X288      3
#define RX_C1_STD_0X288                3
#undef RX_C1_STD_0X288
#define CANTP1_M2_C1_RX_STD_0X285      4
#define RX_C1_STD_0X285                4
#undef RX_C1_STD_0X285
#define CANTP1_M2_C1_RX_STD_0X28A      5
#define RX_C1_STD_0X28A                5
#undef RX_C1_STD_0X28A
#define CANTP1_M2_C1_RX_STD_0X301      6
#define RX_C1_STD_0X301                6
#undef RX_C1_STD_0X301
#define CANTP1_M2_C1_RX_STD_0X306      7
#define RX_C1_STD_0X306                7
#undef RX_C1_STD_0X306
#define CANTP1_M2_C2_TX_STD_0X185      8
#define TX_C2_STD_0X185                8
#undef TX_C2_STD_0X185
#define CANTP1_M2_C2_TX_STD_0X188      9
#define TX_C2_STD_0X188                9
#undef TX_C2_STD_0X188
#define CANTP1_M2_C2_RX_STD_0X287      10
#define RX_C2_STD_0X287                10
#undef RX_C2_STD_0X287
#define CANTP1_M2_C2_RX_STD_0X284      11
#define RX_C2_STD_0X284                11
#undef RX_C2_STD_0X284
#define CANTP1_M2_C2_RX_STD_0X289      12
#define RX_C2_STD_0X289                12
#undef RX_C2_STD_0X289
#define CANTP1_M2_C2_RX_STD_0X286      13
#define RX_C2_STD_0X286                13
#undef RX_C2_STD_0X286
#define CANTP1_M2_C2_RX_STD_0X305      14
#define RX_C2_STD_0X305                14
#undef RX_C2_STD_0X305
#define CANTP1_M2_C2_RX_STD_0X302      15
#define RX_C2_STD_0X302                15
#undef RX_C2_STD_0X302

/* predefine needed TX-definition code to support TX-Custom code */
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M2_C1_STD;
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M2_C1_XTD;

/* predefine needed TX-definition code to support TX-Custom code */
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M2_C2_STD;
extern can_tp1_canMsg* CANTP1_TX_SPMSG_M2_C2_XTD;
extern real_T rt_roundd_snf(real_T u);
extern uint32_T MWDSP_EPH_R_D(real_T evt, uint32_T *sta);
extern uint32_T MWDSP_EPH_R_B(boolean_T evt, uint32_T *sta);
real_T look_SplNBinSZcd(uint32_T numDims, const real_T* u, const
  rt_LUTSplineWork * const SWork);
void rt_Spline2Derivd(const real_T *x, const real_T *y, uint32_T n, real_T *u,
                      real_T *y2);
real_T intrp_NSplcd(uint32_T numDims, const rt_LUTSplineWork * const splWork,
                    uint32_T extrapMethod);
extern real_T look1_binlxpw(real_T u0, const real_T bp0[], const real_T table[],
  uint32_T maxIndex);
extern real_T look2_binlxpw(real_T u0, real_T u1, const real_T bp0[], const
  real_T bp1[], const real_T table[], const uint32_T maxIndex[], uint32_T stride);
extern uint32_T plook_binx(real_T u, const real_T bp[], uint32_T maxIndex,
  real_T *fraction);
extern uint32_T binsearch_u32d(real_T u, const real_T bp[], uint32_T startIndex,
  uint32_T maxIndex);
extern void VCM20_MovingAverage_Init(DW_MovingAverage_VCM20_T *localDW);
extern void VCM20_MovingAverage_Start(DW_MovingAverage_VCM20_T *localDW);
extern void VCM20_MovingAverage(real_T rtu_0, B_MovingAverage_VCM20_T *localB,
  DW_MovingAverage_VCM20_T *localDW);
extern void VCM20_MovingAverage_Term(DW_MovingAverage_VCM20_T *localDW);

/* private model entry point functions */
extern void VCM20_derivatives(void);

#endif                                 /* RTW_HEADER_VCM20_private_h_ */

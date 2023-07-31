/***************************************************************************

   Source file VCM20_trc_ptr.c:

   Definition of function that initializes the global TRC pointers

   RTI1401 7.11 (02-Nov-2018)
   Mon Jul 31 14:00:26 2023

   Copyright 2023, dSPACE GmbH. All rights reserved.

 *****************************************************************************/

/* Include header file. */
#include "VCM20_trc_ptr.h"
#include "VCM20.h"
#include "VCM20_private.h"

/* Compiler options to turn off optimization. */
#if !defined(DS_OPTIMIZE_INIT_TRC_POINTERS)
#ifdef _MCCPPC

#pragma options -nOt -nOr -nOi -nOx

#endif

#ifdef __GNUC__

#pragma GCC optimize ("O0")

#endif

#ifdef _MSC_VER

#pragma optimize ("", off)

#endif
#endif

/* Definition of Global pointers to data type transitions (for TRC-file access) */
volatile real_T *p_0_VCM20_real_T_0 = NULL;
volatile real32_T *p_0_VCM20_real32_T_1 = NULL;
volatile uint8_T *p_0_VCM20_uint8_T_2 = NULL;
volatile boolean_T *p_0_VCM20_boolean_T_3 = NULL;
volatile real_T *p_0_VCM20_real_T_4 = NULL;
volatile real_T *p_0_VCM20_real_T_5 = NULL;
volatile real_T *p_0_VCM20_real_T_6 = NULL;
volatile real_T *p_0_VCM20_real_T_7 = NULL;
volatile real_T *p_0_VCM20_real_T_8 = NULL;
volatile real_T *p_1_VCM20_real_T_0 = NULL;
volatile real32_T *p_1_VCM20_real32_T_1 = NULL;
volatile uint8_T *p_1_VCM20_uint8_T_2 = NULL;
volatile real_T *p_1_VCM20_real_T_3 = NULL;
volatile uint32_T *p_1_VCM20_uint32_T_4 = NULL;
volatile boolean_T *p_1_VCM20_boolean_T_5 = NULL;
volatile real_T *p_2_VCM20_real_T_1 = NULL;
volatile int32_T *p_2_VCM20_int32_T_3 = NULL;
volatile uint32_T *p_2_VCM20_uint32_T_4 = NULL;
volatile int_T *p_2_VCM20_int_T_5 = NULL;
volatile boolean_T *p_2_VCM20_boolean_T_6 = NULL;
volatile uint8_T *p_2_VCM20_uint8_T_7 = NULL;
volatile boolean_T *p_2_VCM20_boolean_T_8 = NULL;
volatile boolean_T *p_2_VCM20_boolean_T_10 = NULL;
volatile boolean_T *p_2_VCM20_boolean_T_12 = NULL;
volatile boolean_T *p_2_VCM20_boolean_T_14 = NULL;
volatile boolean_T *p_2_VCM20_boolean_T_16 = NULL;
volatile boolean_T *p_2_VCM20_boolean_T_18 = NULL;
volatile real_T *p_3_VCM20_real_T_0 = NULL;

/*
 *  Declare the functions, that initially assign TRC pointers
 */
static void rti_init_trc_pointers_0(void);

/* Global pointers to data type transitions are separated in different functions to avoid overloading */
static void rti_init_trc_pointers_0(void)
{
  p_0_VCM20_real_T_0 = &VCM20_B.Saturation;
  p_0_VCM20_real32_T_1 = &VCM20_B.DataTypeConversion_o;
  p_0_VCM20_uint8_T_2 = &VCM20_B.LCDtext[0];
  p_0_VCM20_boolean_T_3 = &VCM20_B.LowerRelop1;
  p_0_VCM20_real_T_4 = &VCM20_B.MovingAverage_pnaev.MovingAverage;
  p_0_VCM20_real_T_5 = &VCM20_B.MovingAverage_pnae.MovingAverage;
  p_0_VCM20_real_T_6 = &VCM20_B.MovingAverage_pna.MovingAverage;
  p_0_VCM20_real_T_7 = &VCM20_B.MovingAverage_pn.MovingAverage;
  p_0_VCM20_real_T_8 = &VCM20_B.MovingAverage_p.MovingAverage;
  p_1_VCM20_real_T_0 = &VCM20_P.Iwf;
  p_1_VCM20_real32_T_1 = &VCM20_P.CompareToConstant_const_k;
  p_1_VCM20_uint8_T_2 = &VCM20_P.Counter_HitValue;
  p_1_VCM20_real_T_3 = &VCM20_P.LeftMax_deg_Gain;
  p_1_VCM20_uint32_T_4 = &VCM20_P.uDLookupTable_maxIndex[0];
  p_1_VCM20_boolean_T_5 = &VCM20_P.Delay_InitialCondition_ge;
  p_2_VCM20_real_T_1 = &VCM20_DW.Delay_DSTATE;
  p_2_VCM20_int32_T_3 = &VCM20_DW.clockTickCounter;
  p_2_VCM20_uint32_T_4 = &VCM20_DW.Counter_ClkEphState;
  p_2_VCM20_int_T_5 = &VCM20_DW.SFunction1_IWORK[0];
  p_2_VCM20_boolean_T_6 = &VCM20_DW.Delay1_DSTATE_c;
  p_2_VCM20_uint8_T_7 = &VCM20_DW.Counter_Count;
  p_2_VCM20_boolean_T_8 = &VCM20_DW.objisempty;
  p_2_VCM20_boolean_T_10 = &VCM20_DW.MovingAverage_pnaev.objisempty;
  p_2_VCM20_boolean_T_12 = &VCM20_DW.MovingAverage_pnae.objisempty;
  p_2_VCM20_boolean_T_14 = &VCM20_DW.MovingAverage_pna.objisempty;
  p_2_VCM20_boolean_T_16 = &VCM20_DW.MovingAverage_pn.objisempty;
  p_2_VCM20_boolean_T_18 = &VCM20_DW.MovingAverage_p.objisempty;
  p_3_VCM20_real_T_0 = &VCM20_X.TransferFcn_CSTATE;
}

void VCM20_rti_init_trc_pointers(void)
{
  rti_init_trc_pointers_0();
}

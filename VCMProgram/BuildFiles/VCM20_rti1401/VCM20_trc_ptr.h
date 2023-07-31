/*********************** dSPACE target specific file *************************

   Header file VCM20_trc_ptr.h:

   Declaration of function that initializes the global TRC pointers

   RTI1401 7.11 (02-Nov-2018)
   Mon Jul 31 14:00:26 2023

   Copyright 2023, dSPACE GmbH. All rights reserved.

 *****************************************************************************/
#ifndef RTI_HEADER_VCM20_trc_ptr_h_
#define RTI_HEADER_VCM20_trc_ptr_h_

/* Include the model header file. */
#include "VCM20.h"
#include "VCM20_private.h"
#ifdef EXTERN_C
#undef EXTERN_C
#endif

#ifdef __cplusplus
#define EXTERN_C                       extern "C"
#else
#define EXTERN_C                       extern
#endif

/*
 *  Declare the global TRC pointers
 */
EXTERN_C volatile real_T *p_0_VCM20_real_T_0;
EXTERN_C volatile real32_T *p_0_VCM20_real32_T_1;
EXTERN_C volatile uint8_T *p_0_VCM20_uint8_T_2;
EXTERN_C volatile boolean_T *p_0_VCM20_boolean_T_3;
EXTERN_C volatile real_T *p_0_VCM20_real_T_4;
EXTERN_C volatile real_T *p_0_VCM20_real_T_5;
EXTERN_C volatile real_T *p_0_VCM20_real_T_6;
EXTERN_C volatile real_T *p_0_VCM20_real_T_7;
EXTERN_C volatile real_T *p_0_VCM20_real_T_8;
EXTERN_C volatile real_T *p_1_VCM20_real_T_0;
EXTERN_C volatile real32_T *p_1_VCM20_real32_T_1;
EXTERN_C volatile uint8_T *p_1_VCM20_uint8_T_2;
EXTERN_C volatile real_T *p_1_VCM20_real_T_3;
EXTERN_C volatile uint32_T *p_1_VCM20_uint32_T_4;
EXTERN_C volatile boolean_T *p_1_VCM20_boolean_T_5;
EXTERN_C volatile real_T *p_2_VCM20_real_T_1;
EXTERN_C volatile int32_T *p_2_VCM20_int32_T_3;
EXTERN_C volatile uint32_T *p_2_VCM20_uint32_T_4;
EXTERN_C volatile int_T *p_2_VCM20_int_T_5;
EXTERN_C volatile boolean_T *p_2_VCM20_boolean_T_6;
EXTERN_C volatile uint8_T *p_2_VCM20_uint8_T_7;
EXTERN_C volatile boolean_T *p_2_VCM20_boolean_T_8;
EXTERN_C volatile boolean_T *p_2_VCM20_boolean_T_10;
EXTERN_C volatile boolean_T *p_2_VCM20_boolean_T_12;
EXTERN_C volatile boolean_T *p_2_VCM20_boolean_T_14;
EXTERN_C volatile boolean_T *p_2_VCM20_boolean_T_16;
EXTERN_C volatile boolean_T *p_2_VCM20_boolean_T_18;
EXTERN_C volatile real_T *p_3_VCM20_real_T_0;

/*
 *  Declare the general function for TRC pointer initialization
 */
EXTERN_C void VCM20_rti_init_trc_pointers(void);

#endif                                 /* RTI_HEADER_VCM20_trc_ptr_h_ */

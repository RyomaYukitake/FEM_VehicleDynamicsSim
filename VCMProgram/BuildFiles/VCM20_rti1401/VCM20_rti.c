/*********************** dSPACE target specific file *************************

   Include file VCM20_rti.c:

   Definition of functions and variables for the system I/O and for
   the hardware and software interrupts used.

   RTI1401 7.11 (02-Nov-2018)
   Mon Jul 31 14:28:36 2023

   Copyright 2023, dSPACE GmbH. All rights reserved.

 *****************************************************************************/

#if !(defined(__RTI_SIMENGINE__) || defined(RTIMP_FRAME))
# error This file may be included only by the RTI(-MP) simulation engine.
#endif

/* Include the model header file. */
#include "VCM20.h"
#include "VCM20_private.h"

/* Defines for block output and parameter structure existence */
#define RTI_rtB_STRUCTURE_EXISTS       1
#define RTI_rtP_STRUCTURE_EXISTS       1
#define RTB_STRUCTURE_NAME             VCM20_B
#define RTP_STRUCTURE_NAME             VCM20_P

/* dSPACE generated includes for header files */
#include <brtenv.h>
#include <rtkernel.h>
#include <rti_assert.h>
#include <rtidefineddatatypes.h>
#include <rtican_ds1401.h>
#ifndef dsRtmGetNumSampleTimes
# define dsRtmGetNumSampleTimes(rtm)   4
#endif

#ifndef dsRtmGetTPtr
# define dsRtmGetTPtr(rtm)             ((rtm)->Timing.t)
#endif

#ifndef dsRtmSetTaskTime
# define dsRtmSetTaskTime(rtm, sti, val) (dsRtmGetTPtr((rtm))[sti] = (val))
#endif

/****** Definitions: task functions for timer tasks *********************/

/* Timer Task 1. (Base rate). */
static void rti_TIMERA(rtk_p_task_control_block task)
{
  /* Task entry code BEGIN */
  /* -- None. -- */
  /* Task entry code END */

  /* Task code. */
  baseRateService(task);

  /* Task exit code BEGIN */
  /* -- None. -- */
  /* Task exit code END */
}

/* ===== Declarations of RTI blocks ======================================== */

/* dSPACE I/O Board DS1_RTICAN #3 */

/* ...  definition of channel struct */
can_tp1_canChannel* can_type1_channel_M3_C1;

/* ...  definition of channel struct */
can_tp1_canChannel* can_type1_channel_M3_C2;

/* declare pointers to CAN message structures for support of TX-Custom code */
can_tp1_canMsg* CANTP1_TX_SPMSG_M3_C1_STD;
can_tp1_canMsg* CANTP1_TX_SPMSG_M3_C1_XTD;

/* declare pointers to CAN message structures for support of TX-Custom code */
can_tp1_canMsg* CANTP1_TX_SPMSG_M3_C2_STD;
can_tp1_canMsg* CANTP1_TX_SPMSG_M3_C2_XTD;

/* ... definition of message variable for the RTICAN blocks */
can_tp1_canMsg* can_type1_msg_M3[CANTP1_M3_NUMMSG];

/* dSPACE I/O Board DS1_RTICAN #1 */

/* ...  definition of channel struct */
can_tp1_canChannel* can_type1_channel_M1_C1;

/* ...  definition of channel struct */
can_tp1_canChannel* can_type1_channel_M1_C2;

/* declare pointers to CAN message structures for support of TX-Custom code */
can_tp1_canMsg* CANTP1_TX_SPMSG_M1_C1_STD;
can_tp1_canMsg* CANTP1_TX_SPMSG_M1_C1_XTD;

/* declare pointers to CAN message structures for support of TX-Custom code */
can_tp1_canMsg* CANTP1_TX_SPMSG_M1_C2_STD;
can_tp1_canMsg* CANTP1_TX_SPMSG_M1_C2_XTD;

/* ... definition of message variable for the RTICAN blocks */
can_tp1_canMsg* can_type1_msg_M1[CANTP1_M1_NUMMSG];

/* dSPACE I/O Board DS1_RTICAN #2 */

/* ...  definition of channel struct */
can_tp1_canChannel* can_type1_channel_M2_C1;

/* ...  definition of channel struct */
can_tp1_canChannel* can_type1_channel_M2_C2;

/* declare pointers to CAN message structures for support of TX-Custom code */
can_tp1_canMsg* CANTP1_TX_SPMSG_M2_C1_STD;
can_tp1_canMsg* CANTP1_TX_SPMSG_M2_C1_XTD;

/* declare pointers to CAN message structures for support of TX-Custom code */
can_tp1_canMsg* CANTP1_TX_SPMSG_M2_C2_STD;
can_tp1_canMsg* CANTP1_TX_SPMSG_M2_C2_XTD;

/* ... definition of message variable for the RTICAN blocks */
can_tp1_canMsg* can_type1_msg_M2[CANTP1_M2_NUMMSG];

/* dSPACE I/O Board RTICAN_GLOBAL #0 */

/* ===== Definition of interface functions for simulation engine =========== */
#if GRTINTERFACE == 1
#ifdef MULTITASKING
# define dsIsSampleHit(RTM,sti)        rtmGetSampleHitPtr(RTM)[sti]
#else
# define dsIsSampleHit(RTM,sti)        rtmIsSampleHit(RTM,sti,0)
#endif

#else
#ifndef rtmStepTask
# define rtmStepTask(rtm, idx)         ((rtm)->Timing.TaskCounters.TID[(idx)] == 0)
#endif

# define dsIsSampleHit(RTM,sti)        rtmStepTask(RTM, sti)
#endif

#undef __INLINE
#if defined(_INLINE)
# define __INLINE                      static inline
#else
# define __INLINE                      static
#endif

/*Define additional variables*/
static time_T dsTFinal = -1.0;

#define dsGetTFinal(rtm)               (dsTFinal)

static time_T dsStepSize = 0.001;

# define dsGetStepSize(rtm)            (dsStepSize)

static void rti_mdl_initialize_host_services(void)
{
  ts_timestamp_type ts = { 0, 0 };

  host_service(1, &ts);
  DsDaq_Init(0, 32, 1);
}

static void rti_mdl_initialize_io_boards(void)
{
  /* Registering of RTI products and modules at VCM */
  {
    vcm_module_register(VCM_MID_RTI1401, (void *) 0,
                        VCM_TXT_RTI1401, 7, 11, 0,
                        VCM_VERSION_RELEASE, 0, 0, 0, VCM_CTRL_NO_ST);

    {
      vcm_module_descriptor_type* msg_mod_ptr;
      msg_mod_ptr = vcm_module_register(VCM_MID_MATLAB, (void *) 0,
        VCM_TXT_MATLAB, 9, 11, 0,
        VCM_VERSION_RELEASE, 0, 0, 0, VCM_CTRL_NO_ST);
      vcm_module_register(VCM_MID_SIMULINK, msg_mod_ptr,
                          VCM_TXT_SIMULINK, 10, 4, 0,
                          VCM_VERSION_RELEASE, 0, 0, 0, VCM_CTRL_NO_ST);
      vcm_module_register(VCM_MID_RTW, msg_mod_ptr,
                          VCM_TXT_RTW, 9, 6, 0,
                          VCM_VERSION_RELEASE, 0, 0, 0, VCM_CTRL_NO_ST);
      vcm_module_register(VCM_MID_STATEFLOW, msg_mod_ptr,
                          VCM_TXT_STATEFLOW, 10, 5, 0,
                          VCM_VERSION_RELEASE, 0, 0, 0, VCM_CTRL_NO_ST);
      vcm_module_register(VCM_MID_STATEFLOW_CODER, msg_mod_ptr,
                          VCM_TXT_STATEFLOW_CODER, 9, 6, 0,
                          VCM_VERSION_RELEASE, 0, 0, 0, VCM_CTRL_NO_ST);
    }

    vcm_module_register(VCM_MID_RTICAN, (void *) 0,
                        VCM_TXT_RTICAN, 3, 4, 7,
                        VCM_VERSION_RELEASE, 0, 0, 0, VCM_CTRL_NO_ST);
  }

  /* dSPACE I/O Board DS1401STDADCT4 #1 */
  adc_tp4_init(ADC_TP4_1_MODULE_ADDR);

  /* dSPACE I/O Board DS1401STDADCT4 #1 Unit:ADC */
  adc_tp4_single_init2(ADC_TP4_1_MODULE_ADDR,
                       ADC_TP4_CH4,
                       ADC_TP4_INPUT_VOLTAGE_RANGE_2,
                       ADC_TP4_TRIG_SW,
                       0.0);
  adc_tp4_single_init2(ADC_TP4_1_MODULE_ADDR,
                       ADC_TP4_CH9,
                       ADC_TP4_INPUT_VOLTAGE_RANGE_2,
                       ADC_TP4_TRIG_SW,
                       0.0);
  adc_tp4_single_init2(ADC_TP4_1_MODULE_ADDR,
                       ADC_TP4_CH10,
                       ADC_TP4_INPUT_VOLTAGE_RANGE_2,
                       ADC_TP4_TRIG_SW,
                       0.0);
  adc_tp4_single_init2(ADC_TP4_1_MODULE_ADDR,
                       ADC_TP4_CH8,
                       ADC_TP4_INPUT_VOLTAGE_RANGE_2,
                       ADC_TP4_TRIG_SW,
                       0.0);

  /* dSPACE I/O Board DS1401STDDIOT4 #1 */
  /* main init for DIO TP4 module 1 */
  dio_tp4_init(DIO_TP4_1_MODULE_ADDR);

  /* dSPACE I/O Board DS1401STDDIOT4 #1 Unit:BITINCH */
  dio_tp4_digin_init(DIO_TP4_1_MODULE_ADDR, 1, DIO_TP4_MASK_CH5);
  dio_tp4_digin_init(DIO_TP4_1_MODULE_ADDR, 1, DIO_TP4_MASK_CH15);
  dio_tp4_digin_init(DIO_TP4_1_MODULE_ADDR, 1, DIO_TP4_MASK_CH10);

  /* dSPACE I/O Board DS1401STDDIOT4 #1 Unit:DIGOUT */
  /* setup digital output mode (supply rails) for channel 4 on port 1 */
  dio_tp4_digout_init(DIO_TP4_1_MODULE_ADDR, 1, DIO_TP4_MASK_CH4, 1, 1);

  /* dSPACE I/O Board DS1_RTICAN #3 */
  /* Initialization of DS1501 board */
  can_tp1_communication_init(can_tp1_address_table[2].module_addr,
    CAN_TP1_INT_ENABLE);

  /* dSPACE RTICAN MASTER SETUP Block */
  /* ... Initialize the CAN communication: 500 kbit/s */
  can_type1_channel_M3_C1 = can_tp1_channel_init(can_tp1_address_table[2].
    module_addr, 0, (500 * 1000), CAN_TP1_STD, CAN_TP1_NO_SUBINT);
  can_tp1_channel_termination_set(can_type1_channel_M3_C1,
    CAN_TP1_TERMINATION_ON);

  /* ... Initialize TX message structs for custom code  */
  CANTP1_TX_SPMSG_M3_C1_STD = can_tp1_msg_tx_register(can_type1_channel_M3_C1, 3,
    1050, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO | CAN_TP1_DELAYCOUNT_INFO |
    CAN_TP1_MSG_INFO, CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG,
    CAN_TP1_TIMEOUT_NORMAL);
  CANTP1_TX_SPMSG_M3_C1_XTD = can_tp1_msg_tx_register(can_type1_channel_M3_C1, 3,
    1100, CAN_TP1_EXT, CAN_TP1_TIMECOUNT_INFO | CAN_TP1_DELAYCOUNT_INFO |
    CAN_TP1_MSG_INFO, CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG,
    CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN MASTER SETUP Block */
  /* ... Initialize the Partial Networking Settings */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:313 */
  /* ... Register message */
  can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X139] = can_tp1_msg_rx_register
    (can_type1_channel_M3_C1, 2, 313, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:289 */
  /* ... Register message */
  can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X121] = can_tp1_msg_rx_register
    (can_type1_channel_M3_C1, 1, 289, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:290 */
  /* ... Register message */
  can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X122] = can_tp1_msg_rx_register
    (can_type1_channel_M3_C1, 1, 290, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:544 */
  /* ... Register message */
  can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X220] = can_tp1_msg_rx_register
    (can_type1_channel_M3_C1, 5, 544, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:308 */
  /* ... Register message */
  can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X134] = can_tp1_msg_rx_register
    (can_type1_channel_M3_C1, 2, 308, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:258 */
  /* ... Register message */
  can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X102] = can_tp1_msg_rx_register
    (can_type1_channel_M3_C1, 0, 258, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:337 */
  /* ... Register message */
  can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X151] = can_tp1_msg_rx_register
    (can_type1_channel_M3_C1, 3, 337, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:369 */
  /* ... Register message */
  can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X171] = can_tp1_msg_rx_register
    (can_type1_channel_M3_C1, 3, 369, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:256 */
  /* ... Register message */
  can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X100] = can_tp1_msg_rx_register
    (can_type1_channel_M3_C1, 0, 256, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN MASTER SETUP Block */
  /* ... Initialize the CAN communication: 500 kbit/s */
  can_type1_channel_M3_C2 = can_tp1_channel_init(can_tp1_address_table[2].
    module_addr, 1, (500 * 1000), CAN_TP1_STD, CAN_TP1_NO_SUBINT);
  can_tp1_channel_termination_set(can_type1_channel_M3_C2,
    CAN_TP1_TERMINATION_ON);

  /* ... Initialize TX message structs for custom code  */
  CANTP1_TX_SPMSG_M3_C2_STD = can_tp1_msg_tx_register(can_type1_channel_M3_C2, 3,
    1050, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO | CAN_TP1_DELAYCOUNT_INFO |
    CAN_TP1_MSG_INFO, CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG,
    CAN_TP1_TIMEOUT_NORMAL);
  CANTP1_TX_SPMSG_M3_C2_XTD = can_tp1_msg_tx_register(can_type1_channel_M3_C2, 3,
    1100, CAN_TP1_EXT, CAN_TP1_TIMECOUNT_INFO | CAN_TP1_DELAYCOUNT_INFO |
    CAN_TP1_MSG_INFO, CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG,
    CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN MASTER SETUP Block */
  /* ... Initialize the Partial Networking Settings */

  /* dSPACE RTICAN RX Message Block: "Brakes" Id:514 */
  /* ... Register message */
  can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X202] = can_tp1_msg_rx_register
    (can_type1_channel_M3_C2, 4, 514, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "APPS_Steer" Id:513 */
  /* ... Register message */
  can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X201] = can_tp1_msg_rx_register
    (can_type1_channel_M3_C2, 4, 513, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_DATA_INFO|CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "APPS_Steer" Id:513 */
  VCM20_B.SFunction1_o5_l = 0;         /* processed - flag */

  /* dSPACE I/O Board DS1_RTICAN #1 */
  /* Initialization of DS1501 board */
  can_tp1_communication_init(can_tp1_address_table[0].module_addr,
    CAN_TP1_INT_ENABLE);

  /* dSPACE RTICAN MASTER SETUP Block */
  /* ... Initialize the CAN communication: 500 kbit/s */
  can_type1_channel_M1_C1 = can_tp1_channel_init(can_tp1_address_table[0].
    module_addr, 0, (500 * 1000), CAN_TP1_STD, CAN_TP1_NO_SUBINT);
  can_tp1_channel_termination_set(can_type1_channel_M1_C1,
    CAN_TP1_TERMINATION_OFF);

  /* ... Initialize TX message structs for custom code  */
  CANTP1_TX_SPMSG_M1_C1_STD = can_tp1_msg_tx_register(can_type1_channel_M1_C1, 3,
    1050, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO | CAN_TP1_DELAYCOUNT_INFO |
    CAN_TP1_MSG_INFO, CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG,
    CAN_TP1_TIMEOUT_NORMAL);
  CANTP1_TX_SPMSG_M1_C1_XTD = can_tp1_msg_tx_register(can_type1_channel_M1_C1, 3,
    1100, CAN_TP1_EXT, CAN_TP1_TIMECOUNT_INFO | CAN_TP1_DELAYCOUNT_INFO |
    CAN_TP1_MSG_INFO, CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG,
    CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN MASTER SETUP Block */
  /* ... Initialize the Partial Networking Settings */

  /* dSPACE RTICAN TX Message Block: "Steer LED" Id:528 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X210] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C1, 0, 528, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "Steer LCD" Id:529 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X211] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C1, 0, 529, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN RX Message Block: "Steer SW" Id:512 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X200] = can_tp1_msg_rx_register
    (can_type1_channel_M1_C1, 0, 512, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:1025 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X401] = can_tp1_msg_rx_register
    (can_type1_channel_M1_C1, 1, 1025, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN MASTER SETUP Block */
  /* ... Initialize the CAN communication: 500 kbit/s */
  can_type1_channel_M1_C2 = can_tp1_channel_init(can_tp1_address_table[0].
    module_addr, 1, (500 * 1000), CAN_TP1_STD, CAN_TP1_NO_SUBINT);
  can_tp1_channel_termination_set(can_type1_channel_M1_C2,
    CAN_TP1_TERMINATION_ON);

  /* ... Initialize TX message structs for custom code  */
  CANTP1_TX_SPMSG_M1_C2_STD = can_tp1_msg_tx_register(can_type1_channel_M1_C2, 3,
    1050, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO | CAN_TP1_DELAYCOUNT_INFO |
    CAN_TP1_MSG_INFO, CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG,
    CAN_TP1_TIMEOUT_NORMAL);
  CANTP1_TX_SPMSG_M1_C2_XTD = can_tp1_msg_tx_register(can_type1_channel_M1_C2, 3,
    1100, CAN_TP1_EXT, CAN_TP1_TIMECOUNT_INFO | CAN_TP1_DELAYCOUNT_INFO |
    CAN_TP1_MSG_INFO, CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG,
    CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN MASTER SETUP Block */
  /* ... Initialize the Partial Networking Settings */

  /* dSPACE RTICAN TX Message Block: "EVO4s004" Id:1284 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X504] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 1, 1284, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s004" Id:1313 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X521] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 4, 1313, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s004" Id:1314 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X522] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 4, 1314, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s004" Id:1793 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X701] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 4, 1793, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1281 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X501] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 1, 1281, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s002" Id:1282 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X502] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 1, 1282, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s005" Id:1285 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X505] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 2, 1285, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s006" Id:1286 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X506] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 2, 1286, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s007" Id:1287 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X507] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 2, 1287, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s008" Id:1288 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X508] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 2, 1288, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s009" Id:1289 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X509] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 2, 1289, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s00A" Id:1290 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50A] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 3, 1290, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1291 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50B] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 3, 1291, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1292 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50C] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 3, 1292, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1293 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50D] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 3, 1293, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1294 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50E] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 3, 1294, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1295 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50F] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 4, 1295, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:1296 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X510] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 4, 1296, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s003" Id:1283 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X503] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 1, 1283, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s003" Id:1797 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X705] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 5, 1797, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s005" Id:1794 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X702] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 5, 1794, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s005" Id:1795 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X703] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 5, 1795, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "EVO4s005" Id:1796 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X704] = can_tp1_msg_tx_register
    (can_type1_channel_M1_C2, 5, 1796, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:514 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X202] = can_tp1_msg_rx_register
    (can_type1_channel_M1_C2, 0, 514, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:513 */
  /* ... Register message */
  can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X201] = can_tp1_msg_rx_register
    (can_type1_channel_M1_C2, 0, 513, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE I/O Board DS1_RTICAN #2 */
  /* Initialization of DS1501 board */
  can_tp1_communication_init(can_tp1_address_table[1].module_addr,
    CAN_TP1_INT_ENABLE);

  /* dSPACE RTICAN MASTER SETUP Block */
  /* ... Initialize the CAN communication: 500 kbit/s */
  can_type1_channel_M2_C1 = can_tp1_channel_init(can_tp1_address_table[1].
    module_addr, 0, (500 * 1000), CAN_TP1_STD, CAN_TP1_NO_SUBINT);
  can_tp1_channel_termination_set(can_type1_channel_M2_C1,
    CAN_TP1_TERMINATION_ON);

  /* ... Initialize TX message structs for custom code  */
  CANTP1_TX_SPMSG_M2_C1_STD = can_tp1_msg_tx_register(can_type1_channel_M2_C1, 3,
    1050, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO | CAN_TP1_DELAYCOUNT_INFO |
    CAN_TP1_MSG_INFO, CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG,
    CAN_TP1_TIMEOUT_NORMAL);
  CANTP1_TX_SPMSG_M2_C1_XTD = can_tp1_msg_tx_register(can_type1_channel_M2_C1, 3,
    1100, CAN_TP1_EXT, CAN_TP1_TIMECOUNT_INFO | CAN_TP1_DELAYCOUNT_INFO |
    CAN_TP1_MSG_INFO, CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG,
    CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN MASTER SETUP Block */
  /* ... Initialize the Partial Networking Settings */

  /* dSPACE RTICAN TX Message Block: "VCM_11" Id:388 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X184] = can_tp1_msg_tx_register
    (can_type1_channel_M2_C1, 0, 388, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "VCM_41" Id:393 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X189] = can_tp1_msg_tx_register
    (can_type1_channel_M2_C1, 1, 393, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN RX Message Block: "MC_11" Id:643 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X283] = can_tp1_msg_rx_register
    (can_type1_channel_M2_C1, 1, 643, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "MC_41" Id:648 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X288] = can_tp1_msg_rx_register
    (can_type1_channel_M2_C1, 3, 648, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "MC_12" Id:645 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X285] = can_tp1_msg_rx_register
    (can_type1_channel_M2_C1, 2, 645, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "MC_42" Id:650 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X28A] = can_tp1_msg_rx_register
    (can_type1_channel_M2_C1, 4, 650, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "MC_13" Id:769 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X301] = can_tp1_msg_rx_register
    (can_type1_channel_M2_C1, 4, 769, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "MC_43" Id:774 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X306] = can_tp1_msg_rx_register
    (can_type1_channel_M2_C1, 5, 774, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN MASTER SETUP Block */
  /* ... Initialize the CAN communication: 500 kbit/s */
  can_type1_channel_M2_C2 = can_tp1_channel_init(can_tp1_address_table[1].
    module_addr, 1, (500 * 1000), CAN_TP1_STD, CAN_TP1_NO_SUBINT);
  can_tp1_channel_termination_set(can_type1_channel_M2_C2,
    CAN_TP1_TERMINATION_ON);

  /* ... Initialize TX message structs for custom code  */
  CANTP1_TX_SPMSG_M2_C2_STD = can_tp1_msg_tx_register(can_type1_channel_M2_C2, 3,
    1050, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO | CAN_TP1_DELAYCOUNT_INFO |
    CAN_TP1_MSG_INFO, CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG,
    CAN_TP1_TIMEOUT_NORMAL);
  CANTP1_TX_SPMSG_M2_C2_XTD = can_tp1_msg_tx_register(can_type1_channel_M2_C2, 3,
    1100, CAN_TP1_EXT, CAN_TP1_TIMECOUNT_INFO | CAN_TP1_DELAYCOUNT_INFO |
    CAN_TP1_MSG_INFO, CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG,
    CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN MASTER SETUP Block */
  /* ... Initialize the Partial Networking Settings */

  /* dSPACE RTICAN TX Message Block: "VCM_21" Id:389 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X185] = can_tp1_msg_tx_register
    (can_type1_channel_M2_C2, 0, 389, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN TX Message Block: "VCM01" Id:392 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X188] = can_tp1_msg_tx_register
    (can_type1_channel_M2_C2, 0, 392, CAN_TP1_STD, CAN_TP1_TIMECOUNT_INFO,
     CAN_TP1_NO_SUBINT, 0, CAN_TP1_TRIGGER_MSG, CAN_TP1_TIMEOUT_NORMAL);

  /* dSPACE RTICAN RX Message Block: "MC_31" Id:647 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X287] = can_tp1_msg_rx_register
    (can_type1_channel_M2_C2, 3, 647, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "MC_21" Id:644 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X284] = can_tp1_msg_rx_register
    (can_type1_channel_M2_C2, 1, 644, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "MC_32" Id:649 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X289] = can_tp1_msg_rx_register
    (can_type1_channel_M2_C2, 3, 649, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "MC_22" Id:646 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X286] = can_tp1_msg_rx_register
    (can_type1_channel_M2_C2, 2, 646, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "MC_33" Id:773 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X305] = can_tp1_msg_rx_register
    (can_type1_channel_M2_C2, 5, 773, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE RTICAN RX Message Block: "MC_23" Id:770 */
  /* ... Register message */
  can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X302] = can_tp1_msg_rx_register
    (can_type1_channel_M2_C2, 4, 770, CAN_TP1_STD, (CAN_TP1_DATA_INFO|
      CAN_TP1_TIMECOUNT_INFO), CAN_TP1_NO_SUBINT);

  /* dSPACE I/O Board RTICAN_GLOBAL #0 */
}

/* Function rti_mdl_slave_load() is empty */
#define rti_mdl_slave_load()

/* Function rti_mdl_rtk_initialize() is empty */
#define rti_mdl_rtk_initialize()

static void rti_mdl_initialize_io_units(void)
{
  /* dSPACE I/O Board DS1401STDDIOT4 #1 Unit:DIGOUT */
  /* re-adjust the supply rails for channel 4 on port 1 */
  dio_tp4_digout_init(DIO_TP4_1_MODULE_ADDR, 1, DIO_TP4_MASK_CH4, 1, 1);

  /* write initialization value to digital output channel 4 on port 1 */
  dio_tp4_digout_write(DIO_TP4_1_MODULE_ADDR , 1, DIO_TP4_MASK_CH4, (UInt16) 0);

  /* dSPACE I/O Board DS1401STDDIOT4 #1 Unit:DIGOUT Group:BITOUTCH */
  /* enable all digital output ports (no port wise Dis-Enabling supported) */
  /* called once per module for all dig output ports */
  dio_tp4_digout_mode_set(DIO_TP4_1_MODULE_ADDR, DIO_TP4_DIGOUT_ENABLE);

  /* dSPACE I/O Board DS1_RTICAN #3 */
  /* Start CAN controller */
  can_tp1_channel_start(can_type1_channel_M3_C1, CAN_TP1_INT_DISABLE);

  /* Start CAN controller */
  can_tp1_channel_start(can_type1_channel_M3_C2, CAN_TP1_INT_DISABLE);

  /* Set the type1CAN error level */
  rtican_type1_error_level = 0;

  /* ... Reset all taskqueue-specific error variables */
  rtican_type1_tq_err_all_reset(2);

  /* ... Clear all message data buffers */
  can_tp1_all_data_clear(can_tp1_address_table[2].module_addr);

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[2][2] = can_tp1_msg_wakeup
              (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X139])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[2][1] = can_tp1_msg_wakeup
              (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X121])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[2][1] = can_tp1_msg_wakeup
              (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X122])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[2][5] = can_tp1_msg_wakeup
              (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X220])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[2][2] = can_tp1_msg_wakeup
              (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X134])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[2][0] = can_tp1_msg_wakeup
              (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X102])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[2][3] = can_tp1_msg_wakeup
              (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X151])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[2][3] = can_tp1_msg_wakeup
              (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X171])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[2][0] = can_tp1_msg_wakeup
              (can_type1_msg_M3[CANTP1_M3_C1_RX_STD_0X100])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[2][4] = can_tp1_msg_wakeup
              (can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X202])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[2][4] = can_tp1_msg_wakeup
              (can_type1_msg_M3[CANTP1_M3_C2_RX_STD_0X201])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  /* dSPACE I/O Board DS1_RTICAN #1 */
  /* Start CAN controller */
  can_tp1_channel_start(can_type1_channel_M1_C1, CAN_TP1_INT_DISABLE);

  /* Start CAN controller */
  can_tp1_channel_start(can_type1_channel_M1_C2, CAN_TP1_INT_DISABLE);

  /* Set the type1CAN error level */
  rtican_type1_error_level = 0;

  /* ... Reset all taskqueue-specific error variables */
  rtican_type1_tq_err_all_reset(0);

  /* ... Clear all message data buffers */
  can_tp1_all_data_clear(can_tp1_address_table[0].module_addr);

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][0] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X210])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][0] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C1_TX_STD_0X211])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][0] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X200])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][1] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X401])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][1] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X504])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][4] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X521])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][4] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X522])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][4] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X701])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][1] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X501])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][1] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X502])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][2] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X505])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][2] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X506])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][2] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X507])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][2] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X508])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][2] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X509])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][3] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50A])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][3] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50B])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][3] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50C])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][3] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50D])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][3] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50E])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][4] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X50F])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][4] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X510])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][1] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X503])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][5] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X705])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][5] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X702])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][5] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X703])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][5] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X704])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][0] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X202])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[0][0] = can_tp1_msg_wakeup
              (can_type1_msg_M1[CANTP1_M1_C2_RX_STD_0X201])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  /* dSPACE I/O Board DS1_RTICAN #2 */
  /* Start CAN controller */
  can_tp1_channel_start(can_type1_channel_M2_C1, CAN_TP1_INT_DISABLE);

  /* Start CAN controller */
  can_tp1_channel_start(can_type1_channel_M2_C2, CAN_TP1_INT_DISABLE);

  /* Set the type1CAN error level */
  rtican_type1_error_level = 0;

  /* ... Reset all taskqueue-specific error variables */
  rtican_type1_tq_err_all_reset(1);

  /* ... Clear all message data buffers */
  can_tp1_all_data_clear(can_tp1_address_table[1].module_addr);

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][0] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X184])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][1] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C1_TX_STD_0X189])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][1] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X283])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][3] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X288])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][2] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X285])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][4] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X28A])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][4] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X301])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][5] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C1_RX_STD_0X306])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][0] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X185])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][0] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C2_TX_STD_0X188])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][3] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X287])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][1] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X284])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][3] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X289])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][2] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X286])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][5] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X305])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }

  {
    static UInt32 numInit = 0;
    if (numInit != 0) {
      /* ... Wake message up */
      while ((rtican_type1_tq_error[1][4] = can_tp1_msg_wakeup
              (can_type1_msg_M2[CANTP1_M2_C2_RX_STD_0X302])) ==
             DSMCOM_BUFFER_OVERFLOW) ;
    }

    ++numInit;
  }
}

/* Function rti_mdl_acknowledge_interrupts() is empty */
#define rti_mdl_acknowledge_interrupts()

/* Function rti_mdl_timetables_register() is empty */
#define rti_mdl_timetables_register()

/* Function rti_mdl_timesync_simstate() is empty */
#define rti_mdl_timesync_simstate()

static void rti_mdl_background(void)
{
  /* DsDaq background call */
  DsDaq_Background(0);

  /* dSPACE I/O Board DS1_RTICAN #3 */
  {
    real_T bg_code_exec_time;
    static real_T bg_code_last_exec_time = 0.0;
    bg_code_exec_time = RTLIB_TIC_READ();
    if ((bg_code_exec_time - bg_code_last_exec_time) > 0.25 ||
        (bg_code_exec_time - bg_code_last_exec_time) < 0) {
      /* ... Check taskqueue-specific error variables */
      rtican_type1_tq_err_all_chk(can_tp1_address_table[2].module_addr, 2);
      bg_code_last_exec_time = bg_code_exec_time;
    }
  }

  /* copy DPMEM - buffers in background */
  {
    /* call update function for CAN Tp1 CAN interface (module number: 3) */
    can_tp1_msg_copy_all_to_mem(can_tp1_address_table[2].module_addr);
  }

  /* dSPACE I/O Board DS1_RTICAN #1 */
  {
    real_T bg_code_exec_time;
    static real_T bg_code_last_exec_time = 0.0;
    bg_code_exec_time = RTLIB_TIC_READ();
    if ((bg_code_exec_time - bg_code_last_exec_time) > 0.25 ||
        (bg_code_exec_time - bg_code_last_exec_time) < 0) {
      /* ... Check taskqueue-specific error variables */
      rtican_type1_tq_err_all_chk(can_tp1_address_table[0].module_addr, 0);
      bg_code_last_exec_time = bg_code_exec_time;
    }
  }

  /* copy DPMEM - buffers in background */
  {
    /* call update function for CAN Tp1 CAN interface (module number: 1) */
    can_tp1_msg_copy_all_to_mem(can_tp1_address_table[0].module_addr);
  }

  /* dSPACE I/O Board DS1_RTICAN #2 */
  {
    real_T bg_code_exec_time;
    static real_T bg_code_last_exec_time = 0.0;
    bg_code_exec_time = RTLIB_TIC_READ();
    if ((bg_code_exec_time - bg_code_last_exec_time) > 0.25 ||
        (bg_code_exec_time - bg_code_last_exec_time) < 0) {
      /* ... Check taskqueue-specific error variables */
      rtican_type1_tq_err_all_chk(can_tp1_address_table[1].module_addr, 1);
      bg_code_last_exec_time = bg_code_exec_time;
    }
  }

  /* copy DPMEM - buffers in background */
  {
    /* call update function for CAN Tp1 CAN interface (module number: 2) */
    can_tp1_msg_copy_all_to_mem(can_tp1_address_table[1].module_addr);
  }
}

__INLINE void rti_mdl_sample_input(void)
{
  /* Calls for base sample time: [0.001, 0] */
  /* dSPACE I/O Board DS1401STDADCT4 #1 Unit:ADC */

  /* dSPACE I/O Board DS1401STDADCT4 #1 Unit:ADC Group:ADC */
  adc_tp4_sw_trigger(ADC_TP4_1_MODULE_ADDR, ADC_TP4_CONV_START_4|
                     ADC_TP4_CONV_START_9|ADC_TP4_CONV_START_10|
                     ADC_TP4_CONV_START_8);

  /* dSPACE I/O Board DS1_RTICAN #3 Unit:DEFAULT */
  /* call update function for CAN Tp1 CAN interface (module number: 3) */
  can_tp1_msg_copy_all_to_mem(can_tp1_address_table[2].module_addr);

  /* dSPACE I/O Board DS1_RTICAN #1 Unit:DEFAULT */
  /* call update function for CAN Tp1 CAN interface (module number: 1) */
  can_tp1_msg_copy_all_to_mem(can_tp1_address_table[0].module_addr);

  /* dSPACE I/O Board DS1_RTICAN #2 Unit:DEFAULT */
  /* call update function for CAN Tp1 CAN interface (module number: 2) */
  can_tp1_msg_copy_all_to_mem(can_tp1_address_table[1].module_addr);

  /* Calls for sample time: [0.01, 0.0] */
  if (dsIsSampleHit(VCM20_M, 2)) {
    /* dSPACE I/O Board DS1_RTICAN #3 Unit:DEFAULT */
    /* call update function for CAN Tp1 CAN interface (module number: 3) */
    can_tp1_msg_copy_all_to_mem(can_tp1_address_table[2].module_addr);
  }

  /* Calls for sample time: [0.1, 0.0] */
  if (dsIsSampleHit(VCM20_M, 3)) {
    /* dSPACE I/O Board DS1401STDDIOT4 #1 Unit:BITINCH */
    {
      /* get digital signal state on channel 5 on port 1 */
      UInt16 inputValue = dio_tp4_digin_read(DIO_TP4_1_MODULE_ADDR, 1,
        DIO_TP4_MASK_CH5);
      VCM20_B.SFunction1_b = (boolean_T) (inputValue >> (5 - 1));
    }

    {
      /* get digital signal state on channel 15 on port 1 */
      UInt16 inputValue = dio_tp4_digin_read(DIO_TP4_1_MODULE_ADDR, 1,
        DIO_TP4_MASK_CH15);
      VCM20_B.SFunction1_o = (boolean_T) (inputValue >> (15 - 1));
    }

    {
      /* get digital signal state on channel 10 on port 1 */
      UInt16 inputValue = dio_tp4_digin_read(DIO_TP4_1_MODULE_ADDR, 1,
        DIO_TP4_MASK_CH10);
      VCM20_B.SFunction1_j = (boolean_T) (inputValue >> (10 - 1));
    }

    /* dSPACE I/O Board DS1_RTICAN #3 Unit:DEFAULT */
    /* call update function for CAN Tp1 CAN interface (module number: 3) */
    can_tp1_msg_copy_all_to_mem(can_tp1_address_table[2].module_addr);

    /* dSPACE I/O Board DS1_RTICAN #1 Unit:DEFAULT */
    /* call update function for CAN Tp1 CAN interface (module number: 1) */
    can_tp1_msg_copy_all_to_mem(can_tp1_address_table[0].module_addr);
  }
}

static void rti_mdl_daq_service()
{
  /* dSPACE Host Service */
  host_service(1, rtk_current_task_absolute_time_ptr_get());
  DsDaq_Service(0, 0, 1, (DsDaqSTimestampStruct *)
                rtk_current_task_absolute_time_ptr_get());
}

#undef __INLINE

/****** [EOF] ****************************************************************/

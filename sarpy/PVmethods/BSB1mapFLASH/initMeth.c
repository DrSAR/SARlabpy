/****************************************************************
 *
 * $Source: /pv/CvsTree/pv/gen/src/prg/methods/FLASH/initMeth.c,v $
 *
 * Copyright (c) 2002 - 2005
 * Bruker BioSpin MRI GmbH
 * D-76275 Ettlingen, Germany
 *
 * All Rights Reserved
 *
 * $Id: initMeth.c,v 1.35.2.3 2009/09/10 16:36:03 mawi Exp $
 *
 ****************************************************************/

static const char resid[] = "$Id: initMeth.c,v 1.35.2.3 2009/09/10 16:36:03 mawi Exp $(C) 2002 Bruker BioSpin MRI GmbH";

#define DEBUG		0
#define DB_MODULE	1
#define DB_LINE_NR	1


#include "method.h"

/*:=MPB=:=======================================================*
 *
 * Global Function: initMeth
 *
 * Description: This procedure is implicitly called when this
 *	method is selected.
 *
 * Error History:
 *
 * Interface:							*/

void initMeth()
/*:=MPE=:=======================================================*/
{
  int dimRange[2] = { 2,3 };
  int lowMat[3]   = { 32, 32, 8 };
  int upMat[3]    = { 2048, 2048, 512 };
  
  
  DB_MSG(( "-->initMeth\n" ));
  
  /* which version of toolboxes should be active */
  
  PTB_VersionRequirement( Yes,20090101,"");
  

  /*  Initialize NA ( see code in parsRelations ) */
  
  Local_NAveragesRange();
  Local_NMovieFramesRange();
  Local_MovieOnOffRange();
  ParxRelsMakeNonEditable("TimeForMovieFrames");

  if(ParxRelsParHasValue("PVM_NRepetitions") == No)
    PVM_NRepetitions = 1;

  /*  Initialise dummy scans */
  dsRange();

  /* 
   * Which parameter classes (see parDefs.h) should be
   * hidden in routine mode 
   */
  
  PTB_SetUserTypeClasses( "Nuclei,"
			  "Sequence_Details,"
			  "RF_Pulses" );
  
  
  /* Initialisation of rf pulse parameters */
  
  /* 1: flip angle in the scan edidor */
  
  if(ParxRelsParHasValue("PVM_ExcPulseAngle") == No)
    PVM_ExcPulseAngle = 30.0;
  ParxRelsShowInEditor("PVM_ExcPulseAngle");
  
  /*
   * 2: pulses declared in parDefinitions.h 
   * in this case - ExcPulse. We initalise it to default name, 
   * 1ms, and the flip angle given in PVM_ExcPulseAngle
   */
  
  if(ParxRelsParHasValue("ExcPulse") == No)
  {
    STB_InitRFPulse(&ExcPulse,
		    CFG_RFPulseDefaultShapename(LIB_EXCITATION),
		    1.0,
		    PVM_ExcPulseAngle);
  }
  
  ExcPulseRange();
  
  /* 3: the corresponding pulse enums */
  
  STB_InitExcPulseEnum("ExcPulseEnum");
  

  /*  Bloch-Siegert parameter initialization	  */  
  
  if(ParxRelsParHasValue("BSFreqOffset") == No)
    BSFreqOffset = 4000;
  
  if(ParxRelsParHasValue("BSPulse") == No)
  {
    STB_InitRFPulse(&BSPulse,
		    CFG_RFPulseDefaultShapename(LIB_EXCITATION),
		    4.0,
		    1.0);
  }
  
  BSPulseRange();
  STB_InitExcPulseEnum("BSPulseEnum");

  if(ParxRelsParHasValue("CollectZeroPowerData") == No)
    CollectZeroPowerData = On;
  
  /* Initialisation of nucleus */  
  STB_InitNuclei(1);

  /* Initialisation of geometry parameters */
  /* A: in-plane */

  STB_InitStandardInplaneGeoPars(2,dimRange,lowMat,upMat,No);

  /* B: slice geometry */

  STB_InitSliceGeoPars(0,0,0);

  /* set gradient limits according to desired obliqueness */
  
  ControlGradientLimits(PVM_MajSliceOri);

  /* Initialisation of atoms */
  
  /* 1: method specific initialisation */
  
  if(ParxRelsParHasValue("PVM_ReadDephaseTime") == No)
    PVM_ReadDephaseTime = 1.5;
  if(ParxRelsParHasValue("PVM_ExSliceRephaseTime") == No)
    PVM_ExSliceRephaseTime = 1.5;

  PVM_AcqStartWaitTime = 0.0;

  if(ParxRelsParHasValue("PVM_RepetitionTime") == No)
    PVM_RepetitionTime = 100.0;
  if(ParxRelsParHasValue("PVM_EchoTime") == No)
    PVM_EchoTime = 6.0;
  if(ParxRelsParHasValue("PVM_DeriveGains") == No)
    PVM_DeriveGains = Yes;
  
  /* 2: call init functions to make sure all atoms have legal values */
  
  STB_InitReadAtoms();  
  STB_InitExSliceAtoms();
  STB_Init2dPhaseAtoms();
  STB_Init3dPhaseAtoms();
  
  /* Initialisation of spoilers */
  if(ParxRelsParHasValue("ReadSpoilerDuration") == No)
    ReadSpoilerDuration = 2;
  if(ParxRelsParHasValue("ReadSpoilerStrength") == No)
    ReadSpoilerStrength = 20;
  if(ParxRelsParHasValue("SliceSpoilerDuration") == No)
    SliceSpoilerDuration = 2;
  if(ParxRelsParHasValue("SliceSpoilerStrength") == No)
    SliceSpoilerStrength = 20;



  EchoTimeModeRange();

  /* Encoding */
  STB_InitEncoding();

  /* initialize digitizer parameter */

  STB_InitDigPars();
  EffSWhRange();

  /* not a csi experiment */
  PTB_SetSpectrocopyDims( 0, 0 );
  
  /* Initialisation of modules */
  STB_InitFatSupModule();
  STB_InitMagTransModule();
  STB_InitSatSlicesModule();
  STB_InitFlowSaturationModule();
  STB_InitTriggerModule();
  STB_InitTaggingModule();
  STB_InitEvolutionModule();
  STB_InitSelIrModule();
  STB_InitBlBloodModule();

  /* initialization of method specific reconstruction */
  if(ParxRelsParHasValue("RecoMethMode") == No)
    RecoMethMode=Default;
  if(ParxRelsParHasValue("WeightingMode") == No)
    WeightingMode=positive_mask;
  GaussBroadRange();
  MaskWeightRange();

  /* Visibility of Scan Editor parameters */
  ParxRelsShowInEditor("PVM_EchoTime1,PVM_NEchoImages");
  ParxRelsHideClassInEditor("ScanEditorInterface");
  ParxRelsMakeEditable("PVM_EchoTime1"); 
  ParxRelsMakeNonEditable("PVM_NEchoImages");
  ParxRelsShowInFile("PVM_EchoTime1,PVM_EchoTime,PVM_NEchoImages");

  ParxRelsHideClassInEditor("PPGparameters"); 
  ParxRelsShowInFile("PPGparameters"); 
  /* 
   * Once all parameters have initial values, the backbone is called
   * to assure they are consistent 
   */
  
  backbone();
  
  
  DB_MSG(("<--initMeth\n"));
  
}



/****************************************************************/
/*		E N D   O F   F I L E				*/
/****************************************************************/










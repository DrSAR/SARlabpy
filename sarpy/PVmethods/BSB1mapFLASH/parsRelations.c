/****************************************************************
 *
 * $Source: /pv/CvsTree/pv/gen/src/prg/methods/FLASH/parsRelations.c,v $
 *
 * Copyright (c) 2002 - 2003
 * Bruker BioSpin MRI GmbH
 * D-76275 Ettlingen, Germany
 *
 * All Rights Reserved
 *
 *
 * $Id: parsRelations.c,v 1.47.2.7 2009/08/14 14:25:05 sako Exp $
 *
 ****************************************************************/

static const char resid[] = "$Id: parsRelations.c,v 1.47.2.7 2009/08/14 14:25:05 sako Exp $ (C) 2002 Bruker BioSpin MRI GmbH";

#define DEBUG		0
#define DB_MODULE	0
#define DB_LINE_NR	0


#include "method.h"


/* ------------------------------------------------------------ 
 * backbone 
 * The main part of method code. The consitency of all parameters 
 * is checked here, relations between them are resolved and, 
 * finally, functions setting the base level parameters are 
 * called.
 */
void backbone( void )
{
  YesNo refAttIsAvailable=No;
  int dim;
  double referenceAttenuation=0,
         minFov[3] = {1e-3, 1e-3, 1e-3},
         minThickness,mindur,fixedTime;
  int    nSlices;
  YesNo RetVal;
 
  
  
  DB_MSG(("-->backbone\n"));


  PVM_UserType  = Expert_User;
  STB_DefaultUserTypeHandler(); /* this controls param visibility
				   for expert/routine */
  
  /* 
   * Nucleus and  PVM_GradCalConst
   * are handled by this funtion: 
   */
  
  
  STB_UpdateNuclei(Yes);
  

  /* coexistance of modules: */
  if(PVM_SelIrOnOff == On)
  {
    PVM_BlBloodOnOff = Off;
    PVM_MovieOnOff = Off;
  }

  if (CollectZeroPowerData == On)
  {
    PVM_NRepetitions = 4;
  }
  else
  {
    PVM_NRepetitions = 2;
  }
  /* handle RF pulse */   
  
  if(PVM_DeriveGains == Yes)
    refAttIsAvailable =
      STB_GetRefAtt(1,PVM_Nucleus1,&referenceAttenuation);         
  else
    refAttIsAvailable = No;
  
  
  
  STB_UpdateRFPulse("ExcPulse",
		    &ExcPulse,
		    refAttIsAvailable,
		    referenceAttenuation);
  
  
  STB_UpdateExcPulseEnum("ExcPulseEnum",
			 &ExcPulseEnum,
			 ExcPulse.Filename,
			 ExcPulse.Classification);
      
  PVM_ExcPulseAngle = ExcPulse.FlipAngle;
  
  /* Bloch-Siegert pulse handling */

  STB_UpdateRFPulse("BSPulse",
		    &BSPulse,
		    refAttIsAvailable,
		    referenceAttenuation);
  
  
  STB_UpdateExcPulseEnum("BSPulseEnum",
			 &BSPulseEnum,
			 BSPulse.Filename,
			 BSPulse.Classification);

  /* sequence atoms: */

  if(EchoTimeMode == Short_TE)
  {
    PVM_ReadDephaseTime = PVM_ExSliceRephaseTime;
    ParxRelsMakeNonEditable("PVM_ReadDephaseTime");
  }  
  else
    ParxRelsMakeEditable("PVM_ReadDephaseTime");
  
  PVM_2dPhaseGradientTime = PVM_ReadDephaseTime;
  PVM_3dPhaseGradientTime = PVM_2dPhaseGradientTime;
  
  /* spoilers */
  ReadSpoilerStrength =  MIN_OF(ReadSpoilerStrength, 100.0);
  ReadSpoilerStrength =  MAX_OF(ReadSpoilerStrength, -100.0);
  SliceSpoilerStrength = MIN_OF(SliceSpoilerStrength, 100.0);
  SliceSpoilerStrength = MAX_OF(SliceSpoilerStrength, -100.0);
  SliceSpoilerDuration = MAX_OF(SliceSpoilerDuration, 2*PVM_RiseTime);
  
  
  
  /* excitation pulse */
  PVM_ExSlicePulseLength   = ExcPulse.Length;
  PVM_ExSliceBandWidth     = ExcPulse.Bandwidth;
  PVM_ExSliceRephaseFactor = ExcPulse.RephaseFactor * 
                             ExcPulse.TrimRephase / 100.0;
  
  
  /* begin Update Geometry: */
  
  /* 
   * 1: in-plane geometry:
   * The STB_StandardInplaneGeoParHandler is called twice:
   * first, with a dummy value of minFov, to make size constraints;
   * then, after the true minFov is found, to do the rest.
   * (because the sizes must be set before we find minFov)
   */
      
  STB_StandardInplaneGeoParHandler(minFov,2.0);

  /* once the image size is decided (after first update of inpl geo), we
     can update the Encoding parameters, to get the acquisition size (PVM_EncMatrix).
  */
  PVM_RareFactor = 1;

  STB_UpdateEncoding(PTB_GetSpatDim(),  /* total dimensions */
		     PVM_Matrix,        /* image size */ 
		     PVM_AntiAlias,     /* a-alias */
		     &PVM_RareFactor,   /* segment size */
                     SEG_SEQUENTIAL,    /* segmenting mode */
		     Yes,               /* ppi in 2nd dim allowed */
		     Yes,               /* ppi ref lines in 2nd dim allowed */
		     Yes);              /* partial ft in 2nd dim allowed */ 

  /* update bandwidth and acquisition time */

  STB_UpdateDigPars(&PVM_EffSWh,
		    PVM_EncMatrix[0],
		    PVM_AntiAlias,
		    &PVM_AcquisitionTime);

  /*
   *  assure a minimum DEOSC that allows full scan shift
   */

  mindur =  MAX_OF(PVM_2dPhaseGradientTime+PVM_RiseTime,PVM_DigEndDelOpt);

  ReadSpoilerDuration = MAX_OF(mindur,ReadSpoilerDuration);


  ControlGradientLimits(PVM_MajSliceOri);

  LocalGeometryMinimaRels(minFov, &minThickness);

  dim = PTB_GetSpatDim();

  /*
   * Constrain minimum fov in 3D phase direction and 
   * minimum slice thickness in case of 3D measurements
   */

  if(dim == 3)
  {
    double min;

    min=MAX_OF(minFov[2],minThickness);
    minFov[2] = minThickness = min;
    
  }

  STB_StandardInplaneGeoParHandler(minFov,2.0);
  
    
  /* 2: slice geometry 
   *   The maximun slice number per packages is set to 1 for 3D
   *   and free for 2D . The total maximum number of slices is 
   *   defined by the return value of CFG_MaxSlices() called by
   *   STB_UpdateSliceGeoPars if no constrain is given by the
   *   arguments (value 0)
   */
 


  if(dim == 3)
  {

    /*
     * Connect slice thickness to FOV in 3rd direction.
     */
    PVM_SliceThick = PVM_Fov[2];

    /* constrain maximum slices per package to 1 */

    STB_UpdateSliceGeoPars(0,0,1,minThickness);
  }
  else
  {
    /* no constrain to slices in 2D mode */

    STB_UpdateSliceGeoPars(0,0,0,minThickness);
    
  }
    /* constrain to 1 slice package for BlackBlood */

  if(PVM_BlBloodOnOff == On)
    STB_UpdateSliceGeoPars(0,1,0,minThickness);
 
    /* end Update Geometry */



  LocalGradientStrengthRels();  
  LocalFrequencyOffsetRels();
  
  PVM_NRepetitions = MAX_OF(1,PVM_NRepetitions);
  
  PVM_NEchoImages = 1;
  
  Local_NAveragesRange();

  if (ParxRelsParHasValue("PVM_MotionSupOnOff") == 0)
     PVM_MotionSupOnOff = Off;


  PARX_hide_pars(NOT_HIDDEN,"PVM_MotionSupOnOff");
 
  /* handling of modules */

  STB_UpdateFatSupModule(PVM_Nucleus1);
  STB_UpdateMagTransModule();
  STB_UpdateSatSlicesModule(PVM_Nucleus1);
  STB_UpdateFlowSaturationModule(PVM_Nucleus1);
  STB_UpdateTriggerModule();
  STB_UpdateTaggingModule(PVM_Nucleus1,PVM_Fov,PVM_Matrix);  
  STB_UpdateEvolutionModule(&EvolutionDuration);
  
  fixedTime = SliceSpoilerDuration +CFG_AmplifierEnable() + PVM_ExSlicePulseLength/2.0 + PVM_TaggingModuleTime;

  /* Update for Black Blood Module   */
  {
  double slabthick[1];
  slabthick[0] = (PVM_SPackArrSliceDistance[0] * (PVM_SPackArrNSlices[0]-1)) + PVM_SliceThick;
  STB_UpdateBlBloodModule(PVM_Nucleus1,slabthick,PVM_SPackArrSliceOffset,1,fixedTime);
  }
  
 
  Local_UpdateMovie();
  echoTimeRels();

  /* update slice-selective inversion */
  SliceSegDurRels();
  nSlices = GTB_NumberOfSlices( PVM_NSPacks, PVM_SPackArrNSlices);
  RetVal = STB_UpdateSelIrModule(PVM_Nucleus1,PVM_SliceThick,PVM_SliceOffset,nSlices,&SliceSegDur,0,fixedTime);
  if(PVM_SelIrOnOff==On)
  ParxRelsCopyPar("PVM_InversionTime","PVM_SelIrInvTime");
  if(PVM_BlBloodOnOff==On)
  ParxRelsCopyPar("PVM_InversionTime","PVM_BlBloodInvTime");
  if(PVM_BlBloodOnOff==Off&&PVM_SelIrOnOff==Off)
  {
  PVM_InversionTime = 0.0;
  ParxRelsHideInEditor("PVM_InversionTime");
  }
  repetitionTimeRels();
  
  /* visibility of method specific reconstruction parameters */
  RecoMethModeVisPar();

  /* set baselevel acquisition parameter */
  SetBaseLevelParam();

  /*setting baselevel reconstruction parameters */
   SetRecoParam(); 
   
  DB_MSG(("<--backbone\n"));
}


/*==============================================================
 * relation of NDummyScans
 *==============================================================*/


void dsRelations(void)
{

  DB_MSG(("-->dsRelations\n"));

  dsRange(); 
  backbone();

  DB_MSG(("<--dsRelations\n"));

}

void dsRange(void)
{
  if(ParxRelsParHasValue("NDummyScans") == No)
    NDummyScans = 0;
  else
    NDummyScans = MAX_OF(0, NDummyScans);
} 



/*=======================================================
 *ExcPulseAngleRelation
 *Redirected relation of PVM_ExcPulseAngle
 *======================================================= */

void ExcPulseAngleRelation(void)
{
  DB_MSG(("-->ExcPulseAngleRelation"));

  ExcPulse.FlipAngle = PVM_ExcPulseAngle;
  ExcPulseRange();
  backbone();

  DB_MSG(("<--ExcPulseAngleRelation"));
}

/*==========================================================
 *
 *  examples for relations concearning special pulses and 
 *  pulselists
 *
 *==========================================================*/



/*===============================================================
 * ExcPulseEnumRelation
 * Relation of ExcPulseEnum (a dynamic enmueration parameter which
 * allows to select one of the existing library exc. pulses)
 * Sets the name and the clasification  of the pulse perameter ExcPulse 
 * according to the selected enum value.
 *===============================================================*/

void ExcPulseEnumRelation(void)
{
  YesNo status;
  DB_MSG(("-->ExcPulsesEnumRelation\n"));
  
  /* set the name and clasification of ExcPulse: */

  status = STB_UpdateExcPulseName("ExcPulseEnum",
				  &ExcPulseEnum,
				  ExcPulse.Filename,
				  &ExcPulse.Classification);
  
  /* call the method relations */
  backbone();
  
  DB_MSG(("<--ExcPulseEnumRelation\n"));                                       
}



/* ===================================================================
 * Relation of ExcPulse
 * 
 * All pulses of type PVM_RF_PULSE_TYPE must have relations like this.
 * However, if you clone this funtion for a different pulse parameter
 * remember to replace the param name in the call to UT_SetRequest!
 *
 * IMPORTANT: this function should not be invoked in the backbone!
 ====================================================================*/

void ExcPulseRelation(void)
{
  DB_MSG(("-->ExcPulseRelation\n"));
  
  /*
   * Tell the request handling system that the parameter
   * ExcPulse has been edited 
   */
  
  UT_SetRequest("ExcPulse");
  
  /* Check the values of ExcPulse */
  
  ExcPulseRange();
  
  /* 
   * call the backbone; further handling will take place there
   * (by means of STB_UpdateRFPulse)  
   */
  
  backbone();
  
  DB_MSG(("<--ExcPulseRelation\n"));
}


/*================================================================
 *	         L O C A L   F U N C T I O N S			
 *==============================================================*/



void ExcPulseRange(void)
{
  DB_MSG(("-->ExcPulseRange\n"));
  
  /* allowed clasification */
  
  switch(ExcPulse.Classification)
  {
    default:
      
      ExcPulse.Classification = LIB_EXCITATION;
      break;
    case LIB_EXCITATION:
    case PVM_EXCITATION:
    case USER_PULSE:
      break;
  }

  /* allowed angle for this pulse */
  
  ExcPulse.FlipAngle = MIN_OF(90.0,ExcPulse.FlipAngle);
  
  
  /* general verifiation of all pulse atributes  */
  
  STB_CheckRFPulse(&ExcPulse);
  
  DB_MSG(("<--ExcPulseRange\n"));
  
}

void echoTimeRels( void )
{
  DB_MSG(("-->echoTimeRels\n"));
  
  PVM_MinEchoTime = 
    PVM_ExSlicePulseLength / 2.0 +
    PVM_RampTime                 +
    PVM_InterGradientWaitTime    +
    PVM_ExSliceRephaseTime       +
    CFG_AmplifierEnable()        +
    BSPulse.Length		 +
    PVM_InterGradientWaitTime    +
    PVM_RiseTime                 +
    PVM_AcquisitionTime * PVM_EchoPosition / 100.0;

  if(EchoTimeMode == Long_TE) 
    PVM_MinEchoTime += (PVM_ReadDephaseTime+PVM_InterGradientWaitTime);
  
  PVM_EchoTime = PVM_EchoTime < PVM_MinEchoTime   ?
                 PVM_MinEchoTime : PVM_EchoTime;
  
  
  /* Set Echo Parameters for Scan Editor  */  
  PVM_EchoTime1 = PVM_EchoTime;

  /* set ppg flag related to echo time mode */
  PVM_ppgFlag1 = (EchoTimeMode == Long_TE)? Yes : No;

  DB_MSG(("<--echoTimeRels\n"));
}


void SliceSegDurRels( void)
{
  double amplifierenable;

  amplifierenable = CFG_AmplifierEnable();

  MinSliceSegDur = minLoopRepetitionTime();

  /* period of the slice loop in the IR mode: */ 
  SliceSegDur = MinSliceSegDur;

}

double minLoopRepetitionTime(void)
/* ---------------------------------------------------------
this function rturns the minimum duration of the innermost
pulse program loop
----------------------------------------------------------*/
{
  double minTr, minD0;
  
  minD0 = 0.04  /* ADC_END_4ch */ + PVM_InterGradientWaitTime;

  minTr = 
	PVM_FatSupModuleTime         +
	PVM_MagTransModuleTime       +
	PVM_FovSatModuleTime         +
	PVM_InFlowSatModuleTime      +
	SliceSpoilerDuration         +
	CFG_AmplifierEnable()        +
	PVM_ExSlicePulseLength / 2.0 +
	PVM_EchoTime                 +
	PVM_AcquisitionTime *(100.0 - PVM_EchoPosition) / 100.0 +
	ReadSpoilerDuration          +
	PVM_RiseTime                 +
        0.02                         +
	minD0;

  return minTr;
}

void repetitionTimeRels( void )
{
  int nSlices,dim;
  double TotalTime, minD0;

  DB_MSG(("-->repetitionTimeRels\n"));

  TotalTime = 0.0;
  nSlices = GTB_NumberOfSlices( PVM_NSPacks, PVM_SPackArrNSlices );

  minD0 = 0.04  /* ADC_END_4ch */ + PVM_InterGradientWaitTime;
  
  if( PVM_ErrorDetected == Yes )
  {
    UT_ReportError("FLASH-backbone: In function call!");
    return;
  }
  
  if(PVM_SelIrOnOff == On)
  {
    PVM_MinRepetitionTime =  
      PVM_SelIrModuleTime +
      PVM_TaggingModuleTime +
      nSlices * SliceSegDur;
  }
  else
  {
    /* min TR in a movie: */
    PVM_MinRepetitionTime = nSlices * minLoopRepetitionTime();
 
    /* if there is no movie, TR also includes some modules: */
    if(PVM_MovieOnOff == Off || PVM_NMovieFrames == 1)
    {
      PVM_MinRepetitionTime += 
	PVM_BlBloodModuleTime +
	PVM_SelIrModuleTime +
	PVM_TaggingModuleTime;
    } 
  }  
   
  PVM_RepetitionTime = ( PVM_RepetitionTime < PVM_MinRepetitionTime ? 
			  PVM_MinRepetitionTime : PVM_RepetitionTime );
 

  /* delay after the slice loop, used only in IR mode to cotrol TR: */
  SliceSegEndDelay = PVM_SelIrOnOff==On? PVM_RepetitionTime - PVM_MinRepetitionTime : 0;
 
  
  /* Calculate Total Scan Time and Set for Scan Editor */ 
  dim = PTB_GetSpatDim();

  if( dim >1 )
  {
    if(PVM_MovieOnOff == On && PVM_NMovieFrames >1 )
    {    
      /* TR is one movie frame, without prep modules: */
      TimeForMovieFrames = PVM_RepetitionTime * PVM_NMovieFrames + PVM_TaggingModuleTime;
      TotalTime = (PVM_BlBloodModuleTime + TimeForMovieFrames)  
	          * PVM_EncMatrix[1] * PVM_NAverages; 
    }
    else
    {
      /* TR includes prep modules */
      TimeForMovieFrames = PVM_RepetitionTime - PVM_BlBloodModuleTime - PVM_TaggingModuleTime;
      TotalTime = PVM_RepetitionTime * PVM_EncMatrix[1] * PVM_NAverages;   
    } 
  }                     

  if( dim >2 )
    TotalTime = TotalTime * PVM_EncMatrix[2];

  /* time for one repetition */
  OneRepTime = TotalTime/1000.0;

  if(PVM_EvolutionOnOff && (PVM_NEvolutionCycles > 1))
      TotalTime = TotalTime *PVM_NEvolutionCycles + EvolutionDuration;
  
  TotalTime = TotalTime * PVM_NRepetitions;
  
  UT_ScanTimeStr(PVM_ScanTimeStr,TotalTime);
  
  
  ParxRelsShowInEditor("PVM_ScanTimeStr");
  ParxRelsMakeNonEditable("PVM_ScanTimeStr");
  
  DB_MSG(("<-- repetitionTimeRels\n"));
}

void LocalGeometryMinimaRels(double *min_fov, double *min_thickness )
{
  /*
   * This function calculates the minima for the fields-of-view in all
   * three directions as well as the minimum slice thickness. It is always
   * assumed that all slices have the same thickness
   * 
   * The results is returned in min_fov[0..2] and min_thickness.
   */

  double readRampInteg;  /* normalised integral rising read gradient ramp   */
  double readDephInteg;  /* normalised integral read dephase gradient       */

  double sliceRampInteg; /* normalised integral falling slice gradient ramp */
  double sliceRephInteg; /* normalised integral slice rephase gradient      */


  DB_MSG(("-->LocalGeometryMinimaRels\n"));
    
  switch( PTB_GetSpatDim() )
    {
    case 3: 
      /* PHASE ENCODING GRADIENT - 3nd DIM 
       *
       * Step #1:
       * Calculate the normalised integral of the phase encoding gradient
       * in the 3rd dimension
       *
       * The variable Phase3dInteg is a parameter defined in the file:
       * parsDefinition.h but it is NOT included in the definition of 
       * MethodClass
       * that appears in parsLayout.h, The value of Phase3dInteg determined 
       * here is used later in "LocalGradientStrengthRels()" 
       */

      Phase3dInteg = MRT_NormGradPulseTime ( PVM_3dPhaseGradientTime,
					     PVM_3dPhaseRampUpTime,
					     PVM_3dPhaseRampUpIntegral,
					     PVM_3dPhaseRampDownTime,
					     PVM_3dPhaseRampDownIntegral );
      /*
       * Step #2:
       * Calculate the resulting minimum field-of-view in that direction
       */

      min_fov[2] = MRT_PhaseFov( Phase3dInteg,
				 PVM_Matrix[2],
				 PVM_Lim3dPhaseGradient,
				 PVM_GradCalConst );
      /* falls through */
    case 2: 
      /* PHASE ENCODING GRADIENT - 2nd DIM 
       *
       * the same for the second dimension
       *
       * The variable Phase2dInteg is a parameter defined in the file:
       * parsDefinition.h but it is NOT included in the definition of 
       * MethodClass that appears in parsLayout.h, The value of 
       * Phase2dInteg determined here is used later in 
       * "LocalGradientStrengthRels()" 
       */
      
      Phase2dInteg = MRT_NormGradPulseTime ( PVM_2dPhaseGradientTime,
					     PVM_2dPhaseRampUpTime,
					     PVM_2dPhaseRampUpIntegral,
					     PVM_2dPhaseRampDownTime,
					     PVM_2dPhaseRampDownIntegral );
      
      min_fov[1] = MRT_PhaseFov( Phase2dInteg,
				 PVM_Matrix[1],
				 PVM_Lim2dPhaseGradient,
				 PVM_GradCalConst );
      /* falls through */
    case 1:
 
      /* FREQUENCY ENCODING GRADIENT
       *
       * Step #1:
       * Calculate the normalised integral of the read dephase gradient
       */

      readDephInteg = MRT_NormGradPulseTime( PVM_ReadDephaseTime,
					     PVM_ReadDephaseRampUpTime,
					     PVM_ReadDephaseRampUpIntegral,
					     PVM_ReadDephaseRampDownTime,
					     PVM_ReadDephaseRampDownIntegral );
      /*
       * Step #2:
       * Calculate the normalised integral of the rising ramp of the read
       * gradient preceding data collection
       */
      readRampInteg = MRT_NormGradRampTime( PVM_ReadRampUpTime,
					    PVM_ReadRampUpIntegral );
      
      /*
       * Step #3
       * Calculate the ratio of the strength of the read gradient to the
       * strength of the read dephase gradient
       *
       * The variable ReadGradRatio is a parameter defined in the file:
       * parsDefinition.h but it is NOT included in the definition of 
       * MethodClass that appears in parsLayout.h, The value of ReadGradRatio 
       * determined here is used later in "LocalGradientStrengthRels()" 
       */
      
      ReadGradRatio = MRT_ReadGradRatio( PVM_AcquisitionTime,
					 PVM_EchoPosition,
					 PVM_AcqStartWaitTime,
					 readDephInteg,
					 readRampInteg );
      /*
       * Step #4
       * Calculate the minimum field of view in the read direction
       */
      
      min_fov[0] = MRT_MinReadFov( PVM_EffSWh,
				   ReadGradRatio,
				   PVM_LimReadGradient,
				   PVM_LimReadDephaseGradient,  
				   PVM_GradCalConst );
      /* falls through */
    default: 
      /* SLICE SELECTION GRADIENT
       *
       * Calculate the normalised integral of the descending gradient 
       * ramp after the RF pulse
       */

      sliceRampInteg = MRT_NormGradRampTime( PVM_ExSliceRampDownTime,
					     PVM_ExSliceRampDownIntegral );
      /*
       * Calculate the normalised integral of the slice selection rephasing
       * gradient
       */

      sliceRephInteg = 
	MRT_NormGradPulseTime( PVM_ExSliceRephaseTime,
			       PVM_ExSliceRephaseRampUpTime,
			       PVM_ExSliceRephaseRampUpIntegral,
			       PVM_ExSliceRephaseRampDownTime,
			       PVM_ExSliceRephaseRampDownIntegral );
      /*
       * Calculate the ratio of the strength of the slice selection 
       * gradient to the strength of the slice selection rephase 
       * gradient
       *
       * The variable SliceGradRatio is a parameter defined in the file:
       * parsDefinition.h but it is NOT included in the definition of 
       * MethodClass that appears in parsLayout.h. 
       * The value of SliceGradRatio determined here is used later in 
       * "LocalGradientStrengthRels()"
       */

      SliceGradRatio = 
	MRT_SliceGradRatio( PVM_ExSlicePulseLength,
			    PVM_ExSliceRephaseFactor,
			    PVM_ExSliceRampDownWaitTime,
			    sliceRampInteg,
			    sliceRephInteg );
      /*
       * Calculate the minimum slice thickness
       */
      
      *min_thickness = MRT_MinSliceThickness( PVM_ExSliceBandWidth,
					      SliceGradRatio,
					      PVM_LimExSliceGradient,
					      PVM_LimExSliceRephaseGradient,
					      PVM_GradCalConst );
      break;
    }
  DB_MSG(("<--LocalGeometryMinimaRels\n"));
}



void LocalGradientStrengthRels( void )
{
  int dim;

  DB_MSG(("-->LocalGradientStrengthRels\n"));
  
  /*
   * This function calculates all the gradient strengths 
   */

  dim = PTB_GetSpatDim();

  /*
   * PHASE ENCODING GRADIENT - 3nd DIM 
   */

  if(dim > 2)
  {
      PVM_3dPhaseGradient = 
	MRT_PhaseGrad( Phase3dInteg,
		       PVM_Matrix[2],
		       PVM_Fov[2],
		       PVM_GradCalConst );
  }  
  else
    PVM_3dPhaseGradient = 0.0;

  /* 
   * PHASE ENCODING GRADIENT - 2nd DIM 
   */
 
  if(dim > 1)
  {     
      PVM_2dPhaseGradient = 
	MRT_PhaseGrad( Phase2dInteg,
		       PVM_Matrix[1],
		       PVM_Fov[1],
		       PVM_GradCalConst );
  }
  else
    PVM_2dPhaseGradient = 0.0;
      
  /* FREQUENCY ENCODING GRADIENT */

  if(dim >0)
  {
  PVM_ReadGradient = 
    MRT_ReadGrad( PVM_EffSWh,
		  PVM_Fov[0],
		  PVM_GradCalConst );
  PVM_ReadDephaseGradient = 
    MRT_ReadDephaseGrad( ReadGradRatio,
			 PVM_ReadGradient );
  }
  else
    PVM_ReadGradient = PVM_ReadDephaseGradient = 0.0;

  /* SLICE SELECTION GRADIENT */

  PVM_ExSliceGradient = 
    MRT_SliceGrad( PVM_ExSliceBandWidth,
		   PVM_SliceThick,
		   PVM_GradCalConst );
  
  PVM_ExSliceRephaseGradient = 
    MRT_SliceRephaseGrad( SliceGradRatio,
			      PVM_ExSliceGradient );
  
  
  DB_MSG(("<--LocalGradientStrengthRels\n"));
  
}

void LocalFrequencyOffsetRels( void )
{
  int spatDim;
  int i,nslices;

  DB_MSG(("-->LocalFrequencyOffsetRels\n"));

  spatDim = PTB_GetSpatDim();

  nslices = GTB_NumberOfSlices(PVM_NSPacks,PVM_SPackArrNSlices);

  /*
   * Calculate offset in read direction
   */
  
  MRT_FrequencyOffsetList(nslices,
			  PVM_EffReadOffset,
			  PVM_ReadGradient,
			  PVM_GradCalConst,
			  PVM_ReadOffsetHz );

  /*
   * Calculate slice offset
   */

  MRT_FrequencyOffsetList(nslices,
			  PVM_EffSliceOffset,
			  PVM_ExSliceGradient,
			  PVM_GradCalConst,
			  PVM_SliceOffsetHz );

  if(spatDim == 3)
  {
    for(i=0;i<nslices;i++)
      PVM_EffPhase2Offset[i] = -PVM_EffSliceOffset[i];
  }
  
  DB_MSG(("<--LocalFrequencyOffsetRels\n"));
  
}



/*===============================================================
 *
 * Range checking routine for parameter PVM_NAverages
 *
 *==============================================================*/


void Local_NAveragesRange(void)
{
  int ival;
  DB_MSG(("-->Local_NAveragesRange\n"));
  
  /* 
   *  Range check of PVM_NAverages: prevent it to be negative or 0
   */
  
  if(ParxRelsParHasValue("PVM_NAverages") == No)
  {
    PVM_NAverages = 1;
  }

  ival = PVM_NAverages;
  PVM_NAverages = MAX_OF(ival,1);
  
  
  DB_MSG(("<--Local_NAveragesRange\n"));

}



void Local_NAveragesHandler(void)
{

  DB_MSG(("-->Local_NAveragesRange\n"));


  Local_NAveragesRange();

  /*
   *   Averages range check is finished, handle the request by
   *   the method:
   */
  
  
  backbone();
  
  
  DB_MSG(("<--Local_NAveragesRange\n"));
  return;
}

/*===========================================================
 * The update of slice geometry is done in backbone.
 * This function connects the Fov in 3rd spatial dimension
 * to the slice thickness.
 *==========================================================*/


void localHandleSliceGeometry(void)
{
  DB_MSG(("-->localHandleSliceGeometry\n"));

  if(PTB_GetSpatDim()==3)
  {

    STB_IsotropicRange();

    switch(PVM_Isotropic)
    {
      case Isotropic_Fov:
      case Isotropic_All:
	PVM_Fov[0] = 
	PVM_Fov[1] = 
	PVM_Fov[2] = PVM_SliceThick;
	break;
      default:
	PVM_Fov[2] = PVM_SliceThick;
	break;
    }
  }

  backbone();


  DB_MSG(("<--localHandleSliceGeometry\n"));
  return;
}

/*==================================================
 *     relation to update Echotime from scan editor
 *==================================================*/

void LocalEchoTime1Relation(void)
{
  DB_MSG(("-->LocalEchoTime1Relation\n"));

  PVM_EchoTime = PVM_EchoTime1;
  backbone();

  DB_MSG(("<--LocalEchoTime1Relation\n"));
}

void Local_NMovieFramesRange(void)
{
  
  /*
   *  Range check of PVM_NMovieFrames: prevent it to be negative or 0
   */
  if(ParxRelsParHasValue("PVM_NMovieFrames") == No)
  {
    PVM_NMovieFrames = 1;
  }
  PVM_NMovieFrames = MAX_OF(PVM_NMovieFrames,1);
  PVM_NMovieFrames = MIN_OF(PVM_NMovieFrames,200);
  
  
}

void Local_MovieOnOffRange(void)
{
  
  /*
   *  Range check of PVM_MovieOnOff
   */
  if(ParxRelsParHasValue("PVM_MovieOnOff") == No)
  {
    PVM_MovieOnOff = Off;
  }
       
}


void Local_NMovieFramesRels(void)
{
  Local_NMovieFramesRange();
  backbone();
  
}

void Local_MovieOnOffRels(void)
{
  Local_MovieOnOffRange();
  backbone();
}


void Local_UpdateMovie(void)
{
  if(PVM_MovieOnOff == Off)
  {
    ParxRelsHideInEditor("PVM_NMovieFrames,TimeForMovieFrames");
  }
  else
  {
    ParxRelsShowInEditor("PVM_NMovieFrames,TimeForMovieFrames");
  }
}

/* relations and range checking of EchoTimeMode */
void EchoTimeModeRels(void)
{
  EchoTimeModeRange();
  backbone();
}

void  EchoTimeModeRange(void)
{
  if(ParxRelsParHasValue("EchoTimeMode") == No)
    EchoTimeMode = Long_TE;

  if(EchoTimeMode != Long_TE && EchoTimeMode != Short_TE)
    EchoTimeMode = Long_TE;
}

/* rangechecking and redirected relations of PVM_EffSWh */

void EffSWhRange(void)
{
  DB_MSG(("-->EffSWhRange"));

  if(!ParxRelsParHasValue("PVM_EffSWh"))
  {
    PVM_EffSWh = 50000;
  }
  else
  {
    PVM_EffSWh = MIN_OF(MAX_OF(PVM_EffSWh,2000.0),1000000);
  }

  DB_MSG(("-->EffSWhRange"));
  return;
}

void EffSWhRel(void)
{
  DB_MSG(("-->EffSWhRel"));

  EffSWhRange();
  backbone();

  DB_MSG(("-->EffSWhRel"));
  return;
}

void localInversionRel(void)
{
    DB_MSG(("-->localInversionRel"));
    if(!ParxRelsParHasValue("PVM_InversionTime"))
       PVM_InversionTime = 0.0;
    PVM_InversionTime = MAX_OF(PVM_InversionTime,0.0);
    if(PVM_SelIrOnOff==On)
    PVM_SelIrInvTime = PVM_InversionTime;
    if(PVM_BlBloodOnOff==On)
    PVM_BlBloodInvTime = PVM_InversionTime;
    backbone();
    DB_MSG(("-->localInversionRel"));

}

void ControlGradientLimits(YesNo NotOblique)
{
  DB_MSG(("-->ControlGradientLimits: Obliqueness forbidden: %s",NotOblique==Yes ? "Yes":"No"));

  if(NotOblique==Yes)
  {
    PVM_LimReadDephaseGradient = 100.0;
    PVM_Lim2dPhaseGradient = 100.0;
    PVM_Lim3dPhaseGradient = 100.0;
  }
  else
  {
    /* Gradient limits in % of max. Value 57 (1/sqrt(3))
       is needed for arbitrary oblique slices. */
    PVM_LimReadDephaseGradient = 57.0;
    PVM_Lim2dPhaseGradient = 57.0;
    PVM_Lim3dPhaseGradient = 57.0;
  }


  DB_MSG(("-->ControlGradientLimits"));
}

/* relations for method specific reconstruction parameters  */
void RecoMethModeRel(void)
{
  DB_MSG(("-->RecoMethModeRel"));

  RecoMethModeVisPar();

  SetNewRecoParam();

  DB_MSG(("<--RecoMethModeRel"));
}

void RecoMethModeVisPar(void)
{
  DB_MSG(("-->RecoMethModeVisPar"));

  if(RecoMethMode==SWI)
  {
   ParxRelsShowInEditor("WeightingMode,GaussBroadening");
   if(WeightingMode==phase_image)
     ParxRelsHideInEditor("MaskWeighting");
   else
     ParxRelsShowInEditor("MaskWeighting"); 
  }
  else
    ParxRelsHideInEditor("WeightingMode,GaussBroadening,MaskWeighting");


  DB_MSG((">--RecoMethModeVisPar"));
}

void MaskModeRel(void)
{
  DB_MSG(("-->MaskModeRel"));
  
  if(WeightingMode==phase_image)
    ParxRelsHideInEditor("MaskWeighting");
  else
    ParxRelsShowInEditor("MaskWeighting"); 

  DB_MSG(("<--MaskModeRel"));
}

void GaussBroadRange(void)
{
  double max;

  DB_MSG(("-->GaussBroadRange"));

  if(ParxRelsParHasValue("GaussBroadening") == No)
    GaussBroadening=1.0;

  max=MAX_OF(PVM_Fov[0],PVM_Fov[1]);
  if(PTB_GetSpatDim()==3)
    max=MAX_OF(max,PVM_Fov[2]);
  
  GaussBroadening = MIN_OF(MAX_OF(0,GaussBroadening),max);

  DB_MSG(("<--GaussBroadRange"));
}

void MaskWeightRange(void)
{
  DB_MSG(("-->MaskWeightRange"));

  if(ParxRelsParHasValue("MaskWeighting") == No)
    MaskWeighting=4.0;

  MaskWeighting = MIN_OF(MAX_OF(0,MaskWeighting),20);

  DB_MSG(("<--MaskWeightRange"));
}


/*===============================================================
 * BSPulseEnumRelation
 * Relation of BSPulseEnum (a dynamic enmueration parameter which
 * allows to select one of the existing library exc. pulses)
 * Sets the name and the clasification  of the pulse perameter BSPulse 
 * according to the selected enum value.
 *===============================================================*/

void BSPulseEnumRelation(void)
{
  YesNo status;
  DB_MSG(("-->BSPulsesEnumRelation\n"));
  
  /* set the name and clasification of BSPulse: */

  status = STB_UpdateExcPulseName("BSPulseEnum",
				  &BSPulseEnum,
				  BSPulse.Filename,
				  &BSPulse.Classification);
  
  /* call the method relations */
  backbone();
  
  DB_MSG(("<--BSPulseEnumRelation\n"));                                       
}



/* ===================================================================
 * Relation of BSPulse
 * 
 * All pulses of type PVM_RF_PULSE_TYPE must have relations like this.
 * However, if you clone this funtion for a different pulse parameter
 * remember to replace the param name in the call to UT_SetRequest!
 *
 * IMPORTANT: this function should not be invoked in the backbone!
 ====================================================================*/

void BSPulseRelation(void)
{
  DB_MSG(("-->BSPulseRelation\n"));
  
  /*
   * Tell the request handling system that the parameter
   * BSPulse has been edited 
   */
  
  UT_SetRequest("BSPulse");
  
  /* Check the values of BSPulse */
  
  BSPulseRange();
  
  /* 
   * call the backbone; further handling will take place there
   * (by means of STB_UpdateRFPulse)  
   */
  
  backbone();
  
  DB_MSG(("<--BSPulseRelation\n"));
}


void BSPulseRange(void)
{
  DB_MSG(("-->BSPulseRange\n"));
  
  /* allowed clasification */
  
  switch(BSPulse.Classification)
  {
    default:
      
      BSPulse.Classification = LIB_EXCITATION;
      break;
    case LIB_EXCITATION:
    case PVM_EXCITATION:
    case USER_PULSE:
      break;
  }

  
  /* general verifiation of all pulse atributes  */
  
  STB_CheckRFPulse(&BSPulse);
  
  DB_MSG(("<--BSPulseRange\n"));
  
}
/****************************************************************/
/*		E N D   O F   F I L E				*/
/****************************************************************/









/****************************************************************
 *
 * $Source: /pv/CvsTree/pv/gen/src/prg/methods/FLASH/parsLayout.h,v $
 *
 * Copyright (c) 1999-2007
 * Bruker BioSpin MRI GmbH
 * D-76275 Ettlingen, Germany
 *
 * All Rights Reserved
 *
 * $Id: parsLayout.h,v 1.21.2.2 2009/08/14 14:25:05 sako Exp $
 *
 ****************************************************************/

/****************************************************************/
/*	PARAMETER CLASSES				       	*/
/****************************************************************/


/*--------------------------------------------------------------*
 * Definition of the PV class...
 *--------------------------------------------------------------*/

parclass
{
  PVM_EffSWh;
  PVM_EchoPosition;
  EchoTimeMode;
  PVM_ReadDephaseTime;
  PVM_ExSliceRephaseTime;
  PVM_MinFov;
  PVM_MinSliceThick;
  ReadSpoilerDuration;
  ReadSpoilerStrength;
  SliceSpoilerDuration;
  SliceSpoilerStrength;
  DigitizerPars;
}
attributes
{
  display_name "Sequence Details";
} Sequence_Details;

parclass
{
  ExcPulseEnum;
  ExcPulse;
} 
attributes
{
  display_name "RF Pulses";
} RF_Pulses;

parclass
{
  NDummyScans;

  PVM_MovieOnOff;
  PVM_NMovieFrames;
  TimeForMovieFrames;

  PVM_EvolutionOnOff;
  Evolution_Parameters;

  PVM_TriggerModule;
  Trigger_Parameters;

  PVM_TaggingOnOff;
  Tagging_Parameters;

  PVM_SelIrOnOff;
  Selective_IR_Parameters;

  PVM_BlBloodOnOff;
  BlackBlood_Parameters;

  PVM_FatSupOnOff;
  Fat_Sup_Parameters;

  PVM_MagTransOnOff;
  Magn_Transfer_Parameters;

  PVM_FovSatOnOff;
  Sat_Slices_Parameters;

  PVM_InFlowSatOnOff;
  Flow_Sat_Parameters;

  PVM_MotionSupOnOff;

} Preparation;

parclass
{
  RecoMethMode;
  WeightingMode;
  MaskWeighting;
  GaussBroadening;
  RECO_wordtype;
  RECO_map_mode;
  RECO_map_percentile;
  RECO_map_error;
  RECO_map_range;
  RECO_map_user_slope;
  RECO_map_user_offset;
}attributes
{
  display_name "Reconstruction Options";
}ReconstructionOptions;

/* The following class is defined to assure that parameters of the 
 * scan editor are properly displayed. */
parclass
{
  PVM_EchoTime1;
  PVM_NEchoImages;
} ScanEditorInterface;

/* Thhe following class is defined to assure that parameters used in 
 * ppg have the correct value for display ppg.*/

parclass
{
  PVM_ppgFlag1;
}PPGparameters;

parclass
{
  Method;
  PVM_EchoTime;
  PVM_RepetitionTime;
  PVM_NAverages;
  PVM_NRepetitions;
  PVM_ScanTimeStr;
  BSFreqOffset;
  BSPulseEnum;
  BSPulse;
  CollectZeroPowerData;
  PVM_ExcPulseAngle;
  PVM_DeriveGains;
  RF_Pulses;
  Nuclei;
  Encoding;
  Sequence_Details;
  StandardInplaneGeometry;
  StandardSliceGeometry;
  Preparation;
  ScanEditorInterface;
  PPGparameters;
  ReconstructionOptions;
} MethodClass;


/****************************************************************/
/*	E N D   O F   F I L E					*/
/****************************************************************/




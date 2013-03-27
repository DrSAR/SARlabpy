/****************************************************************
 *
 * $Source: /pv/CvsTree/pv/gen/src/prg/methods/FLASH/callbackDefs.h,v $
 *
 * Copyright (c) 1999-2002
 * Bruker BioSpin MRI GmbH
 * D-76275 Ettlingen, Germany
 *
 * All Rights Reserved
 *
 *
 * $Id: callbackDefs.h,v 1.20.2.2 2009/08/14 14:25:05 sako Exp $
 *
 ****************************************************************/

/* Encoding */
relations PVM_EncodingHandler backbone;

/* inplane geometry */
relations PVM_InplaneGeometryHandler  backbone;

/* slice geometry: */
relations PVM_SliceGeometryHandler  localHandleSliceGeometry;

/* digitizer parameters and bandwidth */
relations PVM_DigHandler       backbone;
relations PVM_EffSWh           EffSWhRel;

/* sequence atoms */
relations PVM_3dPhaseHandler          backbone;
relations PVM_2dPhaseHandler          backbone;
relations PVM_FreqEncodingHandler     backbone;
relations PVM_ExSliceSelectionHandler backbone;

/* modules */
relations PVM_FatSupHandler     backbone;
relations PVM_MagTransHandler   backbone;
relations PVM_SatSlicesHandler  backbone;
relations PVM_FlowSatHandler    backbone;
relations PVM_TriggerHandler    backbone;
relations PVM_TaggingHandler    backbone;
relations PVM_EvolutionHandler  backbone;
relations PVM_SelIrHandler      backbone;
relations PVM_BlBloodHandler    backbone;

/* other parameters */
relations PVM_NucleiHandler     backbone;
relations PVM_DeriveGains       backbone;
relations PVM_RepetitionTime    backbone;
relations PVM_EchoTime          backbone;
relations PVM_EchoTime1         LocalEchoTime1Relation;
relations PVM_NAverages         Local_NAveragesHandler;
relations PVM_MotionSupOnOff    backbone;
relations PVM_ExcPulseAngle     ExcPulseAngleRelation;
relations PVM_NMovieFrames      Local_NMovieFramesRels;
relations PVM_MovieOnOff        Local_MovieOnOffRels;
relations PVM_NRepetitions      backbone;
relations PVM_InversionTime     localInversionRel;

/* react on parameter adjustments */
relations PVM_AdjResultHandler backbone;
/*
 * Redirect relation for visu creation
 */
relations VisuDerivePars        deriveVisu;

/* redirection of method specific reconstruction method */
relations RecoUserUpdate        RecoUserUpdateRel;

/****************************************************************/
/*	E N D   O F   F I L E					*/
/****************************************************************/







